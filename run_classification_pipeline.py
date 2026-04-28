"""Train and evaluate a reproducible XGBoost cohort classifier.

This script starts from the split-ready table produced by run_preprocess.py.
It avoids leakage by splitting before train-only feature selection and by
excluding patient_id/cohort/cohort_code from the feature matrix.
"""

from __future__ import annotations

import argparse
import json
import os
import pathlib
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "data/classification_pipeline/.matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "data/classification_pipeline/.cache")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
os.makedirs(os.environ["XDG_CACHE_HOME"], exist_ok=True)

import numpy as np
import pandas as pd
import polars as pl
import sklearn
import xgboost
import matplotlib.pyplot as plt
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    classification_report,
    confusion_matrix,
    ConfusionMatrixDisplay,
    f1_score,
    log_loss,
    RocCurveDisplay,
)
from sklearn.model_selection import train_test_split
from sklearn.utils.class_weight import compute_sample_weight
from xgboost import XGBClassifier


RANDOM_STATE = 42
ID_COL = "patient_id"
TARGET_COL = "cohort"
TARGET_CODE_COL = "cohort_code"


def save_evaluation_plots(
    *,
    y_test: pd.Series,
    predictions: np.ndarray,
    probabilities: np.ndarray,
    class_names: list[str],
    importance_frame: pd.DataFrame,
    output_dir: pathlib.Path,
) -> dict[str, str]:
    """Write held-out evaluation plots and return their paths."""
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    paths = {
        "confusion_matrix": plots_dir / "confusion_matrix.png",
        "normalized_confusion_matrix": plots_dir / "confusion_matrix_normalized.png",
        "feature_importance": plots_dir / "feature_importance_top25.png",
        "roc_curves": plots_dir / "roc_curves.png",
    }

    fig, ax = plt.subplots(figsize=(7, 6))
    ConfusionMatrixDisplay.from_predictions(
        y_test,
        predictions,
        display_labels=class_names,
        cmap="Blues",
        colorbar=False,
        ax=ax,
    )
    ax.set_title("Confusion Matrix")
    fig.tight_layout()
    fig.savefig(paths["confusion_matrix"], dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6))
    ConfusionMatrixDisplay.from_predictions(
        y_test,
        predictions,
        display_labels=class_names,
        normalize="true",
        cmap="Blues",
        colorbar=False,
        values_format=".2f",
        ax=ax,
    )
    ax.set_title("Normalized Confusion Matrix")
    fig.tight_layout()
    fig.savefig(paths["normalized_confusion_matrix"], dpi=180)
    plt.close(fig)

    top_importance = importance_frame.head(25).sort_values("importance")
    fig_height = max(6, 0.28 * len(top_importance))
    fig, ax = plt.subplots(figsize=(10, fig_height))
    ax.barh(top_importance["feature"], top_importance["importance"], color="#2f6f9f")
    ax.set_title("Top 25 XGBoost Feature Importances")
    ax.set_xlabel("Importance")
    ax.set_ylabel("")
    fig.tight_layout()
    fig.savefig(paths["feature_importance"], dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 6))
    for idx, class_name in enumerate(class_names):
        binary_truth = (y_test.to_numpy() == idx).astype(int)
        RocCurveDisplay.from_predictions(
            binary_truth,
            probabilities[:, idx],
            name=class_name,
            ax=ax,
        )
    ax.plot([0, 1], [0, 1], linestyle="--", color="0.65", linewidth=1)
    ax.set_title("One-vs-Rest ROC Curves")
    fig.tight_layout()
    fig.savefig(paths["roc_curves"], dpi=180)
    plt.close(fig)

    return {name: str(path) for name, path in paths.items()}


def load_manifest(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def select_top_variance_features(
    X_train: pd.DataFrame,
    feature_cols: list[str],
    top_n_features: int | None,
) -> list[str]:
    """Select top-N variance features using training data only."""
    if top_n_features is None or top_n_features <= 0 or top_n_features >= len(feature_cols):
        return feature_cols

    variances = X_train[feature_cols].var(axis=0, skipna=True)
    variances = variances.replace([np.inf, -np.inf], np.nan).fillna(-np.inf)
    return variances.sort_values(ascending=False).head(top_n_features).index.to_list()


def train_and_evaluate(
    data_path: pathlib.Path,
    manifest_path: pathlib.Path,
    output_dir: pathlib.Path,
    test_size: float,
    top_n_features: int | None,
    random_state: int,
) -> dict[str, Any]:
    manifest = load_manifest(manifest_path)
    label_map = manifest["label_map"]
    code_to_label = {int(code): label for label, code in label_map.items()}
    class_names = [code_to_label[i] for i in sorted(code_to_label)]

    df = pl.read_parquet(data_path).to_pandas()
    feature_cols = manifest["feature_columns"]
    forbidden = {ID_COL, TARGET_COL, TARGET_CODE_COL}
    leaked_features = sorted(forbidden.intersection(feature_cols))
    if leaked_features:
        raise ValueError(f"Forbidden columns found in feature list: {leaked_features}")

    X = df[feature_cols]
    y = df[TARGET_CODE_COL].astype(int)

    (
        X_train,
        X_test,
        y_train,
        y_test,
        train_ids,
        test_ids,
        train_labels,
        test_labels,
    ) = train_test_split(
        X,
        y,
        df[ID_COL],
        df[TARGET_COL],
        test_size=test_size,
        stratify=y,
        random_state=random_state,
    )

    selected_features = select_top_variance_features(X_train, feature_cols, top_n_features)
    X_train_selected = X_train[selected_features]
    X_test_selected = X_test[selected_features]

    sample_weight = compute_sample_weight(class_weight="balanced", y=y_train)
    model = XGBClassifier(
        objective="multi:softprob",
        num_class=len(class_names),
        eval_metric="mlogloss",
        n_estimators=400,
        max_depth=4,
        learning_rate=0.05,
        subsample=0.9,
        colsample_bytree=0.8,
        reg_lambda=1.0,
        random_state=random_state,
        n_jobs=4,
        missing=np.nan,
        tree_method="hist",
    )
    model.fit(X_train_selected, y_train, sample_weight=sample_weight)

    probabilities = model.predict_proba(X_test_selected)
    predictions = np.asarray(probabilities).argmax(axis=1)

    importance_frame = pd.DataFrame(
        {
            "feature": selected_features,
            "importance": model.feature_importances_,
        }
    ).sort_values("importance", ascending=False)
    plot_paths = save_evaluation_plots(
        y_test=y_test,
        predictions=predictions,
        probabilities=probabilities,
        class_names=class_names,
        importance_frame=importance_frame,
        output_dir=output_dir,
    )

    metrics: dict[str, Any] = {
        "data_path": str(data_path),
        "manifest_path": str(manifest_path),
        "package_versions": {
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "polars": pl.__version__,
            "scikit_learn": sklearn.__version__,
            "xgboost": xgboost.__version__,
        },
        "random_state": random_state,
        "test_size": test_size,
        "top_n_features": top_n_features,
        "n_features_available": len(feature_cols),
        "n_features_selected": len(selected_features),
        "label_map": label_map,
        "train_class_counts": {
            code_to_label[int(code)]: int(count)
            for code, count in y_train.value_counts().sort_index().items()
        },
        "test_class_counts": {
            code_to_label[int(code)]: int(count)
            for code, count in y_test.value_counts().sort_index().items()
        },
        "accuracy": float(accuracy_score(y_test, predictions)),
        "balanced_accuracy": float(balanced_accuracy_score(y_test, predictions)),
        "macro_f1": float(f1_score(y_test, predictions, average="macro")),
        "weighted_f1": float(f1_score(y_test, predictions, average="weighted")),
        "log_loss": float(log_loss(y_test, probabilities, labels=sorted(code_to_label))),
        "classification_report": classification_report(
            y_test,
            predictions,
            target_names=class_names,
            output_dict=True,
            zero_division=0,
        ),
        "confusion_matrix": confusion_matrix(y_test, predictions).tolist(),
        "plot_paths": plot_paths,
        "leakage_notes": [
            "Train/test split is stratified and happens before variance feature selection.",
            "Top-variance feature selection is fit on X_train only.",
            "No imputation or scaling is fit before splitting.",
            "XGBoost handles NaN feature values natively.",
            "Balanced sample weights are fit on y_train only.",
            "The input merged table was already globally feature-selected upstream; rebuild from raw modality tables for the strictest final evaluation.",
        ],
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    model_path = output_dir / "xgboost_cohort_model.json"
    metrics_path = output_dir / "classification_metrics.json"
    features_path = output_dir / "selected_features.json"
    predictions_path = output_dir / "test_predictions.csv"
    split_path = output_dir / "split_assignments.csv"
    importances_path = output_dir / "feature_importances.csv"

    model.save_model(model_path)
    metrics["model_path"] = str(model_path)
    metrics["metrics_path"] = str(metrics_path)
    metrics["selected_features_path"] = str(features_path)
    metrics["predictions_path"] = str(predictions_path)
    metrics["split_assignments_path"] = str(split_path)
    metrics["feature_importances_path"] = str(importances_path)

    features_path.write_text(json.dumps(selected_features, indent=2) + "\n")
    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n")

    prediction_frame = pd.DataFrame(
        {
            ID_COL: test_ids.reset_index(drop=True),
            "actual_code": y_test.reset_index(drop=True),
            "actual_cohort": test_labels.reset_index(drop=True),
            "predicted_code": predictions,
            "predicted_cohort": [code_to_label[int(code)] for code in predictions],
        }
    )
    for idx, class_name in enumerate(class_names):
        prediction_frame[f"prob_{class_name}"] = probabilities[:, idx]
    prediction_frame.to_csv(predictions_path, index=False)

    split_frame = pd.concat(
        [
            pd.DataFrame(
                {
                    ID_COL: train_ids.reset_index(drop=True),
                    TARGET_COL: train_labels.reset_index(drop=True),
                    "split": "train",
                }
            ),
            pd.DataFrame(
                {
                    ID_COL: test_ids.reset_index(drop=True),
                    TARGET_COL: test_labels.reset_index(drop=True),
                    "split": "test",
                }
            ),
        ],
        ignore_index=True,
    )
    split_frame.to_csv(split_path, index=False)

    importance_frame.to_csv(importances_path, index=False)

    return metrics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data",
        type=pathlib.Path,
        default=pathlib.Path("data/model_ready_cohort.parquet"),
    )
    parser.add_argument(
        "--manifest",
        type=pathlib.Path,
        default=pathlib.Path("data/model_ready_cohort_manifest.json"),
    )
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=pathlib.Path("data/classification_pipeline"),
    )
    parser.add_argument("--test-size", type=float, default=0.20)
    parser.add_argument(
        "--top-n-features",
        type=int,
        default=300,
        help="Select top-N variance features using training data only. Use 0 to keep all features.",
    )
    parser.add_argument("--random-state", type=int, default=RANDOM_STATE)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    top_n = None if args.top_n_features <= 0 else args.top_n_features
    metrics = train_and_evaluate(
        data_path=args.data,
        manifest_path=args.manifest,
        output_dir=args.output_dir,
        test_size=args.test_size,
        top_n_features=top_n,
        random_state=args.random_state,
    )

    print("Classification pipeline complete")
    print(f"Accuracy: {metrics['accuracy']:.4f}")
    print(f"Balanced accuracy: {metrics['balanced_accuracy']:.4f}")
    print(f"Macro F1: {metrics['macro_f1']:.4f}")
    print(f"Log loss: {metrics['log_loss']:.4f}")
    print(f"Selected features: {metrics['n_features_selected']} of {metrics['n_features_available']}")
    print(f"Metrics written to {metrics['metrics_path']}")


if __name__ == "__main__":
    main()
