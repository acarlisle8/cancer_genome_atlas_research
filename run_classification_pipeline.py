"""Train and evaluate a reproducible XGBoost classifier.

This script starts from the split-ready table produced by run_preprocess.py or
run_preprocess_subtype.py. It uses stratified 5-fold CV and excludes
patient_id and target columns from the feature matrix to prevent leakage.
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
from sklearn.model_selection import StratifiedKFold
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


def train_and_evaluate(
    data_path: pathlib.Path,
    manifest_path: pathlib.Path,
    output_dir: pathlib.Path,
    random_state: int,
) -> dict[str, Any]:
    manifest = load_manifest(manifest_path)
    label_map = manifest["label_map"]
    code_to_label = {int(code): label for label, code in label_map.items()}
    class_names = [code_to_label[i] for i in sorted(code_to_label)]

    target_col = manifest.get("target_column", TARGET_COL)
    target_code_col = manifest.get("target_code_column", TARGET_CODE_COL)

    df = pl.read_parquet(data_path).to_pandas()
    feature_cols = manifest["feature_columns"]
    forbidden = {ID_COL, target_col, target_code_col}
    leaked_features = sorted(forbidden.intersection(feature_cols))
    if leaked_features:
        raise ValueError(f"Forbidden columns found in feature list: {leaked_features}")

    X = df[feature_cols]
    y = df[target_code_col].astype(int)
    patient_ids = df[ID_COL]
    target_labels = df[target_col]

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state)

    fold_accuracies: list[float] = []
    fold_balanced_accuracies: list[float] = []
    fold_macro_f1s: list[float] = []
    fold_log_losses: list[float] = []
    oof_rows: list[dict[str, Any]] = []

    for fold_idx, (train_idx, val_idx) in enumerate(skf.split(X, y)):
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]

        sample_weight = compute_sample_weight(class_weight="balanced", y=y_train)
        fold_model = XGBClassifier(
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
        fold_model.fit(X_train, y_train, sample_weight=sample_weight)

        fold_probs = fold_model.predict_proba(X_val)
        fold_preds = np.asarray(fold_probs).argmax(axis=1)

        fold_accuracies.append(float(accuracy_score(y_val, fold_preds)))
        fold_balanced_accuracies.append(float(balanced_accuracy_score(y_val, fold_preds)))
        fold_macro_f1s.append(float(f1_score(y_val, fold_preds, average="macro", zero_division=0)))
        fold_log_losses.append(float(log_loss(y_val, fold_probs, labels=sorted(code_to_label))))

        for i, orig_idx in enumerate(val_idx):
            prob_dict = {f"prob_{class_names[j]}": float(fold_probs[i, j]) for j in range(len(class_names))}
            oof_rows.append({
                ID_COL: patient_ids.iloc[orig_idx],
                "fold": fold_idx,
                target_col: target_labels.iloc[orig_idx],
                "actual_code": int(y_val.iloc[i]),
                f"actual_{target_col}": target_labels.iloc[orig_idx],
                f"predicted_{target_col}": code_to_label[int(fold_preds[i])],
                "predicted_code": int(fold_preds[i]),
                **prob_dict,
            })

    cv_scores = {
        "accuracy":          {"mean": float(np.mean(fold_accuracies)),          "std": float(np.std(fold_accuracies))},
        "balanced_accuracy": {"mean": float(np.mean(fold_balanced_accuracies)), "std": float(np.std(fold_balanced_accuracies))},
        "macro_f1":          {"mean": float(np.mean(fold_macro_f1s)),           "std": float(np.std(fold_macro_f1s))},
        "log_loss":          {"mean": float(np.mean(fold_log_losses)),          "std": float(np.std(fold_log_losses))},
    }

    full_sample_weight = compute_sample_weight(class_weight="balanced", y=y)
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
    model.fit(X, y, sample_weight=full_sample_weight)

    importance_frame = pd.DataFrame(
        {
            "feature": feature_cols,
            "importance": model.feature_importances_,
        }
    ).sort_values("importance", ascending=False)

    oof_frame = pd.DataFrame(oof_rows).sort_values(ID_COL).reset_index(drop=True)
    oof_y_true = oof_frame["actual_code"]
    oof_y_pred = oof_frame["predicted_code"]
    oof_probs = oof_frame[[f"prob_{c}" for c in class_names]].to_numpy()

    plot_paths = save_evaluation_plots(
        y_test=oof_y_true,
        predictions=oof_y_pred.to_numpy(),
        probabilities=oof_probs,
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
        "n_folds": 5,
        "cv_scores": cv_scores,
        "n_features": len(feature_cols),
        "label_map": label_map,
        "class_counts": {
            code_to_label[int(code)]: int(count)
            for code, count in y.value_counts().sort_index().items()
        },
        "accuracy": cv_scores["accuracy"]["mean"],
        "balanced_accuracy": cv_scores["balanced_accuracy"]["mean"],
        "macro_f1": cv_scores["macro_f1"]["mean"],
        "weighted_f1": float(f1_score(oof_y_true, oof_y_pred, average="weighted", zero_division=0)),
        "log_loss": cv_scores["log_loss"]["mean"],
        "classification_report": classification_report(
            oof_y_true,
            oof_y_pred,
            target_names=class_names,
            output_dict=True,
            zero_division=0,
        ),
        "confusion_matrix": confusion_matrix(oof_y_true, oof_y_pred).tolist(),
        "plot_paths": plot_paths,
        "leakage_notes": [
            "Stratified 5-fold CV — no patient appears in both train and val for their own fold.",
            "No imputation or scaling is fit before splitting.",
            "XGBoost handles NaN feature values natively.",
            "Balanced sample weights are fit on training fold only, never on val.",
            "Final model is refit on full dataset for feature importances and the model artifact.",
            "OOF predictions cover all patients and are used for confusion matrix, ROC curves, and test_predictions.csv.",
        ],
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    model_path = output_dir / f"xgboost_{target_col}_model.json"
    metrics_path = output_dir / "classification_metrics.json"
    predictions_path = output_dir / "test_predictions.csv"
    split_path = output_dir / "split_assignments.csv"
    importances_path = output_dir / "feature_importances.csv"

    model.save_model(model_path)
    metrics["model_path"] = str(model_path)
    metrics["metrics_path"] = str(metrics_path)
    metrics["predictions_path"] = str(predictions_path)
    metrics["split_assignments_path"] = str(split_path)
    metrics["feature_importances_path"] = str(importances_path)

    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n")

    pred_cols = [ID_COL, "fold", f"actual_{target_col}", "actual_code", f"predicted_{target_col}", "predicted_code"] + [f"prob_{c}" for c in class_names]
    oof_frame[pred_cols].to_csv(predictions_path, index=False)

    split_frame = oof_frame[[ID_COL, target_col, "fold"]].copy()
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
    parser.add_argument("--random-state", type=int, default=RANDOM_STATE)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    metrics = train_and_evaluate(
        data_path=args.data,
        manifest_path=args.manifest,
        output_dir=args.output_dir,
        random_state=args.random_state,
    )

    print("Classification pipeline complete")
    print(f"CV accuracy:          {metrics['cv_scores']['accuracy']['mean']:.4f} ± {metrics['cv_scores']['accuracy']['std']:.4f}")
    print(f"CV balanced accuracy: {metrics['cv_scores']['balanced_accuracy']['mean']:.4f} ± {metrics['cv_scores']['balanced_accuracy']['std']:.4f}")
    print(f"CV macro F1:          {metrics['cv_scores']['macro_f1']['mean']:.4f} ± {metrics['cv_scores']['macro_f1']['std']:.4f}")
    print(f"CV log loss:          {metrics['cv_scores']['log_loss']['mean']:.4f} ± {metrics['cv_scores']['log_loss']['std']:.4f}")
    print(f"Features used: {metrics['n_features']}")
    print(f"Metrics written to {metrics['metrics_path']}")


if __name__ == "__main__":
    main()
