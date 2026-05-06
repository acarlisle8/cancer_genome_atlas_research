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
import shap
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

from src.markers import PANELS


RANDOM_STATE = 42
ID_COL = "patient_id"
TARGET_COL = "cohort"
TARGET_CODE_COL = "cohort_code"


def find_marker_columns(feature_cols: list[str], cancer_type: str | None) -> list[str]:
    """Return columns matching the marker panel for this cancer type.

    Marker panels store bare Ensembl IDs ("ENSG00000091831"); the merged
    matrix stores versioned IDs ("ENSG00000091831.16"). Match on the bare
    prefix (split on ".").
    """
    if not cancer_type or cancer_type not in PANELS:
        return []
    panel_ids = set(PANELS[cancer_type].values())
    return [c for c in feature_cols if c.split(".", 1)[0] in panel_ids]


def select_features_for_fold(
    X_train: pd.DataFrame,
    feature_cols: list[str],
    must_keep: list[str],
    top_n: int | None,
) -> list[str]:
    """Pick top-N variance features from the training data, unioned with must_keep.

    Variance is computed on the training fold only — must_keep is domain-defined
    (marker panels) so it carries no leakage. Returns features in the original
    feature_cols order to keep downstream output deterministic.
    """
    if top_n is None or top_n <= 0 or top_n >= len(feature_cols):
        return feature_cols
    variances = X_train[feature_cols].var(axis=0, skipna=True)
    variances = variances.replace([np.inf, -np.inf], np.nan).fillna(-np.inf)
    top_picks = set(variances.sort_values(ascending=False).head(top_n).index.to_list())
    selected = top_picks | set(must_keep)
    # Preserve original feature order so importance frames are stable
    return [c for c in feature_cols if c in selected]


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


def save_shap_artifacts(
    model: XGBClassifier,
    X: pd.DataFrame,
    class_names: list[str],
    output_dir: pathlib.Path,
) -> tuple[dict[str, str], str]:
    """Compute SHAP values for the final model and save plots + a summary CSV.

    For multiclass XGBoost, shap.TreeExplainer returns an array shaped
    (n_samples, n_features, n_classes). We save:
      - shap_summary_top20.png  : beeswarm plot of top-20 features (all classes)
      - shap_per_class/*.png    : per-class beeswarm plots (top 15 per class)
      - shap_values_summary.csv : mean |SHAP| globally + per class

    Returns (plot_paths, summary_csv_path).
    """
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    per_class_dir = plots_dir / "shap_per_class"
    per_class_dir.mkdir(parents=True, exist_ok=True)

    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X)
    if isinstance(shap_values, list):
        shap_values = np.stack(shap_values, axis=-1)  # (samples, features, classes)
    elif shap_values.ndim == 2:
        shap_values = shap_values[..., None]  # binary -> (samples, features, 1)

    feature_names = list(X.columns)
    mean_abs_per_class = np.abs(shap_values).mean(axis=0)            # (features, classes)
    mean_abs_global = mean_abs_per_class.mean(axis=1)                # (features,)

    summary_df = pd.DataFrame({
        "feature": feature_names,
        "mean_abs_shap": mean_abs_global,
    })
    for i, cls in enumerate(class_names):
        summary_df[f"mean_abs_shap_{cls}"] = mean_abs_per_class[:, i]
    summary_df = summary_df.sort_values("mean_abs_shap", ascending=False)
    summary_csv = output_dir / "shap_values_summary.csv"
    summary_df.to_csv(summary_csv, index=False)

    paths: dict[str, str] = {}

    # Top-20 multi-class summary (bar of mean |SHAP| stacked across classes)
    top_20 = summary_df.head(20)["feature"].to_list()
    top_20_idx = [feature_names.index(f) for f in top_20]
    plt.figure()
    shap.summary_plot(
        [shap_values[:, top_20_idx, i] for i in range(shap_values.shape[2])],
        X[top_20].values,
        feature_names=top_20,
        class_names=class_names,
        plot_type="bar",
        show=False,
        plot_size=(10, 8),
    )
    plt.tight_layout()
    summary_path = plots_dir / "shap_summary_top20.png"
    plt.savefig(summary_path, dpi=180, bbox_inches="tight")
    plt.close("all")
    paths["shap_summary_top20"] = str(summary_path)

    # Per-class beeswarm: each class's top 15 features by |SHAP|
    for i, cls in enumerate(class_names):
        cls_importance = pd.Series(mean_abs_per_class[:, i], index=feature_names)
        top_15 = cls_importance.sort_values(ascending=False).head(15).index.to_list()
        top_15_idx = [feature_names.index(f) for f in top_15]
        plt.figure()
        shap.summary_plot(
            shap_values[:, top_15_idx, i],
            X[top_15].values,
            feature_names=top_15,
            show=False,
            plot_size=(9, 6),
        )
        plt.title(f"SHAP — {cls}")
        plt.tight_layout()
        cls_path = per_class_dir / f"shap_{cls}.png"
        plt.savefig(cls_path, dpi=180, bbox_inches="tight")
        plt.close("all")
        paths[f"shap_{cls}"] = str(cls_path)

    return paths, str(summary_csv)


def train_and_evaluate(
    data_path: pathlib.Path,
    manifest_path: pathlib.Path,
    output_dir: pathlib.Path,
    random_state: int,
    top_n_features: int | None,
) -> dict[str, Any]:
    manifest = load_manifest(manifest_path)
    label_map = manifest["label_map"]
    code_to_label = {int(code): label for label, code in label_map.items()}
    class_names = [code_to_label[i] for i in sorted(code_to_label)]

    target_col = manifest.get("target_column", TARGET_COL)
    target_code_col = manifest.get("target_code_column", TARGET_CODE_COL)
    cancer_type = manifest.get("cancer_type")  # only present for subtype tasks

    df = pl.read_parquet(data_path).to_pandas()
    feature_cols = manifest["feature_columns"]
    forbidden = {ID_COL, target_col, target_code_col}
    leaked_features = sorted(forbidden.intersection(feature_cols))
    if leaked_features:
        raise ValueError(f"Forbidden columns found in feature list: {leaked_features}")

    marker_cols = find_marker_columns(feature_cols, cancer_type)

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
    fold_feature_counts: list[int] = []

    for fold_idx, (train_idx, val_idx) in enumerate(skf.split(X, y)):
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]

        # Per-fold feature selection: top-N variance computed on training data
        # only, unioned with the marker panel (which is domain-defined and
        # leakage-free regardless of fold).
        selected = select_features_for_fold(X_train, feature_cols, marker_cols, top_n_features)
        fold_feature_counts.append(len(selected))
        X_train_sel = X_train[selected]
        X_val_sel = X_val[selected]

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
        fold_model.fit(X_train_sel, y_train, sample_weight=sample_weight)

        fold_probs = fold_model.predict_proba(X_val_sel)
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

    # Final model: same selection logic on the full dataset for the artifact +
    # feature importances. Variance computed on all rows here is fine — the CV
    # numbers above are what we report; this fit is for the saved model only.
    final_selected = select_features_for_fold(X, feature_cols, marker_cols, top_n_features)
    X_final = X[final_selected]
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
    model.fit(X_final, y, sample_weight=full_sample_weight)

    importance_frame = pd.DataFrame(
        {
            "feature": final_selected,
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

    shap_plot_paths, shap_csv_path = save_shap_artifacts(
        model=model,
        X=X_final,
        class_names=class_names,
        output_dir=output_dir,
    )
    plot_paths.update(shap_plot_paths)

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
        "n_features_pool": len(feature_cols),
        "n_features_per_fold_mean": float(np.mean(fold_feature_counts)),
        "n_features_final_model": len(final_selected),
        "n_marker_features": len(marker_cols),
        "marker_panel": cancer_type if marker_cols else None,
        "top_n_features_arg": top_n_features,
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
            "Top-N variance feature selection is computed on the training fold only, eliminating the global-feature-selection leak.",
            "Marker panel features are force-included in every fold (domain-defined, leakage-free).",
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
    metrics["shap_summary_csv"] = shap_csv_path

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
    parser.add_argument(
        "--top-n-features",
        type=int,
        default=500,
        help="Top-N variance feature selection per fold (plus marker panel). 0 keeps all features in pool.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    top_n = None if args.top_n_features <= 0 else args.top_n_features
    metrics = train_and_evaluate(
        data_path=args.data,
        manifest_path=args.manifest,
        output_dir=args.output_dir,
        random_state=args.random_state,
        top_n_features=top_n,
    )

    print("Classification pipeline complete")
    print(f"CV accuracy:          {metrics['cv_scores']['accuracy']['mean']:.4f} ± {metrics['cv_scores']['accuracy']['std']:.4f}")
    print(f"CV balanced accuracy: {metrics['cv_scores']['balanced_accuracy']['mean']:.4f} ± {metrics['cv_scores']['balanced_accuracy']['std']:.4f}")
    print(f"CV macro F1:          {metrics['cv_scores']['macro_f1']['mean']:.4f} ± {metrics['cv_scores']['macro_f1']['std']:.4f}")
    print(f"CV log loss:          {metrics['cv_scores']['log_loss']['mean']:.4f} ± {metrics['cv_scores']['log_loss']['std']:.4f}")
    print(f"Features pool: {metrics['n_features_pool']}, per-fold mean: {metrics['n_features_per_fold_mean']:.0f}, marker panel: {metrics['marker_panel']} ({metrics['n_marker_features']} features)")
    print(f"Metrics written to {metrics['metrics_path']}")


if __name__ == "__main__":
    main()
