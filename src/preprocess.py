"""Preprocess the merged TCGA cohort matrix for downstream modeling."""

from __future__ import annotations

import argparse
import json
import pathlib
from collections import Counter
from typing import Any

import polars as pl

from src.utils import get_logger

logger = get_logger(__name__)

ID_COL = "patient_id"
TARGET_COL = "cohort"
TARGET_CODE_COL = "cohort_code"


def _is_cnv_arm(column: str) -> bool:
    return column[-1:] in {"p", "q"} and (column[:-1].isdigit() or column[:-1] == "X")


def _feature_family(column: str) -> str:
    if column.startswith("ENSG"):
        return "rna"
    if column.startswith("cg"):
        return "methylation"
    if _is_cnv_arm(column):
        return "cnv"
    return "other"


def _columns_all_null_in_any_cohort(
    df: pl.DataFrame,
    feature_cols: list[str],
    cohort_col: str = TARGET_COL,
) -> list[str]:
    """Return columns that are completely missing for at least one target cohort.

    These columns are dangerous for cohort prediction because XGBoost can learn
    the cohort from structural missingness instead of molecular signal.
    """
    dropped: list[str] = []
    for column in feature_cols:
        null_summary = (
            df.group_by(cohort_col)
            .agg(
                pl.col(column).null_count().alias("nulls"),
                pl.len().alias("rows"),
            )
            .with_columns((pl.col("nulls") == pl.col("rows")).alias("all_null"))
        )
        if null_summary["all_null"].any():
            dropped.append(column)
    return dropped


def _columns_over_missingness_threshold(
    df: pl.DataFrame,
    feature_cols: list[str],
    max_missing_rate: float,
) -> list[str]:
    n_rows = len(df)
    if n_rows == 0:
        return []

    null_counts = df.select(
        [pl.col(column).null_count().alias(column) for column in feature_cols]
    ).row(0)
    return [
        column
        for column, null_count in zip(feature_cols, null_counts, strict=True)
        if (null_count / n_rows) > max_missing_rate
    ]


def _label_map(cohorts: list[str]) -> dict[str, int]:
    return {cohort: idx for idx, cohort in enumerate(sorted(cohorts))}


def preprocess_for_cohort_model(
    input_path: pathlib.Path,
    output_path: pathlib.Path,
    manifest_path: pathlib.Path,
    max_missing_rate: float = 0.20,
) -> pathlib.Path:
    """Create a split-ready cohort-classification table.

    The output keeps patient_id only for tracking, keeps cohort/cohort_code as
    targets, and records the exact model feature columns in the manifest. No
    train/test split, imputation, scaling, or supervised feature selection is
    performed here; those must be fit after splitting.
    """
    if not 0 <= max_missing_rate < 1:
        raise ValueError("max_missing_rate must be in [0, 1)")

    df = pl.read_parquet(input_path)
    required = {ID_COL, TARGET_COL}
    missing_required = required - set(df.columns)
    if missing_required:
        raise ValueError(f"Input is missing required columns: {sorted(missing_required)}")
    if df[ID_COL].n_unique() != len(df):
        raise ValueError(f"{ID_COL} must be unique before modeling")

    feature_cols = [c for c in df.columns if c not in {ID_COL, TARGET_COL}]
    non_numeric = [
        c
        for c in feature_cols
        if not df.schema[c].is_numeric()
    ]
    if non_numeric:
        raise ValueError(f"Feature columns must be numeric: {non_numeric[:10]}")

    structural_missing_cols = _columns_all_null_in_any_cohort(df, feature_cols)
    remaining_cols = [c for c in feature_cols if c not in set(structural_missing_cols)]
    sparse_cols = _columns_over_missingness_threshold(df, remaining_cols, max_missing_rate)
    model_feature_cols = [c for c in remaining_cols if c not in set(sparse_cols)]

    labels = _label_map(df[TARGET_COL].unique().to_list())
    # Build the label encoding with replace instead of a fitted encoder so it is
    # deterministic and can be reproduced exactly from the manifest.
    out_df = (
        df.select([ID_COL, TARGET_COL] + model_feature_cols)
        .with_columns(
            pl.col(TARGET_COL)
            .replace(labels)
            .cast(pl.Int64)
            .alias(TARGET_CODE_COL)
        )
        .select([ID_COL, TARGET_COL, TARGET_CODE_COL] + model_feature_cols)
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.write_parquet(output_path, compression="snappy")

    family_counts = Counter(_feature_family(c) for c in model_feature_cols)
    dropped_family_counts = Counter(_feature_family(c) for c in structural_missing_cols + sparse_cols)
    manifest: dict[str, Any] = {
        "input_path": str(input_path),
        "output_path": str(output_path),
        "rows": len(out_df),
        "columns": len(out_df.columns),
        "id_column": ID_COL,
        "target_column": TARGET_COL,
        "target_code_column": TARGET_CODE_COL,
        "label_map": labels,
        "feature_columns": model_feature_cols,
        "feature_family_counts": dict(sorted(family_counts.items())),
        "max_missing_rate": max_missing_rate,
        "dropped": {
            "structural_missing_in_any_cohort": structural_missing_cols,
            "over_missingness_threshold": sparse_cols,
            "counts_by_family": dict(sorted(dropped_family_counts.items())),
        },
        "leakage_notes": [
            "Do not include patient_id, cohort, or cohort_code in X.",
            "Split before fitting any imputer, scaler, encoder, selector, or resampler.",
            "XGBoost can handle null feature values directly; no global imputation was applied.",
            "Columns all-null within any cohort were removed to avoid cohort-specific missingness shortcuts.",
        ],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    logger.info(
        "Wrote %s rows x %s columns to %s",
        len(out_df),
        len(out_df.columns),
        output_path,
    )
    logger.info("Kept feature families: %s", dict(sorted(family_counts.items())))
    logger.info(
        "Dropped %d structural-missing and %d sparse columns",
        len(structural_missing_cols),
        len(sparse_cols),
    )
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=pathlib.Path,
        default=pathlib.Path("data/merged_all_cohorts.parquet"),
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=pathlib.Path("data/model_ready_cohort.parquet"),
    )
    parser.add_argument(
        "--manifest",
        type=pathlib.Path,
        default=pathlib.Path("data/model_ready_cohort_manifest.json"),
    )
    parser.add_argument("--max-missing-rate", type=float, default=0.20)
    args = parser.parse_args()

    preprocess_for_cohort_model(
        input_path=args.input,
        output_path=args.output,
        manifest_path=args.manifest,
        max_missing_rate=args.max_missing_rate,
    )


if __name__ == "__main__":
    main()
