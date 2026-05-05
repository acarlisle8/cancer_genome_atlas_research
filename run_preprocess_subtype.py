"""Preprocess per-cancer-type subtype classification tables for BRCA, PRAD, and LUAD."""

from __future__ import annotations

import argparse
import json
import pathlib

from src.preprocess import preprocess_for_subtype_model

CANCER_TYPES = ["BRCA", "PRAD", "LUAD"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=pathlib.Path, default=pathlib.Path("data/model_ready_cohort.parquet"))
    parser.add_argument("--labels", type=pathlib.Path, default=pathlib.Path("data/audit/subtype_label_audit.parquet"))
    parser.add_argument("--label-col", type=str, default="hoadley_subtype_selected")
    parser.add_argument("--output-dir", type=pathlib.Path, default=pathlib.Path("data"))
    parser.add_argument("--max-missing-rate", type=float, default=0.20)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    for cancer_type in CANCER_TYPES:
        ct = cancer_type.lower()
        output_path = args.output_dir / f"model_ready_{ct}.parquet"
        manifest_path = args.output_dir / f"model_ready_{ct}_manifest.json"

        out = preprocess_for_subtype_model(
            features_path=args.features,
            labels_path=args.labels,
            label_col=args.label_col,
            cancer_type=cancer_type,
            output_path=output_path,
            manifest_path=manifest_path,
            max_missing_rate=args.max_missing_rate,
        )

        manifest = json.loads(manifest_path.read_text())
        n_rows = manifest["n_labeled"]
        n_classes = manifest["n_classes"]
        print(f"{cancer_type}: {n_rows} rows, {n_classes} classes → {out}")

    print(f"\nAll three subtype tables written to {args.output_dir}/")


if __name__ == "__main__":
    main()
