from pathlib import Path

from src.preprocess import preprocess_for_cohort_model

INPUT_PATH = Path("data/merged_all_cohorts.parquet")
OUTPUT_PATH = Path("data/model_ready_cohort.parquet")
MANIFEST_PATH = Path("data/model_ready_cohort_manifest.json")

out = preprocess_for_cohort_model(INPUT_PATH, OUTPUT_PATH, MANIFEST_PATH)
print(f"Preprocessed cohort modeling table written to {out}")
print(f"Preprocessing manifest written to {MANIFEST_PATH}")
