from pathlib import Path
from src.merge import merge_all_cohorts

DATA_DIR = Path("data")
OUTPUT_DIR = Path("data")

out = merge_all_cohorts(DATA_DIR, OUTPUT_DIR)
print(f"Merged matrix written to {out}")
