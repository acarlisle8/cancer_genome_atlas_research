"""Multi-cohort multi-modality ingest orchestrator using polars.

Loops 3 cohorts × 3 modalities, writes data/{cohort}/{modality}.parquet,
skips files that already exist so re-runs are resumable.
"""
from pathlib import Path

from src.ingest_cnv_polars import ingest_cnv_polars
from src.ingest_methylation_polars import ingest_methylation_polars
from src.ingest_rnaseq_polars import ingest_rnaseq_polars

DATA_DIR = Path("data")
COHORTS = ["TCGA-BRCA", "TCGA-LUAD", "TCGA-PRAD"]

# (output filename stem, ingest function)
INGESTORS = [
    ("rna_seq",     ingest_rnaseq_polars),
    ("cnv",         ingest_cnv_polars),
    ("methylation", ingest_methylation_polars),
]


for cohort in COHORTS:
    cohort_dir = DATA_DIR / cohort
    for modality, ingest_fn in INGESTORS:
        out_path = cohort_dir / f"{modality}.parquet"
        if out_path.exists():
            print(f"[skip] {cohort} {modality} — already exists at {out_path}")
            continue
        print(f"=== {cohort} {modality} ===")
        try:
            ingest_fn(cohort_dir, project_id=cohort)
        except Exception as exc:
            print(f"[FAIL] {cohort} {modality}: {exc}")
