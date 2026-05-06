"""Multi-cohort multi-modality ingest orchestrator using polars.

Loops 3 cohorts × 3 modalities (rna_seq, cnv, methylation) plus 3 BRCA-only
modalities (rppa, mirna, mutations) for the Phase 4 multi-omic MOFA+ work.
Writes data/{cohort}/{modality}.parquet, skips files that already exist so
re-runs are resumable.
"""
from pathlib import Path

from src.ingest_cnv_polars import ingest_cnv_polars
from src.ingest_methylation_polars import ingest_methylation_polars
from src.ingest_mirna_polars import ingest_mirna_polars
from src.ingest_mutations import ingest_mutations
from src.ingest_rnaseq_polars import ingest_rnaseq_polars
from src.ingest_rppa_polars import ingest_rppa_polars

DATA_DIR = Path("data")
COHORTS = ["TCGA-BRCA", "TCGA-LUAD", "TCGA-PRAD"]

# (output filename stem, ingest function)
INGESTORS = [
    ("rna_seq",     ingest_rnaseq_polars),
    ("cnv",         ingest_cnv_polars),
    ("methylation", ingest_methylation_polars),
]

# Phase 4b modalities: BRCA-only multi-omic extension for MOFA+.
INGESTORS_BRCA_ONLY = [
    ("rppa",      ingest_rppa_polars),
    ("mirna",     ingest_mirna_polars),
    ("mutations", ingest_mutations),
]


for cohort in COHORTS:
    cohort_dir = DATA_DIR / cohort
    cohort_ingestors = list(INGESTORS)
    if cohort == "TCGA-BRCA":
        cohort_ingestors += INGESTORS_BRCA_ONLY
    for modality, ingest_fn in cohort_ingestors:
        out_path = cohort_dir / f"{modality}.parquet"
        if out_path.exists():
            print(f"[skip] {cohort} {modality} — already exists at {out_path}")
            continue
        print(f"=== {cohort} {modality} ===")
        try:
            ingest_fn(cohort_dir, project_id=cohort)
        except Exception as exc:
            print(f"[FAIL] {cohort} {modality}: {exc}")
