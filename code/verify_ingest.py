"""Compare per-cohort patient counts: GDC manifest (access=open) vs Parquet on disk.

For each (cohort, modality), fetches the open-access manifest and counts
unique patient_ids in the corresponding Parquet via DuckDB streaming
COUNT(DISTINCT). Mismatches indicate hidden sync gaps that didn't surface
as ingest crashes (e.g., missing patients silently dropped because the
missing file came late enough in the scan).

Uses DuckDB rather than polars n_unique() because polars materializes the
entire column into a hash set in memory; for methylation parquets with
~580M rows this OOMs even small EC2 instances. DuckDB streams.
"""
import gc
from pathlib import Path

import duckdb

from src.gdc_client import fetch_manifest

DATA_DIR = Path("data")

COHORTS = ["TCGA-BRCA", "TCGA-LUAD", "TCGA-PRAD"]

MODALITIES = [
    ("rna_seq",     "Transcriptome Profiling", "Gene Expression Quantification"),
    ("cnv",         "Copy Number Variation",   "Masked Copy Number Segment"),
    ("methylation", "DNA Methylation",         "Methylation Beta Value"),
]


def count_unique_patients(parquet_path: Path) -> int:
    """Streaming COUNT(DISTINCT patient_id) via DuckDB — no full-column materialization."""
    con = duckdb.connect(":memory:")
    try:
        # Cap memory hard so a small EC2 instance can't OOM. COUNT(DISTINCT) on
        # patient_id needs at most a few hundred unique 12-char strings — KBs.
        con.execute("PRAGMA memory_limit='1GB'")
        con.execute("PRAGMA threads=2")
        n = con.execute(
            f"SELECT COUNT(DISTINCT patient_id) FROM read_parquet('{parquet_path}')"
        ).fetchone()[0]
    finally:
        con.close()
    return n


def main():
    print(f"{'cohort':<11} {'modality':<13} {'manifest':>10} {'parquet':>10} {'diff':>8}")
    print("-" * 60)
    for cohort in COHORTS:
        for modality, dc, dt in MODALITIES:
            manifest = fetch_manifest(cohort, dc, dt)
            manifest_patients = {e["patient_id"] for e in manifest}
            n_manifest = len(manifest_patients)

            pq = DATA_DIR / cohort / f"{modality}.parquet"
            if not pq.exists():
                print(f"{cohort:<11} {modality:<13} {n_manifest:>10} {'MISSING':>10} {'-':>8}")
                continue

            try:
                n_parquet = count_unique_patients(pq)
            except Exception as e:
                print(f"{cohort:<11} {modality:<13} {n_manifest:>10} {'ERROR':>10}  {str(e)[:30]}")
                continue

            diff = n_parquet - n_manifest
            marker = "" if diff == 0 else (" gap" if diff < 0 else " +")
            print(f"{cohort:<11} {modality:<13} {n_manifest:>10} {n_parquet:>10} {diff:>+8}{marker}")

            gc.collect()


if __name__ == "__main__":
    main()
