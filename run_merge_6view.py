"""Phase 4d: BRCA 6-view multi-omic merge.

Reads the 6 BRCA per-modality long-format Parquets from data/TCGA-BRCA/,
applies per-view feature selection (top-5000 variance for RNA + methylation,
all of CNV/RPPA/miRNA, ≥2% recurrence for mutations), pivots each to wide,
6-way inner-joins on patient_id, writes a single merged Parquet.

Output: data/TCGA-BRCA/merged_brca_6view.parquet
"""
from pathlib import Path
from src.merge import merge_brca_6view

out = merge_brca_6view(Path("data"), Path("data/TCGA-BRCA"))
print(f"BRCA 6-view merged matrix written to {out}")
