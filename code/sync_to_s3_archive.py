"""
Copy TCGA files from the public s3://tcga-2-open/ bucket into your own S3 bucket.
Copies are server-side (boto3 copy_object) — no data travels through this machine.
Already-present destination keys are skipped, so re-runs are resumable.

Usage:
    uv run python sync_to_s3.py                                  # default: TCGA-BRCA
    uv run python sync_to_s3.py --project TCGA-LUAD              # single cohort
    uv run python sync_to_s3.py --project TCGA-LUAD TCGA-CHOL    # multiple cohorts
    uv run python sync_to_s3.py --all                            # all 33 TCGA cohorts
"""

import argparse
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import boto3
from botocore import UNSIGNED
from botocore.config import Config

from src.gdc_client import fetch_manifest

SOURCE_BUCKET = "tcga-2-open"
DEST_BUCKET = "g23861422-datsbd-s2026"
DEFAULT_PREFIX = "tcga/"

# Concurrency for per-file HEAD + COPY operations. Each call is a single S3
# request (~135 ms round-trip), so the original serial loop ran at ~7.5
# ops/sec — measured 2026-04-27, ~6,700 ops in 15 min across ACC/BLCA/BRCA.
# Threading parallelizes the per-file requests; with 20 workers, expected
# throughput is ~100-150 ops/sec, ~10-15× faster end-to-end.
MAX_WORKERS = 20

# Primary modalities — pulled into tcga/<name>/ for first-pass analysis.
PRIMARY_MODALITIES = [
    ("Transcriptome Profiling",     "Gene Expression Quantification",    "rnaseq"),
    ("Copy Number Variation",       "Masked Copy Number Segment",        "cnv"),
    ("DNA Methylation",             "Methylation Beta Value",            "methylation"),
    ("Simple Nucleotide Variation", "Masked Somatic Mutation",           "mutations"),
    ("Transcriptome Profiling",     "miRNA Expression Quantification",   "mirna"),
    ("Proteome Profiling",          "Protein Expression Quantification", "rppa"),
    ("Clinical",                    "Clinical Supplement",               "clinical"),
]

# Auxiliary modalities — pulled into tcga/_excluded/<name>/ to mark them as
# "do not include in first-pass analysis" (used for filtering, niche detail).
EXCLUDED_MODALITIES = [
    ("Biospecimen",             "Biospecimen Supplement",            "_excluded/biospecimen"),
    ("Transcriptome Profiling", "Isoform Expression Quantification", "_excluded/mirna_isoforms"),
]

MODALITIES = PRIMARY_MODALITIES + EXCLUDED_MODALITIES

ALL_TCGA_PROJECTS = [
    "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
    "TCGA-COAD", "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC",
    "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG",
    "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",
    "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC",
    "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM",
    "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM",
]


def key_exists(s3_client, bucket: str, key: str) -> bool:
    try:
        s3_client.head_object(Bucket=bucket, Key=key)
        return True
    except Exception:
        return False


def _process_entry(auth_s3, dest_bucket, dest_prefix, modality_name, entry):
    """Worker: HEAD-check then COPY one file. Returns 'copied', 'skipped', or 'failed'."""
    file_id = entry["file_id"]
    file_name = entry["file_name"]
    src_key = f"{file_id}/{file_name}"
    dest_key = f"{dest_prefix}{modality_name}/{file_id}/{file_name}"

    if key_exists(auth_s3, dest_bucket, dest_key):
        return "skipped"

    try:
        auth_s3.copy_object(
            CopySource={"Bucket": SOURCE_BUCKET, "Key": src_key},
            Bucket=dest_bucket,
            Key=dest_key,
        )
        return "copied"
    except Exception as exc:
        print(f"  FAILED {file_name}: {exc}", file=sys.stderr)
        return "failed"


def sync_modality(
    auth_s3,
    dest_bucket: str,
    dest_prefix: str,
    project_id: str,
    data_category: str,
    data_type: str,
    modality_name: str,
) -> tuple[int, int]:
    """Returns (copied, skipped) counts."""
    print(f"\n[{modality_name}] Fetching manifest for {project_id}...")
    entries = fetch_manifest(project_id, data_category, data_type)
    print(f"[{modality_name}] {len(entries)} files in manifest, processing with {MAX_WORKERS} workers")

    copied = skipped = failed = 0
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = [
            pool.submit(_process_entry, auth_s3, dest_bucket, dest_prefix, modality_name, e)
            for e in entries
        ]
        for fut in as_completed(futures):
            result = fut.result()
            if result == "copied":
                copied += 1
            elif result == "skipped":
                skipped += 1
            else:
                failed += 1

    if failed:
        print(f"[{modality_name}] WARNING: {failed} failures", file=sys.stderr)
    return copied, skipped


def main():
    parser = argparse.ArgumentParser(description="Sync TCGA data to s3://g23861422-datsbd-s2026")
    parser.add_argument("--project", nargs="+", help="TCGA project ID(s), space-separated. Default: TCGA-BRCA")
    parser.add_argument("--all", action="store_true", help=f"Sync all {len(ALL_TCGA_PROJECTS)} TCGA projects")
    parser.add_argument("--prefix", default=DEFAULT_PREFIX, help=f"Key prefix in destination bucket (default: {DEFAULT_PREFIX})")
    args = parser.parse_args()

    if args.all:
        projects = ALL_TCGA_PROJECTS
    elif args.project:
        projects = args.project
    else:
        projects = ["TCGA-BRCA"]

    # Authenticated client — reads public source bucket and writes to your bucket
    auth_s3 = boto3.client("s3", region_name="us-east-1")

    print(f"Syncing {len(projects)} project(s) → s3://{DEST_BUCKET}/{args.prefix}")

    grand_copied = grand_skipped = 0
    for proj in projects:
        print(f"\n========== {proj} ==========")
        proj_copied = proj_skipped = 0
        for data_category, data_type, modality_name in MODALITIES:
            copied, skipped = sync_modality(
                auth_s3,
                DEST_BUCKET, args.prefix,
                proj,
                data_category, data_type, modality_name,
            )
            proj_copied += copied
            proj_skipped += skipped
            print(f"[{modality_name}] done — {copied} copied, {skipped} skipped")
        grand_copied += proj_copied
        grand_skipped += proj_skipped
        print(f"== {proj} totals: {proj_copied} copied, {proj_skipped} skipped ==")

    print(f"\n=== ALL DONE === Grand totals: {grand_copied} copied, {grand_skipped} skipped across {len(projects)} project(s)")


if __name__ == "__main__":
    main()
