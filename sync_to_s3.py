"""
Copy TCGA files from the public s3://tcga-2-open/ bucket into your own S3 bucket.
Copies are server-side (boto3 copy_object) — no data travels through this machine.
Already-present destination keys are skipped.

Usage:
    uv run python sync_to_s3.py
    uv run python sync_to_s3.py --project TCGA-LUAD
"""

import argparse
import sys

import boto3
from botocore import UNSIGNED
from botocore.config import Config

from src.gdc_client import fetch_manifest

SOURCE_BUCKET = "tcga-2-open"
DEST_BUCKET = "g23861422-datsbd-s2026"
DEFAULT_PREFIX = "tcga/"

MODALITIES = [
    ("Transcriptome Profiling", "Gene Expression Quantification", "rnaseq"),
    ("Copy Number Variation",   "Masked Copy Number Segment",     "cnv"),
    ("DNA Methylation",        "Methylation Beta Value",          "methylation"),
]


def key_exists(s3_client, bucket: str, key: str) -> bool:
    try:
        s3_client.head_object(Bucket=bucket, Key=key)
        return True
    except Exception:
        return False


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
    print(f"[{modality_name}] {len(entries)} files in manifest")

    copied = skipped = 0
    for i, entry in enumerate(entries, 1):
        file_id = entry["file_id"]
        file_name = entry["file_name"]
        src_key = f"{file_id}/{file_name}"
        dest_key = f"{dest_prefix}{modality_name}/{file_id}/{file_name}"

        if key_exists(auth_s3, dest_bucket, dest_key):
            skipped += 1
            continue

        try:
            auth_s3.copy_object(
                CopySource={"Bucket": SOURCE_BUCKET, "Key": src_key},
                Bucket=dest_bucket,
                Key=dest_key,
            )
            copied += 1
            print(f"  [{i}/{len(entries)}] copied {file_name}")
        except Exception as exc:
            print(f"  [{i}/{len(entries)}] FAILED {file_name}: {exc}", file=sys.stderr)

    return copied, skipped


def main():
    parser = argparse.ArgumentParser(description="Sync TCGA data to s3://g23861422-datsbd-s2026")
    parser.add_argument("--project", default="TCGA-BRCA", help="TCGA project ID (default: TCGA-BRCA)")
    parser.add_argument("--prefix", default=DEFAULT_PREFIX, help=f"Key prefix in destination bucket (default: {DEFAULT_PREFIX})")
    args = parser.parse_args()

    # Authenticated client — reads public source bucket and writes to your bucket
    auth_s3 = boto3.client("s3", region_name="us-east-1")

    print(f"Syncing {args.project} → s3://{DEST_BUCKET}/{args.prefix}")

    total_copied = total_skipped = 0
    for data_category, data_type, modality_name in MODALITIES:
        copied, skipped = sync_modality(
            auth_s3,
            DEST_BUCKET, args.prefix,
            args.project,
            data_category, data_type, modality_name,
        )
        total_copied += copied
        total_skipped += skipped
        print(f"[{modality_name}] done — {copied} copied, {skipped} skipped")

    print(f"\nAll done. Total copied: {total_copied}, skipped: {total_skipped}")


if __name__ == "__main__":
    main()
