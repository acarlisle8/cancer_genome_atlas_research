"""Methylation ingestion (Polars variant): use pl.scan_csv to read all patient
beta-value TSVs from S3 lazily, then sink_parquet for streaming output.

Mirrors src/ingest_rnaseq_polars.py and src/ingest_cnv_polars.py — single
batch scan via object_store, join to manifest by file_id from the source path,
stream to disk via sink_parquet (~850K probes × ~300 patients per cohort).

Methylation files are headerless TSVs with two columns: probe_id, beta_value.
"""
import pathlib
import time

import polars as pl

from src.gdc_client import TCGA_S3_BUCKET, fetch_manifest
from src.utils import get_logger

logger = get_logger(__name__)

_DATA_CATEGORY = "DNA Methylation"
_DATA_TYPE = "Methylation Beta Value"


def ingest_methylation_polars(
    output_dir: pathlib.Path,
    project_id: str = "TCGA-BRCA",
) -> pathlib.Path:
    """
    Fetch GDC methylation manifest, scan all patient beta-value TSVs from S3
    with Polars, join to manifest by file_id, stream to a single Parquet via
    sink_parquet.

    Output schema (long format):
        patient_id: VARCHAR  — 12-char TCGA barcode
        probe_id:   VARCHAR  — methylation probe identifier (e.g. "cg00000029")
        beta_value: FLOAT64  — beta value in [0.0, 1.0]; NA in source becomes null
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Fetching methylation manifest for %s", project_id)
    manifest = fetch_manifest(project_id, _DATA_CATEGORY, _DATA_TYPE)
    logger.info("Manifest contains %d files", len(manifest))

    file_paths = [f"{TCGA_S3_BUCKET}methylation/{e['file_id']}/{e['file_name']}" for e in manifest]
    manifest_lf = pl.LazyFrame(manifest)

    logger.info("Scanning %d files via Polars scan_csv...", len(file_paths))
    t0 = time.time()

    # Headerless TSV; "NA" treated as null so beta_value parses as Float64.
    # Polars scan_csv with has_header=False fails on s3:// URLs (treats them as
    # local paths). Workaround: has_header=True + schema override consumes the
    # first probe row as a "header"; losing 1 of ~850K probes per file is
    # irrelevant when we select top-500 by variance downstream.
    lf = pl.scan_csv(
        file_paths,
        separator="\t",
        has_header=True,
        schema={
            "probe_id":   pl.String,
            "beta_value": pl.Float64,
        },
        null_values=["NA"],
        include_file_paths="_source_file",
        storage_options={"aws_region": "us-east-1"},
    )

    result = (
        lf
        .with_columns(
            pl.col("_source_file").str.extract(r"/([0-9a-f-]{36})/", 1).alias("file_id")
        )
        .join(manifest_lf, on="file_id", how="inner")
        .select([
            "patient_id",
            "probe_id",
            "beta_value",
        ])
    )

    parquet_path = output_dir / "methylation.parquet"
    result.sink_parquet(parquet_path, compression="snappy")

    logger.info("Wrote to %s in %.1fs total", parquet_path, time.time() - t0)
    return parquet_path
