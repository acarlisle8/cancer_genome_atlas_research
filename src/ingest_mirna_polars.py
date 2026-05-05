"""miRNA ingestion (Polars variant): use pl.scan_csv to read all patient
miRBase quantification TSVs from S3 lazily, then sink_parquet for streaming
output.

Mirrors src/ingest_cnv_polars.py — single batch scan via object_store, join
to manifest by file_id extracted from the source path, stream to disk.

miRNA files have header:
    miRNA_ID  read_count  reads_per_million_miRNA_mapped  cross-mapped
1881 miRNAs per patient (full miRBase v21 panel). ~70% of values are 0
(no expression detected) — these are real signal, not missing.
"""
import pathlib
import time

import polars as pl

from src.gdc_client import TCGA_S3_BUCKET, fetch_manifest
from src.utils import get_logger

logger = get_logger(__name__)

_DATA_CATEGORY = "Transcriptome Profiling"
_DATA_TYPE = "miRNA Expression Quantification"


def ingest_mirna_polars(
    output_dir: pathlib.Path,
    project_id: str = "TCGA-BRCA",
) -> pathlib.Path:
    """
    Fetch GDC miRNA manifest, scan all patient miRBase quantification TSVs
    from S3 with Polars, join to manifest by file_id, stream to a single
    Parquet via sink_parquet.

    Output schema (long format):
        patient_id:                     VARCHAR  — 12-char TCGA barcode
        mirna_id:                       VARCHAR  — miRBase ID (e.g. "hsa-let-7a-1")
        reads_per_million_mirna_mapped: FLOAT64  — RPM normalized count
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Fetching miRNA manifest for %s", project_id)
    manifest = fetch_manifest(project_id, _DATA_CATEGORY, _DATA_TYPE)
    logger.info("Manifest contains %d files", len(manifest))

    file_paths = [f"{TCGA_S3_BUCKET}mirna/{e['file_id']}/{e['file_name']}" for e in manifest]
    manifest_lf = pl.LazyFrame(manifest)

    logger.info("Scanning %d files via Polars scan_csv...", len(file_paths))
    t0 = time.time()

    lf = pl.scan_csv(
        file_paths,
        separator="\t",
        has_header=True,
        schema={
            "miRNA_ID":                       pl.String,
            "read_count":                     pl.Int64,
            "reads_per_million_miRNA_mapped": pl.Float64,
            "cross-mapped":                   pl.String,
        },
        include_file_paths="_source_file",
        storage_options={"aws_region": "us-east-1", "skip_signature": "true"},
    )

    result = (
        lf
        .with_columns(
            pl.col("_source_file").str.extract(r"/([0-9a-f-]{36})/", 1).alias("file_id")
        )
        .join(manifest_lf, on="file_id", how="inner")
        .select([
            "patient_id",
            pl.col("miRNA_ID").alias("mirna_id"),
            pl.col("reads_per_million_miRNA_mapped").alias("reads_per_million_mirna_mapped"),
        ])
    )

    parquet_path = output_dir / "mirna.parquet"
    result.sink_parquet(parquet_path, compression="snappy")

    logger.info("Wrote to %s in %.1fs total", parquet_path, time.time() - t0)
    return parquet_path
