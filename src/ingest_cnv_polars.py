"""CNV ingestion (Polars variant): use pl.scan_csv to read all patient
segment TSVs from S3 lazily, then sink_parquet for streaming output.

Mirrors src/ingest_rnaseq_polars.py — single batch scan via object_store,
join to manifest by file_id extracted from the source path, stream to disk.
"""
import pathlib
import time

import polars as pl

from src.gdc_client import TCGA_S3_BUCKET, fetch_manifest
from src.utils import get_logger

logger = get_logger(__name__)

_DATA_CATEGORY = "Copy Number Variation"
_DATA_TYPE = "Masked Copy Number Segment"


def ingest_cnv_polars(
    output_dir: pathlib.Path,
    project_id: str = "TCGA-BRCA",
) -> pathlib.Path:
    """
    Fetch GDC CNV manifest, scan all patient segment TSVs from S3 with Polars,
    join to manifest by file_id, stream to a single Parquet via sink_parquet.

    Output schema (long format):
        patient_id:  VARCHAR  — 12-char TCGA barcode
        chromosome:  VARCHAR
        start:       INT64
        end:         INT64
        copy_number: FLOAT64  — log2 copy ratio (Segment_Mean)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Fetching CNV manifest for %s", project_id)
    manifest = fetch_manifest(project_id, _DATA_CATEGORY, _DATA_TYPE)
    logger.info("Manifest contains %d files", len(manifest))

    file_paths = [f"{TCGA_S3_BUCKET}cnv/{e['file_id']}/{e['file_name']}" for e in manifest]
    manifest_lf = pl.LazyFrame(manifest)

    logger.info("Scanning %d files via Polars scan_csv...", len(file_paths))
    t0 = time.time()

    # TCGA Masked Copy Number Segment TSVs have header:
    #   GDC_Aliquot  Chromosome  Start  End  Num_Probes  Segment_Mean
    lf = pl.scan_csv(
        file_paths,
        separator="\t",
        has_header=True,
        schema={
            "GDC_Aliquot": pl.String,
            "Chromosome":  pl.String,
            "Start":       pl.Int64,
            "End":         pl.Int64,
            "Num_Probes":  pl.Int64,
            "Segment_Mean": pl.Float64,
        },
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
            pl.col("Chromosome").alias("chromosome"),
            pl.col("Start").alias("start"),
            pl.col("End").alias("end"),
            pl.col("Segment_Mean").alias("copy_number"),
        ])
    )

    parquet_path = output_dir / "cnv.parquet"
    result.sink_parquet(parquet_path, compression="snappy")

    logger.info("Wrote to %s in %.1fs total", parquet_path, time.time() - t0)
    return parquet_path
