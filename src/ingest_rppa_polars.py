"""RPPA ingestion (Polars variant): use pl.scan_csv to read all patient
protein-expression TSVs from S3 lazily, then sink_parquet for streaming output.

Mirrors src/ingest_cnv_polars.py — single batch scan via object_store, join
to manifest by file_id extracted from the source path, stream to disk.

RPPA files have header:
    AGID  lab_id  catalog_number  set_id  peptide_target  protein_expression
~5% of antibodies are reported as the literal string "NA" per patient (panel-
level QC failures); null_values=["NA"] coerces those to nulls so the Float64
cast on protein_expression succeeds. CRLF line endings in the source are
handled by Polars without special config.
"""
import pathlib
import time

import polars as pl

from src.gdc_client import TCGA_S3_BUCKET, fetch_manifest
from src.utils import get_logger

logger = get_logger(__name__)

_DATA_CATEGORY = "Proteome Profiling"
_DATA_TYPE = "Protein Expression Quantification"


def ingest_rppa_polars(
    output_dir: pathlib.Path,
    project_id: str = "TCGA-BRCA",
) -> pathlib.Path:
    """
    Fetch GDC RPPA manifest, scan all patient protein-expression TSVs from S3
    with Polars, join to manifest by file_id, stream to a single Parquet via
    sink_parquet.

    Output schema (long format):
        patient_id:         VARCHAR  — 12-char TCGA barcode
        peptide_target:     VARCHAR  — protein/antibody target name (e.g. "1433BETA")
        protein_expression: FLOAT64  — normalized expression; "NA" → null
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Fetching RPPA manifest for %s", project_id)
    manifest = fetch_manifest(project_id, _DATA_CATEGORY, _DATA_TYPE)
    logger.info("Manifest contains %d files", len(manifest))

    file_paths = [f"{TCGA_S3_BUCKET}rppa/{e['file_id']}/{e['file_name']}" for e in manifest]
    manifest_lf = pl.LazyFrame(manifest)

    logger.info("Scanning %d files via Polars scan_csv...", len(file_paths))
    t0 = time.time()

    lf = pl.scan_csv(
        file_paths,
        separator="\t",
        has_header=True,
        schema={
            "AGID":               pl.String,
            "lab_id":             pl.String,
            "catalog_number":     pl.String,
            "set_id":             pl.String,
            "peptide_target":     pl.String,
            "protein_expression": pl.Float64,
        },
        null_values=["NA"],
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
            "peptide_target",
            "protein_expression",
        ])
    )

    parquet_path = output_dir / "rppa.parquet"
    result.sink_parquet(parquet_path, compression="snappy")

    logger.info("Wrote to %s in %.1fs total", parquet_path, time.time() - t0)
    return parquet_path
