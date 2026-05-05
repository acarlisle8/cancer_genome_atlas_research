"""RNA-seq ingestion (Polars variant): use pl.scan_csv to read all patient
TSVs from S3 lazily, then sink_parquet for streaming output.

Polars uses the object_store Rust crate for S3, which historically has better
throughput than DuckDB's httpfs for many small CSV files. sink_parquet streams
to disk without materializing the full dataset in memory — different memory
profile than DuckDB's CREATE TABLE AS SELECT.
"""
import pathlib
import time

import polars as pl

from src.gdc_client import TCGA_S3_BUCKET, fetch_manifest
from src.utils import get_logger

logger = get_logger(__name__)

_DATA_CATEGORY = "Transcriptome Profiling"
_DATA_TYPE = "Gene Expression Quantification"


def ingest_rnaseq_polars(
    output_dir: pathlib.Path,
    project_id: str = "TCGA-BRCA",
) -> pathlib.Path:
    """
    Fetch GDC manifest, scan all patient RNA-seq TSVs from S3 with Polars,
    join to manifest by file_id, stream to a single Parquet via sink_parquet.

    Output schema (long format):
        patient_id:         VARCHAR
        gene_id:            VARCHAR
        gene_name:          VARCHAR
        gene_type:          VARCHAR
        fpkm_unstranded:    DOUBLE
        tpm_unstranded:     DOUBLE
        fpkm_uq_unstranded: DOUBLE
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Fetching RNA-seq manifest for %s", project_id)
    manifest = fetch_manifest(project_id, _DATA_CATEGORY, _DATA_TYPE)
    logger.info("Manifest contains %d files", len(manifest))

    file_paths = [f"{TCGA_S3_BUCKET}rnaseq/{e['file_id']}/{e['file_name']}" for e in manifest]
    manifest_lf = pl.LazyFrame(manifest)

    logger.info("Scanning %d files via Polars scan_csv...", len(file_paths))
    t0 = time.time()

    # TCGA STAR gene-counts files have:
    #   line 1: "# gene-model: GENCODE v36" (comment) → skip via comment_prefix
    #   line 2: real header
    #   lines 3-6: N_unmapped/N_multimapping/N_noFeature/N_ambiguous (non-numeric values
    #             in numeric columns) → keep as String at scan time, filter out, then cast
    lf = pl.scan_csv(
        file_paths,
        separator="\t",
        has_header=True,
        comment_prefix="#",
        schema={
            "gene_id":            pl.String,
            "gene_name":          pl.String,
            "gene_type":          pl.String,
            "unstranded":         pl.String,
            "stranded_first":     pl.String,
            "stranded_second":    pl.String,
            "tpm_unstranded":     pl.String,
            "fpkm_unstranded":    pl.String,
            "fpkm_uq_unstranded": pl.String,
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
        .filter(~pl.col("gene_id").str.starts_with("N_"))
        .select([
            "patient_id",
            "gene_id",
            "gene_name",
            "gene_type",
            pl.col("fpkm_unstranded").cast(pl.Float64),
            pl.col("tpm_unstranded").cast(pl.Float64),
            pl.col("fpkm_uq_unstranded").cast(pl.Float64),
        ])
    )

    parquet_path = output_dir / "rna_seq.parquet"
    result.sink_parquet(parquet_path, compression="snappy")

    logger.info("Wrote to %s in %.1fs total", parquet_path, time.time() - t0)
    return parquet_path
