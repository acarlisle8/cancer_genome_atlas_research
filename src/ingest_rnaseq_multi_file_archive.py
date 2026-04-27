"""RNA-seq ingestion (chunked multi-file): read patient TSVs in chunks of N
files via parallel DuckDB read_csv, with progress logging per chunk.

Each chunk uses the multi-file parallel read pattern (read_csv with a list of
paths). Chunking caps memory pressure and gives progress visibility — without
chunking, DuckDB plans the whole 1,231-file query at once, which can hang.
"""
import pathlib
import time

import pyarrow as pa

from src.gdc_client import TCGA_S3_BUCKET, fetch_manifest, get_duckdb_conn
from src.utils import get_logger

logger = get_logger(__name__)

_DATA_CATEGORY = "Transcriptome Profiling"
_DATA_TYPE = "Gene Expression Quantification"

# Number of files per chunk. Smaller = more frequent progress logs and lower
# memory pressure; larger = fewer round-trips, less Python overhead.
CHUNK_SIZE = 100


def ingest_rnaseq_multi_file(
    output_dir: pathlib.Path,
    project_id: str = "TCGA-BRCA",
) -> pathlib.Path:
    """
    Fetch GDC manifest, read all patient RNA-seq TSVs from S3 in parallel
    chunks via DuckDB, and write a consolidated Parquet file.

    Output schema (long format):
        patient_id:         VARCHAR
        gene_id:            VARCHAR  — Ensembl ID with version
        gene_name:          VARCHAR  — HGNC symbol (e.g. TP53, BRCA1) for SHAP interpretability
        gene_type:          VARCHAR  — protein_coding / lncRNA / miRNA / etc., useful for filtering
        fpkm_unstranded:    DOUBLE   — original normalization
        tpm_unstranded:     DOUBLE   — modern cross-sample standard
        fpkm_uq_unstranded: DOUBLE   — used in published TCGA papers
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Fetching RNA-seq manifest for %s", project_id)
    manifest = fetch_manifest(project_id, _DATA_CATEGORY, _DATA_TYPE)
    logger.info("Manifest contains %d files", len(manifest))

    chunks = [manifest[i:i + CHUNK_SIZE] for i in range(0, len(manifest), CHUNK_SIZE)]
    logger.info("Split into %d chunks of up to %d files each", len(chunks), CHUNK_SIZE)

    con = get_duckdb_conn()
    con.execute("PRAGMA enable_progress_bar")
    con.execute("""
        CREATE TABLE rnaseq_staging (
            patient_id         VARCHAR,
            gene_id            VARCHAR,
            gene_name          VARCHAR,
            gene_type          VARCHAR,
            fpkm_unstranded    DOUBLE,
            tpm_unstranded     DOUBLE,
            fpkm_uq_unstranded DOUBLE
        )
    """)

    overall_start = time.time()
    for i, chunk in enumerate(chunks, 1):
        chunk_paths = [f"{TCGA_S3_BUCKET}rnaseq/{e['file_id']}/{e['file_name']}" for e in chunk]
        chunk_table = pa.Table.from_pylist(chunk)
        con.register("chunk_manifest", chunk_table)

        logger.info("Chunk %d/%d: reading %d files...", i, len(chunks), len(chunk))
        t0 = time.time()
        con.execute("""
            INSERT INTO rnaseq_staging
            SELECT
                m.patient_id,
                r.gene_id,
                r.gene_name,
                r.gene_type,
                r.fpkm_unstranded::DOUBLE     AS fpkm_unstranded,
                r.tpm_unstranded::DOUBLE      AS tpm_unstranded,
                r.fpkm_uq_unstranded::DOUBLE  AS fpkm_uq_unstranded
            FROM read_csv(?, delim='\t', header=true, filename=true) r
            JOIN chunk_manifest m
              ON m.file_id = regexp_extract(r.filename, '/([0-9a-f-]{36})/', 1)
            WHERE NOT starts_with(r.gene_id, 'N_')
        """, [chunk_paths])
        logger.info("Chunk %d/%d done in %.1fs", i, len(chunks), time.time() - t0)

    row_count = con.execute("SELECT COUNT(*) FROM rnaseq_staging").fetchone()[0]
    if row_count == 0:
        raise RuntimeError(f"No rows ingested for {project_id}.")

    parquet_path = output_dir / "rna_seq_multi_file.parquet"
    con.execute(f"COPY rnaseq_staging TO '{parquet_path}' (FORMAT PARQUET, COMPRESSION SNAPPY)")
    logger.info(
        "Wrote %d rows to %s in %.1fs total",
        row_count, parquet_path, time.time() - overall_start,
    )

    return parquet_path
