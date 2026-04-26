"""RNA-seq ingestion: read directly from S3 via DuckDB, write Parquet."""
import csv
import pathlib

from src.gdc_client import TCGA_S3_BUCKET, fetch_manifest, get_duckdb_conn
from src.utils import get_logger

logger = get_logger(__name__)

_DATA_CATEGORY = "Transcriptome Profiling"
_DATA_TYPE = "Gene Expression Quantification"


def ingest_rnaseq(
    output_dir: pathlib.Path,
    project_id: str = "TCGA-BRCA",
) -> pathlib.Path:
    """
    Fetch GDC manifest, read each patient's RNA-seq TSV directly from S3 into
    DuckDB, and write a single consolidated Parquet file.

    Output schema (long format):
        patient_id:      VARCHAR  — 12-char TCGA barcode
        gene_id:         VARCHAR  — Ensembl ID with version
        fpkm_unstranded: DOUBLE

    Args:
        output_dir: Directory where rna_seq.parquet (and errors_rnaseq.csv) are written.
        project_id: GDC project ID (default "TCGA-BRCA").

    Returns:
        pathlib.Path to the written rna_seq.parquet.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Fetching RNA-seq manifest for %s", project_id)
    manifest = fetch_manifest(project_id, _DATA_CATEGORY, _DATA_TYPE)
    logger.info("Manifest contains %d files", len(manifest))

    con = get_duckdb_conn()
    con.execute("""
        CREATE TABLE rnaseq_staging (
            patient_id      VARCHAR,
            gene_id         VARCHAR,
            fpkm_unstranded DOUBLE
        )
    """)

    errors: list[dict] = []

    for entry in manifest:
        s3_path = f"{TCGA_S3_BUCKET}rnaseq/{entry['file_id']}/{entry['file_name']}"
        try:
            con.execute("""
                INSERT INTO rnaseq_staging
                SELECT ? AS patient_id, gene_id, fpkm_unstranded::DOUBLE
                FROM read_csv(?, delim='\t', header=true)
                WHERE NOT starts_with(gene_id, 'N_')
            """, [entry["patient_id"], s3_path])
        except Exception as exc:  # noqa: BLE001
            logger.warning("Skipping %s (file_id=%s): %s", entry["patient_id"], entry["file_id"], exc)
            errors.append({"patient_id": entry["patient_id"], "file_id": entry["file_id"], "reason": str(exc)})

    row_count = con.execute("SELECT COUNT(*) FROM rnaseq_staging").fetchone()[0]
    if row_count == 0:
        raise RuntimeError(f"No rows ingested for {project_id}. Check errors_rnaseq.csv.")

    parquet_path = output_dir / "rna_seq.parquet"
    con.execute(f"COPY rnaseq_staging TO '{parquet_path}' (FORMAT PARQUET, COMPRESSION SNAPPY)")
    logger.info("Wrote %d rows to %s", row_count, parquet_path)

    if errors:
        errors_path = output_dir / "errors_rnaseq.csv"
        with open(errors_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["patient_id", "file_id", "reason"])
            writer.writeheader()
            writer.writerows(errors)
        logger.warning("Wrote %d errors to %s", len(errors), errors_path)

    return parquet_path
