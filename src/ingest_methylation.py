"""Methylation ingestion: read beta value files directly from S3 via DuckDB, write Parquet."""
import csv
import pathlib

from src.gdc_client import TCGA_S3_BUCKET, fetch_manifest, get_duckdb_conn
from src.utils import get_logger

logger = get_logger(__name__)

_DATA_CATEGORY = "DNA Methylation"
_DATA_TYPE = "Methylation Beta Value"


def ingest_methylation(
    output_dir: pathlib.Path,
    project_id: str = "TCGA-BRCA",
) -> pathlib.Path:
    """
    Fetch GDC methylation manifest, read each patient's beta value file directly
    from S3 into DuckDB, and write a consolidated Parquet file.

    Methylation files are headerless TSVs with two columns: probe_id, beta_value.

    Output schema (long format):
        patient_id: VARCHAR  — 12-char TCGA barcode
        probe_id:   VARCHAR  — methylation probe identifier (e.g. "cg00000029")
        beta_value: DOUBLE   — beta value in [0.0, 1.0]

    Args:
        output_dir: Directory where methylation.parquet (and errors_methylation.csv) are written.
        project_id: GDC project ID (default "TCGA-BRCA").

    Returns:
        pathlib.Path to the written methylation.parquet.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Fetching methylation manifest for %s", project_id)
    manifest = fetch_manifest(project_id, _DATA_CATEGORY, _DATA_TYPE)
    logger.info("Manifest contains %d files", len(manifest))

    con = get_duckdb_conn()
    con.execute("""
        CREATE TABLE methylation_staging (
            patient_id VARCHAR,
            probe_id   VARCHAR,
            beta_value DOUBLE
        )
    """)

    errors: list[dict] = []

    for entry in manifest:
        s3_path = f"{TCGA_S3_BUCKET}methylation/{entry['file_id']}/{entry['file_name']}"
        try:
            con.execute("""
                INSERT INTO methylation_staging
                SELECT
                    ? AS patient_id,
                    probe_id,
                    beta_value
                FROM read_csv(
                    ?,
                    delim='\t',
                    header=false,
                    columns={'probe_id': 'VARCHAR', 'beta_value': 'DOUBLE'}
                )
            """, [entry["patient_id"], s3_path])
        except Exception as exc:  # noqa: BLE001
            logger.warning("Skipping %s (file_id=%s): %s", entry["patient_id"], entry["file_id"], exc)
            errors.append({"patient_id": entry["patient_id"], "file_id": entry["file_id"], "error": str(exc)})

    row_count = con.execute("SELECT COUNT(*) FROM methylation_staging").fetchone()[0]

    parquet_path = output_dir / "methylation.parquet"
    con.execute(f"COPY methylation_staging TO '{parquet_path}' (FORMAT PARQUET, COMPRESSION SNAPPY)")
    logger.info("Wrote %d rows to %s", row_count, parquet_path)

    if errors:
        errors_path = output_dir / "errors_methylation.csv"
        with open(errors_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["patient_id", "file_id", "error"])
            writer.writeheader()
            writer.writerows(errors)
        logger.warning("Wrote %d errors to %s", len(errors), errors_path)

    return parquet_path
