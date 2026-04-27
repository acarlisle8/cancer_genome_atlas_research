"""CNV ingestion: read segment files directly from S3 via DuckDB, write Parquet."""
import pathlib

import polars as pl

from src.gdc_client import TCGA_S3_BUCKET, fetch_manifest, get_duckdb_conn
from src.utils import get_logger

logger = get_logger(__name__)


def ingest_cnv(
    output_dir: pathlib.Path,
    project_id: str = "TCGA-BRCA",
) -> pathlib.Path:
    """
    Fetch GDC CNV manifest, read each patient's segment file directly from S3
    into DuckDB, and write a consolidated Parquet file.

    Output schema (long format):
        patient_id:  VARCHAR  — 12-char TCGA barcode
        chromosome:  VARCHAR
        start:       BIGINT
        end:         BIGINT
        copy_number: DOUBLE  — log2 copy ratio (Segment_Mean)

    Args:
        output_dir: Directory to write cnv.parquet (and errors_cnv.csv if needed).
        project_id: TCGA project identifier (default: TCGA-BRCA).

    Returns:
        pathlib.Path to the written cnv.parquet file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Snagging CNV manifest for %s", project_id)
    manifest = fetch_manifest(project_id, "Copy Number Variation", "Masked Copy Number Segment")
    logger.info("Manifest contains %d CNV files", len(manifest))

    con = get_duckdb_conn()
    con.execute("""
        CREATE TABLE cnv_staging (
            patient_id  VARCHAR,
            chromosome  VARCHAR,
            start       BIGINT,
            "end"       BIGINT,
            copy_number DOUBLE
        )
    """)

    errors: list[dict] = []

    for entry in manifest:
        s3_path = f"{TCGA_S3_BUCKET}cnv/{entry['file_id']}/{entry['file_name']}"
        try:
            con.execute("""
                INSERT INTO cnv_staging
                SELECT
                    ? AS patient_id,
                    Chromosome   AS chromosome,
                    Start::BIGINT,
                    "End"::BIGINT,
                    Segment_Mean::DOUBLE AS copy_number
                FROM read_csv(?, delim='\t', header=true)
            """, [entry["patient_id"], s3_path])
        except Exception as exc:  # noqa: BLE001
            logger.warning("Skipping %s (file_id=%s): %s", entry["patient_id"], entry["file_id"], exc)
            errors.append({"patient_id": entry["patient_id"], "file_id": entry["file_id"], "error": str(exc)})

    row_count = con.execute("SELECT COUNT(*) FROM cnv_staging").fetchone()[0]

    if row_count == 0:
        logger.warning("No CNV rows ingested; writing empty Parquet")
        combined = pl.DataFrame({
            "patient_id": pl.Series([], dtype=pl.Utf8),
            "chromosome": pl.Series([], dtype=pl.Utf8),
            "start": pl.Series([], dtype=pl.Int64),
            "end": pl.Series([], dtype=pl.Int64),
            "copy_number": pl.Series([], dtype=pl.Float64),
        })
        parquet_path = output_dir / "cnv.parquet"
        combined.write_parquet(parquet_path, compression="snappy")
    else:
        parquet_path = output_dir / "cnv.parquet"
        con.execute(f"COPY cnv_staging TO '{parquet_path}' (FORMAT PARQUET, COMPRESSION SNAPPY)")

    logger.info("Wrote %d rows to %s", row_count, parquet_path)

    if errors:
        errors_path = output_dir / "errors_cnv.csv"
        pl.DataFrame(errors).write_csv(errors_path)
        logger.warning("Wrote %d errors to %s", len(errors), errors_path)

    return parquet_path
