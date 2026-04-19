"""Methylation beta value aggregator: fetch manifest, download, parse to long-format, write Parquet."""
import pathlib

import polars as pl

from src.gdc_client import download_file, fetch_manifest
from src.utils import get_logger

logger = get_logger(__name__)

METHYLATION_OUTPUT_COLS = ["patient_id", "probe_id", "beta_value"]


def parse_methylation_betas(path: pathlib.Path, patient_id: str) -> pl.DataFrame:
    """
    Read a level3 beta value file (headerless, 2 columns: probe_id and beta_value).

    File format: tab-separated, NO header row. First column = probe ID (e.g., "cg00000029"),
    second column = beta value (float in [0.0, 1.0]).

    Args:
        path: Path to the .methylation_array.sesame.level3betas.txt file
        patient_id: 12-char TCGA patient barcode to attach as column

    Returns:
        DataFrame with columns: [patient_id, probe_id, beta_value]
        Dtypes: [Utf8, Utf8, Float64]
    """
    df = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=["probe_id", "beta_value"],
    )

    # Cast beta_value to Float64 (raises if non-numeric, triggering D-06 error handler)
    df = df.with_columns(pl.col("beta_value").cast(pl.Float64))

    # Add patient_id and select output columns
    df = df.with_columns(pl.lit(patient_id).cast(pl.Utf8).alias("patient_id"))
    df = df.select(METHYLATION_OUTPUT_COLS)

    return df


def ingest_methylation(
    output_dir: pathlib.Path,
    raw_dir: pathlib.Path,
    project_id: str = "TCGA-BRCA",
) -> pathlib.Path:
    """
    Fetch methylation manifest, download per-patient beta value files, aggregate to Parquet.

    Steps:
    1. fetch_manifest for DNA Methylation / Methylation Beta Value
    2. Download each file into raw_dir / "methylation"
    3. parse_methylation_betas per file; on exception: log + record in errors list (D-06)
    4. Concatenate all parsed DataFrames
    5. Write output_dir / "methylation.parquet" with snappy compression
    6. If any errors: write output_dir / "errors_methylation.csv"
    7. Return Parquet path

    Args:
        output_dir: Directory to write methylation.parquet (and errors_methylation.csv if needed)
        raw_dir: Base directory for raw downloads; methylation files go in raw_dir / "methylation"
        project_id: TCGA project identifier (default: TCGA-BRCA)

    Returns:
        pathlib.Path to the written methylation.parquet file
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    methylation_raw_dir = raw_dir / "methylation"

    logger.info("Fetching methylation manifest for %s", project_id)
    manifest = fetch_manifest(
        project_id,
        "DNA Methylation",
        "Methylation Beta Value",
    )
    logger.info("Manifest contains %d methylation files", len(manifest))

    frames: list[pl.DataFrame] = []
    errors: list[dict] = []

    for entry in manifest:
        file_id = entry["file_id"]
        file_name = entry["file_name"]
        patient_id = entry["patient_id"]

        try:
            local_path = download_file(file_id, file_name, methylation_raw_dir)
            df = parse_methylation_betas(local_path, patient_id)
            frames.append(df)
        except Exception as exc:  # noqa: BLE001
            logger.warning(
                "Skipping patient %s (file %s): %s", patient_id, file_id, exc
            )
            errors.append(
                {
                    "patient_id": patient_id,
                    "file_id": file_id,
                    "error": str(exc),
                }
            )

    if not frames:
        logger.warning("No methylation frames parsed; writing empty Parquet")
        combined = pl.DataFrame(
            {
                "patient_id": pl.Series([], dtype=pl.Utf8),
                "probe_id": pl.Series([], dtype=pl.Utf8),
                "beta_value": pl.Series([], dtype=pl.Float64),
            }
        )
    else:
        combined = pl.concat(frames)

    parquet_path = output_dir / "methylation.parquet"
    combined.write_parquet(parquet_path, compression="snappy")
    logger.info("Wrote %d rows to %s", len(combined), parquet_path)

    if errors:
        errors_path = output_dir / "errors_methylation.csv"
        pl.DataFrame(errors).write_csv(errors_path)
        logger.warning("Wrote %d errors to %s", len(errors), errors_path)

    return parquet_path
