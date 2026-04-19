"""CNV segment aggregator: fetch manifest, download, parse to long-format, write Parquet."""
import pathlib

import polars as pl

from src.gdc_client import download_file, fetch_manifest
from src.utils import get_logger

logger = get_logger(__name__)

CNV_RENAME = {
    "Chromosome": "chromosome",
    "Start": "start",
    "End": "end",
    "Segment_Mean": "copy_number",
}
CNV_OUTPUT_COLS = ["patient_id", "chromosome", "start", "end", "copy_number"]
CNV_DTYPES = {
    "chromosome": pl.Utf8,
    "start": pl.Int64,
    "end": pl.Int64,
    "copy_number": pl.Float64,
}


def parse_cnv_seg(path: pathlib.Path, patient_id: str) -> pl.DataFrame:
    """
    Read a masked copy number segment file and return a long-format DataFrame.

    Expects a tab-separated file with header including columns:
    Chromosome, Start, End, Segment_Mean (plus optional GDC_Aliquot, Num_Probes).

    Args:
        path: Path to the .seg.txt or .seg.v2.txt file
        patient_id: 12-char TCGA patient barcode to attach as column

    Returns:
        DataFrame with columns: [patient_id, chromosome, start, end, copy_number]
        Dtypes: [Utf8, Utf8, Int64, Int64, Float64]
    """
    df = pl.read_csv(path, separator="\t")

    # Rename source columns to output names
    df = df.rename({k: v for k, v in CNV_RENAME.items() if k in df.columns})

    # Cast numeric columns
    for col_name, dtype in CNV_DTYPES.items():
        if col_name in df.columns:
            df = df.with_columns(pl.col(col_name).cast(dtype))

    # Add patient_id and select only output columns
    df = df.with_columns(pl.lit(patient_id).cast(pl.Utf8).alias("patient_id"))
    df = df.select(CNV_OUTPUT_COLS)

    return df


def ingest_cnv(
    output_dir: pathlib.Path,
    raw_dir: pathlib.Path,
    project_id: str = "TCGA-BRCA",
) -> pathlib.Path:
    """
    Fetch CNV manifest, download per-patient segment files, aggregate to Parquet.

    Steps:
    1. fetch_manifest for Copy Number Variation / Masked Copy Number Segment
    2. Download each file into raw_dir / "cnv"
    3. parse_cnv_seg per file; on exception: log + record in errors list (D-06)
    4. Concatenate all parsed DataFrames
    5. Write output_dir / "cnv.parquet" with snappy compression
    6. If any errors: write output_dir / "errors_cnv.csv"
    7. Return Parquet path

    Args:
        output_dir: Directory to write cnv.parquet (and errors_cnv.csv if needed)
        raw_dir: Base directory for raw downloads; CNV files go in raw_dir / "cnv"
        project_id: TCGA project identifier (default: TCGA-BRCA)

    Returns:
        pathlib.Path to the written cnv.parquet file
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    cnv_raw_dir = raw_dir / "cnv"

    logger.info("Fetching CNV manifest for %s", project_id)
    manifest = fetch_manifest(
        project_id,
        "Copy Number Variation",
        "Masked Copy Number Segment",
    )
    logger.info("Manifest contains %d CNV files", len(manifest))

    frames: list[pl.DataFrame] = []
    errors: list[dict] = []

    for entry in manifest:
        file_id = entry["file_id"]
        file_name = entry["file_name"]
        patient_id = entry["patient_id"]

        try:
            local_path = download_file(file_id, file_name, cnv_raw_dir)
            df = parse_cnv_seg(local_path, patient_id)
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
        logger.warning("No CNV frames parsed; writing empty Parquet")
        combined = pl.DataFrame(
            {
                "patient_id": pl.Series([], dtype=pl.Utf8),
                "chromosome": pl.Series([], dtype=pl.Utf8),
                "start": pl.Series([], dtype=pl.Int64),
                "end": pl.Series([], dtype=pl.Int64),
                "copy_number": pl.Series([], dtype=pl.Float64),
            }
        )
    else:
        combined = pl.concat(frames)

    parquet_path = output_dir / "cnv.parquet"
    combined.write_parquet(parquet_path, compression="snappy")
    logger.info("Wrote %d rows to %s", len(combined), parquet_path)

    if errors:
        errors_path = output_dir / "errors_cnv.csv"
        pl.DataFrame(errors).write_csv(errors_path)
        logger.warning("Wrote %d errors to %s", len(errors), errors_path)

    return parquet_path
