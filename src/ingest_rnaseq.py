"""RNA-seq ingestion and aggregation for TCGA BRCA cohort.

Fetches the GDC manifest, downloads each patient's augmented_star_gene_counts.tsv
from S3, parses FPKM values into long-format rows, and writes a single consolidated
Parquet file.

Output schema (long format, D-01):
    patient_id:       pl.Utf8   — 12-char TCGA barcode (e.g. "TCGA-BH-A18H")
    gene_id:          pl.Utf8   — Ensembl ID with version (e.g. "ENSG00000000003.15")
    fpkm_unstranded:  pl.Float64
"""

import csv
import pathlib

import polars as pl

from src.gdc_client import fetch_manifest, download_file
from src.utils import get_logger

logger = get_logger(__name__)

# GDC data category / type constants for RNA-seq
_DATA_CATEGORY = "Transcriptome Profiling"
_DATA_TYPE = "Gene Expression Quantification"


def parse_rnaseq_tsv(path: pathlib.Path, patient_id: str) -> pl.DataFrame:
    """
    Read a single augmented_star_gene_counts.tsv file and return long-format rows.

    Steps:
      1. Read TSV with polars (tab separator).
      2. Filter out rows where gene_id starts with "N_" (summary stat rows such as
         N_unmapped, N_multimapping, N_noFeature, N_ambiguous).
      3. Select ["gene_id", "fpkm_unstranded"] and cast fpkm_unstranded to Float64.
      4. Add literal patient_id column.
      5. Return DataFrame with columns in order: ["patient_id", "gene_id", "fpkm_unstranded"].

    Args:
        path:       Path to the .tsv file.
        patient_id: 12-char TCGA barcode to tag all rows with.

    Returns:
        pl.DataFrame with schema: patient_id (Utf8), gene_id (Utf8), fpkm_unstranded (Float64).

    Raises:
        Exception: Re-raises any polars / IO exception so the caller can skip-and-log (D-06).
    """
    df = pl.read_csv(path, separator="\t")

    # Filter out N_ summary rows (T-02-01 mitigation: select only known columns)
    df = df.filter(~pl.col("gene_id").str.starts_with("N_"))

    # Select and cast the two data columns, then add patient literal
    df = (
        df.select([
            pl.col("gene_id").cast(pl.Utf8),
            pl.col("fpkm_unstranded").cast(pl.Float64),
        ])
        .with_columns(pl.lit(patient_id).alias("patient_id"))
        .select(["patient_id", "gene_id", "fpkm_unstranded"])
    )

    return df


def ingest_rnaseq(
    output_dir: pathlib.Path,
    raw_dir: pathlib.Path,
    project_id: str = "TCGA-BRCA",
) -> pathlib.Path:
    """
    Full ingestion pipeline for the RNA-seq modality.

    Steps:
      1. Fetch manifest from GDC for the given project (Transcriptome Profiling /
         Gene Expression Quantification).
      2. For each manifest entry, download the file from S3 to raw_dir/rnaseq/.
      3. Parse each downloaded TSV with parse_rnaseq_tsv. On exception: log a
         warning and record the failure — do NOT abort (D-06 skip-and-log).
      4. Concatenate all successful DataFrames.
      5. Write to output_dir/rna_seq.parquet (snappy compression).
      6. If any errors occurred, write output_dir/errors_rnaseq.csv.
      7. Return the path to the written Parquet file.

    Args:
        output_dir:  Directory where rna_seq.parquet (and errors_rnaseq.csv) are written.
        raw_dir:     Root directory for raw downloads; files land in raw_dir/rnaseq/.
        project_id:  GDC project ID (default "TCGA-BRCA").

    Returns:
        pathlib.Path to the written rna_seq.parquet.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    rnaseq_raw_dir = raw_dir / "rnaseq"

    logger.info("Fetching manifest for project=%s …", project_id)
    entries = fetch_manifest(project_id, _DATA_CATEGORY, _DATA_TYPE)
    logger.info("Manifest contains %d files", len(entries))

    frames: list[pl.DataFrame] = []
    errors: list[dict] = []

    for entry in entries:
        file_id = entry["file_id"]
        file_name = entry["file_name"]
        patient_id = entry["patient_id"]

        try:
            local_path = download_file(file_id, file_name, rnaseq_raw_dir)
            df = parse_rnaseq_tsv(local_path, patient_id)
            frames.append(df)
        except Exception as exc:  # noqa: BLE001
            logger.warning(
                "Skipping patient %s (file_id=%s): %s",
                patient_id,
                file_id,
                exc,
            )
            errors.append({
                "patient_id": patient_id,
                "file_id": file_id,
                "reason": str(exc),
            })

    if not frames:
        raise RuntimeError(
            f"No files parsed successfully for project {project_id}. "
            "Check errors_rnaseq.csv for details."
        )

    combined = pl.concat(frames)
    parquet_path = output_dir / "rna_seq.parquet"
    combined.write_parquet(parquet_path, compression="snappy")
    logger.info("Wrote %d rows to %s", len(combined), parquet_path)

    if errors:
        errors_path = output_dir / "errors_rnaseq.csv"
        with open(errors_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["patient_id", "file_id", "reason"])
            writer.writeheader()
            writer.writerows(errors)
        logger.warning(
            "Wrote %d error entries to %s", len(errors), errors_path
        )

    return parquet_path
