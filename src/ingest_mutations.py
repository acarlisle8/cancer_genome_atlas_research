"""Mutation ingestion: read MAF files directly from S3 via DuckDB, write Parquet.

Output is a per-(patient, gene) binary indicator suitable for a Bernoulli
likelihood view in MOFA+. Silent / non-coding variants are excluded; only
mutations that change protein function are retained.
"""
import csv
import pathlib

from src.gdc_client import TCGA_S3_BUCKET, fetch_manifest, get_duckdb_conn
from src.utils import get_logger

logger = get_logger(__name__)

_DATA_CATEGORY = "Simple Nucleotide Variation"
_DATA_TYPE = "Masked Somatic Mutation"

# Variant_Classification values that change protein sequence or splicing.
# Anything not in this set (Silent, Intron, 5'UTR, 3'UTR, RNA, IGR, Flank) is excluded.
_NON_SILENT_CLASSES = (
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Splice_Site",
    "Splice_Region",
    "Translation_Start_Site",
    "Nonstop_Mutation",
)


def ingest_mutations(
    output_dir: pathlib.Path,
    project_id: str = "TCGA-BRCA",
) -> pathlib.Path:
    """
    Fetch GDC Masked Somatic Mutation manifest, read each patient's MAF file
    directly from S3 into DuckDB, filter to non-silent protein-altering variants,
    deduplicate to (patient, gene), and write a consolidated Parquet file.

    Output schema (long format, Bernoulli-ready):
        patient_id:   VARCHAR — 12-char TCGA barcode
        hugo_symbol:  VARCHAR — gene symbol (e.g. "TP53")
        mutated:      DOUBLE  — always 1.0; absent (patient, gene) pairs implicitly 0

    Args:
        output_dir: Directory where mutations.parquet (and errors_mutations.csv) are written.
        project_id: GDC project ID (default "TCGA-BRCA").

    Returns:
        pathlib.Path to the written mutations.parquet.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Fetching mutations manifest for %s", project_id)
    manifest = fetch_manifest(project_id, _DATA_CATEGORY, _DATA_TYPE)
    logger.info("Manifest contains %d files", len(manifest))

    con = get_duckdb_conn()
    con.execute("""
        CREATE TABLE mutations_staging (
            patient_id  VARCHAR,
            hugo_symbol VARCHAR,
            mutated     DOUBLE
        )
    """)

    classes_sql = ",".join(f"'{c}'" for c in _NON_SILENT_CLASSES)
    errors: list[dict] = []

    for entry in manifest:
        s3_path = f"{TCGA_S3_BUCKET}mutations/{entry['file_id']}/{entry['file_name']}"
        try:
            con.execute(f"""
                INSERT INTO mutations_staging
                SELECT DISTINCT
                    ? AS patient_id,
                    Hugo_Symbol AS hugo_symbol,
                    1.0::DOUBLE AS mutated
                FROM read_csv(
                    ?,
                    delim='\t',
                    header=true,
                    comment='#',
                    compression='gzip',
                    types={{'Hugo_Symbol':'VARCHAR','Variant_Classification':'VARCHAR'}}
                )
                WHERE Variant_Classification IN ({classes_sql})
                  AND Hugo_Symbol IS NOT NULL
                  AND Hugo_Symbol <> ''
                  AND Hugo_Symbol <> 'Unknown'
            """, [entry["patient_id"], s3_path])
        except Exception as exc:  # noqa: BLE001
            logger.warning("Skipping %s (file_id=%s): %s", entry["patient_id"], entry["file_id"], exc)
            errors.append({"patient_id": entry["patient_id"], "file_id": entry["file_id"], "error": str(exc)})

    row_count = con.execute("SELECT COUNT(*) FROM mutations_staging").fetchone()[0]
    if row_count == 0:
        raise RuntimeError(f"No rows ingested for {project_id}. Check errors_mutations.csv.")

    parquet_path = output_dir / "mutations.parquet"
    con.execute(f"COPY mutations_staging TO '{parquet_path}' (FORMAT PARQUET, COMPRESSION SNAPPY)")
    logger.info("Wrote %d rows to %s", row_count, parquet_path)

    if errors:
        errors_path = output_dir / "errors_mutations.csv"
        with open(errors_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["patient_id", "file_id", "error"])
            writer.writeheader()
            writer.writerows(errors)
        logger.warning("Wrote %d errors to %s", len(errors), errors_path)

    return parquet_path
