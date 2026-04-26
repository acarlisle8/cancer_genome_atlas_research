"""GDC API manifest fetcher and DuckDB S3 connection factory for TCGA data."""
import os

import duckdb
import requests

GDC_FILES_URL = "https://api.gdc.cancer.gov/files"
TCGA_S3_BUCKET = "s3://g23861422-datsbd-s2026/tcga/"


def fetch_manifest(
    project_id: str,
    data_category: str,
    data_type: str,
    size: int = 10000,
) -> list[dict]:
    """
    Query GDC API and return manifest entries.

    Returns list of {"file_id": str, "file_name": str, "patient_id": str}.
    patient_id is the 12-char TCGA barcode derived from cases[0].submitter_id
    by taking the first 3 hyphen-delimited segments: "-".join(barcode.split("-")[:3])

    Args:
        project_id: TCGA project ID, e.g. "TCGA-BRCA"
        data_category: GDC data category, e.g. "Transcriptome Profiling"
        data_type: GDC data type, e.g. "Gene Expression Quantification"
        size: Maximum number of results to retrieve (default 10000)

    Returns:
        List of dicts with keys: file_id, file_name, patient_id

    Raises:
        requests.HTTPError: If GDC API returns a non-2xx response
    """
    payload = {
        "filters": {
            "op": "and",
            "content": [
                {
                    "op": "=",
                    "content": {
                        "field": "cases.project.project_id",
                        "value": project_id,
                    },
                },
                {
                    "op": "=",
                    "content": {
                        "field": "data_category",
                        "value": data_category,
                    },
                },
                {
                    "op": "=",
                    "content": {
                        "field": "data_type",
                        "value": data_type,
                    },
                },
            ],
        },
        "fields": "file_id,file_name,cases.submitter_id",
        "format": "json",
        "size": size,
    }

    response = requests.post(GDC_FILES_URL, json=payload)
    response.raise_for_status()

    hits = response.json()["data"]["hits"]

    result = []
    for hit in hits:
        barcode = hit["cases"][0]["submitter_id"]
        patient_id = "-".join(barcode.split("-")[:3])
        result.append(
            {
                "file_id": hit["file_id"],
                "file_name": hit["file_name"],
                "patient_id": patient_id,
            }
        )

    return result


def get_duckdb_conn(db_path: str = ":memory:") -> duckdb.DuckDBPyConnection:
    """
    Create a DuckDB connection configured for authenticated S3 access.

    Reads AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY, AWS_SESSION_TOKEN,
    and AWS_DEFAULT_REGION from environment variables.

    Args:
        db_path: DuckDB database path (default ":memory:")

    Returns:
        Configured DuckDB connection ready to query S3 via read_csv / read_parquet
    """
    con = duckdb.connect(db_path)
    con.execute("INSTALL httpfs; LOAD httpfs;")
    con.execute(f"""
        SET s3_access_key_id     = '{os.environ["AWS_ACCESS_KEY_ID"]}';
        SET s3_secret_access_key = '{os.environ["AWS_SECRET_ACCESS_KEY"]}';
        SET s3_session_token     = '{os.environ.get("AWS_SESSION_TOKEN", "")}';
        SET s3_region            = '{os.environ.get("AWS_DEFAULT_REGION", "us-east-1")}';
    """)
    return con

