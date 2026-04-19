"""GDC API manifest fetcher and S3 downloader for TCGA data."""
import pathlib

import requests
import s3fs

GDC_FILES_URL = "https://api.gdc.cancer.gov/files"
TCGA_S3_BUCKET = "tcga-2-open"


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


def download_file(
    file_id: str,
    file_name: str,
    dest_dir: pathlib.Path,
) -> pathlib.Path:
    """
    Download one file from s3://tcga-2-open/{file_id}/{file_name} to dest_dir/file_name.

    If dest_dir/file_name already exists, skip download and return path (D-07 resume).
    Uses s3fs.S3FileSystem(anon=True) for anonymous S3 access.

    Args:
        file_id: GDC file UUID
        file_name: Name of the file
        dest_dir: Local directory to download into (will be created if missing)

    Returns:
        pathlib.Path to the downloaded (or already-existing) local file
    """
    dest_path = dest_dir / file_name

    if dest_path.exists():
        return dest_path

    dest_dir.mkdir(parents=True, exist_ok=True)

    s3_path = f"s3://{TCGA_S3_BUCKET}/{file_id}/{file_name}"
    fs = s3fs.S3FileSystem(anon=True)

    with fs.open(s3_path, "rb") as remote_file:
        content = remote_file.read()

    with open(dest_path, "wb") as local_file:
        local_file.write(content)

    return dest_path
