"""Tests for src/ingest_methylation.py — DuckDB-direct methylation aggregator."""
import csv
import pathlib
import tempfile
import unittest
from unittest.mock import patch

import duckdb
import polars as pl


# Headerless 2-column TSV fixture (probe_id, beta_value)
METHYLATION_TSV_CONTENT = (
    "cg00000029\t0.123\n"
    "cg00000108\t0.456\n"
    "cg00000165\t0.789\n"
    "cg00000289\t0.012\n"
)


def _write_fake_s3_file(fake_bucket: pathlib.Path, file_id: str, file_name: str, content: str) -> None:
    """Create <fake_bucket>/methylation/<file_id>/<file_name> with given content."""
    file_dir = fake_bucket / "methylation" / file_id
    file_dir.mkdir(parents=True, exist_ok=True)
    (file_dir / file_name).write_text(content)


class TestIngestMethylation(unittest.TestCase):
    """End-to-end tests for ingest_methylation() against a local fake S3 bucket."""

    def _manifest(self):
        return [
            {
                "file_id": "file-id-001",
                "file_name": "p1.methylation_array.sesame.level3betas.txt",
                "patient_id": "TCGA-AA-0001",
            },
            {
                "file_id": "file-id-002",
                "file_name": "p2.methylation_array.sesame.level3betas.txt",
                "patient_id": "TCGA-AA-0002",
            },
        ]

    def test_writes_parquet(self):
        from src.ingest_methylation import ingest_methylation

        manifest = self._manifest()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            for entry in manifest:
                _write_fake_s3_file(fake_bucket, entry["file_id"], entry["file_name"], METHYLATION_TSV_CONTENT)

            with patch("src.ingest_methylation.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_methylation.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_methylation.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                result = ingest_methylation(output_dir, "TCGA-BRCA")

            self.assertTrue((output_dir / "methylation.parquet").exists())
            self.assertEqual(result, output_dir / "methylation.parquet")

    def test_parquet_schema_and_row_count(self):
        from src.ingest_methylation import ingest_methylation

        manifest = self._manifest()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            for entry in manifest:
                _write_fake_s3_file(fake_bucket, entry["file_id"], entry["file_name"], METHYLATION_TSV_CONTENT)

            with patch("src.ingest_methylation.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_methylation.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_methylation.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                ingest_methylation(output_dir, "TCGA-BRCA")

            df = pl.read_parquet(output_dir / "methylation.parquet")
            self.assertEqual(df.columns, ["patient_id", "probe_id", "beta_value"])
            self.assertEqual(df["patient_id"].dtype, pl.Utf8)
            self.assertEqual(df["probe_id"].dtype, pl.Utf8)
            self.assertEqual(df["beta_value"].dtype, pl.Float64)
            # 2 patients * 4 rows each
            self.assertEqual(len(df), 8)

    def test_patient_ids_attached_correctly(self):
        from src.ingest_methylation import ingest_methylation

        manifest = self._manifest()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            for entry in manifest:
                _write_fake_s3_file(fake_bucket, entry["file_id"], entry["file_name"], METHYLATION_TSV_CONTENT)

            with patch("src.ingest_methylation.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_methylation.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_methylation.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                ingest_methylation(output_dir, "TCGA-BRCA")

            df = pl.read_parquet(output_dir / "methylation.parquet")
            patient_ids = set(df["patient_id"].to_list())
            self.assertEqual(patient_ids, {"TCGA-AA-0001", "TCGA-AA-0002"})

    def test_skip_and_log_on_parse_error(self):
        """When one file has non-numeric beta values, write errors_methylation.csv and keep going."""
        from src.ingest_methylation import ingest_methylation

        manifest = [
            {"file_id": "file-good", "file_name": "good.txt", "patient_id": "TCGA-AA-0001"},
            {"file_id": "file-bad",  "file_name": "bad.txt",  "patient_id": "TCGA-AA-XXXX"},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            _write_fake_s3_file(fake_bucket, "file-good", "good.txt", METHYLATION_TSV_CONTENT)
            _write_fake_s3_file(fake_bucket, "file-bad", "bad.txt", "cg001\tNOT_A_NUMBER\ncg002\tALSO_BAD\n")

            with patch("src.ingest_methylation.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_methylation.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_methylation.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                ingest_methylation(output_dir, "TCGA-BRCA")

            # errors_methylation.csv must be written
            errors_csv = output_dir / "errors_methylation.csv"
            self.assertTrue(errors_csv.exists())
            with open(errors_csv, newline="") as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["patient_id"], "TCGA-AA-XXXX")
            self.assertEqual(rows[0]["file_id"], "file-bad")
            self.assertIn("error", rows[0])

            # methylation.parquet still written, with only the good patient's rows
            df = pl.read_parquet(output_dir / "methylation.parquet")
            self.assertEqual(len(df), 4)
            self.assertEqual(set(df["patient_id"].to_list()), {"TCGA-AA-0001"})

    def test_no_error_csv_when_all_succeed(self):
        from src.ingest_methylation import ingest_methylation

        manifest = self._manifest()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            for entry in manifest:
                _write_fake_s3_file(fake_bucket, entry["file_id"], entry["file_name"], METHYLATION_TSV_CONTENT)

            with patch("src.ingest_methylation.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_methylation.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_methylation.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                ingest_methylation(output_dir, "TCGA-BRCA")

            self.assertFalse((output_dir / "errors_methylation.csv").exists())


if __name__ == "__main__":
    unittest.main()
