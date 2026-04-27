"""Tests for src/ingest_cnv.py — DuckDB-direct CNV aggregator."""
import csv
import pathlib
import tempfile
import unittest
from unittest.mock import patch

import duckdb
import polars as pl


CNV_TSV_CONTENT = (
    "GDC_Aliquot\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n"
    "ALIQUOT-01\tchr1\t100\t200\t50\t0.123\n"
    "ALIQUOT-01\tchr2\t300\t400\t75\t-0.456\n"
    "ALIQUOT-01\tchrX\t500\t600\t30\t1.234\n"
)


def _write_fake_s3_file(fake_bucket: pathlib.Path, file_id: str, file_name: str, content: str) -> None:
    """Create <fake_bucket>/cnv/<file_id>/<file_name> with given content."""
    file_dir = fake_bucket / "cnv" / file_id
    file_dir.mkdir(parents=True, exist_ok=True)
    (file_dir / file_name).write_text(content)


class TestIngestCnv(unittest.TestCase):
    """End-to-end tests for ingest_cnv() against a local fake S3 bucket."""

    def _manifest(self):
        return [
            {"file_id": "file-id-001", "file_name": "p1.seg.v2.txt", "patient_id": "TCGA-AA-0001"},
            {"file_id": "file-id-002", "file_name": "p2.seg.v2.txt", "patient_id": "TCGA-AA-0002"},
        ]

    def test_writes_parquet(self):
        from src.ingest_cnv import ingest_cnv

        manifest = self._manifest()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            for entry in manifest:
                _write_fake_s3_file(fake_bucket, entry["file_id"], entry["file_name"], CNV_TSV_CONTENT)

            with patch("src.ingest_cnv.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_cnv.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_cnv.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                result = ingest_cnv(output_dir, "TCGA-BRCA")

            self.assertTrue((output_dir / "cnv.parquet").exists())
            self.assertEqual(result, output_dir / "cnv.parquet")

    def test_parquet_schema_and_row_count(self):
        from src.ingest_cnv import ingest_cnv

        manifest = self._manifest()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            for entry in manifest:
                _write_fake_s3_file(fake_bucket, entry["file_id"], entry["file_name"], CNV_TSV_CONTENT)

            with patch("src.ingest_cnv.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_cnv.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_cnv.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                ingest_cnv(output_dir, "TCGA-BRCA")

            df = pl.read_parquet(output_dir / "cnv.parquet")
            self.assertEqual(df.columns, ["patient_id", "chromosome", "start", "end", "copy_number"])
            self.assertEqual(df["patient_id"].dtype, pl.Utf8)
            self.assertEqual(df["chromosome"].dtype, pl.Utf8)
            self.assertEqual(df["start"].dtype, pl.Int64)
            self.assertEqual(df["end"].dtype, pl.Int64)
            self.assertEqual(df["copy_number"].dtype, pl.Float64)
            # 2 patients * 3 rows each
            self.assertEqual(len(df), 6)

    def test_patient_ids_attached_correctly(self):
        from src.ingest_cnv import ingest_cnv

        manifest = self._manifest()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            for entry in manifest:
                _write_fake_s3_file(fake_bucket, entry["file_id"], entry["file_name"], CNV_TSV_CONTENT)

            with patch("src.ingest_cnv.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_cnv.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_cnv.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                ingest_cnv(output_dir, "TCGA-BRCA")

            df = pl.read_parquet(output_dir / "cnv.parquet")
            self.assertEqual(set(df["patient_id"].to_list()), {"TCGA-AA-0001", "TCGA-AA-0002"})

    def test_skip_and_log_on_parse_error(self):
        """When one file is malformed, write errors_cnv.csv and keep going."""
        from src.ingest_cnv import ingest_cnv

        manifest = [
            {"file_id": "file-good", "file_name": "good.seg.txt", "patient_id": "TCGA-AA-0001"},
            {"file_id": "file-bad",  "file_name": "bad.seg.txt",  "patient_id": "TCGA-AA-XXXX"},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            _write_fake_s3_file(fake_bucket, "file-good", "good.seg.txt", CNV_TSV_CONTENT)
            _write_fake_s3_file(fake_bucket, "file-bad", "bad.seg.txt", "not\ta\tvalid\tfile\n1\t2\t3\t4\n")

            with patch("src.ingest_cnv.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_cnv.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_cnv.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                ingest_cnv(output_dir, "TCGA-BRCA")

            errors_csv = output_dir / "errors_cnv.csv"
            self.assertTrue(errors_csv.exists())
            with open(errors_csv, newline="") as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["patient_id"], "TCGA-AA-XXXX")
            self.assertEqual(rows[0]["file_id"], "file-bad")
            self.assertIn("error", rows[0])

            df = pl.read_parquet(output_dir / "cnv.parquet")
            self.assertEqual(len(df), 3)
            self.assertEqual(set(df["patient_id"].to_list()), {"TCGA-AA-0001"})

    def test_no_error_csv_when_all_succeed(self):
        from src.ingest_cnv import ingest_cnv

        manifest = self._manifest()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            for entry in manifest:
                _write_fake_s3_file(fake_bucket, entry["file_id"], entry["file_name"], CNV_TSV_CONTENT)

            with patch("src.ingest_cnv.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_cnv.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_cnv.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                ingest_cnv(output_dir, "TCGA-BRCA")

            self.assertFalse((output_dir / "errors_cnv.csv").exists())


if __name__ == "__main__":
    unittest.main()
