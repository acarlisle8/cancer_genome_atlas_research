"""Tests for src/ingest_rnaseq.py — DuckDB-direct RNA-seq aggregator."""
import csv
import pathlib
import tempfile
import unittest
from unittest.mock import patch

import duckdb
import polars as pl


RNASEQ_TSV_HEADER = (
    "gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second"
    "\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded\n"
)
RNASEQ_TSV_DATA = (
    "ENSG00000000003.15\tTSPAN6\tprotein_coding\t1200\t600\t600\t5.1\t2.3\t3.1\n"
    "ENSG00000000005.6\tTNMD\tprotein_coding\t0\t0\t0\t0.0\t0.0\t0.0\n"
    "ENSG00000000419.12\tDPM1\tprotein_coding\t800\t400\t400\t3.2\t1.1\t2.0\n"
    "ENSG00000000457.14\tSCYL3\tprotein_coding\t300\t150\t150\t1.5\t0.7\t1.0\n"
    "ENSG00000000460.17\tC1orf112\tprotein_coding\t500\t250\t250\t2.0\t0.9\t1.3\n"
)
# N_ summary rows must be filtered out by the WHERE clause in production SQL
RNASEQ_TSV_N_ROWS = (
    "N_unmapped\t\t\t100\t50\t50\t0.0\t0.0\t0.0\n"
    "N_multimapping\t\t\t200\t100\t100\t0.0\t0.0\t0.0\n"
    "N_noFeature\t\t\t50\t25\t25\t0.0\t0.0\t0.0\n"
    "N_ambiguous\t\t\t30\t15\t15\t0.0\t0.0\t0.0\n"
)
RNASEQ_TSV_CONTENT = RNASEQ_TSV_HEADER + RNASEQ_TSV_DATA + RNASEQ_TSV_N_ROWS


def _write_fake_s3_file(fake_bucket: pathlib.Path, file_id: str, file_name: str, content: str) -> None:
    """Create <fake_bucket>/rnaseq/<file_id>/<file_name> with given content."""
    file_dir = fake_bucket / "rnaseq" / file_id
    file_dir.mkdir(parents=True, exist_ok=True)
    (file_dir / file_name).write_text(content)


class TestIngestRnaseq(unittest.TestCase):
    """End-to-end tests for ingest_rnaseq() against a local fake S3 bucket."""

    def _manifest(self):
        return [
            {"file_id": "file-001", "file_name": "p1.tsv", "patient_id": "TCGA-BH-A18H"},
            {"file_id": "file-002", "file_name": "p2.tsv", "patient_id": "TCGA-E2-A14P"},
        ]

    def test_writes_parquet(self):
        from src.ingest_rnaseq import ingest_rnaseq

        manifest = self._manifest()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            for entry in manifest:
                _write_fake_s3_file(fake_bucket, entry["file_id"], entry["file_name"], RNASEQ_TSV_CONTENT)

            with patch("src.ingest_rnaseq.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_rnaseq.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_rnaseq.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                result = ingest_rnaseq(output_dir, "TCGA-BRCA")

            self.assertTrue((output_dir / "rna_seq.parquet").exists())
            self.assertEqual(result, output_dir / "rna_seq.parquet")

    def test_parquet_schema_and_row_count(self):
        from src.ingest_rnaseq import ingest_rnaseq

        manifest = self._manifest()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            for entry in manifest:
                _write_fake_s3_file(fake_bucket, entry["file_id"], entry["file_name"], RNASEQ_TSV_CONTENT)

            with patch("src.ingest_rnaseq.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_rnaseq.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_rnaseq.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                ingest_rnaseq(output_dir, "TCGA-BRCA")

            df = pl.read_parquet(output_dir / "rna_seq.parquet")
            self.assertEqual(df.columns, ["patient_id", "gene_id", "fpkm_unstranded"])
            self.assertEqual(df["patient_id"].dtype, pl.Utf8)
            self.assertEqual(df["gene_id"].dtype, pl.Utf8)
            self.assertEqual(df["fpkm_unstranded"].dtype, pl.Float64)
            # 2 patients * 5 data rows each (4 N_ summary rows per file are filtered)
            self.assertEqual(len(df), 10)

    def test_n_summary_rows_filtered_out(self):
        from src.ingest_rnaseq import ingest_rnaseq

        manifest = self._manifest()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            for entry in manifest:
                _write_fake_s3_file(fake_bucket, entry["file_id"], entry["file_name"], RNASEQ_TSV_CONTENT)

            with patch("src.ingest_rnaseq.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_rnaseq.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_rnaseq.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                ingest_rnaseq(output_dir, "TCGA-BRCA")

            df = pl.read_parquet(output_dir / "rna_seq.parquet")
            self.assertFalse(any(g.startswith("N_") for g in df["gene_id"].to_list()))

    def test_patient_ids_attached_correctly(self):
        from src.ingest_rnaseq import ingest_rnaseq

        manifest = self._manifest()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            for entry in manifest:
                _write_fake_s3_file(fake_bucket, entry["file_id"], entry["file_name"], RNASEQ_TSV_CONTENT)

            with patch("src.ingest_rnaseq.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_rnaseq.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_rnaseq.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                ingest_rnaseq(output_dir, "TCGA-BRCA")

            df = pl.read_parquet(output_dir / "rna_seq.parquet")
            self.assertEqual(set(df["patient_id"].to_list()), {"TCGA-BH-A18H", "TCGA-E2-A14P"})

    def test_skip_and_log_on_parse_error(self):
        """When one file is malformed, write errors_rnaseq.csv and keep going."""
        from src.ingest_rnaseq import ingest_rnaseq

        manifest = [
            {"file_id": "file-good", "file_name": "good.tsv", "patient_id": "TCGA-BH-A18H"},
            {"file_id": "file-bad",  "file_name": "bad.tsv",  "patient_id": "TCGA-XX-XXXX"},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            _write_fake_s3_file(fake_bucket, "file-good", "good.tsv", RNASEQ_TSV_CONTENT)
            _write_fake_s3_file(fake_bucket, "file-bad", "bad.tsv", "not\ta\tvalid\theader\n1\t2\t3\t4\n")

            with patch("src.ingest_rnaseq.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_rnaseq.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_rnaseq.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                ingest_rnaseq(output_dir, "TCGA-BRCA")

            errors_csv = output_dir / "errors_rnaseq.csv"
            self.assertTrue(errors_csv.exists())
            with open(errors_csv, newline="") as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["patient_id"], "TCGA-XX-XXXX")
            self.assertEqual(rows[0]["file_id"], "file-bad")
            self.assertIn("reason", rows[0])

            df = pl.read_parquet(output_dir / "rna_seq.parquet")
            self.assertEqual(len(df), 5)
            self.assertEqual(set(df["patient_id"].to_list()), {"TCGA-BH-A18H"})

    def test_no_error_csv_when_all_succeed(self):
        from src.ingest_rnaseq import ingest_rnaseq

        manifest = self._manifest()
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / "output"
            fake_bucket = tmp_path / "fake_s3"
            for entry in manifest:
                _write_fake_s3_file(fake_bucket, entry["file_id"], entry["file_name"], RNASEQ_TSV_CONTENT)

            with patch("src.ingest_rnaseq.TCGA_S3_BUCKET", f"{fake_bucket}/"), \
                 patch("src.ingest_rnaseq.fetch_manifest", return_value=manifest), \
                 patch("src.ingest_rnaseq.get_duckdb_conn", return_value=duckdb.connect(":memory:")):
                ingest_rnaseq(output_dir, "TCGA-BRCA")

            self.assertFalse((output_dir / "errors_rnaseq.csv").exists())


if __name__ == "__main__":
    unittest.main()
