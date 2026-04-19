"""Tests for src/ingest_methylation.py — methylation beta value aggregator."""
import pathlib
import tempfile
import unittest
from unittest.mock import patch

import polars as pl


# Headerless 2-column TSV fixture (probe_id, beta_value)
METHYLATION_TSV_CONTENT = (
    "cg00000029\t0.123\n"
    "cg00000108\t0.456\n"
    "cg00000165\t0.789\n"
    "cg00000289\t0.012\n"
)


class TestParseMethylationBetas(unittest.TestCase):
    """Tests for parse_methylation_betas()."""

    def test_output_columns_are_exactly_three(self):
        """parse_methylation_betas output has exactly: patient_id, probe_id, beta_value."""
        from src.ingest_methylation import parse_methylation_betas

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".methylation_array.sesame.level3betas.txt", delete=False
        ) as f:
            f.write(METHYLATION_TSV_CONTENT)
            tmp_path = pathlib.Path(f.name)

        df = parse_methylation_betas(tmp_path, "TCGA-XX-0001")
        self.assertEqual(df.columns, ["patient_id", "probe_id", "beta_value"])
        tmp_path.unlink()

    def test_row_count_matches_data_rows(self):
        """parse_methylation_betas returns 4 rows for a file with 4 data rows."""
        from src.ingest_methylation import parse_methylation_betas

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        ) as f:
            f.write(METHYLATION_TSV_CONTENT)
            tmp_path = pathlib.Path(f.name)

        df = parse_methylation_betas(tmp_path, "TCGA-XX-0001")
        self.assertEqual(len(df), 4)
        tmp_path.unlink()

    def test_beta_value_dtype_is_float64(self):
        """parse_methylation_betas casts beta_value to Float64."""
        from src.ingest_methylation import parse_methylation_betas

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        ) as f:
            f.write(METHYLATION_TSV_CONTENT)
            tmp_path = pathlib.Path(f.name)

        df = parse_methylation_betas(tmp_path, "TCGA-XX-0001")
        self.assertEqual(df["beta_value"].dtype, pl.Float64)
        tmp_path.unlink()

    def test_patient_id_column_value(self):
        """parse_methylation_betas sets patient_id to the provided value."""
        from src.ingest_methylation import parse_methylation_betas

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        ) as f:
            f.write(METHYLATION_TSV_CONTENT)
            tmp_path = pathlib.Path(f.name)

        df = parse_methylation_betas(tmp_path, "TCGA-AB-9999")
        self.assertTrue((df["patient_id"] == "TCGA-AB-9999").all())
        tmp_path.unlink()

    def test_probe_id_dtype_is_utf8(self):
        """parse_methylation_betas produces Utf8 for probe_id and patient_id."""
        from src.ingest_methylation import parse_methylation_betas

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        ) as f:
            f.write(METHYLATION_TSV_CONTENT)
            tmp_path = pathlib.Path(f.name)

        df = parse_methylation_betas(tmp_path, "TCGA-XX-0001")
        self.assertEqual(df["patient_id"].dtype, pl.Utf8)
        self.assertEqual(df["probe_id"].dtype, pl.Utf8)
        tmp_path.unlink()

    def test_beta_values_correct(self):
        """parse_methylation_betas reads beta values correctly."""
        from src.ingest_methylation import parse_methylation_betas

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        ) as f:
            f.write(METHYLATION_TSV_CONTENT)
            tmp_path = pathlib.Path(f.name)

        df = parse_methylation_betas(tmp_path, "TCGA-XX-0001")
        self.assertAlmostEqual(df["beta_value"][0], 0.123)
        self.assertAlmostEqual(df["beta_value"][3], 0.012)
        tmp_path.unlink()


class TestIngestMethylation(unittest.TestCase):
    """Tests for ingest_methylation()."""

    def _make_manifest(self):
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

    def test_ingest_methylation_writes_parquet(self):
        """ingest_methylation writes methylation.parquet to output_dir."""
        from src.ingest_methylation import ingest_methylation

        manifest = self._make_manifest()

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = pathlib.Path(tmpdir) / "output"
            raw_dir = pathlib.Path(tmpdir) / "raw"
            output_dir.mkdir()
            raw_dir.mkdir()

            meth_raw_dir = raw_dir / "methylation"
            meth_raw_dir.mkdir()
            for entry in manifest:
                (meth_raw_dir / entry["file_name"]).write_text(METHYLATION_TSV_CONTENT)

            with patch("src.ingest_methylation.fetch_manifest", return_value=manifest):
                with patch(
                    "src.ingest_methylation.download_file",
                    side_effect=lambda fid, fname, dest: dest / fname,
                ):
                    result = ingest_methylation(output_dir, raw_dir, "TCGA-BRCA")

            self.assertTrue((output_dir / "methylation.parquet").exists())
            self.assertEqual(result, output_dir / "methylation.parquet")

    def test_ingest_methylation_parquet_has_correct_schema(self):
        """ingest_methylation output parquet has exactly 3 columns: patient_id, probe_id, beta_value."""
        from src.ingest_methylation import ingest_methylation

        manifest = self._make_manifest()

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = pathlib.Path(tmpdir) / "output"
            raw_dir = pathlib.Path(tmpdir) / "raw"
            output_dir.mkdir()
            raw_dir.mkdir()

            meth_raw_dir = raw_dir / "methylation"
            meth_raw_dir.mkdir()
            for entry in manifest:
                (meth_raw_dir / entry["file_name"]).write_text(METHYLATION_TSV_CONTENT)

            with patch("src.ingest_methylation.fetch_manifest", return_value=manifest):
                with patch(
                    "src.ingest_methylation.download_file",
                    side_effect=lambda fid, fname, dest: dest / fname,
                ):
                    ingest_methylation(output_dir, raw_dir, "TCGA-BRCA")

            df = pl.read_parquet(output_dir / "methylation.parquet")
            self.assertEqual(df.columns, ["patient_id", "probe_id", "beta_value"])
            # 2 patients * 4 rows each
            self.assertEqual(len(df), 8)

    def test_ingest_methylation_skip_and_log_on_parse_error(self):
        """ingest_methylation writes errors_methylation.csv and continues on parse exception."""
        from src.ingest_methylation import ingest_methylation

        manifest = [
            {
                "file_id": "file-id-001",
                "file_name": "good.methylation_array.sesame.level3betas.txt",
                "patient_id": "TCGA-AA-0001",
            },
            {
                "file_id": "file-id-bad",
                "file_name": "bad.methylation_array.sesame.level3betas.txt",
                "patient_id": "TCGA-AA-XXXX",
            },
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = pathlib.Path(tmpdir) / "output"
            raw_dir = pathlib.Path(tmpdir) / "raw"
            output_dir.mkdir()
            raw_dir.mkdir()

            meth_raw_dir = raw_dir / "methylation"
            meth_raw_dir.mkdir()

            # Good file
            (meth_raw_dir / "good.methylation_array.sesame.level3betas.txt").write_text(
                METHYLATION_TSV_CONTENT
            )
            # Bad file — wrong number of columns (3 instead of 2)
            (meth_raw_dir / "bad.methylation_array.sesame.level3betas.txt").write_text(
                "cg001\t0.5\textra_col\ncg002\t0.3\textra_col\n"
            )

            with patch("src.ingest_methylation.fetch_manifest", return_value=manifest):
                with patch(
                    "src.ingest_methylation.download_file",
                    side_effect=lambda fid, fname, dest: dest / fname,
                ):
                    ingest_methylation(output_dir, raw_dir, "TCGA-BRCA")

            # errors_methylation.csv must be written
            self.assertTrue((output_dir / "errors_methylation.csv").exists())
            # methylation.parquet must still exist (pipeline did not abort)
            self.assertTrue((output_dir / "methylation.parquet").exists())


if __name__ == "__main__":
    unittest.main()
