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


if __name__ == "__main__":
    unittest.main()
