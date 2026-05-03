"""Tests for src/ingest_cnv.py — CNV segment aggregator."""
import io
import pathlib
import tempfile
import unittest
from unittest.mock import MagicMock, patch

import polars as pl


CNV_TSV_CONTENT = (
    "GDC_Aliquot\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n"
    "ALIQUOT-01\tchr1\t100\t200\t50\t0.123\n"
    "ALIQUOT-01\tchr2\t300\t400\t75\t-0.456\n"
    "ALIQUOT-01\tchrX\t500\t600\t30\t1.234\n"
)


class TestParseCnvSeg(unittest.TestCase):
    """Tests for parse_cnv_seg()."""

    def test_output_columns_are_exactly_five(self):
        """parse_cnv_seg output has exactly: patient_id, chromosome, start, end, copy_number."""
        from src.ingest_cnv import parse_cnv_seg

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".seg.txt", delete=False
        ) as f:
            f.write(CNV_TSV_CONTENT)
            tmp_path = pathlib.Path(f.name)

        df = parse_cnv_seg(tmp_path, "TCGA-XX-0001")
        self.assertEqual(
            df.columns, ["patient_id", "chromosome", "start", "end", "copy_number"]
        )
        tmp_path.unlink()

    def test_row_count_matches_data_rows(self):
        """parse_cnv_seg returns 3 rows for a file with 3 data rows."""
        from src.ingest_cnv import parse_cnv_seg

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".seg.txt", delete=False
        ) as f:
            f.write(CNV_TSV_CONTENT)
            tmp_path = pathlib.Path(f.name)

        df = parse_cnv_seg(tmp_path, "TCGA-XX-0001")
        self.assertEqual(len(df), 3)
        tmp_path.unlink()

    def test_no_extra_columns(self):
        """parse_cnv_seg drops GDC_Aliquot and Num_Probes columns."""
        from src.ingest_cnv import parse_cnv_seg

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".seg.txt", delete=False
        ) as f:
            f.write(CNV_TSV_CONTENT)
            tmp_path = pathlib.Path(f.name)

        df = parse_cnv_seg(tmp_path, "TCGA-XX-0001")
        self.assertNotIn("GDC_Aliquot", df.columns)
        self.assertNotIn("Num_Probes", df.columns)
        tmp_path.unlink()

    def test_patient_id_column_value(self):
        """parse_cnv_seg sets patient_id to the provided value."""
        from src.ingest_cnv import parse_cnv_seg

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".seg.txt", delete=False
        ) as f:
            f.write(CNV_TSV_CONTENT)
            tmp_path = pathlib.Path(f.name)

        df = parse_cnv_seg(tmp_path, "TCGA-AB-9999")
        self.assertTrue((df["patient_id"] == "TCGA-AB-9999").all())
        tmp_path.unlink()

    def test_dtypes(self):
        """parse_cnv_seg produces correct dtypes: patient_id/chromosome=Utf8, start/end=Int64, copy_number=Float64."""
        from src.ingest_cnv import parse_cnv_seg

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".seg.txt", delete=False
        ) as f:
            f.write(CNV_TSV_CONTENT)
            tmp_path = pathlib.Path(f.name)

        df = parse_cnv_seg(tmp_path, "TCGA-XX-0001")
        self.assertEqual(df["patient_id"].dtype, pl.Utf8)
        self.assertEqual(df["chromosome"].dtype, pl.Utf8)
        self.assertEqual(df["start"].dtype, pl.Int64)
        self.assertEqual(df["end"].dtype, pl.Int64)
        self.assertEqual(df["copy_number"].dtype, pl.Float64)
        tmp_path.unlink()

    def test_segment_mean_renamed_to_copy_number(self):
        """parse_cnv_seg renames Segment_Mean -> copy_number with correct values."""
        from src.ingest_cnv import parse_cnv_seg

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".seg.txt", delete=False
        ) as f:
            f.write(CNV_TSV_CONTENT)
            tmp_path = pathlib.Path(f.name)

        df = parse_cnv_seg(tmp_path, "TCGA-XX-0001")
        self.assertAlmostEqual(df["copy_number"][0], 0.123)
        self.assertAlmostEqual(df["copy_number"][1], -0.456)
        tmp_path.unlink()


if __name__ == "__main__":
    unittest.main()
