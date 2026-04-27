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


class TestIngestCnv(unittest.TestCase):
    """Tests for ingest_cnv()."""

    def _make_manifest(self):
        return [
            {
                "file_id": "file-id-001",
                "file_name": "patient1.seg.v2.txt",
                "patient_id": "TCGA-AA-0001",
            },
            {
                "file_id": "file-id-002",
                "file_name": "patient2.seg.v2.txt",
                "patient_id": "TCGA-AA-0002",
            },
        ]

    def test_ingest_cnv_writes_parquet(self):
        """ingest_cnv writes cnv.parquet to output_dir."""
        from src.ingest_cnv import ingest_cnv

        manifest = self._make_manifest()

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = pathlib.Path(tmpdir) / "output"
            raw_dir = pathlib.Path(tmpdir) / "raw"
            output_dir.mkdir()
            raw_dir.mkdir()

            # Create fake downloaded files
            cnv_raw_dir = raw_dir / "cnv"
            cnv_raw_dir.mkdir()
            for entry in manifest:
                fpath = cnv_raw_dir / entry["file_name"]
                fpath.write_text(CNV_TSV_CONTENT)

            with patch("src.ingest_cnv.fetch_manifest", return_value=manifest):
                with patch(
                    "src.ingest_cnv.download_file",
                    side_effect=lambda fid, fname, dest: dest / fname,
                ):
                    result = ingest_cnv(output_dir, raw_dir, "TCGA-BRCA")

            self.assertTrue((output_dir / "cnv.parquet").exists())
            self.assertEqual(result, output_dir / "cnv.parquet")

    def test_ingest_cnv_parquet_has_correct_schema(self):
        """ingest_cnv output parquet has correct column names."""
        from src.ingest_cnv import ingest_cnv

        manifest = self._make_manifest()

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = pathlib.Path(tmpdir) / "output"
            raw_dir = pathlib.Path(tmpdir) / "raw"
            output_dir.mkdir()
            raw_dir.mkdir()

            cnv_raw_dir = raw_dir / "cnv"
            cnv_raw_dir.mkdir()
            for entry in manifest:
                fpath = cnv_raw_dir / entry["file_name"]
                fpath.write_text(CNV_TSV_CONTENT)

            with patch("src.ingest_cnv.fetch_manifest", return_value=manifest):
                with patch(
                    "src.ingest_cnv.download_file",
                    side_effect=lambda fid, fname, dest: dest / fname,
                ):
                    ingest_cnv(output_dir, raw_dir, "TCGA-BRCA")

            df = pl.read_parquet(output_dir / "cnv.parquet")
            self.assertEqual(
                df.columns,
                ["patient_id", "chromosome", "start", "end", "copy_number"],
            )
            # 2 patients * 3 rows each
            self.assertEqual(len(df), 6)

    def test_ingest_cnv_skip_and_log_on_parse_error(self):
        """ingest_cnv writes errors_cnv.csv and continues when parse_cnv_seg raises."""
        from src.ingest_cnv import ingest_cnv

        manifest = [
            {
                "file_id": "file-id-001",
                "file_name": "good.seg.txt",
                "patient_id": "TCGA-AA-0001",
            },
            {
                "file_id": "file-id-bad",
                "file_name": "bad.seg.txt",
                "patient_id": "TCGA-AA-XXXX",
            },
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = pathlib.Path(tmpdir) / "output"
            raw_dir = pathlib.Path(tmpdir) / "raw"
            output_dir.mkdir()
            raw_dir.mkdir()

            cnv_raw_dir = raw_dir / "cnv"
            cnv_raw_dir.mkdir()

            # Good file
            (cnv_raw_dir / "good.seg.txt").write_text(CNV_TSV_CONTENT)
            # Bad file (malformed — no valid columns)
            (cnv_raw_dir / "bad.seg.txt").write_text("not\ta\tvalid\tfile\n1\t2\t3\t4\n")

            with patch("src.ingest_cnv.fetch_manifest", return_value=manifest):
                with patch(
                    "src.ingest_cnv.download_file",
                    side_effect=lambda fid, fname, dest: dest / fname,
                ):
                    ingest_cnv(output_dir, raw_dir, "TCGA-BRCA")

            # errors_cnv.csv must be written
            self.assertTrue((output_dir / "errors_cnv.csv").exists())
            # cnv.parquet must still exist (pipeline did not abort)
            self.assertTrue((output_dir / "cnv.parquet").exists())


if __name__ == "__main__":
    unittest.main()
