"""Tests for src/ingest_rnaseq.py — RNA-seq ingestion and parsing."""
import csv
import pathlib
from unittest.mock import patch, MagicMock

import polars as pl
import pytest

from src.ingest_rnaseq import parse_rnaseq_tsv, ingest_rnaseq


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

TSV_HEADER = (
    "gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second"
    "\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded\n"
)

TSV_DATA_ROWS = [
    "ENSG00000000003.15\tTSPAN6\tprotein_coding\t1200\t600\t600\t5.1\t2.3\t3.1\n",
    "ENSG00000000005.6\tTNMD\tprotein_coding\t0\t0\t0\t0.0\t0.0\t0.0\n",
    "ENSG00000000419.12\tDPM1\tprotein_coding\t800\t400\t400\t3.2\t1.1\t2.0\n",
    "ENSG00000000457.14\tSCYL3\tprotein_coding\t300\t150\t150\t1.5\t0.7\t1.0\n",
    "ENSG00000000460.17\tC1orf112\tprotein_coding\t500\t250\t250\t2.0\t0.9\t1.3\n",
]

TSV_N_ROWS = [
    "N_unmapped\t\t\t100\t50\t50\t0.0\t0.0\t0.0\n",
    "N_multimapping\t\t\t200\t100\t100\t0.0\t0.0\t0.0\n",
    "N_noFeature\t\t\t50\t25\t25\t0.0\t0.0\t0.0\n",
    "N_ambiguous\t\t\t30\t15\t15\t0.0\t0.0\t0.0\n",
]


def make_tsv_file(tmp_path: pathlib.Path, filename: str = "sample.tsv") -> pathlib.Path:
    """Write a minimal TSV fixture with data rows and N_ summary rows."""
    content = TSV_HEADER + "".join(TSV_DATA_ROWS) + "".join(TSV_N_ROWS)
    p = tmp_path / filename
    p.write_text(content)
    return p


# ---------------------------------------------------------------------------
# Tests: parse_rnaseq_tsv
# ---------------------------------------------------------------------------


class TestParseRnaseqTsv:
    def test_returns_correct_columns(self, tmp_path):
        tsv = make_tsv_file(tmp_path)
        df = parse_rnaseq_tsv(tsv, "TCGA-BH-A18H")
        assert df.columns == ["patient_id", "gene_id", "fpkm_unstranded"]

    def test_returns_correct_dtypes(self, tmp_path):
        tsv = make_tsv_file(tmp_path)
        df = parse_rnaseq_tsv(tsv, "TCGA-BH-A18H")
        assert df.dtypes[0] == pl.Utf8
        assert df.dtypes[1] == pl.Utf8
        assert df.dtypes[2] == pl.Float64

    def test_filters_n_rows(self, tmp_path):
        """N_ summary rows must be excluded; only 5 data rows remain."""
        tsv = make_tsv_file(tmp_path)
        df = parse_rnaseq_tsv(tsv, "TCGA-BH-A18H")
        assert len(df) == 5

    def test_no_n_prefix_in_gene_id(self, tmp_path):
        tsv = make_tsv_file(tmp_path)
        df = parse_rnaseq_tsv(tsv, "TCGA-BH-A18H")
        assert not any(gid.startswith("N_") for gid in df["gene_id"].to_list())

    def test_patient_id_column_value(self, tmp_path):
        tsv = make_tsv_file(tmp_path)
        df = parse_rnaseq_tsv(tsv, "TCGA-BH-A18H")
        assert all(pid == "TCGA-BH-A18H" for pid in df["patient_id"].to_list())

    def test_fpkm_values_match(self, tmp_path):
        tsv = make_tsv_file(tmp_path)
        df = parse_rnaseq_tsv(tsv, "TCGA-BH-A18H")
        expected = [2.3, 0.0, 1.1, 0.7, 0.9]
        actual = df["fpkm_unstranded"].to_list()
        assert actual == pytest.approx(expected)

    def test_only_three_columns_returned(self, tmp_path):
        """Extra columns from TSV must not appear in output."""
        tsv = make_tsv_file(tmp_path)
        df = parse_rnaseq_tsv(tsv, "TCGA-BH-A18H")
        assert len(df.columns) == 3


# ---------------------------------------------------------------------------
# Tests: ingest_rnaseq
# ---------------------------------------------------------------------------


class TestIngestRnaseq:
    def _make_manifest_entries(self, tmp_path: pathlib.Path):
        """Create two TSV files and return manifest-style entries pointing to them."""
        p1 = make_tsv_file(tmp_path, "patient1.tsv")
        p2 = make_tsv_file(tmp_path, "patient2.tsv")
        return [
            {"file_id": "file-001", "file_name": "patient1.tsv", "patient_id": "TCGA-BH-A18H"},
            {"file_id": "file-002", "file_name": "patient2.tsv", "patient_id": "TCGA-E2-A14P"},
        ], p1, p2

    def test_writes_parquet(self, tmp_path):
        entries, p1, p2 = self._make_manifest_entries(tmp_path)
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        with patch("src.ingest_rnaseq.fetch_manifest", return_value=entries), \
             patch("src.ingest_rnaseq.download_file", side_effect=[p1, p2]):
            result = ingest_rnaseq(output_dir, tmp_path / "raw")

        assert result == output_dir / "rna_seq.parquet"
        assert result.exists()

    def test_parquet_schema(self, tmp_path):
        entries, p1, p2 = self._make_manifest_entries(tmp_path)
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        with patch("src.ingest_rnaseq.fetch_manifest", return_value=entries), \
             patch("src.ingest_rnaseq.download_file", side_effect=[p1, p2]):
            result = ingest_rnaseq(output_dir, tmp_path / "raw")

        df = pl.read_parquet(result)
        assert "patient_id" in df.columns
        assert "gene_id" in df.columns
        assert "fpkm_unstranded" in df.columns
        assert df.dtypes[df.columns.index("patient_id")] == pl.Utf8
        assert df.dtypes[df.columns.index("gene_id")] == pl.Utf8
        assert df.dtypes[df.columns.index("fpkm_unstranded")] == pl.Float64

    def test_parquet_row_count(self, tmp_path):
        """2 patients × 5 data rows each = 10 total rows."""
        entries, p1, p2 = self._make_manifest_entries(tmp_path)
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        with patch("src.ingest_rnaseq.fetch_manifest", return_value=entries), \
             patch("src.ingest_rnaseq.download_file", side_effect=[p1, p2]):
            result = ingest_rnaseq(output_dir, tmp_path / "raw")

        df = pl.read_parquet(result)
        assert len(df) == 10

    def test_fetch_manifest_called_with_correct_args(self, tmp_path):
        entries, p1, p2 = self._make_manifest_entries(tmp_path)
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        with patch("src.ingest_rnaseq.fetch_manifest", return_value=entries) as mock_fm, \
             patch("src.ingest_rnaseq.download_file", side_effect=[p1, p2]):
            ingest_rnaseq(output_dir, tmp_path / "raw", project_id="TCGA-BRCA")

        mock_fm.assert_called_once_with(
            "TCGA-BRCA",
            "Transcriptome Profiling",
            "Gene Expression Quantification",
        )

    def test_error_recorded_in_csv_on_parse_failure(self, tmp_path):
        """When parse_rnaseq_tsv raises, errors_rnaseq.csv is written; pipeline continues."""
        good_tsv = make_tsv_file(tmp_path, "good.tsv")
        bad_tsv = tmp_path / "bad.tsv"
        bad_tsv.write_text("not\ta\tvalid\ttsv\n")

        entries = [
            {"file_id": "file-001", "file_name": "good.tsv", "patient_id": "TCGA-BH-A18H"},
            {"file_id": "file-002", "file_name": "bad.tsv", "patient_id": "TCGA-XX-XXXX"},
        ]
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        with patch("src.ingest_rnaseq.fetch_manifest", return_value=entries), \
             patch("src.ingest_rnaseq.download_file", side_effect=[good_tsv, bad_tsv]):
            result = ingest_rnaseq(output_dir, tmp_path / "raw")

        # Pipeline should still produce a Parquet from the good patient
        assert result.exists()
        df = pl.read_parquet(result)
        assert len(df) == 5  # only the good patient's rows

        # Error CSV must exist and contain the bad patient
        errors_csv = output_dir / "errors_rnaseq.csv"
        assert errors_csv.exists()
        with open(errors_csv, newline="") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == 1
        assert rows[0]["patient_id"] == "TCGA-XX-XXXX"
        assert rows[0]["file_id"] == "file-002"
        assert "reason" in rows[0]

    def test_no_error_csv_when_no_failures(self, tmp_path):
        """errors_rnaseq.csv should NOT be created when all parses succeed."""
        entries, p1, p2 = self._make_manifest_entries(tmp_path)
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        with patch("src.ingest_rnaseq.fetch_manifest", return_value=entries), \
             patch("src.ingest_rnaseq.download_file", side_effect=[p1, p2]):
            ingest_rnaseq(output_dir, tmp_path / "raw")

        assert not (output_dir / "errors_rnaseq.csv").exists()
