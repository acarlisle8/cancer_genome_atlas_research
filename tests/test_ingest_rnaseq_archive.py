"""Tests for src/ingest_rnaseq.py — RNA-seq ingestion and parsing."""
import csv
import pathlib
from unittest.mock import patch, MagicMock

import polars as pl
import pytest

from src.ingest_rnaseq import parse_rnaseq_tsv


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


