"""Tests for src/merge.py — pivot/join orchestrator for TCGA multi-omics merge."""
import pathlib
import tempfile
import unittest

import polars as pl


def _make_data_dir(tmp: pathlib.Path) -> pathlib.Path:
    """Write minimal synthetic Parquet files for all 3 cohorts x 3 modalities.

    Each cohort uses a distinct patient_id prefix to prevent cross-cohort
    collision in the merged output.  CNV rows include both a p-arm and q-arm
    segment per patient so _pivot_cnv produces non-empty output.
    """
    cohorts = ["TCGA-BRCA", "TCGA-LUAD", "TCGA-PRAD"]
    genes = [f"ENSG{i:011d}.1" for i in range(10)]
    probes = [f"cg{i:08d}" for i in range(10)]

    for cohort in cohorts:
        d = tmp / cohort
        d.mkdir(parents=True, exist_ok=True)

        # Distinct patient_id prefixes per cohort: BRCA -> BR, LUAD -> LU, PRAD -> PR
        prefix = cohort[5:7]
        patients = [f"TCGA-{prefix}-P{i:03d}" for i in range(3)]

        # rna_seq.parquet — 3 patients x 10 genes = 30 rows
        # Give each patient different fpkm values so variance > 0 across patients
        rna_rows = [
            {"patient_id": p, "gene_id": g, "fpkm_unstranded": float(i * 10 + j)}
            for i, p in enumerate(patients)
            for j, g in enumerate(genes)
        ]
        pl.DataFrame(rna_rows).write_parquet(d / "rna_seq.parquet")

        # cnv.parquet — 3 patients x 2 segments each = 6 rows
        # Segment 1: midpoint 3000000 < centromere 123400000 -> arm "1p"
        # Segment 2: midpoint 165000000 > centromere 123400000 -> arm "1q"
        cnv_rows = []
        for i, p in enumerate(patients):
            cnv_rows.append({
                "patient_id": p,
                "chromosome": "chr1",
                "start": 1000000,
                "end": 5000000,
                "copy_number": 0.1 * i,
            })
            cnv_rows.append({
                "patient_id": p,
                "chromosome": "chr1",
                "start": 130000000,
                "end": 200000000,
                "copy_number": -0.1 * i,
            })
        pl.DataFrame(cnv_rows).write_parquet(d / "cnv.parquet")

        # methylation.parquet — 3 patients x 10 probes = 30 rows
        meth_rows = [
            {"patient_id": p, "probe_id": pr, "beta_value": 0.1 * (i + j % 5)}
            for i, p in enumerate(patients)
            for j, pr in enumerate(probes)
        ]
        pl.DataFrame(meth_rows).write_parquet(d / "methylation.parquet")

    return tmp


class TestRnaseqPivot(unittest.TestCase):
    """MERGE-01: _pivot_rnaseq produces a wide matrix with patient_id x top gene columns."""

    def test_pivot_shape(self):
        """Output shape is (n_patients, n_selected_genes + 1)."""
        with tempfile.TemporaryDirectory() as tmp:
            d = pathlib.Path(tmp)
            genes = [f"ENSG{i:011d}.1" for i in range(10)]
            rows = [
                {"patient_id": f"P{i:03d}", "gene_id": g, "fpkm_unstranded": float(i + j)}
                for i in range(3)
                for j, g in enumerate(genes)
            ]
            rna_path = d / "rna_seq.parquet"
            pl.DataFrame(rows).write_parquet(rna_path)
            from src.merge import _pivot_rnaseq
            df = _pivot_rnaseq(rna_path, genes[:5])
            # 3 patients, 5 genes + patient_id column = 6 columns
            self.assertEqual(df.shape, (3, 6))

    def test_selected_genes_are_columns(self):
        """The correct gene IDs appear as columns in the wide output."""
        with tempfile.TemporaryDirectory() as tmp:
            d = pathlib.Path(tmp)
            genes = [f"ENSG{i:011d}.1" for i in range(10)]
            rows = [
                {"patient_id": f"P{i:03d}", "gene_id": g, "fpkm_unstranded": float(i + j)}
                for i in range(3)
                for j, g in enumerate(genes)
            ]
            rna_path = d / "rna_seq.parquet"
            pl.DataFrame(rows).write_parquet(rna_path)
            from src.merge import _pivot_rnaseq
            top_genes = genes[:3]
            df = _pivot_rnaseq(rna_path, top_genes)
            for g in top_genes:
                self.assertIn(g, df.columns)

    def test_excluded_genes_not_in_columns(self):
        """Non-selected genes do not appear as columns in the wide output."""
        with tempfile.TemporaryDirectory() as tmp:
            d = pathlib.Path(tmp)
            genes = [f"ENSG{i:011d}.1" for i in range(10)]
            rows = [
                {"patient_id": f"P{i:03d}", "gene_id": g, "fpkm_unstranded": float(i + j)}
                for i in range(3)
                for j, g in enumerate(genes)
            ]
            rna_path = d / "rna_seq.parquet"
            pl.DataFrame(rows).write_parquet(rna_path)
            from src.merge import _pivot_rnaseq
            top_genes = genes[:3]
            excluded = genes[3:]
            df = _pivot_rnaseq(rna_path, top_genes)
            for g in excluded:
                self.assertNotIn(g, df.columns)


class TestCnvPivot(unittest.TestCase):
    """MERGE-02: _pivot_cnv produces chromosome-arm columns named 1p, 1q, etc."""

    def test_arm_column_names(self):
        """Arm columns use NUMBERp/NUMBERq naming, not bare 'p'/'q'."""
        with tempfile.TemporaryDirectory() as tmp:
            cnv_path = pathlib.Path(tmp) / "cnv.parquet"
            rows = []
            for i in range(3):
                # p-arm segment: midpoint 3000000 < centromere 123400000
                rows.append({
                    "patient_id": f"P{i:03d}",
                    "chromosome": "chr1",
                    "start": 1000000,
                    "end": 5000000,
                    "copy_number": 0.1 * i,
                })
                # q-arm segment: midpoint 165000000 > centromere 123400000
                rows.append({
                    "patient_id": f"P{i:03d}",
                    "chromosome": "chr1",
                    "start": 130000000,
                    "end": 200000000,
                    "copy_number": -0.1 * i,
                })
            pl.DataFrame(rows).write_parquet(cnv_path)
            from src.merge import _pivot_cnv
            df = _pivot_cnv(cnv_path)
            self.assertIn("1p", df.columns)
            self.assertIn("1q", df.columns)
            # bare "p" or "q" would indicate wrong column naming
            self.assertNotIn("p", df.columns)
            self.assertNotIn("q", df.columns)

    def test_midpoint_arm_assignment(self):
        """Segments are assigned to the correct arm based on midpoint vs centromere."""
        with tempfile.TemporaryDirectory() as tmp:
            cnv_path = pathlib.Path(tmp) / "cnv.parquet"
            rows = [
                # midpoint = (1000000 + 5000000) / 2 = 3000000 < 123400000 -> "1p"
                {
                    "patient_id": "P000",
                    "chromosome": "chr1",
                    "start": 1000000,
                    "end": 5000000,
                    "copy_number": 0.5,
                },
                # midpoint = (130000000 + 200000000) / 2 = 165000000 > 123400000 -> "1q"
                {
                    "patient_id": "P000",
                    "chromosome": "chr1",
                    "start": 130000000,
                    "end": 200000000,
                    "copy_number": -0.3,
                },
            ]
            pl.DataFrame(rows).write_parquet(cnv_path)
            from src.merge import _pivot_cnv
            df = _pivot_cnv(cnv_path)
            row = df.filter(pl.col("patient_id") == "P000")
            self.assertIsNotNone(row["1p"][0])
            self.assertIsNotNone(row["1q"][0])

    def test_noncanonical_chromosomes_filtered(self):
        """chrM segments produce no arm columns (filtered before centromere join)."""
        with tempfile.TemporaryDirectory() as tmp:
            cnv_path = pathlib.Path(tmp) / "cnv.parquet"
            rows = [
                # canonical segment — provides a real arm so pivot has output
                {
                    "patient_id": "P000",
                    "chromosome": "chr1",
                    "start": 1000000,
                    "end": 5000000,
                    "copy_number": 0.1,
                },
                # non-canonical — should be filtered out
                {
                    "patient_id": "P000",
                    "chromosome": "chrM",
                    "start": 1000,
                    "end": 5000,
                    "copy_number": 0.2,
                },
            ]
            pl.DataFrame(rows).write_parquet(cnv_path)
            from src.merge import _pivot_cnv
            df = _pivot_cnv(cnv_path)
            self.assertNotIn("Mp", df.columns)
            self.assertNotIn("Mq", df.columns)

    def test_one_row_per_patient(self):
        """Pivot produces exactly one row per patient_id."""
        with tempfile.TemporaryDirectory() as tmp:
            cnv_path = pathlib.Path(tmp) / "cnv.parquet"
            rows = []
            for i in range(3):
                rows.append({
                    "patient_id": f"P{i:03d}",
                    "chromosome": "chr1",
                    "start": 1000000,
                    "end": 5000000,
                    "copy_number": 0.1 * i,
                })
                rows.append({
                    "patient_id": f"P{i:03d}",
                    "chromosome": "chr1",
                    "start": 130000000,
                    "end": 200000000,
                    "copy_number": -0.1 * i,
                })
            pl.DataFrame(rows).write_parquet(cnv_path)
            from src.merge import _pivot_cnv
            df = _pivot_cnv(cnv_path)
            self.assertEqual(df["patient_id"].n_unique(), len(df))


class TestMethylationPivot(unittest.TestCase):
    """MERGE-03: _pivot_methylation produces wide matrix with selected probe columns; null beta handled."""

    def test_pivot_shape(self):
        """Output shape is (n_patients, n_selected_probes + 1)."""
        with tempfile.TemporaryDirectory() as tmp:
            d = pathlib.Path(tmp)
            probes = [f"cg{i:08d}" for i in range(10)]
            rows = [
                {"patient_id": f"P{i:03d}", "probe_id": pr, "beta_value": 0.1 * (i + j % 5)}
                for i in range(3)
                for j, pr in enumerate(probes)
            ]
            meth_path = d / "methylation.parquet"
            pl.DataFrame(rows).write_parquet(meth_path)
            from src.merge import _pivot_methylation
            df = _pivot_methylation(meth_path, probes[:5])
            # 3 patients, 5 probes + patient_id column = 6 columns
            self.assertEqual(df.shape, (3, 6))

    def test_null_beta_value_handled(self):
        """null beta_value does not crash the pivot; null cell remains null."""
        with tempfile.TemporaryDirectory() as tmp:
            meth_path = pathlib.Path(tmp) / "methylation.parquet"
            rows = [
                {"patient_id": "P000", "probe_id": "cg00000000", "beta_value": None},
                {"patient_id": "P001", "probe_id": "cg00000000", "beta_value": 0.5},
            ]
            pl.DataFrame(
                rows,
                schema={"patient_id": pl.Utf8, "probe_id": pl.Utf8, "beta_value": pl.Float64},
            ).write_parquet(meth_path)
            from src.merge import _pivot_methylation
            # should not raise
            df = _pivot_methylation(meth_path, ["cg00000000"])
            p000_row = df.filter(pl.col("patient_id") == "P000")
            # null beta remains null (not coerced to 0.0)
            self.assertIsNone(p000_row["cg00000000"][0])

    def test_excluded_probes_absent(self):
        """Probes not in top_probes do not appear as columns."""
        with tempfile.TemporaryDirectory() as tmp:
            meth_path = pathlib.Path(tmp) / "methylation.parquet"
            probes = [f"cg{i:08d}" for i in range(10)]
            rows = [
                {"patient_id": f"P{i:03d}", "probe_id": pr, "beta_value": 0.1 * (i + j % 5)}
                for i in range(3)
                for j, pr in enumerate(probes)
            ]
            pl.DataFrame(rows).write_parquet(meth_path)
            from src.merge import _pivot_methylation
            top_probes = probes[:3]
            excluded = probes[3:]
            df = _pivot_methylation(meth_path, top_probes)
            for pr in excluded:
                self.assertNotIn(pr, df.columns)


class TestMergeAllCohorts(unittest.TestCase):
    """MERGE-04: merge_all_cohorts produces correct shape, schema, and inner-join semantics."""

    def test_output_file_exists(self):
        """Output file exists and is named merged_all_cohorts.parquet."""
        with tempfile.TemporaryDirectory() as tmp:
            data_dir = _make_data_dir(pathlib.Path(tmp))
            from src.merge import merge_all_cohorts
            out = merge_all_cohorts(data_dir, data_dir)
            self.assertTrue(out.exists())
            self.assertEqual(out.name, "merged_all_cohorts.parquet")

    def test_one_row_per_patient(self):
        """Merged output has exactly one row per unique patient_id."""
        with tempfile.TemporaryDirectory() as tmp:
            data_dir = _make_data_dir(pathlib.Path(tmp))
            from src.merge import merge_all_cohorts
            out = merge_all_cohorts(data_dir, data_dir)
            df = pl.read_parquet(out)
            self.assertEqual(df["patient_id"].n_unique(), len(df))

    def test_cohort_column_values(self):
        """cohort column has exactly the three bare cancer type labels."""
        with tempfile.TemporaryDirectory() as tmp:
            data_dir = _make_data_dir(pathlib.Path(tmp))
            from src.merge import merge_all_cohorts
            out = merge_all_cohorts(data_dir, data_dir)
            df = pl.read_parquet(out)
            self.assertIn("cohort", df.columns)
            self.assertEqual(set(df["cohort"].to_list()), {"BRCA", "LUAD", "PRAD"})

    def test_inner_join_excludes_missing_patients(self):
        """A patient present in rna_seq but absent from cnv must not appear in output."""
        with tempfile.TemporaryDirectory() as tmp:
            data_dir = _make_data_dir(pathlib.Path(tmp))
            # Identify one BRCA patient and remove them from the CNV file
            brca_rna = pl.read_parquet(data_dir / "TCGA-BRCA" / "rna_seq.parquet")
            all_brca_patients = brca_rna["patient_id"].unique().to_list()
            excluded = all_brca_patients[-1]
            brca_cnv_full = pl.read_parquet(data_dir / "TCGA-BRCA" / "cnv.parquet")
            brca_cnv_trimmed = brca_cnv_full.filter(pl.col("patient_id") != excluded)
            brca_cnv_trimmed.write_parquet(data_dir / "TCGA-BRCA" / "cnv.parquet")
            from src.merge import merge_all_cohorts
            out = merge_all_cohorts(data_dir, data_dir)
            df = pl.read_parquet(out)
            brca_patients_in_output = (
                df.filter(pl.col("cohort") == "BRCA")["patient_id"].to_list()
            )
            self.assertNotIn(excluded, brca_patients_in_output)


if __name__ == "__main__":
    unittest.main()
