"""Merge: pivot Phase 1 Parquets and join into a single wide feature matrix."""
import pathlib

import duckdb
import polars as pl

from src.markers import ALL_RNA_MARKERS
from src.utils import get_logger

logger = get_logger(__name__)

COHORTS = ["TCGA-BRCA", "TCGA-LUAD", "TCGA-PRAD"]

# Variance-based pre-selection caps. Wider than strictly needed to give the
# training pipeline headroom for in-fold feature selection without re-running
# the merge. RNA has ~60K genes and methylation has ~850K probes upstream;
# capping at 5K each yields a merged matrix of ~10K columns — manageable on
# disk and in memory while preserving most of the informative variance.
N_RNA_GENES = 5000
N_METH_PROBES = 5000

HG38_CENTROMERES = {
    "chr1": 123400000, "chr2": 93900000, "chr3": 90900000,
    "chr4": 50400000, "chr5": 48800000, "chr6": 59800000,
    "chr7": 60100000, "chr8": 45200000, "chr9": 43300000,
    "chr10": 39800000, "chr11": 53400000, "chr12": 35500000,
    "chr13": 17700000, "chr14": 17200000, "chr15": 19000000,
    "chr16": 36800000, "chr17": 25100000, "chr18": 18500000,
    "chr19": 26200000, "chr20": 28100000, "chr21": 12000000,
    "chr22": 15000000, "chrX": 60600000,
}

# Pre-computed set of canonical chromosome names (from HG38_CENTROMERES) for fast lookup
_CANONICAL_CHROMS: set[str] = set(HG38_CENTROMERES.keys())


def _top_n_by_variance(
    parquet_paths: list[pathlib.Path],
    id_col: str,
    value_col: str,
    n: int = 5000,
    must_keep: set[str] | None = None,
) -> list[str]:
    """Return the top-N feature IDs by variance, plus any must_keep features.

    Uses lazy scan to avoid loading all rows into memory simultaneously. For
    versioned Ensembl IDs, must_keep contains bare IDs (e.g., "ENSG00000091831")
    and is matched against the parquet's versioned column (e.g.,
    "ENSG00000091831.16") via split on ".".

    Args:
        parquet_paths: Paths to Parquet files to scan (all cohorts).
        id_col: Column name for the feature identifier (e.g. "gene_id", "probe_id").
        value_col: Column name for the numeric value (e.g. "fpkm_unstranded", "beta_value").
        n: Number of top features to return (default 5000).
        must_keep: Optional set of bare feature IDs to force-include regardless of
            variance rank. Matched against the parquet column via split on ".",
            so versioned IDs in the parquet still match bare IDs in this set.

    Returns:
        List of feature ID strings: top-N by variance, unioned with must_keep
        members that actually exist in the parquets. Length is between n and
        n + |must_keep|.
    """
    must_keep = must_keep or set()
    paths_sql = ", ".join(f"'{p}'" for p in parquet_paths)
    con = duckdb.connect(":memory:")
    try:
        con.execute("PRAGMA memory_limit='4GB'")
        con.execute("PRAGMA threads=2")
        rows = con.execute(f"""
            SELECT {id_col}
            FROM read_parquet([{paths_sql}])
            WHERE {value_col} IS NOT NULL
            GROUP BY {id_col}
            ORDER BY VAR_POP({value_col}) DESC NULLS LAST
            LIMIT {n}
        """).fetchall()
        selected = {r[0] for r in rows}

        if must_keep:
            keeps_sql = ", ".join(f"'{k}'" for k in must_keep)
            keep_rows = con.execute(f"""
                SELECT DISTINCT {id_col}
                FROM read_parquet([{paths_sql}])
                WHERE split_part({id_col}, '.', 1) IN ({keeps_sql})
            """).fetchall()
            keep_versioned = {r[0] for r in keep_rows}
            added = keep_versioned - selected
            if added:
                logger.info(
                    "Force-included %d must_keep features (of %d candidates) "
                    "not already in top-%d variance",
                    len(added), len(must_keep), n,
                )
            selected |= keep_versioned
    finally:
        con.close()
    return list(selected)


def _pivot_rnaseq(parquet_path: pathlib.Path, top_genes: list[str]) -> pl.DataFrame:
    """Pivot long-format RNA-seq Parquet to wide format (patient × gene).

    Filters to the pre-selected global gene list, then pivots so each gene
    becomes a column. Uses aggregate_function='mean' to handle rare duplicate
    (patient_id, gene_id) pairs.

    Args:
        parquet_path: Path to rna_seq.parquet (long format).
        top_genes: List of gene_id strings to keep as columns.

    Returns:
        Wide DataFrame with columns: patient_id, <gene_id_1>, ..., <gene_id_n>.
    """
    # Lab-05 memory pattern: lazy scan + filter pushdown + project only the
    # columns the pivot actually needs. RNA-seq parquets are 0.7–1.5 GB
    # compressed, so eager pl.read_parquet OOMs an 8 GB instance; this
    # collects only ~547K rows (1095 patients × 500 top genes) instead.
    filtered = (
        pl.scan_parquet(parquet_path)
        .filter(pl.col("gene_id").is_in(top_genes))
        .select(["patient_id", "gene_id", "fpkm_uq_unstranded"])
        .collect()
    )
    return filtered.pivot(
        on="gene_id",
        index="patient_id",
        values="fpkm_uq_unstranded",
        aggregate_function="mean",
    )


def _arm_sort_key(arm: str) -> tuple[int, int]:
    """Natural chromosome order: 1p < 1q < 2p < ... < 22q < Xp < Xq."""
    chrom, side = arm[:-1], arm[-1]
    chrom_num = 23 if chrom == "X" else int(chrom)
    side_num = 0 if side == "p" else 1
    return (chrom_num, side_num)


def _pivot_cnv(parquet_path: pathlib.Path) -> pl.DataFrame:
    """Aggregate CNV segments to chromosome-arm means and pivot wide.

    Maps each segment to a chromosome arm (p or q) using hg38 centromere
    positions and the segment midpoint heuristic. Computes mean log2 copy
    ratio per patient per arm, then pivots to wide format.

    Non-canonical chromosomes (chrM, chrUn_*, etc.) are filtered out before
    the centromere join to avoid null chr_arm values.

    Arm column names follow the pattern: 1p, 1q, 2p, ..., Xp, Xq.

    Args:
        parquet_path: Path to cnv.parquet (long segment format).

    Returns:
        Wide DataFrame with columns: patient_id, <arm_1>, ..., <arm_n>.
    """
    # GDC seg files emit bare "1"/"X"; synthetic test data may use "chr1". Normalize to "chr"-prefix form.
    lf = pl.scan_parquet(parquet_path).with_columns(
        pl.when(pl.col("chromosome").str.starts_with("chr"))
        .then(pl.col("chromosome"))
        .otherwise(pl.lit("chr") + pl.col("chromosome"))
        .alias("chromosome")
    )

    df = lf.filter(pl.col("chromosome").is_in(list(_CANONICAL_CHROMS))).collect()

    if df.is_empty():
        raise RuntimeError(
            f"_pivot_cnv: no rows survived chromosome filter on {parquet_path} — "
            f"input parquet has no canonical chromosome values"
        )

    cent_df = pl.DataFrame({
        "chromosome": list(HG38_CENTROMERES.keys()),
        "centromere": list(HG38_CENTROMERES.values()),
    })

    df = (
        df
        .join(cent_df, on="chromosome", how="left")
        .with_columns([
            (
                pl.col("chromosome").str.replace("chr", "")
                + pl.when(
                    ((pl.col("start") + pl.col("end")) / 2) < pl.col("centromere")
                )
                .then(pl.lit("p"))
                .otherwise(pl.lit("q"))
            ).alias("chr_arm")
        ])
    )

    arm_means = (
        df
        .group_by(["patient_id", "chr_arm"])
        .agg(pl.col("copy_number").mean())
    )

    wide = arm_means.pivot(
        on="chr_arm",
        index="patient_id",
        values="copy_number",
        aggregate_function="first",
    )

    arm_cols = sorted((c for c in wide.columns if c != "patient_id"), key=_arm_sort_key)
    return wide.select(["patient_id"] + arm_cols)


def _pivot_methylation(parquet_path: pathlib.Path, top_probes: list[str]) -> pl.DataFrame:
    """Pivot long-format methylation Parquet to wide format (patient × probe).

    Uses lazy scan + filter before collect to avoid loading the full 225M-row
    methylation table into memory. Only the top-N probes are collected eagerly.

    Args:
        parquet_path: Path to methylation.parquet (long format).
        top_probes: List of probe_id strings to keep as columns.

    Returns:
        Wide DataFrame with columns: patient_id, <probe_id_1>, ..., <probe_id_n>.
    """
    filtered = (
        pl.scan_parquet(parquet_path)
        .filter(pl.col("probe_id").is_in(top_probes))
        .collect()
    )
    return filtered.pivot(
        on="probe_id",
        index="patient_id",
        values="beta_value",
        aggregate_function="mean",
    )


def merge_all_cohorts(data_dir: pathlib.Path, output_dir: pathlib.Path) -> pathlib.Path:
    """Merge Phase 1 Parquets across all cohorts into a single wide feature matrix.

    Reads nine Parquet files from data_dir/{cohort}/{modality}.parquet, computes
    global top-500 gene and probe lists by variance, pivots each modality per cohort,
    inner-joins the three modalities on patient_id, adds a cohort label column, and
    stacks all cohorts into a single merged Parquet file.

    Input layout:
        data_dir/TCGA-BRCA/rna_seq.parquet
        data_dir/TCGA-BRCA/cnv.parquet
        data_dir/TCGA-BRCA/methylation.parquet
        data_dir/TCGA-LUAD/...
        data_dir/TCGA-PRAD/...

    Output schema:
        patient_id:  VARCHAR  — 12-char TCGA barcode (one row per patient)
        <gene_ids>:  DOUBLE   — 500 RNA-seq FPKM columns (global top-500 by variance)
        <arm_cols>:  DOUBLE   — ~46 chromosome-arm copy-number columns (1p, 1q, ..., Xq)
        <probe_ids>: DOUBLE   — 500 methylation beta-value columns (global top-500 by variance)
        cohort:      VARCHAR  — bare cancer type: "BRCA", "LUAD", or "PRAD"

    Join semantics: inner join across all three modalities — patients absent from any
    one modality are excluded from the output.

    Args:
        data_dir: Root directory containing per-cohort subdirectories with Phase 1 Parquets.
        output_dir: Directory where merged_all_cohorts.parquet is written.

    Returns:
        pathlib.Path to the written merged_all_cohorts.parquet.

    Raises:
        RuntimeError: If the merged matrix is empty (all patients excluded by inner join
            or Phase 1 Parquets are missing/empty).
        RuntimeError: If duplicate patient_id rows are found in the merged output.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Global feature selection: variance computed across all cohorts combined
    rna_paths = [data_dir / c / "rna_seq.parquet" for c in COHORTS]
    meth_paths = [data_dir / c / "methylation.parquet" for c in COHORTS]

    logger.info(
        "Computing global top-%d genes by variance + %d marker panel genes",
        N_RNA_GENES, len(ALL_RNA_MARKERS),
    )
    top_genes = _top_n_by_variance(
        rna_paths, "gene_id", "fpkm_uq_unstranded",
        N_RNA_GENES, must_keep=ALL_RNA_MARKERS,
    )

    logger.info("Computing global top-%d probes by variance across %d cohorts", N_METH_PROBES, len(COHORTS))
    top_probes = _top_n_by_variance(meth_paths, "probe_id", "beta_value", N_METH_PROBES)

    # Per-cohort transform and join
    cohort_frames = []
    for cohort in COHORTS:
        label = cohort.replace("TCGA-", "")  # "BRCA", "LUAD", "PRAD"
        logger.info("Processing cohort %s", cohort)
        cohort_dir = data_dir / cohort

        rna_wide = _pivot_rnaseq(cohort_dir / "rna_seq.parquet", top_genes)
        logger.info("  rna_wide  %s", rna_wide.shape)
        cnv_wide = _pivot_cnv(cohort_dir / "cnv.parquet")
        logger.info("  cnv_wide  %s", cnv_wide.shape)
        meth_wide = _pivot_methylation(cohort_dir / "methylation.parquet", top_probes)
        logger.info("  meth_wide %s", meth_wide.shape)

        merged_cohort = (
            rna_wide
            .join(cnv_wide, on="patient_id", how="inner")
            .join(meth_wide, on="patient_id", how="inner")
            .with_columns(pl.lit(label).alias("cohort"))
        )
        logger.info("  merged_cohort %s (after 3-way inner join)", merged_cohort.shape)
        cohort_frames.append(merged_cohort)

    # Stack all cohorts — diagonal_relaxed fills null for any arm columns absent in a cohort
    # (e.g. cohort with no chrX segments won't have Xp/Xq; vertical would crash)
    final = pl.concat(cohort_frames, how="diagonal_relaxed")

    # Post-condition guards
    if len(final) == 0:
        raise RuntimeError(
            "Merged matrix is empty — check that Phase 1 Parquets exist and are non-empty."
        )
    if final["patient_id"].n_unique() != len(final):
        raise RuntimeError("Duplicate patient_id rows in merged output")

    # Write output
    out_path = output_dir / "merged_all_cohorts.parquet"
    final.write_parquet(out_path, compression="snappy")
    logger.info("Wrote merged matrix (%d rows x %d cols) to %s", len(final), len(final.columns), out_path)

    return out_path
