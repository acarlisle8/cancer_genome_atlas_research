"""Merge: pivot Phase 1 Parquets and join into a single wide feature matrix."""
import pathlib

import polars as pl

from src.utils import get_logger

logger = get_logger(__name__)

COHORTS = ["TCGA-BRCA", "TCGA-LUAD", "TCGA-PRAD"]
N_RNA_GENES = 500
N_METH_PROBES = 500

# Ensembl base IDs (no version suffix) for canonical subtype marker genes.
# These are pinned into the feature set regardless of variance rank so that
# downstream subtype analyses always have the relevant signal available.
# Sources: PAM50/BRCA (Parker 2009, TCGA 2012), LUAD (Wilkerson 2012,
# TCGA 2014), PRAD iCluster (TCGA Cell 2015).
SUBTYPE_MARKER_GENES: frozenset[str] = frozenset({
    # BRCA PAM50 subset
    "ENSG00000091831",  # ESR1
    "ENSG00000082175",  # PGR
    "ENSG00000141736",  # ERBB2
    "ENSG00000148773",  # MKI67
    "ENSG00000129514",  # FOXA1  (also PRAD)
    "ENSG00000054598",  # FOXC1
    "ENSG00000171791",  # BCL2
    "ENSG00000039068",  # CDH1
    "ENSG00000136997",  # MYC
    "ENSG00000089685",  # BIRC5
    "ENSG00000087586",  # AURKA
    "ENSG00000134057",  # CCNB1
    "ENSG00000101057",  # MYBL2
    "ENSG00000186081",  # KRT5  (also LUAD)
    "ENSG00000186847",  # KRT14
    "ENSG00000128422",  # KRT17
    "ENSG00000146648",  # EGFR
    "ENSG00000160867",  # FGFR4
    "ENSG00000117399",  # CDC20
    "ENSG00000165304",  # MELK
    # LUAD markers
    "ENSG00000136352",  # NKX2-1
    "ENSG00000168878",  # SFTPB
    "ENSG00000168484",  # SFTPC
    "ENSG00000073282",  # TP63
    "ENSG00000205420",  # KRT6A
    "ENSG00000185499",  # MUC1
    # PRAD markers
    "ENSG00000169083",  # AR
    "ENSG00000157554",  # ERG
    "ENSG00000006468",  # ETV1
    "ENSG00000175832",  # ETV4
    "ENSG00000151702",  # FLI1
    "ENSG00000121067",  # SPOP
    "ENSG00000171862",  # PTEN
    "ENSG00000141510",  # TP53
})

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
    n: int = 500,
    pinned_base_ids: frozenset[str] | None = None,
) -> list[str]:
    """Return top-N feature IDs by variance, always including pinned features.

    Pinned features are matched by Ensembl base ID (stripping the version
    suffix, e.g. ENSG00000091831.16 -> ENSG00000091831) and consume slots
    from n, so the total returned is always <= n.

    Args:
        parquet_paths: Paths to Parquet files to scan (all cohorts).
        id_col: Column name for the feature identifier.
        value_col: Column name for the numeric value.
        n: Total number of features to return (default 500).
        pinned_base_ids: Ensembl base IDs to force-include regardless of rank.

    Returns:
        List of feature ID strings (versioned if present in data).
    """
    lf = pl.scan_parquet([str(p) for p in parquet_paths])
    all_stats = (
        lf
        .group_by(id_col)
        .agg(pl.col(value_col).var().alias("variance"))
        .collect()
    )

    if pinned_base_ids:
        all_stats = all_stats.with_columns(
            pl.col(id_col).str.split(".").list.first().alias("_base_id")
        )
        pinned = all_stats.filter(pl.col("_base_id").is_in(pinned_base_ids))
        pinned_versioned: set[str] = set(pinned[id_col].to_list())
        n_fill = max(0, n - len(pinned_versioned))
        top_var = (
            all_stats
            .filter(~pl.col(id_col).is_in(pinned_versioned))
            .sort("variance", descending=True, nulls_last=True)
            .head(n_fill)
        )
        return pinned[id_col].to_list() + top_var[id_col].to_list()

    return (
        all_stats
        .sort("variance", descending=True, nulls_last=True)
        .head(n)[id_col]
        .to_list()
    )


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
    df = pl.read_parquet(parquet_path)
    filtered = df.filter(pl.col("gene_id").is_in(top_genes))
    return filtered.pivot(
        on="gene_id",
        index="patient_id",
        values="fpkm_unstranded",
        aggregate_function="mean",
    )


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
    df = pl.read_parquet(parquet_path)

    # Drop non-canonical chromosomes (chrM, chrUn_*, etc.) before centromere join
    df = df.filter(pl.col("chromosome").is_in(list(HG38_CENTROMERES.keys())))

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
    # Sort arm columns alphabetically so pl.concat across cohorts produces consistent schemas.
    # group_by output order is non-deterministic; pivot inherits that order.
    arm_cols = sorted(c for c in wide.columns if c != "patient_id")
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

    logger.info("Computing global top-%d genes by variance across %d cohorts", N_RNA_GENES, len(COHORTS))
    top_genes = _top_n_by_variance(
        rna_paths, "gene_id", "fpkm_unstranded", N_RNA_GENES, SUBTYPE_MARKER_GENES
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
        cnv_wide = _pivot_cnv(cohort_dir / "cnv.parquet")
        meth_wide = _pivot_methylation(cohort_dir / "methylation.parquet", top_probes)

        merged_cohort = (
            rna_wide
            .join(cnv_wide, on="patient_id", how="inner")
            .join(meth_wide, on="patient_id", how="inner")
            .with_columns(pl.lit(label).alias("cohort"))
        )
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
