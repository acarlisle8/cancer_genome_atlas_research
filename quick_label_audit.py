"""
subtype_label_audit.py

Audits TCGA subtype labels and subtype-marker gene presence for a model-ready
multi-omic cohort (BRCA, LUAD, PRAD).

Inputs:
  --cohort-parquet : the existing model_ready_cohort.parquet
  --hoadley-tsv    : Hoadley pan-cancer subtype table, downloaded manually from
                     https://gdc.cancer.gov/about-data/publications/panimmune
                     (filename: TCGASubtype.20170308.tsv)
  --out-dir        : where to write audit artifacts
  --cbio-cache     : local cache dir for cBioPortal downloads (default: data/cbio_cache)

Outputs:
  subtype_label_audit.parquet      one row per patient with all subtype calls joined
  subtype_marker_coverage.json     which marker genes are present/missing in RNA features
  subtype_counts.json              per-cohort x subtype-source sample counts
  audit_report.txt                 human-readable summary

Usage:
  python subtype_label_audit.py \\
      --cohort-parquet  data/model_ready_cohort.parquet \\
      --hoadley-tsv     data/TCGASubtype.20170308.tsv \\
      --out-dir         data/audit/
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List

import polars as pl
import requests

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
log = logging.getLogger("subtype_audit")

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

CBIOPORTAL_API = "https://www.cbioportal.org/api"

# cBioPortal pan_can_atlas_2018 studies carry per-marker-paper SUBTYPE columns
STUDIES: Dict[str, str] = {
    "BRCA": "brca_tcga_pan_can_atlas_2018",
    "LUAD": "luad_tcga_pan_can_atlas_2018",
    "PRAD": "prad_tcga_pan_can_atlas_2018",
}

# ---------------------------------------------------------------------------
# Subtype-marker gene panels (gene symbol -> Ensembl ID, no version)
#
# These are *subsets* of the canonical panels chosen for high-confidence
# mapping. Verify against your annotation version (GENCODE/Ensembl release)
# before treating absence as a real gap. Extend per your needs.
#
# BRCA:  PAM50 subset (Parker et al. 2009; TCGA Nature 2012)
# LUAD:  TRU/PP/PI markers (Wilkerson 2012; TCGA Nature 2014)
# PRAD:  iCluster markers (TCGA Cell 2015) - mostly mutation-driven, but
#        AR / ERG / FOXA1 expression also stratifies subtypes.
# ---------------------------------------------------------------------------

PAM50_BRCA: Dict[str, str] = {
    "ESR1":  "ENSG00000091831",
    "PGR":   "ENSG00000082175",
    "ERBB2": "ENSG00000141736",
    "MKI67": "ENSG00000148773",
    "FOXA1": "ENSG00000129514",
    "FOXC1": "ENSG00000054598",
    "BCL2":  "ENSG00000171791",
    "CDH1":  "ENSG00000039068",
    "MYC":   "ENSG00000136997",
    "BIRC5": "ENSG00000089685",
    "AURKA": "ENSG00000087586",
    "CCNB1": "ENSG00000134057",
    "MYBL2": "ENSG00000101057",
    "KRT5":  "ENSG00000186081",
    "KRT14": "ENSG00000186847",
    "KRT17": "ENSG00000128422",
    "EGFR":  "ENSG00000146648",
    "FGFR4": "ENSG00000160867",
    "CDC20": "ENSG00000117399",
    "MELK":  "ENSG00000165304",
}

LUAD_MARKERS: Dict[str, str] = {
    "NKX2-1": "ENSG00000136352",
    "SFTPB":  "ENSG00000168878",
    "SFTPC":  "ENSG00000168484",
    "TP63":   "ENSG00000073282",
    "KRT5":   "ENSG00000186081",
    "KRT6A":  "ENSG00000205420",
    "MUC1":   "ENSG00000185499",
}

PRAD_MARKERS: Dict[str, str] = {
    "AR":    "ENSG00000169083",
    "ERG":   "ENSG00000157554",
    "ETV1":  "ENSG00000006468",
    "ETV4":  "ENSG00000175832",
    "FLI1":  "ENSG00000151702",
    "SPOP":  "ENSG00000121067",
    "FOXA1": "ENSG00000129514",
    "PTEN":  "ENSG00000171862",
    "TP53":  "ENSG00000141510",
}

PANELS: Dict[str, Dict[str, str]] = {
    "BRCA": PAM50_BRCA,
    "LUAD": LUAD_MARKERS,
    "PRAD": PRAD_MARKERS,
}

# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_model_ready(parquet_path: Path) -> pl.DataFrame:
    """Read just patient_id + cohort + RNA feature columns from the cohort parquet."""
    log.info("Loading model-ready cohort: %s", parquet_path)
    df = pl.read_parquet(parquet_path)
    log.info("  shape=%s  cohorts=%s", df.shape, df["cohort"].unique().to_list())
    if df["patient_id"].n_unique() != df.height:
        log.warning("  patient_id is NOT unique (%d unique / %d rows) - dedupe upstream",
                    df["patient_id"].n_unique(), df.height)
    return df


def load_hoadley(tsv_path: Path) -> pl.DataFrame:
    """Load TCGASubtype.20170308.tsv and normalize to a 12-char patient_id key.

    The file has roughly these columns:
      sampleID  cancer.type  Subtype_mRNA  Subtype_DNAmeth  Subtype_protein
      Subtype_miRNA  Subtype_CNA  Subtype_Integrative  Subtype_other  Subtype_Selected

    Subtype_Selected is the canonical per-cancer call. We keep the rest for
    sensitivity analyses (e.g., does mRNA-based and methylation-based agree?).
    """
    log.info("Loading Hoadley subtypes: %s", tsv_path)
    df = pl.read_csv(
        tsv_path,
        separator="\t",
        null_values=["", "NA", "Not_Applicable", "[Not Applicable]", "NotApplicable"],
        infer_schema_length=10000,
    )
    log.info("  rows=%d  columns=%s", df.height, df.columns)

    # Locate sample-ID column (header naming has drifted across releases)
    id_col = next(
        (c for c in df.columns if c.lower() in ("sampleid", "pan.samplesid", "sample_id")),
        None,
    )
    if id_col is None:
        raise ValueError(f"No sample ID column in {tsv_path}; columns={df.columns}")

    # Rename Subtype_* cols with hoadley_ prefix to keep them distinct downstream
    rename = {
        c: f"hoadley_{c.lower()}"
        for c in df.columns
        if c.lower().startswith("subtype")
    }
    df = df.rename(rename)

    df = df.with_columns(
        pl.col(id_col).str.slice(0, 12).alias("patient_id"),
    )

    # Hoadley file is at sample/aliquot level. Dedupe to patient.
    # If a patient has multiple aliquots, keep the first (Subtype_Selected
    # is generally consistent across aliquots of the same primary tumor).
    before = df.height
    df = df.unique(subset=["patient_id"], keep="first")
    log.info("  patients (post-dedupe): %d (dropped %d duplicates)", df.height, before - df.height)

    keep_cols = ["patient_id"] + [c for c in df.columns if c.startswith("hoadley_")]
    return df.select(keep_cols)


def load_cbioportal_subtypes(cohort: str, cache_dir: Path) -> pl.DataFrame:
    """Pull per-cancer SUBTYPE from the cBioPortal REST API (patient-level).

    The GitHub datahub raw URLs return Git LFS pointers, so we use the public
    cBioPortal API: /studies/{studyId}/clinical-data with
    clinicalDataType=PATIENT and attributeId=SUBTYPE.
    """
    study = STUDIES[cohort]
    cache_path = cache_dir / f"{study}__subtype_api.json"

    if not cache_path.exists():
        cache_dir.mkdir(parents=True, exist_ok=True)
        url = (
            f"{CBIOPORTAL_API}/studies/{study}/clinical-data"
            f"?clinicalDataType=PATIENT&attributeId=SUBTYPE&pageSize=10000"
        )
        log.info("Downloading %s", url)
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        cache_path.write_bytes(r.content)
    else:
        log.info("Using cached %s", cache_path)

    records = json.loads(cache_path.read_text())
    out_col = f"{cohort.lower()}_marker_subtype"

    if not records:
        log.warning("  %s: no SUBTYPE records returned", cohort)
        return pl.DataFrame({"patient_id": pl.Series([], dtype=pl.Utf8),
                             out_col: pl.Series([], dtype=pl.Utf8)})

    df = pl.DataFrame({
        "patient_id": [rec["patientId"] for rec in records],
        out_col:      [rec["value"]     for rec in records],
    }).unique(subset=["patient_id"], keep="first")

    n_labeled = df.filter(pl.col(out_col).is_not_null()).height
    log.info("  %s: %d patients, %d with non-null SUBTYPE", cohort, df.height, n_labeled)
    return df


# ---------------------------------------------------------------------------
# Audit logic
# ---------------------------------------------------------------------------

def strip_ensembl_version(feat: str) -> str:
    """ENSG00000141736.14 -> ENSG00000141736"""
    return feat.split(".", 1)[0]


def check_marker_coverage(
    rna_features: List[str],
    panels: Dict[str, Dict[str, str]],
) -> Dict:
    """For each cohort's marker panel, report which genes are in RNA features."""
    rna_base = {strip_ensembl_version(f) for f in rna_features if f.startswith("ENSG")}
    out: Dict = {}
    for cohort, panel in panels.items():
        present = {sym: ensg for sym, ensg in panel.items() if ensg in rna_base}
        missing = {sym: ensg for sym, ensg in panel.items() if ensg not in rna_base}
        out[cohort] = {
            "panel_size": len(panel),
            "n_present": len(present),
            "n_missing": len(missing),
            "present": present,
            "missing": missing,
        }
        log.info("  %s panel: %d/%d markers present", cohort, len(present), len(panel))
    return out


def join_subtype_calls(
    cohort_df: pl.DataFrame,
    hoadley_df: pl.DataFrame,
    cbio_dfs: Dict[str, pl.DataFrame],
) -> pl.DataFrame:
    """Build the per-patient audit table: id + cohort + all subtype calls."""
    base = cohort_df.select(["patient_id", "cohort"])

    out = base.join(hoadley_df, on="patient_id", how="left")
    for _, cbio_df in cbio_dfs.items():
        out = out.join(cbio_df, on="patient_id", how="left")

    # Coalesce per-cohort marker subtypes into a single column for convenience
    out = out.with_columns(
        pl.when(pl.col("cohort") == "BRCA").then(pl.col("brca_marker_subtype"))
        .when(pl.col("cohort") == "LUAD").then(pl.col("luad_marker_subtype"))
        .when(pl.col("cohort") == "PRAD").then(pl.col("prad_marker_subtype"))
        .otherwise(None)
        .alias("marker_subtype")
    )
    return out


def per_cohort_subtype_counts(audit: pl.DataFrame) -> Dict:
    """Sample counts per cohort x subtype source."""
    summary: Dict = {}
    subtype_cols = [c for c in audit.columns if "subtype" in c.lower() and c != "cohort"]

    for cohort in audit["cohort"].unique().to_list():
        cohort_df = audit.filter(pl.col("cohort") == cohort)
        n_total = cohort_df.height
        cohort_summary = {"n_patients": n_total, "by_source": {}}

        for col in subtype_cols:
            n_labeled = cohort_df.filter(pl.col(col).is_not_null()).height
            counts = (
                cohort_df.group_by(col)
                         .agg(pl.len().alias("n"))
                         .sort("n", descending=True)
            )
            cohort_summary["by_source"][col] = {
                "n_labeled": n_labeled,
                "n_unlabeled": n_total - n_labeled,
                "counts": {
                    (k if k is not None else "__null__"): v
                    for k, v in zip(counts[col].to_list(), counts["n"].to_list())
                },
            }
        summary[cohort] = cohort_summary
    return summary


def write_report(
    path: Path,
    audit: pl.DataFrame,
    coverage: Dict,
    counts: Dict,
) -> None:
    lines: List[str] = []
    lines.append("=" * 72)
    lines.append("TCGA SUBTYPE LABEL AUDIT")
    lines.append("=" * 72)
    lines.append(f"Total patients: {audit.height}")
    lines.append("")

    lines.append("--- SUBTYPE-MARKER GENE COVERAGE IN RNA FEATURES ---")
    for cohort, info in coverage.items():
        lines.append(f"\n{cohort}: {info['n_present']}/{info['panel_size']} markers present")
        if info["missing"]:
            lines.append("  MISSING from RNA features:")
            for sym, ensg in info["missing"].items():
                lines.append(f"    {sym} ({ensg})")
    lines.append("")

    lines.append("--- PER-COHORT SUBTYPE COUNTS ---")
    for cohort, info in counts.items():
        lines.append(f"\n{cohort} (n={info['n_patients']}):")
        for source, src_info in info["by_source"].items():
            lines.append(
                f"  [{source}]  labeled={src_info['n_labeled']}  "
                f"unlabeled={src_info['n_unlabeled']}"
            )
            for sub, n in src_info["counts"].items():
                if sub != "__null__":
                    lines.append(f"      {sub}: {n}")

    path.write_text("\n".join(lines))
    log.info("Wrote report: %s", path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: List[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--cohort-parquet", type=Path, required=True)
    p.add_argument("--hoadley-tsv", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--cbio-cache", type=Path, default=Path("data/cbio_cache"))
    args = p.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    cohort_df = load_model_ready(args.cohort_parquet)
    rna_features = [c for c in cohort_df.columns if c.startswith("ENSG")]
    log.info("RNA features in cohort: %d", len(rna_features))

    hoadley_df = load_hoadley(args.hoadley_tsv)

    cbio_dfs = {
        cohort: load_cbioportal_subtypes(cohort, args.cbio_cache)
        for cohort in STUDIES
    }

    coverage = check_marker_coverage(rna_features, PANELS)
    audit = join_subtype_calls(cohort_df, hoadley_df, cbio_dfs)
    counts = per_cohort_subtype_counts(audit)

    audit.write_parquet(args.out_dir / "subtype_label_audit.parquet")
    (args.out_dir / "subtype_marker_coverage.json").write_text(json.dumps(coverage, indent=2))
    (args.out_dir / "subtype_counts.json").write_text(json.dumps(counts, indent=2))
    write_report(args.out_dir / "audit_report.txt", audit, coverage, counts)

    log.info("Done. Outputs in %s", args.out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
