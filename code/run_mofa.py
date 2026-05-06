"""Run MOFA+ multi-omic factor model on a single cohort.

Loads data/merged_all_cohorts.parquet, filters to the requested cohort, splits
features into RNA / methylation / CNV views by name pattern, applies
modality-appropriate preprocessing, and fits MOFA+ via mofapy2.

Outputs (under --out-dir):
  - mofa_model.hdf5                     full model (samples + loadings + factors)
  - factor_scores.csv                   patient_id × factor_<i> (latent Z)
  - variance_explained.csv              factor × modality (R² per view)
  - top_loadings_top25_<view>.csv       per view, factor × top-25 features
  - run_manifest.json                   args + per-modality shapes + run timing

Usage:
  python run_mofa.py --cohort BRCA
  python run_mofa.py --cohort BRCA --n-factors 20 --max-iter 1500
  python run_mofa.py --cohort BRCA --max-iter 50 --convergence fast   # smoke

Preprocessing rationale (recorded here so the writeup can cite it):
  RNA:  values are FPKM_UQ (right-skewed, range 0..hundreds).
        log2(x + 1) compresses the tail and makes the distribution
        approximately gaussian per gene.
  meth: values are beta (0..1, bimodal). M-value = log2(β/(1-β)) is the
        standard transform — variance-stabilizing and well-modeled by
        a gaussian likelihood.
  CNV:  values are already segment-mean log2(tumor/normal) ratios, mean
        approximately 0 and bounded ±2.5. No transform; MOFA+ scales views
        internally so the small CNV view doesn't dominate or get drowned.

NaNs are preserved through preprocessing. mofapy2 handles missing values
natively — they're treated as latent in the likelihood.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time
from typing import Callable, Iterable

import h5py
import numpy as np
import pandas as pd
import polars as pl
from mofapy2.run.entry_point import entry_point

ID_COL = "patient_id"
COHORT_COL = "cohort"
COHORT_CODE_COL = "cohort_code"
META_COLS = {ID_COL, COHORT_COL, COHORT_CODE_COL}
COHORTS = ["BRCA", "LUAD", "PRAD"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cohort", choices=COHORTS, required=True)
    p.add_argument("--input", type=pathlib.Path,
                   default=pathlib.Path("data/TCGA-BRCA/merged_brca_6view.parquet"),
                   help="Wide merged parquet. Default is the BRCA 6-view file from "
                        "Phase 4d. Pass merged_all_cohorts.parquet for the legacy 3-view path.")
    p.add_argument("--out-dir", type=pathlib.Path, default=None,
                   help="Default: data/mofa_<COHORT>")
    p.add_argument("--n-factors", type=int, default=15,
                   help="Upper bound; ARD prunes inactive factors during training.")
    p.add_argument("--max-iter", type=int, default=1000)
    p.add_argument("--convergence", choices=["fast", "medium", "slow"], default="medium")
    p.add_argument("--max-missing-rate", type=float, default=None,
                   help="Drop features missing in >this fraction of patients. "
                        "When unset, uses per-view defaults from MAX_MISSING_RATES "
                        "(methylation=0.40, others=0.20). When set, overrides "
                        "uniformly for all views.")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def split_views(df: pd.DataFrame) -> dict[str, list[str]]:
    """Group non-meta columns into views by column-name prefix.

    Prefixes are emitted by src/merge.merge_brca_6view (Phase 4d):
        ENSG*  → RNA            cg*    → methylation
        RPPA_* → RPPA           MIR_*  → miRNA
        MUT_*  → mutations      else   → CNV (chr-arm names: 1p, 1q, ..., Xq)

    Returns a dict in fixed view order (mutations last so the Bernoulli view
    sits at the end of likelihoods/matrices for readability). Empty views
    are kept in the dict and skipped at training time — that way the legacy
    3-view merged_all_cohorts.parquet still works (RPPA/miRNA/mutations
    just come back empty).
    """
    views: dict[str, list[str]] = {
        "RNA": [], "methylation": [], "CNV": [],
        "RPPA": [], "miRNA": [], "mutations": [],
    }
    for c in df.columns:
        if c in META_COLS:
            continue
        if c.startswith("ENSG"):
            views["RNA"].append(c)
        elif c.startswith("cg"):
            views["methylation"].append(c)
        elif c.startswith("RPPA_"):
            views["RPPA"].append(c)
        elif c.startswith("MIR_"):
            views["miRNA"].append(c)
        elif c.startswith("MUT_"):
            views["mutations"].append(c)
        else:
            views["CNV"].append(c)
    return views


def filter_by_missingness(arr: np.ndarray, cols: list[str], max_rate: float) -> tuple[np.ndarray, list[str]]:
    """Drop columns whose missing rate exceeds max_rate."""
    miss = np.isnan(arr).mean(axis=0)
    keep_idx = np.where(miss <= max_rate)[0]
    kept_cols = [cols[i] for i in keep_idx]
    return arr[:, keep_idx], kept_cols


def preprocess_rna(X: np.ndarray) -> np.ndarray:
    """log2(x + 1) for FPKM_UQ. NaNs preserved."""
    return np.log2(X + 1.0)


def preprocess_methylation(X: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    """beta -> M-value with edge clipping. NaNs preserved."""
    Xc = np.clip(X, eps, 1.0 - eps)
    return np.log2(Xc / (1.0 - Xc))


def preprocess_cnv(X: np.ndarray) -> np.ndarray:
    """CNV is already log2-ratio. No transform; mofapy2 scale_views handles scale."""
    return X.astype(float)


def preprocess_rppa(X: np.ndarray) -> np.ndarray:
    """RPPA values are already lab-normalized (median-centered, log2-ratio).
    No transform; scale_views handles cross-view scale parity."""
    return X.astype(float)


def preprocess_mirna(X: np.ndarray) -> np.ndarray:
    """miRNA values are RPM (reads per million, library-normalized at ingest).
    log2(x + 1) compresses the heavy right tail (many low-expression miRNAs,
    a few very high) and stabilizes variance for the Gaussian likelihood.
    NaNs preserved (per-modality parquet is 0% missing today, but keep the
    same defensive handling as RNA)."""
    return np.log2(X + 1.0)


def preprocess_mutations(X: np.ndarray) -> np.ndarray:
    """Mutations are 0/1 indicators by construction (merge fills absences with 0,
    pivots non-silent variants to 1). Just cast to float so mofapy2's Bernoulli
    likelihood gets a clean dtype."""
    return X.astype(float)


PREPROCESSORS: dict[str, "Callable[[np.ndarray], np.ndarray]"] = {
    "RNA": preprocess_rna,
    "methylation": preprocess_methylation,
    "CNV": preprocess_cnv,
    "RPPA": preprocess_rppa,
    "miRNA": preprocess_mirna,
    "mutations": preprocess_mutations,
}

LIKELIHOODS: dict[str, str] = {
    "RNA": "gaussian",
    "methylation": "gaussian",
    "CNV": "gaussian",
    "RPPA": "gaussian",
    "miRNA": "gaussian",
    "mutations": "bernoulli",
}

# Per-view per-column missingness thresholds. A column is dropped if its NaN
# rate within the view exceeds this. Methylation needs a higher cap because
# the merged 6-view's mean meth missingness is ~27% (vs ~18% on the source
# parquet — the 6-way inner join shifts the kept-patient population toward
# worse meth coverage). At 0.20 only ~8% of meth probes survive; at 0.40
# we keep ~half. Other views are well below 20% missing so they pass cleanly.
# --max-missing-rate on the CLI overrides this dict uniformly when set.
MAX_MISSING_RATES: dict[str, float] = {
    "RNA": 0.20,
    "methylation": 0.40,
    "CNV": 0.20,
    "RPPA": 0.20,
    "miRNA": 0.20,
    "mutations": 0.20,
}


def extract_factor_scores(model_path: pathlib.Path, sample_names: list[str]) -> pd.DataFrame:
    """Read the latent Z matrix from the trained HDF5 model."""
    with h5py.File(model_path, "r") as f:
        # mofapy2 writes expectations under /expectations/Z/<group>
        groups = list(f["expectations/Z"].keys())
        Z_blocks = [np.array(f[f"expectations/Z/{g}"]) for g in groups]
        Z = np.concatenate(Z_blocks, axis=1) if Z_blocks[0].ndim == 2 else None
        if Z is None:
            raise RuntimeError("Unexpected Z shape in HDF5")
    n_factors = Z.shape[0]
    df = pd.DataFrame(Z.T, index=sample_names, columns=[f"factor_{i+1}" for i in range(n_factors)])
    df.index.name = ID_COL
    return df


def extract_variance_explained(model_path: pathlib.Path, view_names: list[str]) -> pd.DataFrame:
    """R² per factor per view, from the trained model. Stored as (views, factors)
    in HDF5; we return (factors, views) for downstream readability."""
    with h5py.File(model_path, "r") as f:
        groups = list(f["variance_explained/r2_per_factor"].keys())
        per_group = [np.array(f[f"variance_explained/r2_per_factor/{g}"]) for g in groups]
    r2 = per_group[0]  # one group per cohort, shape = (n_views, n_factors)
    n_views, n_factors = r2.shape
    df = pd.DataFrame(
        r2.T,
        index=[f"factor_{i+1}" for i in range(n_factors)],
        columns=view_names[:n_views],
    )
    return df


def extract_top_loadings(model_path: pathlib.Path, view_idx: int, view_name: str,
                          feature_names: list[str], top_k: int = 25) -> pd.DataFrame:
    """For one view: top-K features per factor by |loading|.
    HDF5 W is stored as (factors, features) for each view."""
    with h5py.File(model_path, "r") as f:
        W = np.array(f[f"expectations/W/{view_name}"])
    if W.shape[1] != len(feature_names):
        W = W.T
    n_factors, n_features = W.shape
    rows = []
    for k in range(n_factors):
        row = W[k, :]
        order = np.argsort(-np.abs(row))[:top_k]
        for rank, idx in enumerate(order, 1):
            rows.append({
                "factor": f"factor_{k+1}",
                "rank": rank,
                "feature": feature_names[idx],
                "loading": float(row[idx]),
            })
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir or pathlib.Path(f"data/mofa_{args.cohort}")
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] {args.input}", flush=True)
    # Lazy scan + filter pushdown: defensive against --input being pointed at a
    # long-format modality parquet (methylation = 444M rows, OOM-kills a 7.6 GB
    # box if read eagerly). For the wide merged_* parquet this is essentially
    # free — the filter prunes nothing — but the convention matters.
    lf = pl.scan_parquet(args.input)
    schema_names = lf.collect_schema().names()
    if COHORT_COL in schema_names:
        lf = lf.filter(pl.col(COHORT_COL) == args.cohort)
        print(f"[load] filtering to cohort={args.cohort}", flush=True)
    else:
        # The Phase-4d BRCA 6-view file has no cohort column — every row is
        # BRCA by construction. --cohort is still required because it labels
        # the output dir, run manifest, and mofapy2 groups_names.
        print(f"[load] no '{COHORT_COL}' column in input — single-cohort file, "
              f"using all rows (--cohort={args.cohort} used for output labeling)",
              flush=True)
    df = lf.collect(engine="streaming").to_pandas()
    print(f"[load] rows={len(df)}", flush=True)

    views = split_views(df)
    sample_names = df[ID_COL].astype(str).tolist()

    # Effective per-view missingness threshold: CLI override (uniform) wins
    # over per-view defaults from MAX_MISSING_RATES.
    if args.max_missing_rate is not None:
        eff_rates = {v: args.max_missing_rate for v in MAX_MISSING_RATES}
        print(f"[filter] --max-missing-rate={args.max_missing_rate:.0%} "
              f"overrides per-view defaults uniformly", flush=True)
    else:
        eff_rates = dict(MAX_MISSING_RATES)

    matrices: list[np.ndarray] = []
    feature_lists: list[list[str]] = []
    pre_shapes: dict[str, dict] = {}

    for view_name, cols in views.items():
        if not cols:
            print(f"[skip] view {view_name}: no columns", flush=True)
            continue
        rate = eff_rates[view_name]
        X_raw = df[cols].to_numpy(dtype=float)
        X_kept, kept_cols = filter_by_missingness(X_raw, cols, rate)
        X_pre = PREPROCESSORS[view_name](X_kept)
        matrices.append(X_pre)
        feature_lists.append(kept_cols)
        pre_shapes[view_name] = {
            "n_features_in": len(cols),
            "n_features_kept": len(kept_cols),
            "n_dropped_by_missing": len(cols) - len(kept_cols),
            "likelihood": LIKELIHOODS[view_name],
            "max_missing_rate": rate,
        }
        print(f"[view] {view_name} ({LIKELIHOODS[view_name]}): kept "
              f"{len(kept_cols)}/{len(cols)} features after missingness filter "
              f"({rate:.0%})", flush=True)

    view_names = [v for v in views if views[v]]

    print(f"[mofa] building entry point: {args.n_factors} factors, "
          f"{args.max_iter} max iter, convergence={args.convergence}", flush=True)
    ent = entry_point()
    ent.set_data_options(scale_views=True)
    ent.set_data_matrix(
        data=[[m] for m in matrices],
        likelihoods=[LIKELIHOODS[v] for v in view_names],
        views_names=view_names,
        groups_names=[args.cohort],
        samples_names=[sample_names],
        features_names=feature_lists,
    )
    ent.set_model_options(
        factors=args.n_factors,
        ard_factors=True,
        ard_weights=True,
        spikeslab_weights=False,
    )
    ent.set_train_options(
        iter=args.max_iter,
        convergence_mode=args.convergence,
        dropR2=0.001,
        gpu_mode=False,
        verbose=False,
        seed=args.seed,
    )
    ent.build()

    t0 = time.time()
    ent.run()
    elapsed = time.time() - t0
    print(f"[mofa] training done in {elapsed:.1f}s", flush=True)

    model_path = out_dir / "mofa_model.hdf5"
    ent.save(str(model_path), save_data=True)
    print(f"[save] {model_path}", flush=True)

    factor_df = extract_factor_scores(model_path, sample_names)
    factor_df.to_csv(out_dir / "factor_scores.csv")

    var_df = extract_variance_explained(model_path, view_names)
    var_df.to_csv(out_dir / "variance_explained.csv")

    for i, vn in enumerate(view_names):
        loadings_df = extract_top_loadings(model_path, i, vn, feature_lists[i])
        loadings_df.to_csv(out_dir / f"top_loadings_top25_{vn}.csv", index=False)

    manifest = {
        "cohort": args.cohort,
        "input": str(args.input),
        "out_dir": str(out_dir),
        "n_samples": len(sample_names),
        "n_factors_requested": args.n_factors,
        "n_factors_active": int(var_df.shape[0]),
        "max_iter": args.max_iter,
        "convergence": args.convergence,
        "max_missing_rate_cli": args.max_missing_rate,
        "max_missing_rates_per_view": eff_rates,
        "seed": args.seed,
        "training_seconds": round(elapsed, 1),
        "views": pre_shapes,
        "variance_explained_top": {
            f"factor_{i+1}": {vn: float(var_df.iloc[i, j]) for j, vn in enumerate(view_names)}
            for i in range(min(5, var_df.shape[0]))
        },
    }
    (out_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"[done] {args.cohort}: {var_df.shape[0]} active factors", flush=True)
    # Note: values are already R² in percent (mofapy2 convention), do not re-scale.
    print("Variance explained (first 8 factors, %):")
    print(var_df.head(8).round(2).to_string())


if __name__ == "__main__":
    main()
