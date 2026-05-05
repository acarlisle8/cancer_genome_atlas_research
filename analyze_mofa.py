"""Post-training analysis of a MOFA+ model.

Reads the HDF5 model + CSVs produced by run_mofa.py and produces:
  - variance_explained_heatmap.png   factors × modalities, R² heatmap
  - kmeans_silhouette.csv            silhouette score for k=2..10
  - cluster_assignments.csv          patient_id × cluster (best k by silhouette)
  - cluster_vs_subtype.json          ARI / NMI vs published-subtype labels
                                     (only if --known-labels provided)

Usage:
  python analyze_mofa.py --model-dir data/mofa_BRCA \\
      --known-labels data/audit/subtype_label_audit.parquet \\
      --label-col hoadley_subtype_selected \\
      --cancer-type BRCA
"""
from __future__ import annotations

import argparse
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model-dir", type=pathlib.Path, required=True,
                   help="Directory produced by run_mofa.py (contains factor_scores.csv etc.)")
    p.add_argument("--k-min", type=int, default=2)
    p.add_argument("--k-max", type=int, default=10)
    p.add_argument("--known-labels", type=pathlib.Path, default=None,
                   help="Parquet with patient_id + label column (subtype audit table).")
    p.add_argument("--label-col", type=str, default="hoadley_subtype_selected")
    p.add_argument("--cancer-type", type=str, default=None,
                   help="Filter labels to rows starting with this prefix (e.g. BRCA).")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def variance_heatmap(var_df: pd.DataFrame, out_path: pathlib.Path) -> None:
    """Heatmap of factor × view R² (percent)."""
    pct = (var_df * 100).round(2)
    fig, ax = plt.subplots(figsize=(max(4, 0.5 * pct.shape[1] + 2), max(5, 0.35 * pct.shape[0] + 1)))
    im = ax.imshow(pct.values, aspect="auto", cmap="viridis")
    ax.set_xticks(range(pct.shape[1]))
    ax.set_xticklabels(pct.columns)
    ax.set_yticks(range(pct.shape[0]))
    ax.set_yticklabels(pct.index)
    for i in range(pct.shape[0]):
        for j in range(pct.shape[1]):
            v = pct.values[i, j]
            ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                    color="white" if v < pct.values.max() * 0.6 else "black", fontsize=8)
    ax.set_title("Variance explained per factor per view (%)")
    fig.colorbar(im, ax=ax, label="R² (%)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def silhouette_sweep(Z: np.ndarray, k_min: int, k_max: int, seed: int) -> pd.DataFrame:
    rows = []
    for k in range(k_min, k_max + 1):
        km = KMeans(n_clusters=k, random_state=seed, n_init=10).fit(Z)
        sil = silhouette_score(Z, km.labels_)
        rows.append({"k": k, "silhouette": float(sil), "inertia": float(km.inertia_)})
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    md = args.model_dir
    if not md.exists():
        raise SystemExit(f"model dir not found: {md}")

    factor_df = pd.read_csv(md / "factor_scores.csv", index_col=0)
    var_df = pd.read_csv(md / "variance_explained.csv", index_col=0)
    print(f"[load] {factor_df.shape[0]} samples × {factor_df.shape[1]} factors")

    Z = factor_df.to_numpy()
    Z = np.nan_to_num(Z, nan=0.0)

    # 1) variance heatmap
    variance_heatmap(var_df, md / "variance_explained_heatmap.png")
    print(f"[plot] {md / 'variance_explained_heatmap.png'}")

    # 2) silhouette sweep on factor scores
    sil_df = silhouette_sweep(Z, args.k_min, args.k_max, args.seed)
    sil_df.to_csv(md / "kmeans_silhouette.csv", index=False)
    best_k = int(sil_df.loc[sil_df["silhouette"].idxmax(), "k"])
    best_sil = float(sil_df["silhouette"].max())
    print(f"[silhouette] best k={best_k} (silhouette={best_sil:.3f})")

    # 3) cluster assignments at best k
    km = KMeans(n_clusters=best_k, random_state=args.seed, n_init=10).fit(Z)
    cluster_df = pd.DataFrame({"patient_id": factor_df.index, "cluster": km.labels_})
    cluster_df.to_csv(md / "cluster_assignments.csv", index=False)
    print(f"[clusters] sizes: {dict(zip(*np.unique(km.labels_, return_counts=True)))}")

    # 4) compare to known subtype labels if provided
    summary = {
        "best_k": best_k,
        "best_silhouette": best_sil,
        "cluster_sizes": {int(k): int(v) for k, v in zip(*np.unique(km.labels_, return_counts=True))},
        "agreement_vs_known": None,
    }
    if args.known_labels is not None:
        labels_df = pl.read_parquet(args.known_labels).select(["patient_id", args.label_col]).to_pandas()
        if args.cancer_type:
            labels_df = labels_df[labels_df[args.label_col].astype(str).str.startswith(args.cancer_type)]
        merged = cluster_df.merge(labels_df, on="patient_id", how="inner").dropna(subset=[args.label_col])
        if len(merged) > 0:
            ari = adjusted_rand_score(merged[args.label_col], merged["cluster"])
            nmi = normalized_mutual_info_score(merged[args.label_col], merged["cluster"])
            summary["agreement_vs_known"] = {
                "n_labeled": int(len(merged)),
                "label_col": args.label_col,
                "ari": round(float(ari), 4),
                "nmi": round(float(nmi), 4),
            }
            print(f"[agreement] vs {args.label_col} (n={len(merged)}): ARI={ari:.3f}, NMI={nmi:.3f}")
        else:
            summary["agreement_vs_known"] = {"error": "no overlap between clusters and labels"}

    (md / "cluster_vs_subtype.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"[done] summary -> {md / 'cluster_vs_subtype.json'}")


if __name__ == "__main__":
    main()
