"""Generate report figures for the Phase 4e/4f BRCA MOFA+ analysis.

Reads artifacts produced by run_mofa.py + analyze_mofa.py from data/mofa_BRCA/
plus the PAM50 labels from the audit parquet, writes PNG figures to
data/mofa_BRCA/figures/.

Usage:
    uv run python make_report_figures.py

Optional inputs:
- /tmp/mofa_full.log — if present, an ELBO convergence figure is produced.
  (Capture this on the next run with: `... | tee /tmp/mofa_full.log`)
"""
from __future__ import annotations

import pathlib
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl

MODEL_DIR = pathlib.Path("data/mofa_BRCA")
FIGURES_DIR = MODEL_DIR / "figures"
LABELS_PATH = pathlib.Path("data/audit/subtype_label_audit.parquet")
LABEL_COL = "hoadley_subtype_selected"
LOG_PATH = pathlib.Path("/tmp/mofa_full.log")

VIEWS = ["RNA", "methylation", "CNV", "RPPA", "miRNA", "mutations"]

PAM50_ORDER = ["BRCA.Basal", "BRCA.Her2", "BRCA.LumA", "BRCA.LumB", "BRCA.Normal"]
PAM50_COLORS = {
    "BRCA.Basal":  "#d62728",  # red
    "BRCA.Her2":   "#ff7f0e",  # orange
    "BRCA.LumA":   "#1f77b4",  # blue
    "BRCA.LumB":   "#9467bd",  # purple
    "BRCA.Normal": "#2ca02c",  # green
}


def fig_variance_heatmap(out_path: pathlib.Path) -> None:
    """Variance explained per factor per view, in percent.

    mofapy2 stores R² values already in percent in variance_explained.csv —
    do NOT multiply by 100 again.
    """
    var_df = pd.read_csv(MODEL_DIR / "variance_explained.csv", index_col=0)
    pct = var_df.round(2)
    fig, ax = plt.subplots(figsize=(7.5, max(5, 0.4 * pct.shape[0] + 1)))
    im = ax.imshow(pct.values, aspect="auto", cmap="viridis")
    ax.set_xticks(range(pct.shape[1]))
    ax.set_xticklabels(pct.columns, rotation=15)
    ax.set_yticks(range(pct.shape[0]))
    ax.set_yticklabels(pct.index)
    vmax = pct.values.max()
    for i in range(pct.shape[0]):
        for j in range(pct.shape[1]):
            v = pct.values[i, j]
            ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                    color="white" if v < vmax * 0.5 else "black", fontsize=8)
    ax.set_title("MOFA+ variance explained per factor per view (%)")
    fig.colorbar(im, ax=ax, label="R² (%)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def fig_f1_f2_scatter(out_path: pathlib.Path) -> None:
    """Factor 1 vs Factor 2 scatter, colored by PAM50 subtype."""
    factor_df = pd.read_csv(MODEL_DIR / "factor_scores.csv", index_col=0)
    factor_df.index.name = "patient_id"
    labels_df = (
        pl.read_parquet(LABELS_PATH)
        .select(["patient_id", LABEL_COL])
        .to_pandas()
        .dropna(subset=[LABEL_COL])
    )
    merged = factor_df[["factor_1", "factor_2"]].reset_index().merge(
        labels_df, on="patient_id", how="inner"
    )

    fig, ax = plt.subplots(figsize=(7, 6))
    for label in PAM50_ORDER:
        sub = merged[merged[LABEL_COL] == label]
        if sub.empty:
            continue
        ax.scatter(sub["factor_1"], sub["factor_2"],
                   c=PAM50_COLORS[label],
                   label=f"{label.replace('BRCA.', '')} (n={len(sub)})",
                   alpha=0.65, s=22, edgecolor="white", linewidth=0.3)
    ax.axhline(0, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.axvline(0, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.set_xlabel("Factor 1  (basal-vs-luminal axis)")
    ax.set_ylabel("Factor 2  (proliferation / driver-mutation axis)")
    ax.set_title(f"MOFA+ Factor 1 vs Factor 2 by PAM50 subtype  (n={len(merged)})")
    ax.legend(title="PAM50", loc="best", framealpha=0.9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def fig_cluster_contingency(out_path: pathlib.Path) -> None:
    """Cluster × PAM50 contingency heatmap (counts annotated, colored by row %)."""
    clusters_df = pd.read_csv(MODEL_DIR / "cluster_assignments.csv")
    labels_df = (
        pl.read_parquet(LABELS_PATH)
        .select(["patient_id", LABEL_COL])
        .to_pandas()
        .dropna(subset=[LABEL_COL])
    )
    merged = clusters_df.merge(labels_df, on="patient_id", how="inner")
    ct = pd.crosstab(merged["cluster"], merged[LABEL_COL])
    cols_present = [c for c in PAM50_ORDER if c in ct.columns]
    ct = ct[cols_present]
    ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100

    n_clusters = ct.shape[0]
    fig, ax = plt.subplots(figsize=(7.5, max(3.5, 0.7 * n_clusters + 1.5)))
    im = ax.imshow(ct_pct.values, aspect="auto", cmap="Blues", vmin=0, vmax=100)
    ax.set_xticks(range(ct.shape[1]))
    ax.set_xticklabels([c.replace("BRCA.", "") for c in ct.columns])
    ax.set_yticks(range(n_clusters))
    ax.set_yticklabels([f"cluster {c}\n(n={ct.loc[c].sum()})" for c in ct.index])
    for i in range(n_clusters):
        for j in range(ct.shape[1]):
            n = ct.values[i, j]
            pct = ct_pct.values[i, j]
            ax.text(j, i, f"{n}\n({pct:.0f}%)", ha="center", va="center",
                    color="white" if pct > 50 else "black", fontsize=9)
    ax.set_title(f"Cluster × PAM50 contingency  (k={n_clusters} k-means on factor scores)")
    ax.set_xlabel("PAM50 subtype")
    fig.colorbar(im, ax=ax, label="% of cluster")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def fig_top_loadings(factor_label: str, top_k: int, out_path: pathlib.Path) -> None:
    """Top-K loadings (by |loading|) for a factor, one panel per view."""
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), squeeze=False)
    for ax, view in zip(axes.flatten(), VIEWS):
        df = pd.read_csv(MODEL_DIR / f"top_loadings_top25_{view}.csv")
        sub = df[df["factor"] == factor_label].copy()
        if sub.empty:
            ax.set_title(f"{view}  (no loadings — view absent)")
            ax.axis("off")
            continue
        sub["abs_loading"] = sub["loading"].abs()
        sub = sub.nlargest(top_k, "abs_loading").sort_values("loading")
        colors = ["#d62728" if v < 0 else "#1f77b4" for v in sub["loading"]]
        ax.barh(range(len(sub)), sub["loading"], color=colors)
        ax.set_yticks(range(len(sub)))
        labels = [s if len(s) <= 22 else s[:20] + ".." for s in sub["feature"]]
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_title(view, fontsize=11)
        ax.axvline(0, color="black", linewidth=0.5)
        ax.grid(True, axis="x", alpha=0.3)
    fig.suptitle(f"Top-{top_k} loadings on {factor_label}  (red = negative, blue = positive)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def fig_elbo_convergence(log_path: pathlib.Path, out_path: pathlib.Path) -> bool:
    """Parse iter / ELBO / deltaELBO / Factors lines from the mofapy2 stdout
    log; emit a 2-panel figure (ELBO + log-scale ΔELBO%). Returns False if
    the log isn't found or doesn't contain iteration lines."""
    if not log_path.exists():
        return False
    pat = re.compile(
        r"^Iteration (\d+):\s+time=([\d.]+),\s+ELBO=(-?[\d.]+),\s+"
        r"deltaELBO=[\d.]+\s+\(([\d.]+)%\),\s+Factors=(\d+)"
    )
    iters: list[int] = []
    elbos: list[float] = []
    deltas: list[float] = []
    factors: list[int] = []
    for line in log_path.read_text().splitlines():
        m = pat.match(line)
        if m:
            iters.append(int(m.group(1)))
            elbos.append(float(m.group(3)))
            deltas.append(float(m.group(4)))
            factors.append(int(m.group(5)))
    if not iters:
        return False

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ax1, ax2 = axes
    ax1.plot(iters, elbos, "o-", color="tab:blue", markersize=2.5, linewidth=1)
    ax1.set_xlabel("Iteration")
    ax1.set_ylabel("ELBO")
    ax1.set_title(f"ELBO over {len(iters)} iterations  (final = {elbos[-1]:,.0f})")
    ax1.grid(True, alpha=0.3)

    ax2.plot(iters, deltas, "o-", color="tab:orange", markersize=2.5, linewidth=1)
    ax2.set_yscale("log")
    ax2.set_xlabel("Iteration")
    ax2.set_ylabel("ΔELBO (%)  [log scale]")
    ax2.set_title("Convergence")
    ax2.grid(True, alpha=0.3, which="both")
    fig.suptitle("MOFA+ training convergence (full run, --convergence medium)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return True


def main() -> None:
    if not MODEL_DIR.exists():
        raise SystemExit(f"model dir not found: {MODEL_DIR} — run run_mofa.py first")
    if not LABELS_PATH.exists():
        raise SystemExit(f"labels parquet not found: {LABELS_PATH}")
    FIGURES_DIR.mkdir(exist_ok=True)

    fig_variance_heatmap(FIGURES_DIR / "01_variance_explained_heatmap.png")
    print(f"[fig 1] {FIGURES_DIR / '01_variance_explained_heatmap.png'}")

    fig_f1_f2_scatter(FIGURES_DIR / "02_factor1_vs_factor2_pam50.png")
    print(f"[fig 2] {FIGURES_DIR / '02_factor1_vs_factor2_pam50.png'}")

    fig_cluster_contingency(FIGURES_DIR / "03_cluster_pam50_contingency.png")
    print(f"[fig 3] {FIGURES_DIR / '03_cluster_pam50_contingency.png'}")

    fig_top_loadings("factor_1", top_k=10,
                     out_path=FIGURES_DIR / "04_top_loadings_factor1.png")
    print(f"[fig 4] {FIGURES_DIR / '04_top_loadings_factor1.png'}")

    fig_top_loadings("factor_2", top_k=10,
                     out_path=FIGURES_DIR / "05_top_loadings_factor2.png")
    print(f"[fig 5] {FIGURES_DIR / '05_top_loadings_factor2.png'}")

    if fig_elbo_convergence(LOG_PATH, FIGURES_DIR / "06_elbo_convergence.png"):
        print(f"[fig 6] {FIGURES_DIR / '06_elbo_convergence.png'}")
    else:
        print(f"[fig 6] skipped — {LOG_PATH} not found")

    print(f"\nDone. Figures in {FIGURES_DIR}/")


if __name__ == "__main__":
    main()
