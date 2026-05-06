"""Generate project-wide overview figures for the report.

Reads existing artifacts from data/{classification_pipeline*,mofa_BRCA}/ and
produces cross-phase summary figures into data/figures/. Companion to
make_report_figures.py (which produces Phase-4-specific figures into
data/mofa_BRCA/figures/).

Figures produced:
  01  BRCA modality coverage funnel (per-modality patient counts → 6-way intersection)
  02  Cross-task XGBoost accuracy comparison (cohort + 3 subtype tasks)
  03  SHAP-vs-published-panel hits per cohort
  04  Silhouette score over k for MOFA+ factor scores
  05  Phase 3 vs Phase 4 BRCA: supervised XGBoost vs unsupervised MOFA+ clustering

Usage:
    uv run python make_overview_figures.py
"""
from __future__ import annotations

import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DATA = pathlib.Path("data")
FIGURES_DIR = DATA / "figures"

# Per-modality BRCA patient counts. Sourced from session logs (verified
# via DuckDB earlier in the project) — re-running n_unique on the 444M-row
# methylation parquet via Polars OOMs the 7.6 GB box; see known-issues.md
# "Polars memory landmines". DuckDB streaming would work but these counts
# are stable and documented, so hardcoding is the right call here.
BRCA_MODALITY_COUNTS = {
    "RNA-seq":     1095,
    "methylation": 1097,
    "CNV":         1098,
    "RPPA":         881,
    "miRNA":       1079,
    "mutations":    967,
}
BRCA_INTERSECTION = 744  # post 6-way inner-join in merged_brca_6view


def fig_modality_coverage(out_path: pathlib.Path) -> None:
    """Bar chart of BRCA per-modality patient counts vs the 6-way intersection.

    Visualizes the data funnel: each modality has 800–1100 patients, but
    the 6-way intersection (patients with all six modalities present)
    drops to 744 — the actual MOFA+ training cohort.
    """
    modalities = list(BRCA_MODALITY_COUNTS.keys())
    counts = list(BRCA_MODALITY_COUNTS.values())

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(modalities, counts, color="#4c78a8", edgecolor="white")
    for bar, c in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, c + 15, f"{c}",
                ha="center", va="bottom", fontsize=10)

    ax.axhline(BRCA_INTERSECTION, color="#d62728", linestyle="--", linewidth=2,
               label=f"6-way intersection (n={BRCA_INTERSECTION})")
    ax.text(len(modalities) - 0.5, BRCA_INTERSECTION + 15,
            f"6-way intersection: {BRCA_INTERSECTION}",
            color="#d62728", ha="right", va="bottom", fontsize=10, fontweight="bold")

    ax.set_ylabel("Patients (n)")
    ax.set_title("TCGA-BRCA modality coverage funnel\n"
                 "Per-modality patient counts → multi-omic training cohort")
    ax.set_ylim(0, max(counts) * 1.15)
    ax.legend(loc="lower left")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def fig_classifier_accuracy(out_path: pathlib.Path) -> None:
    """Cross-task XGBoost accuracy comparison.

    Reads classification_metrics.json from each of the 4 Phase 3 runs
    (cohort + per-cohort subtype). Plots accuracy / balanced accuracy /
    macro F1 with std-error caps from 5-fold CV.
    """
    runs = [
        ("Cohort\n(BRCA/LUAD/PRAD)", "classification_pipeline", 3),
        ("BRCA subtype\n(PAM50, 5-class)", "classification_pipeline_brca", 5),
        ("PRAD subtype\n(4-class)", "classification_pipeline_prad", 4),
        ("LUAD subtype\n(6-class)", "classification_pipeline_luad", 6),
    ]
    metrics_to_plot = [("accuracy", "Accuracy"),
                       ("balanced_accuracy", "Balanced acc."),
                       ("macro_f1", "Macro F1")]

    rows = []
    for label, run_dir, n_classes in runs:
        with open(DATA / run_dir / "classification_metrics.json") as f:
            m = json.load(f)
        rows.append({
            "task": label,
            "n_classes": n_classes,
            **{
                k: m["cv_scores"][k]["mean"] for k, _ in metrics_to_plot
            },
            **{
                f"{k}_std": m["cv_scores"][k]["std"] for k, _ in metrics_to_plot
            },
        })

    fig, ax = plt.subplots(figsize=(10, 5.5))
    x = np.arange(len(runs))
    bar_width = 0.25
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    for i, (key, display) in enumerate(metrics_to_plot):
        means = [r[key] for r in rows]
        stds = [r[f"{key}_std"] for r in rows]
        bars = ax.bar(x + i * bar_width, means, bar_width, yerr=stds,
                      capsize=4, label=display, color=colors[i],
                      edgecolor="white")
        for bar, m in zip(bars, means):
            ax.text(bar.get_x() + bar.get_width() / 2, m + 0.015, f"{m:.3f}",
                    ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x + bar_width)
    ax.set_xticklabels([r["task"] for r in rows], fontsize=9)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("Score (5-fold CV mean ± std)")
    ax.set_title("XGBoost classifier performance across tasks\n"
                 "5-fold CV with per-fold variance feature selection + marker-panel force-include")
    ax.axhline(1.0, color="black", linewidth=0.5, linestyle=":", alpha=0.4)
    ax.legend(loc="lower left")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def fig_shap_panel_hits(out_path: pathlib.Path) -> None:
    """SHAP-vs-published-panel hit-rate per cohort.

    Reads shap_panel_comparison.json from each per-cohort subtype run.
    Shows hits in top-10 / 25 / 50 vs total panel size. BRCA / PAM50 has
    the strongest panel recovery; LUAD has near-zero (consistent with the
    LUAD model's lower accuracy).
    """
    cohorts = [
        ("BRCA / PAM50", "classification_pipeline_brca"),
        ("PRAD / TCGA-Cell-2015", "classification_pipeline_prad"),
        ("LUAD / Wilkerson", "classification_pipeline_luad"),
    ]

    rows = []
    for label, run_dir in cohorts:
        with open(DATA / run_dir / "shap_panel_comparison.json") as f:
            m = json.load(f)
        rows.append({
            "cohort": label,
            "panel_size": m["panel_size"],
            "top10": m["hits_top10"]["hits"],
            "top25": m["hits_top25"]["hits"],
            "top50": m["hits_top50"]["hits"],
        })

    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = np.arange(len(rows))
    bar_width = 0.25
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    cuts = [("top10", "top-10"), ("top25", "top-25"), ("top50", "top-50")]
    for i, (key, display) in enumerate(cuts):
        hits = [r[key] for r in rows]
        bars = ax.bar(x + i * bar_width, hits, bar_width,
                      label=f"{display} SHAP", color=colors[i], edgecolor="white")
        for bar, h, total in zip(bars, hits, [r["panel_size"] for r in rows]):
            pct = 100 * h / total if total else 0
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.3,
                    f"{h}/{total}\n({pct:.0f}%)",
                    ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x + bar_width)
    ax.set_xticklabels([f"{r['cohort']}\n(panel n={r['panel_size']})" for r in rows],
                       fontsize=9)
    ax.set_ylabel("Panel-gene hits in top-K SHAP features")
    ax.set_title("SHAP top features vs published marker panels\n"
                 "BRCA/PAM50 strongly recovers the panel; LUAD does not "
                 "(consistent with low LUAD subtype accuracy)")
    ax.set_ylim(0, max([r["panel_size"] for r in rows]) * 1.2)
    ax.legend(loc="upper right")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def fig_silhouette_curve(out_path: pathlib.Path) -> None:
    """Silhouette score over k for k-means on MOFA+ factor scores.

    Reads the k=2..10 sweep saved during this session at
    kmeans_silhouette_sweep_k2-10.csv. Best k is 2 (silhouette 0.194), with
    k=3..10 all in 0.11–0.13 — confirms the data has weak globular
    structure beyond the basal-vs-rest split.
    """
    df = pd.read_csv(DATA / "mofa_BRCA" / "kmeans_silhouette_sweep_k2-10.csv")

    fig, ax1 = plt.subplots(figsize=(8.5, 5))
    ax1.plot(df["k"], df["silhouette"], "o-", color="#1f77b4",
             markersize=8, linewidth=2, label="silhouette")
    ax1.set_xlabel("k (number of clusters)")
    ax1.set_ylabel("Silhouette score", color="#1f77b4")
    ax1.tick_params(axis="y", labelcolor="#1f77b4")
    ax1.set_xticks(df["k"])
    ax1.grid(True, alpha=0.3)

    # Mark best-k
    best_idx = df["silhouette"].idxmax()
    best_k = int(df.loc[best_idx, "k"])
    best_sil = float(df.loc[best_idx, "silhouette"])
    ax1.axvline(best_k, color="#d62728", linestyle="--", linewidth=1.5, alpha=0.7,
                label=f"best k={best_k} (silhouette={best_sil:.3f})")

    # Mark PAM50-cardinality reference
    ax1.axvline(5, color="gray", linestyle=":", linewidth=1.5, alpha=0.7,
                label="k=5 (PAM50 cardinality)")

    ax2 = ax1.twinx()
    ax2.plot(df["k"], df["inertia"], "s-", color="#ff7f0e",
             markersize=6, linewidth=1.5, alpha=0.7, label="inertia")
    ax2.set_ylabel("Inertia (WCSS)", color="#ff7f0e")
    ax2.tick_params(axis="y", labelcolor="#ff7f0e")

    ax1.set_title("k-means on MOFA+ factor scores: silhouette + inertia over k\n"
                  "Weak globular structure beyond basal-vs-rest split (k=2)")
    # Combined legend
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def fig_supervised_vs_unsupervised_brca(out_path: pathlib.Path) -> None:
    """Phase 3 (supervised XGBoost) vs Phase 4 (unsupervised MOFA+) on BRCA.

    Two complementary stories:
      (a) Supervised: train on PAM50 labels, get 0.87 ± 0.01 5-fold CV
          accuracy with strong PAM50 panel recovery (15/20 in top-25).
      (b) Unsupervised: no labels used, recover 97%-pure basal cluster +
          76%-pure LumA cluster from multi-omic factor scores.

    Different paradigms, different metrics — both succeed in their own
    framing.
    """
    with open(DATA / "classification_pipeline_brca" / "classification_metrics.json") as f:
        ph3 = json.load(f)
    with open(DATA / "mofa_BRCA" / "cluster_vs_subtype.json") as f:
        ph4 = json.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # Left panel: Phase 3 supervised
    ax_l = axes[0]
    metrics = [("accuracy", "Accuracy"),
               ("balanced_accuracy", "Balanced acc."),
               ("macro_f1", "Macro F1")]
    means = [ph3["cv_scores"][k]["mean"] for k, _ in metrics]
    stds = [ph3["cv_scores"][k]["std"] for k, _ in metrics]
    bars = ax_l.bar(range(len(metrics)), means, yerr=stds, capsize=4,
                    color="#1f77b4", edgecolor="white")
    for bar, m in zip(bars, means):
        ax_l.text(bar.get_x() + bar.get_width() / 2, m + 0.02, f"{m:.3f}",
                  ha="center", va="bottom", fontsize=10)
    ax_l.set_xticks(range(len(metrics)))
    ax_l.set_xticklabels([d for _, d in metrics])
    ax_l.set_ylim(0, 1.1)
    ax_l.set_ylabel("5-fold CV score (mean ± std)")
    ax_l.set_title("Phase 3: supervised XGBoost on PAM50\n"
                   "(5 classes, RNA only)")
    ax_l.axhline(1.0, color="black", linewidth=0.5, linestyle=":", alpha=0.4)
    ax_l.grid(True, axis="y", alpha=0.3)

    # Right panel: Phase 4 unsupervised
    ax_r = axes[1]
    # Pull cluster purities from cluster assignments + label join (re-derive)
    clusters_df = pd.read_csv(DATA / "mofa_BRCA" / "cluster_assignments.csv")
    import polars as pl
    labels_df = (
        pl.read_parquet(DATA / "audit" / "subtype_label_audit.parquet")
        .select(["patient_id", "hoadley_subtype_selected"])
        .to_pandas()
        .dropna(subset=["hoadley_subtype_selected"])
    )
    merged = clusters_df.merge(labels_df, on="patient_id", how="inner")
    ct = pd.crosstab(merged["cluster"], merged["hoadley_subtype_selected"])
    purities = []
    cluster_n = []
    cluster_labels = []
    for c in sorted(ct.index):
        row = ct.loc[c]
        top_class = row.idxmax().replace("BRCA.", "")
        pct = 100 * row.max() / row.sum()
        purities.append(pct)
        cluster_n.append(row.sum())
        cluster_labels.append(f"cluster {c}\nn={row.sum()}\nmajority: {top_class}")

    bars = ax_r.bar(range(len(purities)), purities,
                    color=["#d62728" if p >= 90 else "#1f77b4" if p >= 70 else "#aaaaaa"
                           for p in purities],
                    edgecolor="white")
    for bar, p in zip(bars, purities):
        ax_r.text(bar.get_x() + bar.get_width() / 2, p + 1.5, f"{p:.0f}%",
                  ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax_r.set_xticks(range(len(purities)))
    ax_r.set_xticklabels(cluster_labels, fontsize=9)
    ax_r.set_ylim(0, 110)
    ax_r.set_ylabel("Majority-class purity within cluster (%)")
    title_l2 = (f"ARI={ph4['agreement_vs_known']['ari']:.2f}, "
                f"NMI={ph4['agreement_vs_known']['nmi']:.2f} vs PAM50")
    ax_r.set_title(f"Phase 4: unsupervised MOFA+ k=4 clustering\n"
                   f"(6 modalities, no labels used) — {title_l2}")
    ax_r.axhline(100, color="black", linewidth=0.5, linestyle=":", alpha=0.4)
    ax_r.grid(True, axis="y", alpha=0.3)

    fig.suptitle("BRCA subtyping: supervised vs unsupervised paradigms",
                 fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIGURES_DIR.mkdir(exist_ok=True)
    fig_modality_coverage(FIGURES_DIR / "01_brca_modality_coverage.png")
    print(f"[fig 1] {FIGURES_DIR / '01_brca_modality_coverage.png'}")

    fig_classifier_accuracy(FIGURES_DIR / "02_classifier_accuracy_comparison.png")
    print(f"[fig 2] {FIGURES_DIR / '02_classifier_accuracy_comparison.png'}")

    fig_shap_panel_hits(FIGURES_DIR / "03_shap_panel_hits.png")
    print(f"[fig 3] {FIGURES_DIR / '03_shap_panel_hits.png'}")

    fig_silhouette_curve(FIGURES_DIR / "04_silhouette_curve.png")
    print(f"[fig 4] {FIGURES_DIR / '04_silhouette_curve.png'}")

    fig_supervised_vs_unsupervised_brca(FIGURES_DIR / "05_supervised_vs_unsupervised_brca.png")
    print(f"[fig 5] {FIGURES_DIR / '05_supervised_vs_unsupervised_brca.png'}")

    print(f"\nDone. Project-overview figures in {FIGURES_DIR}/")
    print(f"Phase 4-specific figures still at data/mofa_BRCA/figures/")


if __name__ == "__main__":
    main()
