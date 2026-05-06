"""Compare top SHAP features per subtype model against published marker panels.

For each cancer (BRCA, PRAD, LUAD):
  - Reads shap_values_summary.csv
  - Classifies each feature as RNA / methylation / CNV by name pattern
  - For RNA features, matches Ensembl ID prefix against the panel in src/markers
  - Reports overall top-K, per-class top-K, panel rank table, and hit rates

Outputs both a printed summary and a JSON report under
data/classification_pipeline_<task>/shap_panel_comparison.json
"""
from __future__ import annotations

import json
import pathlib
from typing import Iterable

import polars as pl

from src.markers import PANELS

TASKS = ["BRCA", "PRAD", "LUAD"]

PUBLISHED_NOTES = {
    "BRCA": "PAM50 (Parker 2009; TCGA Nature 2012). Subtypes: LumA, LumB, Her2, Basal, Normal-like.",
    "PRAD": "TCGA Cell 2015 (iCluster + AR/ERG/FOXA1/SPOP markers).",
    "LUAD": "Wilkerson 2012 (TRU / proximal-inflammatory / proximal-proliferative).",
}


def feature_kind(name: str) -> str:
    if name.startswith("ENSG"):
        return "rna"
    if name.startswith("cg"):
        return "methylation"
    return "cnv"  # chromosome arms (1p, 1q, ...)


def ensembl_prefix(feat: str) -> str:
    """Strip the version suffix from an Ensembl ID."""
    return feat.split(".", 1)[0]


def panel_for(cancer: str) -> dict[str, str]:
    return PANELS.get(cancer, {})


def panel_lookup(cancer: str) -> dict[str, str]:
    """Reverse: ENSG_id -> gene_symbol."""
    return {v: k for k, v in panel_for(cancer).items()}


def annotate(feat: str, panel_rev: dict[str, str]) -> tuple[str, str]:
    """Return (kind, label). Label is gene symbol if RNA-and-in-panel, else feat."""
    kind = feature_kind(feat)
    if kind == "rna":
        sym = panel_rev.get(ensembl_prefix(feat))
        return kind, sym if sym else feat
    return kind, feat


def topk_table(df: pl.DataFrame, k: int, panel_rev: dict[str, str]) -> list[dict]:
    rows = []
    for r in df.head(k).iter_rows(named=True):
        kind, label = annotate(r["feature"], panel_rev)
        in_panel = kind == "rna" and ensembl_prefix(r["feature"]) in panel_rev
        rows.append({
            "rank": len(rows) + 1,
            "feature": r["feature"],
            "label": label,
            "kind": kind,
            "in_panel": in_panel,
            "mean_abs_shap": float(r["mean_abs_shap"]),
        })
    return rows


def panel_rank_table(df: pl.DataFrame, panel_rev: dict[str, str]) -> list[dict]:
    """For each panel gene, find its rank in the SHAP ranking (or None if absent)."""
    feats = df["feature"].to_list()
    feat_idx: dict[str, int] = {}
    for i, f in enumerate(feats):
        if f.startswith("ENSG"):
            feat_idx.setdefault(ensembl_prefix(f), i)

    rows = []
    for ensg, sym in panel_rev.items():
        idx = feat_idx.get(ensg)
        rows.append({
            "gene": sym,
            "ensg": ensg,
            "rank": idx + 1 if idx is not None else None,
            "in_features": idx is not None,
            "mean_abs_shap": float(df["mean_abs_shap"][idx]) if idx is not None else None,
        })
    rows.sort(key=lambda x: (x["rank"] is None, x["rank"] or 10**9))
    return rows


def hit_rate(panel_rev: dict[str, str], top_rows: list[dict]) -> dict:
    panel_size = len(panel_rev)
    hits = sum(1 for r in top_rows if r["in_panel"])
    return {"panel_size": panel_size, "hits": hits, "rate": hits / panel_size if panel_size else 0.0}


def per_class_top(df: pl.DataFrame, panel_rev: dict[str, str], top_n: int) -> dict[str, list[dict]]:
    """For each per-class column, return the top-N features by that class's |SHAP|."""
    per_class = {}
    for col in df.columns:
        if not col.startswith("mean_abs_shap_"):
            continue
        cls = col[len("mean_abs_shap_"):]
        ranked = df.sort(col, descending=True).head(top_n)
        rows = []
        for r in ranked.iter_rows(named=True):
            kind, label = annotate(r["feature"], panel_rev)
            in_panel = kind == "rna" and ensembl_prefix(r["feature"]) in panel_rev
            rows.append({"feature": r["feature"], "label": label, "kind": kind,
                         "in_panel": in_panel, "mean_abs_shap_class": float(r[col])})
        per_class[cls] = rows
    return per_class


def fmt_top_block(title: str, rows: Iterable[dict]) -> str:
    out = [f"  {title}"]
    for r in rows:
        marker = "*" if r.get("in_panel") else " "
        out.append(f"    {marker} {r['rank']:>2}. [{r['kind']:>11}] {r['label']:<22} shap={r['mean_abs_shap']:.4f}")
    return "\n".join(out)


def fmt_panel_rank_block(rows: list[dict]) -> str:
    out = ["  panel-gene ranks (in SHAP overall)"]
    for r in rows:
        rank = r["rank"]
        rank_s = f"{rank:>4}" if rank is not None else "  --"
        shap_s = f"{r['mean_abs_shap']:.4f}" if r["mean_abs_shap"] is not None else "   -- "
        miss = "" if r["in_features"] else "  (not in selected pool)"
        out.append(f"    rank {rank_s}  {r['gene']:<8} shap={shap_s}{miss}")
    return "\n".join(out)


def main() -> None:
    out_lines: list[str] = []
    out_lines.append("=" * 72)
    out_lines.append("Top SHAP features vs published marker panels")
    out_lines.append("=" * 72)

    for cancer in TASKS:
        ct = cancer.lower()
        path = pathlib.Path(f"data/classification_pipeline_{ct}/shap_values_summary.csv")
        if not path.exists():
            out_lines.append(f"\n[skip] {cancer}: {path} not found")
            continue

        df = pl.read_csv(path)
        panel_rev = panel_lookup(cancer)
        top10 = topk_table(df, 10, panel_rev)
        top25 = topk_table(df, 25, panel_rev)
        top50 = topk_table(df, 50, panel_rev)

        rate10 = hit_rate(panel_rev, top10)
        rate25 = hit_rate(panel_rev, top25)
        rate50 = hit_rate(panel_rev, top50)
        panel_ranks = panel_rank_table(df, panel_rev)
        per_class = per_class_top(df, panel_rev, top_n=8)

        out_lines.append(f"\n{'-'*72}")
        out_lines.append(f"{cancer}  — {PUBLISHED_NOTES[cancer]}")
        out_lines.append(f"  panel size: {len(panel_rev)} genes")
        out_lines.append(f"  pool: {df.height} features ranked")
        out_lines.append(f"  panel hits in top10/25/50: "
                         f"{rate10['hits']}/{rate10['panel_size']} | "
                         f"{rate25['hits']}/{rate25['panel_size']} | "
                         f"{rate50['hits']}/{rate50['panel_size']}")
        out_lines.append(fmt_top_block("top 10 overall (by mean |SHAP|, * = in panel)", top10))
        out_lines.append(fmt_panel_rank_block(panel_ranks))

        out_lines.append("  per-class top 8 (by class mean |SHAP|, * = in panel):")
        for cls, rows in per_class.items():
            out_lines.append(f"    [{cls}]")
            for i, r in enumerate(rows, 1):
                marker = "*" if r["in_panel"] else " "
                out_lines.append(f"      {marker} {i:>2}. [{r['kind']:>11}] {r['label']:<22} shap={r['mean_abs_shap_class']:.4f}")

        report = {
            "cancer": cancer,
            "published_reference": PUBLISHED_NOTES[cancer],
            "panel_size": len(panel_rev),
            "panel_genes": list(panel_rev.values()),
            "n_features_ranked": df.height,
            "hits_top10": rate10,
            "hits_top25": rate25,
            "hits_top50": rate50,
            "top25_overall": top25,
            "panel_gene_ranks": panel_ranks,
            "per_class_top8": per_class,
        }
        report_path = pathlib.Path(f"data/classification_pipeline_{ct}/shap_panel_comparison.json")
        report_path.write_text(json.dumps(report, indent=2) + "\n")
        out_lines.append(f"  report: {report_path}")

    text = "\n".join(out_lines)
    print(text)
    pathlib.Path("data/shap_panel_comparison_summary.txt").write_text(text + "\n")


if __name__ == "__main__":
    main()
