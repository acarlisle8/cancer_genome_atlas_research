# Phase 02 Merge — Verification Checklist

Practical acceptance checks to run **after** `python run_merge.py` finishes,
to confirm `data/merged_all_cohorts.parquet` is correct end-to-end (not just
"exit code 0 with garbage").

Run the script below in a shell from project root. Every assertion should
print `OK`. Anything else is a real issue worth investigating before moving
on to Phase 03 (XGBoost + SHAP).

---

## One-shot verification script

```python
from pathlib import Path
import polars as pl

OUT = Path("data/merged_all_cohorts.parquet")

# 1. File exists and opens
assert OUT.exists(), f"missing {OUT}"
df = pl.read_parquet(OUT)
print(f"OK file opens — shape {df.shape}")

# 2. Reasonable row count (~1,000 patients across BRCA + LUAD + PRAD)
assert 500 < len(df) < 3000, f"row count looks wrong: {len(df)}"
print(f"OK row count {len(df)}")

# 3. Reasonable column count (500 genes + ~46 arms + 500 probes + patient_id + cohort)
assert 1000 < len(df.columns) < 1100, f"column count looks wrong: {len(df.columns)}"
print(f"OK column count {len(df.columns)}")

# 4. Cohort column has exactly the 3 expected values
cohorts = set(df["cohort"].unique().to_list())
assert cohorts == {"BRCA", "LUAD", "PRAD"}, f"cohorts wrong: {cohorts}"
print(f"OK cohort values {cohorts}")

# 5. patient_id is unique (one row per patient)
assert df["patient_id"].n_unique() == len(df), "duplicate patient_ids"
print(f"OK patient_id uniqueness")

# 6. Per-cohort patient counts not wildly imbalanced
counts = df.group_by("cohort").len().sort("len")
print(f"   cohort counts:\n{counts}")
mn, mx = counts["len"].min(), counts["len"].max()
assert mx / mn < 5, f"cohort imbalance suspicious: min={mn} max={mx}"
print(f"OK cohort balance (ratio {mx/mn:.1f})")

# 7. Gene FPKM columns: nonneg, mostly small positives
gene_cols = [c for c in df.columns if c.startswith("ENSG")]
assert len(gene_cols) > 100, f"too few gene columns: {len(gene_cols)}"
gene_min = df[gene_cols].min().min_horizontal().item()
assert gene_min is None or gene_min >= 0, f"negative FPKM found: {gene_min}"
print(f"OK {len(gene_cols)} gene columns, all >=0")

# 8. Methylation beta columns: in [0, 1]
probe_cols = [c for c in df.columns if c.startswith("cg")]
assert len(probe_cols) > 100, f"too few probe columns: {len(probe_cols)}"
probe_min = df[probe_cols].min().min_horizontal().item()
probe_max = df[probe_cols].max().max_horizontal().item()
assert (probe_min is None or probe_min >= 0) and (probe_max is None or probe_max <= 1), \
    f"beta out of [0,1]: min={probe_min} max={probe_max}"
print(f"OK {len(probe_cols)} probe columns in [0,1]")

# 9. CNV arm columns: log2 copy ratios, typically -3..3
arm_cols = [c for c in df.columns if c[0].isdigit() and (c.endswith("p") or c.endswith("q"))]
assert 30 <= len(arm_cols) <= 50, f"unexpected arm column count: {len(arm_cols)}"
arm_max = df[arm_cols].max().max_horizontal().item()
assert arm_max is None or abs(arm_max) < 10, f"CNV log-ratio crazy: max={arm_max}"
print(f"OK {len(arm_cols)} arm columns")

# 10. No column entirely null
all_null = [c for c in df.columns if df[c].null_count() == len(df)]
assert not all_null, f"all-null columns: {all_null}"
print(f"OK no all-null columns")

# 11. File size sanity
size_mb = OUT.stat().st_size / 1e6
assert 50 < size_mb < 5000, f"file size suspicious: {size_mb:.0f} MB"
print(f"OK file size {size_mb:.0f} MB")

print("\n=== ALL 11 CHECKS PASSED ===")
```

---

## What to investigate if a check fails

| Check | If it fails, suspect |
|-------|---------------------|
| 1 (exists/opens) | merge crashed — check log or re-run with verbose |
| 2 (row count) | inner-join too restrictive (small) or no dedup (large) |
| 3 (col count) | top-N variance computation didn't pick 500 of each |
| 4 (cohort values) | cohort labelling broken in `merge_all_cohorts` |
| 5 (patient_id unique) | duplicate rows from join — Aidan's RuntimeError guard should catch |
| 6 (cohort balance) | one cohort lost most patients in the inner join — look at per-cohort patient counts in raw Parquets |
| 7 (FPKM nonneg) | parsing went wrong somewhere upstream |
| 8 (beta in [0,1]) | methylation parsing or scaling issue |
| 9 (arm column count) | `_pivot_cnv` filtered too aggressively or chrM leaked in |
| 10 (no all-null) | a feature got selected but the column ended up empty after pivot — bug in top-N selection or filtering |
| 11 (file size) | obvious truncation or empty Parquet |

---

_Owner: Zachary_
_Created: 2026-04-28 in support of running merge end-to-end on EC2 after the polars ingest landed_
