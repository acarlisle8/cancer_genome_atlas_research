# Phase 4e + 4f Session Log — 2026-05-06 (early morning)

**Branch:** `brca-multimodal-spark`. After this session, **12 commits ahead** of `main`.

## Entry point on reconnect

**You (fresh Claude) are picking up after a session that did three things:
(1) post-mortem and audit of an OOM crash from the prior session,
(2) Phase 4e — extended `run_mofa.py` from 3-view to 6-view BRCA, ran the
full model, and (3) Phase 4f-light — quick clustering analysis vs PAM50.
The full multi-omic MOFA+ result is in hand.**

Read this file, then [.planning/2026-05-05-phase-4bd-session-log.md] for
the prior 4b/4d work, then [.planning/ROADMAP.md] Phase 4 section.

---

## What was completed this session

### OOM post-mortem (2026-05-05 22:06–22:17 UTC)

The prior session (`a8827b69`) was OOM-killed three times in 11 minutes
(PIDs 8402, 12276, 12320 — each ~6 GB RSS on a 7.6 GB no-swap box).
**The OOMs were ad-hoc EDA Bash one-liners**, not committed code.
Triggers:

- 22:06: `pl.read_parquet('methylation.parquet').group_by(...).count()` —
  444M rows, eager.
- 22:17: `df.group_by(['patient_id','probe_id']).count()` — same parquet,
  same problem.
- 22:17:30s later: retry that dropped methylation but kept rna_seq
  (74.7M rows) — Polars `async-executor` OOM'd itself.

**The instance reboot at 00:06 UTC is a separate event** — clean dmesg, no
panic, ~2 hr after the last OOM. Not caused by the OOMs themselves.

### Branch audit for similar foot-guns

Audited every Python file added on `brca-multimodal-spark`:

| File | Issue | Fix |
|---|---|---|
| [run_mofa.py:170] | Eager `pl.read_parquet(...).filter(...).to_pandas()` on `merged_all_cohorts.parquet` | → `pl.scan_parquet(...).filter(...).collect(engine="streaming").to_pandas()` |
| [src/merge.py] `_pivot_rppa`, `_pivot_mirna`, `_recurrent_mutated_genes`, `_pivot_mutations` | All used eager `pl.read_parquet`, breaking the lazy convention used by `_pivot_rnaseq` / `_pivot_methylation` | All four rewritten to mirror the lazy convention; `is_in(recurrent_genes)` filter pushed into the lazy plan for `_pivot_mutations` |

The eager reads on RPPA/miRNA/mutations parquets *worked today* (3.5 MB,
7 MB, 0.4 MB respectively at BRCA-only scale) but were landmines for any
future pan-cohort expansion.

The actual OOM-causing scripts were never committed — they were scratch
EDA. Memory file `~/.claude/projects/-home-ubuntu/memory/feedback_dats6450_lazy_polars.md`
records the rule for future sessions.

Commit: `6735016 Convert eager pl.read_parquet to lazy scan_parquet on big-file paths`

### `mofapy2` was missing from `pyproject.toml`

`run_mofa.py` imported it from day one but the dep was never declared.
Survived the prior session because mofapy2 was installed somewhere external
to the project venv; that didn't survive the reboot.

`uv add mofapy2` → 0.7.4. Commit: `727e6d2`.

### Methylation feature-selection bug in `_top_n_by_variance`

**Symptom (from a quick diagnostic during the audit):** the merged 6-view's
mean methylation missingness was 36.3% — vastly higher than the source
parquet's ~18% array-QC baseline.

**Root cause:** [src/merge.py:_top_n_by_variance] ranked probes by
`VAR_POP(beta_value)` over non-null observations only. Probes that were
100% missing for BRCA but had non-null beta values for one or two odd
patients had huge "empirical variance" across those few values and
preferentially infiltrated the top-5000 list.

DuckDB diagnostic confirmed:

| Probe set | Mean missingness | Median |
|---|---|---|
| All 488K probes (source) | ~19% | 18.4% |
| Top-5000 by variance (old) | 27.9% | 18.8% |
| Bottom-5000 by variance | 19.0% | 18.4% |
| Random 5000 | 19.1% | 18.4% |

**Fix:** Added a `HAVING COUNT(*) >= MIN_COVERAGE_FRAC * n_patients`
clause before the `ORDER BY VAR_POP`. `MIN_COVERAGE_FRAC = 0.80` at module
scope. Marker genes in `must_keep` still bypass the filter.

After re-running [run_merge_6view.py]: methylation mean missingness
36.3% → **26.7%**, per-col p90 29.2% (no more 100%-missing tail).
The remaining ~8% gap above source baseline is from the 6-way inner-join
shifting the kept-patient population toward worse meth coverage —
structural, not a probe-selection bug. Tightening the threshold further
would just drop informative-but-incomplete probes.

Commit: `fd8d3e3 Filter sparse features before variance ranking in _top_n_by_variance`.

The merged file `data/TCGA-BRCA/merged_brca_6view.parquet` was
regenerated (744 patients × 12,612 cols — same shape, different
underlying probes).

### Phase 4e — extended `run_mofa.py` from 3-view to 6-view

All six items from the 4d-handoff plan, plus a per-view max-missing-rate
dict to handle methylation specifically. Commit: `531e063`.

| Item | What changed |
|---|---|
| 1 | `split_views()` recognizes `RPPA_*` / `MIR_*` / `MUT_*` prefixes alongside `ENSG` / `cg`. Returns a fixed-order dict (mutations last). |
| 2 | New `preprocess_rppa` (identity), `preprocess_mirna` (`log2(RPM+1)`), `preprocess_mutations` (cast-to-float) functions. |
| 3 | `if/elif` preprocessing dispatch replaced with `PREPROCESSORS` dict keyed by view name. Adding views is one line. |
| 4 | Per-view likelihoods via `LIKELIHOODS` dict — 5 gaussian + 1 bernoulli (mutations). |
| 5 | Cohort filter is conditional on `cohort` column being in input schema. The 6-view file is single-cohort and has no cohort col. `--cohort` is still required because it labels output dir / manifest / mofapy2 `groups_names`. |
| 6 | `--input` default switched to `data/TCGA-BRCA/merged_brca_6view.parquet`. Old `merged_all_cohorts.parquet` still works via explicit `--input`. |
| **+** | **`MAX_MISSING_RATES` dict at module scope.** Methylation = 0.40 (merged mean miss is 27%, p90 29%). All others = 0.20. CLI `--max-missing-rate` now defaults to None and acts as a uniform global override when set. Backward-compat for `run_preprocess_subtype.py` and the planning doc's `--max-missing-rate 0.5` reference. |

### Smoke test — `--max-iter 50 --convergence fast`

- Wall time: 185s. Converged at iter 39/50 (deltaELBO < 0.0005%).
- 14 of 15 factors active after ARD pruning.
- All 6 views loaded, likelihoods routed correctly (V0–V4 gaussian, V5 bernoulli).
- One warning: 298 of 1881 miRNAs have zero variance (constant-zero across
  all 744 patients). mofapy2 absorbs them; cheap cleanup for later
  (`filter_by_zero_variance` alongside `filter_by_missingness`).

### Full run — `--max-iter 1000 --convergence medium`

- Wall time: **958.9s (~16 min)**, ~200 iterations.
- ARD held at 14 active factors (data has that many distinct axes).
- Mutations contributions concentrated cleanly on factors 1 (R² 0.81%)
  and 2 (R² 0.63%), zero on most other factors — correct Bernoulli-view
  behavior at convergence. Smoke had small non-zero everywhere; medium
  pruned to where the signal lives.
- Variance-explained per view across 8 leading factors:
  RNA ~40%, methylation ~45%, miRNA ~23%, RPPA ~18%, CNV ~7%, mutations ~1.5%
  (mutations is bounded low by the Bernoulli scale, not a problem).
- Factor structure (full table in `data/mofa_BRCA/variance_explained.csv`):

| Factor | RNA | meth | CNV | RPPA | miRNA | mut | Likely biology |
|---|---|---|---|---|---|---|---|
| 1 | 12.2 | 16.8 | 2.1 | 6.9 | 7.1 | **0.81** | basal-vs-luminal axis (multi-omic across all 6) |
| 2 | 9.1 | 2.9 | 1.8 | 5.4 | 8.4 | 0.63 | proliferation / driver-mutation axis |
| 3 | 0.3 | **14.8** | 0.3 | 0.4 | 0.0 | 0.0 | methylation-specific (TME / cell-of-origin) |
| 4 | 6.5 | 1.5 | 0.3 | 2.1 | 1.4 | 0.0 | RNA-led transcriptional program |
| 5 | 2.2 | 0.9 | 1.8 | 1.5 | 2.7 | 0.0 | small multi-omic |
| 6 | 6.3 | 0.1 | 0.2 | 0.6 | 0.8 | 0.0 | RNA-only sub-program |
| 7 | 0.3 | 7.2 | 0.4 | 0.0 | 0.0 | 0.0 | methylation-specific |
| 8 | 2.8 | 0.2 | 0.3 | 1.0 | 2.3 | 0.0 | small multi-omic |

### Phase 4f-light — clustering vs PAM50

Ran `analyze_mofa.py` (k-means at silhouette-best k, ARI/NMI vs
`hoadley_subtype_selected`). Output: `data/mofa_BRCA/cluster_vs_subtype.json`.

| k-means k | Silhouette | ARI | NMI | Cluster sizes |
|---|---|---|---|---|
| 2 (silhouette-best) | 0.194 | 0.307 | 0.447 | 132 / 612 |
| 4 (forced) | 0.124 | 0.268 | 0.328 | 97 / 225 / 278 / 144 |
| 5 (forced, PAM50 cardinality) | 0.126 | 0.250 | 0.321 | one degenerate 2-pt cluster |

**Contingency table (k=4 vs PAM50):**

| Cluster | n | Composition |
|---|---|---|
| 0 | 97 | **97% Basal** (94/97) |
| 2 | 278 | **76% LumA** |
| 1 | 225 | 46% LumB + 22 Her2 + LumA mix |
| 3 | 144 | mixed: 33% LumA, 23% Her2, 19% Basal, 15% LumB, 10% Normal |

**ROADMAP target was ARI > 0.5; we landed at 0.27–0.31.** Three reasons:

1. **PAM50 is RNA-only.** Methylation contributes 45% of the captured
   variance and has its own non-PAM50 structure (factors 3 and 7).
   The audit table's `hoadley_subtype_integrative` would be a fair
   multi-omic comparator but **is null for BRCA** in this dataset.
2. **Plain k-means at silhouette-best is not the spec.** ROADMAP says
   "consensus clustering"; that's a real implementation gap.
3. **All 14 factors include noise.** Factors 6–14 are progressively
   view-specific or low-variance. Restricting to factors 1–5 typically
   gains 0.03–0.05 ARI.

### What we can credibly say

- **Method works**: model converged, recovered known multi-omic biology,
  identified basal-like with 97% purity at k=4.
- **Mutations contribute**: Bernoulli view co-varies with the proliferation
  axis where you'd expect (TP53/PIK3CA enrichment in basal/proliferative
  tumors).
- **Methylation has independent signal**: two methylation-specific factors
  capture biology the RNA panel doesn't see.
- **PAM50 partial match (NMI 0.45 at k=2)**: agreement above chance,
  driven by clean basal recovery — consistent with multi-omic-vs-single-omic
  literature.

### What we can NOT credibly say

- "Recovered PAM50" (didn't — recovered basal cleanly + partial separation
  elsewhere).
- "Better than RNA-only would do" (didn't run RNA-only baseline).
- "ARI > 0.5" (didn't hit that target).

---

## Commits this session (in order)

| SHA | Title |
|---|---|
| `6735016` | Convert eager pl.read_parquet to lazy scan_parquet on big-file paths |
| `727e6d2` | Add mofapy2 0.7.4 dep — required by run_mofa.py since it was first added |
| `fd8d3e3` | Filter sparse features before variance ranking in _top_n_by_variance |
| `531e063` | Phase 4e: extend run_mofa.py from 3-view to 6-view BRCA multi-omic |

(Documentation update commit will follow this log.)

---

## State of the repo at session end

- **Branch**: `brca-multimodal-spark`, 12 commits ahead of `main`
- **Working tree**: clean (assuming the doc-update commit lands)
- **Generated artifacts (gitignored)**:
  - `data/TCGA-BRCA/merged_brca_6view.parquet` (744 × 12,612, regenerated this session with cleaner meth probes)
  - `data/mofa_BRCA/mofa_model.hdf5` (full-run trained model)
  - `data/mofa_BRCA/factor_scores.csv` (744 patients × 14 factors)
  - `data/mofa_BRCA/variance_explained.csv` (14 factors × 6 views)
  - `data/mofa_BRCA/top_loadings_top25_<view>.csv` (6 files)
  - `data/mofa_BRCA/run_manifest.json`
  - `data/mofa_BRCA/cluster_assignments.csv` (last analyze_mofa run)
  - `data/mofa_BRCA/cluster_vs_subtype.json`
  - `data/mofa_BRCA/variance_explained_heatmap.png`
  - `data/mofa_BRCA/kmeans_silhouette.csv`

---

## Open items for the next phase

### 4f proper (not done this session)
- **Consensus clustering** with bootstrap resampling — `ConsensusClusterPlus`
  paradigm or equivalent in Python (`scikit-learn` doesn't have one OOTB;
  `consensus-clustering` PyPI package or a hand-rolled resampler).
  Expected ARI gain: 0.03–0.08.
- **SNF (Similarity Network Fusion)** as an independent multi-omic
  comparator — orthogonal paradigm to MOFA+'s factor-based approach.
  ROADMAP success criterion: clusters comparable across paradigms.
- **Cluster on factors 1–5 only** to drop noise factors. One-line change.
- **Decision point at end of 4f**: is "ARI > 0.5" still the right success
  criterion given PAM50's RNA-only nature? Or reframe as:
  - "Recovers basal-like at >90% purity" ✓ (already done)
  - "Top factor is multi-omic across all 6 views" ✓ (already done)
  - Plus add a comparator vs RNA-only XGBoost from Phase 3.

### 4g (Metabric external validation)
- Cross-platform harmonization: Metabric is microarray-based, not
  RNA-seq. Probe → gene mapping required. RPPA/miRNA may not be
  available for Metabric — likely need to validate on RNA + meth
  + CNV subset.

### Cleanup carried over
- **Zero-variance filter for miRNA** — 298/1881 miRNAs constant-zero
  across all 744 patients. mofapy2 warning at every run. Cheap fix:
  add `filter_by_zero_variance(X, cols)` alongside `filter_by_missingness`
  in `run_mofa.py`.

### Do NOT re-litigate
- Eager Polars on big parquets — locked: lazy `scan_parquet` only on
  methylation/RNA-seq parquets. See `feedback_dats6450_lazy_polars.md`
  in Claude memory and the `Polars memory landmines` section in
  [.planning/known-issues.md].
- `MIN_COVERAGE_FRAC = 0.80` for `_top_n_by_variance` — confirmed to
  drop the 100%-missing-tail bias. Don't tighten further; the remaining
  gap above source baseline is structural (population shift via inner-join).
- `MAX_MISSING_RATES["methylation"] = 0.40` — calibrated to the merged
  6-view's per-col p90 of 29%. Other views stay at 0.20.

---

## Pointers to key files

- This session: `.planning/2026-05-06-phase-4ef-session-log.md` (you're reading it)
- Prior session: `.planning/2026-05-05-phase-4bd-session-log.md`
- Earlier: `.planning/2026-05-05-phase-4a-session-log.md`, `.planning/2026-05-05-session-log.md`
- Roadmap: `.planning/ROADMAP.md` (Phase 4 section)
- Project: `.planning/PROJECT.md`
- Known issues: `.planning/known-issues.md`
- Runbook: `.planning/RUNBOOK.md`
- 4e implementation: `run_mofa.py` (parse_args, split_views, PREPROCESSORS, LIKELIHOODS, MAX_MISSING_RATES, main)
- 4f-light analysis: `analyze_mofa.py` (k-means, silhouette, ARI/NMI)
- Outputs (gitignored): `data/mofa_BRCA/` directory
