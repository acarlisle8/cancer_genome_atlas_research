# 2026-05-05 — Session log

Subtype classification rewrite + SHAP-vs-panel comparison + MOFA+ prep.
Long evening session. Picks up from a prior session that disconnected mid-merge.

---

## Arc of the night

1. **Resumed mid-merge after disconnect.** Prior session was running `run_merge.py`
   with the new variance pool (500→5000) and force-marker pipeline. Process had
   died mid-PRAD methylation pivot. Restarted from scratch (~14 min).

2. **Re-ran the preprocessing chain.** `run_preprocess.py` (cohort) → 2104 × 7419,
   then `run_preprocess_subtype.py` → BRCA / PRAD / LUAD model-ready parquets.

3. **Trained 4 XGBoost models** with the new 5-fold CV + per-fold variance +
   force-marker pipeline. Sequential, ~8.5 min total on the 2-core box.

4. **Wrote `compare_shap_to_panels.py`** to map top SHAP features against
   published marker panels (PAM50 / Wilkerson / TCGA Cell PRAD).

5. **Cleaned up git history** into 3 logical commits, opened PR #6, squash-merged
   to main. Hit GitHub's 100 MB file limit; reverted parquets out of tracking via
   gitignore amendment.

6. **Squash-merge accidentally deleted both regenerated parquets** during pull
   (mechanism explained below). Regenerated via run_merge.py + run_preprocess.py.

7. **MOFA+ prep for tomorrow.** Installed mofapy2/mofax, wrote `run_mofa.py` +
   `analyze_mofa.py`, validated end-to-end with a smoke run.

---

## Results

### XGBoost classifiers (5-fold stratified CV, per-fold variance + force-marker, 500 features per fold)

| Task | Accuracy | Balanced Acc | Macro F1 | Log Loss | Classes | Marker panel |
|---|---|---|---|---|---|---|
| Cohort | 0.9990 ± 0.0012 | 0.9987 ± 0.0016 | 0.9990 ± 0.0012 | 0.0036 ± 0.0016 | 3 | none (cohort task) |
| BRCA subtype | 0.8675 ± 0.0094 | 0.7927 ± 0.0344 | 0.7982 ± 0.0274 | 0.3846 ± 0.0413 | 5 | BRCA (20 PAM50 genes) |
| PRAD subtype | 0.9279 ± 0.0117 | 0.8909 ± 0.0323 | 0.8993 ± 0.0233 | 0.2119 ± 0.0119 | 4 | PRAD (9 genes) |
| LUAD subtype | 0.5913 ± 0.0733 | 0.5747 ± 0.0824 | 0.5746 ± 0.0773 | 1.1247 ± 0.1274 | 6 | LUAD (7 genes) |

### Cohort vs Aidan's v1 baseline (commit `2d65eb1`, single train/test split, 300 features global)

| Metric | v1 | new | delta |
|---|---|---|---|
| Accuracy | 0.9929 | 0.9990 | +0.61 pp |
| Macro F1 | 0.9928 | 0.9990 | +0.62 pp |
| Log loss | 0.0138 | 0.0036 | −74% |

Caveat: not apples-to-apples (different split, feature count, *and* data version
since v1 was on the pre-CNV-bug-fix data). Real story is "no regression and a
proper CV-based number to report" rather than the +0.6 pp.

### SHAP vs published panels

**BRCA / PAM50:** strong validation.
- 17/20 panel genes in top-50, 9/20 in top-10.
- Top-10 SHAP is essentially PAM50: ESR1, FOXA1, FOXC1, BIRC5, EGFR, ERBB2,
  AURKA, KRT14, CDC20.
- Per-class drivers match published intrinsic-subtype biology
  (Basal: FOXA1↓ + FOXC1↑; Her2: ERBB2 #1; LumA: AURKA/CDC20/ESR1/BCL2;
  LumB: ESR1/BIRC5/EGFR/KRT17; Normal-like: KRT14/KRT5/EGFR — flagging the
  TCGA "Normal-like is debated" caveat).

**PRAD / TCGA Cell 2015:** partial, but two interesting non-hits worth keeping
for the writeup.
- ERG (#1) and ETV1 (#2) lead overall — these are the labels themselves
  (PRAD.1-ERG, PRAD.2-ETV1).
- **AR ranks 407/507 with ~0 mean |SHAP|.** Not a model failure: AR is
  universally expressed in PRAD, so it doesn't discriminate between subtypes
  even though it's the defining lineage marker. SHAP measures discrimination,
  not biological centrality.
- **SPOP ranks 191** despite being a labeled subtype: SPOP defines a
  *mutational* subtype, not an expression subtype, so SPOP expression itself
  doesn't shift between classes.

**LUAD / Wilkerson:** weak.
- Only NKX2-1 (TTF-1, terminal respiratory unit marker) in top-50 of the
  7-gene panel; dominates LUAD.6 SHAP.
- Other Wilkerson genes (SFTPB, SFTPC, TP63, KRT5/6A) at ranks 226-430 with
  near-zero |SHAP|. Consistent with the model's overall 0.59 acc — 230
  patients × 6 classes is too small.

---

## What got committed (PR #6, squash-merged as `3c17b30`)

Three logical commits before squash:

| Commit | What |
|---|---|
| `28d21a2` | merge: force-include subtype marker panels + widen variance pool 500→5000. Adds `src/markers.py` (PAM50 + Wilkerson + TCGA Cell PRAD as Ensembl IDs), threads `must_keep` through `_top_n_by_variance`, fixes anon-S3 reads in ingest scripts via `skip_signature=true`. |
| `640b0b2` (re-made as `2ec7f53` after parquet revert) | classification: 5-fold CV + per-fold variance + force-marker, add subtype task. Rewrites `run_classification_pipeline.py` (CV + per-fold variance + leakage fix + SHAP per-class beeswarms + OOF predictions). Adds `run_preprocess_subtype.py`, `preprocess_for_subtype_model`, `quick_label_audit.py`, cBioPortal cache, TCGASubtype.tsv. New deps: pandas / numpy / scikit-learn / xgboost / matplotlib / shap. Adds regenerated `data/classification_pipeline/*` cohort outputs. |
| `4c717b7` | analysis: SHAP-vs-published-panel comparison. Adds `compare_shap_to_panels.py` and headline finding summaries. |

Squash commit on main: `3c17b30` — "Subtype classification: 5-fold CV +
per-fold variance + marker-panel force-include (#6)".

### What's *not* yet committed (on disk, ready for tomorrow)

- `run_mofa.py` — train MOFA+ on one cohort. Splits merged matrix into
  RNA / methylation / CNV views, applies modality-appropriate preprocessing
  (log2(x+1) on RNA, M-values on methylation, no transform on CNV), drops
  features missing in >max-missing-rate of patients, fits via
  `mofapy2.run.entry_point`. Saves HDF5 + factor scores CSV + variance-
  explained CSV + per-view top-25 loadings.
- `analyze_mofa.py` — post-training: variance heatmap, k-means silhouette
  sweep (k=2..10), ARI/NMI vs published subtype labels.
- mofapy2 0.7.4 / mofax 0.3.7 / h5py / seaborn installed via
  `uv pip` (not yet in pyproject — should be `uv add` once results land).

---

## Storage decisions

- `data/merged_all_cohorts.parquet` (now 117 MB, was 11 MB pre-bug-fix at the
  variance pool of 500) is **gitignored** — it broke GitHub's 100 MB hard cap.
- `data/model_ready_cohort.parquet` (100 MB) likewise gitignored.
- Both regenerate via `run_merge.py` + `run_preprocess.py` (~15 min on 2 cores).
- `data/model_ready_{brca,luad,prad}.parquet` were never tracked (gitignored
  via the `data/*` rule); no change.

---

## Lessons from the night (durable)

1. **Back up working-tree files before checkout/pull when path tracking
   differs.** Lost the freshly-regenerated parquets when checking out main
   after a feature branch had `git rm --cached`'d them. `git checkout main`
   silently overwrote the working-tree 117 MB version with main's old
   tracked 11 MB version, then `git pull` deleted the old version because the
   squash-merge commit recorded a deletion. Both steps were predictable;
   should have `cp` to /tmp first. Saved as memory.

2. **No Co-Authored-By Claude in commits.** Saved as memory.

3. **Sequential vs parallel ML training: check CPU count first.** This box is
   2 cores, each XGBoost configured with `n_jobs=4`. Sequential was the right
   call but I justified it after launching, not before. Should `nproc` /
   `free -h` first.

---

## Open items going into tomorrow

**Active:**
- [ ] Run MOFA+ for BRCA / PRAD / LUAD (~20-30 min each)
- [ ] Run `analyze_mofa.py` for each (silhouette + ARI/NMI vs known subtypes)
- [ ] Decide whether to re-run any cohort with `--max-missing-rate 0.5` to
      keep more methylation features (per-cohort drop is 2368 → 82 for BRCA
      at the default 0.20)
- [ ] Commit the MOFA+ scripts + results once they look right

**Deferred (raised tonight, not actioned):**
- Update `.planning/ROADMAP.md` and `.planning/RUNBOOK.md` to reflect today's
  work (still references the pre-subtype, pre-CV state).
- Decide LUAD subtype treatment for the writeup — collapse rare subtypes,
  drop subtype claim, or report as honest limitation.
- Reproducibility check (fresh-clone-and-run end-to-end).
- Tests for `preprocess_for_subtype_model` and the marker force-include path.
- The writeup itself.

---

## Smoke run results (BRCA, 5 factors, 50 iters, fast convergence)

For sanity-check before the morning runs:
- Trained in 84.4s.
- 4 active factors after ARD pruning.
- Variance explained per view (R², %): factor_1 RNA=10.8 / meth=11.6 / CNV=1.2;
  factor_2 RNA=9.2 / meth=2.8 / CNV=1.6; factor_3 4.5 / 3.4 / 4.6; factor_4 6.5 / 1.6 / 0.3.
- Top loadings on factor_1 are keratins (KRT17, KRT5, KRT8, etc.) — basal-vs-
  luminal axis, exactly the right kind of signal.
- k-means on factor scores: best k=5 by silhouette (0.333). ARI=0.225,
  NMI=0.334 vs `hoadley_subtype_selected`.

Both numbers should improve substantially with full training (15-20 factors,
1000 iter, medium convergence) — the smoke run was just to confirm
end-to-end plumbing.
