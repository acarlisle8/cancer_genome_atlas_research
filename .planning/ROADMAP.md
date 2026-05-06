# Roadmap: TCGA Cancer Type Classifier

## Overview

Phases 1-3 built a single-node DuckDB/Polars pipeline that ingests three modalities for three cohorts, joins them into a feature matrix, and trains XGBoost classifiers (cohort 3-class + per-cohort subtype) with SHAP analysis. Phase 4 extends to multi-omic integration on a real Spark-on-EC2 cluster, scoped to BRCA, with MOFA+ as the integration framework and Metabric as the external validation set.

## Phases

- [x] **Phase 1: Ingest** - All three modalities stream from S3 via DuckDB, filtered to tumor samples, written to Parquet
- [x] **Phase 2: Merge** - Cohort Parquets pivoted and joined into a single wide feature matrix ready for ML
- [x] **Phase 3: Classify + Deploy** - XGBoost trained on cohort + per-cohort subtype tasks with 5-fold CV, SHAP per-class, panel comparison vs PAM50/Wilkerson/iCluster
- [ ] **Phase 4: BRCA multi-omic on Spark** - Three new modality readers (RPPA, mutations, miRNA) + Spark-native ingest + MOFA+ 6-view integration + consensus clustering + SNF comparator + Metabric external validation, scoped to BRCA

## Phase Details

### Phase 1: Ingest
**Goal**: Three modality Parquets (rna_seq, cnv, methylation) exist on disk for all three cohorts, containing only tumor samples streamed from S3 via DuckDB
**Depends on**: Nothing (first phase)
**Requirements**: INGEST-01, INGEST-02, INGEST-03, INGEST-04
**Success Criteria** (what must be TRUE):
  1. Running any ingest script produces a Parquet file without downloading raw files locally
  2. Each output Parquet contains only rows where the sample type is tumor (no normal tissue contamination)
  3. Methylation ingestor uses DuckDB streaming (not the old pandas-based approach) and writes methylation.parquet
  4. All three modalities produce Parquet files for TCGA-BRCA, TCGA-LUAD, and TCGA-PRAD
**Plans**: 5 plans

Plans:
- [x] 01-01-PLAN.md — Fix gdc_client.py: tumor filter, public S3 bucket, anonymous credentials
- [x] 01-02-PLAN.md — Fix S3 path in ingest_rnaseq.py and ingest_cnv.py (remove modality prefix)
- [x] 01-03-PLAN.md — Rewrite ingest_methylation.py as DuckDB S3 streaming ingestor
- [x] 01-04-PLAN.md — Create run_ingest.py orchestrator (loops 3 cohorts x 3 modalities)
- [x] 01-05-PLAN.md — Rewrite broken test suite; add test_run_ingest.py

### Phase 2: Merge
**Goal**: A single merged_all_cohorts.parquet exists with ~1,000 features × ~1,000 patients, stacking all three cohorts, with a cohort label column
**Depends on**: Phase 1
**Requirements**: MERGE-01, MERGE-02, MERGE-03, MERGE-04
**Success Criteria** (what must be TRUE):
  1. merged_all_cohorts.parquet loads without error and has one row per patient
  2. The feature matrix contains RNA-seq columns (top 500 genes), CNV arm columns, and methylation columns (top 500 probes)
  3. Only patients present in all three modalities appear (inner join enforced)
  4. A cohort column identifies each patient as BRCA, LUAD, or PRAD
**Plans**: 2 plans

Plans:
- [x] 02-01-PLAN.md — Implement src/merge.py: all transform helpers and merge_all_cohorts() orchestrator
- [x] 02-02-PLAN.md — Write tests/test_merge.py (four test classes) and create run_merge.py runner

### Phase 3: Classify + Deploy
**Goal**: XGBoost produces a test accuracy number and SHAP summary plot from the merged feature matrix, and the entire pipeline runs start-to-finish on EC2
**Depends on**: Phase 2
**Requirements**: CLASS-01, CLASS-02, CLASS-03, DEPLOY-01
**Success Criteria** (what must be TRUE):
  1. XGBoost trains on 3-class cancer type labels (BRCA / LUAD / PRAD) without error — done (cohort acc 0.999 ± 0.001)
  2. A test accuracy score prints to stdout after training — done
  3. A SHAP summary plot is saved to disk showing which features drive predictions — done (per-class SHAP for cohort + each subtype task)
  4. The full pipeline (ingest → merge → classify) completes without manual intervention on an EC2 instance — done
**Plans**: shipped via PR #4 (polars ingest), PR #5 (preprocessing/model), PR #6 (5-fold CV + marker force-include + SHAP-vs-panel)
**Outputs**: cohort + BRCA (0.87) + LUAD (0.59) + PRAD (0.93) subtype classifiers, SHAP panel comparison report

### Phase 4: BRCA multi-omic on Spark
**Goal**: Six-modality MOFA+ on BRCA recovers PAM50 subtypes via unsupervised multi-omic factor analysis, run on a real Spark-on-EC2 cluster. Consensus clustering on factor scores + SNF as comparator. External validation on Metabric BRCA. Decision point at end on whether to scale to more cohorts.
**Depends on**: Phase 3 (existing BRCA classification baseline)
**Requirements**: TBD (to be specified by `/gsd-plan-phase 4`)
**Success Criteria** (what must be TRUE):
  1. Spark cluster spins up via setup-spark-cluster.sh and runs an end-to-end Spark job on BRCA methylation
  2. Three new modality readers (RPPA, mutations, miRNA) produce per-patient Parquet for BRCA matching the existing Phase 1 pattern
  3. MOFA+ produces a 6-view factor model on BRCA with active-factor count via ARD pruning
  4. Consensus clustering on MOFA+ factors recovers PAM50 subtypes with ARI > 0.5
  5. SNF run on the same 6 views produces clusters comparable to MOFA+ via convergence analysis
  6. BRCA model and MOFA+ factors generalize to Metabric BRCA cohort (cross-platform harmonization handled)

Sub-phases:
  4a. [x] Spark cluster scaffolding — `spark/setup-spark-cluster.sh` ported from midterm-02. Cluster brought up + torn down without smoke test (rubric demonstration deprioritized; Phase 4 doesn't need the cluster downstream).
  4b. [x] New modality ingest readers — RPPA / miRNA (Polars `scan_csv` batch) + mutations (DuckDB; MAF format quirks). All BRCA-scope. Outputs at `data/TCGA-BRCA/{rppa,mirna,mutations}.parquet`.
  4c. [cut] Spark-native ingest port — would re-do existing Phase 1 parquets in PySpark with no scientific value; cut once 4a's "Spark on EC2" demonstration was deprioritized.
  4d. [x] Multi-omic merge for BRCA — `merge_brca_6view` in `src/merge.py` produces `data/TCGA-BRCA/merged_brca_6view.parquet` (744 patients × 12,612 cols). Re-run on 2026-05-06 after fixing a sparse-feature selection bug in `_top_n_by_variance` (see known-issues.md).
  4e. [x] MOFA+ multi-omic with modality-appropriate likelihoods (Gaussian × 5 + Bernoulli × 1). `run_mofa.py` extended; full run 2026-05-06 in 16 min, 14 active factors after ARD pruning. Factor 1 is multi-omic across all 6 views (basal-vs-luminal axis); factor 2 is proliferation/driver-mutation. Methylation-specific factors 3 and 7.
  4f. [partial] k-means + silhouette + ARI/NMI vs PAM50 done via `analyze_mofa.py`. ARI 0.27–0.31 (below the > 0.5 target — see known-issues.md "PAM50 vs multi-omic comparison ceiling"). Cluster 0 at k=4 is 97% basal-like. **Not yet done**: consensus clustering with bootstrap resampling; SNF comparator.
  4g. [ ] Metabric BRCA external validation — cross-platform harmonization (microarray ↔ RNA-seq). Likely RNA + meth + CNV only; RPPA/miRNA Metabric availability TBD.
  4h. [ ] Phase results writeup + decision: scale to more cohorts (Phase 5) or stop here.

**Out of scope for Phase 4**:
  - Cohort scaling beyond BRCA (deferred to Phase 5 if Phase 4 shows multi-omic value)
  - Survival analysis (clinical modality) — separate optional ML task
  - LUAD / PRAD multi-omic — same reason as cohort scaling
**UI hint**: no

## Progress

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Ingest | 5/5 | Complete | 2026-04-27 |
| 2. Merge | 2/2 | Complete | 2026-04-28 |
| 3. Classify + Deploy | shipped via PRs #4-#6 | Complete | 2026-05-05 |
| 4. BRCA multi-omic on Spark | 5/8 (4a/4b/4d/4e done; 4f partial; 4c cut) | In progress | - |
