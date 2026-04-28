# Roadmap: TCGA Cancer Type Classifier

## Overview

Three phases mirror the three-step pipeline: get the data out of S3 and into Parquet, join and pivot it into a feature matrix, then train the classifier and run the whole thing on EC2 to produce the course deliverable. Each phase delivers one complete, independently verifiable capability.

## Phases

- [x] **Phase 1: Ingest** - All three modalities stream from S3 via DuckDB, filtered to tumor samples, written to Parquet
- [x] **Phase 2: Merge** - Cohort Parquets pivoted and joined into a single wide feature matrix ready for ML
- [ ] **Phase 3: Classify + Deploy** - XGBoost trains on the merged matrix, SHAP plots produced, full pipeline runs end-to-end on EC2

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
  1. XGBoost trains on 3-class cancer type labels (BRCA / LUAD / PRAD) without error
  2. A test accuracy score prints to stdout after training
  3. A SHAP summary plot is saved to disk showing which features drive predictions
  4. The full pipeline (ingest → merge → classify) completes without manual intervention on an EC2 instance
**Plans**: TBD
**UI hint**: no

## Progress

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Ingest | 5/5 | Complete | 2026-04-27 |
| 2. Merge | 2/2 | Complete | 2026-04-28 |
| 3. Classify + Deploy | 0/? | Not started | - |
