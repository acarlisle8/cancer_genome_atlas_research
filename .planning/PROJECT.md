# TCGA Cancer Genomics Pipeline

## What This Is

A reproducible ETL pipeline that ingests TCGA data from AWS public S3, aggregates per-patient files into Parquet, and joins RNA-seq, copy number variation, and DNA methylation data across shared patient identifiers — covering several cancer cohorts. The integrated dataset feeds XGBoost classifiers for cancer subtype prediction and pan-cancer type classification, with SHAP analysis to identify driving genes.

Built for DATS 6540 Big Data Analytics (Group 2: Aidan Carlisle and Zachary Cardell).

## Core Value

A working, reproducible pipeline that produces a clean integrated Parquet dataset — the ML results prove it works, but the pipeline is the deliverable.

## Requirements

### Validated

- [x] Ingest TCGA data from public AWS S3 bucket (Phase 1)
- [x] Aggregate per-patient files into consolidated Parquet files per modality (Phase 1)
- [x] Join RNA-seq, CNV, and methylation data on shared patient identifiers (Phase 2)
- [x] Train XGBoost classifier for molecular subtype prediction (Phase 3 — BRCA 0.87, PRAD 0.93, LUAD 0.59 with 5-fold CV)
- [x] Train XGBoost classifier for cancer type classification across cohorts (Phase 3 — 0.999 ± 0.001 on BRCA/LUAD/PRAD)
- [x] SHAP analysis to identify gene-level feature importance (Phase 3 — per-class SHAP for cohort + each subtype task)
- [x] Compare results against published classifiers — internal panel comparison (Phase 3 — PAM50/Wilkerson/iCluster panel hits in top SHAP features)
- [x] Multi-omic integration via MOFA+ (RNA + methylation + CNV + RPPA + mutations + miRNA) — Phase 4 (BRCA 6-view factor model, 14 active factors, factor 1 multi-omic across all 6 views)
- [x] Spark-on-EC2 cluster scaffolding — Phase 4a (cluster scaffolded; smoke test deprioritized; cluster not used for downstream Phase 4 work)

### Active

- [ ] Cover several cancer cohorts (5-10 TCGA cancer types) — currently 3 of 33; deferred to Phase 5 pending Phase 4 results
- [ ] External validation against independent cohort (Metabric for BRCA) — Phase 4g
- [ ] Unsupervised subtype recovery via consensus clustering on MOFA+ factors + SNF comparator — Phase 4f (k-means + silhouette done with ARI 0.27–0.31 vs PAM50; consensus clustering and SNF still pending)

### Out of Scope

- Novel ML methodology — results just need to be produced and reproducible
- Real-time or streaming ingestion — batch ETL only
- Web UI or dashboard — analysis outputs are enough
- Mobile or external user-facing interface

## Context

- Dataset: The Cancer Genome Atlas (TCGA) via GDC, hosted on public AWS S3 (`s3://tcga-2-open/`)
- Data modalities: RNA-seq (gene expression), copy number variation, DNA methylation
- Per-patient file format: inconsistent across patients and cancer types — aggregation is non-trivial
- Compute strategy: develop locally with DuckDB + Polars (fits for single-cohort work), scale to AWS EMR/Spark only if full multi-cohort volume demands it (personal account — cost-conscious)
- Team of 2, ~2-3 week timeline
- Validated labels for subtype classification exist in TCGA clinical metadata

## Constraints

- **Timeline**: ~2-3 weeks to completion — must scope tightly
- **Cost**: Personal AWS account — EMR should only be used if local compute is genuinely insufficient
- **Data format**: TCGA files are inconsistent per-patient; aggregation logic must handle missing files and format variations
- **Reproducibility**: Pipeline must be reproducible (scripted, not ad-hoc notebook steps)

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| DuckDB + Polars as primary compute (Phases 1-3) | Fast local development, handles 3-cohort sizes without cluster cost | Validated — shipped Phases 1-3 |
| Spark on EC2 cluster (Phase 4+) | Multi-omic at BRCA scale + course rubric requires distributed compute on cluster | Deprioritized — cluster scaffolded in 4a but rubric demonstration was cut as low-value; downstream 4d–4f run on local DuckDB+Polars+mofapy2 |
| Lazy Polars (`pl.scan_parquet`) for any pass over the methylation/RNA-seq parquets | 7.6 GB EC2 instance with no swap; eager `pl.read_parquet` on the 444M-row methylation parquet OOM-kills consistently. Three OOMs paid 2026-05-05. | Locked 2026-05-06 — see [known-issues.md](known-issues.md) and [feedback memory](.../memory/feedback_dats6450_lazy_polars.md) |
| Per-view max-missing-rate dict in `run_mofa.py` (methylation 0.40, others 0.20) | Merged 6-view's mean meth missingness is ~27% (vs ~18% source baseline) due to 6-way inner-join shifting kept patients toward worse meth coverage; uniform 0.20 drops 92% of meth probes | Active — Phase 4e |
| `MIN_COVERAGE_FRAC = 0.80` in `_top_n_by_variance` | Variance ranking on probes/genes preferentially admitted 100%-missing-but-bimodal probes; HAVING clause filters before ORDER BY | Active — Phase 4d (re-run on 2026-05-06) |
| Parquet as intermediate format | Columnar, efficient for downstream ML feature extraction | Validated |
| XGBoost for supervised classification | Standard baseline; results just need to be produced, not novel | Validated |
| SHAP for interpretability | Satisfies requirement to identify driving genes; published comparison point | Validated |
| MOFA+ as primary multi-omic integration method | Multi-likelihood (Gaussian + Bernoulli), handles missing data, ARD sparsity | Active — Phase 4 |
| SNF as multi-omic comparator | Network-based paradigm distinct from MOFA+'s factor-based — convergence across paradigms strengthens unsupervised story | Active — Phase 4 |
| Consensus clustering on MOFA+ factor scores (replaces plain k-means) | Stability assessment via resampling; matches published TCGA subtype methodology | Active — Phase 4 |
| BRCA-deep before cohort-wide | Test multi-omic hypothesis on best-instrumented cohort (Metabric exists) before scaling effort | Active — Phase 4 |

## Evolution

This document evolves at phase transitions and milestone boundaries.

**After each phase transition** (via `/gsd-transition`):
1. Requirements invalidated? → Move to Out of Scope with reason
2. Requirements validated? → Move to Validated with phase reference
3. New requirements emerged? → Add to Active
4. Decisions to log? → Add to Key Decisions
5. "What This Is" still accurate? → Update if drifted

**After each milestone** (via `/gsd-complete-milestone`):
1. Full review of all sections
2. Core Value check — still the right priority?
3. Audit Out of Scope — reasons still valid?
4. Update Context with current state

---
*Last updated: 2026-04-18 after initialization*
