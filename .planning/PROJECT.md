# TCGA Cancer Genomics Pipeline

## What This Is

A reproducible ETL pipeline that ingests TCGA data from AWS public S3, aggregates per-patient files into Parquet, and joins RNA-seq, copy number variation, and DNA methylation data across shared patient identifiers — covering several cancer cohorts. The integrated dataset feeds XGBoost classifiers for cancer subtype prediction and pan-cancer type classification, with SHAP analysis to identify driving genes.

Built for DATS 6540 Big Data Analytics (Group 2: Aidan Carlisle and Zachary Cardell).

## Core Value

A working, reproducible pipeline that produces a clean integrated Parquet dataset — the ML results prove it works, but the pipeline is the deliverable.

## Requirements

### Validated

(None yet — ship to validate)

### Active

- [ ] Ingest TCGA data from public AWS S3 bucket
- [ ] Aggregate per-patient files into consolidated Parquet files per modality
- [ ] Join RNA-seq, CNV, and methylation data on shared patient identifiers
- [ ] Cover several cancer cohorts (5-10 TCGA cancer types)
- [ ] Train XGBoost classifier for molecular subtype prediction
- [ ] Train XGBoost classifier for cancer type classification across cohorts
- [ ] SHAP analysis to identify gene-level feature importance
- [ ] Compare results against published classifiers for validation

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
| DuckDB + Polars as primary compute | Fast local development, handles medium cohort sizes without EMR cost | — Pending |
| EMR/Spark as fallback only | Cost-conscious; personal AWS account | — Pending |
| Parquet as intermediate format | Columnar, efficient for downstream ML feature extraction | — Pending |
| XGBoost for classification | Standard baseline; results just need to be produced, not novel | — Pending |
| SHAP for interpretability | Satisfies requirement to identify driving genes; published comparison point | — Pending |

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
