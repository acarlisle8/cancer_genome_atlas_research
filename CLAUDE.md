# TCGA Cancer Genomics Pipeline

## Project Overview

TCGA ETL pipeline ingesting from public AWS S3, aggregating per-patient RNA-seq, CNV, and methylation files into Parquet, joining on shared patient identifiers across multiple cancer cohorts, and running XGBoost classifiers with SHAP analysis.

**Core deliverable:** Reproducible pipeline. ML results prove it works.

## GSD Workflow

This project uses GSD for structured execution. Planning artifacts live in `.planning/`.

- **Roadmap:** `.planning/ROADMAP.md`
- **Requirements:** `.planning/REQUIREMENTS.md`
- **State:** `.planning/STATE.md`

### Commands

```
/gsd-plan-phase 1        # Plan the next phase
/gsd-execute-phase 1     # Execute a planned phase
/gsd-progress            # Check project status
```

## Tech Stack

- **Data ingestion:** `boto3` or `s3fs` for S3 access
- **Processing:** DuckDB + Polars (primary), EMR/Spark (fallback for scale)
- **Output format:** Parquet
- **ML:** XGBoost, SHAP
- **Compute:** Local-first; AWS (personal account) for scale if needed

## Key Constraints

- Timeline: ~2-3 weeks
- Cost-conscious: avoid EMR unless local compute is genuinely insufficient
- Pipeline must be reproducible (scriptable, not ad-hoc notebook steps)
- TCGA data is public: `s3://tcga-2-open/`
