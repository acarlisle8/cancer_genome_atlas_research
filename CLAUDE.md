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

### Branch & merge workflow

All phase work happens in Claude Code worktrees on `claude/*` branches — never commit directly to `main` (enforced by `.githooks/pre-commit`).

**After each phase is complete, remind the user to:**
1. Open a PR from the worktree branch into `main`
2. Use **Squash and merge** on GitHub so `main` gets one clean commit per phase
3. Delete the worktree branch after merge

## Tech Stack

- **Data ingestion:** `boto3` or `s3fs` for S3 access
- **Processing:** DuckDB + Polars (Phases 1-3, single-cohort scope), Spark on EC2 cluster (Phase 4+, multi-omic + scale)
- **Output format:** Parquet
- **ML:** XGBoost, SHAP, MOFA+ (multi-omic factor analysis), SNF (similarity network fusion)
- **Compute:** Local single-node for Phases 1-3; multi-node Spark cluster on EC2 for Phase 4+

## Key Constraints

- Cost-conscious: prefer standalone Spark cluster on EC2 over managed EMR
- Pipeline must be reproducible (scriptable, not ad-hoc notebook steps)
- TCGA data is public: `s3://tcga-2-open/` (raw); personal sync at `s3://g23861422-datsbd-s2026/tcga/`
