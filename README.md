# DATS 6450 Final Project: Multi-omic TCGA-BRCA Subtyping at Scale

> **Question:** On a six-modality TCGA-BRCA dataset (~744 patients), do
> supervised and unsupervised analyses recover the same molecular-subtype
> structure?

End-to-end ETL and modeling pipeline ingesting TCGA from public AWS S3,
joining six modalities (RNA-seq, methylation, CNV, RPPA, miRNA, mutations)
on shared patient identifiers, and running both supervised XGBoost on
PAM50 labels and unsupervised MOFA+ multi-omic factor analysis.

## Deliverables

- [final_report.qmd](final_report.qmd): source (Quarto)
- [final_report.pdf](final_report.pdf): rendered PDF (12 pages, ~9 main + 3 appendix)
- [final_report.html](final_report.html): rendered self-contained HTML
- [references.bib](references.bib): bibliography

## Repo layout

```
.
├── final_report.qmd / final_report.pdf / final_report.html   ← deliverables
├── references.bib                                             ← bibliography
├── README.md                                                  ← this file
├── pyproject.toml / uv.lock / requirements.txt                ← dependencies
├── code/                  ← entry-point scripts (run_*.py, analyze_*.py, etc.)
├── src/                   ← library modules (ingest_*, merge.py, preprocess.py, ...)
├── tests/                 ← pytest suite
├── spark/                 ← Spark-on-EC2 cluster scaffolding
├── data/                  ← parquet outputs + figures (mostly gitignored)
└── .planning/             ← planning artifacts, runbook, session logs
```

## Reproducing the results

Pre-requisite: Python 3.11+ and [uv](https://docs.astral.sh/uv/). Then:

```bash
uv sync
```

| Stage              | Command                                                | Wall time |
|--------------------|--------------------------------------------------------|-----------|
| Ingest             | `uv run python code/run_ingest.py`                     | ~30–60 min |
| Merge (6-view BRCA)| `uv run python code/run_merge_6view.py`                | ~8 min    |
| Supervised + SHAP  | `uv run python code/run_classification_pipeline.py`    | ~10 min   |
| MOFA+ training     | `uv run python code/run_mofa.py --cohort BRCA`         | ~16 min   |
| MOFA+ analysis     | `uv run python code/analyze_mofa.py --model-dir data/mofa_BRCA --known-labels data/audit/subtype_label_audit.parquet --label-col hoadley_subtype_selected --cancer-type BRCA` | <1 min |
| Render report      | `quarto render final_report.qmd`                       | <30 s     |

Detailed runbook: [.planning/RUNBOOK.md](.planning/RUNBOOK.md).

## Compute

Single AWS EC2 `t3` (2 vCPU, 7.6 GB RAM, no swap). The 444 M-row
methylation parquet (~5.5 GB compressed, ~30 GB decompressed) exceeds
memory if loaded eagerly, so the pipeline uses DuckDB streaming and
Polars lazy evaluation throughout. A multi-node Spark cluster is
scaffolded under [spark/](spark/) but not on the critical path for the
BRCA results.

## Authors

Aidan Carlisle, Zachary Cardell. DATS 6450, Spring 2026.
