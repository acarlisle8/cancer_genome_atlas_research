---
phase: 01-data-ingestion-modality-aggregation
plan: "02"
subsystem: data-ingestion
tags: [rnaseq, gdc-api, s3, parquet, polars, tdd, long-format]
dependency_graph:
  requires:
    - src/gdc_client.py (fetch_manifest, download_file) — from 01-01
    - src/utils.py (get_logger) — from 01-01
  provides:
    - src/ingest_rnaseq.py (parse_rnaseq_tsv, ingest_rnaseq)
    - data/brca/rna_seq.parquet (produced at runtime)
  affects:
    - Phase 2 join logic consumes rna_seq.parquet on patient_id key
tech_stack:
  added:
    - polars (read_csv, filter, select, cast, write_parquet)
  patterns:
    - TDD (RED/GREEN per plan)
    - Long-format schema: patient_id × gene_id rows (D-01)
    - D-06 skip-and-log: failed patients written to errors_rnaseq.csv
    - D-07 resume: download_file skips existing files
    - T-02-01 mitigation: explicit column selection
    - T-02-02 mitigation: per-patient try/except
key_files:
  created:
    - src/ingest_rnaseq.py
    - tests/test_ingest_rnaseq.py
  modified: []
decisions:
  - key: N_ row filter applied at parse level
    detail: "polars filter(~col('gene_id').str.starts_with('N_')) removes unmapped/multimapping/noFeature/ambiguous summary rows before any downstream logic"
  - key: fpkm_unstranded chosen over tpm_unstranded
    detail: "FPKM retained as-is per CONTEXT.md Claude's discretion — use normalized values as provided; can be revisited if ML results poor in Phase 3"
  - key: errors_rnaseq.csv only created on failure
    detail: "File is NOT written when all parses succeed — reduces clutter in clean runs"
metrics:
  duration_seconds: 180
  completed_date: "2026-04-19"
  tasks_completed: 1
  files_created: 2
---

# Phase 1 Plan 02: RNA-seq Ingestion and Aggregation Summary

**One-liner:** RNA-seq modality aggregator using GDC manifest fetch + anonymous S3 download + polars TSV parsing into long-format Parquet [patient_id, gene_id, fpkm_unstranded] with D-06 skip-and-log error CSV and snappy compression.

## Tasks Completed

| Task | Name | Type | Commits | Status |
|------|------|------|---------|--------|
| 1 | Implement src/ingest_rnaseq.py — parse single TSV to long-format DataFrame | auto/tdd | 72fa8a1 (RED), 3304ca8 (GREEN) | Done |

## What Was Built

**src/ingest_rnaseq.py** — Two public functions:

- `parse_rnaseq_tsv(path, patient_id)`: Reads an augmented_star_gene_counts.tsv with polars, filters rows where `gene_id.starts_with("N_")` (summary stats), selects `[gene_id, fpkm_unstranded]`, casts to `[Utf8, Float64]`, adds literal `patient_id` column. Returns 3-column DataFrame: `[patient_id, gene_id, fpkm_unstranded]`.

- `ingest_rnaseq(output_dir, raw_dir, project_id="TCGA-BRCA")`: Calls `fetch_manifest("TCGA-BRCA", "Transcriptome Profiling", "Gene Expression Quantification")`, downloads each entry to `raw_dir/rnaseq/`, parses each TSV via `parse_rnaseq_tsv`. Failed patients are caught, logged, and recorded to `errors_rnaseq.csv` (D-06). Successful DataFrames are concatenated and written to `output_dir/rna_seq.parquet` with snappy compression. Returns the parquet path.

**tests/test_ingest_rnaseq.py** — 13 unit tests:
- 7 tests for `parse_rnaseq_tsv`: column names, dtypes, N_ filtering, patient_id literal, fpkm values, only 3 columns
- 6 tests for `ingest_rnaseq`: writes Parquet, schema, row count, manifest call args, error CSV on failure, no error CSV when clean

## Verification

```
python -m pytest tests/test_ingest_rnaseq.py -v   # 13 passed
python -m pytest tests/ -v                        # 34 passed (full suite)
python -c "from src.ingest_rnaseq import ingest_rnaseq, parse_rnaseq_tsv; print('OK')"
```

## TDD Gate Compliance

- RED gate: commit 72fa8a1 (`test(01-02): add failing tests for ingest_rnaseq`)
- GREEN gate: commit 3304ca8 (`feat(01-02): implement ingest_rnaseq.py — RNA-seq modality aggregator`)

Both TDD gates satisfied.

## Deviations from Plan

None — plan executed exactly as written.

## Threat Surface Scan

All threats from the plan's threat model addressed:

| Threat ID | Mitigation Applied |
|-----------|-------------------|
| T-02-01 | `pl.read_csv` with `separator="\t"`, then explicit `.select(["gene_id", "fpkm_unstranded"])` — unknown columns silently ignored |
| T-02-02 | `try/except` wraps `download_file` + `parse_rnaseq_tsv` per patient; errors logged and written to CSV, pipeline continues |
| T-02-03 | Accepted — errors CSV in local output dir under user control |
| T-02-04 | Accepted — output path is `output_dir / "rna_seq.parquet"` with no interpolation |

No new threat surface introduced beyond the plan's model.

## Known Stubs

None.

## Self-Check: PASSED

- `src/ingest_rnaseq.py` exists
- `tests/test_ingest_rnaseq.py` exists
- RED commit 72fa8a1 verified in git log
- GREEN commit 3304ca8 verified in git log
- All 13 plan-02 tests pass; all 34 total tests pass
