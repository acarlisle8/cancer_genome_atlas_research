---
phase: 01-data-ingestion-modality-aggregation
plan: 03
subsystem: ingest-cnv-methylation
status: complete
tags: [ingest, cnv, methylation, parquet, polars]
key-files:
  created:
    - src/ingest_cnv.py
    - src/ingest_methylation.py
    - tests/test_ingest_cnv.py
    - tests/test_ingest_methylation.py
metrics:
  tasks_completed: 2
  tests_added: 18
  tests_passing: 18
---

## What Was Built

Implemented `src/ingest_cnv.py` and `src/ingest_methylation.py` — the CNV and methylation modality aggregators. Both follow the same fetch → download → parse → Parquet pattern established by ingest_rnaseq.py.

### src/ingest_cnv.py
- `parse_cnv_seg(path, patient_id)` — reads masked copy number segment TSV; renames `Chromosome`→`chromosome`, `Start`→`start`, `End`→`end`, `Segment_Mean`→`copy_number`; drops extra columns (GDC_Aliquot, Num_Probes); returns `[patient_id, chromosome, start, end, copy_number]` DataFrame with dtypes `[Utf8, Utf8, Int64, Int64, Float64]`
- `ingest_cnv(output_dir, raw_dir, project_id)` — calls `fetch_manifest` with `data_category="Copy Number Variation"`, `data_type="Masked Copy Number Segment"`; downloads to `raw_dir/cnv/`; writes `output_dir/cnv.parquet` (snappy); writes `errors_cnv.csv` on failures (D-06)

### src/ingest_methylation.py
- `parse_methylation_betas(path, patient_id)` — reads headerless 2-column TSV with `has_header=False, new_columns=["probe_id", "beta_value"]`; casts beta_value to Float64; returns `[patient_id, probe_id, beta_value]` DataFrame
- `ingest_methylation(output_dir, raw_dir, project_id)` — calls `fetch_manifest` with `data_category="DNA Methylation"`, `data_type="Methylation Beta Value"`; downloads; writes `output_dir/methylation.parquet` (snappy); writes `errors_methylation.csv` on failures (D-06)

## Commits

| Hash | Phase | Message |
|------|-------|---------|
| 1e3e26c | test | test(01-03): add failing tests for ingest_cnv (RED) |
| d8021d8 | feat | feat(01-03): implement ingest_cnv.py — CNV segment aggregator (GREEN) |
| ff9c374 | test | test(01-03): add failing tests for ingest_methylation (RED) |
| 996f7c0 | feat | feat(01-03): implement ingest_methylation.py — methylation beta value aggregator (GREEN) |

## Deviations

None — all D-01 through D-07 context decisions honored. CNV column renames exact as specified. Methylation headerless read exact as specified.

## Self-Check: PASSED

- `grep "def parse_cnv_seg" src/ingest_cnv.py` ✓
- `grep "def ingest_cnv" src/ingest_cnv.py` ✓
- `grep "Segment_Mean" src/ingest_cnv.py` ✓
- `grep "errors_cnv.csv" src/ingest_cnv.py` ✓
- `grep "def parse_methylation_betas" src/ingest_methylation.py` ✓
- `grep "has_header=False" src/ingest_methylation.py` ✓
- `grep "errors_methylation.csv" src/ingest_methylation.py` ✓
- `python -m pytest tests/test_ingest_cnv.py tests/test_ingest_methylation.py -q` → 18 passed ✓
