---
status: partial
phase: 02-merge
source: [02-VERIFICATION.md]
started: 2026-04-27T21:35:00Z
updated: 2026-04-27T21:35:00Z
---

## Current Test

[awaiting human testing]

## Tests

### 1. End-to-end run on EC2 with real Phase 1 Parquets

expected: `python run_merge.py` completes without error; `data/merged_all_cohorts.parquet` is written with ~1,050 columns and ~1,000 rows; cohort column contains {BRCA, LUAD, PRAD}; one row per patient; no OOM error during methylation processing
result: [pending]

## Summary

total: 1
passed: 0
issues: 0
pending: 1
skipped: 0
blocked: 0

## Gaps
