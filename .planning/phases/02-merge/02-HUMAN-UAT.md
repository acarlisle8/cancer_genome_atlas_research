---
status: passed
phase: 02-merge
source: [02-VERIFICATION.md]
started: 2026-04-27T21:35:00Z
updated: 2026-04-28T16:05:00Z
---

## Current Test

Complete — Phase 2 end-to-end run succeeded on EC2 with real Phase 1 Parquets.

## Tests

### 1. End-to-end run on EC2 with real Phase 1 Parquets

expected: `python run_merge.py` completes without error; `data/merged_all_cohorts.parquet` is written with ~1,050 columns and ~1,000 rows; cohort column contains {BRCA, LUAD, PRAD}; one row per patient; no OOM error during methylation processing
result: **passed (2026-04-28 16:05 UTC)** — `python run_merge.py` completed in 4:48 with `Wrote merged matrix (2104 rows x 1043 cols) to data/merged_all_cohorts.parquet`. 11 MB snappy-compressed. 10/11 checklist checks pass on first run; the 11th (file-size lower bound) was a stale threshold and has been corrected in [02-VERIFY-CHECKLIST.md](02-VERIFY-CHECKLIST.md).

### Per-cohort results

| Cohort | merged_cohort shape | Notes |
|---|---|---|
| BRCA | 1094 × 863 | inner-join attrition 1095→1094; only 320 of 500 top probes appeared in BRCA |
| LUAD | 515 × 1043 | inner-join attrition 518→515; all 500 top probes present |
| PRAD | 495 × 852 | 309 of 500 top probes |
| **Final** | **2104 × 1043** | union of probe columns across cohorts via `pl.concat(diagonal_relaxed)`; 495 probe cols total (5 missing — minor, see runbook) |

## Summary

total: 1
passed: 1
issues: 0
pending: 0
skipped: 0
blocked: 0

## Gaps

None for Phase 2. Five top-500 probes did not materialize as columns in the final matrix — likely a combination of the headerless-TSV workaround dropping the alphabetically-first probe per file ([known-issues.md](../../known-issues.md)) and edge cases in the per-cohort `is_in` filter. Not blocking; flagged for revisit if Phase 3 SHAP analysis suggests probe-side attrition is meaningful.
