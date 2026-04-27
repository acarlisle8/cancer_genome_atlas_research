---
phase: 02-merge
plan: 02
subsystem: test
tags: [polars, parquet, pivot, join, cnv, methylation, rnaseq, unittest, tempfile]

# Dependency graph
requires:
  - phase: 02-merge
    plan: 01
    provides: "src/merge.py with _pivot_rnaseq, _pivot_cnv, _pivot_methylation, merge_all_cohorts"
provides:
  - "tests/test_merge.py with TestRnaseqPivot, TestCnvPivot, TestMethylationPivot, TestMergeAllCohorts"
  - "run_merge.py merge pipeline entry point"
affects: [03-ml]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "deferred import inside each test method: from src.merge import X"
    - "self-contained test with tempfile.TemporaryDirectory + pl.write_parquet synthetic data"
    - "_make_data_dir module-level helper writing 3 cohorts x 3 modalities of synthetic Parquet"
    - "pl.DataFrame(rows, schema={...}) for typed null columns in test data"

key-files:
  created:
    - tests/test_merge.py
    - run_merge.py
  modified:
    - src/merge.py

key-decisions:
  - "All pl.read_parquet() calls inside TemporaryDirectory with-block to avoid FileNotFoundError after cleanup"
  - "Distinct patient_id prefix per cohort in _make_data_dir (TCGA-BR-, TCGA-LU-, TCGA-PR-) to guarantee no cross-cohort collisions"
  - "Both p-arm and q-arm CNV segments in _make_data_dir so _pivot_cnv produces non-empty output for each cohort"
  - "Fixed _pivot_cnv arm column sort in src/merge.py (Rule 1 bug) to make pl.concat schema-consistent"

# Metrics
duration: 15min
completed: 2026-04-27
---

# Phase 02 Plan 02: Tests and Runner Summary

**14-test suite covering all four MERGE-XX requirements with synthetic Parquet in TemporaryDirectory, plus run_merge.py runner mirroring run_pipeline.py, and a Rule 1 bug fix in _pivot_cnv for non-deterministic arm column ordering**

## Performance

- **Duration:** ~15 min
- **Started:** 2026-04-27T21:05:00Z
- **Completed:** 2026-04-27T21:21:00Z
- **Tasks:** 2
- **Files modified:** 3 (tests/test_merge.py created, run_merge.py created, src/merge.py fixed)

## Accomplishments

- `tests/test_merge.py`: 14 tests across 4 classes, all passing
  - `TestRnaseqPivot` (3 tests): shape, selected-gene columns, excluded-gene absence
  - `TestCnvPivot` (4 tests): arm naming (1p/1q not bare p/q), midpoint assignment, chrM filtered, one row per patient
  - `TestMethylationPivot` (3 tests): shape, null beta_value preserved, excluded probes absent
  - `TestMergeAllCohorts` (4 tests): output file exists, one row per patient, cohort values {BRCA,LUAD,PRAD}, inner-join excludes missing patients
- `run_merge.py`: minimal runner exactly mirroring `run_pipeline.py` — no argparse, no `__main__` guard, direct call to `merge_all_cohorts(DATA_DIR, OUTPUT_DIR)`
- `src/merge.py` (bug fix): added deterministic arm column sort in `_pivot_cnv` to prevent `ShapeError` in `pl.concat`

## Task Commits

1. **Task 1: Write tests/test_merge.py + fix _pivot_cnv** - `07a3f4b` (feat)
2. **Task 2: Create run_merge.py** - `e629019` (feat)

## Files Created/Modified

- `/Users/aidancarlisle/Documents/cancer_genome_atlas_research/tests/test_merge.py` — 14 tests, 4 classes, `_make_data_dir` helper, deferred imports, no mocking; 360 lines
- `/Users/aidancarlisle/Documents/cancer_genome_atlas_research/run_merge.py` — 8-line runner
- `/Users/aidancarlisle/Documents/cancer_genome_atlas_research/src/merge.py` — `_pivot_cnv` arm column sort fix (5-line change)

## Decisions Made

- **All read_parquet inside with-block:** `pl.read_parquet(out)` in `TestMergeAllCohorts` tests is inside the `with tempfile.TemporaryDirectory()` block — file is on disk while it's being read.
- **Distinct patient_id prefixes in _make_data_dir:** Cohort prefixes `BR`/`LU`/`PR` ensure no cross-cohort patient_id collision, satisfying `merge_all_cohorts`'s post-condition assert on uniqueness.
- **Two-segment CNV rows per patient:** Each patient in `_make_data_dir` gets one p-arm and one q-arm segment so `_pivot_cnv` produces a non-empty wide frame for all three cohorts and `pl.concat` works.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed non-deterministic arm column order in _pivot_cnv**
- **Found during:** Task 1 execution — `TestMergeAllCohorts.test_cohort_column_values` failed with `polars.exceptions.ShapeError: unable to vstack, column names don't match: "1q" and "1p"`
- **Issue:** `_pivot_cnv` uses `group_by(["patient_id", "chr_arm"])` which has non-deterministic output order. `pl.pivot` inherits that order for the new columns. When BRCA produces `[patient_id, 1q, 1p]` and LUAD produces `[patient_id, 1p, 1q]`, `pl.concat(how='vertical')` raises `ShapeError`.
- **Fix:** After the `pivot`, select columns with arm columns sorted alphabetically: `arm_cols = sorted(c for c in wide.columns if c != "patient_id")` → `return wide.select(["patient_id"] + arm_cols)`
- **Files modified:** src/merge.py (lines 133-142, 5 lines added/changed)
- **Commit:** `07a3f4b`

---

**Total deviations:** 1 auto-fixed (Rule 1 — non-deterministic column order bug)
**Impact on plan:** Acceptance criteria satisfied; fix is also necessary for correctness in production use (real cohorts would have the same ordering instability).

## Issues Encountered

None beyond the auto-fixed Rule 1 bug above.

## User Setup Required

None — no external services, no environment variables.

## Next Phase Readiness

- Phase 2 is fully complete: `src/merge.py` and `tests/test_merge.py` both verified passing
- `run_merge.py` ready to run once Phase 1 Parquets exist on EC2
- Phase 3 ML can import from `src/merge.py` via `merge_all_cohorts()` and assume `merged_all_cohorts.parquet` at `data/merged_all_cohorts.parquet`

## Known Stubs

None — all test assertions use real computed values from synthetic Parquet; `run_merge.py` makes a real call to `merge_all_cohorts`.

---
*Phase: 02-merge*
*Completed: 2026-04-27*

## Self-Check: PASSED

- `tests/test_merge.py` exists: FOUND
- `run_merge.py` exists: FOUND
- `02-02-SUMMARY.md` exists: FOUND
- Commit `07a3f4b` exists: FOUND
- Commit `e629019` exists: FOUND
- No unexpected file deletions in either commit
- `python -m pytest tests/test_merge.py -x -q` exits 0 (14 passed)
