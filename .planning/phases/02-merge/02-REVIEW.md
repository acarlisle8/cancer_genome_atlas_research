---
phase: 02-merge
reviewed: 2026-04-27T21:25:53Z
depth: standard
files_reviewed: 3
files_reviewed_list:
  - src/merge.py
  - tests/test_merge.py
  - run_merge.py
findings:
  critical: 0
  warning: 3
  info: 3
  total: 6
status: issues_found
---

# Phase 02: Code Review Report

**Reviewed:** 2026-04-27T21:25:53Z
**Depth:** standard
**Files Reviewed:** 3
**Status:** issues_found

## Summary

Three files were reviewed: the core merge module (`src/merge.py`), its test suite (`tests/test_merge.py`), and the CLI runner (`run_merge.py`). The overall structure is solid — the pivot/join orchestration is well-documented, the inner-join semantics are correct, and the lazy scan for methylation is a good call on memory grounds.

Three warnings require fixes before this can run reliably on real TCGA data. The most significant is a runtime crash risk in `pl.concat` when CNV arm column sets differ across cohorts (e.g., one cohort has no chrX segments). The second is a logic bug in `_top_n_by_variance` where null variances sort to the top under Polars' default behavior, incorrectly selecting uninformative features. The third is an `assert` statement that is stripped by Python's optimizer flag and cannot be relied on as a data integrity guard. Three informational items cover missing test coverage, missing `__main__` guard, and a hardcoded path.

## Warnings

### WR-01: `pl.concat(..., how="vertical")` crashes when CNV arm columns differ across cohorts

**File:** `src/merge.py:240`
**Issue:** `_pivot_cnv` produces arm columns only for chromosomes that actually appear in a given cohort's segment file. If LUAD has no chrX segments, its wide frame has no `Xp`/`Xq` columns; if BRCA does, the two frames have different schemas. `pl.concat(..., how="vertical")` requires identical schemas and raises `unable to append to a DataFrame of width N with a DataFrame of width M`. This will crash in production on any cohort whose CNV data lacks coverage on one or more chromosomes. RNA-seq and methylation columns are controlled by the global top-N feature list and are not affected.

**Fix:** Use `how="diagonal_relaxed"` for the final concat. This fills missing arm columns with null for cohorts that lack them, which is the correct domain semantics (no segment data = missing, not zero).

```python
# line 240
final = pl.concat(cohort_frames, how="diagonal_relaxed")
```

---

### WR-02: `_top_n_by_variance` selects null-variance features first due to Polars null sort order

**File:** `src/merge.py:52`
**Issue:** In Polars, `.sort("variance", descending=True)` on a `LazyFrame` places nulls **first** by default (Polars treats null as the largest value in a descending sort). A gene or probe that appears in only one row across all cohorts yields `var() = null`. Those null-variance features sort to the top of the list and get selected into `top_genes` / `top_probes`, displacing genuinely high-variance features. This is a logic bug: a feature with undefined variance should rank last, not first. The impact is limited in full TCGA data (most genes appear in thousands of samples) but can corrupt feature selection in partial datasets or new cohorts with limited samples.

**Fix:** Add `nulls_last=True` to the sort:

```python
# line 52
.sort("variance", descending=True, nulls_last=True)
```

---

### WR-03: Duplicate patient_id guard uses `assert` (stripped by `python -O`)

**File:** `src/merge.py:247`
**Issue:** The post-condition check `assert final["patient_id"].n_unique() == len(final), "..."` is a data integrity guard against silent data corruption (duplicate rows). Python's `-O` (optimize) flag and tools like `pyinstaller` strip all `assert` statements. If this pipeline is ever run optimized, duplicate rows would pass through silently to the downstream ML model, potentially inflating apparent accuracy. The docstring even documents this as raising `AssertionError`, which is inconsistent with the `RuntimeError` raised three lines earlier for the empty-matrix case.

**Fix:** Replace with an explicit `RuntimeError` using the same pattern as the empty-matrix check:

```python
# line 247
n_unique = final["patient_id"].n_unique()
if n_unique != len(final):
    raise RuntimeError(
        f"Duplicate patient_id rows in merged output: "
        f"{len(final)} rows but only {n_unique} unique patient IDs."
    )
```

Also update the docstring `Raises:` section to say `RuntimeError` (not `AssertionError`) for this condition.

---

## Info

### IN-01: `_top_n_by_variance` has no direct unit tests

**File:** `tests/test_merge.py`
**Issue:** The test suite covers `_pivot_rnaseq`, `_pivot_cnv`, `_pivot_methylation`, and `merge_all_cohorts`, but `_top_n_by_variance` has no dedicated test class. There is no test verifying (a) that the top-N genes are actually the highest-variance ones, (b) that `nulls_last` behavior is correct (before the fix for WR-02), or (c) that fewer-than-N features in the data returns all available features without error. The function is the entry point for global feature selection, so a bug there affects every downstream column.

**Fix:** Add a `TestTopNByVariance` class with at least:
- a test that verifies the returned IDs are ordered by descending variance
- a test where the dataset has fewer than N features (returns all of them, no error)
- a test where one feature has a single data point (null variance sorts last after WR-02 fix)

---

### IN-02: `run_merge.py` lacks a `__main__` guard and error handling

**File:** `run_merge.py:1-8`
**Issue:** The script executes `merge_all_cohorts(...)` at module import time — there is no `if __name__ == "__main__":` guard. This means importing `run_merge` anywhere (e.g., in a test, notebook, or REPL) triggers a full pipeline run against `data/` in the current working directory. Additionally there is no `try/except` around the call, so a missing `data/` directory produces an unformatted `FileNotFoundError` with no actionable message for the user.

**Fix:**
```python
from pathlib import Path
from src.merge import merge_all_cohorts

DATA_DIR = Path("data")
OUTPUT_DIR = Path("data")

if __name__ == "__main__":
    out = merge_all_cohorts(DATA_DIR, OUTPUT_DIR)
    print(f"Merged matrix written to {out}")
```

---

### IN-03: Hardcoded relative path in `run_merge.py`

**File:** `run_merge.py:4-5`
**Issue:** `DATA_DIR = Path("data")` and `OUTPUT_DIR = Path("data")` are hardcoded relative paths. The script must be run from the project root for these to resolve correctly; running from any other directory silently targets the wrong data directory (or creates a new one). There is no validation that `DATA_DIR` exists before calling `merge_all_cohorts`.

**Fix:** Accept paths as CLI arguments (or at minimum use `Path(__file__).parent / "data"` to anchor relative to the script location):

```python
import argparse
from pathlib import Path
from src.merge import merge_all_cohorts

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run TCGA multi-omics merge.")
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--output-dir", type=Path, default=Path("data"))
    args = parser.parse_args()

    if not args.data_dir.exists():
        raise SystemExit(f"data-dir does not exist: {args.data_dir.resolve()}")

    out = merge_all_cohorts(args.data_dir, args.output_dir)
    print(f"Merged matrix written to {out}")
```

---

_Reviewed: 2026-04-27T21:25:53Z_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
