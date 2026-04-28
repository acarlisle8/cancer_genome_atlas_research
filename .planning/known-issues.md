# Known Issues

Living document — append new issues as we hit them. Each entry has a short
description, the workaround, and the rationale so future-us can revisit if a
library update fixes the underlying bug.

---

## Polars `scan_csv` fails on S3 URLs when `has_header=False`

**First seen:** 2026-04-28, during initial multi-cohort ingest run.

**Symptom:** Calling `pl.scan_csv("s3://...", has_header=False, ...)` raises
`FileNotFoundError: No such file or directory (os error 2)` on the S3 URL,
even when the object exists. RNA-seq and CNV ingestors (which use
`has_header=True`) work fine on identical bucket paths.

**Diagnosis:** Polars takes a different internal code path when
`has_header=False`. That path appears to skip the URL-scheme detection that
normally routes `s3://` URLs through the `object_store` Rust crate, and
instead treats the path as a local filesystem path — which obviously doesn't
resolve. Confirmed by isolating: the *only* difference between a working
single-file scan and the broken one was the `has_header` flag.

**Workaround:** [src/ingest_methylation_polars.py](../src/ingest_methylation_polars.py)
sets `has_header=True` and provides a `schema={"probe_id": pl.String,
"beta_value": pl.Float64}` to override the column names polars would derive
from the first row. The first row of each headerless TSV is consumed as the
"header" and lost.

**Cost:** 1 probe lost per patient file. Each methylation file has ~850,000
probes; downstream feature selection picks the top-500 by variance across
all three cohorts pooled. The probability that the lost probe is one of the
top-500 is ~0.06%, and the variance computation is across pooled cohorts so
even one unlucky cohort doesn't drop a feature globally.

**Revisit when:** Polars releases a fix for `scan_csv` headerless + remote
URI handling. Current pin is whatever `uv.lock` resolves on this branch.

---

## GDC `Masked Copy Number Segment` files use bare chromosome names ("1", "X")

**First seen:** 2026-04-28, debugging "Merged matrix is empty" error in run_merge.py.

**Symptom:** [src/merge.py](../src/merge.py) `_pivot_cnv` filtered every row out
of the CNV parquet, then pivoted to an empty wide frame, then the three-way
inner join in `merge_all_cohorts` produced zero rows for every cohort. Final
`pl.concat` raised `RuntimeError: Merged matrix is empty`.

**Diagnosis:** The GDC `Masked Copy Number Segment` TSV files store the
`Chromosome` column as bare values like `"1"`, `"2"`, ..., `"22"`, `"X"` —
*without* a `"chr"` prefix. Both [src/ingest_cnv.py](../src/ingest_cnv.py) and
[src/ingest_cnv_polars.py](../src/ingest_cnv_polars.py) preserve this format
verbatim. Aidan's `_pivot_cnv` filter, however, used
`pl.col("chromosome").is_in(list(HG38_CENTROMERES.keys()))` where the dict
keys are `"chr1"..."chr22","chrX"` — so the filter never matched anything.
Confirmed via DuckDB: `WHERE chromosome IN ('chr1',...) → 0 rows`,
`WHERE chromosome IN ('1',...) → 463,441 rows`.

The unit tests in `tests/test_merge.py` masked the bug because the synthetic
fixtures wrote `"chr1"` directly, which matched the filter.

**Workaround:** [src/merge.py](../src/merge.py) `_pivot_cnv` now normalizes the
chromosome column at the top of the function — adds the `"chr"` prefix
unconditionally if it's missing. Robust to either input format. The arm-name
logic later (`str.replace("chr", "")`) was already format-agnostic, so the
rest of the function is unchanged.

**Cost:** None — the normalization is a one-pass transform on a column that's
already loaded.

**Revisit when:** GDC standardizes their seg-file format. Don't hold breath.

---
