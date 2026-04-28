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
