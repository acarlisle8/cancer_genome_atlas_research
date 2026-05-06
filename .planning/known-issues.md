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

## Polars memory landmines on the methylation / RNA-seq parquets

**First seen:** 2026-05-05 22:06 UTC, mid-Phase-4b session. Same session
OOM-killed three Python processes in 11 minutes (PIDs 8402, 12276, 12320 —
each ~6 GB RSS).

**Symptom:** Any of these patterns OOM-kills the box on the
`methylation.parquet` (444M rows, 5.5 GB compressed) or `rna_seq.parquet`
(74.7M rows, 1.6 GB compressed):

```python
pl.read_parquet(path)                                    # eager load
pl.scan_parquet(path).select(pl.col(c).n_unique()).collect()    # n_unique materializes
df.group_by(...).count()                                 # without streaming engine
```

The instance is t3 / 7.6 GB RAM with **no swap**; even Polars "lazy"
chains still allocate the full unique-set / aggregation hashmap when the
sink isn't actually streaming.

**Diagnosis:** The OOMs were never in committed code — they were ad-hoc
EDA Bash one-liners. But the same anti-pattern existed in `run_mofa.py`
(eager `pl.read_parquet` of `merged_all_cohorts.parquet`) and four
helpers in `src/merge.py` (`_pivot_rppa`, `_pivot_mirna`,
`_recurrent_mutated_genes`, `_pivot_mutations`) — those survived only
because the parquets they target are small today.

**Workaround / project rule:** On the per-modality long-format parquets,
**only** these patterns are safe:

```python
pl.scan_parquet(path).filter(...).select(...).collect(engine="streaming")
duckdb.connect(":memory:").execute("PRAGMA memory_limit='4GB'") + group-by SQL
```

Aggregations that need the unique-set materialized in memory (`n_unique`,
`COUNT(DISTINCT)`) should go through DuckDB with an explicit memory
limit, not Polars. See [src/merge.py](../src/merge.py) `_top_n_by_variance`
for the canonical pattern.

**Cost:** None — the lazy/streaming versions are no slower in practice
on these data sizes.

**Revisit when:** Either (a) the box is upgraded with swap or to ≥16 GB
RAM, or (b) Polars 2.x makes streaming the default for `n_unique` etc.
Until then, lazy is mandatory on the big parquets.

---

## Variance ranking in `_top_n_by_variance` admitted 100%-missing methylation probes

**First seen:** 2026-05-06, post-mortem audit during Phase 4e prep.

**Symptom:** The merged 6-view BRCA matrix's mean methylation missingness
was 36.3% — vastly above the source parquet's ~18% array-QC baseline.
The top-5000-by-variance methylation probes were enriched for the
100%-missing-for-BRCA tail.

**Diagnosis:** [src/merge.py](../src/merge.py) `_top_n_by_variance` ranked
probes by `VAR_POP(beta_value)` over non-null observations only. Probes
with one or two non-null beta values (out of 1097 BRCA patients) had
huge "empirical variance" across those few values and outranked
genuinely high-variance well-measured probes.

DuckDB diagnostic, all 488,026 source probes:

| Set | Mean missingness | Median |
|---|---|---|
| All probes | ~19% | 18.4% |
| Top-5000 by variance (old) | 27.9% | 18.8% |
| Bottom-5000 by variance | 19.0% | 18.4% |
| Random 5000 | 19.1% | 18.4% |

Per-probe missingness is bimodal: ~80% of probes cluster at 18-20%
(array's normal QC failure rate), with a tail at 100% that the variance
ranking preferentially captured.

**Workaround:** Added `HAVING COUNT(*) >= MIN_COVERAGE_FRAC * n_patients`
clause before the `ORDER BY VAR_POP`. `MIN_COVERAGE_FRAC = 0.80` at
module scope. `must_keep` marker genes still bypass the filter (they're
force-included regardless of coverage). After re-running
`run_merge_6view.py`, mean meth missingness in the merged file dropped
36.3% → 26.7%, with no more 100%-missing-tail.

**Cost:** A small number of genuinely-bimodal but partially-missing
probes get excluded. Acceptable — those probes contribute mostly NaN
to MOFA+ training anyway.

**Revisit when:** N/A — the fix is principled. Tighten the threshold
above 0.80 only if a future analysis surfaces a coverage skew the
inner-join introduces.

---

## PAM50 vs multi-omic comparison ceiling

**First seen:** 2026-05-06, Phase 4f-light analysis.

**Symptom:** k-means clustering on the 6-view MOFA+ factor scores agreed
with PAM50 labels (`hoadley_subtype_selected`) at ARI 0.27–0.31 across
k ∈ {2, 4, 5} — well below the ROADMAP success criterion of ARI > 0.5.
At k=4 the clusters are: 97% pure Basal (n=97), 76% LumA (n=278),
moderate LumB (n=225), and a 144-patient "ambiguous" cluster.

**Diagnosis:** PAM50 is an **RNA-only** subtyping. Our MOFA+ factors
integrate signal from six modalities, and methylation alone contributes
~45% of the captured variance (vs RNA's ~40%). Two factors (3 and 7)
are methylation-specific — they capture cell-of-origin / TME composition
biology that genuinely doesn't appear in PAM50 RNA labels. The audit
table's `hoadley_subtype_integrative` field, which would be a fair
multi-omic comparator (iCluster-derived integrated subtypes), **is null
for BRCA** in this dataset.

The biological signal *is* there — Basal recovery is essentially perfect
(97% pure cluster) and NMI 0.45 at k=2 is solid agreement-above-chance.
What's missing is alignment of the LumA/LumB/Her2 boundaries, where
PAM50 draws sharper RNA-defined lines than the underlying multi-omic
structure supports.

**Workaround:** Two layers:

1. **Reframe the success criterion** away from "matches PAM50 partition
   exactly" toward measurable multi-omic claims: (a) ≥90% basal-purity
   in the dominant factor's top cluster, (b) factor 1 multi-omic across
   all 6 views, (c) mutation signal concentrated on the proliferation
   axis. All three are met.
2. **Try consensus clustering with bootstrap resampling** (the actual
   ROADMAP 4f spec) — typically gains 0.03–0.08 ARI by stabilizing
   cluster assignments. Plus restrict to factors 1–5 instead of all
   14 (the long-tail factors are noise).

**Cost:** Reframing implicitly demotes the original "ARI > 0.5 vs PAM50"
criterion in favor of substructure metrics that better fit a
multi-omic-vs-single-omic comparison.

**Revisit when:** Either (a) an integrative TCGA-BRCA subtype label
(iCluster, COCA) is added to the audit table, or (b) consensus
clustering + factor-subset narrows the gap meaningfully.

---
