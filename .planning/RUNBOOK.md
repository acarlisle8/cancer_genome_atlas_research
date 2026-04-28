# Runbook — How to run this end-to-end

Operational guide for executing Phase 1 (ingest) → Phase 2 (merge) on a fresh
EC2 instance. Written 2026-04-28 after the first successful end-to-end run.

For project context, read [PROJECT.md](PROJECT.md). For phase status, read
[ROADMAP.md](ROADMAP.md). For non-obvious bugs and workarounds, read
[known-issues.md](known-issues.md).

---

## TL;DR

```bash
# One-time setup
git clone <repo> && cd cancer_genome_atlas_research
uv sync
source .venv/bin/activate

# AWS credentials must be configured — see "S3 access" below
aws configure   # or attach EC2 instance profile

# Run pipeline
python run_ingest.py    # ~30-60 min, writes ~17 GB to data/
python run_merge.py     # ~5 min, writes data/merged_all_cohorts.parquet (~11 MB)

# Verify (paste from .planning/phases/02-merge/02-VERIFY-CHECKLIST.md)
```

---

## What you'll get

After both scripts run, `data/` contains:

```
data/
├── TCGA-BRCA/
│   ├── rna_seq.parquet         (~1.5 GB,  74.6M rows)
│   ├── cnv.parquet             (~7 MB,    463K rows)
│   └── methylation.parquet     (~5.5 GB,  444M rows)
├── TCGA-LUAD/                  (~4.6 GB total)
├── TCGA-PRAD/                  (~4.1 GB total)
└── merged_all_cohorts.parquet  (11 MB, 2104 rows × 1043 cols)
```

The **`merged_all_cohorts.parquet`** is the Phase 2 deliverable that Phase 3
(XGBoost + SHAP) will consume. It's committed to the repo (carved out as an
exception in `.gitignore`) so Phase 3 can run against it without re-syncing
the 17 GB of cohort parquets.

---

## Prerequisites

### Compute
- **EC2 instance**: `t3.large` or larger (8 GB RAM, 2 vCPU). The 8 GB minimum
  is real — we hit OOM repeatedly on smaller and had to switch portions of the
  pipeline from polars to DuckDB streaming aggregates.
- **Disk**: 25 GB free for `data/` + venv + cache.
- **Python 3.12** (3.10+ should work; we developed on 3.12.3).
- **`uv`** for dependency management (`curl -LsSf https://astral.sh/uv/install.sh | sh`).

### S3 access — important caveat

The polars ingestors in `src/ingest_*_polars.py` read from
`s3://g23861422-datsbd-s2026/tcga/` — that's **Zachary's personal AWS S3
bucket**, populated via `sync_to_s3.py` from the public TCGA bucket. You have
three options:

1. **Get IAM read access to Zachary's bucket** (talk to Zachary). Cheapest,
   fastest.
2. **Re-sync to your own bucket**:
   ```bash
   # Edit sync_to_s3.py: change DEST_BUCKET to your bucket
   # Edit src/gdc_client.py: change TCGA_S3_BUCKET to match
   python sync_to_s3.py --project TCGA-BRCA TCGA-LUAD TCGA-PRAD
   ```
   This pulls from `s3://tcga-2-open/` (public) and writes to your bucket.
   Storage cost is on you (~25 GB for 3 cohorts).
3. **Skip Phase 1 entirely**: `data/merged_all_cohorts.parquet` (11 MB) is
   already committed to the repo. If all you need is to run Phase 3, just
   `git pull` and it's there.

The original DuckDB ingestors (`src/ingest_*.py` without `_polars`) read
directly from `s3://tcga-2-open/` (public) but were superseded by the polars
versions for throughput reasons. They still work if you'd rather skip the
personal-bucket detour — `src/gdc_client.py:TCGA_S3_BUCKET` controls the
source.

### AWS auth resolution
Both DuckDB and polars resolve credentials via the AWS SDK default chain:
env vars → `~/.aws/credentials` profile → SSO → EC2 instance profile.
Configured via `aws configure` or by attaching an IAM role to the EC2
instance.

---

## Run sequence

### Step 1: Ingest (Phase 1)
```bash
python run_ingest.py
```

What it does: orchestrates 3 cohorts × 3 modalities = 9 ingestion jobs. Each
job uses polars `scan_csv` to lazily stream patient TSVs from S3 and
`sink_parquet` to write a consolidated long-format parquet per
(cohort, modality). Resumable — skips outputs that already exist.

Expected duration: **30–60 min** depending on instance and S3 throughput.
The methylation files are the long pole (~12 GB to scan).

Outputs: 9 files at `data/{cohort}/{modality}.parquet`.

### Step 2: Merge (Phase 2)
```bash
python run_merge.py
```

What it does:
1. DuckDB streaming aggregate over the 3 RNA-seq parquets to pick the top-500
   gene_ids by variance (~25s).
2. Same for the 3 methylation parquets to pick the top-500 probe_ids by
   variance (~3 min — heaviest single step).
3. For each cohort:
   - `_pivot_rnaseq` filters to top genes, pivots wide (~10s)
   - `_pivot_cnv` aggregates segments to chromosome-arm means, pivots (<1s)
   - `_pivot_methylation` filters to top probes, pivots (~30s)
   - 3-way inner join on `patient_id`
4. `pl.concat(diagonal_relaxed)` stacks the three cohort frames; writes the
   merged matrix.

Expected duration: **~5 min total**.

Output: `data/merged_all_cohorts.parquet` — 2104 rows × 1043 cols.

### Step 3: Verify
Paste the script from
[.planning/phases/02-merge/02-VERIFY-CHECKLIST.md](phases/02-merge/02-VERIFY-CHECKLIST.md)
into a Python REPL or run it as a script. All 11 checks should print `OK`.

---

## What was changed since `9d75d37` (Aidan's original Phase 2 commit)

These commits land between `9d75d37` (Aidan's Phase 2 implementation) and
HEAD on the `claude/polars-full-pipeline` branch. Listed in the order they
were made:

| Commit | Reason |
|---|---|
| `6ccade1` | Lab-05 streaming pattern in `_top_n_by_variance` (column projection + `engine="streaming"`); RNA-seq switched to `fpkm_uq_unstranded` for TCGA-standard normalization. |
| `9f699f9` | Switched `_top_n_by_variance` from polars streaming to DuckDB — polars streaming engine fell back to eager on `var() + sort + head` and OOM'd 8 GB instance. |
| `88287d1` | Rewrote `_pivot_rnaseq` with lazy scan + filter pushdown + column projection — eager `pl.read_parquet` on the 1.5 GB BRCA RNA-seq file OOM'd. |
| `695e96a` | **The merge bug fix.** `_pivot_cnv` now normalizes chromosome to "chr"-prefix — GDC seg files emit bare "1"/"X", so the canonical filter dropped every row and produced the empty merge. Also: empty-frame guard, lazy scan, natural arm sort, `_CANONICAL_CHROMS` is now actually used. |
| `8898233` | Logging in `merge_all_cohorts`: per-cohort/per-pivot shape lines so the next silent zero-row failure is debuggable in seconds. |
| `f2b295b` | Test fixtures aligned to real ingest format: `fpkm_uq_unstranded` column name + bare chromosome `"1"`/`"M"` (not `"chr1"`/`"chrM"`). The synthetic-data mismatch is what hid the chromosome bug from the test suite. |

Plus several earlier ingest/sync commits (`5a3cb3c`, `785eaa2`, `dcb18c6`,
`a382582`, etc.) that landed Phase 1.

---

## Troubleshooting

### `MemoryError` or instance hangs

We OOM'd this 8 GB instance multiple times during development. Patterns:

- **Polars streaming engine silently falls back to eager** on certain
  operation chains (`var() + sort + head`, `unique()` on high-cardinality
  columns). If you see RAM climbing instead of staying flat, suspect this.
  The fix is usually to push that step through DuckDB —
  `_top_n_by_variance` is the canonical example.
- **Eager `pl.read_parquet` on multi-GB files**: replace with
  `pl.scan_parquet(...).filter(...).select([...]).collect()` to push
  predicates and projections into the parquet reader.
- **DuckDB without memory cap**: always set `PRAGMA memory_limit='4GB'`
  (or similar) when running streaming aggregates so it spills to disk
  rather than OOMs. See [src/merge.py](../src/merge.py) `_top_n_by_variance`.

### Merge succeeds but `merged_all_cohorts.parquet` is empty
Should be impossible after `695e96a` — but if you see it, run the per-pivot
shape lines now logged inside `merge_all_cohorts` (Step 2 above) to identify
which pivot is producing 0 rows.

### S3 access denied
`gdc_client.py:get_duckdb_conn` resolves AWS creds via the default chain.
Verify with `aws s3 ls s3://g23861422-datsbd-s2026/tcga/` — if that errors,
your creds aren't reaching DuckDB/polars either.

### Methylation ingest produces "1 probe lost per file" warnings
Expected. See [known-issues.md](known-issues.md) — polars `scan_csv` has a
bug with headerless CSVs over S3 and we work around it by treating the first
row as a header. The lost probe (one per ~850K) is statistically irrelevant
to top-500 feature selection.

---

## After Phase 2 → before Phase 3

The merged matrix is the input to XGBoost. Expected shape:
`2104 rows × 1043 cols`, with columns:

- `patient_id` (str)
- 500 RNA-seq gene columns (`ENSG...`)
- 41 CNV chromosome-arm columns (`1p`, `1q`, ..., `Xp`, `Xq`)
- 495 methylation probe columns (`cg...`) — note: 5 less than the 500 selected,
  see [02-HUMAN-UAT.md](phases/02-merge/02-HUMAN-UAT.md) "Gaps"
- `cohort` (str: BRCA / LUAD / PRAD)

Phase 3 entry point will be a new `src/classify.py` with an XGBoost
multiclass classifier on the `cohort` label, plus SHAP summary plots.
