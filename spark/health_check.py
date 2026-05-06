"""
Cancer Genome Atlas Research — Spark Cluster Health Check

Run after `setup-spark-cluster.sh` to confirm the cluster is up,
has the expected parallelism, and can do basic DataFrame ops.

Optional S3 read probe: set the S3_PROBE_PATH env var to a known
TCGA file path to verify Spark can read from S3 with the configured
IAM-instance-profile credentials.

Usage:
    uv run python spark/health_check.py spark://<MASTER_PRIVATE_IP>:7077

    # With S3 probe:
    S3_PROBE_PATH=s3a://g23861422-datsbd-s2026/tcga/rppa/<uuid>/<file>.tsv \\
        uv run python spark/health_check.py spark://<MASTER_PRIVATE_IP>:7077

The default expected parallelism is 6 (3 workers x 2 cores on t3.large).
Override with EXPECTED_PARALLELISM env var if cluster shape differs.
"""

import os
import sys
import time

from pyspark.sql import SparkSession

MASTER_URL = sys.argv[1] if len(sys.argv) > 1 else "local[*]"
EXPECTED_PARALLELISM = int(os.environ.get("EXPECTED_PARALLELISM", "6"))
S3_PROBE_PATH = os.environ.get("S3_PROBE_PATH")

print(f"Connecting to {MASTER_URL} ...")
t0 = time.time()

spark = (
    SparkSession.builder
    .master(MASTER_URL)
    .appName("CGAR-HealthCheck")
    .config("spark.jars.packages",
            "org.apache.hadoop:hadoop-aws:3.3.4,"
            "com.amazonaws:aws-java-sdk-bundle:1.12.262")
    .config("spark.hadoop.fs.s3a.aws.credentials.provider",
            "com.amazonaws.auth.InstanceProfileCredentialsProvider")
    .config("spark.hadoop.fs.s3a.endpoint", "s3.amazonaws.com")
    .config("spark.hadoop.fs.s3a.path.style.access", "false")
    .getOrCreate()
)
spark.sparkContext.setLogLevel("WARN")

# 1. Parallelism
parallelism = spark.sparkContext.defaultParallelism
print(f"[1] defaultParallelism = {parallelism}")
assert parallelism >= EXPECTED_PARALLELISM, (
    f"Expected >= {EXPECTED_PARALLELISM}, got {parallelism}. "
    "Workers may not all be connected. Check Spark UI: http://<MASTER_PUBLIC_IP>:8080"
)
print(f"    OK")

# 2. Basic DataFrame
print("[2] Basic DataFrame test ...")
df = spark.range(1_000_000)
even_count = df.filter(df.id % 2 == 0).count()
assert even_count == 500_000, f"Expected 500000 even rows, got {even_count}"
print(f"    OK — {even_count} rows after filter")

# 3. Optional S3 probe
if S3_PROBE_PATH:
    print(f"[3] S3 read probe: {S3_PROBE_PATH} ...")
    if S3_PROBE_PATH.endswith(".parquet") or "/parquet" in S3_PROBE_PATH:
        sample = spark.read.parquet(S3_PROBE_PATH).limit(5).collect()
    else:
        sample = spark.read.text(S3_PROBE_PATH).limit(5).collect()
    assert len(sample) > 0, f"Read 0 rows from {S3_PROBE_PATH}"
    print(f"    OK — read {len(sample)} rows from S3")
else:
    print("[3] S3 probe skipped (set S3_PROBE_PATH env var to enable)")

elapsed = time.time() - t0
print(f"\n[4] Total elapsed: {elapsed:.1f}s")
print("\nHealth check PASSED.")

spark.stop()
