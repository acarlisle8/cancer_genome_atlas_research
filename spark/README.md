# Spark Cluster Scaffolding

Automated setup for a 4-node Apache Spark cluster on AWS EC2 (1 master + 3 workers, t3.large, Ubuntu 22.04). Pulled and adapted from `midterm-02-zcardell-gwu/`. Used by Phase 4 for BRCA multi-omic ingest, MOFA+, and SNF.

## Layout

```
spark/
├── setup-spark-cluster.sh        Bring up the cluster (run from this dir)
├── cleanup-spark-cluster.sh      Tear it down
├── health_check.py               Post-setup sanity check from dev machine
├── cluster-files/                Files SCPed to every cluster node
│   ├── pyproject.toml            Cluster-side Python deps (uv sync on each node)
│   └── spark_installation_test.py   Spark + AWS jars sanity test (run on master)
└── README.md                     This file
```

The following are **generated** by `setup-spark-cluster.sh` and gitignored — do not commit:

- `cluster-config.txt`, `cluster-ips.txt` — cluster IPs and instance IDs
- `spark-cluster-key-*.pem` — SSH key (secret)
- `cluster-setup-*.log` — setup logs
- `ssh_to_master_node.sh` — auto-generated SSH helper

## Usage

### Bring up the cluster

```bash
# Get your laptop's public IP from https://ipchicken.com/
cd spark/
./setup-spark-cluster.sh <YOUR_LAPTOP_IP>
```

The script auto-detects the IAM instance profile of the dev EC2 machine and propagates it to the cluster nodes for S3 access.

Setup output goes to `spark/cluster-setup-YYYYMMDD-HHMMSS.log`. On success, you'll see the master URL `spark://<MASTER_PRIVATE_IP>:7077`.

### Verify the cluster

From the dev machine:

```bash
uv run python spark/health_check.py spark://<MASTER_PRIVATE_IP>:7077
```

Optional S3 read probe:

```bash
S3_PROBE_PATH=s3a://g23861422-datsbd-s2026/tcga/rppa/<uuid>/<file>.tsv \
    uv run python spark/health_check.py spark://<MASTER_PRIVATE_IP>:7077
```

### Run installation test on the master node

```bash
./ssh_to_master_node.sh
cd ~/spark-cluster
uv run python spark_installation_test.py
```

### Tear down

```bash
cd spark/
./cleanup-spark-cluster.sh
```

This terminates all 4 EC2 instances, deletes the security group + key pair, and removes the local `cluster-config.txt`, `cluster-ips.txt`, and `*.pem` files.

## Cluster shape

| Component | Count | Type | Disk |
|---|---|---|---|
| Master | 1 | t3.large (2 vCPU, 8 GB) | 100 GB gp3 |
| Worker | 3 | t3.large (2 vCPU, 8 GB) | 100 GB gp3 |

Total: 8 vCPU, 32 GB RAM across the cluster.

Spark 3.4.4 on Hadoop 3.3.4 with `hadoop-aws:3.3.4` + `aws-java-sdk-bundle:1.12.262` jars for S3a access.

## Notes

- The `cluster-files/pyproject.toml` deps install on **every node** via `uv sync` during setup. Includes `mofapy2` so the master can run MOFA+ orchestration jobs distributed via Spark scheduler.
- Region is hardcoded to `us-east-1` (where TCGA data is hosted).
- Always run `setup-spark-cluster.sh` and `cleanup-spark-cluster.sh` from the `spark/` directory — they look for `cluster-files/` and `cluster-config.txt` in the current working directory.
