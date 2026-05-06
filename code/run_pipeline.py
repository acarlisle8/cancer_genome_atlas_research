from pathlib import Path
from src.ingest_rnaseq import ingest_rnaseq
from src.ingest_cnv import ingest_cnv
from src.ingest_methylation import ingest_methylation

OUT = Path("data/processed")

ingest_rnaseq(OUT, project_id="TCGA-BRCA")
ingest_cnv(OUT, project_id="TCGA-BRCA")
ingest_methylation(OUT, project_id="TCGA-BRCA")
