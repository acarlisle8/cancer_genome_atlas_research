"""Canonical subtype-marker gene panels.

These are domain-defined gene panels used for molecular subtype classification.
Force-including them in the merged feature matrix prevents variance-based
filtering from accidentally dropping the genes that actually define subtype
in published classifiers (PAM50 for BRCA, etc.).

Ensembl IDs are stored *without* version suffix because GENCODE annotation
versions drift over time. Match against the parquet's versioned `gene_id`
column via prefix (e.g., "ENSG00000091831." matches "ENSG00000091831.16").

References:
    PAM50 BRCA: Parker et al. 2009; TCGA Nature 2012
    LUAD     : Wilkerson 2012; TCGA Nature 2014
    PRAD     : TCGA Cell 2015 (iCluster + AR/ERG/FOXA1 expression markers)
"""
from __future__ import annotations


PAM50_BRCA: dict[str, str] = {
    "ESR1":  "ENSG00000091831",
    "PGR":   "ENSG00000082175",
    "ERBB2": "ENSG00000141736",
    "MKI67": "ENSG00000148773",
    "FOXA1": "ENSG00000129514",
    "FOXC1": "ENSG00000054598",
    "BCL2":  "ENSG00000171791",
    "CDH1":  "ENSG00000039068",
    "MYC":   "ENSG00000136997",
    "BIRC5": "ENSG00000089685",
    "AURKA": "ENSG00000087586",
    "CCNB1": "ENSG00000134057",
    "MYBL2": "ENSG00000101057",
    "KRT5":  "ENSG00000186081",
    "KRT14": "ENSG00000186847",
    "KRT17": "ENSG00000128422",
    "EGFR":  "ENSG00000146648",
    "FGFR4": "ENSG00000160867",
    "CDC20": "ENSG00000117399",
    "MELK":  "ENSG00000165304",
}

LUAD_MARKERS: dict[str, str] = {
    "NKX2-1": "ENSG00000136352",
    "SFTPB":  "ENSG00000168878",
    "SFTPC":  "ENSG00000168484",
    "TP63":   "ENSG00000073282",
    "KRT5":   "ENSG00000186081",
    "KRT6A":  "ENSG00000205420",
    "MUC1":   "ENSG00000185499",
}

PRAD_MARKERS: dict[str, str] = {
    "AR":    "ENSG00000169083",
    "ERG":   "ENSG00000157554",
    "ETV1":  "ENSG00000006468",
    "ETV4":  "ENSG00000175832",
    "FLI1":  "ENSG00000151702",
    "SPOP":  "ENSG00000121067",
    "FOXA1": "ENSG00000129514",
    "PTEN":  "ENSG00000171862",
    "TP53":  "ENSG00000141510",
}

PANELS: dict[str, dict[str, str]] = {
    "BRCA": PAM50_BRCA,
    "LUAD": LUAD_MARKERS,
    "PRAD": PRAD_MARKERS,
}

# Union of all marker Ensembl IDs (no version suffix), used by the merge step
# to force-include marker genes alongside variance-picked genes.
ALL_RNA_MARKERS: set[str] = (
    set(PAM50_BRCA.values())
    | set(LUAD_MARKERS.values())
    | set(PRAD_MARKERS.values())
)
