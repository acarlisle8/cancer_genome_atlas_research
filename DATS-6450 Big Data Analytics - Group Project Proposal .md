DATS 6540 \- Big Data Analytics  
Group 2 \- Project Proposal  
Team Members: Aidan Carlisle and Zachary Cardell

Dataset: 

1) GDC’s The Cancer Genome Atlas (TCGA): [https://registry.opendata.aws/tcga/](https://registry.opendata.aws/tcga/)

The growth of data volume and computational demand has impacted many professional fields, including finance, marketing and biology. One such field, cancer genomics demonstrates the need for improved data engineering techniques, including distributed cloud computing. The Cancer Genome Atlas (TCGA) covers several molecular data types and methodologies, including gene expression, copy number variation, and DNA methylation. This atlas also spans 33 cancer types, while being spread across public repositories in inconsistent per-patient file formats. Distributed cloud computing provides the infrastructure to build a stable, reproducible ETL workflow that can unify these sources into a single integrated dataset, making downstream tasks like molecular cancer subtype classification feasible.   
This project builds a pipeline that ingests TCGA data directly from its public S3 bucket, aggregates per-patient files into Parquet files, and joins patients across RNA-seq, copy number, and methylation on shared identifiers. We plan to use lightweight tools such as DuckDB and Polars, for exploratory analysis and single cancer type work, while scaling up to distributed processing on “EMR”/Spark only when the full pan-cancer data volume requires it. We will then train XGBoost classifiers on two tasks: molecular subtype prediction using validated labels, and cancer type classification across all 33 TCGA cohorts. SHAP analysis can be used to identify which genes drive predictions, with results compared against published classifiers for confirmation.