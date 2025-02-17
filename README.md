# ProteoDeconv Manuscript Code

This repository contains the code and analysis pipeline for the proteoDeconv manuscript (currently under review).

## Overview

The pipeline leverages the R packages {targets, renv} to ensure reproducibility and effective environment management.

## Getting Started

### Prerequisites

Ensure the following system dependencies are installed:
- libcurl4-openssl-dev
- libssl-dev
- libpng-dev
- libxml2-dev
- libmagick++-dev
- cmake
- libmbedtls-dev
- libharfbuzz-dev
- libfribidi-dev

Additionally, install Quarto and Docker if required.

### Data

Required datasets:
- **PXD004352** – proteinGroups.txt  
  Place in: `data/raw/rieckmann/total`
- **PXD056050** (currently password-protected) – report.pg_matrix.tsv and report.unique_genes_matrix.tsv  
  Place in: `data/raw/immune_cells`
- **scRNA-seq Data** – Download "NSCLC PBMCs Single Cell RNA-Seq (Fig. 2ab)" from the CIBERSORTx website and place `Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt` in `data/raw/NSCLC_PBMCs/`

Note: The pipeline is designed to tolerate missing datasets; individual targets will fail when data is unavailable.

## Setup

Restore the R environment by running:

```{R}
renv::restore()
```

## Running the Pipeline

Execute the pipeline with:
```{R}
targets::tar_make()
```

## License

The code is licensed under the MIT License.