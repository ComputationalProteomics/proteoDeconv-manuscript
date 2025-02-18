# Proteomics Cell-Type Deconvolution Pipeline
This codebase contains the complete analysis pipeline used in our manuscript (currently under review) for proteomics deconvolution.

## Overview

The pipeline is built using these key R packages:
- {targets} for pipeline management and reproducibility
- {renv} for consistent package versioning and environment control
- {proteoDeconv}, our package for facilitating proteomics cell-type deconvolution


### System Requirements

First, install these system libraries:

```{bash}
libcurl4-openssl-dev
libssl-dev
libpng-dev
libxml2-dev
libmagick++-dev
cmake
libmbedtls-dev
libharfbuzz-dev
libfribidi-dev
```

You'll also need:
- Quarto for report generation
- Docker (for running the CIBERSORTx container)

### Required Datasets

Our pipeline works with three key datasets:

1. **Proteomics Data (PXD004352)**
   - File: `proteinGroups.txt`
   - Place in: `data/raw/rieckmann/total`

2. **Immune Cell Data (PXD056050)**
   - Currently password-protected
   - Files: 
     - `report.pg_matrix.tsv`
     - `report.unique_genes_matrix.tsv`
   - Place in: `data/raw/immune_cells`

3. **scRNA-seq Reference Data**
   - Download "NSCLC PBMCs Single Cell RNA-Seq (Fig. 2ab)" from CIBERSORTx
   - Place `Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt` in `data/raw/NSCLC_PBMCs/`

The pipeline can run with incomplete data - for example, if the scRNA-seq reference data is missing, only the steps requiring that dataset will be skipped.

### CIBERSORTx Setup

1. Download the `CIBERSORT.R` script from CIBERSORTx and place it in `R/cibersort`
2. Request a token from CIBERSORTx
3. Create an `.Renviron` file with your credentials:
   ```R
   CIBERSORTX_TOKEN = your_token_here
   CIBERSORTX_EMAIL = your_email_here
   ```

## Running the Analysis

1. Set up your R environment:
   ```R
   renv::restore()
   ```

2. Launch the pipeline:
   ```R
   targets::tar_make()
   ```

