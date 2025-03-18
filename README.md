# Code accompanying the manuscript "Considerations and Software for Successful Immune Cell Deconvolution using Proteomics Data"

This codebase contains the complete analysis pipeline used in our manuscript (currently under review) for proteomics deconvolution.

## Overview

The pipeline is built using these key R packages:

- `targets` for pipeline management and reproducibility
- `renv` for consistent package versioning and environment control
- `proteoDeconv`, our package for facilitating proteomics cell-type deconvolution

For more information about the proteoDeconv package, please visit the [proteoDeconv repository](https://github.com/ComputationalProteomics/proteoDeconv).

## Getting Started

First, clone this repository to your local machine:

```bash
git clone https://github.com/ComputationalProteomics/proteoDeconv-manuscript.git
cd proteoDeconv-manuscript
```

After cloning, you'll need to add the required data files, CIBERSORT.R script, and .Renviron file as described in the sections below.

## Installation & Setup

You can run this pipeline using one of two approaches:

### Option 1: Using Docker

1. **Install Docker** on your system.

2. **Pull and run the Docker image** from within the cloned repository
   directory:

   ```bash
   docker run -d \
   -p 127.0.0.1:8787:8787 \
   -v "$(pwd):/home/rstudio/proteoDeconv-manuscript" \
   -v /home/rstudio/proteoDeconv-manuscript/renv \
   -v /var/run/docker.sock:/var/run/docker.sock \
   manszamore/proteodeconv-manuscript:latest
   ```

   Note: The Docker socket is mounted (via -v /var/run/docker.sock:/var/run/docker.sock) to enable Docker-in-Docker functionality, which is required to run the CIBERSORTx Docker image from within this container.


3. **Access RStudio Server**:
   - Open your browser and navigate to `http://localhost:8787/`
   - The Docker image comes with all required dependencies pre-installed

### Option 2: Local Installation

1. **System Requirements**:
   - Tested on Ubuntu and macOS
   - R version 4.4.2 (we recommend using [rig](https://github.com/r-lib/rig) for R version management)
   - Quarto for report generation
   - Docker (for running the CIBERSORTx container)

2. **System Libraries** (Ubuntu/Debian):

   ```bash
   apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev make libgsl0-dev libglpk-dev \
    libicu-dev libpng-dev libxml2-dev python3 libx11-dev cmake xz-utils \
    pandoc zlib1g-dev libfontconfig1-dev libfreetype6-dev libfribidi-dev \
    libharfbuzz-dev libjpeg-dev libtiff-dev curl
   ```

3. **R Environment Setup**:

   ```R
   renv::restore()
   ```

## Required Datasets

Regardless of installation method, you'll need to add these datasets:

1. **Rieckmann et al. DDA Data (PXD004352)**
   - File: `proteinGroups.txt`
   - Place in: `data/raw/rieckmann/total`

2. **Immune Cell DIA Data (PXD056050)**
   - Currently password-protected, requiring a reviewer account
   - Files:
     - `report.pg_matrix.tsv`
     - `report.unique_genes_matrix.tsv`
   - Place in: `data/raw/immune_cells`

3. **scRNA-seq Reference Data**
   - Download "NSCLC PBMCs Single Cell RNA-Seq (Fig. 2ab)" from CIBERSORTx
   - Place `Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt` in `data/raw/NSCLC_PBMCs/`

The pipeline can run with incomplete data - for example, if the scRNA-seq reference data is missing, only the steps requiring that dataset will be skipped.

## CIBERSORT Setup

Due to licensing restrictions, we cannot include the CIBERSORT code in this repository. To run the CIBERSORT analysis, you will need to:

1. Download the `CIBERSORT.R` script (version 1.04) from the CIBERSORT website and place it in `R`

### CIBERSORTx Credentials

In order to run CIBERSORTx, you will need to:

1. Request a token from CIBERSORTx (free for academic use)
2. Create an `.Renviron` file in this folder with your credentials:

   ```R
   CIBERSORTX_TOKEN = your_token_here
   CIBERSORTX_EMAIL = your_email_here
   ```

### CIBERSORT setup

1. Download the CIBERSORT.R script from the CIBERSORT website
2. Place it in the `R` directory


## Running the Analysis

Launch the pipeline:

```R
targets::tar_make()
```

The complete pipeline takes approximately 5 hours to run in its entirety.

You can also run specific parts of the pipeline using pattern matching:

```R
# Only run targets related to preprocessing (e.g. merge vs slice, etc.)
targets::tar_make(matches("preprocessing"))
```

## Expected Output

After successful completion, the pipeline will generate:

- Manuscript figures in the `manuscript_figures/` directory
- A supplementary table in `supplementary_table.xlsx` containing simulation results

Furthermore, individual targets can be read using `tar_read` to access the results of each step.

## Issues

If you encounter any issues, please open an issue in the GitHub repository.

## Citation

Please cite our work if you use this pipeline or proteoDeconv in your research. The citation details will be added once the manuscript is published.
