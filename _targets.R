library(targets)
library(tarchetypes)
library(here)
library(crew)
library(autometric)
library(ggplot2)

here::i_am("_targets.R")
tar_config_set(
  store = here::here("_targets"),
  script = here::here("_targets.R")
)


here <- function(...) fs::path_rel(here::here(...))

tar_option_set(
  packages = c("proteoDeconv", "tibble", "tidyverse", "ggplot2", "readxl", "here", "patchwork", "limma", "vsn"),
  imports = "proteoDeconv",
  format = "qs",
  error = "null",
  workspace_on_error = TRUE,
  storage = "worker",
  retrieval = "worker",
  memory = "transient",
  garbage_collection = 10
)

tar_source()


tar_plan(
  # File targets
  tar_file(dda_pg, here("data", "raw", "rieckmann", "total", "proteinGroups.txt")),

  tar_file(dia_pg, here("data", "raw", "immune_cells", "report.pg_matrix.tsv")),
  tar_file(dia_ug, here("data", "raw", "immune_cells", "report.unique_genes_matrix.tsv")),
  tar_file(sc_rna_seq_pbmc, here("data", "raw", "NSCLC_PBMCs", "Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt")),
  file = list(
    dia_pg = dia_pg,
    dia_ug = dia_ug,
    dda_pg = dda_pg,
    sc_rna_seq_pbmc = sc_rna_seq_pbmc
  ),

  # Preprocessed data
  dia_data = preprocess_protein_data(
    input = dia_pg,
    method = "diann",
    imputation_mode = "lowest_value",
    duplicate_mode = "slice",
    remove_samples_pattern = "PBMC|Mix|NK"
  ),
  dia_data_full = preprocess_protein_data(
    input = dia_pg,
    method = "diann",
    imputation_mode = "lowest_value",
    duplicate_mode = "slice",
    remove_samples_pattern = "NK"
  ),
  dia_data_reduced = preprocess_protein_data(
    input = dia_pg,
    method = "diann",
    imputation_mode = "lowest_value",
    duplicate_mode = "slice",
    remove_samples_pattern = "Bnaive|Mix|NKactive|PBMC|Tcell"
  ),
  dda_data = preprocess_protein_data(
    input = dda_pg,
    method = "maxquant",
    imputation_mode = "lowest_value",
    duplicate_mode = "slice",
    remove_samples_pattern = "Erythrocyte|Thrombocyte|Library",
    update_symbols = TRUE
  ),


  sig_dda = create_signature_matrix(
    refsample = dda_data,
    phenoclasses = create_general_phenoclasses(dda_data)
  ),

  signature_params_analysis,
  normalization_analysis,
  algorithm_analysis,
  preprocessing_analysis,
  sim_validation_analysis,
  sig_source_analysis

)
