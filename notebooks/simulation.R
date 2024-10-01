get_simulation_metrics_mixes <- function(source_data = "pg", signature) {
  if (source_data == "pg") {
    data <-
      read_tsv(here("data/raw_internal/bigsort/report.pg_matrix.tsv")) |>
      prepare_diann_matrix() |>
      update_gene_symbols() |>
      handle_gene_groups() |>
      handle_missing_values(imputation_mode = "lowest_value") |>
      handle_duplicate_genes(duplicate_mode = "slice") |>
      select(-matches("PBMC|Mix"))
  } else if (source_data == "ug") {
    data <-
      read_tsv(here("data/raw_internal/bigsort/report.unique_genes_matrix.tsv")) |>
      prepare_diann_matrix() |>
      update_gene_symbols() |>
      handle_gene_groups() |>
      handle_missing_values(imputation_mode = "lowest_value") |>
      handle_duplicate_genes(duplicate_mode = "slice") |>
      select(-matches("PBMC|Mix"))
  }

  annotation <- data.frame(
    "ID" = colnames(data |> select(-Genes)),
    "cell_type" = stringr::str_extract(colnames(data |> select(-Genes)), "(?<=_)[^_]+(?=-)")
  )

  counts <- data |>
    column_to_rownames("Genes") |>
    as.matrix()
  tpm <- data |>
    handle_scaling(unlog = FALSE, tpm = TRUE) |>
    column_to_rownames("Genes") |>
    as.matrix()

  counts_sparse <- Matrix::Matrix(counts, sparse = TRUE)
  tpm_sparse <- Matrix::Matrix(tpm, sparse = TRUE)

  dataset <- SimBu::dataset(
    annotation = annotation,
    count_matrix = counts_sparse,
    tpm_matrix = tpm_sparse,
    name = "A"
  )


  simulation <- SimBu::simulate_bulk(dataset,
    scenario = "random",
    scaling_factor = "NONE",
    seed = 4
  )

  simulated_data <- SummarizedExperiment::assays(simulation$bulk)[["bulk_tpm"]]

  signature <- read_tsv(
    here("data/signature_matrices/LM7c.txt")
  )
  simulated_data <- simulated_data |>
    as.matrix() |>
    as_tibble(rownames = "Genes")
  simulated_deconvoluted <- deconvolute("cibersortx", simulated_data, signature, perm = 1) |>
    pivot_longer(cols = -cell_type, names_to = "sample", values_to = "cell_count")


  cell_fractions <- simulation$cell_fractions |>
    as_tibble(rownames = "sample") |>
    pivot_longer(cols = -sample, names_to = "cell_type", values_to = "cell_count")

  cell_type_mapping <- read_csv(here("data/signature_matrices/source/cell_type_mapping.csv"))
  cell_type_mapping <- setNames(cell_type_mapping$group, cell_type_mapping$cell_type)

  cell_fractions$cell_type <- dplyr::recode(cell_fractions$cell_type, !!!cell_type_mapping)

  joined_df <- dplyr::inner_join(cell_fractions, simulated_deconvoluted, by = c("cell_type", "sample"), suffix = c("_cf", "_sd"))

  rmse <- sqrt(mean((joined_df$cell_count_cf - joined_df$cell_count_sd)^2))

  correlation <- cor(joined_df$cell_count_cf, joined_df$cell_count_sd, method = "pearson")
  correlation2 <- cor(joined_df$cell_count_cf, joined_df$cell_count_sd, method = "spearman")

  metrics <- list(
    rmse = rmse,
    correlation_pearson = correlation,
    correlation_spearman = correlation2
  )
  metrics
}


get_simulation_metrics <- function(source_data = "rieckmann", signature) {
  immune_cells_lowest_value_imputed <-
    read_excel(here("data/raw_external/Rieckmann_2017_S3.xlsx"), na = "NA") |>
    rename(Genes = `Gene names`) |>
    select(Genes, starts_with("LFQ.intensity_")) |>
    update_gene_symbols(verbose = TRUE) |>
    handle_gene_groups() |>
    handle_missing_values(imputation_mode = "lowest_value") |>
    handle_duplicate_genes(duplicate_mode = "slice") |>
    select(-contains(c("Erythrocyte", "Thrombocyte", "Th", "Tregs")))

  annotation <- data.frame(
    "ID" = colnames(immune_cells_lowest_value_imputed |> select(-Genes)),
    "cell_type" = gsub("LFQ.intensity(.imputed)?_(.*)_\\d+_(activated|steady-state)", "\\2_\\3", colnames(immune_cells_lowest_value_imputed |> select(-Genes)))
  )

  cell_type_mapping <- read_csv(here("data/signature_matrices/source/cell_type_mapping.csv"))
  cell_type_mapping <- setNames(cell_type_mapping$group, cell_type_mapping$cell_type)


  annotation$cell_type <- dplyr::recode(annotation$cell_type, !!!cell_type_mapping)


  counts <- immune_cells_lowest_value_imputed |>
    column_to_rownames("Genes") |>
    as.matrix()
  tpm <- immune_cells_lowest_value_imputed |>
    handle_scaling(unlog = FALSE, tpm = TRUE) |>
    column_to_rownames("Genes") |>
    as.matrix()

  counts_sparse <- Matrix::Matrix(counts, sparse = TRUE)
  tpm_sparse <- Matrix::Matrix(tpm, sparse = TRUE)

  dataset <- SimBu::dataset(
    annotation = annotation,
    count_matrix = counts_sparse,
    tpm_matrix = tpm_sparse,
    name = "immune_cells_lowest_value_imputed"
  )


  simulation <- SimBu::simulate_bulk(dataset,
    scenario = "random",
    scaling_factor = "NONE",
    seed = 4
  )

  simulated_data <- SummarizedExperiment::assays(simulation$bulk)[["bulk_tpm"]]


  simulated_data <- simulated_data |>
    as.matrix() |>
    as_tibble(rownames = "Genes")
  simulated_deconvoluted <- deconvolute("cibersortx", simulated_data, signature, perm = 1) |>
    pivot_longer(cols = -cell_type, names_to = "sample", values_to = "cell_count")


  cell_fractions <- simulation$cell_fractions |>
    as_tibble(rownames = "sample") |>
    pivot_longer(cols = -sample, names_to = "cell_type", values_to = "cell_count")

  cell_type_mapping <- read_csv(here("data/signature_matrices/source/cell_type_mapping.csv"))
  cell_type_mapping <- setNames(cell_type_mapping$group, cell_type_mapping$cell_type)


  cell_fractions$cell_type <- dplyr::recode(cell_fractions$cell_type, !!!cell_type_mapping)

  joined_df <- dplyr::inner_join(cell_fractions, simulated_deconvoluted, by = c("cell_type", "sample"), suffix = c("_cf", "_sd"))

  rmse <- sqrt(mean((joined_df$cell_count_cf - joined_df$cell_count_sd)^2))

  correlation <- cor(joined_df$cell_count_cf, joined_df$cell_count_sd, method = "pearson")
  correlation2 <- cor(joined_df$cell_count_cf, joined_df$cell_count_sd, method = "spearman")

  metrics <- list(
    rmse = rmse,
    correlation_pearson = correlation,
    correlation_spearman = correlation2
  )
  metrics
}






get_simulation_metrics2 <- function(data, cell_types, signature, seed) {
  annotation <- data.frame(
    "ID" = colnames(data |> select(-Genes)),
    "cell_type" = cell_types
  )

  counts <- data |>
    column_to_rownames("Genes") |>
    as.matrix()
  tpm <- data |>
    handle_scaling(unlog = FALSE, tpm = TRUE) |>
    column_to_rownames("Genes") |>
    as.matrix()

  counts_sparse <- Matrix::Matrix(counts, sparse = TRUE)
  tpm_sparse <- Matrix::Matrix(tpm, sparse = TRUE)

  dataset <- SimBu::dataset(
    annotation = annotation,
    count_matrix = counts_sparse,
    tpm_matrix = tpm_sparse,
    name = "dataset"
  )

  simulation <- SimBu::simulate_bulk(dataset,
    scenario = "random",
    scaling_factor = "NONE",
    seed = seed
  )

  simulated_data <- SummarizedExperiment::assays(simulation$bulk)[["bulk_tpm"]]

  simulated_data <- simulated_data |>
    as.matrix() |>
    as_tibble(rownames = "Genes")
  simulated_deconvoluted <- deconvolute("cibersortx", simulated_data, signature, perm = 1) |>
    pivot_longer(cols = -cell_type, names_to = "sample", values_to = "cell_count")

  cell_fractions <- simulation$cell_fractions |>
    as_tibble(rownames = "sample") |>
    pivot_longer(cols = -sample, names_to = "cell_type", values_to = "cell_count")

  joined_df <- dplyr::inner_join(cell_fractions, simulated_deconvoluted, by = c("cell_type", "sample"), suffix = c("_cf", "_sd"))
  joined_df
}

calculate_simulation_metrics <- function(joined_df) {
  rmse <- sqrt(mean((joined_df$cell_count_cf - joined_df$cell_count_sd)^2))

  correlation <- cor(joined_df$cell_count_cf, joined_df$cell_count_sd, method = "pearson")
  correlation2 <- cor(joined_df$cell_count_cf, joined_df$cell_count_sd, method = "spearman")

  metrics <- list(
    rmse = rmse,
    correlation_pearson = correlation,
    correlation_spearman = correlation2
  )
  metrics
}

run_sim2 <- function() {
  immune_cells_lowest_value_imputed <-
    read_excel(here("data/raw_external/Rieckmann_2017_S3.xlsx"), na = "NA") |>
    rename(Genes = `Gene names`) |>
    select(Genes, starts_with("LFQ.intensity_")) |> # also available: LFQ.Intensity.Imputed
    update_gene_symbols(verbose = TRUE) |>
    handle_gene_groups() |>
    handle_missing_values(imputation_mode = "lowest_value") |>
    handle_duplicate_genes(duplicate_mode = "slice") |>
    select(-contains(c("Erythrocyte", "Thrombocyte", "Th", "Tregs")))

  cell_types <- gsub("LFQ.intensity(.imputed)?_(.*)_\\d+_(activated|steady-state)", "\\2_\\3", colnames(immune_cells_lowest_value_imputed |> select(-Genes)))
  cell_type_mapping <- read_csv(here("data/signature_matrices/source/cell_type_mapping.csv"))
  cell_type_mapping <- setNames(cell_type_mapping$group, cell_type_mapping$cell_type)
  cell_types <- dplyr::recode(cell_types, !!!cell_type_mapping)

  signature <- read_tsv(
    here("data/signature_matrices/LM7c.txt")
  )

  sim_df <- get_simulation_metrics2(immune_cells_lowest_value_imputed, cell_types, signature, seed = 4)
  calculate_simulation_metrics(sim_df)
}

db <- memoise::cache_filesystem(tools::R_user_dir("proteoDeconv", "cache"))
get_simulation_metrics <<- memoise::memoise(get_simulation_metrics, cache = db)
get_simulation_metrics2 <<- memoise::memoise(get_simulation_metrics2, cache = db)
