preprocess_protein_data <- function(input,
                                    method = "diann",
                                    imputation_mode = "lowest_value",
                                    duplicate_mode = "slice",
                                    remove_samples_pattern = NULL,
                                    unlog = FALSE,
                                    tpm = TRUE,
                                    update_symbols = FALSE,
                                    normalize = "none") {
  if (is.character(input) && length(input) == 1) {
    df <- read_tsv(input, na = c("", "NA", "0"))
  } else if (is.data.frame(input)) {
    df <- input
  } else {
    stop("`input` must be either a file path (string) or a data frame.")
  }

  if (!is.null(remove_samples_pattern)) {
    df <- df %>% select(-matches(remove_samples_pattern))
  }


  if (method == "diann") {
    df <- prepare_diann_matrix(df)
  } else if (method == "maxquant") {
    colnames(df) <- make.names(colnames(df))
    df <- df |>
      rename(Genes = `Gene.names`) |>
      select(Genes, starts_with("LFQ"))
  }
  df <- handle_gene_groups(df)

  if (update_symbols) {
    df <- update_gene_symbols(df, verbose = FALSE)
  }

  if (imputation_mode == "RF") {
    df <- handle_missing_values(df, imputation_mode = imputation_mode, ntree = 10)
  } else {
    df <- handle_missing_values(df, imputation_mode = imputation_mode)
  }
  df <- handle_duplicate_genes(df, duplicate_mode = duplicate_mode)

  if (normalize != "none") {
    expr_mat <- as.matrix(df %>% select(-Genes))
    expr_mat <- log2(expr_mat)
    normalized_mat <- switch(normalize,
      "cycloess" = {
        limma::normalizeCyclicLoess(expr_mat, method = "fast")
      },
      "quantile" = {
        limma::normalizeQuantiles(expr_mat)
      },
      "vsn" = {
        vsn::justvsn(expr_mat)
      },
      stop("Unknown normalization method: ", normalize)
    )
    df <- bind_cols(
      Genes = df$Genes,
      as_tibble(normalized_mat, .name_repair = "minimal")
    )
    colnames(df)[-1] <- colnames(expr_mat)
  }

  df <- handle_scaling(df, unlog = unlog, tpm = tpm)

  return(df)
}
