preprocess_protein_data <- function(
  input,
  method = "diann",
  imputation_mode = "lowest_value",
  duplicate_mode = "slice",
  remove_samples_pattern = NULL,
  unlog = FALSE,
  tpm = TRUE,
  update_symbols = FALSE,
  normalize = "none"
) {
  if (is.character(input) && length(input) == 1) {
    df <- readr::read_tsv(input, na = c("", "NA", "0"))
  } else if (is.data.frame(input)) {
    df <- input
  } else {
    stop("`input` must be either a file path (string) or a data frame.")
  }

  if (!is.null(remove_samples_pattern)) {
    df <- dplyr::select(df, -dplyr::matches(remove_samples_pattern))
  }

  if (method == "diann") {
    df <- prepare_diann_matrix(df)
  } else if (method == "maxquant") {
    colnames(df) <- make.names(colnames(df))
    df <- df %>%
      dplyr::rename(Genes = `Gene.names`) %>%
      dplyr::select(Genes, dplyr::starts_with("LFQ"))
  }

  mat <- tibble_to_matrix(df)

  mat <- extract_identifiers(mat)

  if (update_symbols) {
    mat <- update_gene_symbols(mat, verbose = FALSE)
  }

  if (imputation_mode == "RF") {
    mat <- handle_missing_values(
      mat,
      imputation_mode = imputation_mode,
      ntree = 100
    )
  } else {
    mat <- handle_missing_values(mat, imputation_mode = imputation_mode)
  }

  mat <- handle_duplicates(mat, duplicate_mode = duplicate_mode)

  if (normalize != "none") {
    expr_mat <- log2(mat)
    normalized_mat <- switch(
      normalize,
      cycloess = limma::normalizeCyclicLoess(expr_mat, method = "fast"),
      quantile = limma::normalizeQuantiles(expr_mat),
      vsn = vsn::justvsn(expr_mat),
      stop("Unknown normalization method: ", normalize)
    )
    rownames(normalized_mat) <- rownames(expr_mat)
    colnames(normalized_mat) <- colnames(expr_mat)
    mat <- normalized_mat
  }

  mat <- handle_scaling(mat, unlog = unlog, tpm = tpm)

  return(mat)
}

tibble_to_matrix <- function(data) {
  if (!is.data.frame(data)) {
    stop("Input must be a data frame or tibble")
  }
  char_columns <- vapply(data, is.character, logical(1))
  if (sum(char_columns) != 1) {
    stop(
      "Data must contain exactly one character column, but found ",
      sum(char_columns)
    )
  }
  rowname_col <- names(data)[char_columns]
  row_names <- data[[rowname_col]]
  data_rest <- data[setdiff(names(data), rowname_col)]
  m <- as.matrix(data_rest)
  rownames(m) <- row_names
  m
}
