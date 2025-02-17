run_epic <- function(preprocessed_data, signature_df, with_other_cells, method_label, ...) {
  signature_df_epic <- handle_scaling(signature_df, tpm = TRUE, unlog = FALSE)
  preprocessed_data_epic <- handle_scaling(preprocessed_data, tpm = TRUE, unlog = FALSE)
  
  preprocessed_data <- proteoDeconv::handle_input_data(preprocessed_data, as_tibble = FALSE)
  signature_df <- proteoDeconv::handle_input_data(signature_df, as_tibble = FALSE)
  preprocessed_data_epic <- proteoDeconv::handle_input_data(preprocessed_data_epic, as_tibble = FALSE)
  signature_df_epic <- proteoDeconv::handle_input_data(signature_df_epic, as_tibble = FALSE)
  
  epic_signature <- list(
    refProfiles = as.matrix(signature_df_epic),
    sigGenes = rownames(signature_df_epic)
  )
  
  epic_res <- EPIC::EPIC(
    preprocessed_data_epic,
    withOtherCells = with_other_cells,
    reference = epic_signature
  )
  
  result <- epic_res$mRNAProportions %>%
    tibble::as_tibble(rownames = "sample") %>%
    tidyr::pivot_longer(
      cols = -sample,
      names_to = "cell_type",
      values_to = "cell_count"
    ) %>%
    dplyr::mutate(method = method_label) %>%
    dplyr::mutate(!!!list(...))
  
  return(result)
}

deconvolute_data <- function(method, preprocessed_data, signature_df, ...) {
  method_lower <- tolower(method)
  
  if(method_lower %in% c("cibersort", "cibersortx")) {
    deconv_result <- deconvolute(method, preprocessed_data, signature_df) %>%
      tidyr::pivot_longer(
        cols = -cell_type,
        names_to = "sample",
        values_to = "cell_count"
      ) %>%
      dplyr::mutate(method = method) %>%
      dplyr::mutate(!!!list(...))
    return(deconv_result)
  } else if(method_lower == "bayesdebulk") {
    signature_df <- read_tsv("data/signature_matrices/LM7c.txt")
    signature_df_epic <- handle_scaling(signature_df, tpm = TRUE, unlog = FALSE)
    preprocessed_data_epic <- handle_scaling(preprocessed_data, tpm = TRUE, unlog = FALSE)
    
    preprocessed_data <- proteoDeconv::handle_input_data(preprocessed_data, as_tibble = FALSE)
    signature_df <- proteoDeconv::handle_input_data(signature_df, as_tibble = FALSE)
    preprocessed_data_epic <- proteoDeconv::handle_input_data(preprocessed_data_epic, as_tibble = FALSE)
    signature_df_epic <- proteoDeconv::handle_input_data(signature_df_epic, as_tibble = FALSE)
    
    return(run_bayesdebulk(preprocessed_data, signature_df, ...))
  } else if(method_lower == "epic") {
    return(run_epic(preprocessed_data, signature_df, TRUE, "EPIC", ...))
  } else if(method_lower == "epic_noothercells") {
    return(run_epic(preprocessed_data, signature_df, FALSE, "EPIC, without uncharacterized cells", ...))
  } else {
    stop("Method not supported")
  }
}

run_bayesdebulk <- function(data,
                            signature_matrix,
                            n_iter = 1000,
                            burn_in = 100,
                            ...) {
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(BayesDeBulk)
  library(stringr)


  get_signature_markers <- function(Y, signature_matrix) {
    if (length(Y) > 1) {
      genes <- unique(c(rownames(Y[[1]]), rownames(Y[[2]])))
    } else {
      genes <- rownames(Y[[1]])
    }

    cell_types <- colnames(signature_matrix)
    index_matrix <- NULL

    for (s in seq_len(ncol(signature_matrix))) {
      for (k in seq_len(ncol(signature_matrix))) {
        if (k != s) {
          i <- (signature_matrix[, s] > 1000 &
            signature_matrix[, s] > 3 * signature_matrix[, k])
          marker_unique <- rownames(signature_matrix)[i]
          marker_unique <- marker_unique[!is.na(match(marker_unique, genes))]

          if (length(marker_unique) >= 1) {
            new_rows <- cbind(
              rep(cell_types[s], length(marker_unique)),
              rep(cell_types[k], length(marker_unique)),
              marker_unique
            )
            index_matrix <- if (is.null(index_matrix)) new_rows else rbind(index_matrix, new_rows)
          }
        }
      }
    }
    return(index_matrix)
  }

  markers <- get_signature_markers(list(data), signature_matrix)

  bayes <- BayesDeBulk(
    n.iter = n_iter,
    burn.in = burn_in,
    Y = list(data),
    markers = markers
  )

  bayes_cell_fractions <- bayes$cell.fraction

  transform_bayesdebulk_output <- function(bayes_cell_fractions, ...) {
    df <- as.data.frame(bayes_cell_fractions) %>%
      tibble::rownames_to_column("sample")

    df_long <- tidyr::pivot_longer(
      df,
      cols = -sample,
      names_to = "cell_type",
      values_to = "cell_count"
    ) %>%
      dplyr::mutate(method = "BayesDeBulk") %>%
      dplyr::mutate(!!!list(...))

    return(df_long)
  }

  result <- transform_bayesdebulk_output(bayes_cell_fractions, ...)
  return(result)
}
