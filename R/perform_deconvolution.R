deconvolute_data <- function(
  method,
  preprocessed_data,
  signature_df,
  return_df = FALSE,
  ...
) {
  method_lower <- tolower(method)

  if (method_lower == "epic_noothercells") {
    result <- deconvolute(
      "epic",
      preprocessed_data,
      signature_df,
      with_other_cells = FALSE
    )
  } else {
    result <- deconvolute(
      method,
      preprocessed_data,
      signature_df
    )
  }
  if (!return_df) {
    return(result)
  } else {
    result <- result |>
      tibble::as_tibble(rownames = "sample") |>
      tidyr::pivot_longer(
        cols = -sample,
        names_to = "cell_type",
        values_to = "cell_count"
      ) %>%
      dplyr::mutate(method = method) %>%
      dplyr::mutate(!!!list(...))
    return(result)
  }
}
