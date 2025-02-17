calculate_metrics <- function(cell_fractions, deconvoluted_data, by_cell = TRUE, ...) {
  cell_fractions_long <- cell_fractions |>
    tidyr::pivot_longer(
      cols = -sample,
      names_to = "cell_type",
      values_to = "cell_count"
    )

  joined_df <- dplyr::inner_join(
    cell_fractions_long,
    deconvoluted_data,
    by = c("cell_type", "sample"),
    suffix = c("_cf", "_sd")
  )

  if (!by_cell) {
    rmse <- sqrt(mean((joined_df$cell_count_cf - joined_df$cell_count_sd)^2, na.rm = TRUE))
    correlation_pearson <- cor(
      joined_df$cell_count_cf,
      joined_df$cell_count_sd,
      method = "pearson",
      use = "complete.obs"
    )
    correlation_spearman <- cor(
      joined_df$cell_count_cf,
      joined_df$cell_count_sd,
      method = "spearman",
      use = "complete.obs"
    )
    metrics <- tibble::tibble(
      rmse = rmse,
      correlation_pearson = correlation_pearson,
      correlation_spearman = correlation_spearman
    )
  } else {
    metrics <- joined_df |>
      dplyr::group_by(cell_type) |>
      dplyr::summarise(
        rmse = sqrt(mean((cell_count_cf - cell_count_sd)^2, na.rm = TRUE)),
        correlation_pearson = cor(cell_count_cf, cell_count_sd, method = "pearson", use = "complete.obs"),
        correlation_spearman = cor(cell_count_cf, cell_count_sd, method = "spearman", use = "complete.obs"),
        .groups = "drop"
      )
  }

  metrics <- metrics %>%
    dplyr::mutate(!!!list(...))

  return(metrics)
}
