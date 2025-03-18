calculate_metrics <- function(
  cell_fractions,
  deconvoluted_data,
  by_cell = TRUE,
  ...
) {
  if (is.matrix(cell_fractions)) {
    cell_fractions <- as_tibble(cell_fractions, rownames = "sample")
  }

  if (is.matrix(deconvoluted_data)) {
    deconvoluted_data <- as_tibble(deconvoluted_data, rownames = "sample")
  }

  if (!"cell_type" %in% colnames(cell_fractions)) {
    cell_fractions_long <- cell_fractions %>%
      pivot_longer(
        cols = -sample,
        names_to = "cell_type",
        values_to = "cell_count"
      )
  } else {
    cell_fractions_long <- cell_fractions
  }

  if (!"cell_type" %in% colnames(deconvoluted_data)) {
    deconvoluted_data <- deconvoluted_data %>%
      pivot_longer(
        cols = -sample,
        names_to = "cell_type",
        values_to = "cell_count"
      )
  }

  joined_df <- inner_join(
    cell_fractions_long,
    deconvoluted_data,
    by = c("sample", "cell_type"),
    suffix = c("_cf", "_sd")
  )

  if (!by_cell) {
    rmse <- sqrt(mean(
      (joined_df$cell_count_cf - joined_df$cell_count_sd)^2,
      na.rm = TRUE
    ))
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
    metrics <- tibble(
      rmse = rmse,
      correlation_pearson = correlation_pearson,
      correlation_spearman = correlation_spearman
    )
  } else {
    metrics <- joined_df %>%
      group_by(cell_type) %>%
      summarise(
        rmse = sqrt(mean((cell_count_cf - cell_count_sd)^2, na.rm = TRUE)),
        correlation_pearson = cor(
          cell_count_cf,
          cell_count_sd,
          method = "pearson",
          use = "complete.obs"
        ),
        correlation_spearman = cor(
          cell_count_cf,
          cell_count_sd,
          method = "spearman",
          use = "complete.obs"
        ),
        .groups = "drop"
      )
  }

  metrics <- metrics %>%
    mutate(!!!list(...))

  return(metrics)
}
