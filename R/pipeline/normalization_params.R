normalization_params <- tibble::tibble(
  normalization = c("none", "cycloess", "quantile", "vsn")
)

norm_param_map_obj <- tar_map(
  values = normalization_params,
  tar_target(
    norm_dda_data,
    preprocess_protein_data(
      input = dda_pg,
      method = "maxquant",
      imputation_mode = "lowest_value",
      duplicate_mode = "slice",
      remove_samples_pattern = "Erythrocyte|Thrombocyte|Library",
      update_symbols = TRUE,
      unlog = ifelse(normalization != "none", TRUE, FALSE),
      tpm = TRUE,
      normalize = normalization
    )
  ),
  tar_target(
    norm_sig_dda,
    create_signature_matrix(
      refsample = norm_dda_data,
      phenoclasses = create_general_phenoclasses(norm_dda_data)
    )
  ),
  tar_target(
    norm_simulation_dda_random,
    simulate_data(
      data = norm_dda_data,
      cell_types = map_cell_groups(colnames(norm_dda_data |> dplyr::select(-Genes))),
      seed = 4,
      scenario = "random"
    )
  ),
  tar_target(
    norm_deconvoluted_dda_random,
    deconvolute_data(
      method = "cibersort",
      preprocessed_data = norm_simulation_dda_random$simulated_data,
      signature_df = norm_sig_dda
    )
  ),
  tar_target(
    norm_metrics_dda_random,
    calculate_metrics(
      norm_simulation_dda_random$cell_fractions,
      norm_deconvoluted_dda_random,
      normalization = normalization
    )
  ),
  tar_target(
    norm_deconvoluted_pure_dda,
    deconvolute_data(
      method = "cibersort",
      preprocessed_data = norm_dda_data,
      signature_df = norm_sig_dda,
      normalization = normalization
    )
  ),
  tar_target(
    norm_dia_data,
    preprocess_protein_data(
      input = dia_pg,
      method = "diann",
      imputation_mode = "lowest_value",
      duplicate_mode = "slice",
      remove_samples_pattern = "NK",
      unlog = FALSE,
      tpm = TRUE,
      normalize = normalization
    )
  ),
  tar_target(
    norm_simulation_dia_random,
    simulate_data(
      data = norm_dia_data,
      cell_types = map_cell_groups(colnames(norm_dia_data |> dplyr::select(-Genes))),
      seed = 4,
      scenario = "random"
    )
  ),
  tar_target(
    norm_deconvoluted_dia_random,
    deconvolute_data(
      method = "cibersort",
      preprocessed_data = norm_simulation_dia_random$simulated_data,
      signature_df = norm_sig_dda
    )
  ),
  tar_target(
    norm_metrics_dia_random,
    calculate_metrics(
      norm_simulation_dia_random$cell_fractions,
      norm_deconvoluted_dia_random,
      normalization = normalization
    )
  ),
  tar_target(
    norm_deconvoluted_pure_dia,
    deconvolute_data(
      method = "cibersort",
      preprocessed_data = norm_dia_data,
      signature_df = norm_sig_dda,
      normalization = normalization
    )
  )
)


normalization_analysis <- tar_plan(
  norm_param_map_obj,
  tar_combine(
    combined_norm_metrics_dda_random,
    norm_param_map_obj[["norm_metrics_dda_random"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  ),
  tar_combine(
    combined_norm_metrics_dia_random,
    norm_param_map_obj[["norm_metrics_dia_random"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  ),
  tar_combine(
    combined_norm_deconvoluted_pure_dda,
    norm_param_map_obj[["norm_deconvoluted_pure_dda"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  ),
  tar_combine(
    combined_norm_deconvoluted_pure_dia,
    norm_param_map_obj[["norm_deconvoluted_pure_dia"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  )
)
