signature_params_dda <- tibble::tibble(
  method_id = c(
    "50-200",
    "100-250",
    "200-400",
    "350-550",
    "550-750",
    "750-950"
  ),
  g_min = c(50, 100, 200, 350, 550, 750),
  g_max = c(200, 250, 400, 550, 750, 950),
  q_value = rep(0.01, 6),
  filter = rep(FALSE, 6)
)

sig_param_map_obj <- tar_map(
  values = signature_params_dda,
  names = method_id,
  tar_target(
    signature_params_sig_dda,
    create_signature_matrix(
      refsample = dda_data,
      phenoclasses = create_general_phenoclasses(dda_data),
      g_min = g_min,
      g_max = g_max,
      q_value = q_value,
      filter = filter,
      verbose = TRUE
    )
  ),
  tar_target(
    signature_params_deconvoluted_dda,
    deconvolute_data(
      method = "cibersort",
      preprocessed_data = signature_params_simulation_dda$simulated_data,
      signature_df = signature_params_sig_dda
    )
  ),
  tar_target(
    signature_params_metrics_dda,
    calculate_metrics(
      signature_params_simulation_dda$cell_fractions,
      signature_params_deconvoluted_dda,
      method_id = method_id,
      g_min = g_min,
      g_max = g_max,
      q_value = q_value,
      filter = filter
    )
  ),
  tar_target(
    signature_params_deconvoluted_pure_dda,
    deconvolute_data(
      method = "cibersort",
      preprocessed_data = dda_data,
      signature_df = signature_params_sig_dda,
      return_df = TRUE,
      method_id = method_id,
      g_min = g_min,
      g_max = g_max,
      q_value = q_value,
      filter = filter
    )
  ),
  tar_target(
    signature_params_deconvoluted_dia,
    deconvolute_data(
      method = "cibersort",
      preprocessed_data = signature_params_simulation_dia$simulated_data,
      signature_df = signature_params_sig_dda
    )
  ),
  tar_target(
    signature_params_metrics_dia,
    calculate_metrics(
      signature_params_simulation_dia$cell_fractions,
      signature_params_deconvoluted_dia,
      method_id = method_id,
      g_min = g_min,
      g_max = g_max,
      q_value = q_value,
      filter = filter
    )
  ),
  tar_target(
    signature_params_deconvoluted_pure_dia,
    deconvolute_data(
      method = "cibersort",
      preprocessed_data = dia_data_full,
      signature_df = signature_params_sig_dda,
      return_df = TRUE,
      method_id = method_id,
      g_min = g_min,
      g_max = g_max,
      q_value = q_value,
      filter = filter
    )
  ),
  tar_target(
    signature_params_deconvoluted_dia_reduced,
    deconvolute_data(
      method = "cibersort",
      preprocessed_data = signature_params_simulation_dia_reduced$simulated_data,
      signature_df = signature_params_sig_dda
    )
  ),
  tar_target(
    signature_params_metrics_dia_reduced,
    calculate_metrics(
      signature_params_simulation_dia_reduced$cell_fractions,
      signature_params_deconvoluted_dia_reduced,
      method_id = method_id,
      g_min = g_min,
      g_max = g_max,
      q_value = q_value,
      filter = filter
    )
  ),
  tar_target(
    signature_params_deconvoluted_pure_dia_reduced,
    deconvolute_data(
      method = "cibersort",
      preprocessed_data = dia_data_reduced,
      signature_df = signature_params_sig_dda,
      return_df = TRUE,
      method_id = method_id,
      g_min = g_min,
      g_max = g_max,
      q_value = q_value,
      filter = filter
    )
  )
)

signature_params_analysis <- tar_plan(
  tar_target(
    signature_params_simulation_dda,
    simulate_data(
      data = dda_data,
      cell_types = map_cell_groups(colnames(dda_data)),
      seed = 4
    )
  ),
  tar_target(
    signature_params_simulation_dia,
    simulate_data(
      data = dia_data,
      cell_types = map_cell_groups(colnames(dia_data)),
      seed = 4
    )
  ),
  tar_target(
    signature_params_simulation_dia_reduced,
    simulate_data(
      data = dia_data_reduced,
      cell_types = map_cell_groups(colnames(dia_data_reduced)),
      seed = 4
    )
  ),
  sig_param_map_obj,
  tar_combine(
    combined_sig_params_metrics_dda,
    sig_param_map_obj[["signature_params_metrics_dda"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  ),
  tar_combine(
    combined_sig_params_metrics_dia,
    sig_param_map_obj[["signature_params_metrics_dia"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  ),
  tar_combine(
    combined_sig_params_deconvoluted_pure_dda,
    sig_param_map_obj[["signature_params_deconvoluted_pure_dda"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  ),
  tar_combine(
    combined_sig_params_deconvoluted_pure_dia,
    sig_param_map_obj[["signature_params_deconvoluted_pure_dia"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  ),
  tar_combine(
    combined_sig_params_metrics_dia_reduced,
    sig_param_map_obj[["signature_params_metrics_dia_reduced"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  ),
  tar_combine(
    combined_sig_params_deconvoluted_pure_dia_reduced,
    sig_param_map_obj[["signature_params_deconvoluted_pure_dia_reduced"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  )
)
