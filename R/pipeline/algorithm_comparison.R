algorithm_params <- tibble::tibble(
  algorithm_id = c(
    "cibersort",
    "cibersortx",
    "epic",
    "epic_noOtherCells",
    "bayesdebulk"
  ),
  algorithm = c(
    "cibersort",
    "cibersortx",
    "epic",
    "epic_noOtherCells",
    "bayesdebulk"
  )
)
algo_param_map_dda <- tar_map(
  values = algorithm_params,
  names = algorithm_id,
  tar_target(
    algo_deconvoluted_dda,
    deconvolute_data(
      method = algorithm,
      preprocessed_data = algo_simulation_dda$simulated_data,
      signature_df = algo_sig_dda
    )
  ),
  tar_target(
    algo_metrics_dda,
    calculate_metrics(
      algo_simulation_dda$cell_fractions,
      algo_deconvoluted_dda,
      algorithm_id = algorithm_id
    )
  ),
  tar_target(
    algo_deconvoluted_pure_dda,
    deconvolute_data(
      method = algorithm,
      preprocessed_data = dda_data,
      signature_df = algo_sig_dda,
      return_df = TRUE,
      algorithm_id = algorithm_id
    )
  )
)

algo_param_map_dia <- tar_map(
  values = algorithm_params,
  names = algorithm_id,
  tar_target(
    algo_deconvoluted_dia,
    deconvolute_data(
      method = algorithm,
      preprocessed_data = algo_simulation_dia$simulated_data,
      signature_df = algo_sig_dda
    )
  ),
  tar_target(
    algo_metrics_dia,
    calculate_metrics(
      algo_simulation_dia$cell_fractions,
      algo_deconvoluted_dia,
      algorithm_id = algorithm_id
    )
  ),
  tar_target(
    algo_deconvoluted_pure_dia,
    deconvolute_data(
      method = algorithm,
      preprocessed_data = dia_data_full,
      signature_df = algo_sig_dda,
      return_df = TRUE,
      algorithm_id = algorithm_id
    )
  )
)

algorithm_analysis <- tar_plan(
  tar_target(
    algo_simulation_dda,
    simulate_data(
      data = dda_data,
      cell_types = map_cell_groups(colnames(dda_data)),
      seed = 4
    )
  ),
  tar_target(
    algo_sig_dda,
    create_signature_matrix(
      refsample = dda_data,
      phenoclasses = create_general_phenoclasses(dda_data)
    )
  ),
  tar_target(
    algo_simulation_dia,
    simulate_data(
      data = dia_data,
      cell_types = map_cell_groups(colnames(dia_data)),
      seed = 4
    )
  ),
  algo_param_map_dda,
  algo_param_map_dia,
  tar_combine(
    combined_algo_metrics_dda,
    algo_param_map_dda[["algo_metrics_dda"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  ),
  tar_combine(
    combined_algo_metrics_dia,
    algo_param_map_dia[["algo_metrics_dia"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  ),
  tar_combine(
    combined_algo_deconvoluted_pure_dda,
    algo_param_map_dda[["algo_deconvoluted_pure_dda"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  ),
  tar_combine(
    combined_algo_deconvoluted_pure_dia,
    algo_param_map_dia[["algo_deconvoluted_pure_dia"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  )
)
