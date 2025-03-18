preprocessing_datasets <- tidyr::expand_grid(
  dataset_id = c(
    "dda_pg",
    "dia_pg",
    "dia_ug",
    "dia_pg_reduced",
    "dia_ug_reduced"
  ), # "dda_rerun_pg", "dda_table",
  imp_mode = c("lowest_value"),
  dup_mode = c("slice", "merge"),
  upd_symbols = c(FALSE, TRUE),
  unlog = FALSE,
  tpm = TRUE
)

additional_dda_pg <- tidyr::expand_grid(
  dataset_id = c("dda_pg", "dia_pg", "dia_pg_reduced"),
  imp_mode = c("knn", "zero", "min", "MinDet", "MinProb", "RF"),
  dup_mode = "slice",
  upd_symbols = TRUE,
  unlog = FALSE,
  tpm = TRUE
)

preprocessing_datasets <- dplyr::bind_rows(
  preprocessing_datasets,
  additional_dda_pg
)


preprocessing_obj <- tar_map(
  values = preprocessing_datasets,
  tar_target(
    preprocessed_data,
    preprocess_protein_data(
      input = file[[gsub("_reduced", "", dataset_id)]],
      method = ifelse(str_detect(dataset_id, "dia"), "diann", "maxquant"),
      imputation_mode = imp_mode,
      duplicate_mode = dup_mode,
      update_symbols = upd_symbols,
      unlog = unlog,
      tpm = tpm,
      remove_samples_pattern = case_when(
        str_detect(dataset_id, "reduced") ~ "Bnaive|Mix|NKactive|PBMC|Tcell",
        str_detect(dataset_id, "dia") ~ "NKactive",
        TRUE ~ "Erythrocyte|Thrombocyte|Library"
      )
    )
  ),
  tar_target(
    simulation,
    simulate_data(
      data = preprocessed_data,
      cell_types = map_cell_groups(colnames(preprocessed_data)),
      seed = 4
    )
  ),
  tar_target(
    deconvoluted,
    deconvolute_data(
      method = "cibersort",
      preprocessed_data = simulation$simulated_data,
      signature_df = sig_dda
    )
  ),
  tar_target(
    simulation_metrics,
    calculate_metrics(
      simulation$cell_fractions,
      deconvoluted,
      imp_mode = imp_mode,
      dup_mode = dup_mode,
      upd_symbols = upd_symbols,
      dataset_id = dataset_id
    )
  ),
  tar_target(
    deconvoluted_pure,
    deconvolute_data(
      method = "cibersort",
      preprocessed_data = preprocessed_data,
      signature_df = sig_dda,
      return_df = TRUE,
      imp_mode = imp_mode,
      dup_mode = dup_mode,
      upd_symbols = upd_symbols,
      dataset_id = dataset_id
    )
  )
)

preprocessing_analysis <- tar_plan(
  preprocessing_obj,
  tar_combine(
    combined_preprocessing_sim_metrics,
    preprocessing_obj[["simulation_metrics"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  ),
  tar_combine(
    combined_preprocessing_deconvoluted_pure,
    preprocessing_obj[["deconvoluted_pure"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  )
)
