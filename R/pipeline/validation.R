sim_validation_datasets <- tidyr::expand_grid(
  dataset_id = c("dia_pg", "dia_ug"),
  imp_mode = c("lowest_value"),
  dup_mode = c("slice", "merge"),
  upd_symbols = c(FALSE, TRUE),
  unlog = FALSE,
  tpm = TRUE
)


sim_validation_obj <- tar_map(
  values = sim_validation_datasets,
  tar_target(
    sim_val_preprocessed_for_simulation,
    preprocess_protein_data(
      input = file[[dataset_id]],
      method = ifelse(str_detect(dataset_id, "dia"), "diann", "maxquant"),
      imputation_mode = imp_mode,
      duplicate_mode = dup_mode,
      update_symbols = upd_symbols,
      unlog = unlog,
      tpm = tpm,
      remove_samples_pattern = if (str_detect(dataset_id, "dia"))
        "PBMC|Mix|NK" else "Erythrocyte|Thrombocyte|Library"
    )
  ),
  tar_target(
    sim_val_preprocessed_for_deconvolution,
    preprocess_protein_data(
      input = file[[dataset_id]],
      method = ifelse(str_detect(dataset_id, "dia"), "diann", "maxquant"),
      imputation_mode = imp_mode,
      duplicate_mode = dup_mode,
      update_symbols = upd_symbols,
      unlog = unlog,
      tpm = tpm,
      remove_samples_pattern = if (str_detect(dataset_id, "dia")) "NK" else
        "Erythrocyte|Thrombocyte|Library"
    )
  ),
  tar_target(
    sim_val_simulation,
    simulate_data(
      data = sim_val_preprocessed_for_simulation,
      cell_types = map_cell_groups(colnames(
        sim_val_preprocessed_for_simulation
      )),
      seed = 4,
      scenario = "even",
      whitelist = c("Monocytes", "CD8+ T cells")
    )
  ),
  tar_target(
    sim_val_deconvoluted_sim,
    deconvolute_data(
      method = "cibersort",
      preprocessed_data = sim_val_simulation$simulated_data,
      signature_df = sig_dda,
      return_df = TRUE,
      imp_mode = imp_mode,
      dup_mode = dup_mode,
      dataset_id = dataset_id
    )
  ),
  tar_target(
    sim_val_deconvoluted_pure,
    deconvolute_data(
      method = "cibersort",
      preprocessed_data = sim_val_preprocessed_for_deconvolution,
      signature_df = sig_dda,
      return_df = TRUE,
      imp_mode = imp_mode,
      dup_mode = dup_mode,
      dataset_id = dataset_id
    )
  ),
  tar_target(
    sim_val_deconvoluted_both,
    dplyr::bind_rows(
      sim_val_deconvoluted_pure,
      sim_val_deconvoluted_sim,
      .id = "id"
    )
  )
)

sim_validation_analysis <- tar_plan(
  sim_validation_obj,
  tar_combine(
    sim_val_combined_deconv_both,
    sim_validation_obj[["sim_val_deconvoluted_both"]],
    command = dplyr::bind_rows(!!!.x, .id = "id")
  )
)
