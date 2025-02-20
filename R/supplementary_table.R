make_supplementary_table <- function() {
  mean_metrics2 <- function(data, group_var = id) {
    data %>%
      group_by({{ group_var }}) %>%
      summarise(
        across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
        across(where(negate(is.numeric)), ~ first(.x)),
        .groups = "drop"
      ) |>
      select(-cell_type)
  }

  supplementary_table <- list(
    `Algorithms, per CT` = tar_read(combined_algo_metrics_dda),
    `Algorithms, mean` = tar_read(combined_algo_metrics_dda) |> 
      mean_metrics2() |> 
      arrange(desc(correlation_pearson)),
    
    `Preprocessing, per CT` = tar_read(combined_preprocessing_sim_metrics) |> 
      filter(imp_mode == "lowest_value") |> 
      select(-imp_mode) |> 
      filter(!str_detect(dataset_id, "dia")),
    `Preprocessing, mean` = tar_read(combined_preprocessing_sim_metrics) |> 
      filter(imp_mode == "lowest_value") |> 
      select(-imp_mode) |> 
      filter(!str_detect(dataset_id, "dia")) |> 
      mean_metrics2() |> 
      arrange(desc(correlation_pearson)),
    
    `Normalization, per CT` = tar_read(combined_norm_metrics_dda_random),
    `Normalization, mean` = tar_read(combined_norm_metrics_dda_random) |> 
      mean_metrics2() |> 
      arrange(desc(correlation_pearson)),
    
    `Imputation, per CT` = tar_read(combined_preprocessing_sim_metrics) |> 
      filter(dup_mode == "slice", upd_symbols == TRUE) |> 
      select(-c(dup_mode, upd_symbols)) |> 
      filter(!str_detect(dataset_id, "dia")),
    `Imputation, mean` = tar_read(combined_preprocessing_sim_metrics) |> 
      filter(dup_mode == "slice", upd_symbols == TRUE) |> 
      select(-c(dup_mode, upd_symbols)) |> 
      filter(!str_detect(dataset_id, "dia")) |> 
      mean_metrics2() |> 
      arrange(desc(correlation_pearson)),
    
    `Signature params, per CT` = tar_read(combined_sig_params_metrics_dda) |> 
      filter(str_detect(method_id, "-")),
    `Signature params, mean` = tar_read(combined_sig_params_metrics_dda) |> 
      filter(str_detect(method_id, "-")) |> 
      mean_metrics2() |> 
      arrange(desc(correlation_pearson))
  )

  file_path <- here("supplementary_table.xlsx")
  write_xlsx(supplementary_table, file_path)

  return(file_path)
}