sig_source_params <- tidyr::expand_grid(
    data_source       = c("dda_pg", "dia_pg"),
    signature_source  = c("dda", "sc_rna_seq_pbmc")
)

mapping_sig_comparison <- list(
    "B cells" = c("\\bB\\b", "Bcell", "Bnaive", "_B\\.(memory|naive|plasma)"),
    "CD4+ T cells" = c("\\bT4\\b", "CD4", "T4\\.(naive|CM|EM|EMRA)"),
    "CD8+ T cells" = c("\\bT8\\b", "CD8", "T8\\.(naive|CM|EM|EMRA)"),
    "Monocytes" = c("MO\\.classical", "MO\\.intermediate", "MO\\.nonclassical", "Mono"),
    "NK cells" = c("NK\\.bright", "NK\\.dim", "NK\\b", "NKactive")
)

sig_source_obj <- tar_map(
    values = sig_source_params,
    tar_target(
        sig_source_preprocessed_data,
        preprocess_protein_data(
            input = file[[data_source]],
            method = ifelse(str_detect(data_source, "dia"), "diann", "maxquant"),
            imputation_mode = "lowest_value",
            duplicate_mode = "slice",
            update_symbols = TRUE,
            unlog = FALSE,
            tpm = TRUE,
            remove_samples_pattern = if (str_detect(data_source, "dia")) {
                "NK"
            } else {
                "Erythrocyte|Thrombocyte|Library"
            }
        )
    ),
    tar_target(
        sig_source_simulation,
        simulate_data(
            data = sig_source_preprocessed_data,
            cell_types = map_cell_groups(
                colnames(sig_source_preprocessed_data %>% select(-Genes)),
                mapping_rules = mapping_sig_comparison
            ),
            blacklist = "Unknown",
            seed = 4
        )
    ),
    tar_target(
        signature_mat,
        if (signature_source == "dda") {
            create_signature_matrix(
                refsample = dda_data,
                phenoclasses = create_general_phenoclasses(
                    dda_data,
                    mapping_rules = mapping_sig_comparison
                ),
                g_min = 300,
                g_max = 500,
                q_value = 0.01,
                filter = TRUE
            )
        } else {
            create_signature_matrix(
                refsample = sc_rna_seq_pbmc_data,
                phenoclasses = NULL,
                single_cell = TRUE,
                g_min = 300,
                g_max = 500,
                q_value = 0.01,
                filter = TRUE,
                verbose = TRUE
            )
        }
    ),
    tar_target(
        sig_source_deconvoluted_sim,
        deconvolute_data(
            method            = "cibersort",
            preprocessed_data = sig_source_simulation$simulated_data,
            signature_df      = signature_mat
        )
    ),
    tar_target(
        sig_source_metrics,
        calculate_metrics(
            sig_source_simulation$cell_fractions,
            sig_source_deconvoluted_sim,
            data_source = data_source,
            signature_source = signature_source,
            by_cell = TRUE
        )
    ),
    tar_target(
        sig_source_deconvoluted_pure,
        deconvolute_data(
            method = "cibersort",
            preprocessed_data = sig_source_preprocessed_data,
            signature_df = signature_mat,
            data_source = data_source,
            signature_source = signature_source
        )
    )
)

sig_source_analysis <- tar_plan(
    tar_target(
        sc_rna_seq_pbmc_data,
        {
            df <- read_tsv(sc_rna_seq_pbmc, name_repair = "minimal")
            df[, !grepl("NKT cells", names(df))]
        }
    ),
    sig_source_obj,
    tar_combine(
        sig_source_combined_metrics,
        sig_source_obj[["sig_source_metrics"]],
        command = bind_rows(!!!.x, .id = "id")
    ),
    tar_combine(
        sig_source_combined_deconvoluted_pure,
        sig_source_obj[["sig_source_deconvoluted_pure"]],
        command = bind_rows(!!!.x, .id = "id")
    )
)
