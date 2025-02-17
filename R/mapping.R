prepare_diann_matrix <- function(data) {
  data %>%
    select(Genes, starts_with("D:")) %>%
    rename_with(~ str_extract(.x, "(?<=\\\\)[^\\\\]+(?=\\.(mzML|raw))|Genes"), everything()) %>%
    select(Genes, na.omit(colnames(.)))
}

mapping_list <- list(
  "B cells" = c("\\bB\\b", "Bcell", "Bnaive", "_B\\.(memory|naive|plasma)"),
  "CD4+ T cells" = c("\\bT4\\b", "T4\\.(naive|CM|EM|EMRA)"), # "Tregs", "nTregs", "mTregs", "Th1", "Th2", "Th17",
  "CD8+ T cells" = c("\\bT8\\b", "CD8T", "T8\\.(naive|CM|EM|EMRA)"),
  "Dendritic cells" = c("mDC", "pDC", "DC"),
  "Monocytes" = c("MO\\.classical", "MO\\.intermediate", "MO\\.nonclassical", "Mono"),
  "NK cells" = c("NK\\.bright", "NK\\.dim", "NK", "NKactive"),
  "Granulocytes" = c("Basophil", "Eosinophil", "Neutrophil"),
  "Erythrocytes" = c("Erythrocyte"),
  "Thrombocytes" = c("Thrombocyte")
)

map_cell_groups <- function(column_names, mapping_rules = NULL, default_group = "Unknown", verbose = FALSE) {
  if (is.null(mapping_rules)) {
    mapping_rules <- mapping_list
  }
    proteoDeconv::map_cell_groups(column_names, mapping_rules, default_group = default_group, verbose = verbose)
}


create_general_phenoclasses <- function(immune_cells,
                                        mapping_rules = NULL,
                                        verbose = FALSE) {
  if (is.null(mapping_rules)) {
    mapping_rules <- mapping_list
  }
  proteoDeconv::create_phenoclasses(immune_cells, mapping_rules, verbose = verbose)
}
