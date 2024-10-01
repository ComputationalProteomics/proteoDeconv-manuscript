library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(cowplot)
library(patchwork)

pkgload::load_all("~/Documents/proteoDeconv")
library(here)
library(tidyverse)
library(readxl)
library(SimBu)

here::i_am("proteoDeconv-manuscript.Rproj")

source(here("notebooks/simulation.R"))

scale_colour_brewer_d <- function(..., palette = "Dark2") {
  scale_colour_brewer(..., palette = palette)
}

scale_fill_brewer_d <- function(..., palette = "Dark2") {
  scale_fill_brewer(..., palette = palette)
}

Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/usr/local/bin", sep = ":")) # This appears to be needed when running in renv

options(
  ggplot2.discrete.colour = scale_colour_brewer_d,
  ggplot2.discrete.fill = scale_fill_brewer_d
)
theme_set(
  theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "bottom") +
    theme(
      strip.text = element_text(size = rel(1)),
      strip.background = element_rect(fill = "white", colour = "white", linewidth = 0)
    )
)

immunedeconv::set_cibersort_binary(here("data/CIBERSORT.R"))
immunedeconv::set_cibersort_mat(here("data/signature_matrices/LM22.txt"))



rename_cell_types <- function(cell_type) {
  case_match(cell_type,
    "MO" ~ "Monocytes",
    "DC" ~ "Dendritic cells",
    "T cells CD4" ~ "CD4+ T cells",
    "T cells CD8" ~ "CD8+ T cells",
    "T4 cells" ~ "CD4+ T cells",
    "T8 cells" ~ "CD8+ T cells",
    "NK T cells" ~ "NK cells",
    "otherCells" ~ "Uncharacterized cells",
    .default = cell_type
  )
}

rename_samples <- function(sample_part) {
  case_match(sample_part,
    "Bcell" ~ "B cells",
    "Bnaive" ~ "B naive cells",
    "CD8T" ~ "CD8+ T cells",
    "Mix" ~ "Mix",
    "Mono" ~ "Monocytes",
    "NKactive" ~ "NK cells",
    "Tcell" ~ "T cells",
    "DC" ~ "Dendritic cells",
    .default = sample_part
  )
}

prepare_diann_matrix <- function(data) {
  data %>%
    select(Genes, starts_with("D:")) %>%
    rename_with(~ str_extract(.x, "(?<=\\\\)[^\\\\]+(?=\\.(mzML|raw))|Genes"), everything()) %>%
    select(Genes, na.omit(colnames(.)))
}
