# Shared constants and helpers for SCZ pathway enrichment scripts.
# Source from any pipeline script with: source("R/shared.R")

# ---- Cell-class groupings (SEA-AD subclass level) --------------------------
EXC <- c("L2_3 IT","L4 IT","L5 IT","L5 ET","L5_6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")
INH <- c("Lamp5","Pax6","Pvalb","Sncg","Sst","Sst Chodl","Vip","Chandelier")
GLI <- c("Astro","Oligo","OPC","Micro-PVM","Endo","VLMC")

# ---- Palettes (Okabe-Ito + custom class colours) ---------------------------
COL_UP     <- "#D55E00"   # vermillion
COL_DOWN   <- "#0072B2"   # blue
COL_NS     <- "grey75"
CLASS_COL  <- c(Excitatory = "#117733",
                Inhibitory = "#882255",
                Glia       = "#DDCC77")

# ---- Helpers ---------------------------------------------------------------
order_ct <- function(cts) {
  ord <- c(intersect(EXC, cts), intersect(INH, cts), intersect(GLI, cts))
  factor(cts, levels = ord)
}

class_of <- function(ct) dplyr::case_when(
  ct %in% EXC ~ "Excitatory",
  ct %in% INH ~ "Inhibitory",
  ct %in% GLI ~ "Glia",
  TRUE        ~ NA_character_
)

# Significance marker from padj (used widely in heatmaps)
sig_mark <- function(p) dplyr::case_when(
  p < 0.001 ~ "***",
  p < 0.01  ~ "**",
  p < 0.05  ~ "*",
  p < 0.1   ~ "+",
  TRUE      ~ ""
)

# ---- Default file paths (override at top of scripts if needed) -------------
INPUT_DE_CSV <- "data/DE_genes_all_cells_scz.csv"
GWAS_NO_MHC  <- "~/Github/scz_cell_type_enrichment/data/gwas/scz_gwas_gene_set_no_mhc.csv"
GWAS_WITH_MHC <- "~/Github/scz_cell_type_enrichment/data/gwas/scz_gwas_gene_set.csv"

RESULTS_TABLES  <- "results/tables"
RESULTS_FIGURES <- "results/figures"
