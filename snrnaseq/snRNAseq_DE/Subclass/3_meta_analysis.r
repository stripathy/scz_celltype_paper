setwd("P1_SCZ_DE_fresh")
library(dplyr)
library(ggplot2)
library(tibble)
library(ggrepel)


library(dplyr)

# Get subclass names only
meta <- read.csv("Files/Pseudobulk_metadata_subclass_Mclean.csv")

colnames(meta) <- gsub("\\.", " ", colnames(meta))
cell_types <- colnames(meta)[1:24]
cell_types <- gsub("^L2 3", "L2_3", cell_types)
cell_types <- gsub("^L5 6", "L5_6", cell_types)
cell_types <- gsub("^Micro PVM", "Micro-PVM", cell_types)

# Read meta result files
files <- list.files(
  "/scratch/nendresz/P1_SCZ_DE_fresh/Files",
  pattern = "^meta_results_.*\\.csv$",
  full.names = TRUE
)

results <- list()

for (f in files) {
  ct <- basename(f) %>%
    sub("^meta_results_", "", x = .) %>%
    sub("\\.csv$", "", x = .)

  results[[ct]] <- read.csv(f)
}

# Keep subclasses only
results_subclass <- results[names(results) %in% cell_types]

names(results_subclass)


all_sig_genes <- list()

for (ct in names(results_subclass)) {
  df <- results_subclass[[ct]]

  sig_df <- df %>%
    filter(padj < 0.05)

  all_sig_genes[[ct]] <- sig_df$genes

  cat(ct, ":", nrow(sig_df), "sig genes\n")
}

sig_genes_unique <- sort(unique(unlist(all_sig_genes)))


write.csv(sig_genes_unique, "Files/sig_genes_all_cells_scz.csv")

all_sig_genes_w_estimates <- list()

for (ct in names(results_subclass)) {
  df <- results_subclass[[ct]]

  sig_df <- df %>%
    filter(padj < 0.05) %>%
    select(genes, estimate)

  all_sig_genes_w_estimates[[ct]] <- sig_df
}


saveRDS(all_sig_genes_w_estimates, "Files/sig_genes_all_cells_scz_w_estimates.RDS")



#ALL GENES ALL CELL TYPES


all_genes_df <- bind_rows(
  results_subclass,
  .id = "CellType"
)


write.csv(all_genes_df, "Files/DE_genes_all_cells_scz.csv", row.names = FALSE)


all_genes_df <- read.csv("Files/DE_genes_all_cells_scz.csv")

library(dplyr)

all_genes_df <- all_genes_df %>%
  mutate(
    CellType = gsub("Lamp5Lhx6", "Lamp5_Lhx6", CellType),
    cell_type = gsub("Lamp5Lhx6", "Lamp5_Lhx6", cell_type)
  )

write.csv(all_genes_df, "Files/DE_genes_all_cells_scz.csv", row.names = FALSE)


write.csv(all_genes_df,"/scratch/nendresz/scz_celltype_paper/transcriptomic/data/figure_inputs/DE_genes_all_cells_scz.csv", row.names = FALSE)