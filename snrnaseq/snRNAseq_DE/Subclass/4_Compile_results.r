setwd("scz_celltype_paper/snrnaseq/snRNAseq_DE")
library(dplyr)
library(ggplot2)
library(tibble)
library(ggrepel)


# Get subclass names
meta <- read.csv("Files/Pseudobulk_metadata_subclass_McLean.csv")
colnames(meta) <- gsub("\\.", " ", colnames(meta))

# every column that is not donor metadata (see 2_DE.r)
cell_types <- setdiff(colnames(meta), c("Donor", "Age", "Sex", "Diagnosis", "PMI"))
cell_types <- gsub("^L2 3", "L2_3", cell_types)
cell_types <- gsub("^L5 6", "L5_6", cell_types)
cell_types <- gsub("^Micro PVM", "Micro-PVM", cell_types)


#Cohort-specific results
cohorts <- c("Batiuk", "Frohlich", "MSSM1", "McLean", "MSSM2", "HBCC", "Multiome")
cohort_results <- list()

for (ct in cell_types) {
  safe_type <- gsub("/", "_", ct)

  for (cohort in cohorts) {
    file <- paste0("Files/DE_results_", cohort, "_", safe_type, ".rds")
    if (!file.exists(file)) next

    df <- readRDS(file) %>%
      mutate(
        cell_type = ct,
        cohort = cohort
      )

    cohort_results[[paste(cohort, ct, sep = "_")]] <- df
  }
}

cohort_genes_df <- bind_rows(cohort_results) %>%
  mutate(cell_type = gsub("Lamp5Lhx6", "Lamp5_Lhx6", cell_type))

write.csv(cohort_genes_df, "Files/DE_genes_all_cohorts_subclass.csv", row.names = FALSE)


# 7-cohort meta-analysis
files <- list.files(
  "/scratch/nendresz/scz_celltype_paper/snrnaseq/snRNAseq_DE/Files",
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

# All genes, all subclass cell types
all_genes_df <- bind_rows(
  results_subclass,
  .id = "CellType"
) %>%
  mutate(
    CellType = gsub("Lamp5Lhx6", "Lamp5_Lhx6", CellType),
    cell_type = gsub("Lamp5Lhx6", "Lamp5_Lhx6", cell_type)
  )

  all_genes_df <- all_genes_df %>% select(-X)

write.csv(all_genes_df, "Files/DE_genes_all_cells_scz.csv", row.names = FALSE)

# Write for DE figure
write.csv(all_genes_df, "/scratch/nendresz/scz_celltype_paper/transcriptomic/data/figure_inputs/DE_genes_all_cells_scz.csv", row.names = FALSE)

#Forest plot dataframe
sst_pvalb_mrna <- cohort_genes_df %>%
  filter(toupper(genes) %in% c("SST", "PVALB"))

write.csv(sst_pvalb_mrna,"/scratch/nendresz/scz_celltype_paper/transcriptomic/data/figure_inputs/meta_results_cohorts_subclass_forest.csv",row.names = FALSE)




