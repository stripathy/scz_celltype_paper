#conda activate r_env_meta_analysis
setwd("scz_celltype_paper/snrnaseq/snRNAseq_DE")
library(dplyr)
library(ggplot2)
library(tibble)
library(ggrepel)
library(metafor)
#List cohorts
cohorts <- c("Batiuk", "Frohlich", "MSSM1", "McLean", "MSSM2", "HBCC", "Multiome")
#Just to get cell types
meta <- read.csv("Files/Pseudobulk_metadata_subclass_McLean.csv")
colnames(meta)  <- gsub("\\.", " ", colnames(meta))
# every column that is not donor metadata (see 2_DE.r)
cell_types <- setdiff(colnames(meta), c("Donor", "Age", "Sex", "Diagnosis", "PMI"))
cell_types <- gsub("^L2 3", "L2_3", cell_types)
cell_types <- gsub("^L5 6", "L5_6", cell_types)
cell_types <- gsub("^Micro PVM", "Micro-PVM", cell_types)
print(cell_types)

# Get one list of all cohorts per cell types
de_list <- list()
for (ct in cell_types) {
  de_list[[ct]] <- list()
  for (cohort in cohorts) {
    file <- paste0(
      "Files/DE_results_",
      cohort,
      "_",
      (ct),
      ".rds"
    )
    de_list[[ct]][[cohort]] <- readRDS(file)
  }
}


meta_all <- list()

for (ct in cell_types) {

  
  message("Meta-analysis for ", ct)
  
  all_data <- bind_rows(
    lapply(names(de_list[[ct]]), function(cohort) {
      de_list[[ct]][[cohort]] %>%
        mutate(
          cohort = cohort,
          SE = abs(logFC / t)
        )
    })
  )
  
  all_genes <- unique(all_data$genes)
  meta_results <- list()
  
  for (g in all_genes) {
    
    df <- all_data %>%
      filter(genes == g, !is.na(SE), SE != 0, !is.na(logFC))
    
    if (nrow(df) > 4) {
      
      res <- tryCatch({
        model <- rma(yi = logFC, sei = SE, data = df, method = "REML")
        
        data.frame(
          cell_type = ct,
          genes = g,
          estimate = as.numeric(model$b),
          se = model$se,
          pval = model$pval,
          ci.lb = model$ci.lb,
          ci.ub = model$ci.ub,
          k = model$k,
          tau2 = model$tau2,
          I2 = model$I2
        )
      }, error = function(e) NULL)
      
      if (!is.null(res)) meta_results[[g]] <- res
    }
  }
  
  final_results <- bind_rows(meta_results)
  final_results$padj <- p.adjust(final_results$pval, method = "fdr")
  
  meta_all[[ct]] <- final_results
  
write.csv(final_results, paste0("Files/meta_results_", ct, ".csv"))}
