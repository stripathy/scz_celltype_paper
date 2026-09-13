# conda activate r_env_meta_analysis
# ============================================================================
# 3a_meta_per_gene.r — per-gene meta-analysis across the 7 datasets, SUBCLASS level
# ============================================================================
# Per-gene meta-analysis for the 23-24 SEA-AD subclasses, filling the step
# between the two scripts either side of it:
#
#   2_DE.r             writes Files/DE_results_{cohort}_{celltype}.rds
#   >>> THIS FILE <<<  writes Files/meta_results_{celltype}.csv
#   3_meta_analysis.r  reads  Files/meta_results_*.csv
#
# The body is `Supertypes/3_meta_analysis.r` lines 19-53 with the cell-type
# selection inverted: that script keeps setdiff(cell_types, subclasses); this one
# keeps the subclasses themselves. Everything else is identical — same
# rma(method = "REML"), same nrow(df) > 4 threshold, same SE = |logFC/t|, same BH
# correction, same output columns and filenames. The one necessary difference is
# where the cell-type list comes from: the supertype script globs
# meta_results_*.csv, which cannot work here because this script writes them, so
# the list is taken from the pseudobulk metadata as 3_meta_analysis.r does.
#
# Added 2026-09-13. Re-running it from independently rebuilt pseudobulks recovers
# 12,490 of the 12,492 Sst genes in the committed DE_genes_all_cells_scz.csv,
# with effect sizes at r = 0.9968 and 161 of the 165 FDR < 0.05 genes; SST itself
# comes back at beta -0.459 / padj 0.048 against the committed -0.458 / 0.049.
# The small residual is TMM normalisation drift from a different per-dataset gene
# universe in the rebuild, not a difference in method — see issue 19 in
# ../../KNOWN_ISSUES.md if you are reconciling numbers to the last decimal.
# ============================================================================
setwd("P1_SCZ_DE_fresh")
library(dplyr)
library(metafor)

# Cell types come from the pseudobulk metadata, exactly as 3_meta_analysis.r
# derives them — the supertype script instead globs meta_results_*.csv, which
# cannot work here because this script is what creates those files.
meta <- read.csv("Files/Pseudobulk_metadata_subclass_Mclean.csv")
colnames(meta) <- gsub("\\.", " ", colnames(meta))
cell_types <- colnames(meta)[1:24]
cell_types <- gsub("^L2 3", "L2_3", cell_types)
cell_types <- gsub("^L5 6", "L5_6", cell_types)
cell_types <- gsub("^Micro PVM", "Micro-PVM", cell_types)

cohorts <- c("Bat", "OFC", "MtSinai", "Mclean", "MSSM", "HBCC", "Multi")

# 2_DE.r saves with safe_type <- gsub("/", "_", type); the names above already
# carry the underscore form, so they are used verbatim.
de_list <- list()
for (ct in cell_types) {
  de_list[[ct]] <- list()
  for (cohort in cohorts) {
    f <- paste0("Files/DE_results_", cohort, "_", ct, ".rds")
    if (!file.exists(f)) { message("Skipping missing: ", f); next }
    de_list[[ct]][[cohort]] <- readRDS(f)
  }
}

for (ct in cell_types) {
  if (is.null(de_list[[ct]]) || length(de_list[[ct]]) == 0) next
  message("Meta-analysis for ", ct)

  all_data <- bind_rows(lapply(names(de_list[[ct]]), function(cohort)
    de_list[[ct]][[cohort]] %>% mutate(cohort = cohort, SE = abs(logFC / t))))

  meta_results <- list()
  for (g in unique(all_data$genes)) {
    df <- all_data %>% filter(genes == g, !is.na(SE), SE != 0, !is.na(logFC))
    if (nrow(df) > 4) {
      res <- tryCatch({
        m <- rma(yi = logFC, sei = SE, data = df, method = "REML")
        data.frame(cell_type = ct, genes = g, estimate = as.numeric(m$b), se = m$se,
                   pval = m$pval, ci.lb = m$ci.lb, ci.ub = m$ci.ub,
                   k = m$k, tau2 = m$tau2, I2 = m$I2)
      }, error = function(e) NULL)
      if (!is.null(res)) meta_results[[g]] <- res
    }
  }

  final_results <- bind_rows(meta_results)
  if (nrow(final_results) > 0) final_results$padj <- p.adjust(final_results$pval, method = "fdr")
  write.csv(final_results, paste0("Files/meta_results_", ct, ".csv"), row.names = FALSE)
}
