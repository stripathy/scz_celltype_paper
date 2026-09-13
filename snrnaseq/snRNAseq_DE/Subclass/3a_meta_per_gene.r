# conda activate r_env_meta_analysis
# ============================================================================
# 3a_meta_per_gene.r — per-gene meta-analysis across the 7 datasets, SUBCLASS level
# ============================================================================
# ⚠️ RECONSTRUCTED, 2026-09-13 — this is NOT a recovered copy of the script that
# produced the published numbers. That script was never in the repo (issue 2 in
# ../../KNOWN_ISSUES.md); the original ran in Nicole's cluster working directory
# and has not been retrieved. This file was written from its recipe so that the
# chain behind Figure 2 can be executed end to end. Replace it with the original
# if that is ever recovered, and delete this notice.
#
# It fills the gap between the two committed scripts:
#
#   2_DE.r             writes Files/DE_results_{cohort}_{celltype}.rds
#   >>> THIS FILE <<<  writes Files/meta_results_{celltype}.csv
#   3_meta_analysis.r  reads  Files/meta_results_*.csv
#
# The body is `Supertypes/3_meta_analysis.r` lines 19-53 with the cell-type
# selection inverted: that script keeps setdiff(cell_types, subclasses); the
# subclass run keeps the subclasses themselves. Everything else is unchanged —
# same rma(method = "REML"), same nrow(df) > 4 threshold, same SE = |logFC/t|,
# same BH correction, same output columns and filenames.
#
# Evidence that this is the right recipe (all from the committed output,
# transcriptomic/data/figure_inputs/DE_genes_all_cells_scz.csv):
#   - its columns are exactly those built at Supertypes/3_meta_analysis.r:43-44
#     (cell_type, genes, estimate, se, pval, ci.lb, ci.ub, k, tau2, I2) + padj;
#   - k takes values {5, 6, 7} — consistent with 7 datasets and nrow(df) > 4;
#   - tau2 and I2 are non-zero, so the pooling was REML, not FE.
# Re-running this recipe from independently rebuilt pseudobulks recovers 12,490
# of the 12,492 committed Sst genes with effect sizes r = 0.9968 and 161 of 165
# FDR < 0.05 genes; it is not bit-exact because the rebuild used a different
# per-dataset gene universe. See issue 19.
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
