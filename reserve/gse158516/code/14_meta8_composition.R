#!/usr/bin/env Rscript
# Add GSE158516 to the snRNA-seq compositional meta-analysis: 7 datasets -> 8.
#
# THE `_fe` COLUMNS ARE EXACT, NOT AN APPROXIMATION. This is a two-stage update -- the
# published 7-dataset pooled summary combined with GSE158516 as an eighth study -- because
# the per-cohort estimates for the original seven are not in this repo. That would
# normally be an approximation. It is not, for two reasons established in
# 16_validate_pooling_method.R:
#
#   1. The published compositional pool is FIXED-EFFECT (FE reproduces it to ~1e-16;
#      REML/DL/ML/HE/SJ/EB/PM are all off by ~2e-2), despite the Methods describing
#      random-effects REML. `final_results_crumblr_7_cohorts.csv` carries no tau2/I2/k
#      columns, consistent with FE, while the DE table does carry them.
#   2. Inverse-variance weighting is associative, so FE-pooling (7-dataset summary +
#      GSE158516) is algebraically identical to FE-pooling all eight cohorts at once.
#      Verified directly against per-cohort data for the 16 SST supertypes: max |dbeta|
#      = 7.5e-16.
#
# So `beta_fe`/`fdr_fe` are the exact 8-dataset meta-analysis under the paper's own
# pooling model. `beta_dl`/`fdr_dl` answer a different question -- what a random-effects
# pool would say -- and are reported alongside the per-cell-type Q test so that any
# heterogeneity the FE model ignores stays visible.
#
# Run with the counts prefix to use as the eighth study:
#   Rscript 14_meta8_composition.R           # all 26 donors (age covaried)
#   Rscript 14_meta8_composition.R under70_  # 21 donors, matching the paper's <70 rule

suppressPackageStartupMessages({library(metafor); library(data.table)})

BASE   <- "/Users/shreejoy/Github/scz_celltype_paper/reserve/gse158516"
OUTD   <- file.path(BASE, "output")
SHARED <- "/Users/shreejoy/Github/scz_celltype_paper/shared/snrnaseq_de/nicole_scz_snrnaseq_betas"

PREFIX <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else ""
LABEL  <- if (PREFIX == "") "all 26 donors" else "21 donors (age <= 70)"
cat(sprintf("Eighth study: GSE158516, %s\n\n", LABEL))

HEADLINE <- c("Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25", "L6b_1", "L6b_2", "L6b_4")

# ---- the existing 7-dataset pooled result (neuronal + non-neuronal) ----------
meta7 <- rbind(
  data.table(fread(file.path(SHARED, "final_results_crumblr_7_cohorts.csv")),
             compartment = "neuronal"),
  data.table(fread(file.path(SHARED, "final_results_crumblr_7_nonN_cohorts.csv")),
             compartment = "non-neuronal"))
setnames(meta7, c("CellType", "estimate", "se", "pval", "padj"),
         c("cell_type", "beta7", "se7", "p7", "fdr7"))
meta7 <- meta7[, .(cell_type, compartment, beta7, se7, p7, fdr7)]

# ---- GSE158516 as the eighth study -------------------------------------------
g <- fread(file.path(OUTD, paste0(PREFIX, "composition_crumblr.csv")))[level == "supertype"]
# crumblr/dream's topTable reports t but not SE; SE = logFC / t
g[, se8 := logFC / t]
g <- g[, .(cell_type, beta8 = logFC, se8, p8 = P.Value)]

m <- merge(meta7, g, by = "cell_type")
cat(sprintf("cell types in both: %d of %d (meta-7) and %d (GSE158516)\n",
            nrow(m), nrow(meta7), nrow(g)))
cat(sprintf("median SE: meta-7 %.4f, GSE158516 %.4f  ->  GSE158516 carries a median "
            , median(m$se7), median(m$se8)))
cat(sprintf("%.1f%% of the inverse-variance weight\n\n",
            100 * median(m$se7^2 / (m$se7^2 + m$se8^2))))

# ---- combine, per cell type ---------------------------------------------------
res <- rbindlist(lapply(seq_len(nrow(m)), function(i) {
  yi  <- c(m$beta7[i], m$beta8[i]); sei <- c(m$se7[i], m$se8[i])
  fe  <- rma(yi = yi, sei = sei, method = "FE")
  dl  <- rma(yi = yi, sei = sei, method = "DL")
  data.table(cell_type = m$cell_type[i], compartment = m$compartment[i],
             beta7 = m$beta7[i], se7 = m$se7[i], p7 = m$p7[i], fdr7 = m$fdr7[i],
             beta8_new = m$beta8[i], se8_new = m$se8[i], p8_new = m$p8[i],
             weight8 = m$se7[i]^2 / (m$se7[i]^2 + m$se8[i]^2),
             beta_fe = fe$beta[1], se_fe = fe$se, p_fe = fe$pval,
             beta_dl = dl$beta[1], se_dl = dl$se, p_dl = dl$pval,
             Q = fe$QE, Q_p = fe$QEp, I2 = dl$I2)
}))

# FDR within compartment, matching the paper's stratified correction
res[, fdr_fe := p.adjust(p_fe, "BH"), by = compartment]
res[, fdr_dl := p.adjust(p_dl, "BH"), by = compartment]
res[, concordant := sign(beta7) == sign(beta8_new)]
res[, shift := beta_fe - beta7]

# ---- what actually changes ----------------------------------------------------
cat("=== Do the two studies agree in direction? ===\n")
cat(sprintf("  all cell types:   %d/%d (%.0f%%) concordant\n",
            sum(res$concordant), nrow(res), 100 * mean(res$concordant)))
sig7 <- res[fdr7 < 0.10]
cat(sprintf("  meta-7 FDR<0.10:  %d/%d (%.0f%%) concordant\n",
            sum(sig7$concordant), nrow(sig7), 100 * mean(sig7$concordant)))
cat(sprintf("  heterogeneity:    %d/%d cell types with Q p<0.05\n\n",
            sum(res$Q_p < 0.05), nrow(res)))

cat("=== The paper's headline supertypes ===\n")
h <- res[cell_type %in% HEADLINE][order(match(cell_type, HEADLINE))]
print(h[, .(cell_type, beta7 = round(beta7, 3), fdr7 = signif(fdr7, 2),
            beta8_new = round(beta8_new, 3), p8_new = signif(p8_new, 2),
            beta_fe = round(beta_fe, 3), fdr_fe = signif(fdr_fe, 2),
            fdr_dl = signif(fdr_dl, 2), Q_p = signif(Q_p, 2))])

cat("\n=== Significance status changes (FDR 0.05 / 0.10 / 0.20) ===\n")
for (thr in c(0.05, 0.10, 0.20)) {
  gained <- res[fdr7 >= thr & fdr_fe < thr, cell_type]
  lost   <- res[fdr7 < thr & fdr_fe >= thr, cell_type]
  cat(sprintf("  FDR<%.2f: %d -> %d types", thr, sum(res$fdr7 < thr),
              sum(res$fdr_fe < thr)))
  if (length(gained)) cat(sprintf("  | gained: %s", paste(gained, collapse = ", ")))
  if (length(lost))   cat(sprintf("  | lost: %s", paste(lost, collapse = ", ")))
  cat("\n")
}

fwrite(res, file.path(OUTD, paste0(PREFIX, "meta8_composition.csv")))

cat("\n=== Largest shifts in the pooled estimate ===\n")
print(head(res[order(-abs(shift)), .(cell_type, beta7 = round(beta7, 3),
                                     beta8_new = round(beta8_new, 3),
                                     beta_fe = round(beta_fe, 3),
                                     shift = round(shift, 4),
                                     weight8 = round(weight8, 3))], 10))

cat(sprintf("\nmedian |shift| in beta: %.4f (%.1f%% of the median |beta7| of %.3f)\n",
            median(abs(res$shift)),
            100 * median(abs(res$shift)) / median(abs(res$beta7)),
            median(abs(res$beta7))))
cat(sprintf("wrote %s\n", file.path(OUTD, paste0(PREFIX, "meta8_composition.csv"))))
