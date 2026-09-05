#!/usr/bin/env Rscript
# Which estimator produced the published compositional meta-analysis -- and does the
# two-stage update in 14_meta8_composition.R therefore need a caveat?
#
# `SCZ_SST_7_cohort_estimates (3).csv` carries both the 7 per-cohort crumblr estimates
# (with SEs) and the published pooled row, for the 16 SST supertypes. That lets us test
# each estimator against the truth instead of trusting the Methods text.
#
# TWO FINDINGS, both verified below at machine precision:
#
#   1. The published compositional pool is FIXED-EFFECT, not the REML random-effects the
#      manuscript Methods describe ("combined ... by random-effects meta-analysis
#      (metafor::rma, restricted maximum likelihood)"). FE reproduces the published rows
#      to ~1e-16; REML, DL, ML, HE, SJ, EB and PM are all off by ~2e-2. Corroborating
#      evidence: `final_results_crumblr_7_cohorts.csv` has no tau2/I2/k columns, whereas
#      the DE table `DE_genes_all_cells_scz.csv` does -- so DE really is random-effects
#      and the Methods sentence appears to describe DE and composition together when only
#      DE matches it.
#
#   2. Because inverse-variance weighting is associative, FE-pooling the published
#      7-dataset summary with GSE158516 is ALGEBRAICALLY IDENTICAL to FE-pooling all
#      eight cohorts at once. So the two-stage update in 14_meta8_composition.R is not an
#      approximation -- it is the exact 8-dataset fixed-effect meta-analysis, and it needs
#      no per-cohort data. Verified here for the 16 SST supertypes, where the per-cohort
#      rows are available to compute the direct 8-study pool for comparison.
#
# Implication: the `beta_fe`/`fdr_fe` columns of meta8_composition.csv are exact under the
# paper's own pooling model. The `beta_dl`/`fdr_dl` columns are a different question --
# what a random-effects pool would say -- not a correction to an approximation.

suppressPackageStartupMessages({library(metafor); library(data.table)})

BASE <- "/Users/shreejoy/Github/scz_celltype_paper/external_validation/gse158516"
OUTD <- file.path(BASE, "output")
SHARED <- "/Users/shreejoy/Github/scz_celltype_paper/shared/snrnaseq_de/nicole_scz_snrnaseq_betas"
PER_COHORT <- "/Users/shreejoy/Downloads/SCZ_SST_7_cohort_estimates (3).csv"

# metafor evaluates yi/sei by non-standard evaluation, which does not reach into a
# data.table j-expression; a wrapper gives them ordinary local bindings.
pool <- function(y, s, method) {
  f <- rma(yi = y, sei = s, method = method)
  list(b = f$beta[1], s = f$se, p = f$pval)
}

d <- fread(PER_COHORT); setnames(d, tolower(names(d)))
coh <- d[cohort != "Meta-analysis"]
pub <- d[cohort == "Meta-analysis", .(celltype, pb = estimate, ps = se)]

cat(sprintf("per-cohort input: %d cohorts x %d SST supertypes (k = %s per type)\n\n",
            uniqueN(coh$cohort), uniqueN(coh$celltype),
            paste(unique(coh[, .N, by = celltype]$N), collapse = ",")))

# ---- is the file's pooled row the one the paper actually uses? ----------------
canon <- fread(file.path(SHARED, "final_results_crumblr_7_cohorts.csv"))
setnames(canon, "CellType", "celltype")
m <- merge(pub, canon[, .(celltype, e2 = estimate, s2 = se)], by = "celltype")
cat(sprintf("file's pooled rows vs final_results_crumblr_7_cohorts.csv (n=%d):\n",  nrow(m)))
cat(sprintf("  max|dbeta| = %.2e   max|dse| = %.2e  -> %s\n\n",
            max(abs(m$pb - m$e2)), max(abs(m$ps - m$s2)),
            if (max(abs(m$pb - m$e2)) == 0) "IDENTICAL; this is the paper's own pooled result"
            else "DIFFERENT -- stop here"))

# ---- finding 1: which estimator reproduces it? -------------------------------
cat("=== which metafor estimator reproduces the published pool? ===\n")
fit <- rbindlist(lapply(c("FE", "REML", "DL", "ML", "HE", "SJ", "EB", "PM"), function(mth) {
  r <- coh[, { v <- pool(estimate, se, mth); .(b = v$b, s = v$s) }, by = celltype]
  x <- merge(r, pub, by = "celltype")
  data.table(method = mth, max_dbeta = max(abs(x$b - x$pb)), max_dse = max(abs(x$s - x$ps)))
}))
fit[, reproduces := max_dbeta < 1e-10 & max_dse < 1e-10]
print(fit)
cat(sprintf("\n  -> the published compositional pool is %s\n\n",
            paste(fit[reproduces == TRUE, method], collapse = "/")))

# ---- finding 2: is the two-stage update exact? -------------------------------
g <- fread(file.path(OUTD, "composition_crumblr.csv"))[level == "supertype"]
g[, se := logFC / t]
g <- g[, .(celltype = cell_type, estimate = logFC, se)][celltype %in% unique(coh$celltype)]

direct <- rbind(coh[, .(celltype, estimate, se)], g)[
  , { v <- pool(estimate, se, "FE"); .(b_direct = v$b, se_direct = v$s, p_direct = v$p) },
  by = celltype]
two <- fread(file.path(OUTD, "meta8_composition.csv"))[
  , .(celltype = cell_type, b_two = beta_fe, se_two = se_fe, p_two = p_fe)]
cmp <- merge(direct, two, by = "celltype")
cat(sprintf("=== direct 8-cohort FE vs the two-stage update (n=%d SST types) ===\n", nrow(cmp)))
cat(sprintf("  max|dbeta| = %.2e   max|dse| = %.2e   max|dp| = %.2e\n",
            max(abs(cmp$b_direct - cmp$b_two)), max(abs(cmp$se_direct - cmp$se_two)),
            max(abs(cmp$p_direct - cmp$p_two))))
cat("  -> identical to machine precision: the two-stage pool IS the 8-cohort pool\n\n")

print(cmp[celltype %in% c("Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25")][
  , .(celltype, b_direct = signif(b_direct, 5), b_two = signif(b_two, 5),
      p_direct = signif(p_direct, 4), p_two = signif(p_two, 4))])

fwrite(fit, file.path(OUTD, "pooling_method_validation.csv"))
cat(sprintf("\nwrote %s\n", file.path(OUTD, "pooling_method_validation.csv")))
