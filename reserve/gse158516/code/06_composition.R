#!/usr/bin/env Rscript
# Compositional analysis of GSE158516 with crumblr -- the same estimator the paper's
# 7-dataset meta-analysis uses, so the betas are directly comparable.
#
# Paper model:  composition ~ diagnosis + age at death + sex + PMI
# Here sex drops out: all 26 donors are male.
#
# Neuronal and non-neuronal compartments are analysed separately, as in the paper, so
# proportions reflect within-compartment composition and are not driven by grey/white
# matter sampling differences between sections.

suppressPackageStartupMessages({
  library(crumblr); library(variancePartition); library(limma); library(data.table)
})

BASE <- "/Users/shreejoy/Github/scz_celltype_paper/reserve/gse158516"
OUTD <- file.path(BASE, "output")

# Optional prefix selects a donor subset written by 12_sensitivity_under70.py
# (e.g. "under70_"); empty means the full 26-donor cohort.
PREFIX <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else ""
cat(sprintf("counts prefix: '%s'\n", PREFIX))

meta <- fread(file.path(BASE, "data/sample_metadata.csv"))
setnames(meta, c("Deidentified ID", "Diagnosis", "Age", "Sequencing Batch", "PFC_pH"),
         c("sample_id", "diagnosis", "age", "batch", "pH"))
meta[, diagnosis := factor(diagnosis, levels = c("Control", "Schizophrenia"))]

run_crumblr <- function(level, compartment_of) {
  cnt <- fread(file.path(OUTD, sprintf("%scounts_%s.csv", PREFIX, level)))
  ids <- cnt[[1]]; cnt <- as.matrix(cnt[, -1, drop = FALSE]); rownames(cnt) <- ids
  cnt <- cnt[ids %in% meta$sample_id, , drop = FALSE]
  md  <- as.data.frame(meta[match(rownames(cnt), meta$sample_id)])

  res <- list()
  for (comp in c("neuronal", "non-neuronal")) {
    keep <- colnames(cnt)[compartment_of(colnames(cnt)) == comp]
    keep <- keep[colSums(cnt[, keep, drop = FALSE]) > 0]
    if (length(keep) < 2) next
    sub <- cnt[, keep, drop = FALSE]
    # crumblr needs every donor to contribute cells in the compartment
    ok <- rowSums(sub) > 0
    sub <- sub[ok, , drop = FALSE]; mdc <- md[ok, ]

    cobj <- crumblr(sub)
    form <- ~ diagnosis + age + PMI
    fit  <- dream(cobj, form, mdc)
    fit  <- eBayes(fit)
    tt   <- topTable(fit, coef = "diagnosisSchizophrenia", number = Inf, sort.by = "none")
    tt$cell_type   <- rownames(tt)
    tt$compartment <- comp
    tt$level       <- level
    tt$n_donors    <- nrow(sub)
    tt$mean_prop_control <- colMeans(prop.table(sub[mdc$diagnosis == "Control", , drop = FALSE], 1))
    tt$mean_prop_scz     <- colMeans(prop.table(sub[mdc$diagnosis == "Schizophrenia", , drop = FALSE], 1))
    res[[comp]] <- tt
  }
  out <- rbindlist(res, fill = TRUE)
  # FDR within compartment, matching the paper's stratified correction
  out[, FDR := p.adjust(P.Value, method = "BH"), by = compartment]
  out[order(P.Value)]
}

# Compartment assignment is read off the SEA-AD Class label exported alongside counts.
cls <- fread(file.path(OUTD, "celltype_compartment.csv"))
comp_lookup <- setNames(cls$compartment, cls$cell_type)
compartment_of <- function(x) unname(comp_lookup[x])

all_res <- rbindlist(lapply(c("subclass", "supertype"), run_crumblr, compartment_of),
                     fill = TRUE)
fwrite(all_res, file.path(OUTD, paste0(PREFIX, "composition_crumblr.csv")))

cat("\n=== SUBCLASS, neuronal compartment (top 12 by P) ===\n")
print(head(all_res[level == "subclass" & compartment == "neuronal",
                   .(cell_type, logFC, P.Value, FDR, mean_prop_control, mean_prop_scz)], 12))

cat("\n=== Sst and L6b supertypes (the paper's headline types) ===\n")
targets <- c("Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25", "L6b_1", "L6b_2", "L6b_4")
print(all_res[level == "supertype" & cell_type %in% targets,
              .(cell_type, logFC, P.Value, FDR, mean_prop_control, mean_prop_scz)][
                order(match(cell_type, targets))])

cat("\n=== Sst / Pvalb / L6b at subclass level ===\n")
print(all_res[level == "subclass" & cell_type %in% c("Sst", "Pvalb", "Vip", "Lamp5", "L6b"),
              .(cell_type, logFC, P.Value, FDR, mean_prop_control, mean_prop_scz)])

cat(sprintf("\nwrote %s\n", file.path(OUTD, paste0(PREFIX, "composition_crumblr.csv"))))
