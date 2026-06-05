#!/usr/bin/env Rscript
# ============================================================================
# Cell-type-specific differential expression (SCZ vs Control) for the Xenium
# data, native-R edgeR. Pseudobulk QL F-test, the DE analogue of the crumblr
# compositional analysis. Replaces the archived edgepython script.
#
# Same cell/donor inclusion as crumblr (set in build_de_input.py): cortical
# cells passing corr_qc_pass, typed by corr_subclass / corr_supertype, 24 donors.
#
# Per cell type: DGEList -> filterByExpr -> calcNormFactors(TMM) ->
#   estimateDisp(robust) -> glmQLFit(robust) -> glmQLFTest(coef = dxSCZ),
#   design ~ diagnosis + sex + age(centered).
#
# Usage:  Rscript run_de.R [subclass|supertype]    (default: subclass)
# Input:  output/de/pseudobulk_{level}.csv  + pseudobulk_{level}_samples.csv
# Output: output/de/de_results_{level}_v2.csv   (review vs committed before swap)
# ============================================================================
suppressPackageStartupMessages({
  library(edgeR); library(dplyr); library(tidyr); library(readr)
})

args  <- commandArgs(trailingOnly = TRUE)
level <- if (length(args) >= 1) args[1] else "subclass"
IN    <- "output/de"
MIN_PER_GROUP <- 3                     # min donors per group for a cell type to be tested

cat(sprintf("DE level: %s\n", level))
long <- read_csv(sprintf("%s/pseudobulk_%s.csv", IN, level), show_col_types = FALSE)
samp <- read_csv(sprintf("%s/pseudobulk_%s_samples.csv", IN, level), show_col_types = FALSE)

run_ct <- function(ct) {
  sm <- samp[samp$celltype == ct, ]
  n_scz <- sum(sm$diagnosis == "SCZ"); n_ctrl <- sum(sm$diagnosis == "Control")
  if (n_scz < MIN_PER_GROUP || n_ctrl < MIN_PER_GROUP) return(NULL)

  # genes x donors count matrix
  w <- long[long$celltype == ct, c("gene", "donor", "count")] %>%
       pivot_wider(names_from = donor, values_from = count)
  genes <- w$gene
  mat <- as.matrix(w[, -1]); rownames(mat) <- genes
  sm <- sm[match(colnames(mat), sm$donor), ]

  # design ~ diagnosis + sex + age_c (drop a covariate if constant in this subset)
  dx  <- factor(sm$diagnosis, levels = c("Control", "SCZ"))
  terms <- "dx"
  if (length(unique(sm$sex)) > 1) { sex <- factor(sm$sex); terms <- c(terms, "sex") }
  age_c <- sm$age - mean(sm$age)
  if (sd(age_c) > 0) terms <- c(terms, "age_c")
  design <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))))
  if (!"dxSCZ" %in% colnames(design)) return(NULL)

  d <- DGEList(counts = mat)
  keep <- filterByExpr(d, design)
  if (sum(keep) < 5) return(NULL)
  d <- d[keep, , keep.lib.sizes = FALSE]
  d <- calcNormFactors(d, method = "TMM")
  d <- estimateDisp(d, design, robust = TRUE)
  fit <- glmQLFit(d, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "dxSCZ")
  tt  <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  tt$gene <- rownames(tt)
  tt$n_scz <- n_scz; tt$n_ctrl <- n_ctrl; tt$n_samples <- nrow(sm)
  tt$n_genes_tested <- nrow(tt)
  tt$celltype <- ct
  tt$cell_class <- sm$cell_class[1]
  tt
}

cts <- unique(long$celltype)
res <- list()
for (ct in cts) {
  r <- tryCatch(run_ct(ct), error = function(e) { cat(sprintf("  ERROR %s: %s\n", ct, conditionMessage(e))); NULL })
  if (!is.null(r)) { res[[ct]] <- r; cat(sprintf("  %-16s %3d genes, %d donors (%d SCZ / %d Ctrl)\n",
                                                  ct, nrow(r), r$n_samples[1], r$n_scz[1], r$n_ctrl[1])) }
  else cat(sprintf("  %-16s skipped\n", ct))
}

out <- bind_rows(res) %>%
  select(gene, logFC, logCPM, PValue, F = `F`, FDR,
         n_scz, n_ctrl, n_samples, n_genes_tested, celltype, cell_class)
outpath <- sprintf("%s/de_results_%s_v2.csv", IN, level)
write_csv(out, outpath)
cat(sprintf("\nSaved %s  (%d rows, %d cell types)\n", outpath, nrow(out), length(res)))
