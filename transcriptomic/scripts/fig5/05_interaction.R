#!/usr/bin/env Rscript
# Step 5 | Diagnosis x stratum interaction -- the test behind "graded by
# depletion status".
#
# Every donor contributes one pseudobulk per stratum (step 4), so the contrast is
# taken WITHIN donor and donor-level confounds (medication, age, sex, PMI, batch)
# cancel by construction. Repeated measures on the same donor are handled with
# duplicateCorrelation(block = donor).
#   model A: ~ dx * stratum + covars      (stratum reference = non_depleted)
#            dxSCZ:stratumdepleted     = depleted vs non-depleted difference
#            dxSCZ:stratumintermediate = intermediate vs non-depleted difference
#   model B: ~ dx * score + covars        (score: non = 0, int = 1, dep = 2)
#            dxSCZ:score = single-df linear trend across the three strata
# Meta: metafor REML, k >= 5; BH within coefficient. Then GSEA on interaction z.
#
# Runtime ~35 min (voom is fit twice per model for duplicateCorrelation).
# Enrichment on these ranks is a separate step (05b) so it can be re-run cheaply.
source("transcriptomic/scripts/fig5/_common.R")
suppressPackageStartupMessages({
  library(edgeR); library(limma); library(parallel)
})
set.seed(42)
N_CORES <- max(1, detectCores() - 2)

cohorts <- donor_cohorts()
say("cohorts: %s", paste(cohorts, collapse = ", "))
covars <- load_covars()

run_cohort <- function(co) {
  d <- load_donor_stratum(co)
  m <- d$meta |>
    left_join(covars |> filter(cohort == co) |> select(-cohort), by = "donor") |>
    mutate(key = paste(donor, stratum, sep = "|"),
           # non_depleted first: it is the reference level for model A
           stratum = factor(stratum, c("non_depleted", "intermediate", "depleted")),
           score = as.integer(stratum) - 1L) |>
    prep_model_frame()
  use_pmi <- attr(m, "use_pmi")
  M <- d$counts[, m$key]
  y <- calcNormFactors(DGEList(M[rowSums(M >= 1) >= 0.8 * nrow(m), ]))
  say("%s: %d samples (%d donors), %d genes, pmi=%s",
      co, nrow(m), n_distinct(m$donor), nrow(y), use_pmi)

  fit_int <- function(f, coefs) {
    des <- model.matrix(f, data = m)
    stopifnot(nrow(des) == ncol(y))
    v <- voom(y, des)
    dc <- duplicateCorrelation(v, des, block = m$donor)
    v <- voom(y, des, block = m$donor, correlation = dc$consensus.correlation)
    fit <- eBayes(lmFit(v, des, block = m$donor, correlation = dc$consensus.correlation))
    map_dfr(coefs, function(cf) {
      tt <- topTable(fit, coef = cf, number = Inf, sort.by = "none")
      tibble(cohort = co, coef = cf, gene = rownames(tt), logFC = tt$logFC,
             SE = tt$logFC / tt$t, P.Value = tt$P.Value)
    })
  }
  fA <- if (use_pmi) ~ dx * stratum + age_s + sex + pmi_s else ~ dx * stratum + age_s + sex
  fB <- if (use_pmi) ~ dx * score   + age_s + sex + pmi_s else ~ dx * score   + age_s + sex
  bind_rows(fit_int(fA, c("dxSCZ:stratumdepleted", "dxSCZ:stratumintermediate")),
            fit_int(fB, "dxSCZ:score"))
}
percohort <- map_dfr(cohorts, run_cohort)
write_csv(percohort, file.path(P$pb, "interaction_percohort.csv"))

say("meta-analysis of interaction coefficients ...")
meta <- map_dfr(unique(percohort$coef), function(cf) {
  dd <- percohort |> filter(coef == cf)
  genes <- dd |> count(gene) |> filter(n >= 5) |> pull(gene)
  say("  %s: %d genes with k >= 5", cf, length(genes))
  res <- mclapply(genes, function(gn) {
    d <- dd[dd$gene == gn, ]; fit_meta(d$logFC, d$SE)
  }, mc.cores = N_CORES)
  bind_rows(setNames(res, genes), .id = "gene") |> mutate(coef = cf)
}) |> group_by(coef) |> mutate(padj = p.adjust(pval, "BH")) |> ungroup()
write_csv(meta, file.path(P$pb, "interaction_meta.csv"))

for (cf in unique(meta$coef)) {
  s <- meta |> filter(coef == cf)
  say("%s: %d genes | FDR<0.05 %d | FDR<0.10 %d", cf, nrow(s),
      sum(s$padj < 0.05), sum(s$padj < 0.10))
}
say("FDR<0.05 interaction genes (depleted vs non-depleted):")
print(meta |> filter(coef == "dxSCZ:stratumdepleted", padj < 0.05) |> arrange(padj) |>
        transmute(gene, z = round(zval, 2), padj = signif(padj, 2)))
