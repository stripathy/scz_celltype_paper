#!/usr/bin/env Rscript
# S11: formal diagnosis x stratum interaction, the test behind Fig 5's
# "graded by depletion status" claim.
# Per cohort (donor x stratum pseudobulks from script 28, joint TMM):
#   voom -> duplicateCorrelation(block = donor) -> voom -> lmFit(block, cor)
#   model A: ~ dx * stratum + covars   (stratum ref = non_depleted)
#     -> coefficient dxSCZ:stratumdepleted = depleted-vs-non-depleted
#        difference in the diagnosis effect (also the intermediate one)
#   model B: ~ dx * score + covars     (score: non=0, int=1, dep=2)
#     -> dxSCZ:score = single-df linear trend of the diagnosis effect across
#        strata ("graded" as one coefficient)
# Meta: metafor REML, k >= 5; BH within coefficient. Then fgsea on the
# interaction z ranks. Donor-level confounds cancel in these coefficients.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(tibble); library(arrow); library(edgeR); library(limma)
  library(metafor); library(parallel); library(fgsea); library(msigdbr)
})
set.seed(42)
EXP  <- "transcriptomic/data/stratum_pseudobulks_export"
POUT <- "transcriptomic/results/sst_strata_gsea/pseudobulk"
DOUT <- file.path(POUT, "donor_stratum")
N_CORES <- max(1, detectCores() - 2)
t0 <- Sys.time()
say <- function(...) cat(sprintf("[%6.1f min] ", as.numeric(difftime(Sys.time(), t0, units = "mins"))), sprintf(...), "\n", sep = "")

cohorts <- str_remove(basename(Sys.glob(file.path(DOUT, "*_donor_stratum_meta.csv"))),
                      "_donor_stratum_meta\\.csv$")
say("cohorts: %s", paste(cohorts, collapse = ", "))

covars <- map_dfr(cohorts, function(co)
  read_csv(file.path(EXP, paste0(co, "_groups.csv")),
           col_types = cols(donor = col_character(), .default = col_guess())) |>
    distinct(donor, sex, age, pmi) |>
    mutate(cohort = co,
           age = suppressWarnings(as.numeric(as.character(age))),
           pmi = suppressWarnings(as.numeric(as.character(pmi)))))

run_cohort <- function(co) {
  pb <- as.data.frame(read_parquet(file.path(DOUT, paste0(co, "_donor_stratum_counts.parquet"))))
  rn <- pb$gene; pb$gene <- NULL
  M <- vapply(pb, as.numeric, numeric(nrow(pb))); rownames(M) <- rn
  m <- read_csv(file.path(DOUT, paste0(co, "_donor_stratum_meta.csv")),
                col_types = cols(donor = col_character(), .default = col_guess())) |>
    left_join(covars |> filter(cohort == co) |> select(-cohort), by = "donor") |>
    mutate(key = paste(donor, stratum, sep = "|"),
           dx = factor(diagnosis, c("Control", "SCZ")),
           stratum = factor(stratum, c("non_depleted", "intermediate", "depleted")),
           score = as.integer(stratum) - 1L,
           sex = factor(sex),
           age_s = scale(ifelse(is.na(age), median(age, na.rm = TRUE), age))[, 1])
  use_pmi <- mean(is.na(m$pmi)) < 0.5
  if (use_pmi) m$pmi_s <- scale(ifelse(is.na(m$pmi), median(m$pmi, na.rm = TRUE), m$pmi))[, 1]
  M <- M[, m$key]
  keep <- rowSums(M >= 1) >= 0.8 * nrow(m)
  y <- calcNormFactors(DGEList(M[keep, ]))
  say("%s: %d samples (%d donors), %d genes, pmi=%s",
      co, nrow(m), n_distinct(m$donor), sum(keep), use_pmi)

  fit_int <- function(f, coefs) {
    des <- model.matrix(f, data = m)
    v <- voom(y, des)
    dc <- duplicateCorrelation(v, des, block = m$donor)
    v <- voom(y, des, block = m$donor, correlation = dc$consensus.correlation)
    fit <- eBayes(lmFit(v, des, block = m$donor,
                        correlation = dc$consensus.correlation))
    map_dfr(coefs, function(cf) {
      tt <- topTable(fit, coef = cf, number = Inf, sort.by = "none")
      tibble(cohort = co, coef = cf, gene = rownames(tt),
             logFC = tt$logFC, SE = tt$logFC / tt$t, P.Value = tt$P.Value)
    })
  }
  fA <- if (use_pmi) ~ dx * stratum + age_s + sex + pmi_s else ~ dx * stratum + age_s + sex
  fB <- if (use_pmi) ~ dx * score + age_s + sex + pmi_s else ~ dx * score + age_s + sex
  bind_rows(fit_int(fA, c("dxSCZ:stratumdepleted", "dxSCZ:stratumintermediate")),
            fit_int(fB, "dxSCZ:score"))
}
percohort <- map_dfr(cohorts, run_cohort)
write_csv(percohort, file.path(POUT, "interaction_percohort.csv"))

say("meta-analysis of interaction coefficients ...")
meta_one <- function(d) {
  fit <- tryCatch(rma(yi = d$logFC, sei = d$SE, method = "REML",
                      control = list(maxiter = 1000)),
                  error = function(e) tryCatch(rma(yi = d$logFC, sei = d$SE, method = "DL"),
                                               error = function(e2) NULL))
  if (is.null(fit)) return(NULL)
  tibble(estimate = fit$beta[1], se = fit$se, zval = fit$zval, pval = fit$pval, k = fit$k)
}
meta <- map_dfr(unique(percohort$coef), function(cf) {
  dd <- percohort |> filter(coef == cf)
  genes <- dd |> count(gene) |> filter(n >= 5) |> pull(gene)
  say("  %s: %d genes with k >= 5", cf, length(genes))
  res <- mclapply(genes, function(gn) meta_one(dd[dd$gene == gn, ]), mc.cores = N_CORES)
  bind_rows(setNames(res, genes), .id = "gene") |> mutate(coef = cf)
}) |> group_by(coef) |> mutate(padj = p.adjust(pval, "BH")) |> ungroup()
write_csv(meta, file.path(POUT, "interaction_meta.csv"))

for (cf in unique(meta$coef)) {
  s <- meta |> filter(coef == cf)
  say("%s: %d genes | FDR<0.05: %d | FDR<0.10: %d", cf, nrow(s),
      sum(s$padj < 0.05), sum(s$padj < 0.10))
}
say("key genes (interaction z):")
print(meta |> filter(gene %in% c("SST", "VGF", "RPL36", "NDUFS8", "RPL10", "SCO1")) |>
        select(gene, coef, z = zval, padj) |>
        mutate(z = round(z, 2), padj = signif(padj, 2)) |>
        pivot_wider(names_from = coef, values_from = c(z, padj)), width = 220)
say("top 15 trend-interaction genes (dx:score, by pval):")
print(meta |> filter(coef == "dxSCZ:score") |> arrange(pval) |>
        transmute(gene, z = round(zval, 2), padj = signif(padj, 2)) |> head(15), n = 15)

# ---- interaction GSEA -----------------------------------------------------------
say("interaction GSEA ...")
collect <- function(cat, subcat = NULL) {
  args <- list(species = "Homo sapiens", collection = cat)
  if (!is.null(subcat)) args$subcollection <- subcat
  do.call(msigdbr, args) |> transmute(gs_name, gene = gene_symbol)
}
gs <- bind_rows(collect("H"), collect("C5", "GO:BP"), collect("C5", "GO:CC"),
                collect("C5", "GO:MF"), collect("C2", "CP:REACTOME")) |> distinct()
gsl <- split(gs$gene, gs$gs_name)
gsea <- map_dfr(c("dxSCZ:stratumdepleted", "dxSCZ:score"), function(cf) {
  d <- meta |> filter(coef == cf, !is.na(zval)) |>
    arrange(desc(abs(zval))) |> distinct(gene, .keep_all = TRUE)
  ranks <- sort(setNames(d$zval, d$gene), decreasing = TRUE)
  fgsea(gsl, ranks, minSize = 10, maxSize = 500, nPermSimple = 10000) |>
    as_tibble() |> mutate(coef = cf)
})
write_csv(gsea |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = "|")),
          file.path(POUT, "interaction_gsea.csv"))
KEY <- c("REACTOME_TRANSLATION", "GOCC_RIBOSOMAL_SUBUNIT", "GOBP_CYTOPLASMIC_TRANSLATION",
         "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
         "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "GOBP_OXIDATIVE_PHOSPHORYLATION",
         "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
         "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING", "GOCC_SYNAPTIC_MEMBRANE",
         "GOBP_NEUROTRANSMITTER_SECRETION")
fmt <- function(NES, padj) sprintf("%.1f%s", NES, ifelse(padj < 0.01, "**",
                 ifelse(padj < 0.05, "*", ifelse(padj < 0.1, "+", ""))))
say("=== interaction GSEA, key blocks (negative NES = diagnosis effect more")
say("=== negative in depleted than non-depleted) ===")
print(gsea |> filter(pathway %in% KEY) |> mutate(cell = fmt(NES, padj)) |>
        select(pathway, coef, cell) |>
        pivot_wider(names_from = coef, values_from = cell), n = 12, width = 160)
say("=== top 10 interaction gene sets by padj (dep-vs-non coefficient) ===")
print(gsea |> filter(coef == "dxSCZ:stratumdepleted") |> arrange(padj) |>
        transmute(pathway = str_trunc(pathway, 55), NES = round(NES, 2),
                  padj = signif(padj, 2)) |> head(10), n = 10)
