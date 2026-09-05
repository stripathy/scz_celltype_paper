#!/usr/bin/env Rscript
# Early-epoch sensitivity for the A9 CPS analysis: donors with CPS <= 0.6
# (our operationalization of Gabitto's early epoch — their main text places the
# early/late transition at CPS 0.4-0.6; their exact two-epoch effect sizes come
# from a model specified in their Supplementary Note). Severely affected donors
# excluded as before (all are late-epoch anyway). Everything else mirrors 21b.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(tibble)
  library(arrow); library(edgeR); library(limma); library(fgsea); library(msigdbr)
})
set.seed(42)
OUT  <- "transcriptomic/results/sst_strata_gsea"
AOUT <- file.path(OUT, "seaad_a9")
MIN_CELLS <- 10
CPS_MAX <- 0.6
STRATA <- c("depleted", "intermediate", "non_depleted")
KEY <- c("GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING", "GOBP_NEUROTRANSMITTER_SECRETION",
         "REACTOME_TRANSLATION", "GOCC_RIBOSOMAL_SUBUNIT", "GOBP_CYTOPLASMIC_TRANSLATION",
         "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
         "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "GOBP_OXIDATIVE_PHOSPHORYLATION",
         "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
         "GOMF_UBIQUITIN_PROTEIN_TRANSFERASE_ACTIVITY",
         "GOMF_PROTEIN_TYROSINE_KINASE_ACTIVITY", "GOBP_PROTEIN_UBIQUITINATION")
t0 <- Sys.time()
say <- function(...) cat(sprintf("[%5.1f min] ", as.numeric(difftime(Sys.time(), t0, units="mins"))), sprintf(...), "\n", sep = "")

pb <- as.data.frame(read_parquet(file.path(AOUT, "a9_stratum_pseudobulk.parquet")))
idx <- grep("__index_level_0__", colnames(pb), value = TRUE)
rownames(pb) <- pb[[idx]]; pb[[idx]] <- NULL; pb <- as.matrix(pb)
meta <- read_csv(file.path(AOUT, "a9_stratum_pseudobulk_meta.csv"), show_col_types = FALSE) |>
  mutate(severe = as.character(severe) %in% c("Y", "True", "TRUE", "1"),
         age = as.numeric(str_replace(as.character(age), "\\+", "")),
         pmi = suppressWarnings(as.numeric(pmi))) |>
  filter(!severe, CPS <= CPS_MAX)
say("early-epoch donors (CPS <= %.1f): %d", CPS_MAX, n_distinct(meta$donor))

de_cps <- function(keys) {
  m <- meta |> filter(key %in% keys, n_cells >= MIN_CELLS, !is.na(CPS)) |>
    mutate(sex = factor(sex),
           age_s = scale(ifelse(is.na(age), median(age, na.rm = TRUE), age))[, 1],
           pmi_s = scale(ifelse(is.na(pmi), median(pmi, na.rm = TRUE), pmi))[, 1],
           cps_s = scale(CPS)[, 1])
  say("  n groups: %d", nrow(m))
  y <- DGEList(pb[, m$key])
  y <- calcNormFactors(y[filterByExpr(y, group = m$CPS > median(m$CPS)), ])
  des <- model.matrix(~ cps_s + sex + age_s + pmi_s, data = m)
  stopifnot(nrow(des) == ncol(y))
  fit <- eBayes(lmFit(voom(y, des), des))
  topTable(fit, coef = "cps_s", number = Inf, sort.by = "none") |>
    rownames_to_column("gene") |> as_tibble()
}
collect <- function(cat, subcat = NULL) {
  args <- list(species = "Homo sapiens", collection = cat)
  if (!is.null(subcat)) args$subcollection <- subcat
  do.call(msigdbr, args) |> transmute(gs_name, gene = gene_symbol)
}
gs <- bind_rows(collect("H"), collect("C5", "GO:BP"), collect("C5", "GO:CC"),
                collect("C5", "GO:MF"), collect("C2", "CP:REACTOME")) |> distinct()
gsl <- split(gs$gene, gs$gs_name)
run_g <- function(de, lab) {
  d <- de |> filter(!is.na(t)) |> arrange(desc(abs(t))) |> distinct(gene, .keep_all = TRUE)
  say("fgsea: %s (%d genes)", lab, nrow(d))
  fgsea(gsl, sort(setNames(d$t, d$gene), decreasing = TRUE),
        minSize = 10, maxSize = 500, nPermSimple = 10000) |>
    as_tibble() |> mutate(signature = lab)
}
res <- list(); des_l <- list()
for (st in STRATA) {
  keys <- meta |> filter(stratum == st) |> pull(key)
  de <- de_cps(keys)
  des_l[[st]] <- de |> mutate(signature = paste0("a9early_", st))
  res[[st]] <- run_g(de, paste0("a9early_", st))
}
g <- bind_rows(res)
write_csv(g |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = "|")),
          file.path(AOUT, "gsea_a9_strata_early.csv"))
write_csv(bind_rows(des_l), file.path(AOUT, "de_a9_strata_early.csv"))

say("=== EARLY-epoch burden (padj<0.05) ===")
print(g |> filter(padj < 0.05) |> count(signature, dir = ifelse(NES > 0, "up", "down")) |>
        pivot_wider(names_from = dir, values_from = n, values_fill = 0))
fmt <- function(NES, padj) sprintf("%.1f%s", NES, ifelse(padj < 0.01, "**",
                 ifelse(padj < 0.05, "*", ifelse(padj < 0.1, "+", ""))))
say("=== EARLY-epoch key blocks ===")
print(g |> filter(pathway %in% KEY) |> mutate(cell = fmt(NES, padj)) |>
        select(pathway, signature, cell) |>
        pivot_wider(names_from = signature, values_from = cell), n = 14, width = 200)
say("=== EARLY-epoch marquee genes (t) ===")
print(bind_rows(des_l) |>
        filter(gene %in% c("SST", "VGF", "CALB1", "NGF", "MME", "HCN1", "MAPK8", "LNX2")) |>
        select(gene, signature, t) |> mutate(t = round(t, 2)) |>
        pivot_wider(names_from = signature, values_from = t))
