#!/usr/bin/env Rscript
# Swap-in: replace the IVW stratum signatures with the pseudobulk meta-DE
# (from 27a), rerun GSEA, quantify IVW-vs-pseudobulk concordance (S8), splice
# the canonical files (archiving the IVW versions), and print the headline
# tables. Then script 18 can be rerun unchanged to refresh the figure.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(tibble); library(fgsea); library(msigdbr)
})
set.seed(42)
OUT <- "transcriptomic/results/sst_strata_gsea"
PB  <- file.path(OUT, "pseudobulk")
ARCH <- file.path(OUT, "archive_ivw")
dir.create(ARCH, showWarnings = FALSE)
STRATA <- c("depleted", "intermediate", "non_depleted")
t0 <- Sys.time()
say <- function(...) cat(sprintf("[%5.1f min] ", as.numeric(difftime(Sys.time(), t0, units="mins"))), sprintf(...), "\n", sep = "")

meta <- read_csv(file.path(PB, "stratum_meta_de.csv"), show_col_types = FALSE)

# ---- archive IVW canonical files (once) ---------------------------------------
for (f in c("stratum_gene_signatures.csv", "gsea_all_signatures.csv")) {
  if (!file.exists(file.path(ARCH, f)))
    file.copy(file.path(OUT, f), file.path(ARCH, f))
}
sig_ivw <- read_csv(file.path(ARCH, "stratum_gene_signatures.csv"), show_col_types = FALSE)
g_ivw   <- read_csv(file.path(ARCH, "gsea_all_signatures.csv"), show_col_types = FALSE)

# ---- new canonical gene signatures ---------------------------------------------
sig_pb <- meta |> filter(stratum %in% STRATA) |>
  transmute(stratum, gene, z = zval, estimate, se, k, padj)
write_csv(sig_pb, file.path(OUT, "stratum_gene_signatures.csv"))

# S8: gene-level IVW vs pseudobulk concordance
say("S8 gene-level concordance (IVW vs pseudobulk z):")
for (st in STRATA) {
  m <- sig_ivw |> filter(stratum == st) |> select(gene, z_ivw = z) |>
    inner_join(sig_pb |> filter(stratum == st) |> select(gene, z_pb = z), by = "gene")
  say("  %s: r = %.2f, rho = %.2f (n = %d)", st,
      cor(m$z_ivw, m$z_pb), cor(m$z_ivw, m$z_pb, method = "spearman"), nrow(m))
}

# ---- GSEA on the new signatures -------------------------------------------------
collect <- function(cat, subcat = NULL) {
  args <- list(species = "Homo sapiens", collection = cat)
  if (!is.null(subcat)) args$subcollection <- subcat
  do.call(msigdbr, args) |> transmute(gs_name, gene = gene_symbol)
}
all_gs <- bind_rows(collect("H"), collect("C5", "GO:BP"), collect("C5", "GO:CC"),
                    collect("C5", "GO:MF"), collect("C2", "CP:REACTOME")) |>
  distinct(gs_name, gene)
gs_source <- g_ivw |> distinct(pathway, source) |> rename(gs_name = pathway)
gene_sets <- split(all_gs$gene, all_gs$gs_name)
run_gsea <- function(genes, zs, label) {
  d <- tibble(gene = genes, z = zs) |> filter(!is.na(z)) |>
    arrange(desc(abs(z))) |> distinct(gene, .keep_all = TRUE)
  say("fgsea: %s (%d genes)", label, nrow(d))
  fgsea(gene_sets, sort(setNames(d$z, d$gene), decreasing = TRUE),
        minSize = 10, maxSize = 500, nPermSimple = 10000) |>
    as_tibble() |> mutate(signature = label)
}
g_new <- map_dfr(c(STRATA, "all_sst"), function(st) {
  d <- meta |> filter(stratum == st)
  run_gsea(d$gene, d$zval, ifelse(st == "all_sst", "all_sst_pb", st))
})
g_new <- g_new |> left_join(gs_source, by = c(pathway = "gs_name"))

# splice: keep every non-strata signature from the IVW file (Sst_subclass anchor,
# nicole_* groups); replace the three strata; add all_sst_pb
g_keep <- g_ivw |> filter(!signature %in% c(STRATA, "all_sst_pb"))
g_out <- bind_rows(g_keep,
                   g_new |> mutate(leadingEdge = sapply(leadingEdge, pa