#!/usr/bin/env Rscript
# GSEA on the pseudobulk stratum meta-DE (script 27a output), replacing the IVW
# shortcut. Same house fgsea settings as script 17 / the subclass pipeline.
# Also runs the S8 sensitivity: IVW-vs-pseudobulk concordance (gene z and NES).
# Writes gsea_all_signatures.csv + stratum_gene_signatures.csv in the SAME schema
# script 18 (the figure) expects, so the figure rebuild is a path switch.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(tibble); library(fgsea); library(msigdbr); library(ggplot2); library(cowplot)
})
set.seed(42)
OUT  <- "transcriptomic/results/sst_strata_gsea"
POUT <- file.path(OUT, "pseudobulk")
STRATA <- c("depleted", "intermediate", "non_depleted")
t0 <- Sys.time()
say <- function(...) cat(sprintf("[%5.1f min] ", as.numeric(difftime(Sys.time(), t0, units="mins"))), sprintf(...), "\n", sep = "")

meta <- read_csv(file.path(POUT, "stratum_meta_de.csv"), show_col_types = FALSE)

# ---- signatures in script-18 schema -------------------------------------------
sig <- meta |> filter(stratum %in% STRATA) |>
  transmute(stratum, gene, n_st = k, est_ivw = estimate, se_ivw = se, z = zval)
write_csv(sig, file.path(POUT, "stratum_gene_signatures.csv"))

# ---- gene sets ----------------------------------------------------------------
say("building gene sets ...")
collect <- function(cat, subcat = NULL) {
  args <- list(species = "Homo sapiens", collection = cat)
  if (!is.null(subcat)) args$subcollection <- subcat
  do.call(msigdbr, args) |>
    transmute(gs_name, gene = gene_symbol, source = paste(cat, subcat, sep = ":"))
}
all_gs <- bind_rows(collect("H"), collect("C5", "GO:BP"), collect("C5", "GO:CC"),
                    collect("C5", "GO:MF"), collect("C2", "CP:REACTOME")) |>
  distinct(gs_name, gene, .keep_all = TRUE)
gs_source <- all_gs |> distinct(gs_name, source)
gene_sets <- split(all_gs$gene, all_gs$gs_name)

run_gsea <- function(genes, zs, label) {
  d <- tibble(gene = genes, z = zs) |> filter(!is.na(z)) |>
    arrange(desc(abs(z))) |> distinct(gene, .keep_all = TRUE)
  ranks <- sort(setNames(d$z, d$gene), decreasing = TRUE)
  say("fgsea: %s (%d genes)", label, length(ranks))
  fgsea(gene_sets, ranks, minSize = 10, maxSize = 500, nPermSimple = 10000) |>
    as_tibble() |> mutate(signature = label)
}

gsea <- map_dfr(c(STRATA, "all_sst"), function(st) {
  d <- meta |> filter(stratum == st)
  run_gsea(d$gene, d$zval, if (st == "all_sst") "Sst_subclass" else st)
}) |> left_join(gs_source, by = c(pathway = "gs_name"))
write_csv(gsea |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = "|")),
          file.path(POUT, "gsea_all_signatures.csv"))

say("=== burden (padj<0.05) — pseudobulk meta ===")
print(gsea |> filter(padj < 0.05) |>
        count(signature, dir = ifelse(NES > 0, "up", "down")) |>
        pivot_wider(names_from = dir, values_from = n, values_fill = 0))

KEY <- c("GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING", "GOCC_PRESYNAPSE",
         "GOBP_NEUROTRANSMITTER_SECRETION", "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES",
         "REACTOME_TRANSLATION", "GOCC_RIBOSOMAL_SUBUNIT", "GOBP_CYTOPLASMIC_TRANSLATION",
         "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
         "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "GOBP_OXIDATIVE_PHOSPHORYLATION",
         "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT", "GOCC_MOTILE_CILIUM")
fmt <- function(NES, padj) sprintf("%.1f%s", NES, ifelse(padj < 0.01, "**",
                 ifelse(padj < 0.05, "*", ifelse(padj < 0.1, "+", ""))))
say("=== key blocks — pseudobulk meta ===")
print(gsea |> filter(pathway %in% KEY) |> mutate(cell = fmt(NES, padj)) |>
        select(pathway, signature, cell) |>
        pivot_wider(names_from = signature, values_from = cell) |>
        select(pathway, Sst_subclass, depleted, intermediate, non_depleted),
      n = 13, width = 200)

# ---- S8: IVW vs pseudobulk concordance -----------------------------------------
say("=== S8: IVW vs pseudobulk concordance ===")
ivw_sig <- read_csv(file.path(OUT, "stratum_gene_signatures.csv"), show_col_types = FALSE)
ivw_g   <- read_csv(file.path(OUT, "gsea_all_signatures.csv"), show_col_types = FALSE)
gz <- inner_join(ivw_sig |> select(stratum, gene, z_ivw = z),
                 sig |> select(stratum, gene, z_pb = z), by = c("stratum", "gene"))
for (st in STRATA) {
  s <- gz |> filter(stratum == st)
  say("  %s genes: r = %.2f, rho = %.2f (n = %d)", st,
      cor(s$z_ivw, s$z_pb), cor(s$z_ivw, s$z_pb, method = "spearman"), nrow(s))
}
nz <- inner_join(ivw_g |> filter(signature %in% STRATA) |> select(signature, pathway, NES_ivw = NES),
                 gsea |> filter(signature %in% STRATA) |> select(signature, pathway, NES_pb = NES),
                 by = c("signature", "pathway"))
for (st in STRATA) {
  s <- nz |> filter(signature == st)
  say("  %s NES: r = %.2f (n = %d sets)", st, cor(s$NES_ivw, s$NES_pb), nrow(s))
}
write_csv(nz, file.path(POUT, "s8_ivw_vs_pseudobulk_NES.csv"))
p <- ggplot(nz |> mutate(stratum = factor(signature, STRATA)),
            aes(NES_ivw, NES_pb)) +
  geom_point(size = 0.5, alpha = 0.25) +
  geom_abline(linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  facet_wrap(~stratum, nrow = 1) +
  labs(x = "NES, IVW shortcut", y = "NES, pseudobulk meta",
       title = "S8: GSEA concordance between the IVW shortcut and the pseudobulk meta-analysis") +
  theme_cowplot(font_size = 14) +
  theme(strip.background = element_rect(fill = "grey92"))
ggsave(file.path(POUT, "figS8_ivw_vs_pseudobulk.png"), p, width = 13, height = 4.6,
       dpi = 200, bg = "white")
say("wrote %s", file.path(POUT, "gsea_all_signatures.csv"))
