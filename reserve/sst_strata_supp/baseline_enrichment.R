#!/usr/bin/env Rscript
# Pathway enrichment on the Fig. 4g baseline signature: genes distinguishing the
# five SCZ-depleted Sst supertypes from the not-depleted ones in the SEA-AD
# neurotypical DLPFC reference (580 genes at Wilcoxon FDR < 0.05 & |log2FC| > 0.25).
# (1) ORA (hypergeometric) per direction against an expressed-gene background
#     (detected in >=10% of nuclei in either group), gene sets restricted to the
#     background, 10-500 genes, BH within direction.
# (2) Preranked fgsea on baseline log2FC over the background (identity contrast).
# (3) Cross-check: baseline log2FC of the Fig. S8 disease-state modules
#     (leading-edge unions) - do depleted types express these programs more at
#     baseline?
source("transcriptomic/scripts/sst_strata/_common.R")
suppressPackageStartupMessages({ library(msigdbr); library(fgsea) })
OUT <- file.path(P$supp, "baseline_marker_enrichment")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

v <- read_csv("genetics/results/figures/r_panels/panel_volcano_vulnerable_vs_notdepleted.csv",
              show_col_types = FALSE)
bg <- v |> filter(pct_vulnerable >= 0.10 | pct_not_depleted >= 0.10)
sig <- bg |> filter(fdr_wilcoxon < 0.05, abs(log2FC) > 0.25)
up <- sig$gene[sig$log2FC > 0]; dn <- sig$gene[sig$log2FC < 0]
cat(sprintf("background %d | signature %d (up-in-depleted %d, up-in-spared %d)\n",
            nrow(bg), nrow(sig), length(up), length(dn)))

gs <- bind_rows(msigdbr(species = "Homo sapiens", collection = "H"),
                msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"),
                msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"),
                msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:CC"),
                msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:MF")) |>
  distinct(gs_name, gene_symbol) |> filter(gene_symbol %in% bg$gene)
sets <- split(gs$gene_symbol, gs$gs_name)
sets <- sets[lengths(sets) >= 10 & lengths(sets) <= 500]
cat(sprintf("gene sets after background restriction + size filter: %d\n", length(sets)))

ora <- function(hits, label) {
  N <- nrow(bg); K <- length(hits)
  map_dfr(names(sets), function(p) {
    s <- sets[[p]]; k <- sum(hits %in% s)
    tibble(pathway = p, size = length(s), overlap = k,
           expected = K * length(s) / N,
           p = phyper(k - 1, length(s), N - length(s), K, lower.tail = FALSE),
           genes = paste(intersect(hits, s), collapse = "|"))
  }) |> mutate(padj = p.adjust(p, "BH"), direction = label) |> arrange(p)
}
res_ora <- bind_rows(ora(up, "higher_in_depleted"), ora(dn, "higher_in_spared"))
write_csv(res_ora, file.path(OUT, "ora_580_baseline_signature.csv"))
for (d in unique(res_ora$direction)) {
  cat(sprintf("\n=== ORA %s (top 12, padj<0.05 marked) ===\n", d))
  print(res_ora |> filter(direction == d) |> head(12) |>
        transmute(pathway = substr(pathway, 1, 58), overlap, size,
                  enr = round(overlap/expected, 1), padj = signif(padj, 2)), n = 12)
  cat(sprintf("  sets at padj<0.05: %d\n", sum(res_ora$padj < 0.05 & res_ora$direction == d)))
}

ranks <- setNames(bg$log2FC, bg$gene)
ranks <- ranks[!duplicated(names(ranks))]
gsea <- fgsea(sets, ranks, minSize = 10, maxSize = 500, eps = 0) |>
  as_tibble() |> arrange(padj) |>
  mutate(leadingEdge = vapply(leadingEdge, paste, "", collapse = "|"))
write_csv(gsea, file.path(OUT, "gsea_baseline_log2fc.csv"))
cat("\n=== preranked GSEA on baseline log2FC: top 12 by padj, each direction ===\n")
print(gsea |> filter(NES > 0) |> head(12) |> transmute(pathway = substr(pathway,1,58), NES = round(NES,2), padj = signif(padj,2)), n=12)
print(gsea |> filter(NES < 0) |> head(12) |> transmute(pathway = substr(pathway,1,58), NES = round(NES,2), padj = signif(padj,2)), n=12)
cat(sprintf("GSEA padj<0.05: %d positive (higher-in-depleted), %d negative\n",
            sum(gsea$padj < 0.05 & gsea$NES > 0), sum(gsea$padj < 0.05 & gsea$NES < 0)))

g5 <- read_csv(file.path(P$pb, "gsea_all_signatures.csv"), show_col_types = FALSE)
le_union <- function(ps) g5 |> filter(pathway %in% ps, signature == "depleted") |>
  pull(leadingEdge) |> strsplit("|", fixed = TRUE) |> unlist() |> unique()
MODS <- list(
  translation = le_union(c("REACTOME_TRANSLATION","GOCC_RIBOSOMAL_SUBUNIT",
                           "GOBP_CYTOPLASMIC_TRANSLATION","REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION")),
  oxphos = le_union(c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","GOBP_OXIDATIVE_PHOSPHORYLATION",
                      "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT")),
  synaptic = le_union(c("GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING","GOCC_SYNAPTIC_MEMBRANE",
                        "GOBP_NEUROTRANSMITTER_SECRETION","REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES")))
cat("\n=== baseline log2FC (depleted vs spared, neurotypical) of the disease-state modules ===\n")
cat(sprintf("background median log2FC = %.3f\n", median(bg$log2FC)))
mod_stats <- map_dfr(names(MODS), function(m) {
  x <- bg |> filter(gene %in% MODS[[m]])
  wt <- wilcox.test(x$log2FC, bg$log2FC[!bg$gene %in% MODS[[m]]])
  tibble(module = m, n = nrow(x), median_log2FC = round(median(x$log2FC), 3),
         pct_positive = round(100 * mean(x$log2FC > 0)), p_vs_background = signif(wt$p.value, 2))
})
print(mod_stats)
write_csv(mod_stats, file.path(OUT, "modules_baseline_log2fc.csv"))
cat("\nwrote", OUT, "\n")
