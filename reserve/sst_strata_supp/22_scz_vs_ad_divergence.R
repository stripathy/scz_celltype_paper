#!/usr/bin/env Rscript
# Three jobs, no new DE:
# 1) SCZ-depleted vs AD-depleted divergence figure (pathway NES + gene level)
# 2) Audit of how the highlighted pathway blocks were selected: classify ALL
#    SCZ depleted-stratum significant sets by leading-edge family and report
#    whether the highlighted blocks are representative; write the full
#    annotated table (supplementary-table-ready).
# 3) Are Gabitto's marquee genes (MME, NGF, MAPK8, LNX2, RPL4, ATP5MPL...)
#    present/tested in our SCZ analyses?
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(tibble)
  library(ggplot2); library(cowplot); library(ggrepel)
})
OUT  <- "transcriptomic/results/sst_strata_gsea"
AOUT <- file.path(OUT, "seaad_a9")
FS <- 14

g_scz <- read_csv(file.path(OUT, "gsea_all_signatures.csv"), show_col_types = FALSE)
g_a9  <- read_csv(file.path(AOUT, "gsea_a9_strata.csv"), show_col_types = FALSE)
sig   <- read_csv(file.path(OUT, "stratum_gene_signatures.csv"), show_col_types = FALSE)
de_a9 <- read_csv(file.path(AOUT, "de_a9_strata.csv"), show_col_types = FALSE)

# ---- module gene lists (leading edges of the block sets, SCZ depleted) --------
blocks <- list(
  Synaptic = c("GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING", "GOCC_PRESYNAPSE",
               "GOBP_NEUROTRANSMITTER_SECRETION", "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES"),
  Translation = c("REACTOME_TRANSLATION", "GOCC_RIBOSOMAL_SUBUNIT",
                  "GOBP_CYTOPLASMIC_TRANSLATION", "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION"),
  `OxPhos/mito` = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "GOBP_OXIDATIVE_PHOSPHORYLATION",
                    "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"))
le_union <- function(pws) g_scz |> filter(pathway %in% pws, signature == "depleted") |>
  pull(leadingEdge) |> strsplit("|", fixed = TRUE) |> unlist() |> unique()
mods <- lapply(blocks, le_union)
fam_of <- function(le, thr = 0.5) {
  gg <- strsplit(le, "|", fixed = TRUE)[[1]]
  fr <- vapply(mods, function(m) mean(gg %in% m), numeric(1))
  if (max(fr) >= thr) names(which.max(fr)) else "Other"
}
FAM_COLS <- c(Synaptic = "#0072B2", Translation = "#CC79A7",
              `OxPhos/mito` = "#D55E00", Other = "grey78")

# ---- 2) selection audit -------------------------------------------------------
dep <- g_scz |> filter(signature == "depleted") |>
  mutate(family = vapply(leadingEdge, fam_of, character(1)))
sig_dn <- dep |> filter(padj < 0.05, NES < 0) |> arrange(padj)
cat("=== Selection audit: SCZ depleted-stratum significant DOWN sets (n=",
    nrow(sig_dn), ") by leading-edge family (>=50% overlap) ===\n", sep = "")
print(count(sig_dn, family) |> mutate(pct = round(100 * n / sum(n))))
cat("\nTop 20 by padj, with family:\n")
print(sig_dn |> transmute(rank = row_number(), pathway = str_trunc(pathway, 55),
                          NES = round(NES, 2), padj = signif(padj, 2), family) |>
        head(20), n = 20)
cat(sprintf("\nTop-20 coverage by the three highlighted families: %d/20\n",
            sum(head(sig_dn$family, 20) != "Other")))
cat(sprintf("Top-30 coverage: %d/30\n", sum(head(sig_dn$family, 30) != "Other")))
# full annotated table across strata (supp-table-ready)
wide <- g_scz |> filter(signature %in% c("depleted", "intermediate", "non_depleted")) |>
  select(pathway, source, signature, NES, padj) |>
  pivot_wider(names_from = signature, values_from = c(NES, padj)) |>
  left_join(dep |> select(pathway, family), by = "pathway") |>
  arrange(padj_depleted)
write_csv(wide, file.path(OUT, "suppl_table_strata_gsea_annotated.csv"))

# ---- 3) Gabitto marquee genes in our analyses ---------------------------------
marquee <- c("MME", "NGF", "NGFR", "MAPK8", "LNX2", "RPL4", "ATP5MPL", "HCN1")
sub_de <- read_csv("shared/snrnaseq_de/DE_genes_all_cells_scz.csv", show_col_types = FALSE) |>
  filter(cell_type == "Sst") |> mutate(z = estimate / se) |>
  select(gene = genes, z_scz_subclass = z)
xen_panel <- unique(read_csv("spatial/output/de/de_results_supertype.csv",
                             show_col_types = FALSE)$gene)
tab_g <- tibble(gene = marquee) |>
  left_join(sig |> select(gene, stratum, z) |>
              pivot_wider(names_from = stratum, values_from = z,
                          names_prefix = "z_scz_"), by = "gene") |>
  left_join(sub_de, by = "gene") |>
  left_join(de_a9 |> filter(signature == "a9_depleted") |>
              select(gene, t_a9_depleted = t), by = "gene") |>
  mutate(on_xenium_panel = gene %in% xen_panel,
         across(where(is.numeric), ~round(.x, 2)))
cat("\n=== Gabitto marquee genes in our analyses ===\n")
print(tab_g, n = 10, width = 200)
write_csv(tab_g, file.path(AOUT, "gabitto_marquee_genes_in_scz.csv"))

# ---- 1) divergence figure -----------------------------------------------------
pw <- dep |> select(pathway, family, NES_scz = NES, padj_scz = padj) |>
  inner_join(g_a9 |> filter(signature == "a9_depleted") |>
               select(pathway, NES_a9 = NES, padj_a9 = padj), by = "pathway")
r_pw <- cor(pw$NES_scz, pw$NES_a9)
lab_pw <- pw |> filter(family != "Other") |> group_by(family) |>
  slice_min(padj_scz, n = 2) |> ungroup() |>
  mutate(lab = str_trunc(str_to_sentence(str_replace_all(
    str_remove(pathway, "^(GOBP|GOCC|GOMF|REACTOME|HALLMARK)_"), "_", " ")), 34))
p1 <- ggplot(pw, aes(NES_scz, NES_a9, colour = family)) +
  geom_point(data = ~ filter(.x, family == "Other"), size = 0.7, alpha = 0.3) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_abline(linetype = "dashed", colour = "grey60") +
  geom_point(data = ~ filter(.x, family != "Other"), size = 1.9) +
  geom_text_repel(data = lab_pw, aes(label = lab), size = 3.6, seed = 2,
                  max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = FAM_COLS, name = NULL) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 4.8,
           label = sprintf("r = %.2f (n = %s gene sets)", r_pw, format(nrow(pw), big.mark = ","))) +
  labs(x = "GSEA NES, SCZ depleted stratum (case vs control)",
       y = "GSEA NES, AD depleted stratum (along CPS)") +
  theme_cowplot(font_size = FS) + theme(legend.position = "bottom")

gn <- sig |> filter(stratum == "depleted") |> select(gene, z_scz = z) |>
  inner_join(de_a9 |> filter(signature == "a9_depleted") |> select(gene, t_a9 = t),
             by = "gene") |>
  mutate(module = case_when(gene %in% mods$Translation ~ "Translation",
                            gene %in% mods$`OxPhos/mito` ~ "OxPhos/mito",
                            gene %in% mods$Synaptic ~ "Synaptic",
                            TRUE ~ "Other"),
         module = factor(module, names(FAM_COLS)))
r_gn <- cor(gn$z_scz, gn$t_a9)
hl <- gn |> filter(gene %in% c("SST", "VGF", "CALB1", "MME", "NGF", "HCN1",
                               "RPL36", "RPL10", "NDUFS8", "UQCRH", "NTRK2"))
p2 <- ggplot(gn, aes(z_scz, t_a9, colour = module)) +
  geom_point(data = ~ filter(.x, module == "Other"), size = 0.5, alpha = 0.2) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(data = ~ filter(.x, module != "Other"), size = 1.4, alpha = 0.8) +
  geom_point(data = hl, colour = "black", size = 2.2, shape = 21, stroke = 0.8) +
  geom_text_repel(data = hl, aes(label = gene), colour = "black", size = 4,
                  fontface = "italic", seed = 3, max.overlaps = 25) +
  scale_colour_manual(values = FAM_COLS, name = NULL) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 4.8,
           label = sprintf("r = %.2f (n = %s genes)", r_gn, format(nrow(gn), big.mark = ","))) +
  labs(x = "Gene z, SCZ depleted stratum (case vs control)",
       y = "Gene t, AD depleted stratum (along CPS)") +
  theme_cowplot(font_size = FS) + theme(legend.position = "bottom")

fig <- plot_grid(p1, p2, nrow = 1, labels = c("a", "b"), label_size = 22)
ggsave(file.path(AOUT, "figS_scz_vs_ad_depleted_divergence.png"), fig,
       width = 14.5, height = 6.8, dpi = 200, bg = "white")
write_csv(pw, file.path(AOUT, "scz_vs_ad_depleted_pathways.csv"))
cat("\nwrote figS_scz_vs_ad_depleted_divergence.png\n")
