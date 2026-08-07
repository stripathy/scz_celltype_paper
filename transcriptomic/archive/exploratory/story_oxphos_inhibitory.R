#!/usr/bin/env Rscript
# Story 1: OxPhos collapse in inhibitory neurons
# 5-panel exploratory figure designed to spark hypotheses.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(cowplot); library(fgsea); library(msigdbr); library(ggrepel)
})

OUT_DIR <- "results/exploratory/story_oxphos"
dir.create(OUT_DIR, showWarnings = FALSE)
GWAS_FILE <- "~/Github/scz_cell_type_enrichment/data/gwas/scz_gwas_gene_set_no_mhc.csv"

EXC <- c("L2_3 IT","L4 IT","L5 IT","L5 ET","L5_6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")
INH <- c("Lamp5","Pax6","Pvalb","Sncg","Sst","Sst Chodl","Vip","Chandelier")
GLI <- c("Astro","Oligo","OPC","Micro-PVM","Endo","VLMC")
order_ct <- function(cts) {
  ord <- c(intersect(EXC, cts), intersect(INH, cts), intersect(GLI, cts))
  factor(cts, levels = ord)
}
class_of <- function(ct) case_when(ct %in% EXC ~ "Excitatory",
                                    ct %in% INH ~ "Inhibitory",
                                    ct %in% GLI ~ "Glia",
                                    TRUE        ~ NA_character_)
CLASS_COL <- c(Excitatory = "#117733", Inhibitory = "#882255", Glia = "#DDCC77")

# Load DE data
df <- read_csv("data/DE_genes_all_cells_scz.csv", show_col_types = FALSE) |>
  filter(!is.na(estimate), !is.na(se), se > 0) |>
  mutate(z = estimate / se)
gwas_set <- read_csv(GWAS_FILE, show_col_types = FALSE)$gene

# Existing GSEA cache (broad analysis)
prior <- readRDS("results/gsea_cache.rds")
prior_long <- bind_rows(prior)

# ----- Build fine-grained mitochondrial gene sets ---------------------------
cat("Pulling fine-grained mitochondrial gene sets...\n")
collect <- function(cat, subcat = NULL) {
  args <- list(species = "Homo sapiens", collection = cat)
  if (!is.null(subcat)) args$subcollection <- subcat
  do.call(msigdbr, args) |> transmute(gs_name, gene = gene_symbol)
}
all_msigdb <- bind_rows(
  collect("H"),
  collect("C5","GO:BP"), collect("C5","GO:CC"), collect("C5","GO:MF"),
  collect("C2","CP:REACTOME"), collect("C2","CP:KEGG_LEGACY"),
  collect("C2","CP:KEGG_MEDICUS")
) |> distinct(gs_name, gene)

mito_sets <- list(
  "Complex I (NADH:UQ)"   = c("GOCC_NADH_DEHYDROGENASE_COMPLEX",
                              "REACTOME_COMPLEX_I_BIOGENESIS",
                              "GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE",
                              "GOBP_NADH_DEHYDROGENASE_COMPLEX_ASSEMBLY"),
  "Complex II (SDH)"      = c("GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_SUCCINATE_TO_UBIQUINONE"),
  "Complex III (UQ:Cyt c)"= c("GOCC_RESPIRATORY_CHAIN_COMPLEX_III",
                              "GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_UBIQUINOL_TO_CYTOCHROME_C",
                              "REACTOME_COMPLEX_III_ASSEMBLY"),
  "Complex IV (CytOx)"    = c("GOCC_RESPIRATORY_CHAIN_COMPLEX_IV",
                              "REACTOME_COMPLEX_IV_ASSEMBLY",
                              "GOBP_RESPIRATORY_CHAIN_COMPLEX_IV_ASSEMBLY"),
  "Complex V (ATP syn)"   = c("GOCC_PROTON_TRANSPORTING_ATP_SYNTHASE_COMPLEX",
                              "GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",
                              "GOBP_MITOCHONDRIAL_PROTON_TRANSPORTING_ATP_SYNTHASE_COMPLEX_ASSEMBLY",
                              "GOBP_PROTON_MOTIVE_FORCE_DRIVEN_ATP_SYNTHESIS"),
  "TCA cycle"           = c("REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
                            "GOBP_TRICARBOXYLIC_ACID_CYCLE",
                            "KEGG_CITRATE_CYCLE_TCA_CYCLE"),
  "Mito translation"    = c("GOBP_MITOCHONDRIAL_TRANSLATION",
                            "REACTOME_MITOCHONDRIAL_TRANSLATION",
                            "GOCC_MITOCHONDRIAL_RIBOSOME"),
  "Mito import / fold"  = c("REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT",
                            "GOBP_PROTEIN_TARGETING_TO_MITOCHONDRION",
                            "REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION"),
  "Inner membrane"      = c("GOCC_MITOCHONDRIAL_INNER_MEMBRANE",
                            "GOCC_ORGANELLE_INNER_MEMBRANE"),
  "Aerobic resp / OxPhos summary" = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                                       "GOBP_AEROBIC_RESPIRATION",
                                       "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT")
)

# For each label, build a union gene set from available MSigDB pathways
build_union <- function(names) {
  hits <- all_msigdb |> filter(gs_name %in% names)
  if (nrow(hits) == 0) return(character(0))
  unique(hits$gene)
}
gs_for_mito <- map(mito_sets, build_union)
gs_for_mito <- gs_for_mito[map_int(gs_for_mito, length) >= 8]
cat("Mitochondrial subprocess gene set sizes:\n")
print(map_int(gs_for_mito, length))

# Run fgsea per cell type for these focused sets
cts <- sort(unique(df$cell_type))
mito_results <- list()
set.seed(7)
for (ct in cts) {
  sub <- df |> filter(cell_type == ct) |>
    arrange(desc(abs(z))) |> distinct(genes, .keep_all = TRUE)
  if (nrow(sub) < 500) next
  ranks <- setNames(sub$z, sub$genes)
  ranks <- sort(ranks, decreasing = TRUE)
  res <- fgsea(pathways = gs_for_mito, stats = ranks,
               minSize = 5, maxSize = 2000, nPermSimple = 20000)
  res$cell_type <- ct
  mito_results[[ct]] <- as_tibble(res)
}
mito_long <- bind_rows(mito_results) |>
  mutate(pathway = factor(pathway, levels = names(gs_for_mito)))
write_csv(mito_long |> select(-leadingEdge),
          file.path(OUT_DIR, "mito_subprocess_gsea.csv"))

# =========================================================================
# Panel A — Forest plot of HALLMARK_OXPHOS NES across all cell types
# =========================================================================
oxphos_df <- prior_long |>
  filter(pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") |>
  mutate(class = class_of(cell_type),
         cell_type_ord = reorder(cell_type, NES))

p_A <- ggplot(oxphos_df, aes(x = NES, y = cell_type_ord, colour = class)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = NES, yend = cell_type_ord), linewidth = 0.6) +
  geom_point(aes(size = -log10(padj + 1e-6))) +
  geom_text(aes(label = ifelse(padj < 0.001, "***",
                        ifelse(padj < 0.01,  "**",
                        ifelse(padj < 0.05,  "*",
                        ifelse(padj < 0.1,   "+", ""))))),
            colour = "black", hjust = -0.5, vjust = 0.4, size = 4) +
  scale_colour_manual(values = CLASS_COL, name = "Cell class") +
  scale_size_continuous(name = expression(-log[10]~"FDR"), range = c(2, 7)) +
  scale_x_continuous(limits = c(-2.5, 2.0)) +
  labs(x = "GSEA NES (HALLMARK_OXIDATIVE_PHOSPHORYLATION)", y = NULL,
       subtitle = "A. OxPhos enrichment ranks inhibitory neurons (Sst, Vip, Lamp5, Pvalb) at the most-down end") +
  theme_cowplot(font_size = 12) +
  theme(plot.subtitle  = element_text(size = 12, face = "bold"),
        axis.text.y    = element_text(size = 11),
        legend.position = "right")

# =========================================================================
# Panel B — Mitochondrial subprocess heatmap
# =========================================================================
mito_long$cell_type <- order_ct(mito_long$cell_type)

p_B <- ggplot(mito_long, aes(x = cell_type, y = pathway, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(padj < 0.001, "***",
                        ifelse(padj < 0.01,  "**",
                        ifelse(padj < 0.05,  "*",
                        ifelse(padj < 0.1,   "+", ""))))),
            size = 3.5, vjust = 0.7) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, name = "NES",
                       limits = c(-3, 3), oob = scales::squish) +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL,
       subtitle = "B. Mitochondrial subprocess decomposition - which part of mito biology is hit?") +
  theme_cowplot(font_size = 12) +
  theme(plot.subtitle = element_text(size = 12, face = "bold"),
        axis.text.x   = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y   = element_text(size = 10),
        panel.grid    = element_blank(),
        axis.line     = element_blank(), axis.ticks = element_blank())

# =========================================================================
# Panel C — Leading-edge gene heatmap (z-scores)
# =========================================================================
# Get leading-edge OxPhos genes from prior cache for top inhibitory cell types
key_inh <- c("Sst","Vip","Lamp5","Pvalb","Sncg","Pax6","Chandelier")
le_per_ct <- prior_long |>
  filter(pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
         cell_type %in% key_inh) |>
  rowwise() |>
  mutate(le = list(leadingEdge)) |>
  ungroup() |>
  select(cell_type, le) |>
  unnest(le) |>
  count(le, sort = TRUE) |>
  rename(gene = le, n_inh_le = n) |>
  filter(n_inh_le >= 2)  # in leading edge of >= 2 inhibitory cell types

top_le <- le_per_ct |> head(30) |> pull(gene)
cat(sprintf("\nTop %d leading-edge OxPhos genes (>=2 inhibitory cell types)\n", length(top_le)))

# Build z-score heatmap for these genes across all neurons + glia
le_z <- df |> filter(genes %in% top_le) |>
  select(cell_type, genes, z, padj) |>
  mutate(cell_type = order_ct(cell_type),
         genes     = factor(genes, levels = top_le),
         in_gwas   = genes %in% gwas_set)

le_z$genes <- factor(le_z$genes, levels = rev(top_le))
p_C <- ggplot(le_z, aes(x = cell_type, y = genes, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(padj < 0.05, "*",
                        ifelse(padj < 0.1,  "+", ""))),
            size = 3, vjust = 0.7) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, name = "DE z",
                       limits = c(-4, 4), oob = scales::squish) +
  labs(x = NULL, y = NULL,
       subtitle = sprintf("C. Top %d OxPhos leading-edge genes (in >=2 inhibitory cell types) - z-scores",
                          length(top_le))) +
  theme_cowplot(font_size = 11) +
  theme(plot.subtitle = element_text(size = 12, face = "bold"),
        axis.text.x   = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y   = element_text(size = 9, face = ifelse(rev(top_le) %in% gwas_set,
                                                              "bold.italic", "plain")),
        panel.grid    = element_blank(),
        axis.line     = element_blank(), axis.ticks = element_blank())

# =========================================================================
# Panel D — Focal gene forest plots (GOT2 + GWAS-overlapping OxPhos genes)
# =========================================================================
focal_genes <- top_le[1:6]
if (!"GOT2" %in% focal_genes) focal_genes <- c("GOT2", focal_genes[1:5])
focal_genes <- unique(focal_genes)
cat(sprintf("Focal genes for panel D: %s\n", paste(focal_genes, collapse = ", ")))
cat(sprintf("  (GWAS-overlapping: %s)\n",
            paste(intersect(focal_genes, gwas_set), collapse = ", ")))

focal_df <- df |> filter(genes %in% focal_genes) |>
  mutate(cell_type = order_ct(cell_type),
         class     = class_of(as.character(cell_type)),
         genes     = factor(genes, levels = focal_genes),
         in_gwas   = genes %in% gwas_set,
         ci_lo     = estimate - 1.96*se,
         ci_hi     = estimate + 1.96*se,
         sig       = padj < 0.1)

p_D <- ggplot(focal_df, aes(x = estimate, y = cell_type, colour = class)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, yend = cell_type), linewidth = 0.4) +
  geom_point(aes(shape = sig), size = 2) +
  facet_wrap(~ genes, ncol = 3, scales = "free_x") +
  scale_colour_manual(values = CLASS_COL, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(x = "Effect estimate (SCZ vs control)", y = NULL,
       subtitle = sprintf("D. Top OxPhos leading-edge genes (filled = padj<0.1). %s in PGC3 GWAS.",
                           paste(intersect(focal_genes, gwas_set), collapse = ", "))) +
  theme_cowplot(font_size = 11) +
  theme(plot.subtitle = element_text(size = 12, face = "bold"),
        axis.text.y   = element_text(size = 9),
        strip.text    = element_text(face = "bold.italic", size = 12),
        strip.background = element_blank())

# =========================================================================
# Panel E — Cross-pathway correlation: OxPhos NES vs Synaptic NES across cells
# =========================================================================
xy <- prior_long |>
  filter(pathway %in% c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                        "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING",
                        "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
                        "REACTOME_TRANSLATION")) |>
  select(cell_type, pathway, NES) |>
  pivot_wider(names_from = pathway, values_from = NES) |>
  mutate(class = class_of(cell_type)) |>
  filter(!is.na(HALLMARK_OXIDATIVE_PHOSPHORYLATION),
         !is.na(GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING))

cor_test <- cor.test(xy$HALLMARK_OXIDATIVE_PHOSPHORYLATION,
                     xy$GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING,
                     method = "spearman")
cat(sprintf("\nSpearman rho(OxPhos vs trans-syn-signaling) across cell types: %.2f (p=%.3g)\n",
            cor_test$estimate, cor_test$p.value))

p_E <- ggplot(xy, aes(x = HALLMARK_OXIDATIVE_PHOSPHORYLATION,
                       y = GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING,
                       colour = class)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(aes(label = cell_type), size = 3.2,
                  max.overlaps = Inf, seed = 1,
                  box.padding = 0.3, point.padding = 0.2,
                  segment.colour = "grey60", show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey40",
              linetype = "dashed", linewidth = 0.4, inherit.aes = FALSE,
              aes(x = HALLMARK_OXIDATIVE_PHOSPHORYLATION,
                  y = GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING)) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Spearman rho = %.2f  (p = %.2g)",
                           cor_test$estimate, cor_test$p.value),
           hjust = -0.05, vjust = 1.5, size = 4, colour = "black") +
  scale_colour_manual(values = CLASS_COL, name = "Class") +
  labs(x = "NES (HALLMARK OxPhos)",
       y = "NES (trans-synaptic signaling regulation)",
       subtitle = "E. Are OxPhos and synaptic collapse coupled across cell types?") +
  theme_cowplot(font_size = 12) +
  theme(plot.subtitle = element_text(size = 12, face = "bold"))

# =========================================================================
# Compose figure
# =========================================================================
row1 <- plot_grid(p_A, p_E, nrow = 1, rel_widths = c(1.3, 1))
row2 <- p_B
row3 <- plot_grid(p_C, p_D, nrow = 1, rel_widths = c(1.4, 1.1))

full <- plot_grid(row1, row2, row3, ncol = 1,
                  rel_heights = c(1, 1.05, 1.65))
ggsave(file.path(OUT_DIR, "figure_oxphos_inhibitory.png"),
       full, width = 18, height = 21, dpi = 200, bg = "white",
       limitsize = FALSE)
ggsave(file.path(OUT_DIR, "figure_oxphos_inhibitory.pdf"),
       full, width = 18, height = 21, bg = "white",
       limitsize = FALSE)
cat(sprintf("\nWrote %s/figure_oxphos_inhibitory.{png,pdf}\n", OUT_DIR))
