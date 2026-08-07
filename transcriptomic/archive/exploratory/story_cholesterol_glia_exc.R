#!/usr/bin/env Rscript
# Story 2: Cholesterol biosynthesis collapse in glia + excitatory neurons
# 5-panel exploratory figure with mevalonate pathway breakdown + lipid context

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(cowplot); library(fgsea); library(msigdbr); library(ggrepel)
})

OUT_DIR <- "results/exploratory/story_cholesterol"
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

df <- read_csv("data/DE_genes_all_cells_scz.csv", show_col_types = FALSE) |>
  filter(!is.na(estimate), !is.na(se), se > 0) |>
  mutate(z = estimate / se)
gwas_set <- read_csv(GWAS_FILE, show_col_types = FALSE)$gene

prior <- readRDS("results/gsea_cache.rds")
prior_long <- bind_rows(prior)

# ----- Build fine-grained lipid / sterol gene sets --------------------------
cat("Pulling lipid metabolism gene sets...\n")
collect <- function(cat, subcat = NULL) {
  args <- list(species = "Homo sapiens", collection = cat)
  if (!is.null(subcat)) args$subcollection <- subcat
  do.call(msigdbr, args) |> transmute(gs_name, gene = gene_symbol)
}
all_msigdb <- bind_rows(
  collect("H"),
  collect("C5","GO:BP"), collect("C5","GO:CC"), collect("C5","GO:MF"),
  collect("C2","CP:REACTOME"), collect("C2","CP:KEGG_LEGACY"),
  collect("C2","CP:WIKIPATHWAYS")
) |> distinct(gs_name, gene)

lipid_sets <- list(
  "Mevalonate / isoprenoid"   = c("REACTOME_CHOLESTEROL_BIOSYNTHESIS",
                                  "GOBP_ISOPRENOID_BIOSYNTHETIC_PROCESS",
                                  "WP_CHOLESTEROL_BIOSYNTHESIS_PATHWAY",
                                  "KEGG_STEROID_BIOSYNTHESIS"),
  "Cholesterol biosynthesis"  = c("REACTOME_CHOLESTEROL_BIOSYNTHESIS",
                                  "GOBP_CHOLESTEROL_BIOSYNTHETIC_PROCESS"),
  "Sterol biosynthesis"       = c("GOBP_STEROL_BIOSYNTHETIC_PROCESS",
                                  "GOBP_STEROL_METABOLIC_PROCESS"),
  "SREBP regulation"          = c("REACTOME_REGULATION_OF_CHOLESTEROL_BIOSYNTHESIS_BY_SREBP_SREBF",
                                  "GOBP_REGULATION_OF_CHOLESTEROL_BIOSYNTHETIC_PROCESS"),
  "Cholesterol homeostasis"   = c("HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                                  "GOBP_CHOLESTEROL_HOMEOSTASIS"),
  "Steroid metabolism"        = c("REACTOME_METABOLISM_OF_STEROIDS",
                                  "GOBP_STEROID_METABOLIC_PROCESS"),
  "Fatty acid metabolism"     = c("HALLMARK_FATTY_ACID_METABOLISM",
                                  "GOBP_FATTY_ACID_METABOLIC_PROCESS"),
  "Sphingolipid metabolism"   = c("REACTOME_SPHINGOLIPID_METABOLISM",
                                  "GOBP_SPHINGOLIPID_METABOLIC_PROCESS"),
  "Phospholipid metabolism"   = c("REACTOME_PHOSPHOLIPID_METABOLISM",
                                  "GOBP_PHOSPHOLIPID_METABOLIC_PROCESS"),
  "Lipid transport"           = c("GOBP_LIPID_TRANSPORT",
                                  "REACTOME_PLASMA_LIPOPROTEIN_REMODELING",
                                  "REACTOME_LIPOPROTEIN_METABOLISM"),
  "Myelination"               = c("GOBP_MYELINATION",
                                  "GOBP_AXON_ENSHEATHMENT")
)

build_union <- function(names) {
  hits <- all_msigdb |> filter(gs_name %in% names)
  if (nrow(hits) == 0) return(character(0))
  unique(hits$gene)
}
gs_for_lipid <- map(lipid_sets, build_union)
gs_for_lipid <- gs_for_lipid[map_int(gs_for_lipid, length) >= 8]
cat("Lipid subprocess gene set sizes:\n")
print(map_int(gs_for_lipid, length))

cts <- sort(unique(df$cell_type))
lipid_results <- list()
set.seed(11)
for (ct in cts) {
  sub <- df |> filter(cell_type == ct) |>
    arrange(desc(abs(z))) |> distinct(genes, .keep_all = TRUE)
  if (nrow(sub) < 500) next
  ranks <- setNames(sub$z, sub$genes)
  ranks <- sort(ranks, decreasing = TRUE)
  res <- fgsea(pathways = gs_for_lipid, stats = ranks,
               minSize = 5, maxSize = 2000, nPermSimple = 20000)
  res$cell_type <- ct
  lipid_results[[ct]] <- as_tibble(res)
}
lipid_long <- bind_rows(lipid_results) |>
  mutate(pathway = factor(pathway, levels = names(gs_for_lipid)))
write_csv(lipid_long |> select(-leadingEdge),
          file.path(OUT_DIR, "lipid_subprocess_gsea.csv"))

# =========================================================================
# Panel A — Forest plot of REACTOME_CHOLESTEROL_BIOSYNTHESIS NES
# =========================================================================
chol_df <- prior_long |>
  filter(pathway == "REACTOME_CHOLESTEROL_BIOSYNTHESIS") |>
  mutate(class = class_of(cell_type),
         cell_type_ord = reorder(cell_type, NES))

p_A <- ggplot(chol_df, aes(x = NES, y = cell_type_ord, colour = class)) +
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
  scale_x_continuous(limits = c(-3, 2)) +
  labs(x = "GSEA NES (REACTOME_CHOLESTEROL_BIOSYNTHESIS)", y = NULL,
       subtitle = "A. Cholesterol biosynthesis strongest-down in upper/mid excitatory (L5 IT, L5_6 NP, L6 IT) + Astro/Oligo") +
  theme_cowplot(font_size = 12) +
  theme(plot.subtitle  = element_text(size = 12, face = "bold"),
        axis.text.y    = element_text(size = 11),
        legend.position = "right")

# =========================================================================
# Panel B — Lipid subprocess heatmap
# =========================================================================
lipid_long$cell_type <- order_ct(lipid_long$cell_type)

p_B <- ggplot(lipid_long, aes(x = cell_type, y = pathway, fill = NES)) +
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
       subtitle = "B. Lipid subprocess decomposition - is it specifically cholesterol biosynthesis, or broader lipid metabolism?") +
  theme_cowplot(font_size = 12) +
  theme(plot.subtitle = element_text(size = 12, face = "bold"),
        axis.text.x   = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y   = element_text(size = 10),
        panel.grid    = element_blank(),
        axis.line     = element_blank(), axis.ticks = element_blank())

# =========================================================================
# Panel C — Leading-edge gene heatmap (mevalonate pathway, z-scores)
# =========================================================================
key_affected <- c("Astro","Oligo","L5 IT","L5_6 NP","L6 IT","L6 CT","L4 IT","L2_3 IT")
le_chol <- prior_long |>
  filter(pathway == "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
         cell_type %in% key_affected) |>
  rowwise() |>
  mutate(le = list(leadingEdge)) |>
  ungroup() |>
  select(cell_type, le) |>
  unnest(le) |>
  count(le, sort = TRUE) |>
  rename(gene = le, n_le = n) |>
  filter(n_le >= 2)

top_chol_le <- le_chol |> head(25) |> pull(gene)
cat(sprintf("\nTop %d cholesterol leading-edge genes (>=2 affected cell types)\n",
            length(top_chol_le)))

# Order genes by mevalonate pathway step (manual annotation of canonical enzymes)
pathway_order <- c(
  "ACAT2","HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","GGPS1",
  "FDFT1","SQLE","LSS","CYP51A1","TM7SF2","MSMO1","NSDHL","HSD17B7",
  "EBP","SC5D","DHCR24","DHCR7",
  "SREBF1","SREBF2","INSIG1","INSIG2","HMGCR","SCAP","NPC1","NPC2",
  "LDLR","ABCA1","SOAT1"
)
ordered_top <- c(intersect(pathway_order, top_chol_le),
                 setdiff(top_chol_le, pathway_order))

le_z <- df |> filter(genes %in% ordered_top) |>
  select(cell_type, genes, z, padj) |>
  mutate(cell_type = order_ct(cell_type),
         genes     = factor(genes, levels = rev(ordered_top)))

p_C <- ggplot(le_z, aes(x = cell_type, y = genes, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(padj < 0.05, "*",
                        ifelse(padj < 0.1,  "+", ""))),
            size = 3, vjust = 0.7) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, name = "DE z",
                       limits = c(-4, 4), oob = scales::squish) +
  labs(x = NULL, y = NULL,
       subtitle = sprintf("C. Top %d cholesterol-biosynthesis leading-edge genes (ordered roughly along mevalonate pathway)",
                          length(ordered_top))) +
  theme_cowplot(font_size = 11) +
  theme(plot.subtitle = element_text(size = 12, face = "bold"),
        axis.text.x   = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y   = element_text(size = 10,
                                     face = ifelse(rev(ordered_top) %in% gwas_set,
                                                   "bold.italic", "plain")),
        panel.grid    = element_blank(),
        axis.line     = element_blank(), axis.ticks = element_blank())

# =========================================================================
# Panel D — SREBF2 spotlight + downstream targets
# =========================================================================
# SREBF2 itself + its canonical targets
srebf2_targets <- c("SREBF2","HMGCR","HMGCS1","FDFT1","LSS","SQLE","MVD",
                    "DHCR7","DHCR24","LDLR","INSIG1")
focal_df <- df |> filter(genes %in% srebf2_targets) |>
  mutate(cell_type = order_ct(cell_type),
         class     = class_of(as.character(cell_type)),
         genes     = factor(genes, levels = srebf2_targets),
         ci_lo     = estimate - 1.96*se,
         ci_hi     = estimate + 1.96*se,
         sig       = padj < 0.1)

p_D <- ggplot(focal_df, aes(x = estimate, y = cell_type, colour = class)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, yend = cell_type), linewidth = 0.4) +
  geom_point(aes(shape = sig), size = 1.8) +
  facet_wrap(~ genes, ncol = 4, scales = "free_x") +
  scale_colour_manual(values = CLASS_COL, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(x = "Effect estimate", y = NULL,
       subtitle = "D. SREBF2 (PGC3 GWAS gene) + canonical downstream targets: coordinated suppression in same cell types") +
  theme_cowplot(font_size = 11) +
  theme(plot.subtitle = element_text(size = 12, face = "bold"),
        axis.text.y   = element_text(size = 8),
        strip.text    = element_text(face = "bold.italic", size = 11),
        strip.background = element_blank())

# =========================================================================
# Panel E — Cross-pathway correlation: cholesterol vs myelin / phospholipid
# =========================================================================
xy <- prior_long |>
  filter(pathway %in% c("REACTOME_CHOLESTEROL_BIOSYNTHESIS",
                        "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                        "GOBP_MYELINATION",
                        "REACTOME_PHOSPHOLIPID_METABOLISM",
                        "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING")) |>
  select(cell_type, pathway, NES) |>
  pivot_wider(names_from = pathway, values_from = NES) |>
  mutate(class = class_of(cell_type))

cor1 <- cor.test(xy$REACTOME_CHOLESTEROL_BIOSYNTHESIS,
                 xy$HALLMARK_OXIDATIVE_PHOSPHORYLATION, method = "spearman")
cor2 <- cor.test(xy$REACTOME_CHOLESTEROL_BIOSYNTHESIS,
                 xy$GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING, method = "spearman")
cat(sprintf("\nSpearman rho(cholesterol vs oxphos):    %.2f (p=%.3g)\n",
            cor1$estimate, cor1$p.value))
cat(sprintf("Spearman rho(cholesterol vs synaptic):  %.2f (p=%.3g)\n",
            cor2$estimate, cor2$p.value))

p_E <- ggplot(xy, aes(x = REACTOME_CHOLESTEROL_BIOSYNTHESIS,
                       y = HALLMARK_OXIDATIVE_PHOSPHORYLATION,
                       colour = class)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(aes(label = cell_type), size = 3.2,
                  max.overlaps = Inf, seed = 1,
                  box.padding = 0.3, point.padding = 0.2,
                  segment.colour = "grey60", show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Spearman rho = %.2f  (p = %.2g)",
                           cor1$estimate, cor1$p.value),
           hjust = -0.05, vjust = 1.5, size = 4, colour = "black") +
  scale_colour_manual(values = CLASS_COL, name = "Class") +
  labs(x = "NES (cholesterol biosynthesis)",
       y = "NES (HALLMARK OxPhos)",
       subtitle = "E. Is cholesterol collapse coupled to OxPhos collapse across cell types?") +
  theme_cowplot(font_size = 12) +
  theme(plot.subtitle = element_text(size = 12, face = "bold"))

# Compose figure
row1 <- plot_grid(p_A, p_E, nrow = 1, rel_widths = c(1.3, 1))
row2 <- p_B
row3 <- plot_grid(p_C, p_D, nrow = 1, rel_widths = c(1.1, 1.4))

full <- plot_grid(row1, row2, row3, ncol = 1,
                  rel_heights = c(1, 1.05, 1.55))
ggsave(file.path(OUT_DIR, "figure_cholesterol_glia_exc.png"),
       full, width = 18, height = 21, dpi = 200, bg = "white",
       limitsize = FALSE)
ggsave(file.path(OUT_DIR, "figure_cholesterol_glia_exc.pdf"),
       full, width = 18, height = 21, bg = "white",
       limitsize = FALSE)
cat(sprintf("\nWrote %s/figure_cholesterol_glia_exc.{png,pdf}\n", OUT_DIR))
