#!/usr/bin/env Rscript
# Mechanistic hypothesis-testing for OxPhos suppression in inhibitory neurons.
# Tests: PGC-1a axis, glycolytic compensation, activity-dependent change,
# antioxidant response, mitochondrial dynamics / quality control.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(cowplot); library(fgsea); library(msigdbr); library(ggrepel)
})

OUT_DIR <- "results/exploratory/story_oxphos_mech"
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

# ---- Build gene sets for each mechanism ------------------------------------
# 1. PGC-1a transcriptional axis (mito biogenesis master regulators + co-activators)
pgc1_axis <- c(
  "PPARGC1A","PPARGC1B","PPRC1",       # PGC-1 family
  "NRF1","GABPA","GABPB1","GABPB2",    # NRF1/GABP
  "NFE2L2","NFE2L1",                   # NRF2/NRF3
  "TFAM","TFB1M","TFB2M","POLRMT",     # mtDNA TFs
  "ESRRA","ESRRB","ESRRG",             # ERR receptors
  "YY1","MEF2A","MEF2C","MEF2D",       # accessory TFs
  "PPARA","PPARD","PPARG"              # PPARs
)

# 2. Glycolysis (HALLMARK glycolysis + key enzymes)
collect <- function(cat, subcat = NULL) {
  args <- list(species = "Homo sapiens", collection = cat)
  if (!is.null(subcat)) args$subcollection <- subcat
  do.call(msigdbr, args) |> transmute(gs_name, gene = gene_symbol)
}
H <- collect("H")
glycolysis <- H |> filter(gs_name == "HALLMARK_GLYCOLYSIS") |> pull(gene) |> unique()

# 3. Antioxidant / ROS response (HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY + NRF2 targets)
ros <- H |> filter(gs_name == "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY") |>
  pull(gene) |> unique()
nrf2_targets <- c(
  "HMOX1","NQO1","GCLC","GCLM","TXN","TXNRD1","GSTP1","GSTM1","GSTM2",
  "GSR","PRDX1","PRDX6","FTH1","FTL","SQSTM1","MAFG","KEAP1","SOD1","SOD2","CAT","GPX1","GPX4"
)
antioxidant <- unique(c(ros, nrf2_targets))

# 4. Immediate early genes (activity-dependent transcription)
ieg <- c("FOS","FOSB","JUN","JUNB","JUND","EGR1","EGR2","EGR3","EGR4",
         "ARC","NR4A1","NR4A2","NR4A3","NPAS4","BDNF","HOMER1","DUSP1",
         "IER2","IER3","ATF3","BTG2","GADD45B","CYR61")

# 5. Mitochondrial dynamics (fission/fusion)
mito_dyn <- c("DNM1L","MFF","FIS1","MIEF1","MIEF2","MTFP1",     # fission
              "MFN1","MFN2","OPA1","INF2",                       # fusion
              "STOML2","SLC25A46","CHCHD3","CHCHD6")             # cristae

# 6. Mitophagy / mito quality control
mitophagy <- c("PINK1","PRKN","BNIP3","BNIP3L","FUNDC1","NIPSNAP1","NIPSNAP2",
               "CALCOCO2","OPTN","TBK1","TAX1BP1","SQSTM1","NBR1",
               "MAP1LC3A","MAP1LC3B","MAP1LC3C","GABARAPL1","GABARAPL2","ATG7","ATG12",
               "LONP1","CLPP","CLPX","YME1L1","AFG3L2","SPG7","HSPD1","HSPE1","HSPA9")

mech_sets <- list(
  "PGC-1a / mito biogenesis"    = pgc1_axis,
  "Glycolysis (HALLMARK)"       = glycolysis,
  "Antioxidant / ROS / NRF2"    = antioxidant,
  "Activity-dependent (IEGs)"   = ieg,
  "Mito dynamics (fission/fusion)" = mito_dyn,
  "Mito quality control / mitophagy" = mitophagy
)
cat("Mechanism gene set sizes:\n")
print(map_int(mech_sets, length))

# ---- Run GSEA per cell type for these custom sets --------------------------
cts <- sort(unique(df$cell_type))
mech_results <- list()
set.seed(99)
for (ct in cts) {
  sub <- df |> filter(cell_type == ct) |>
    arrange(desc(abs(z))) |> distinct(genes, .keep_all = TRUE)
  if (nrow(sub) < 500) next
  ranks <- setNames(sub$z, sub$genes)
  ranks <- sort(ranks, decreasing = TRUE)
  res <- fgsea(pathways = mech_sets, stats = ranks,
               minSize = 5, maxSize = 500, nPermSimple = 30000)
  res$cell_type <- ct
  mech_results[[ct]] <- as_tibble(res)
}
mech_long <- bind_rows(mech_results) |>
  mutate(pathway = factor(pathway, levels = names(mech_sets)))
write_csv(mech_long |> select(-leadingEdge),
          file.path(OUT_DIR, "mechanism_gsea.csv"))

# =========================================================================
# Panel A — PGC-1a axis gene-level forest plot
# =========================================================================
# Focus on inhibitory cell types (plus Astro as comparison)
focus_ct <- c("Sst","Pvalb","Vip","Lamp5","Sncg","Chandelier","Astro")

# Filter to TFs that have data in at least 3 of these cell types
tf_data <- df |> filter(genes %in% pgc1_axis, cell_type %in% focus_ct) |>
  group_by(genes) |> filter(n() >= 3) |> ungroup()
present_tfs <- intersect(pgc1_axis, unique(tf_data$genes))
cat(sprintf("\nPGC-1a axis TFs present in >=3 of %d inhibitory cell types: %d / %d\n",
            length(focus_ct), length(present_tfs), length(pgc1_axis)))
print(present_tfs)

tf_df <- df |> filter(genes %in% present_tfs, cell_type %in% focus_ct) |>
  mutate(cell_type = factor(cell_type, levels = focus_ct),
         genes = factor(genes, levels = present_tfs),
         ci_lo = estimate - 1.96*se,
         ci_hi = estimate + 1.96*se,
         sig   = padj < 0.1,
         class = class_of(as.character(cell_type)))

p_A <- ggplot(tf_df, aes(x = estimate, y = cell_type, colour = class)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, yend = cell_type), linewidth = 0.4) +
  geom_point(aes(shape = sig), size = 2) +
  facet_wrap(~ genes, ncol = 5, scales = "free_x") +
  scale_colour_manual(values = CLASS_COL, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(x = "Effect estimate (SCZ vs control)", y = NULL,
       subtitle = "A. Mitochondrial biogenesis TFs (PGC-1a axis) — filled = padj<0.1") +
  theme_cowplot(font_size = 10) +
  theme(plot.subtitle = element_text(size = 12, face = "bold"),
        axis.text.y   = element_text(size = 9),
        strip.text    = element_text(face = "bold.italic", size = 10),
        strip.background = element_blank())

# =========================================================================
# Panel B — Mechanism heatmap (GSEA NES per gene set per cell type)
# =========================================================================
mech_long$cell_type <- order_ct(mech_long$cell_type)

p_B <- ggplot(mech_long, aes(x = cell_type, y = pathway, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(padj < 0.001, "***",
                        ifelse(padj < 0.01,  "**",
                        ifelse(padj < 0.05,  "*",
                        ifelse(padj < 0.1,   "+", ""))))),
            size = 3.6, vjust = 0.7) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, name = "NES",
                       limits = c(-3, 3), oob = scales::squish) +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL,
       subtitle = "B. Mechanism gene-set GSEA across all cell types — which hypothesis is supported where?") +
  theme_cowplot(font_size = 12) +
  theme(plot.subtitle = element_text(size = 12, face = "bold"),
        axis.text.x   = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y   = element_text(size = 10),
        panel.grid    = element_blank(),
        axis.line     = element_blank(), axis.ticks = element_blank())

# =========================================================================
# Panel C — Compensation test: OxPhos NES vs Glycolysis NES
# =========================================================================
xy <- mech_long |> filter(pathway == "Glycolysis (HALLMARK)") |>
  select(cell_type, glyco_NES = NES) |>
  left_join(prior_long |> filter(pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") |>
              select(cell_type, oxphos_NES = NES),
            by = "cell_type") |>
  mutate(class = class_of(as.character(cell_type)))

cor_test_C <- cor.test(xy$oxphos_NES, xy$glyco_NES, method = "spearman")
cat(sprintf("\nSpearman rho(OxPhos NES vs Glycolysis NES): %.2f (p=%.3g)\n",
            cor_test_C$estimate, cor_test_C$p.value))

p_C <- ggplot(xy, aes(x = oxphos_NES, y = glyco_NES, colour = class)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_abline(slope = -1, intercept = 0, linetype = "dotted",
              colour = "grey50", linewidth = 0.3) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(aes(label = cell_type), size = 3.2,
                  max.overlaps = Inf, seed = 1, box.padding = 0.3,
                  segment.colour = "grey60", show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Spearman rho = %.2f  (p = %.2g)\nDotted line = perfect Warburg compensation",
                           cor_test_C$estimate, cor_test_C$p.value),
           hjust = -0.05, vjust = 1.3, size = 3.8, colour = "black") +
  scale_colour_manual(values = CLASS_COL, name = "Class") +
  labs(x = "NES (OxPhos)",
       y = "NES (Glycolysis)",
       subtitle = "C. Glycolytic compensation? If OxPhos down forces a Warburg switch, we'd expect negative slope.") +
  theme_cowplot(font_size = 12) +
  theme(plot.subtitle = element_text(size = 12, face = "bold"))

# =========================================================================
# Panel D — Activity hypothesis: IEG NES vs OxPhos NES
# =========================================================================
xy2 <- mech_long |> filter(pathway == "Activity-dependent (IEGs)") |>
  select(cell_type, ieg_NES = NES, ieg_padj = padj) |>
  left_join(prior_long |> filter(pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") |>
              select(cell_type, oxphos_NES = NES, oxphos_padj = padj),
            by = "cell_type") |>
  mutate(class = class_of(as.character(cell_type)))

cor_test_D <- cor.test(xy2$oxphos_NES, xy2$ieg_NES, method = "spearman")
cat(sprintf("Spearman rho(OxPhos NES vs IEG NES): %.2f (p=%.3g)\n",
            cor_test_D$estimate, cor_test_D$p.value))

p_D <- ggplot(xy2, aes(x = oxphos_NES, y = ieg_NES, colour = class)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              colour = "grey50", linewidth = 0.3) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(aes(label = cell_type), size = 3.2,
                  max.overlaps = Inf, seed = 2, box.padding = 0.3,
                  segment.colour = "grey60", show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Spearman rho = %.2f  (p = %.2g)\nDotted line = activity-driven OxPhos hypothesis",
                           cor_test_D$estimate, cor_test_D$p.value),
           hjust = -0.05, vjust = 1.3, size = 3.8, colour = "black") +
  scale_colour_manual(values = CLASS_COL, guide = "none") +
  labs(x = "NES (OxPhos)",
       y = "NES (Activity-dependent IEGs)",
       subtitle = "D. Activity hypothesis: if firing demand drives OxPhos, IEGs should track OxPhos.") +
  theme_cowplot(font_size = 12) +
  theme(plot.subtitle = element_text(size = 12, face = "bold"))

# =========================================================================
# Panel E — Within-inhibitory ranking + key IEG forest plot
# =========================================================================
key_iegs <- c("FOS","ARC","BDNF","NPAS4","EGR1","NR4A1","HOMER1")
ieg_data <- df |> filter(genes %in% key_iegs, cell_type %in% INH) |>
  mutate(cell_type = factor(cell_type, levels = INH),
         genes = factor(genes, levels = key_iegs),
         ci_lo = estimate - 1.96*se,
         ci_hi = estimate + 1.96*se,
         sig   = padj < 0.1)

p_E <- ggplot(ieg_data, aes(x = estimate, y = cell_type)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, yend = cell_type),
               colour = "#882255", linewidth = 0.4) +
  geom_point(aes(shape = sig), colour = "#882255", size = 2) +
  facet_wrap(~ genes, ncol = 4, scales = "free_x") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(x = "Effect estimate (SCZ vs control)", y = NULL,
       subtitle = "E. Activity-dependent IEGs across inhibitory cells (filled = padj<0.1)") +
  theme_cowplot(font_size = 10) +
  theme(plot.subtitle = element_text(size = 12, face = "bold"),
        axis.text.y   = element_text(size = 9),
        strip.text    = element_text(face = "bold.italic", size = 11),
        strip.background = element_blank())

# =========================================================================
# Compose figure
# =========================================================================
row1 <- p_A
row2 <- p_B
row3 <- plot_grid(p_C, p_D, nrow = 1, rel_widths = c(1, 1))
row4 <- p_E

full <- plot_grid(row1, row2, row3, row4, ncol = 1,
                  rel_heights = c(1.2, 1.0, 1.0, 1.1))
ggsave(file.path(OUT_DIR, "figure_oxphos_mechanism_tests.png"),
       full, width = 16, height = 22, dpi = 200, bg = "white",
       limitsize = FALSE)
ggsave(file.path(OUT_DIR, "figure_oxphos_mechanism_tests.pdf"),
       full, width = 16, height = 22, bg = "white",
       limitsize = FALSE)
cat(sprintf("\nWrote %s/figure_oxphos_mechanism_tests.{png,pdf}\n", OUT_DIR))

# =========================================================================
# Brief text summary of which hypothesis wins per cell type
# =========================================================================
cat("\n=== Mechanism-by-mechanism NES across inhibitory cell types ===\n")
inh_summary <- mech_long |> filter(cell_type %in% INH) |>
  select(cell_type, pathway, NES, padj) |>
  pivot_wider(names_from = pathway, values_from = c(NES, padj))
print(mech_long |> filter(cell_type %in% INH) |>
        mutate(sig = case_when(padj < 0.001 ~ "***",
                                padj < 0.01  ~ "**",
                                padj < 0.05  ~ "*",
                                padj < 0.1   ~ "+",  TRUE ~ "")) |>
        select(cell_type, pathway, NES, padj, sig) |>
        arrange(cell_type, padj), n = 50)
