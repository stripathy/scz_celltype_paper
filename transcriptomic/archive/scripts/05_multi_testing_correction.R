#!/usr/bin/env Rscript
# Cross-cell-type multiple-testing correction for GSEA results.
# Three corrections:
#   (A) Global BH across all (cell type x pathway) p-values
#   (B) Per-pathway Fisher's combined p across cell types, then BH
#   (C) Per-pathway Stouffer's signed Z combination (preserves direction)
# Then re-make the curated heatmap with all three correction layers shown.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(cowplot); library(scales)
})

OUT_DIR <- "results"

EXC <- c("L2_3 IT","L4 IT","L5 IT","L5 ET","L5_6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")
INH <- c("Lamp5","Pax6","Pvalb","Sncg","Sst","Sst Chodl","Vip","Chandelier")
GLI <- c("Astro","Oligo","OPC","Micro-PVM","Endo","VLMC")
order_ct <- function(cts) {
  ord <- c(intersect(EXC, cts), intersect(INH, cts), intersect(GLI, cts))
  factor(cts, levels = ord)
}

# ---- Load GSEA results -----------------------------------------------------
g <- read_csv(file.path(OUT_DIR, "gsea_all_celltypes.csv"), show_col_types = FALSE)
cat(sprintf("Loaded %d GSEA results (%d pathways x %d cell types)\n",
            nrow(g), n_distinct(g$pathway), n_distinct(g$cell_type)))

# ---- (A) Global BH across all (pathway x cell type) -----------------------
g <- g |>
  mutate(padj_global = p.adjust(pval, method = "BH"))

cat("\n=== (A) Global BH correction ===\n")
cat(sprintf("Tests with padj_global < 0.05: %d (was %d at per-cell-type FDR)\n",
            sum(g$padj_global < 0.05, na.rm = TRUE),
            sum(g$padj < 0.05, na.rm = TRUE)))
cat(sprintf("Tests with padj_global < 0.01: %d\n",
            sum(g$padj_global < 0.01, na.rm = TRUE)))

# Distribution of "which cell types still significant" per pathway
global_sig_by_pw <- g |> filter(padj_global < 0.05) |>
  count(pathway, sort = TRUE)
cat(sprintf("Pathways with >=1 global-BH hit: %d\n", nrow(global_sig_by_pw)))

# ---- (B) Fisher's combined p per pathway across cell types -----------------
fisher_combine <- function(p) {
  p <- p[!is.na(p) & p > 0]
  if (length(p) < 2) return(NA_real_)
  k <- length(p)
  chi <- -2 * sum(log(p))
  pchisq(chi, df = 2 * k, lower.tail = FALSE)
}

pw_fisher <- g |>
  group_by(pathway, source) |>
  summarise(k_cells       = n(),
            mean_NES      = mean(NES, na.rm = TRUE),
            sd_NES        = sd(NES, na.rm = TRUE),
            p_fisher      = fisher_combine(pval),
            n_pos_NES     = sum(NES > 0, na.rm = TRUE),
            n_neg_NES     = sum(NES < 0, na.rm = TRUE),
            n_sig_per_ct  = sum(padj < 0.05, na.rm = TRUE),
            .groups = "drop") |>
  mutate(padj_fisher = p.adjust(p_fisher, method = "BH"))

cat("\n=== (B) Per-pathway Fisher's combined p, BH-corrected ===\n")
cat(sprintf("Pathways with padj_fisher < 0.05: %d / %d\n",
            sum(pw_fisher$padj_fisher < 0.05, na.rm = TRUE), nrow(pw_fisher)))
cat(sprintf("Pathways with padj_fisher < 0.01: %d\n",
            sum(pw_fisher$padj_fisher < 0.01, na.rm = TRUE)))

# ---- (C) Stouffer's signed-Z combination per pathway -----------------------
# Convert per-test p to a signed Z using sign(NES); combine with equal weights;
# this preserves direction (collapses bidirectional signals).
stouffer_signed <- function(p, nes) {
  ok <- !is.na(p) & p > 0 & p < 1 & !is.na(nes)
  p <- p[ok]; nes <- nes[ok]
  if (length(p) < 2) return(c(z = NA_real_, p = NA_real_))
  z_i <- sign(nes) * qnorm(p / 2, lower.tail = FALSE)
  z   <- sum(z_i) / sqrt(length(z_i))
  p_two <- 2 * pnorm(-abs(z))
  c(z = z, p = p_two)
}

pw_stouffer <- g |>
  group_by(pathway, source) |>
  summarise(stouffer = list(stouffer_signed(pval, NES)),
            mean_NES = mean(NES, na.rm = TRUE),
            n_cells  = n(),
            .groups = "drop") |>
  mutate(stouffer_Z = sapply(stouffer, `[[`, "z"),
         p_stouffer = sapply(stouffer, `[[`, "p"),
         padj_stouffer = p.adjust(p_stouffer, method = "BH")) |>
  select(-stouffer)

cat("\n=== (C) Per-pathway Stouffer's signed-Z, BH-corrected ===\n")
cat(sprintf("Pathways with padj_stouffer < 0.05: %d / %d\n",
            sum(pw_stouffer$padj_stouffer < 0.05, na.rm = TRUE),
            nrow(pw_stouffer)))
cat(sprintf("Of those, |Z| > 3:  %d\n",
            sum(pw_stouffer$padj_stouffer < 0.05 &
                abs(pw_stouffer$stouffer_Z) > 3, na.rm = TRUE)))

# Save combined pathway-level summary
pw_summary <- pw_fisher |>
  left_join(pw_stouffer |> select(pathway, stouffer_Z, p_stouffer, padj_stouffer),
            by = "pathway")
write_csv(pw_summary, file.path(OUT_DIR, "pathway_multitest_summary.csv"))
cat(sprintf("\nWrote %s\n",
            file.path(OUT_DIR, "pathway_multitest_summary.csv")))

# ---- Top pathways by each correction --------------------------------------
cat("\n=== Top 20 pathways by Stouffer signed-Z (most directionally coherent) ===\n")
print(pw_summary |> arrange(desc(abs(stouffer_Z))) |> head(20) |>
        select(pathway, source, k_cells, mean_NES, stouffer_Z,
               padj_stouffer, padj_fisher, n_sig_per_ct), n = 20)

cat("\n=== Top 20 pathways by Fisher's combined p ===\n")
print(pw_summary |> arrange(p_fisher) |> head(20) |>
        select(pathway, source, k_cells, mean_NES, p_fisher,
               padj_fisher, padj_stouffer, n_sig_per_ct), n = 20)

# ---- Re-make curated heatmap with multi-correction annotation --------------
curated <- tribble(
  ~theme,                          ~display,                                 ~pathway,
  "Synaptic / vesicle (DOWN)",     "trans-synaptic signaling (reg)",         "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING",
  "Synaptic / vesicle (DOWN)",     "presynapse",                             "GOCC_PRESYNAPSE",
  "Synaptic / vesicle (DOWN)",     "postsynaptic membrane",                  "GOCC_POSTSYNAPTIC_MEMBRANE",
  "Synaptic / vesicle (DOWN)",     "synaptic vesicle cycle",                 "GOBP_SYNAPTIC_VESICLE_CYCLE",
  "Synaptic / vesicle (DOWN)",     "neurotransmitter secretion",             "GOBP_NEUROTRANSMITTER_SECRETION",
  "Synaptic / vesicle (DOWN)",     "transmission across chemical synapses",  "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES",
  "Ion channels / excitability",   "regulation of membrane potential",       "GOBP_REGULATION_OF_MEMBRANE_POTENTIAL",
  "Ion channels / excitability",   "gated channel activity",                 "GOMF_GATED_CHANNEL_ACTIVITY",
  "Ion channels / excitability",   "monoatomic cation channel activity",     "GOMF_MONOATOMIC_CATION_CHANNEL_ACTIVITY",
  "Mitochondrial / OxPhos",        "HALLMARK oxidative phosphorylation",     "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "Mitochondrial / OxPhos",        "respiratory electron transport",         "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "Mitochondrial / OxPhos",        "aerobic respiration",                    "GOBP_AEROBIC_RESPIRATION",
  "Cholesterol / lipid",           "HALLMARK cholesterol homeostasis",       "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "Cholesterol / lipid",           "cholesterol biosynthesis (Reactome)",    "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
  "Cholesterol / lipid",           "sterol biosynthetic process",            "GOBP_STEROL_BIOSYNTHETIC_PROCESS",
  "Axon guidance",                 "signaling by ROBO receptors",            "REACTOME_SIGNALING_BY_ROBO_RECEPTORS",
  "Axon guidance",                 "SLIT/ROBO regulation",                   "REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS",
  "Translation / ribosome",        "REACTOME translation",                   "REACTOME_TRANSLATION",
  "Translation / ribosome",        "cytoplasmic translation",                "GOBP_CYTOPLASMIC_TRANSLATION",
  "Translation / ribosome",        "ribosomal subunit",                      "GOCC_RIBOSOMAL_SUBUNIT",
  "Growth / stress (UP)",          "HALLMARK MYC targets V1",                "HALLMARK_MYC_TARGETS_V1",
  "Growth / stress (UP)",          "HALLMARK mTORC1 signaling",              "HALLMARK_MTORC1_SIGNALING",
  "Growth / stress (UP)",          "HSF1 activation (chaperones)",           "REACTOME_HSF1_ACTIVATION",
  "Growth / stress (UP)",          "chaperone-mediated protein folding",     "GOBP_CHAPERONE_MEDIATED_PROTEIN_FOLDING",
  "BMP / TGFb signaling (UP)",     "signaling by BMP",                       "REACTOME_SIGNALING_BY_BMP",
  "BMP / TGFb signaling (UP)",     "I-SMAD binding",                         "GOMF_I_SMAD_BINDING",
  "BMP / TGFb signaling (UP)",     "HALLMARK TGF-beta signaling",            "HALLMARK_TGF_BETA_SIGNALING"
)

# Annotate curated pathways with pathway-level corrections
curated_ann <- curated |>
  left_join(pw_summary |>
              select(pathway, padj_fisher, stouffer_Z, padj_stouffer, n_sig_per_ct),
            by = "pathway")

cat("\n=== Curated pathways: pathway-level correction summary ===\n")
print(curated_ann |> select(theme, display, padj_fisher, stouffer_Z,
                            padj_stouffer, n_sig_per_ct), n = 30)
write_csv(curated_ann, file.path(OUT_DIR, "curated_pathway_corrections.csv"))

# ---- Heatmap with both per-cell + global-BH significance markers ----------
d <- g |> inner_join(curated, by = "pathway") |>
  select(cell_type, theme, display, pathway, NES, padj, padj_global)
d$cell_type <- order_ct(d$cell_type)
d$display   <- factor(d$display,  levels = rev(unique(curated$display)))
d$theme     <- factor(d$theme,    levels = unique(curated$theme))
# Marker text reflects per-cell-type FDR + an overlay for global-BH survival
d <- d |> mutate(
  marker_local = case_when(padj < 0.001 ~ "***",
                            padj < 0.01  ~ "**",
                            padj < 0.05  ~ "*",
                            padj < 0.1   ~ "+", TRUE ~ ""),
  marker_global = if_else(padj_global < 0.05, "●", "")
)

# Pathway-level row strip labels showing Fisher BH
strip_labels <- curated_ann |>
  mutate(label = case_when(
    padj_fisher < 0.001 ~ sprintf("%s ***", display),
    padj_fisher < 0.01  ~ sprintf("%s **",  display),
    padj_fisher < 0.05  ~ sprintf("%s *",   display),
    padj_fisher < 0.1   ~ sprintf("%s +",   display),
    TRUE                ~ as.character(display)
  ))
d$display_label <- factor(strip_labels$label[match(d$display, strip_labels$display)],
                          levels = rev(strip_labels$label))

p_heat <- ggplot(d, aes(x = cell_type, y = display_label, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = marker_local), size = 3.4, vjust = 0.65) +
  geom_text(aes(label = marker_global), size = 2.8, vjust = -1.4, hjust = -0.95,
            colour = "black", fontface = "bold") +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y",
             switch = "y", labeller = labeller(theme = label_wrap_gen(18))) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, name = "NES",
                       limits = c(-3, 3), oob = scales::squish,
                       breaks = c(-3,-2,-1,0,1,2,3)) +
  labs(x = NULL, y = NULL,
       title = "Pan-cell-type GSEA with multi-testing correction layers",
       subtitle = paste("Cells: per-cell-type FDR (***<0.001, **<0.01, *<0.05, +<0.1)",
                        " - Row labels: pathway-level Fisher BH",
                        " - Black dot: global BH across all 23 x 12K tests < 0.05",
                        sep = "\n"),
       caption = "Pathway-level (Fisher BH) tests 'is this pathway broadly affected across cell types'") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x       = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y       = element_text(size = 10),
        plot.title        = element_text(size = 14, face = "bold"),
        plot.subtitle     = element_text(size = 10, colour = "grey20"),
        plot.caption      = element_text(size = 9, colour = "grey40"),
        panel.grid        = element_blank(),
        axis.line         = element_blank(),
        axis.ticks        = element_blank(),
        strip.placement   = "outside",
        strip.background  = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1,
                                         face = "bold", size = 11),
        panel.spacing.y   = unit(0.25, "lines"))

ggsave(file.path(OUT_DIR, "GSEA_themes_multitest.png"),
       p_heat, width = 15, height = 12, dpi = 220, bg = "white")
ggsave(file.path(OUT_DIR, "GSEA_themes_multitest.pdf"),
       p_heat, width = 15, height = 12, bg = "white")
cat("\nWrote GSEA_themes_multitest.{png,pdf}\n")

# ---- Summary table of curated theme survival ------------------------------
cat("\n=== Survival of curated themes under stringent multi-testing ===\n")
summary_tbl <- curated_ann |>
  mutate(survives_fisher_05   = padj_fisher   < 0.05,
         survives_stouffer_05 = padj_stouffer < 0.05) |>
  group_by(theme) |>
  summarise(n_pathways   = n(),
            n_fisher_sig = sum(survives_fisher_05, na.rm = TRUE),
            n_stouffer_sig = sum(survives_stouffer_05, na.rm = TRUE),
            mean_Z       = mean(stouffer_Z, na.rm = TRUE),
            .groups = "drop")
print(summary_tbl, n = 30)
write_csv(summary_tbl, file.path(OUT_DIR, "theme_survival_summary.csv"))

cat("\nDone.\n")
