#!/usr/bin/env Rscript
# Clean readable summary figure:
#   Panel A: Sst ORA vs GSEA for the same set of pathways (the rescue)
#   Panel B: Cross-cell-type GSEA NES heatmap, curated pathway set

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(scales)
  library(ggplot2); library(cowplot)
})

OUT_DIR <- "results"

EXC <- c("L2_3 IT","L4 IT","L5 IT","L5 ET","L5_6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")
INH <- c("Lamp5","Pax6","Pvalb","Sncg","Sst","Sst Chodl","Vip","Chandelier")
GLI <- c("Astro","Oligo","OPC","Micro-PVM","Endo","VLMC")
order_ct <- function(cts) {
  ord <- c(intersect(EXC, cts), intersect(INH, cts), intersect(GLI, cts))
  factor(cts, levels = ord)
}

g <- read_csv(file.path(OUT_DIR, "gsea_all_celltypes.csv"), show_col_types = FALSE)

# ----- Curated set of biologically interpretable pathways -------------------
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

d <- g |> inner_join(curated, by = "pathway") |>
  select(cell_type, theme, display, pathway, NES, padj)
d$cell_type <- order_ct(d$cell_type)
# Lock row order to the curated table
d$display   <- factor(d$display,
                      levels = rev(unique(curated$display)))
d$theme     <- factor(d$theme,
                      levels = unique(curated$theme))

p_heat <- ggplot(d, aes(x = cell_type, y = display, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(padj < 0.001, "***",
                       ifelse(padj < 0.01,  "**",
                       ifelse(padj < 0.05,  "*",
                       ifelse(padj < 0.1,   "+",  ""))))),
            size = 3.6, vjust = 0.65) +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y",
             switch = "y", labeller = labeller(theme = label_wrap_gen(18))) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, name = "NES",
                       limits = c(-3, 3), oob = scales::squish,
                       breaks = c(-3,-2,-1,0,1,2,3)) +
  labs(x = NULL, y = NULL,
       title = "Pan-cell-type GSEA NES of curated pathway themes in schizophrenia",
       caption = "***padj<0.001  **padj<0.01  *padj<0.05  +padj<0.1   (orange=up, blue=down)") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x       = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y       = element_text(size = 10),
        plot.title        = element_text(size = 14, face = "bold"),
        plot.caption      = element_text(size = 10, colour = "grey30"),
        panel.grid        = element_blank(),
        axis.line         = element_blank(),
        axis.ticks        = element_blank(),
        strip.placement   = "outside",
        strip.background  = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1,
                                         face = "bold", size = 11),
        panel.spacing.y   = unit(0.25, "lines"))

ggsave(file.path(OUT_DIR, "GSEA_themes_heatmap.png"),
       p_heat, width = 14, height = 11, dpi = 220, bg = "white")
ggsave(file.path(OUT_DIR, "GSEA_themes_heatmap.pdf"),
       p_heat, width = 14, height = 11, bg = "white")
cat("Wrote GSEA_themes_heatmap.{png,pdf}\n")

# ----- Sst ORA vs GSEA bar chart ---------------------------------------------
sst_g <- g |> filter(cell_type == "Sst") |>
  select(pathway, NES, padj_gsea = padj)

# Read ORA results from previous run
ora <- read_csv("results/exploratory/archive_ora_sst/enrich_fdr_full.csv", show_col_types = FALSE)
ora_sst_down <- ora |> filter(query_set == "down") |>
  transmute(term_name, padj_ora = p_value)
ora_sst_up   <- ora |> filter(query_set == "up") |>
  transmute(term_name, padj_ora = p_value)

# Map MSigDB pathway names to human-readable terms used in ORA
map_ora <- tribble(
  ~pathway,                                            ~ora_term,
  "GOCC_PRESYNAPSE",                                   "presynapse",
  "GOCC_SYNAPSE",                                      "synapse",
  "GOCC_SYNAPTIC_VESICLE",                             "synaptic vesicle",
  "GOCC_EXOCYTIC_VESICLE",                             "exocytic vesicle",
  "GOCC_SECRETORY_VESICLE",                            "secretory vesicle",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",                NA,
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",           NA,
  "GOBP_AEROBIC_RESPIRATION",                          NA,
  "GOMF_I_SMAD_BINDING",                               "I-SMAD binding",
  "REACTOME_SIGNALING_BY_BMP",                         NA,
  "REACTOME_CHOLESTEROL_BIOSYNTHESIS",                 NA,
  "GOMF_STRUCTURAL_CONSTITUENT_OF_RIBOSOME",           NA,
  "REACTOME_TRANSLATION",                              NA,
  "REACTOME_SIGNALING_BY_ROBO_RECEPTORS",              NA,
  "REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS", NA
)
compare <- map_ora |>
  left_join(sst_g, by = "pathway") |>
  left_join(bind_rows(ora_sst_down, ora_sst_up) |> distinct(term_name, .keep_all = TRUE),
            by = c("ora_term" = "term_name")) |>
  mutate(direction = ifelse(NES > 0, "up", "down"),
         display = str_replace_all(pathway, c("HALLMARK_" = "Hallmark: ",
                                              "REACTOME_" = "Reactome: ",
                                              "GOCC_"     = "GO:CC: ",
                                              "GOBP_"     = "GO:BP: ",
                                              "GOMF_"     = "GO:MF: ",
                                              "_"         = " ")) |> str_to_lower(),
         display = str_replace(display, "^(hallmark|reactome|go:cc|go:bp|go:mf): (.*)",
                               "\\1: \\2"))

cmp_long <- compare |>
  transmute(display, GSEA = -log10(padj_gsea), ORA = -log10(padj_ora),
            direction) |>
  pivot_longer(c(GSEA, ORA), names_to = "method", values_to = "neglog10p") |>
  mutate(neglog10p = replace_na(neglog10p, 0))

p_cmp <- ggplot(cmp_long,
                aes(x = neglog10p, y = reorder(display, neglog10p),
                    fill = method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed",
             colour = "grey30", linewidth = 0.4) +
  scale_fill_manual(values = c(GSEA = "#117733", ORA = "#888888"),
                    name = "Method") +
  labs(x = expression(-log[10]~"(adjusted p)"), y = NULL,
       title = "Sst — same pathways, ORA vs GSEA",
       subtitle = "GSEA recovers signals that ORA missed",
       caption = "0 = pathway not returned by that method (above gprofiler's enrichment requirements)") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.y    = element_text(size = 10),
        plot.title     = element_text(size = 13, face = "bold"),
        plot.subtitle  = element_text(size = 11, colour = "grey20"),
        plot.caption   = element_text(size = 9, colour = "grey40"),
        legend.position = c(0.85, 0.15))

ggsave(file.path(OUT_DIR, "Sst_ORA_vs_GSEA.png"),
       p_cmp, width = 11, height = 7, dpi = 220, bg = "white")
ggsave(file.path(OUT_DIR, "Sst_ORA_vs_GSEA.pdf"),
       p_cmp, width = 11, height = 7, bg = "white")
cat("Wrote Sst_ORA_vs_GSEA.{png,pdf}\n")
cat("Done.\n")
