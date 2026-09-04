#!/usr/bin/env Rscript
# Concordance panel v2: foreground the MODULE-level SCZ/AD divergence in the
# depleted stratum. Left: gene scatter with full module coloring + centroid
# markers (translation/OxPhos displaced on the SCZ axis, not the AD axis).
# Right: module x disorder slope panel (standardized mean effect ± 95% CI) —
# the interaction in one glance; SST and VGF drawn as the shared-core exception.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(cowplot); library(ggrepel)
})
OUT  <- "transcriptomic/results/sst_strata_gsea"
AOUT <- file.path(OUT, "seaad_a9")
FS <- 14
FAM_COLS <- c(Synaptic = "#0072B2", Translation = "#CC79A7",
              `OxPhos/mito` = "#D55E00", `SST / VGF` = "#009E73", Other = "grey85")

g   <- read_csv(file.path(OUT, "gsea_all_signatures.csv"), show_col_types = FALSE)
sig <- read_csv(file.path(OUT, "stratum_gene_signatures.csv"), show_col_types = FALSE)
a9  <- read_csv(file.path(AOUT, "de_a9_strata.csv"), show_col_types = FALSE) |>
  filter(signature == "a9_depleted") |> select(gene, t_ad = t)

blocks <- list(
  Synaptic = c("GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING", "GOCC_PRESYNAPSE",
               "GOBP_NEUROTRANSMITTER_SECRETION", "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES"),
  Translation = c("REACTOME_TRANSLATION", "GOCC_RIBOSOMAL_SUBUNIT",
                  "GOBP_CYTOPLASMIC_TRANSLATION", "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION"),
  `OxPhos/mito` = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "GOBP_OXIDATIVE_PHOSPHORYLATION",
                    "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"))
le_union <- function(pws) g |> filter(pathway %in% pws, signature == "depleted") |>
  pull(leadingEdge) |> strsplit("|", fixed = TRUE) |> unlist() |> unique()
mods <- lapply(blocks, le_union)

m <- sig |> filter(stratum == "depleted") |> select(gene, z_scz = z) |>
  inner_join(a9, by = "gene") |>
  mutate(module = case_when(gene %in% c("SST", "VGF") ~ "SST / VGF",
                            gene %in% mods$Translation ~ "Translation",
                            gene %in% mods$`OxPhos/mito` ~ "OxPhos/mito",
                            gene %in% mods$Synaptic ~ "Synaptic",
                            TRUE ~ "Other"),
         module = factor(module, names(FAM_COLS)))

cent <- m |> filter(module != "Other") |>
  group_by(module) |>
  summarise(n = n(), x = mean(z_scz), y = mean(t_ad),
            sex_ = sd(z_scz)/sqrt(n()), sey = sd(t_ad)/sqrt(n()), .groups = "drop")
cat("module centroids (depleted stratum):\n"); print(cent |> mutate(across(where(is.numeric), ~round(.x, 2))))

# ---- left: scatter with module clouds + centroids ------------------------------
hl <- m |> filter(gene %in% c("SST", "VGF"))
p1 <- ggplot(m, aes(z_scz, t_ad, colour = module)) +
  geom_point(data = ~ filter(.x, module == "Other"), size = 0.45, alpha = 0.25) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_abline(linetype = "dashed", colour = "grey60") +
  geom_point(data = ~ filter(.x, !module %in% c("Other", "SST / VGF")),
             size = 1.3, alpha = 0.65) +
  geom_point(data = cent |> filter(module != "SST / VGF"),
             aes(x, y, fill = module), colour = "black", shape = 23, size = 5,
             stroke = 0.9, show.legend = FALSE) +
  geom_point(data = hl, colour = "black", size = 2.6, shape = 21, stroke = 0.9,
             aes(fill = module), show.legend = FALSE) +
  geom_text_repel(data = hl, aes(label = gene), colour = "black", size = 4.4,
                  fontface = "italic", seed = 11, nudge_y = -0.7) +
  geom_text_repel(data = cent |> filter(module != "SST / VGF"),
                  aes(x, y, label = module), colour = "black", size = 4.4,
                  fontface = "bold", seed = 12, nudge_y = 0.9,
                  show.legend = FALSE) +
  scale_colour_manual(values = FAM_COLS, name = NULL,
                      breaks = c("Synaptic", "Translation", "OxPhos/mito")) +
  scale_fill_manual(values = FAM_COLS, guide = "none") +
  labs(x = "Gene z, SCZ depleted stratum (case vs control)",
       y = "Gene t, AD depleted stratum (along CPS)") +
  theme_cowplot(font_size = FS) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3, alpha = 1)))

# ---- right: module x disorder slope panel (standardized) -----------------------
sdz <- sd(m$z_scz); sdt <- sd(m$t_ad)
slope <- bind_rows(
  m |> filter(module != "Other") |> group_by(module) |>
    summarise(disorder = "SCZ", mean = mean(z_scz)/sdz,
              se = sd(z_scz)/sdz/sqrt(n()), .groups = "drop"),
  m |> filter(module != "Other") |> group_by(module) |>
    summarise(disorder = "AD", mean = mean(t_ad)/sdt,
              se = sd(t_ad)/sdt/sqrt(n()), .groups = "drop")) |>
  mutate(disorder = factor(disorder, c("SCZ", "AD")))
p2 <- ggplot(slope, aes(disorder, mean, colour = module, group = module)) +
  geom_hline(yintercept = 0, colour = "grey75") +
  geom_line(linewidth = 1.1) +
  geom_pointrange(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se),
                  size = 0.7, linewidth = 0.9) +
  geom_text_repel(data = slope |> filter(disorder == "AD"),
                  aes(label = module), size = 4.2, fontface = "bold",
                  nudge_x = 0.18, direction = "y", hjust = 0, seed = 13,
                  show.legend = FALSE) +
  scale_colour_manual(values = FAM_COLS, guide = "none") +
  scale_x_discrete(expand = expansion(mult = c(0.15, 0.55))) +
  labs(x = NULL, y = "Standardized mean effect in depleted stratum\n(mean / gene-wide SD, ± 95% CI)") +
  theme_cowplot(font_size = FS)

fig <- plot_grid(p1, p2, nrow = 1, rel_widths = c(0.62, 0.38),
                 labels = c("", ""), label_size = 15)
ggsave(file.path(AOUT, "figS_gene_concordance_v2.png"), fig,
       width = 13.5, height = 6.4, dpi = 200, bg = "white")
cat("wrote figS_gene_concordance_v2.png\n")
