#!/usr/bin/env Rscript
# Supplementary figure | Diagnosis x depletion-group interaction.
#
# Supports the results sentence "Formal models testing the interaction between
# depletion group status and SCZ diagnosis supported these findings, with effects
# concentrated in protein synthesis ... and oxidative phosphorylation ...", which
# currently cites no figure. Every donor contributes one pseudobulk per group, so
# the contrast is taken WITHIN donor and donor-level factors cancel (Methods).
#
#   a  interaction NES for the gene sets shown in Fig. 5d, plus the two sets
#      named in the text, which are NOT both in Fig. 5d
#   b  the 15 gene sets with the strongest interaction, showing that the effect
#      is concentrated in translation rather than spread across the transcriptome
#
# A negative NES means the diagnosis effect is more negative in the depleted
# group than in the non-depleted group.
source("transcriptomic/scripts/fig5/_common.R")
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })

ig <- read_csv(file.path(P$pb, "interaction_gsea.csv"), show_col_types = FALSE) |>
  filter(coef == "dxSCZ:stratumdepleted")
gs <- msigdb_sets()$sets

CITED <- c(REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE = "translation",
           GOBP_ELECTRON_TRANSPORT_CHAIN = "oxphos")
BL <- c(synaptic = "Synaptic", ubiquitin = "Ubiquitin",
        translation = "Cytosolic translation", oxphos = "Oxidative phosphorylation")
pretty_set <- function(x) {
  y <- x |> str_remove("^GOBP_|^GOCC_|^GOMF_|^REACTOME_") |>
    str_replace_all("_", " ") |> str_to_sentence()
  y <- str_replace(y, "^Srp dependent cotranslational protein targeting to membrane",
                   "SRP-dependent cotranslational protein targeting")
  y <- str_replace(y, "^Response of eif2ak4 gcn2 to amino acid deficiency",
                   "GCN2 response to amino-acid deficiency")
  y <- str_replace(y, " nmd$", " (NMD)")
  y <- str_replace(y, "k63", "K63")
  y <- str_replace(y, "slits and robos", "SLITs and ROBOs")
  y
}

# ---- a | the Fig. 5d blocks, plus the two sets the text names ----------------
d_a <- bind_rows(
  BLOCKS |> select(pathway, block) |> mutate(cited = FALSE),
  tibble(pathway = names(CITED), block = unname(CITED), cited = TRUE)) |>
  inner_join(ig, by = "pathway") |>
  mutate(block = factor(BL[block], unname(BL)),
         lab = pretty_set(pathway),
         stars = stars_of(padj))
say("panel a values:")
print(d_a |> transmute(block, lab = str_trunc(lab, 46), NES = round(NES, 2),
                       padj = signif(padj, 2), cited))

p_a <- ggplot(d_a, aes(NES, reorder(lab, NES), colour = block)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_segment(aes(x = 0, xend = NES, yend = reorder(lab, NES)), linewidth = 0.7 * GEO) +
  geom_point(aes(shape = cited), size = 2.4 * GEO) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  geom_text(aes(label = stars), hjust = ifelse(d_a$NES < 0, 1.4, -0.4),
            size = LBL_SMALL, colour = "grey20") +
  scale_colour_manual(values = setNames(unname(FAM_COLS[names(BL)]), unname(BL)), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.16, 0.16))) +
  labs(x = "Interaction NES (depleted vs non-depleted)", y = NULL) +
  theme_fig5() +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key.size = unit(7, "pt"),
        axis.text.y = element_text(size = LBL_SMALL * 2.845))

# ---- b | the strongest interactions overall ---------------------------------
# Several high-ranking Reactome sets carry names that do not suggest the ribosome
# (influenza infection, SLIT/ROBO expression, selenoamino acid metabolism). They
# are dominated by the same ribosomal-protein genes, so the panel marks the
# proportion of each set that is ribosomal rather than letting the names mislead.
top <- ig |> arrange(padj) |> head(15) |>
  mutate(rp_frac = map_dbl(pathway, ~ mean(str_detect(gs[[.x]], "^RP[LS][0-9]"))),
         lab = pretty_set(pathway))
say("panel b, ribosomal-protein content of the top sets:")
print(top |> transmute(lab = str_trunc(lab, 50), padj = signif(padj, 2),
                       rp_pct = round(100 * rp_frac)))

p_b <- ggplot(top, aes(-log10(padj), reorder(lab, -log10(padj)), fill = 100 * rp_frac)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * rp_frac)), hjust = -0.25,
            size = LBL_SMALL, colour = "grey25") +
  scale_fill_gradient(low = "grey80", high = FAM_COLS[["translation"]],
                      limits = c(0, 100), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = expression("Interaction "*-log[10]~"FDR"), y = NULL,
       subtitle = "% ribosomal protein genes") +
  theme_fig5() +
  theme(plot.subtitle = element_text(size = LEGEND_TEXT, colour = "grey35"),
        axis.text.y = element_text(size = LBL_SMALL * 2.845))

fig <- plot_grid(p_a, p_b, nrow = 1, rel_widths = c(0.5, 0.5), align = "h", axis = "tb",
                 labels = c("a", "b"), label_size = PANEL_LABEL)
ggsave(file.path(P$supp, "figS_interaction.png"), fig, width = FIG_W, height = 3.4,
       dpi = FIG_DPI, bg = "white")
ggsave(file.path(P$supp, "figS_interaction.pdf"), fig, width = FIG_W, height = 3.4, bg = "white")
say("wrote figS_interaction.(png|pdf)")
