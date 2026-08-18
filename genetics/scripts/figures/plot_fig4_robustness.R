# Supplementary figure: is Figure 4a robust to the two discretionary choices
# behind it -- which SCZ GWAS, and which cortical region supplies the normotypic
# expression?
#
# Each facet re-runs the panel end to end under one configuration. The taxonomy,
# the specificity recipe, the MAGMA settings and the compositional data are held
# fixed throughout; only the named input changes. The y-axis is therefore
# identical across facets and only x moves.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(ggrepel); library(cowplot)
})

TABLES <- "/Users/shreejoy/Github/scz_celltype_paper/genetics/results/tables"
FIGDIR <- "/Users/shreejoy/Github/scz_celltype_paper/genetics/results/figures"
BASE <- 12

runs <- c("DLPFC + Bigdeli (main figure)", "DLPFC + PGC3", "MTG + Bigdeli")
d <- read_csv(file.path(TABLES, "fig4_robustness_sst16.csv"), show_col_types = FALSE) |>
  pivot_longer(all_of(runs), names_to = "run", values_to = "nlp") |>
  mutate(run = factor(run, levels = runs),
         depletion = -comp_beta,
         status = ifelse(depleted_fdr20, "Depleted in SCZ (FDR < 0.20)", "Not depleted"))

# plotmath-safe P label: 3 decimals down to 0.001, scientific below that
fmt_p <- function(p) ifelse(p >= 1e-3, sprintf("italic(P) == %.3f", p),
  ifelse(p < 1e-4, "italic(P) < 10^-4",     # cor.test's exact Spearman P underflows here
    sprintf("italic(P) == %.1f %%*%% 10^%.0f", p / 10^floor(log10(p)), floor(log10(p)))))

# one inset per facet, pinned to the lower-right corner -- the only region
# empty in all three facets, since depletion rises with enrichment
stats_df <- d |> group_by(run) |>
  summarise(rho = cor(nlp, depletion, method = "spearman"),
            p = suppressWarnings(cor.test(nlp, depletion, method = "spearman")$p.value),
            x = max(nlp), y = min(depletion), .groups = "drop") |>
  mutate(label = sprintf("rho == %.2f * ',' ~ %s", rho, fmt_p(p)))

# label only the supertypes a reader tracks between facets
KEEP <- c("Sst_23", "Sst_2", "Sst_3", "Sst_20", "Sst_25", "Sst_22")

p <- ggplot(d, aes(nlp, depletion)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70", linewidth = 0.25) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#444", fill = "#888",
              alpha = 0.13, linewidth = 0.35, linetype = "dashed") +
  geom_point(aes(fill = color, colour = status, stroke = status),
             shape = 21, size = 3.1, alpha = 0.95) +
  geom_text_repel(aes(label = ifelse(supertype %in% KEEP, supertype, "")),
                  size = 3.0, colour = "grey15", box.padding = 0.3,
                  point.padding = 0.25, max.overlaps = 30, force = 6,
                  min.segment.length = 0, segment.size = 0.2,
                  segment.colour = "grey65", seed = 1) +
  geom_text(data = stats_df, aes(x = x, y = y, label = label),
            parse = TRUE, hjust = 1.0, vjust = -0.4, size = 3.4, inherit.aes = FALSE) +
  facet_wrap(~run, nrow = 1, scales = "free_x") +
  scale_fill_identity() +
  scale_colour_manual(values = c("Depleted in SCZ (FDR < 0.20)" = "black",
                                 "Not depleted" = "grey60")) +
  scale_discrete_manual("stroke",
                        values = c("Depleted in SCZ (FDR < 0.20)" = 0.9,
                                   "Not depleted" = 0.25)) +
  guides(colour = guide_legend(override.aes = list(fill = "grey85", size = 3)),
         stroke = "none") +
  labs(x = expression("SCZ GWAS enrichment ("*-log[10]~italic(P)*")"),
       y = expression("Cell depletion in SCZ ("*-beta*")")) +
  theme_cowplot(font_size = BASE) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = BASE, face = "plain"),
        legend.position = "bottom", legend.title = element_blank(),
        legend.justification = "center",
        panel.spacing = unit(14, "pt"),
        plot.margin = margin(6, 10, 4, 6))

# ---------------------------------------------------------------- row 2
# The AD axis of Figure 4i in the SEA-AD DLPFC cohort (as in the main figure)
# and recomputed on the SEA-AD MTG cohort. Same donors, same crumblr model and
# covariates; only the dissected region differs, so this asks whether the
# cross-disorder concordance is a property of the supertypes or of the region
# they were counted in. Drawn as two facets with the same geometry as row 1.
ad_regions <- c("SEA-AD DLPFC (main figure)", "SEA-AD MTG")
ad <- read_csv(file.path(TABLES, "fig4_robustness_ad_mtg.csv"), show_col_types = FALSE) |>
  left_join(distinct(select(d, supertype, color, status)), by = "supertype") |>
  transmute(supertype, color, status, scz = -scz_beta,
            `SEA-AD DLPFC (main figure)` = -ad_slope, `SEA-AD MTG` = -ad_slope_mtg) |>
  pivot_longer(all_of(ad_regions), names_to = "region", values_to = "ad") |>
  mutate(region = factor(region, levels = ad_regions))

ad_stats <- ad |> group_by(region) |>
  summarise(rho = cor(scz, ad, method = "spearman"),
            p = suppressWarnings(cor.test(scz, ad, method = "spearman")$p.value),
            x = max(scz), y = min(ad), .groups = "drop") |>
  mutate(label = sprintf("rho == %.2f * ',' ~ %s", rho, fmt_p(p)))

p2 <- ggplot(ad, aes(scz, ad)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70", linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70", linewidth = 0.25) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#444", fill = "#888",
              alpha = 0.13, linewidth = 0.35, linetype = "dashed") +
  geom_point(aes(fill = color, colour = status, stroke = status),
             shape = 21, size = 3.1, alpha = 0.95) +
  geom_text_repel(aes(label = supertype), size = 3.0, colour = "grey15",
                  box.padding = 0.3, point.padding = 0.25, max.overlaps = 30,
                  force = 6, min.segment.length = 0, segment.size = 0.2,
                  segment.colour = "grey65", seed = 1) +
  geom_text(data = ad_stats, aes(x = x, y = y, label = label),
            parse = TRUE, hjust = 1.0, vjust = -0.4, size = 3.4, inherit.aes = FALSE) +
  facet_wrap(~region, nrow = 1, scales = "free_y") +
  scale_fill_identity() +
  scale_colour_manual(values = c("Depleted in SCZ (FDR < 0.20)" = "black",
                                 "Not depleted" = "grey60"), guide = "none") +
  scale_discrete_manual("stroke",
                        values = c("Depleted in SCZ (FDR < 0.20)" = 0.9,
                                   "Not depleted" = 0.25), guide = "none") +
  labs(x = expression("Cell depletion in SCZ ("*-beta*")"),
       y = expression("Cell depletion in AD ("*-beta*"/SD CPS)")) +
  theme_cowplot(font_size = BASE) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = BASE, face = "plain"),
        panel.spacing = unit(14, "pt"),
        plot.margin = margin(6, 10, 4, 6))

fig <- plot_grid(
  plot_grid(p, labels = "a", label_size = 15),
  plot_grid(p2, NULL, nrow = 1, rel_widths = c(2.08, 0.92), labels = c("b", ""),
            label_size = 15),
  ncol = 1, rel_heights = c(1, 0.9))

for (ext in c("png", "pdf")) {
  ggsave(file.path(FIGDIR, paste0("supp_fig4a_robustness.", ext)), fig,
         width = 10.6, height = 8.0, bg = "white", dpi = 220)
}
print(as.data.frame(ad_stats[, c("region", "rho", "p")]), digits = 3)
cat("wrote supp_fig4a_robustness.{png,pdf}\n")
print(as.data.frame(stats_df[, c("run", "rho", "p")]), digits = 3)
