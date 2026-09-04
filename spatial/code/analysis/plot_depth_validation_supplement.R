#!/usr/bin/env Rscript
# Depth-model cross-validation supplement (held-out SEA-AD MERFISH donors).
# Purpose: show the depth assignments are NOT circular -- the model generalizes
# to held-out donors AND carries depth signal beyond the cell-type label.
# Reads output/depth_validation/*.csv from validate_depth_model_cv.py.
#
# House style: transcriptomic Fig 09 / genetics composite -- theme_cowplot(7),
# axis titles 6.5 pt / text 6 pt, lowercase bold panel labels (8 pt),
# cowplot::plot_grid, plain-text plotmath stats (no box), 7.1 in width.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(cowplot); library(ggrepel); library(scales)
})
cowplot::set_null_device("agg")

BASE <- 7
CLASS_COL <- c("Glutamatergic" = "#117733", "GABAergic" = "#882255", "Non-neuronal" = "#DDCC77")
ACCENT <- "#0072B2"; HILITE <- "#D55E00"; GREY <- "grey70"

theme_house <- function(base_size = BASE, grid = FALSE) {
  theme_cowplot(font_size = base_size) +
    theme(plot.subtitle = element_text(size = base_size, hjust = 0.5, margin = margin(b = 1)),
          axis.title  = element_text(size = base_size - 0.5),
          axis.text   = element_text(size = base_size - 1, color = "black"),
          axis.line   = element_line(linewidth = 0.3, color = "grey20"),
          axis.ticks  = element_line(linewidth = 0.25, color = "grey20"),
          legend.text = element_text(size = base_size - 1.5),
          legend.title = element_text(size = base_size - 1),
          legend.key.size = unit(8, "pt"),
          panel.grid.major = if (grid) element_line(color = "grey92", linewidth = 0.15) else element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(3, 5, 3, 3))
}

.args <- commandArgs(trailingOnly = TRUE)
DD <- if (length(.args) > 0) .args[1] else "output/depth_validation"
summ <- read_csv(file.path(DD, "cv_summary.csv"), show_col_types = FALSE)
gv <- function(k) summ$value[summ$metric == k]
oof_r2 <- gv("oof_r2"); oof_mae <- gv("oof_mae"); oof_r <- gv("oof_pearson_r"); n_don <- gv("n_donors")

per_donor <- read_csv(file.path(DD, "cv_per_donor.csv"), show_col_types = FALSE)
dec       <- read_csv(file.path(DD, "cv_decomposition.csv"), show_col_types = FALSE)
wsub      <- read_csv(file.path(DD, "cv_within_subclass.csv"), show_col_types = FALSE)
samp      <- read_csv(file.path(DD, "cv_predictions_sample.csv"), show_col_types = FALSE)
# pooled within-subclass r (subclass-mean removed) computed directly from the OOF sample
mres <- samp$manual - ave(samp$manual, samp$subclass)
pres <- samp$pred   - ave(samp$pred,   samp$subclass)
pooled_within_r <- cor(mres, pres)

# ---- (a) out-of-fold per-cell density: predicted vs manual depth ----
pa <- ggplot(samp, aes(manual, pred)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_bin2d(bins = 80) +
  scale_fill_viridis_c(option = "mako", trans = "log10", direction = -1, name = "cells",
                       guide = guide_colorbar(barwidth = unit(4, "pt"), barheight = unit(26, "pt"))) +
  annotate("text", x = 0.02, y = 1.10, hjust = 0, vjust = 1, size = 2.4, lineheight = 0.95,
           label = sprintf("italic(R)^2 == '%.3f'", oof_r2), parse = TRUE) +
  annotate("text", x = 0.02, y = 1.00, hjust = 0, vjust = 1, size = 2.4,
           label = sprintf("MAE == '%.3f'", oof_mae), parse = TRUE) +
  scale_x_continuous(limits = c(-0.02, 1.14), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(-0.05, 1.16), breaks = seq(0, 1, 0.25)) +
  labs(subtitle = "Held-out MERFISH cells",
       x = "Manual SEA-AD depth", y = "Predicted depth (held-out)") +
  theme_house() + theme(aspect.ratio = 1, legend.position = "right")

# ---- (b) per-donor held-out R^2 ----
per_donor <- per_donor %>% arrange(r2) %>% mutate(donor = factor(donor, levels = donor))
pb <- ggplot(per_donor, aes(r2, donor)) +
  geom_vline(xintercept = oof_r2, linetype = "dashed", color = HILITE, linewidth = 0.3) +
  geom_segment(aes(x = min(per_donor$r2) - 0.02, xend = r2, yend = donor),
               color = "grey82", linewidth = 0.3) +
  geom_point(color = ACCENT, size = 1.0) +
  annotate("text", x = oof_r2, y = 1.5, label = sprintf(" overall %.3f", oof_r2),
           hjust = 0, vjust = 0, size = 2.1, color = HILITE) +
  labs(subtitle = sprintf("Held-out test donors (n = %d)", n_don),
       x = expression("Held-out"~italic(R)^2), y = NULL) +
  theme_house() + theme(axis.text.y = element_text(size = BASE - 2.6))

# ---- (c) variance decomposition: subclass-only vs neighborhood ----
dec <- dec %>% mutate(short = ifelse(grepl("Subclass", model), "Subclass\nidentity", "Neighborhood\nmodel"),
                      short = factor(short, levels = c("Subclass\nidentity", "Neighborhood\nmodel")))
pc <- ggplot(dec, aes(short, r2, fill = short)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", r2)), vjust = -0.4, size = 2.5) +
  scale_fill_manual(values = c("Subclass\nidentity" = GREY, "Neighborhood\nmodel" = ACCENT), guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.1))) +
  labs(subtitle = "Depth variance explained", x = NULL, y = expression("Out-of-fold"~italic(R)^2)) +
  theme_house() + theme(axis.text.x = element_text(size = BASE - 1))

# ---- (d) within-subclass depth resolution (subclass-mean removed) ----
wsub <- wsub %>% arrange(within_r) %>% mutate(subclass = factor(subclass, levels = subclass))
pd_ <- ggplot(wsub, aes(within_r, subclass, fill = klass)) +
  geom_col(width = 0.72) +
  annotate("text", x = 0.02, y = nrow(wsub) - 0.3, hjust = 0, vjust = 1, size = 2.3,
           label = sprintf("pooled~italic(r) == '%.2f'", pooled_within_r), parse = TRUE) +
  scale_fill_manual(values = CLASS_COL, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.03))) +
  labs(subtitle = "Within-subclass depth resolution",
       x = "Within-subclass Pearson r", y = NULL) +
  theme_house() + theme(axis.text.y = element_text(size = BASE - 2.6), legend.position = "none")

leg <- get_legend(pd_ + theme(legend.position = "bottom") +
                  guides(fill = guide_legend(override.aes = list(size = 3))))
grid <- plot_grid(pa, pb, pc, pd_, ncol = 2, labels = c("a", "b", "c", "d"),
                  label_size = 8, label_fontface = "bold", align = "hv", axis = "tblr")
fig <- plot_grid(grid, leg, ncol = 1, rel_heights = c(1, 0.05))

ggsave(file.path(DD, "depth_validation_supplement.png"), fig,
       width = 7.1, height = 6.4, dpi = 400, bg = "white", device = ragg::agg_png)
ggsave(file.path(DD, "depth_validation_supplement.pdf"), fig,
       width = 7.1, height = 6.4, bg = "white", device = pdf)
cat(sprintf("OOF R2=%.3f MAE=%.3f | subclass-only R2=%.2f | pooled within-r=%.2f\n",
            oof_r2, oof_mae, dec$r2[dec$short=="Subclass\nidentity"], pooled_within_r))
