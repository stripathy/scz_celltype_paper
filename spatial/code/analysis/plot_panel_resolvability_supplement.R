#!/usr/bin/env Rscript
# Supplementary figure: panel-choice sensitivity — cell-type resolvability on the
# 300-gene Xenium panel vs the full transcriptome.
#   (a) subclass identity (24-way), faceted by class
#   (b) Sst supertype resolution (16 types), Sst_25 highlighted
# Dumbbell: classification F1 on the panel (blue) vs the transcriptome (dark),
# nearest-centroid Pearson classifier, leave-one-donor-out CV on the SEA-AD snRNA
# reference. House style (ggplot2 + cowplot), large legible text.
# Data: build via resolvability_report.py -> data/resolvability_{subclass,sst_supertype}.csv
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(cowplot)
})
cowplot::set_null_device("agg")

DD <- "output/celltyping_supplement/data"; OUT <- "supplemental_figures"   # consolidated supplemental-figure home
dir.create(OUT, showWarnings = FALSE)
PANEL_COL <- "#3182BD"; TX_COL <- "#252525"; SEG <- "#C7C7C7"; HI <- "Sst_25"; HI_COL <- "#B8860B"
BASE <- 13                                   # large base font (legible supplement)
SETLAB <- c(panel = "Xenium panel (300 genes)", all = "full transcriptome")

mklong <- function(d) d |>
  select(cell_type, F1_panel, F1_all, any_of("class")) |>
  pivot_longer(c(F1_panel, F1_all), names_to = "set", values_to = "F1") |>
  mutate(set = factor(SETLAB[sub("F1_", "", set)], levels = SETLAB))

build <- function(d, order_vec, facet = FALSE, hi = NULL) {
  ord <- rev(order_vec)                                       # factor-level order; order_vec[1] lands at top
  d$cell_type <- factor(d$cell_type, levels = ord)
  L <- mklong(d); L$cell_type <- factor(L$cell_type, levels = ord)
  g <- ggplot() +
    geom_segment(data = d, aes(y = cell_type, yend = cell_type, x = F1_panel, xend = F1_all),
                 color = SEG, linewidth = 1.3) +
    geom_point(data = L, aes(x = F1, y = cell_type, fill = set),
               shape = 21, color = "white", stroke = 0.45, size = 3.3) +
    scale_fill_manual(values = setNames(c(PANEL_COL, TX_COL), SETLAB), name = NULL) +
    scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = c(0.02, 0.04))) +
    labs(x = "classification F1  (leave-one-donor-out CV)", y = NULL) +
    theme_cowplot(font_size = BASE) +
    theme(axis.text.y = element_text(size = BASE - 1.5),
          axis.text.x = element_text(size = BASE - 1),
          axis.title.x = element_text(size = BASE),
          axis.line = element_line(linewidth = 0.4),
          panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
          legend.position = "bottom", legend.text = element_text(size = BASE - 1))
  if (facet)
    g <- g + facet_grid(rows = vars(class), scales = "free_y", space = "free_y") +
      theme(strip.background = element_rect(fill = "grey92", color = NA),
            strip.text = element_text(size = BASE - 1, face = "bold"))
  if (!is.null(hi)) {
    yc <- ifelse(ord == hi, HI_COL, "black"); yf <- ifelse(ord == hi, "bold", "plain")
    g <- g + theme(axis.text.y = element_text(color = yc, face = yf, size = BASE - 1.5))
  }
  g
}

source("code/analysis/_xenium_exclusions.R")   # types too rare to analyse here (SM1)

sc <- read_csv(file.path(DD, "resolvability_subclass.csv"), show_col_types = FALSE) |>
  drop_too_rare("cell_type", "subclass")
sc$class <- factor(sc$class, levels = c("GABAergic", "Glutamatergic", "Non-neuronal"))
ss <- read_csv(file.path(DD, "resolvability_sst_supertype.csv"), show_col_types = FALSE) |>
  drop_too_rare("cell_type", "supertype")

# --- match the marker dot-plot cell-type ordering (single source of truth = dot-plot long CSVs) ---
DOT <- "output/depth_validation/lieber_layers"
sc_order <- read_csv(file.path(DOT, "subclass_dotplot_long.csv"), show_col_types = FALSE) |>
  distinct(cell_type, ct_idx) |> arrange(ct_idx) |> pull(cell_type)   # inhibitory -> excitatory -> non-neuronal
# Lamp5 Lhx6 is deliberately absent from both the dot plot and this panel — see
# _xenium_exclusions.R and SM1. It is NOT restored into the ordering here.
sst_order <- read_csv(file.path(DOT, "sst_dotplot_long.csv"), show_col_types = FALSE) |>
  distinct(cell_type, ct_idx) |> arrange(ct_idx) |> pull(cell_type)   # pia -> WM (median predicted depth)

pa <- build(sc, order_vec = sc_order, facet = TRUE)
pb <- build(ss, order_vec = sst_order, hi = HI)
leg <- get_legend(pa + theme(legend.position = "bottom"))
grid <- plot_grid(pa + theme(legend.position = "none"),
                  pb + theme(legend.position = "none"),
                  ncol = 1, labels = c("a", "b"), label_size = 18, label_fontface = "bold",
                  rel_heights = c(1.55, 1), align = "v", axis = "lr")
fig <- plot_grid(grid, leg, ncol = 1, rel_heights = c(1, 0.04))

ggsave(file.path(OUT, "panel_resolvability_supplement.png"), fig,
       width = 7.2, height = 10.6, dpi = 400, bg = "white", device = ragg::agg_png)
ggsave(file.path(OUT, "panel_resolvability_supplement.pdf"), fig,
       width = 7.2, height = 10.6, bg = "white", device = pdf)
cat("saved panel_resolvability_supplement.png/.pdf\n")
