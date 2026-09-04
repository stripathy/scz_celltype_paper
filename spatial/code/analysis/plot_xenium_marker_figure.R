#!/usr/bin/env Rscript
# Supplemental Xenium marker figure (full-page width, 7.09"):
#   (a) subclass-level marker dot plot   (inhibitory -> excitatory -> non-neuronal)
#   (b) SST supertype-level marker dot plot (ordered pia->WM; Sst Chodl at right)
# ggplot2 + cowplot house style; data from {subclass,sst}_dotplot_long.csv.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(cowplot); library(scales)
})
cowplot::set_null_device("agg")
DD <- "output/depth_validation/lieber_layers"            # data inputs (dot-plot long CSVs)
FIGDIR <- "supplemental_figures"                          # consolidated cell-typing supplemental-figure home
dir.create(FIGDIR, showWarnings = FALSE)
REDS <- c("#f7f7f7", "#fdd0a2", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15")
W <- 7.09; BASE <- 7.8

dot_panel <- function(csv, sec_levels, sec_relabel = NULL, ylab_fun = function(d) d$cell_type) {
  d <- read_csv(file.path(DD, csv), show_col_types = FALSE)
  if (!is.null(sec_relabel)) d$section <- dplyr::recode(d$section, !!!sec_relabel)
  d$section <- factor(d$section, levels = sec_levels)
  d$ylab <- ylab_fun(d)
  d$ylab <- factor(d$ylab, levels = unique(d$ylab[order(d$ct_idx)]))
  d$gene <- factor(d$gene, levels = unique(d$gene[order(d$gene_idx)]))
  ggplot(d, aes(gene, ylab)) +
    geom_point(aes(size = pct, colour = scaled_mean)) +
    facet_grid(~ section, scales = "free_x", space = "free_x") +
    scale_colour_gradientn(colours = REDS, limits = c(0, 1),
                           name = "scaled mean\nexpression",
                           breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                           guide = guide_colourbar(barheight = unit(3.5, "pt"),
                                                   barwidth = unit(60, "pt"),
                                                   title.vjust = 1)) +
    scale_size_area(max_size = 2.8, limits = c(0, 1), breaks = c(.1, .25, .5, .9),
                    labels = percent_format(accuracy = 1), name = "% expressing") +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = NULL) +
    theme_cowplot(font_size = BASE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = BASE - 1.5),
          axis.text.y = element_text(size = BASE - 0.5),
          axis.line = element_line(linewidth = 0.3),
          strip.clip = "off",
          strip.background = element_rect(fill = "grey90", colour = NA),
          strip.text = element_text(size = BASE - 1, face = "bold", margin = margin(2, 1, 2, 1)),
          panel.spacing = unit(2.5, "pt"),
          panel.border = element_rect(colour = "grey85", fill = NA, linewidth = 0.3),
          legend.position = "bottom", legend.box = "horizontal",
          legend.justification = "center", legend.box.just = "center",
          legend.box.spacing = unit(4, "pt"), legend.spacing.x = unit(10, "pt"),
          legend.title = element_text(size = BASE - 1), legend.text = element_text(size = BASE - 2),
          plot.margin = margin(4, 4, 2, 4))
}

pa <- dot_panel("subclass_dotplot_long.csv",
                c("pan-Inh", "Inhibitory", "Excitatory", "Non-neuronal"))
pb <- dot_panel("sst_dotplot_long.csv",
                c("GABA/SST", "Sst"),
                sec_relabel = c(anchor = "GABA/SST"))   # supertypes stay depth-ordered; depth numbers dropped from labels

leg <- get_legend(pa)
grid <- plot_grid(pa + theme(legend.position = "none"),
                  pb + theme(legend.position = "none"),
                  ncol = 1, labels = c("a", "b"), label_size = 12, align = "v", axis = "lr",
                  rel_heights = c(1.18, 1))
fig <- plot_grid(grid, leg, ncol = 1, rel_heights = c(1, 0.065))

ggsave(file.path(FIGDIR, "xenium_marker_figure.png"), fig, width = W, height = 9.2, dpi = 400, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIGDIR, "xenium_marker_figure.pdf"), fig, width = W, height = 9.2, bg = "white", device = pdf)
cat("saved", file.path(FIGDIR, "xenium_marker_figure.png/.pdf"), "\n")
