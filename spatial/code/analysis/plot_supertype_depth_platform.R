#!/usr/bin/env Rscript
# Xenium vs SEA-AD MERFISH cortical depth per NEURONAL supertype.
# Same layout as plot_supertype_depth_casecontrol.R — (a) Glutamatergic,
# (b) GABAergic, supertypes grouped by subclass and ordered pia->WM — but the
# violins are split by PLATFORM rather than diagnosis.
#
# The two depths are different kinds of measurement and the figure is meant to
# show that: MERFISH depth is SEA-AD's hand annotation, Xenium depth is the
# neighbourhood model's prediction, and the model was trained on the MERFISH
# annotation. Agreement in position is therefore expected; what the figure adds
# over the median-depth scatter in S2 is where that agreement breaks down, and
# the systematic narrowing of the Xenium distributions (a model predicting a
# conditional mean shrinks spread).
#
# Deliberately no per-supertype significance test: with 10^5-10^6 cells every
# supertype would reach any p threshold, which would say nothing. Effect size
# (difference in median depth) is reported instead.
#
# Data: build_supertype_depth_platform_data.py -> output/depth_platform/
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(cowplot)
})
cowplot::set_null_device("agg")

DD     <- "output/depth_platform"
FIGDIR <- "supplemental_figures"; dir.create(FIGDIR, showWarnings = FALSE)

# ColorBrewer PuOr — deliberately not the blue/red of the Control-vs-SCZ depth
# figure, so the two are never mistaken for each other at a glance
PLAT_COL <- c(MERFISH = "#5E3C99", Xenium = "#E66101")
BIG_DELTA <- 0.10                                   # |median difference| worth flagging
LAYER_B <- c(.10, .40, .55, .70, .90)
LAYER_T <- c(.05, .25, .475, .625, .80, .95); LAYER_N <- c("L1","L2/3","L4","L5","L6","WM")
GLUT_ORDER <- c("L2/3 IT","L4 IT","L5 IT","L5 ET","L5/6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")
GABA_ORDER <- c("Lamp5","Sncg","Vip","Pax6","Chandelier","Pvalb","Sst","Sst Chodl")
SUBORDER <- list(Glutamatergic = GLUT_ORDER, GABAergic = GABA_ORDER)
BASE <- 7           # 7.1 in canvas: matches the other Xenium supplemental figures

df  <- read_csv(file.path(DD, "supertype_depth_platform.csv"), show_col_types = FALSE)
sm  <- read_csv(file.path(DD, "supertype_depth_platform_summary.csv"), show_col_types = FALSE)
df$platform <- factor(df$platform, levels = c("MERFISH", "Xenium"))

r_med <- cor(sm$median_MERFISH, sm$median_Xenium)
cat(sprintf("%d supertypes | median-depth r = %.3f | median |delta| = %.3f | |delta| > %.2f: %d\n",
            nrow(sm), r_med, median(abs(sm$delta_median)), BIG_DELTA,
            sum(abs(sm$delta_median) > BIG_DELTA)))
cat(sprintf("  IQR: MERFISH %.3f vs Xenium %.3f (Xenium %.0f%% narrower)\n",
            median(sm$iqr_MERFISH), median(sm$iqr_Xenium),
            100 * (1 - median(sm$iqr_Xenium) / median(sm$iqr_MERFISH))))

shorten <- function(sup, sub) sub(paste0("^", sub, "[_ ]"), "", sup)

plot_class <- function(cls) {
  ro <- sm |> filter(class == cls) |>
    mutate(subclass = factor(subclass, levels = SUBORDER[[cls]])) |>
    arrange(subclass, median_MERFISH) |>
    mutate(xpos = row_number(), short = mapply(shorten, supertype, as.character(subclass)))
  d <- df |> filter(supertype %in% ro$supertype) |>
    left_join(ro |> select(supertype, xpos), by = "supertype")
  # inner-95% clip per supertype x platform (display only; summaries used all cells)
  dc <- d |> group_by(supertype, platform) |>
    filter(depth >= quantile(depth, .025), depth <= quantile(depth, .975)) |> ungroup()
  bnds <- ro |> group_by(subclass) |>
    summarise(mx = max(xpos), mid = mean(xpos), .groups = "drop")
  vlines <- head(sort(bnds$mx), -1) + 0.5
  lab_y <- -0.025; ylim_top <- -0.215

  # No in-plot statistics. The per-supertype median differences are printed to the
  # console and written to supertype_depth_platform_summary.csv, for the figure
  # legend to draw on; the figure itself carries only the data.
  ggplot(dc, aes(x = xpos, y = depth, fill = platform,
                 group = interaction(xpos, platform))) +
    geom_violin(position = position_dodge(0.85), width = 0.92, alpha = 0.9,
                linewidth = 0.07, color = "grey35", scale = "width") +
    stat_summary(fun = median, geom = "point", position = position_dodge(0.85),
                 size = 0.18, color = "white") +
    geom_hline(yintercept = LAYER_B, linetype = "dashed", color = "grey80", linewidth = 0.2) +
    geom_vline(xintercept = vlines, color = "grey55", linewidth = 0.2) +
    annotate("text", x = bnds$mid, y = lab_y, label = as.character(bnds$subclass),
             fontface = "bold", size = 1.75, vjust = 0.5, hjust = 0, angle = 90) +
    scale_y_reverse(breaks = LAYER_T, labels = LAYER_N, expand = expansion(mult = c(0.02, 0))) +
    scale_x_continuous(breaks = ro$xpos, labels = ro$short, expand = expansion(add = 0.7)) +
    scale_fill_manual(values = PLAT_COL, name = NULL) +
    coord_cartesian(ylim = c(1.02, ylim_top)) +
    labs(x = NULL, y = "cortical depth (pia at top, WM at bottom)", title = cls) +
    theme_cowplot(font_size = BASE) +
    theme(axis.text.x = element_text(size = BASE - 2.5, angle = 90, vjust = 0.5, hjust = 1,
                                     margin = margin(t = 0.5)),
          axis.text.y = element_text(size = BASE - 1),
          axis.title.y = element_text(size = BASE - 0.5),
          axis.line = element_line(linewidth = 0.25),
          axis.ticks = element_line(linewidth = 0.2),
          axis.ticks.length = unit(1.5, "pt"),
          plot.title = element_text(size = BASE + 0.5, face = "bold"),
          legend.position = "none", plot.margin = margin(3, 3, 2, 3))
}

pa <- plot_class("Glutamatergic"); pb <- plot_class("GABAergic")
leg <- get_legend(pa + theme(legend.position = "top", legend.justification = "center",
                             legend.text = element_text(size = BASE - 1),
                             legend.key.size = unit(8, "pt"),
                             legend.box.spacing = unit(2, "pt")))
grid <- plot_grid(pa, pb, ncol = 1, labels = c("a", "b"), label_size = 8,
                  label_fontface = "bold", hjust = -0.3, vjust = 1.3,
                  align = "v", axis = "lr")

# Only the platform colour key sits inside the figure — everything describing what
# the figure shows belongs in the manuscript figure legend, not baked into the PDF.
fig <- plot_grid(leg, grid, ncol = 1, rel_heights = c(0.035, 1))

ggsave(file.path(FIGDIR, "supertype_depth_platform.png"), fig, width = 7.1, height = 8.6,
       dpi = 400, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIGDIR, "supertype_depth_platform.pdf"), fig, width = 7.1, height = 8.6,
       bg = "white", device = pdf)
cat("saved", file.path(FIGDIR, "supertype_depth_platform.png/.pdf"), "\n")

# numbers for the figure legend, printed rather than drawn into the figure
cat("\nfor the legend:\n")
cat(sprintf("  %d neuronal supertypes, >= 50 cells on both platforms\n", nrow(sm)))
cat(sprintf("  median depth r = %.3f (Spearman %.3f), median |difference| = %.3f\n",
            r_med, cor(sm$median_MERFISH, sm$median_Xenium, method = "spearman"),
            median(abs(sm$delta_median))))
cat(sprintf("  Xenium IQR %.0f%% narrower (%.3f vs %.3f)\n",
            100 * (1 - median(sm$iqr_Xenium) / median(sm$iqr_MERFISH)),
            median(sm$iqr_Xenium), median(sm$iqr_MERFISH)))
big <- abs(sm$delta_median) > BIG_DELTA
cat(sprintf("  %d supertypes differ by > %.2f: %s\n", sum(big), BIG_DELTA,
            paste(sprintf("%s (%+.2f)", sm$supertype[big], sm$delta_median[big]),
                  collapse = ", ")))
