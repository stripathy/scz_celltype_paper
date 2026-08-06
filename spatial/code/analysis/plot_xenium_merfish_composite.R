#!/usr/bin/env Rscript
# Xenium vs SEA-AD MERFISH concordance, summary and distribution, in one figure.
#   (a) subclass proportions            (b) neuronal supertype proportions
#   (c) subclass median cortical depth
#   (d) per-supertype depth distributions, Glutamatergic
#   (e) per-supertype depth distributions, GABAergic
#
# The supertype median-depth scatter that used to sit alongside (c) is dropped:
# panels d and e show the same comparison per supertype with the full
# distribution rather than one point, so the scatter added nothing.
#
# Panels d/e carry a caveat that a/b do not, and it belongs in the legend: MERFISH
# depth is SEA-AD's manual annotation and the Xenium depth model was trained on
# it, so agreement in position is expected by construction. The proportion panels
# have no such dependency — the classifier never saw MERFISH labels.
#
# No in-plot statistics anywhere; values are printed for the figure legend.
#
# Data: build_celltyping_validation_data.py   -> output/celltyping_supplement/data
#       build_supertype_depth_platform_data.py -> output/depth_platform
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(cowplot); library(ggrepel); library(scales)
})
cowplot::set_null_device("agg")

SD     <- "output/celltyping_supplement/data"
DD     <- "output/depth_platform"
FIGDIR <- "../manuscript/figures/supplementary"                  # the single home for submission figures
FIGSTEM <- "S03_xenium_merfish_concordance"

BASE <- 7                                      # one size across all five panels
CLASS_COLORS <- c(Glutamatergic = "#E65100", GABAergic = "#2E7D32", `Non-neuronal` = "#1565C0")
PLAT_COL <- c(MERFISH = "#5E3C99", Xenium = "#E66101")   # PuOr; not the Ctrl/SCZ blue-red
PT_SIZE <- 1.5; LAB_SIZE <- 1.6
LAYER_B <- c(.10, .40, .55, .70, .90)
LAYER_T <- c(.05, .25, .475, .625, .80, .95); LAYER_N <- c("L1","L2/3","L4","L5","L6","WM")
GLUT_ORDER <- c("L2/3 IT","L4 IT","L5 IT","L5 ET","L5/6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")
GABA_ORDER <- c("Lamp5","Sncg","Vip","Pax6","Chandelier","Pvalb","Sst","Sst Chodl")
SUBORDER <- list(Glutamatergic = GLUT_ORDER, GABAergic = GABA_ORDER)

# ---------------- (a-c) concordance scatters ----------------
theme_scatter <- function() {
  theme_classic(base_size = BASE) +
    theme(axis.title = element_text(size = BASE - 0.5),
          axis.text = element_text(size = BASE - 1.5, color = "black"),
          axis.line = element_line(linewidth = 0.25),
          axis.ticks = element_line(linewidth = 0.2),
          plot.margin = margin(3, 4, 2, 3), legend.position = "none", aspect.ratio = 1)
}
stat_lab <- function(r, rho, n, tag = "")
  sprintf("r = %.2f%s\nrho = %.2f\nn = %d", r, tag, rho, n)

prop_panel <- function(df, label_all = FALSE, n_label = 8) {
  df$lr  <- abs(log10(df$xenium_prop / df$merfish_prop))
  df$lab <- if (label_all) df$celltype else ifelse(rank(-df$lr) <= n_label, df$celltype, "")
  r   <- cor(log10(df$merfish_prop), log10(df$xenium_prop))
  rho <- cor(df$merfish_prop, df$xenium_prop, method = "spearman")
  lims <- range(c(df$merfish_prop, df$xenium_prop))
  ggplot(df, aes(merfish_prop, xenium_prop, fill = klass)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    geom_point(shape = 21, size = PT_SIZE, color = "white", stroke = 0.2, alpha = 0.95) +
    geom_text_repel(aes(label = lab), size = LAB_SIZE, color = "grey20", max.overlaps = 30,
                    min.segment.length = 0.1, segment.size = 0.15,
                    segment.color = "grey70", seed = 1) +
    annotate("text", x = lims[1], y = lims[2], hjust = 0, vjust = 1,
             label = stat_lab(r, rho, nrow(df), " (log)"), size = 1.7, lineheight = 0.95) +
    scale_x_log10(labels = label_percent(accuracy = 0.1)) +
    scale_y_log10(labels = label_percent(accuracy = 0.1)) +
    scale_fill_manual(values = CLASS_COLORS, name = NULL) +
    labs(x = "MERFISH proportion", y = "Xenium proportion") +
    theme_scatter()
}
depth_scatter <- function(df, label_all = TRUE, n_label = 8) {
  df$dev <- abs(df$xenium_depth - df$merfish_depth)
  df$lab <- if (label_all) df$celltype else ifelse(rank(-df$dev) <= n_label, df$celltype, "")
  r   <- cor(df$merfish_depth, df$xenium_depth)
  rho <- cor(df$merfish_depth, df$xenium_depth, method = "spearman")
  ggplot(df, aes(merfish_depth, xenium_depth, fill = klass)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    geom_point(shape = 21, size = PT_SIZE, color = "white", stroke = 0.2, alpha = 0.95) +
    geom_text_repel(aes(label = lab), size = LAB_SIZE, color = "grey20", max.overlaps = 30,
                    min.segment.length = 0.1, segment.size = 0.15,
                    segment.color = "grey70", seed = 1) +
    annotate("text", x = 0.97, y = 0.03, hjust = 1, vjust = 0,
             label = stat_lab(r, rho, nrow(df)), size = 1.7, lineheight = 0.95) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_fill_manual(values = CLASS_COLORS, name = NULL) +
    labs(x = "MERFISH median depth", y = "Xenium median depth") +
    theme_scatter()
}

ps <- read.csv(file.path(SD, "prop_subclass.csv"))
pt <- read.csv(file.path(SD, "prop_supertype.csv"))
ds <- read.csv(file.path(SD, "depth_subclass.csv"))
pa <- prop_panel(ps, label_all = TRUE)
pb <- prop_panel(pt[pt$klass != "Non-neuronal", ], label_all = FALSE, n_label = 8)
pc <- depth_scatter(ds, label_all = TRUE)
class_leg <- get_legend(pa + theme(legend.position = "bottom",
                                   legend.text = element_text(size = BASE - 1),
                                   legend.key.size = unit(8, "pt"),
                                   legend.box.spacing = unit(2, "pt")) +
                        guides(fill = guide_legend(override.aes = list(size = 2))))

# ---------------- (d, e) per-supertype depth distributions ----------------
df <- read_csv(file.path(DD, "supertype_depth_platform.csv.gz"), show_col_types = FALSE)
sm <- read_csv(file.path(DD, "supertype_depth_platform_summary.csv"), show_col_types = FALSE)
df$platform <- factor(df$platform, levels = c("MERFISH", "Xenium"))

shorten <- function(sup, sub) sub(paste0("^", sub, "[_ ]"), "", sup)
violin_class <- function(cls) {
  ro <- sm |> filter(class == cls) |>
    mutate(subclass = factor(subclass, levels = SUBORDER[[cls]])) |>
    arrange(subclass, median_MERFISH) |>
    mutate(xpos = row_number(), short = mapply(shorten, supertype, as.character(subclass)))
  d <- df |> filter(supertype %in% ro$supertype) |>
    left_join(ro |> select(supertype, xpos), by = "supertype")
  dc <- d |> group_by(supertype, platform) |>
    filter(depth >= quantile(depth, .025), depth <= quantile(depth, .975)) |> ungroup()
  bnds <- ro |> group_by(subclass) |> summarise(mx = max(xpos), mid = mean(xpos), .groups = "drop")
  ggplot(dc, aes(x = xpos, y = depth, fill = platform,
                 group = interaction(xpos, platform))) +
    geom_violin(position = position_dodge(0.85), width = 0.92, alpha = 0.9,
                linewidth = 0.07, color = "grey35", scale = "width") +
    stat_summary(fun = median, geom = "point", position = position_dodge(0.85),
                 size = 0.18, color = "white") +
    geom_hline(yintercept = LAYER_B, linetype = "dashed", color = "grey80", linewidth = 0.2) +
    geom_vline(xintercept = head(sort(bnds$mx), -1) + 0.5, color = "grey55", linewidth = 0.2) +
    annotate("text", x = bnds$mid, y = -0.025, label = as.character(bnds$subclass),
             fontface = "bold", size = 1.75, vjust = 0.5, hjust = 0, angle = 90) +
    scale_y_reverse(breaks = LAYER_T, labels = LAYER_N, expand = expansion(mult = c(0.02, 0))) +
    scale_x_continuous(breaks = ro$xpos, labels = ro$short, expand = expansion(add = 0.7)) +
    scale_fill_manual(values = PLAT_COL, name = NULL) +
    coord_cartesian(ylim = c(1.02, -0.215)) +
    labs(x = NULL, y = "cortical depth (pia at top, WM at bottom)", title = cls) +
    theme_cowplot(font_size = BASE) +
    theme(axis.text.x = element_text(size = BASE - 2.5, angle = 90, vjust = 0.5, hjust = 1,
                                     margin = margin(t = 0.5)),
          axis.text.y = element_text(size = BASE - 1),
          axis.title.y = element_text(size = BASE - 0.5),
          axis.line = element_line(linewidth = 0.25),
          axis.ticks = element_line(linewidth = 0.2), axis.ticks.length = unit(1.5, "pt"),
          plot.title = element_text(size = BASE + 0.5, face = "bold"),
          legend.position = "none", plot.margin = margin(3, 3, 2, 3))
}
pd_ <- violin_class("Glutamatergic"); pe <- violin_class("GABAergic")
plat_leg <- get_legend(pd_ + theme(legend.position = "top", legend.justification = "center",
                                   legend.text = element_text(size = BASE - 1),
                                   legend.key.size = unit(8, "pt"),
                                   legend.box.spacing = unit(2, "pt")))

# ---------------- assemble ----------------
kw <- list(label_size = 8, label_fontface = "bold", hjust = -0.3, vjust = 1.3)
row_scat <- do.call(plot_grid, c(list(pa, pb, pc, nrow = 1, labels = c("a", "b", "c")), kw))
row_d <- do.call(plot_grid, c(list(pd_, nrow = 1, labels = "d"), kw))
row_e <- do.call(plot_grid, c(list(pe,  nrow = 1, labels = "e"), kw))
fig <- plot_grid(row_scat, class_leg, plat_leg, row_d, row_e, ncol = 1,
                 rel_heights = c(0.72, 0.05, 0.04, 1, 1))

W <- 7.1; H <- 10.2
ggsave(file.path(FIGDIR, paste0(FIGSTEM, ".png")), fig, width = W, height = H,
       dpi = 400, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIGDIR, paste0(FIGSTEM, ".pdf")), fig, width = W, height = H,
       bg = "white", device = pdf)
cat(sprintf("saved %s/%s.png/.pdf (%.1f x %.1f in)\n", FIGDIR, FIGSTEM, W, H))

# ---- numbers for the figure legend ----
r_med <- cor(sm$median_MERFISH, sm$median_Xenium)
g <- sm |> group_by(subclass)
cen <- sm |> mutate(m = median_MERFISH - g$median_MERFISH[match(subclass, g$subclass)])
within <- sm |> group_by(subclass) |>
  mutate(m = median_MERFISH - mean(median_MERFISH), x = median_Xenium - mean(median_Xenium)) |>
  ungroup() |> filter(subclass %in% (sm |> count(subclass) |> filter(n >= 3))$subclass)
cat("\nfor the legend:\n")
cat(sprintf("  a subclass proportions   r(log10) = %.2f, rho = %.2f, n = %d\n",
            cor(log10(ps$merfish_prop), log10(ps$xenium_prop)),
            cor(ps$merfish_prop, ps$xenium_prop, method = "spearman"), nrow(ps)))
ptn <- pt[pt$klass != "Non-neuronal", ]
cat(sprintf("  b supertype proportions  r(log10) = %.2f, rho = %.2f, n = %d\n",
            cor(log10(ptn$merfish_prop), log10(ptn$xenium_prop)),
            cor(ptn$merfish_prop, ptn$xenium_prop, method = "spearman"), nrow(ptn)))
cat(sprintf("  c subclass median depth  r = %.2f, rho = %.2f, n = %d\n",
            cor(ds$merfish_depth, ds$xenium_depth),
            cor(ds$merfish_depth, ds$xenium_depth, method = "spearman"), nrow(ds)))
cat(sprintf("  d,e %d supertypes (%d glut / %d gaba); median depth r = %.2f, within-subclass r = %.2f\n",
            nrow(sm), sum(sm$class == "Glutamatergic"), sum(sm$class == "GABAergic"),
            r_med, cor(within$m, within$x)))
cat(sprintf("      median |difference| = %.3f; Xenium IQR %.0f%% narrower (%.3f vs %.3f)\n",
            median(abs(sm$delta_median)),
            100 * (1 - median(sm$iqr_Xenium) / median(sm$iqr_MERFISH)),
            median(sm$iqr_Xenium), median(sm$iqr_MERFISH)))
