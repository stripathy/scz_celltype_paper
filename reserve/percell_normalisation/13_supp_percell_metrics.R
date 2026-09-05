#!/usr/bin/env Rscript
# ============================================================================
# SUPPLEMENTAL — per-cell SST / PVALB reduction in SCZ across normalisations,
# on the canonical cells (corr_subclass + cortical + corr_qc, 24 donors). Four
# metrics per cell type:
#   1. Raw dots / cell           (count)
#   2. Dots / cell area          (grains/100 um^2; the Dienel 2023 per-neuron measure)
#   3. Library-size-normalised   (count / total counts, to median library)
#   4. Library size              (total counts; context)
# Boxplots = per-donor means (12 vs 12). p = negative-binomial mixed-effects
# model (donor random intercept; offset log cell area [col 2] / log library size
# [col 3]). This is the per-cell companion to the composite's edgeR panels.
#
# Input : results/tables/percell_grain_density.csv  (scripts/11_grain_density.py)
# Output: results/figures/S_percell_metrics.{png,pdf}
# ============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(ggsignif); library(cowplot); library(glmmTMB)
})
d <- read_csv("results/tables/percell_grain_density.csv", show_col_types = FALSE)
d$dx <- factor(d$dx, levels = c("Control", "SCZ"))
COL  <- c(Control = "#0072B2", SCZ = "#D55E00"); BASE <- 12
GENES <- c(SST = "SST in Sst", PVALB = "PVALB in Pvalb")

d <- d %>% group_by(gene) %>% mutate(
  m_raw     = count,
  m_area    = count / cell_area_um2 * 100,
  m_libnorm = count / total_counts * median(total_counts),
  m_lib     = total_counts) %>% ungroup()

# NB mixed model -> SCZ/Control rate ratio (exp of the dxSCZ coefficient) + p.
# Returns BOTH so every cited number is verifiable from a CSV (precision rule).
nb_fit <- function(df, offset = NULL, response = "count") {
  f <- if (is.null(offset)) sprintf("%s ~ dx + (1|sample)", response)
       else sprintf("%s ~ dx + offset(log(%s)) + (1|sample)", response, offset)
  co <- summary(glmmTMB(as.formula(f), family = nbinom2, data = df))$coefficients$cond["dxSCZ", ]
  c(ratio = unname(exp(co[1])), p = unname(co[4]))   # exp(Estimate), Pr(>|z|)
}
METRICS <- list(
  list(key = "m_raw",     lab = "Raw dots / cell",            fit = function(df) nb_fit(df)),
  list(key = "m_area",    lab = "Dots / cell area (grains/100 um2)", fit = function(df) nb_fit(df, offset = "cell_area_um2")),
  list(key = "m_libnorm", lab = "Lib-size-normalised",        fit = function(df) nb_fit(df, offset = "total_counts")),
  list(key = "m_lib",     lab = "Library size",               fit = function(df) nb_fit(df, response = "total_counts")))

# Stats table FIRST (one fit per gene x metric) -> CSV; panels look p up from it
# so the on-figure p and the cited p come from the same source of truth.
stats <- do.call(rbind, lapply(names(GENES), function(g)
  do.call(rbind, lapply(METRICS, function(m) {
    fr <- m$fit(d[d$gene == g, ])
    data.frame(gene = g, metric = m$key, label = m$lab,
               scz_over_ctrl = unname(fr["ratio"]), p = unname(fr["p"]))
  }))))
write_csv(stats, "results/tables/S_percell_stats.csv")
cat("Saved results/tables/S_percell_stats.csv\n"); print(stats, digits = 3)

pb <- d %>% group_by(gene, sample, dx) %>% summarise(across(starts_with("m_"), mean), .groups = "drop")
theme_panel <- function(show_x) theme_cowplot(font_size = BASE) +
  theme(legend.position = "none", axis.title = element_blank(),
        axis.text.x = if (show_x) element_text(size = BASE) else element_blank(),
        axis.ticks.x = if (show_x) element_line() else element_blank(),
        plot.margin = margin(3, 4, 3, 4))

panel <- function(g, m, show_x) {
  df <- pb[pb$gene == g, ]; df$y <- df[[m$key]]
  p <- stats$p[stats$gene == g & stats$metric == m$key]
  yr <- range(df$y)
  ggplot(df, aes(dx, y)) +
    geom_boxplot(aes(fill = dx), width = 0.55, outlier.shape = NA, alpha = 0.55, linewidth = 0.4) +
    geom_jitter(aes(fill = dx), shape = 21, colour = "grey25", stroke = 0.3, width = 0.13, height = 0,
                size = 1.6, alpha = 0.9) +
    geom_signif(comparisons = list(c("Control", "SCZ")), annotations = sprintf("italic(p)=='%.3f'", p),
                parse = TRUE, y_position = yr[2] + 0.05 * diff(yr), tip_length = 0.015, textsize = 3.0) +
    scale_fill_manual(values = COL) + scale_y_continuous(expand = expansion(mult = c(0.06, 0.20))) +
    theme_panel(show_x)
}
hdr  <- function(t) ggdraw() + draw_label(t, size = BASE - 1, lineheight = 0.85)
rlab <- function(t) ggdraw() + draw_label(t, size = BASE, angle = 90, fontface = "italic")
RW <- c(0.18, 1, 1, 1, 1)
header <- plot_grid(NULL, plotlist = lapply(METRICS, function(m) hdr(m$lab)), ncol = 5, rel_widths = RW)
rows <- lapply(names(GENES), function(g)
  plot_grid(rlab(GENES[g]), plotlist = lapply(METRICS, function(m) panel(g, m, show_x = (g == "PVALB"))),
            ncol = 5, rel_widths = RW))
body <- plot_grid(header, rows[[1]], rows[[2]], ncol = 1, rel_heights = c(0.16, 1, 1))
title <- ggdraw() + draw_label("Supplementary: per-cell SST / PVALB reduction in SCZ across normalisations (canonical cells)",
                               size = BASE, x = 0.5)
cap <- ggdraw() + draw_label("Points = per-donor means (12 vs 12). p: negative-binomial mixed model (donor random intercept). Col 2 = Dienel grain density.",
                             size = BASE - 4.5, x = 0.5, colour = "grey30")
final <- plot_grid(title, body, cap, ncol = 1, rel_heights = c(0.05, 1, 0.05))
ggsave("results/figures/S_percell_metrics.png", final, width = 11.5, height = 6.2, dpi = 200, bg = "white")
ggsave("results/figures/S_percell_metrics.pdf", final, width = 11.5, height = 6.2, bg = "white")
cat("Saved results/figures/S_percell_metrics.{png,pdf}\n")
