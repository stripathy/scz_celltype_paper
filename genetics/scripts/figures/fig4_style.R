# Shared Figure 4 aesthetic vocabulary -- the single source of truth.
#
# Sourced by scz_sst_hcn1_story.R (panels a-g), fig4_new_panels.R (panels h-j)
# and render_fig4_new_panels.R, so every panel is styled from one place. This
# file used to be a hand-maintained mirror of an identical block inside
# scz_sst_hcn1_story.R; the copies drifted, and a font change applied to one of
# them silently left the other seven panels at the old size. Change sizes,
# colours and helpers here only.

suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })

# Figure 2 is drawn at 7.1 in wide with 7 pt base text. This figure is drawn on
# a larger 8.0 in canvas -- which gives the 10 panels more room for labels during
# layout -- and every text size is Figure 2's value multiplied by 8.0/7.1. Once
# the figure is scaled to a 7.1 in column, its text matches Figure 2 exactly.
FIG_SCALE        <- 8.0 / 7.1
# Base was 7 pt at print scale, matching Figure 2 (09_composite_figure.R).
# Raised to 8 pt: Figure 4 packs 10 panels into the same canvas, so its axes and
# insets read smaller than Figure 2's at the same nominal size. Resulting print
# sizes: axis titles 7.5 pt, tick labels 7.0 pt (Nature guidance is 5-7 pt).
BASE_FONT_SIZE   <- 8 * FIG_SCALE
PANEL_LABEL_SIZE <- 9 * FIG_SCALE

# HCN1 is blue everywhere it is called out (gene-driver scatter, locus track,
# and the marker volcano), so the reader tracks one gene by colour across panels.
HCN1_COLOR  <- "#1565C0"
# CALB1 takes Sst_25's own SEA-AD colour: it is the marker of that population.
CALB1_COLOR <- "#693d07"

# Binary depletion status is carried by the point OUTLINE (colour + thickness),
# leaving fill free for the SEA-AD supertype identity colour.
DEPLETION_OUTLINE <- c(`Depleted (FDR < 0.20)` = "black",
                       `Not depleted`          = "grey60")
DEPLETION_STROKE  <- c(`Depleted (FDR < 0.20)` = 0.7,
                       `Not depleted`          = 0.2)

# Direction colours for the depleted-vs-not-depleted volcano, taken from the
# Sst family palette so the panel reads as part of the Sst story.
COL_UP_VULN <- "#7b4c10"   # darker Sst tone  (higher in the depleted types)
COL_UP_NOTD <- "#f2ad49"   # lighter Sst tone (higher in the not-depleted types)
COL_NS      <- "grey80"
# In-panel text sizes, expressed as multiples of BASE exactly as Figure 2 does
# (09_composite_figure.R): gene/point labels BASE*0.32, bold callouts BASE*0.34,
# secondary labels BASE*0.28, inset statistics BASE*0.26.
LBL_GENE  <- BASE_FONT_SIZE * 0.32   # = 2.24 at 7.1 in print scale
LBL_CALL  <- BASE_FONT_SIZE * 0.34   # = 2.38 at print scale
LBL_SMALL <- BASE_FONT_SIZE * 0.28   # = 1.96 at print scale
# Inset statistics were the smallest text in the figure at BASE*0.26 (5.2 pt at
# print scale); raised to match the gene/callout labels.
LBL_STAT  <- BASE_FONT_SIZE * 0.32


# Exact copy of theme_panel() in scz_sst_hcn1_story.R -- keep byte-identical so
# panels built here are indistinguishable from the existing ones.
theme_panel <- function(base_size = BASE_FONT_SIZE) {
  theme_cowplot(font_size = base_size) +
    theme(
      plot.title        = element_blank(),
      plot.subtitle     = element_blank(),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.title        = element_text(size = base_size - 0.5 * FIG_SCALE),
      axis.text         = element_text(size = base_size - 1.0 * FIG_SCALE),
      axis.line         = element_line(linewidth = 0.25 * FIG_SCALE, color = "grey20"),
      axis.ticks        = element_line(linewidth = 0.25 * FIG_SCALE, color = "grey20"),
      legend.title      = element_text(size = base_size - 1.0 * FIG_SCALE),
      legend.text       = element_text(size = base_size - 1.5 * FIG_SCALE),
      legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
      legend.key.size   = unit(9 * FIG_SCALE, "pt"),
      plot.margin       = margin(3, 4, 3, 4))
}

# Plain-language axis labels, shared across panels for a consistent vocabulary
# (mirrors scz_sst_hcn1_story.R).
LAB_GWAS <- expression("SCZ GWAS enrichment ("*-log[10]~italic(P)*")")
LAB_DEPL <- expression("Cell depletion in SCZ ("*-beta*")")
LAB_HCN1 <- expression(italic("HCN1")*" expression ("*log[2]*" CP10K+1)")

# Shared with scz_sst_hcn1_story.R: cell depletion (-beta) axis range.
DEPLETION_LIM <- c(-0.13, 0.34)

panel_label_kwargs <- list(
  label_size = PANEL_LABEL_SIZE, label_fontface = "bold",
  label_x = 0.0, label_y = 1.0, hjust = -0.3, vjust = 1.3)

add_depletion_status <- function(df) {
  df |> dplyr::mutate(
    depleted_status = factor(
      ifelse(depleted_fdr20, "Depleted (FDR < 0.20)", "Not depleted"),
      levels = c("Depleted (FDR < 0.20)", "Not depleted")))
}

depletion_scales <- function(show_legend = TRUE) {
  color_guide <- if (show_legend) {
    guide_legend(override.aes = list(shape = 21, fill = "grey85", size = 2.5,
                                     stroke = unname(DEPLETION_STROKE)))
  } else "none"
  list(
    scale_color_manual(values = DEPLETION_OUTLINE,
                       breaks = names(DEPLETION_OUTLINE),
                       name = NULL, guide = color_guide),
    scale_discrete_manual("stroke", values = DEPLETION_STROKE, guide = "none"))
}

spearman_rp <- function(x, y) {
  list(rho = cor(x, y, method = "spearman"),
       p   = suppressWarnings(cor.test(x, y, method = "spearman")$p.value))
}

inset_spearman <- function(rho, p,
                           corner = c("top-right", "bottom-right",
                                      "top-left", "bottom-left"),
                           size = LBL_STAT) {
  corner <- match.arg(corner)
  xval  <- if (grepl("right", corner)) Inf  else -Inf
  yval  <- if (grepl("top",   corner)) Inf  else -Inf
  hjust <- if (grepl("right", corner)) 1.05 else -0.05
  vjust <- if (grepl("top",   corner)) 1.5  else -1.0
  # R's exact Spearman test (AS 89) underflows to p == 0 in the extreme tail --
  # panel j (rho = 0.87, n = 16) hits this, where the true value is ~1e-5
  # (asymptotic 1.3e-5; 2e6-draw permutation 3.5e-5). Never print a false zero:
  # below 1e-4 fall back to a threshold that holds under exact, asymptotic and
  # permutation alike.
  lab <- if (!is.finite(p) || p < 1e-4) {
    sprintf("rho == %.2f * ',' ~ italic(p) < 10^-4", rho)
  } else {
    sprintf("rho == %.2f * ',' ~ italic(p) == %.3g", rho, p)
  }
  annotate("text", x = xval, y = yval, hjust = hjust, vjust = vjust,
           label = lab, parse = TRUE, size = size, color = "#222")
}

# The standard supertype-scatter point + label geoms used by panels B/E.
supertype_points <- function(size = 2.8) {
  geom_point(aes(fill = color, color = depleted_status, stroke = depleted_status),
             shape = 21, size = size, alpha = 0.92)
}

supertype_labels <- function(size = LBL_GENE) {
  ggrepel::geom_text_repel(aes(label = supertype), size = size,
                           # padding/force raised with the larger base font so
                           # neighbouring supertype labels stay separated
                           box.padding = 0.32, point.padding = 0.20,
                           max.overlaps = 30, segment.color = "#999",
                           segment.size = 0.2, force = 4.5,
                           min.segment.length = 0, seed = 1)
}

trend_line <- function() {
  geom_smooth(method = "lm", formula = y ~ x, color = "#444", fill = "#888",
              alpha = 0.18, linewidth = 0.3, linetype = "dashed", se = TRUE)
}
