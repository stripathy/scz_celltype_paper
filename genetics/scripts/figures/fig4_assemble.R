# Assembly of Figure 4 from its ten finished panels.
#
# Kept separate from the panel builders so that "what the figure looks like"
# (row grouping, relative widths, canvas size) lives in one place, and so the
# whole figure can be produced by a single call:
#
#   fig <- build_figure4(panels)   # panels = named list a..j of ggplot objects
#   render_figure4(fig, stem)
#
# Requires fig4_style.R (panel_label_kwargs) and cowplot.

suppressPackageStartupMessages({ library(cowplot) })

FIG4_PANEL_IDS <- letters[1:10]   # a..j

# Canvas: same aspect as the Figure 2 composite (7.1 x 6.625 in) so the two
# figures share text sizes once both are scaled to a 7.1 in column. Drawn wider
# (8.0 in) to give ten panels room during layout; fig4_style.R scales every font
# by 8.0/7.1 to compensate.
FIG4_WIDTH_IN  <- 8.0
FIG4_HEIGHT_IN <- 6.625 * (8.0 / 7.1)

#' Assemble the ten Figure 4 panels into the final multipanel.
#'
#' @param panels Named list of ggplot objects with names "a" through "j":
#'   a  genetics vs depletion      f  exemplar morphologies
#'   b  Sst_25 gene drivers        g  exemplar voltage traces
#'   c  HCN1 locus zoom            h  depleted vs not-depleted marker volcano
#'   d  HCN1 expression vs depl.   i  CALB1 violin
#'   e  HCN1 expression vs sag     j  SCZ vs AD (DLPFC) concordance
#' @return A cowplot object.
build_figure4 <- function(panels) {
  missing <- setdiff(FIG4_PANEL_IDS, names(panels))
  if (length(missing))
    stop("build_figure4(): missing panel(s): ", paste(missing, collapse = ", "))

  # Vertical alignment across rows: pad the three first-column panels so a, d
  # and h share a left edge (their y-axis titles/tick labels line up down the
  # figure). Equal *widths* are not possible without discarding the per-row
  # width tuning -- a occupies 33% of its row, d 26%, h 44% -- so only the left
  # edge is matched.
  col1 <- align_plots(panels$a, panels$d, panels$h, align = "v", axis = "l")

  row_abc <- do.call(plot_grid,
    c(list(col1[[1]], panels$b, panels$c,
           nrow = 1, rel_widths = c(1, 1, 1),
           labels = c("a", "b", "c")),
      panel_label_kwargs))

  row_defg <- do.call(plot_grid,
    c(list(col1[[2]], panels$e, panels$f, panels$g,
           nrow = 1, rel_widths = c(1.15, 1.15, 0.76, 1.03),
           align = "h", axis = "tb",
           labels = c("d", "e", "f", "g")),
      panel_label_kwargs))

  # Unequal widths by information density: the volcano carries ~20 gene labels
  # and the AD scatter 16 supertype labels, whereas the violin is two categories
  # and would only gain whitespace from an equal share.
  row_hij <- do.call(plot_grid,
    c(list(col1[[3]], panels$i, panels$j,
           nrow = 1, rel_widths = c(1.25, 0.52, 1.10),
           align = "h", axis = "tb",
           labels = c("h", "i", "j")),
      panel_label_kwargs))

  plot_grid(row_abc, row_defg, row_hij, ncol = 1, rel_heights = c(1, 1, 1.06))
}

#' Write Figure 4 to PNG + PDF + SVG.
#'
#' @param fig    Object returned by build_figure4().
#' @param stem   Output path WITHOUT an extension.
render_figure4 <- function(fig, stem,
                           width = FIG4_WIDTH_IN, height = FIG4_HEIGHT_IN) {
  save_figure(fig, stem, width, height)
}
