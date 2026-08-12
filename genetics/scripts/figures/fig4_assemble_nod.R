# Nine-panel variant of the Figure 4 assembly: the HCN1-expression-vs-depletion
# panel (old d) is dropped and its width handed to the morphology panel, which
# was the narrowest slot in its row (0.76 of 4.09) and could not fit three cells
# legibly.
#
# Panels are relabelled sequentially, so the old e..j become d..i.
#
#   a  genetics vs depletion     e  exemplar voltage traces  (was g)
#   b  Sst gene drivers          f  marker volcano           (was h)
#   c  HCN1 locus zoom           g  CALB1 violin             (was i)
#   d  HCN1 expr vs sag (was e)  h  SCZ vs AD concordance    (was j)
#   e  exemplar morphologies (was f)
suppressPackageStartupMessages({ library(cowplot) })

FIG4_WIDTH_IN  <- 8.0
FIG4_HEIGHT_IN <- 6.625 * (8.0 / 7.1)

build_figure4_nod <- function(panels) {
  need <- c("a", "b", "c", "e", "f", "g", "h", "i", "j")
  missing <- setdiff(need, names(panels))
  if (length(missing))
    stop("build_figure4_nod(): missing panel(s): ", paste(missing, collapse = ", "))

  # a and h keep a shared left edge; the old d is gone so it drops out here.
  col1 <- align_plots(panels$a, panels$h, align = "v", axis = "l")

  row_abc <- do.call(plot_grid,
    c(list(col1[[1]], panels$b, panels$c,
           nrow = 1, rel_widths = c(1, 1, 1),
           labels = c("a", "b", "c")),
      panel_label_kwargs))

  # Row total is unchanged at 4.09 so the figure keeps its proportions; the
  # 1.15 freed by dropping d goes almost entirely to the morphologies, taking
  # them from 0.76 to 1.80 (2.4x wider).
  row_efg <- do.call(plot_grid,
    c(list(panels$e, panels$f, panels$g,
           nrow = 1, rel_widths = c(1.14, 1.80, 1.15),
           align = "h", axis = "tb",
           labels = c("d", "e", "f")),
      panel_label_kwargs))

  row_hij <- do.call(plot_grid,
    c(list(col1[[2]], panels$i, panels$j,
           nrow = 1, rel_widths = c(1.25, 0.52, 1.10),
           align = "h", axis = "tb",
           labels = c("g", "h", "i")),
      panel_label_kwargs))

  plot_grid(row_abc, row_efg, row_hij, ncol = 1, rel_heights = c(1, 1, 1.06))
}
