# ============================================================================
# _exemplar_panels.R — Figure 2 panels d and h: one representative Xenium cell
# per group, drawn from the committed exemplar_*.csv tables.
# ============================================================================
# Shared by 09_composite_figure.R (the figure) and any preview script, so the
# panel is defined once. Sourced, not run: the caller supplies BASE and TAB.
#
# What a reader has to be able to tell from the panel without the legend:
#   - the solid outline is the cell, the dashed one the nucleus,
#   - each red dot is one detected transcript of the marker gene,
#   - how many there are, and
#   - how big the cell is (5 um bar).
# The count label therefore says what it counts ("35 SST transcripts"), the
# first cell of the block names its two outlines (or a one-line key does), and
# the scale bar is labelled once at a size that can be read at print scale.
# Everything else is deliberately unlabelled so the two cells stay the subject.
#
# Inputs (written by scripts/10_xenium_exemplar_cells.py):
#   <TAB>/exemplar_<GENE>_<DX>_boundary.csv   cell polygon, microns, centroid at 0
#   <TAB>/exemplar_<GENE>_<DX>_nucleus.csv    nucleus polygon
#   <TAB>/exemplar_<GENE>_<DX>_dots.csv       marker transcripts inside the cell
# ============================================================================
suppressPackageStartupMessages({ library(readr); library(ggplot2); library(cowplot) })

DOT_COL  <- "#B2182B"   # marker-transcript dots; the same red carries the count label
SB_UM    <- 5           # scale-bar length, microns
TOP_BAND <- 0.26        # extra headroom above the cell, as a fraction of lim, for the
                        # one-line count label. Reserved for every cell so the label
                        # can never sit on a tall polygon (the Pvalb cells reach higher
                        # than the Sst cells) and all four viewports stay identical.

.ex_read <- function(what, gene, dx, tab)
  read_csv(sprintf("%s/exemplar_%s_%s_%s.csv", tab, gene, dx, what), show_col_types = FALSE)

# One coordinate limit for every cell so the zoom, and hence the physical length
# of the 5 um bar, is identical across all four panels. `expand` leaves room
# around the largest cell for outline labels.
exemplar_lim <- function(pairs, tab, expand = 1.05) {
  expand * max(vapply(pairs, function(p) {
    bd <- .ex_read("boundary", p[1], p[2], tab); max(abs(c(bd$x, bd$y))) }, numeric(1)))
}

# gene, dx        which exemplar
# lim             shared half-width of the square viewport (microns)
# scalebar_lab    draw "5 um" above this cell's bar (once per figure is enough)
# label_outlines  name the solid/dashed outlines with short leaders (first cell only)
# count_style     "word": "35 / SST transcripts";  "number": bare count (the old panel)
# `scale` is the caller's FIG_SCALE: the figure is drawn on a canvas wider than
# the print column, so every absolute size (pt text, mm linewidths) is multiplied
# by it. Data-space quantities -- lim, SB_UM, TOP_BAND -- must not be.
build_exemplar <- function(gene, dx, lim, tab, base,
                           scalebar_lab = FALSE, label_outlines = FALSE,
                           count_style = c("word", "number"), scale = 1) {
  count_style <- match.arg(count_style)
  bd <- .ex_read("boundary", gene, dx, tab)
  nu <- .ex_read("nucleus",  gene, dx, tab)
  dt <- .ex_read("dots",     gene, dx, tab)
  nd <- nrow(dt)
  top <- lim * (1 + TOP_BAND)

  # Text sizes follow the composite's hierarchy (09_composite_figure.R: panel
  # subtitle BASE > axis title BASE-0.5 > tick labels BASE-1, nothing above 7 pt),
  # converted from pt to ggplot's mm, so every label here sits on the same tier as
  # its counterpart in panels a-c. The count is the datum: axis-title tier, bold.
  mm <- function(pt) pt / .pt
  sz_count <- mm(base - 0.5 * scale)   # axis-title tier
  sz_bar   <- mm(base - 0.5 * scale)   # axis-title tier (was 0.26*BASE = 1.8 pt: illegible)
  sz_lab   <- mm(base - 1.0 * scale)   # tick-label tier

  p <- ggplot() +
    geom_polygon(data = bd, aes(x, y), fill = "grey93", colour = "grey40", linewidth = 0.35 * scale) +
    geom_polygon(data = nu, aes(x, y), fill = NA, colour = "grey25", linewidth = 0.3 * scale,
                 linetype = "22") +
    { if (nd > 0) geom_point(data = dt, aes(x, y), colour = DOT_COL, size = 0.55 * scale, alpha = 0.9) } +
    # 5 um bar, bottom right of every cell
    annotate("segment", x = lim * 0.96 - SB_UM, xend = lim * 0.96,
             y = -lim * 0.94, yend = -lim * 0.94, linewidth = 0.9 * scale, lineend = "square") +
    # Count, top left, in the reserved band: "35 SST transcripts" on one line. The
    # count is bold because it is the datum; the gene symbol is set upright, as
    # everywhere else in Figure 2; one text size keeps the line inside the panel.
    annotate("text", x = -lim * 0.98, y = top * 0.98, parse = TRUE, hjust = 0, vjust = 1,
             size = sz_count, colour = DOT_COL,
             label = if (count_style == "word") sprintf('bold("%d")~plain("%s transcripts")', nd, gene)
                     else sprintf('bold("%d")', nd)) +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, top)) +
    theme_void(base_size = base) +
    theme(plot.margin = margin(scale, 2 * scale, scale, 2 * scale))

  # plotmath "mu" rather than a Unicode micro sign: the composite's .pdf is written
  # by the base pdf() device, which cannot encode U+03BC/U+00B5 and drops the glyph;
  # plotmath draws mu through the symbol font on every device.
  if (scalebar_lab)
    p <- p + annotate("text", x = lim * 0.96 - SB_UM / 2, y = -lim * 0.94 + 0.06 * lim,
                      label = "5~mu*m", parse = TRUE, size = sz_bar, vjust = 0)

  # Name the outlines on one cell. "cell" leads from the rightmost boundary
  # vertex out to the right margin at mid-height; "nucleus" leads from the lowest
  # nucleus vertex down to the bottom-left corner, which the scale bar (bottom
  # right) and the count (top left) leave free.
  if (label_outlines) {
    br <- bd[which.max(bd$x), ]; nb <- nu[which.min(nu$y), ]
    lab_y <- -lim * 0.90
    p <- p +
      annotate("segment", x = br$x, xend = lim * 0.80, y = br$y, yend = br$y,
               linewidth = 0.3 * scale, colour = "grey30") +
      annotate("text", x = lim * 0.82, y = br$y, label = "cell", hjust = 0, vjust = 0.5,
               size = sz_lab, colour = "grey20") +
      annotate("segment", x = nb$x, xend = nb$x, y = nb$y, yend = lab_y + 0.10 * lim,
               linewidth = 0.3 * scale, colour = "grey30") +
      annotate("text", x = nb$x, y = lab_y, label = "nucleus", hjust = 0.5, vjust = 1,
               size = sz_lab, colour = "grey20")
  }
  p
}

# A one-line key: solid = cell, dashed = nucleus, dot = marker transcript. An
# alternative to label_outlines when the leaders make the first cell too busy;
# drawn once beneath panel d and shared with h.
exemplar_key <- function(gene, base, scale = 1) {
  sz <- (base - 1.0 * scale) / .pt   # tick-label tier, matching the outline labels
  ggplot() +
    annotate("segment", x = 0.02, xend = 0.09, y = 0.5, yend = 0.5, colour = "grey40", linewidth = 0.5 * scale) +
    annotate("text", x = 0.11, y = 0.5, label = "cell", hjust = 0, size = sz, colour = "grey20") +
    annotate("segment", x = 0.27, xend = 0.34, y = 0.5, yend = 0.5, colour = "grey25",
             linewidth = 0.5 * scale, linetype = "22") +
    annotate("text", x = 0.36, y = 0.5, label = "nucleus", hjust = 0, size = sz, colour = "grey20") +
    annotate("point", x = 0.60, y = 0.5, colour = DOT_COL, size = 0.9 * scale) +
    annotate("text", x = 0.63, y = 0.5, label = sprintf("%s transcript", gene),
             hjust = 0, size = sz, colour = "grey20") +
    xlim(0, 1) + ylim(0, 1) + theme_void()
}

# Control | SCZ pair for one marker.
#   show_hdr  draw the Control/SCZ strip (row 1); row 2 gets a blank strip of the
#             same height so its cells are the same size and the 5 um bar is the
#             same physical length in both rows.
#   key       "none" | "draw" (the one-line key under this pair) | "spacer" (a
#             blank of the key's height, for the row that does not carry it — same
#             reason as the header spacer).
ex_pair <- function(gene, lim, tab, base, show_hdr,
                    scalebar_lab_on = c("none", "Control", "SCZ"),
                    label_outlines_on = "none", key = c("none", "draw", "spacer"),
                    count_style = "word", scale = 1) {
  scalebar_lab_on <- match.arg(scalebar_lab_on); key <- match.arg(key)
  hdr  <- function(t) ggdraw() + draw_label(t, size = base)   # subtitle tier, like "SST / Sst" in b
  cell <- function(dx) build_exemplar(gene, dx, lim, tab, base,
                                      scalebar_lab   = identical(scalebar_lab_on, dx),
                                      label_outlines = identical(label_outlines_on, dx),
                                      count_style    = count_style, scale = scale)
  rows <- list(if (show_hdr) plot_grid(hdr("Control"), hdr("SCZ"), ncol = 2) else ggdraw(),
               plot_grid(cell("Control"), cell("SCZ"), ncol = 2))
  hts  <- c(0.16, 1)
  if (key != "none") { rows <- c(rows, list(if (key == "draw") exemplar_key(gene, base, scale) else ggdraw()))
                       hts  <- c(hts, 0.14) }
  plot_grid(plotlist = rows, ncol = 1, rel_heights = hts)
}
