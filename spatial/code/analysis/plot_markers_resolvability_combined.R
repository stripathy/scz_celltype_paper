#!/usr/bin/env Rscript
# Combined supplemental figure: marker genes + panel resolvability, 2 x 2.
#   row 1 = subclass        (a) marker dot plot   (b) panel-vs-transcriptome F1
#   row 2 = Sst supertype   (c) marker dot plot   (d) panel-vs-transcriptome F1
# Left column ~62% of width (dot plots carry many gene columns), right ~38%.
# Cell-type ordering is shared with the dot-plot long CSVs, so rows correspond.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(cowplot); library(scales)
})
cowplot::set_null_device("agg")
source("code/analysis/_xenium_exclusions.R")   # types too rare to analyse here (SM1)

DOT <- "output/depth_validation/lieber_layers"
RES <- "output/celltyping_supplement/data"
OUT <- "../manuscript/figures/supplementary"                     # the single home for submission figures
FIGSTEM <- "S02_xenium_celltype_annotation"

BASE <- 8.0
REDS <- c("#f7f7f7", "#fdd0a2", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15")
PANEL_COL <- "#3182BD"; TX_COL <- "#252525"; SEG <- "#C7C7C7"
HI <- "Sst_25"; HI_COL <- "#B8860B"
SETLAB <- c(panel = "Xenium panel (300 genes)", all = "full transcriptome")

# ---------------- Kwon/Lieber label agreement heatmap ----------------
# Ported from plot_lieber_celltype_diagonal.py so the panel shares this
# figure's fonts and sizing; the ordering logic is reproduced exactly.
KWON_ORDER <- c("L2/3 Ex", "L4/5 Ex", "L5 Ex", "L6 Ex", "MGE", "CGE",
                "Ast", "Oligo", "Ambig/In/Endo", "Mic", "Endo", "Ambig/Oligo", "NA")
KWON_GROUPS <- list(Excitatory = c("L2/3 Ex", "L4/5 Ex", "L5 Ex", "L6 Ex"),
                    Inhibitory = c("MGE", "CGE"),
                    `Glia / vascular` = c("Ast", "Oligo", "Ambig/In/Endo", "Mic", "Endo"),
                    Ambiguous = c("Ambig/Oligo", "NA"))
# ColorBrewer Greens-9, deliberately NOT the Reds used by the marker dot plots:
# within one figure a colour must carry one meaning, and red already encodes
# scaled mean expression in panels b/d while blue encodes "Xenium panel" in
# panels c/e. Green is unused elsewhere here, so it cannot be misread.
# Lowest stop is pure white so a zero cell is indistinguishable from the page.
HEAT_COLS <- c("#ffffff", "#e5f5e0", "#c7e9c0", "#a1d99b", "#74c476",
               "#41ab5d", "#238b45", "#006d2c", "#00441b")

build_lieber_heatmap <- function() {
  raw <- read.csv(file.path(DOT, "subclass_vs_celltype.csv"),
                  row.names = 1, check.names = FALSE)
  raw <- raw[rownames(raw) != "Unassigned", , drop = FALSE]
  raw <- raw[, colSums(raw, na.rm = TRUE) > 0, drop = FALSE]
  raw <- raw[, KWON_ORDER[KWON_ORDER %in% colnames(raw)], drop = FALSE]
  m <- as.matrix(raw)
  pct <- 100 * m / rowSums(m)                       # % of each of our subclasses
  pct <- pct[!rownames(pct) %in% XENIUM_TOO_RARE_SUBCLASS, , drop = FALSE]
  if (any(rownames(as.matrix(raw)) %in% XENIUM_TOO_RARE_SUBCLASS))
    message("  excluded as too rare in Xenium (SM1): ",
            paste(intersect(rownames(raw), XENIUM_TOO_RARE_SUBCLASS), collapse = ", "))

  # order our subclasses by dominant Kwon type, then by share -> block diagonal
  dom <- max.col(pct, ties.method = "first")
  ord <- order(dom, -apply(pct, 1, max))
  pct <- pct[ord, , drop = FALSE]; dom <- dom[ord]

  d <- as.data.frame(as.table(pct))
  names(d) <- c("subclass", "kwon", "value")
  d$subclass <- factor(d$subclass, levels = rownames(pct))
  d$kwon <- factor(d$kwon, levels = rev(colnames(pct)))   # first group at top

  # separators: between Kwon groups (horizontal), and where the dominant Kwon
  # group of our subclasses changes (vertical)
  grp_of <- setNames(rep(names(KWON_GROUPS), lengths(KWON_GROUPS)), unlist(KWON_GROUPS))
  present <- colnames(pct)
  hb <- cumsum(sapply(KWON_GROUPS, function(g) sum(present %in% g)))
  hb <- hb[hb > 0 & hb < length(present)]
  vb <- which(diff(match(grp_of[present[dom]], names(KWON_GROUPS))) != 0) + 0.5

  ggplot(d, aes(subclass, kwon, fill = value)) +
    geom_tile(colour = NA) +
    geom_text(data = subset(d, value >= 10),
              # 80: greens are darker in luminance than greys at the same value,
              # so black text stays the better choice further up the ramp than it
              # would on a grey scale, and white only wins from ~80% on
              aes(label = sprintf("%.0f", value),
                  colour = value >= 80), size = (BASE - 2.5) * 0.35, show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white")) +
    { if (length(hb)) geom_hline(yintercept = length(present) - hb + 0.5,
                                 colour = "#444", linewidth = 0.3) } +
    { if (length(vb)) geom_vline(xintercept = vb, colour = "#444", linewidth = 0.3) } +
    scale_fill_gradientn(colours = HEAT_COLS, limits = c(0, 100),
                         breaks = c(0, 50, 100), name = "% of subclass",
                         guide = guide_colourbar(barheight = unit(3.5, "pt"),
                                                 barwidth = unit(46, "pt"),
                                                 title.position = "top",
                                                 title.hjust = 0.5)) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Reannotated SEA-AD subclass", y = "Kwon et al. cell type") +
    theme_cowplot(font_size = BASE) +
    # axis text matched to the dot plots' y labels (BASE - 1) so the three
    # panel types read at one size; the colourbar sits in the empty
    # excitatory x non-neuronal corner rather than costing a row of height
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = BASE - 1),
          axis.text.y = element_text(size = BASE - 1),
          axis.title.y = element_text(size = BASE - 0.5),
          # the 45-degree labels leave a tall bounding box, which pushes the x
          # title so far down it reads as belonging to the panel below; pull it
          # back up against the labels
          axis.title.x = element_text(size = BASE - 0.5, margin = margin(t = -7)),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.25),
          panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.4),
          legend.position = c(0.865, 0.87),
          legend.justification = c(0.5, 0.5),
          legend.direction = "horizontal",
          legend.background = element_rect(fill = "white", colour = NA),
          legend.margin = margin(1, 3, 1, 3),
          legend.title = element_text(size = BASE - 1.5),
          legend.text = element_text(size = BASE - 2),
          plot.margin = margin(3, 3, 0, 3))
}

# ---------------- marker dot plot ----------------
dot_panel <- function(csv, sec_levels, sec_relabel = NULL, hi = NULL) {
  d <- read_csv(file.path(DOT, csv), show_col_types = FALSE)
  if (!is.null(sec_relabel)) d$section <- dplyr::recode(d$section, !!!sec_relabel)
  d$section <- factor(d$section, levels = sec_levels)
  d$ylab <- factor(d$cell_type, levels = unique(d$cell_type[order(d$ct_idx)]))
  d$gene <- factor(d$gene, levels = unique(d$gene[order(d$gene_idx)]))
  # this column carries the shared y labels for its row, so the Sst_25 callout
  # lives here rather than on the resolvability panel
  ylv <- rev(levels(d$ylab))
  y_el <- if (is.null(hi)) element_text(size = BASE - 1) else
    element_text(size = BASE - 1, colour = ifelse(ylv == hi, HI_COL, "black"),
                 face = ifelse(ylv == hi, "bold", "plain"))
  ggplot(d, aes(gene, ylab)) +
    geom_point(aes(size = pct, colour = scaled_mean)) +
    facet_grid(~ section, scales = "free_x", space = "free_x") +
    scale_colour_gradientn(colours = REDS, limits = c(0, 1),
                           name = "scaled mean expression",
                           breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                           guide = guide_colourbar(barheight = unit(3.5, "pt"),
                                                   barwidth = unit(50, "pt"),
                                                   title.vjust = 1)) +
    scale_size_area(max_size = 2.1, limits = c(0, 1), breaks = c(.1, .25, .5, .9),
                    labels = percent_format(accuracy = 1), name = "% expressing") +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = NULL) +
    theme_cowplot(font_size = BASE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     face = "italic", size = BASE - 2),
          axis.text.y = y_el,
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.25),
          strip.clip = "off",
          strip.background = element_rect(fill = "grey90", colour = NA),
          strip.text = element_text(size = BASE - 1.5, margin = margin(1.5, 1, 1.5, 1)),
          panel.spacing = unit(2, "pt"),
          panel.border = element_rect(colour = "grey85", fill = NA, linewidth = 0.3),
          legend.position = "bottom", legend.box = "horizontal",
          legend.justification = "center",
          legend.title = element_text(size = BASE - 1),
          legend.text = element_text(size = BASE - 1.5),
          plot.margin = margin(3, 3, 2, 3))
}

# ---------------- resolvability dumbbell ----------------
# y labels come from the dot plot in the same row, so this panel draws none;
# that requires the two panels to carry exactly the same cell types.
res_panel <- function(d, order_vec) {
  ord <- rev(order_vec)
  # anything still absent after the too-rare exclusion is a genuine mismatch
  # between the two data sources and should be reported, not quietly dropped
  dropped <- setdiff(d$cell_type, ord)
  if (length(dropped))
    message("  note: not shown (absent from the marker dot plot): ",
            paste(dropped, collapse = ", "))
  d <- d[d$cell_type %in% ord, ]
  d$cell_type <- factor(d$cell_type, levels = ord)
  L <- d |> select(cell_type, F1_panel, F1_all) |>
    pivot_longer(c(F1_panel, F1_all), names_to = "set", values_to = "F1") |>
    mutate(set = factor(SETLAB[sub("F1_", "", set)], levels = SETLAB))
  g <- ggplot() +
    geom_segment(data = d, aes(y = cell_type, yend = cell_type, x = F1_panel, xend = F1_all),
                 color = SEG, linewidth = 0.8) +
    geom_point(data = L, aes(x = F1, y = cell_type, fill = set),
               shape = 21, color = "white", stroke = 0.3, size = 1.7) +
    scale_fill_manual(values = setNames(c(PANEL_COL, TX_COL), SETLAB), name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5),
                       expand = expansion(mult = c(0.03, 0.05))) +
    labs(x = "classification F1", y = NULL) +
    theme_cowplot(font_size = BASE) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = BASE - 1.5),
          axis.title.x = element_text(size = BASE - 0.5),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.25),
          panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
          legend.position = "bottom",
          legend.text = element_text(size = BASE - 1.5),
          plot.margin = margin(3, 2, 2, 1))
  g
}

# ---------------- data ----------------
# "pan-Inh" is wider than its 3-gene facet at this column width and collides
# with the axis line; "GABA" fits and matches the row-2 "GABA/SST" strip
pa <- dot_panel("subclass_dotplot_long.csv",
                c("GABA", "Inhibitory", "Excitatory", "Non-neuronal"),
                sec_relabel = c(`pan-Inh` = "GABA"))
pc <- dot_panel("sst_dotplot_long.csv", c("GABA/SST", "Sst"),
                sec_relabel = c(anchor = "GABA/SST"), hi = HI)

sc_order <- read_csv(file.path(DOT, "subclass_dotplot_long.csv"), show_col_types = FALSE) |>
  distinct(cell_type, ct_idx) |> arrange(ct_idx) |> pull(cell_type)
sst_order <- read_csv(file.path(DOT, "sst_dotplot_long.csv"), show_col_types = FALSE) |>
  distinct(cell_type, ct_idx) |> arrange(ct_idx) |> pull(cell_type)

sc <- read_csv(file.path(RES, "resolvability_subclass.csv"), show_col_types = FALSE) |>
  drop_too_rare("cell_type", "subclass")
ss <- read_csv(file.path(RES, "resolvability_sst_supertype.csv"), show_col_types = FALSE) |>
  drop_too_rare("cell_type", "supertype")
pb <- res_panel(sc, sc_order)
pd <- res_panel(ss, sst_order)

# ---------------- assemble ----------------
ph <- build_lieber_heatmap()

dot_leg  <- get_legend(pa)
res_leg  <- get_legend(pb)
kw <- list(label_size = 8, label_fontface = "bold", label_x = 0, label_y = 1,
           hjust = -0.3, vjust = 1.3)
strip <- function(p) p + theme(legend.position = "none")

# align="h"/axis="tb" makes the two plotting regions in a row share top and
# bottom edges, so the dot plot's y labels index the dumbbell rows correctly
r2 <- align_plots(strip(pa), strip(pb), align = "h", axis = "tb")
r3 <- align_plots(strip(pc), strip(pd), align = "h", axis = "tb")
row1 <- do.call(plot_grid, c(list(ph, nrow = 1, labels = "a"), kw))
row2 <- do.call(plot_grid, c(r2, list(nrow = 1, rel_widths = c(0.66, 0.34),
                                      labels = c("b", "c")), kw))
row3 <- do.call(plot_grid, c(r3, list(nrow = 1, rel_widths = c(0.66, 0.34),
                                      labels = c("d", "e")), kw))
legs <- plot_grid(dot_leg, res_leg, nrow = 1, rel_widths = c(0.66, 0.34))
fig  <- plot_grid(row1, row2, row3, legs, ncol = 1,
                  rel_heights = c(1.22, 1.30, 1, 0.11))

W <- 7.1; H <- 9.9
ggsave(file.path(OUT, paste0(FIGSTEM, ".png")), fig,
       width = W, height = H, dpi = 400, bg = "white", device = ragg::agg_png)
ggsave(file.path(OUT, paste0(FIGSTEM, ".pdf")), fig,
       width = W, height = H, bg = "white", device = pdf)
cat("saved", file.path(OUT, paste0(FIGSTEM, ".png/.pdf")),
    sprintf("(%.1f x %.1f in)\n", W, H))
