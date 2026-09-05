#!/usr/bin/env Rscript
# ===================================================================
# scz_sst_hcn1_story.R
#
# Renders Figure 4: the 9-panel SCZ-genetics / HCN1 / Sst figure.
# Data are pre-extracted into results/figures/r_panels/ by export_panels_abc.py,
# export_panel_d_genetrack.py, export_panels_dehi.py, export_panels_ef.py and
# export_panel_ad_concordance.py. This script reads those CSVs, builds each
# panel with ggplot2, and hands them to build_figure4_nod() in
# fig4_assemble_nod.R for the final layout.
#
# Shared style (fonts, colours, theme_panel, inset helpers) comes from
# fig4_style.R, which also styles panels f-h in fig4_new_panels.R.
#
# Panel guide (top->bottom, left->right):
#   a  genetics_vs_depletion_c  SCZ common-variant enrichment x depletion
#   b  gene_driver_plot         Sst_2 gene drivers (specificity x MAGMA p)
#   c  hcn1_locus_plot          HCN1 locus zoom + fine-mapping + gene track
#   d  sag_c                    HCN1 expression x patch-seq voltage sag
#   e  morphology_plot          Five exemplar Sst reconstructions
#   f  ephys_traces_plot        Voltage responses for the same five cells
#   g  marker volcano           Depleted vs not-depleted Sst markers
#   h  CALB1 violin             CALB1 by depletion group
#   i  AD concordance           SCZ depletion x SEA-AD DLPFC CPS slope
#
# The builders still pass panels under their historical letters a,b,c,e,f,g,h,i,j;
# build_figure4_nod() drops the old d (HCN1 expression x depletion) and relabels
# the rest sequentially to a-i. Convergence across a and d (both on the same 16
# Sst supertypes, sharing the SEA-AD fill + depletion-outline encoding):
# genetics <-> depletion, HCN1 expression <-> intrinsic sag.
# ===================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(glue)
  library(ggplot2); library(cowplot); library(ggrepel); library(scales)
})

# ────────────────────────────────────────────────────────────────────
# Paths
# ────────────────────────────────────────────────────────────────────
# Repo-relative: resolve this script's own location, so the render runs from
# any clone. Rscript exposes it via --file=; fall back to the cwd convention.
# Repo-relative: resolve this script's own location, so it runs from any clone.
# Rscript exposes the path via --file=; the fallback assumes the repo root cwd.
.script_dir <- function() {
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else NA_character_
}
.sd <- .script_dir()
REPO     <- normalizePath(if (!is.na(.sd)) file.path(.sd, "..", "..") else "genetics")
DATA_DIR <- file.path(REPO, "results", "figures", "r_panels")
FIGDIR   <- file.path(REPO, "results", "figures")


# ────────────────────────────────────────────────────────────────────
# Output helper — every figure is written as PNG (raster preview),
# PDF (vector, for submission) and SVG (vector, for hand-editing in
# Illustrator/Inkscape). `stem` is a path WITHOUT an extension.
# svglite is used for SVG: it needs no X11/cairo and embeds text as
# real text, so labels stay editable.
# ────────────────────────────────────────────────────────────────────
FIG_DPI <- 400

save_figure <- function(plot, stem, width, height, dpi = FIG_DPI) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
         dpi = dpi, bg = "white")
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
         bg = "white")
  ggsave(paste0(stem, ".svg"), plot, width = width, height = height,
         bg = "white", device = svglite::svglite)
  invisible(paste0(stem, c(".png", ".pdf", ".svg")))
}

# The panel CSVs are snapshots of upstream analyses, so they can go stale
# without any error. r_panels/MANIFEST.tsv records the source each one came
# from; this stops the render if any source has been rerun since.
source(file.path(dirname(REPO), "shared", "figure_inputs.R"))
fi_check(DATA_DIR,
         refresh_cmd = paste("re-run the export_panels_*.py exporters for the",
                             "stale panels, then re-run this script"))

read_panel <- function(filename) {
  read_csv(file.path(DATA_DIR, filename), show_col_types = FALSE)
}

# ────────────────────────────────────────────────────────────────────
# Visual constants
# ────────────────────────────────────────────────────────────────────
# Sized for a Nature Neuroscience-style double-column figure: 6.5" wide.
# Text scales proportionally LESS than spatial elements so it stays
# legible at print scale (~7 pt body, ~12 pt panel labels).
# Figure 2 is drawn at 7.1 in wide with 7 pt base text. This figure is drawn on
# a larger 8.0 in canvas -- which gives the 10 panels more room for labels during
# layout -- and every text size is Figure 2's value multiplied by 8.0/7.1. Once
# the figure is scaled to a 7.1 in column, its text matches Figure 2 exactly.
# Shared aesthetic vocabulary -- font sizes, depletion encoding, theme_panel(),
# inset_spearman() and friends -- lives in fig4_style.R and is sourced by both
# this script (panels a-g) and fig4_new_panels.R (panels h-j). The two files
# previously carried duplicate copies, so a font change applied to one of them
# silently left the other seven panels at the old size.
source(file.path(REPO, "scripts", "figures", "fig4_style.R"))

# ---- constants specific to this script, not part of the shared vocabulary ----
FIG_WIDTH_IN  <- 6.5
FIG_HEIGHT_IN <- 6.36   # extra 0.75 in over the aspect-preserving 5.61; the
                        # extra height is routed to rows 2/3 so the middle and
                        # bottom panels are less vertically crowded.
HCN1_COLOR    <- CALB1_COLOR  # was #1565C0; blue tied HCN1 to nothing else in
                             # the figure, so it now shares CALB1's Sst_25 tone
LAB_HCN1      <- expression(italic("HCN1")*" expression ("*log[2]*" CP10K+1)")
# Shared axis range so d-x and e-x line up (DEPLETION_LIM comes from fig4_style).
HCN1_LIM      <- c(1.35, 3.65)


# ====================================================================
# B) Genetic risk vs SCZ cell-abundance depletion (16 Sst supertypes)
# ====================================================================
build_genetics_vs_depletion <- function(compact = FALSE) {
  df <- read_panel("panel_B_genetics_vs_depletion.csv") |> add_depletion_status() |>
    mutate(depletion = -comp_beta)   # plot depletion (−β): depleted types point up
  rp <- spearman_rp(df$scz_neg_log10_p, df$depletion)

  # compact: uniform dot size and no legend (used in the trimmed-down variant).
  pts <- if (compact)
    geom_point(aes(fill = color, color = depleted_status, stroke = depleted_status),
               shape = 21, size = 2.8, alpha = 0.92)
  else
    geom_point(aes(fill = color, size = pt_size,
                   color = depleted_status, stroke = depleted_status),
               shape = 21, alpha = 0.92)

  ggplot(df, aes(x = scz_neg_log10_p, y = depletion)) +
    geom_smooth(method = "lm", formula = y ~ x, color = "#444",
                fill = "#888", alpha = 0.18, linewidth = 0.3,
                linetype = "dashed", se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "#aaa", linewidth = 0.2) +
    pts +
    geom_text_repel(aes(label = supertype), size = LBL_GENE,
                    box.padding = 0.25, point.padding = 0.20,
                    max.overlaps = 30, segment.color = "#999",
                    segment.size = 0.2, force = 3, min.segment.length = 0) +
    scale_fill_identity() +
    depletion_scales(show_legend = !compact) +
    (if (!compact) scale_size(range = c(1.2, 3), guide = "none")) +
    inset_spearman(rp$rho, rp$p, corner = "top-left", size = LBL_STAT) +
    coord_cartesian(ylim = DEPLETION_LIM) +
    labs(x = LAB_GWAS, y = LAB_DEPL) +
    theme_panel() +
    (if (compact) theme(legend.position = "none")
     else theme(legend.position = c(0.98, 0.02),
                legend.justification = c("right", "bottom")))
}


# ====================================================================
# C) HCN1 expression vs SCZ MAGMA enrichment across 16 Sst supertypes
# ====================================================================
# The direct convergence statement: per-Sst-supertype mean HCN1 expression
# tracks SCZ genetic enrichment. Same x-axis units as Panel E (log2 CP10K+1)


# ====================================================================
# D) HCN1 locus zoom (Manhattan + gene track, hg38)
# ====================================================================
build_hcn1_locus_plot <- function() {
  snps    <- read_panel("panel_D_snps.csv")
  cs      <- read_panel("panel_D_credible_set.csv")
  genes   <- read_panel("panel_D_genes.csv")
  exons   <- read_panel("panel_D_exons.csv")

  # Keep the track focused on HCN1: drop the neighbouring uncharacterised LOC
  # RNA gene and collapse to a single lane (it was the only other occupant).
  genes <- genes |> filter(!grepl("^LOC", symbol)) |>
    mutate(lane = 0L, n_lanes = 1L)
  exons <- exons |> filter(!grepl("^LOC", symbol)) |> mutate(lane = 0L)
  meta    <- read_panel("panel_D_meta.csv")

  # --- top: Manhattan with the SuSiE-R credible set ---
  # PIPs come from Bigdeli 2026 Supplementary Table 13 (SuSiE-R,
  # European-ancestry), NOT FINEMAP -- the legend used to say FINEMAP
  # because the panel was built on PGC3.
  manhattan <- ggplot(snps, aes(x = pos_mb, y = neg_log10_p)) +
    annotate("rect", xmin = meta$hcn1_start_mb, xmax = meta$hcn1_stop_mb,
             ymin = -Inf, ymax = Inf, fill = HCN1_COLOR, alpha = 0.10) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed",
               color = "purple", linewidth = 0.2, alpha = 0.6) +
    geom_point(data = filter(snps, !in_credible_set),
               color = "#aaa", alpha = 0.45, size = 0.4, stroke = 0) +
    geom_point(data = filter(snps, !in_credible_set & is_gwas_sig),
               color = "#F4A6A6", alpha = 0.75, size = 0.5, stroke = 0) +
    geom_point(data = filter(cs, !is_lead),
               aes(x = pos_mb, y = neg_log10_p, fill = pip, size = pip),
               shape = 21, color = "black", stroke = 0.3) +
    geom_point(data = filter(cs, is_lead),
               aes(x = pos_mb, y = neg_log10_p, fill = pip, size = pip),
               shape = 23, color = "black", stroke = 0.4) +
    scale_fill_gradientn(
      colours = c("#FFEDA0", "#FEB24C", "#F03B20", "#BD0026"),
      limits  = c(0, max(cs$pip)),
      breaks  = c(0.05, 0.5),       # fewer ticks → more compact legend
      trans   = "sqrt",
      name    = "SuSiE-R PIP") +
    scale_size(range = c(1.0, 3.2), guide = "none") +
    # HCN1 is on the minus strand: reversing the coordinate axis makes
    # transcription read left-to-right with the 3' end on the right, which is
    # what readers expect. Applied to both stacked subplots so they stay
    # aligned; the direction chevrons follow the scale and now point right.
    scale_x_reverse() +
    coord_cartesian(xlim = c(meta$win_lo_mb, meta$win_hi_mb)) +
    labs(x = NULL, y = expression("SCZ GWAS ("*-log[10]~italic(P)*")")) +
    theme_panel() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          # Reversing the x-axis moved the credible set to the upper RIGHT,
          # which is where this legend used to sit. Anchor it left instead.
          legend.position = c(0.02, 0.98),
          legend.justification = c("left", "top"),
          legend.key.width = unit(0.15, "cm"),
          legend.key.height = unit(0.18, "cm"),
          plot.margin = margin(2, 4, 0, 4))

  # --- bottom: hg38 RefSeq gene track ---
  n_lanes <- max(genes$n_lanes)
  genes_with_y <- genes |> mutate(y = n_lanes - 1 - lane)
  exons_with_y <- exons |>
    left_join(select(genes_with_y, symbol, y, strand),
              by = c("symbol", "strand"))

  # HCN1 has 10 exon segments with median width ~260 bp; at this panel's
  # render width (a 700 kb window in ~2 in of horizontal space), real-scale
  # exons are <0.001 in wide and effectively invisible. Enforce a minimum
  # displayed width of ~0.5% of the window (≈ 3.5 kb here) so each exon
  # renders as a visible shape — standard genome-browser convention.
  window_mb         <- meta$win_hi_mb - meta$win_lo_mb
  min_exon_width_mb <- window_mb * 0.005
  exons_with_y <- exons_with_y |>
    mutate(width_mb     = stop_mb - start_mb,
           center_mb    = (start_mb + stop_mb) / 2,
           display_w    = pmax(width_mb, min_exon_width_mb),
           start_mb_dsp = center_mb - display_w / 2,
           stop_mb_dsp  = center_mb + display_w / 2)

  # Transcription-direction chevrons spaced along each gene body, pointing in
  # the strand direction (HCN1 is on the minus strand, so they point left).
  # Standard genome-browser convention; drawn under the exon boxes.
  arrow_marks <- do.call(rbind, lapply(seq_len(nrow(genes_with_y)), function(i) {
    g   <- genes_with_y[i, ]
    pos <- seq(g$start_mb, g$stop_mb, length.out = 9)[2:8]
    dx  <- window_mb * 0.020 * ifelse(g$strand == "-", -1, 1)
    data.frame(x = pos, xend = pos + dx, y = g$y, is_hcn1 = g$is_hcn1)
  }))

  gene_track <- ggplot() +
    geom_segment(data = genes_with_y,
                 aes(x = start_mb, xend = stop_mb, y = y, yend = y,
                     color = factor(is_hcn1)),
                 linewidth = 0.5, alpha = 0.75) +
    geom_segment(data = arrow_marks,
                 aes(x = x, xend = xend, y = y, yend = y,
                     color = factor(is_hcn1)),
                 linewidth = 0.45, alpha = 0.95,
                 arrow = arrow(length = unit(0.09, "cm"), type = "open")) +
    geom_rect(data = filter(exons_with_y, seg_type == "CDS"),
              aes(xmin = start_mb_dsp, xmax = stop_mb_dsp,
                  ymin = y - 0.40, ymax = y + 0.40,
                  fill = factor(is_hcn1)),
              color = "black", linewidth = 0.15, alpha = 0.95) +
    geom_rect(data = filter(exons_with_y, seg_type == "UTR"),
              aes(xmin = start_mb_dsp, xmax = stop_mb_dsp,
                  ymin = y - 0.24, ymax = y + 0.24,
                  fill = factor(is_hcn1)),
              color = "black", linewidth = 0.12, alpha = 0.6) +
    geom_text(data = filter(genes_with_y, !is_hcn1),
              aes(x = (pmax(start_mb, meta$win_lo_mb) +
                       pmin(stop_mb,  meta$win_hi_mb)) / 2,
                  y = y + 0.55, label = symbol),
              color = "#333", fontface = "italic", size = 1.8, hjust = 0.5) +
    geom_text(data = filter(genes_with_y, is_hcn1),
              aes(x = (pmax(start_mb, meta$win_lo_mb) +
                       pmin(stop_mb,  meta$win_hi_mb)) / 2,
                  y = y + 0.80, label = symbol),
              color = HCN1_COLOR, fontface = "bold.italic", size = 3.6, hjust = 0.5) +
    scale_color_manual(values = c("FALSE" = "#333", "TRUE" = HCN1_COLOR),
                       guide = "none") +
    scale_fill_manual(values = c("FALSE" = "#2166ac", "TRUE" = HCN1_COLOR),
                      guide = "none") +
    scale_x_reverse() +
    coord_cartesian(xlim = c(meta$win_lo_mb, meta$win_hi_mb),
                    ylim = c(-0.5, n_lanes - 0.5 + 0.9)) +
    labs(x = "Position on chr5 (hg38, Mb)", y = NULL) +
    theme_panel() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = "#fafafa", color = NA),
          plot.margin = margin(0, 4, 2, 4))

  # Give the gene track ~50% more vertical real estate than before
  # (rel_heights 2.0:1.0 instead of 3:1) so HCN1 exon structure is legible.
  plot_grid(manhattan, gene_track, ncol = 1, align = "v", axis = "lr",
            rel_heights = c(2.0, 1.0))
}


# ====================================================================
# E) HCN1 normalized expression vs sag (16 Sst supertypes)
# ====================================================================
build_hcn1_expression_vs_sag <- function(compact = FALSE) {
  # Convert per-supertype HCN1 mean from log1p(CP10K) → log2(CP10K + 1):
  # the underlying argument is the same; division by ln(2) just changes
  # the log base for readability (log2 is conventional in genomics).
  #
  # Point size encodes SCZ MAGMA enrichment (−log10 P, joined from panel B)
  # — same field as the y-axis of Panel A and Panel C — so the reader can
  # see at a glance that the high-HCN1 / high-sag supertypes are also the
  # genetically most-enriched ones. (n_cells is still used to weight the
  # OLS fit, since mean_sag is a within-supertype average that should be
  # less trusted for supertypes with few Patch-seq cells.)
  panel_e <- read_panel("panel_E_hcn1_vs_sag.csv")
  panel_b <- read_panel("panel_B_genetics_vs_depletion.csv")

  df <- panel_e |>
    inner_join(select(panel_b, supertype, scz_neg_log10_p),
               by = "supertype") |>
    add_depletion_status() |>
    mutate(HCN1_expr_log2 = HCN1_expr / log(2))

  rp <- spearman_rp(df$HCN1_expr_log2, df$mean_sag)

  # compact: uniform dot size and no size legend (used in the trimmed variant).
  pts <- if (compact)
    geom_point(aes(fill = color, color = depleted_status, stroke = depleted_status),
               shape = 21, size = 2.8, alpha = 0.92)
  else
    geom_point(aes(fill = color, size = scz_neg_log10_p,
                   color = depleted_status, stroke = depleted_status),
               shape = 21, alpha = 0.92)

  ggplot(df, aes(x = HCN1_expr_log2, y = mean_sag)) +
    geom_smooth(method = "lm", formula = y ~ x, color = "#444",
                fill = "#888", alpha = 0.15, linewidth = 0.35,
                linetype = "dashed",
                aes(weight = sqrt(n_cells))) +
    pts +
    geom_text_repel(aes(label = supertype), size = LBL_GENE,
                    box.padding = 0.25, point.padding = 0.20,
                    max.overlaps = 30, segment.color = "#999",
                    segment.size = 0.2, force = 3, min.segment.length = 0) +
    scale_fill_identity() +
    # Depletion legend suppressed here — Panel B carries it for the figure.
    depletion_scales(show_legend = FALSE) +
    (if (!compact) scale_size_continuous(
      range = c(1.2, 3.5),
      breaks = c(2, 5, 10),
      name = expression(-log[10](italic(P)[MAGMA])))) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
    inset_spearman(rp$rho, rp$p, corner = "bottom-right", size = LBL_STAT) +
    labs(x = LAB_HCN1, y = "Sag ratio (Patch-seq)") +
    theme_panel() +
    (if (compact) theme(legend.position = "none")
     else theme(legend.position = c(0.02, 0.98),
                legend.justification = c("left", "top"),
                legend.key.height = unit(0.25, "cm"),
                legend.spacing.y = unit(0.02, "cm")))
}


# ====================================================================
# F) Cell morphologies (Sst_25 and Sst_5)
# ====================================================================
# Cortical layer model (canonical MTG fractions used in patchseq_builder).
CORTEX_THICKNESS_UM <- 3006
LAYER_FRACTIONS <- list(L1_L2 = 0.089, L2_L3 = 0.149, L3_L4 = 0.427)

build_morphology_plot <- function() {
  segs <- read_panel("panel_F_morphology.csv")
  meta <- read_panel("panel_F_meta.csv")

  # Per-component styling: dendrites + soma drawn thick at full opacity;
  # axons drawn thin and semi-transparent so they don't dominate.
  segs <- segs |> mutate(
    is_axon = comp_type == 2,
    lw      = ifelse(is_axon, 0.12, 0.35),
    alpha_v = ifelse(is_axon, 0.55, 1.0))

  # Place the two cells side-by-side at common µm scale.
  half_w  <- max(abs(range(c(segs$x0, segs$x1))))
  x_gap   <- 170
  # Generalised from the original two-cell layout. The original spaced both
  # cells by the GLOBAL half-width; with three cells that wastes horizontal
  # room and, under coord_fixed(), shrinks the whole panel. Pack each cell by
  # its OWN half-width instead, edge to edge.
  ord_cells <- meta$supertype[order(meta$cell_order)]
  hw_tbl <- segs |> group_by(cell) |>
    summarise(hw = max(abs(c(x0, x1)), na.rm = TRUE), .groups = "drop")
  hw <- hw_tbl$hw[match(ord_cells, hw_tbl$cell)]
  # Even centre-to-centre spacing (edge-to-edge packing left visibly uneven
  # gaps, because the cells differ a lot in width). Pitch is set by the widest
  # adjacent PAIR so nothing can overlap, then centres are distributed evenly.
  n_cells <- length(ord_cells)
  pitch <- max(hw[-n_cells] + hw[-1]) + x_gap
  shifts <- (seq_len(n_cells) - (n_cells + 1) / 2) * pitch
  names(shifts) <- meta$supertype[order(meta$cell_order)]
  segs <- segs |> mutate(x0_sh = x0 + shifts[cell],
                         x1_sh = x1 + shifts[cell])
  meta <- meta |> mutate(soma_x_um_sh = soma_x_um + shifts[supertype])

  # Per-cell label-y: just above each cell's topmost dendrite/axon.
  cell_tops <- segs |> group_by(cell) |>
    summarise(top_y = min(pmin(y0, y1), na.rm = TRUE)) |>
    rename(supertype = cell)
  # Both labels sit at a common height in the headroom ABOVE the pia line, so
  # they cannot collide with the pia/L1-L2/L2-L3/L3-L4 boundaries.
  meta <- meta |> left_join(cell_tops, by = "supertype") |>
    mutate(label_y = -140)

  # Vertical zoom: from pia (with a bit of headroom for "Pia" label) to
  # just below the deepest dendrite + scale-bar room.
  max_y    <- max(segs$y0, segs$y1)
  y_bottom <- max_y + 330    # scale bar now sits inside the panel (lower left)
  y_top    <- -290   # headroom for the cell labels, which sit above the
                     # pia line and above the "Pia" annotation at y = -50
  # extra left margin: the L1/L2/L3 and Pia annotations are anchored at the
  # left edge and would otherwise sit on top of the first cell.
  x_extent <- c(min(shifts - hw) - 260, max(shifts + hw) + 40)

  # Layer-boundary lines (Pia is solid; L1/L2, L2/L3, L3/L4 are dashed).
  layer_lines <- tibble(
    y         = c(0,
                   LAYER_FRACTIONS$L1_L2,
                   LAYER_FRACTIONS$L2_L3,
                   LAYER_FRACTIONS$L3_L4) * c(1, rep(CORTEX_THICKNESS_UM, 3)),
    linewidth = c(0.25, 0.18, 0.18, 0.18),
    linetype  = c("solid", "dashed", "dashed", "dashed"))

  # Layer text positions (LEFT side of panel).
  layer_label_x <- x_extent[1] + 25
  layer_labels  <- tibble(
    text = c("L1", "L2", "L3"),
    y    = c( LAYER_FRACTIONS$L1_L2 / 2,
              (LAYER_FRACTIONS$L1_L2 + LAYER_FRACTIONS$L2_L3) / 2,
              (LAYER_FRACTIONS$L2_L3 + LAYER_FRACTIONS$L3_L4) / 2) *
             CORTEX_THICKNESS_UM)

  # Depleted set is the Fig-3 vulnerable group; everything else is not depleted.
  DEPLETED_SUPERTYPES <- c("Sst_25", "Sst_22", "Sst_2", "Sst_20", "Sst_3")
  bar_y_val <- max_y + 200
  status_bar <- meta |>
    mutate(status = ifelse(supertype %in% DEPLETED_SUPERTYPES,
                           "Depleted in SCZ", "Not depleted"),
           hw_i = hw[match(supertype, ord_cells)],
           lo = soma_x_um_sh - hw_i, hi = soma_x_um_sh + hw_i) |>
    group_by(status) |>
    summarise(x_start = min(lo), x_end = max(hi), .groups = "drop") |>
    mutate(bar_y = bar_y_val,
           x_start = x_start + 40, x_end = x_end - 40,   # gap between groups
           status_col = ifelse(status == "Depleted in SCZ", "#3d3d3d", "#8c8c8c"))

  ggplot() +
    geom_hline(data = layer_lines,
               aes(yintercept = y, linewidth = linewidth, linetype = linetype),
               color = "#999", alpha = 0.55, show.legend = FALSE) +
    geom_segment(data = arrange(segs, is_axon),  # dendrites painted on top
                 aes(x = x0_sh, xend = x1_sh, y = y0, yend = y1,
                     color = cell, alpha = alpha_v, linewidth = lw)) +
    geom_point(data = meta,
               aes(x = soma_x_um_sh, y = soma_y_um),
               shape = 21, color = "black", fill = "black", size = 1.2) +
    # Cell labels (just above each cell's topmost dendrite, colored by supertype).
    geom_text(data = meta,
              aes(x = soma_x_um_sh, y = label_y,
                  label = supertype,
                  color = supertype),
              size = LBL_CALL, hjust = 0.5, vjust = 0) +
    # Layer text on the LEFT.
    geom_text(data = layer_labels,
              aes(x = layer_label_x, y = y, label = text),
              color = "#777", size = LBL_SMALL, hjust = 0, vjust = 0.5,
              fontface = "italic") +
    annotate("text", x = layer_label_x, y = -50, label = "Pia",
             color = "#555", size = LBL_SMALL, hjust = 0) +
    # Compositional status, bracketed under the cells it applies to. Groups
    # come from the Fig-3 compositional analysis, not from anything drawn here.
    geom_segment(data = status_bar,
                 aes(x = x_start, xend = x_end, y = bar_y, yend = bar_y,
                     color = status_col),
                 linewidth = 0.4, inherit.aes = FALSE) +
    geom_text(data = status_bar,
              aes(x = (x_start + x_end) / 2, y = bar_y + 55,
                  label = status, color = status_col),
              size = LBL_SMALL, vjust = 1, inherit.aes = FALSE) +
    # 100 µm scale bar (label offset clearly below the bar so they don't overlap).
    annotate("segment",
             x    = x_extent[1] + 45,
             xend = x_extent[1] + 145,
             y    = max_y - 190, yend = max_y - 190,
             linewidth = 0.8) +
    # vjust = 1 anchors the label's top edge, so with scale_y_reverse() the text
    # is pushed away from the bar rather than back onto it (it collided once the
    # base font was raised).
    annotate("text", x = x_extent[1] + 95,
             y = max_y - 150, label = "100 µm",
             size = LBL_SMALL, hjust = 0.5, vjust = 1) +
    scale_color_manual(values = setNames(meta$color, meta$supertype),
                       guide = "none") +
    scale_linewidth_identity() +
    scale_alpha_identity() +
    scale_linetype_identity() +
    scale_y_reverse() +     # depth increases downward; pia on top
    coord_fixed(ratio = 1, xlim = x_extent, ylim = c(y_bottom, y_top),
                clip = "off") +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = margin(8, 10, 8, 22))   # shifted right within its cell
}


# ====================================================================
# G) Patch-seq voltage traces
# ====================================================================
build_ephys_traces_plot <- function() {
  traces <- read_panel("panel_G_traces.csv")
  meta   <- read_panel("panel_G_meta.csv")

  # Sst_25 (the high-sag exemplar) drives the sag-annotation glyphs.
  sst25       <- filter(meta, supertype == "Sst_25")
  sst25_color <- sst25$color

  v_range <- range(traces$v_mV)
  v_pad   <- diff(v_range) * 0.10

  # Per-trace supertype labels, placed mid-step where the two traces are
  # most separated. Sit ~5 mV above each cell's steady-state V — comfortably
  # above the trace, away from the sag-annotation glyphs at the front of
  # the step (V_peak_t_ms ≈ 50 ms for Sst_25).
  # Labels sit just above each trace's own plateau and are staggered in time,
  # so cells with similar steady-state voltage (Sst_20 and Sst_3 plateaus
  # differ by ~1.5 mV) cannot collide. A short leader drops from the label to
  # the trace.
  LABEL_POS <- tibble::tribble(
    ~supertype, ~label_t_ms, ~label_v_mV,
    "Sst_25",         845,       -60.8,   # own clear space above its plateau
    # "Sag" now sits at LBL_CALL and occupies ~160-340 ms in the same
    # -65..-74 band, so these two start after it.
    "Sst_3",          560,       -69.2,
    "Sst_20",         850,       -71.5,   # later and lower than Sst_3's label
    "Sst_1",          880,       -82.3,   # below, so its leader crosses nothing
    "Sst_5",          430,       -82.9)   # empty band above the deepest trace
  trace_labels <- meta |>
    left_join(LABEL_POS, by = "supertype") |>
    mutate(leader_v = V_steady + 0.4)

  ggplot(traces, aes(x = t_ms, y = v_mV, color = cell)) +
    geom_line(linewidth = 0.45, alpha = 0.95) +
    geom_segment(data = trace_labels,
                 aes(x = label_t_ms, xend = label_t_ms,
                     y = label_v_mV, yend = leader_v, color = supertype),
                 linewidth = 0.25, alpha = 0.8, inherit.aes = FALSE) +
    geom_text(data = trace_labels,
              aes(x = label_t_ms, y = label_v_mV,
                  label = supertype, color = supertype),
              size = LBL_CALL, hjust = 0.5, vjust = 0,
              inherit.aes = FALSE) +
    # Sag annotation on Sst_25: peak (circle) + steady-state (square) +
    # double-headed arrow labeled "Sag".
    geom_point(data = sst25, aes(x = V_peak_t_ms, y = V_peak),
               shape = 21, color = sst25_color, fill = "white",
               size = 2, stroke = 0.7, inherit.aes = FALSE) +
    geom_point(data = sst25, aes(x = 950, y = V_steady),
               shape = 22, color = sst25_color, fill = "white",
               size = 2, stroke = 0.7, inherit.aes = FALSE) +
    geom_segment(data = sst25,
                 aes(x = V_peak_t_ms + 60, xend = V_peak_t_ms + 60,
                     y = V_peak, yend = V_steady),
                 color = sst25_color, linewidth = 0.35,
                 arrow = arrow(ends = "both", length = unit(0.08, "cm")),
                 inherit.aes = FALSE) +
    geom_text(data = sst25,
              aes(x = V_peak_t_ms + 110, y = (V_peak + V_steady) / 2,
                  label = "Sag"),
              color = sst25_color, size = LBL_CALL,
              hjust = 0, inherit.aes = FALSE) +
    scale_color_manual(values = setNames(meta$color, meta$supertype),
                       guide = "none") +
    coord_cartesian(ylim = c(v_range[1] - v_pad * 0.15, v_range[2] + v_pad * 0.1)) +
    labs(x = "Time (ms)",
         y = "Membrane potential (mV)") +
    theme_panel() +
    theme(legend.position = "none")
}


# ====================================================================
# C (restored)) Sst_25 gene-driver scatter
# ====================================================================
# Panel b: which genes drive Sst_2's SCZ enrichment. X = gene specificity in
# Sst_2 (log10); Y = -log10 SCZ MAGMA gene p. Drivers (Sst_2-specific AND
# GWAS-significant) are filled, other genes grey, HCN1 called out. Guide lines:
# Sst_2 specificity 90th percentile (vertical), MAGMA FDR 0.05 and genome-wide
# 5e-8 (horizontal). Data: panel_C_gene_drivers / _top_labels / _meta, from
# export_panels_abc.py. Sst_2 is both the most SCZ-enriched Sst supertype and
# the most significantly depleted, so this panel names the type that sits at
# the corner of panel a.
GENOME_WIDE_NEG_LOG10P <- -log10(5e-8)   # 7.30

build_gene_driver_plot <- function() {
  g    <- read_panel("panel_C_gene_drivers.csv")
  lab  <- read_panel("panel_C_top_labels.csv")
  meta <- read_panel("panel_C_meta.csv")

  # clip="off" (needed so repelled labels may sit outside the panel) also lets
  # points outside ylim render below the axis -- 11,457 genes fall under
  # ylim_lo. Restrict the point layers to the displayed range instead.
  in_view <- function(d) filter(d, neg_log10_p >= meta$ylim_lo,
                                  neg_log10_p <= meta$ylim_hi)
  g_grey <- g |> filter(!is_driver, !is_hcn1) |> in_view()
  g_red  <- g |> filter(is_driver,  !is_hcn1) |> in_view()
  g_hcn1 <- g |> filter(is_hcn1)
  lab_g  <- lab |> filter(!is_hcn1)

  ggplot(mapping = aes(specificity, neg_log10_p)) +
    geom_vline(xintercept = meta$spec_q90, linetype = "dashed",
               color = "#c4c4c4", linewidth = 0.25) +
    geom_hline(yintercept = meta$fdr_neg_log10p_cutoff, linetype = "dashed",
               color = "#c4c4c4", linewidth = 0.25) +
    geom_hline(yintercept = GENOME_WIDE_NEG_LOG10P, linetype = "dotted",
               color = "purple", linewidth = 0.3, alpha = 0.7) +
    geom_point(data = g_grey, color = "grey80", size = 0.35, alpha = 0.5) +
    geom_point(data = g_red,  color = "#d1483d", size = 0.7,  alpha = 0.8) +
    geom_point(data = g_hcn1, fill = HCN1_COLOR, color = "white",
               shape = 21, size = 2.8, stroke = 0.3) +
    # The HCN1 marker is drawn in its own layer, so the gene-label repel cannot
    # see it and parked RBFOX1 underneath it. Carry HCN1 into this layer with an
    # empty label: it repels as an obstacle but draws nothing.
    geom_text_repel(data = bind_rows(lab_g,
                                     transmute(g_hcn1, specificity, neg_log10_p,
                                               symbol = "")),
                    aes(label = symbol), size = LBL_GENE,
                    fontface = "italic", color = "#333",
                    box.padding = 0.30, max.overlaps = 40,
                    segment.color = "#bbb", segment.size = 0.2,
                    min.segment.length = 0, force = 4, seed = 1) +
    # callout goes up and to the LEFT, over the unlabelled sub-threshold cloud;
    # to the right it lands on the driver labels (RBMS3 for the Sst_20 target)
    geom_text_repel(data = g_hcn1, aes(label = symbol), size = LBL_CALL,
                    fontface = "bold.italic", color = HCN1_COLOR,
                    nudge_y = 3.2, nudge_x = -0.30, box.padding = 0.6,
                    segment.color = HCN1_COLOR, segment.size = 0.4,
                    min.segment.length = 0, seed = 1) +
    # The two threshold lines carry no in-plot text: there is no empty region to
    # park it in (the right edge collides with the driver-gene labels, the left
    # edge sits on the point cloud), and the figure legend already identifies
    # them -- "horizontal lines, MAGMA FDR 0.05 (grey) and genome-wide
    # significance (5e-8, purple)". Descriptive text belongs in the legend.
    scale_x_log10(limits = c(meta$xlim_lo, meta$xlim_hi),
                  breaks = c(0.003, 0.01, 0.03, 0.1, 0.3)) +
    coord_cartesian(ylim = c(meta$ylim_lo, meta$ylim_hi), clip = "off") +
    labs(x = sprintf("Gene specificity in %s (log scale)", meta$target_type),
         y = expression(-log[10]*"(SCZ MAGMA gene "*italic(p)*")")) +
    theme_panel()
}


# ====================================================================
# Assemble & save
# ====================================================================
hcn1_locus_plot         <- build_hcn1_locus_plot()
morphology_plot         <- build_morphology_plot()
ephys_traces_plot       <- build_ephys_traces_plot()

panel_label_kwargs <- list(
  label_size = PANEL_LABEL_SIZE, label_fontface = "bold",
  label_x = 0.0, label_y = 1.0,  hjust = -0.3, vjust = 1.3)



# Panel b of Figure 4. (Built here, where the gene-driver panel was
# originally introduced for the retired 7-panel variant.)
gene_driver_plot         <- build_gene_driver_plot()

# ====================================================================
# Compact-form panels a, d and e. `compact = TRUE` drops the in-panel
# legend and uses a uniform dot size; the depleted/not-depleted outline
# encoding is identical in every scatter, so it is described once in the
# figure legend rather than repeated in each panel.
# ====================================================================
genetics_vs_depletion_c <- build_genetics_vs_depletion(compact = TRUE)
sag_c                   <- build_hcn1_expression_vs_sag(compact = TRUE)


# ====================================================================
# Variant: full Fig 4 — the compact layout plus a third row that
# characterises the vulnerable population molecularly (marker volcano,
# CALB1) and shows the same Sst supertypes decline in Alzheimer's disease.
# Builders are shared with the standalone renderer via fig4_new_panels.R.
# ====================================================================
source(file.path(REPO, "scripts", "figures", "fig4_new_panels.R"))
source(file.path(REPO, "scripts", "figures", "fig4_assemble_nod.R"))

figure_full <- build_figure4_nod(list(
  a = genetics_vs_depletion_c,             # SCZ GWAS enrichment vs depletion
  b = gene_driver_plot,                    # Sst_2 gene drivers
  c = hcn1_locus_plot,                     # HCN1 locus zoom + fine-mapping
  e = sag_c,                               # HCN1 expression vs patch-seq sag
  f = morphology_plot,                     # exemplar reconstructions
  g = ephys_traces_plot,                   # exemplar voltage traces
  h = build_marker_volcano("cell"),              # depleted vs not-depleted markers
  i = build_violin("CALB1"),               # CALB1 by depletion group
  j = build_ad_concordance(show_legend = FALSE)))   # SCZ vs AD (DLPFC)

# Written straight into the submission folder, as Figure 2 and the supplementary
# figures are; the figure number lives in the stem.
MAINFIG <- file.path(dirname(REPO), "manuscript", "figures", "main")
OUT_FULL_STEM <- file.path(MAINFIG, "Fig4_scz_genetics_hcn1_sst")
render_figure4(figure_full, OUT_FULL_STEM)
message(glue("\u2192 Saved {OUT_FULL_STEM}.{{png,pdf,svg}}  (Figure 4, 9 panels)"))