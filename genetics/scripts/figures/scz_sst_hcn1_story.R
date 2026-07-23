#!/usr/bin/env Rscript
# ===================================================================
# scz_sst_hcn1_story.R
#
# Renders the 7-panel SCZ-cell-type-enrichment / HCN1 / Sst figure.
# Data are pre-extracted by scripts/figures/export_for_R.py into CSVs
# in results/figures/r_panels/. This script reads those CSVs, builds
# each panel with ggplot2, and assembles them with cowplot::plot_grid.
#
# Panel guide (top→bottom, left→right):
#   A  enrichment_landscape    SCZ enrichment across 137 SEA-AD supertypes
#   B  genetics_vs_depletion   Genetic risk × SCZ cell-abundance depletion
#   C  hcn1_vs_scz_plot        HCN1 expression × SCZ enrichment (Sst only)
#   D  hcn1_locus_plot         HCN1 locus zoom + gene track (hg38)
#   E  hcn1_expression_vs_sag  HCN1 expression × intrinsic sag
#   F  morphology_plot         Two example Sst-cell morphologies
#   G  ephys_traces_plot       Voltage responses for the same two cells
#
# Convergence triangle across panels B, C, E (all on the same 16 Sst
# supertypes, sharing the SEA-AD fill + depletion-outline encoding):
#   B: genetics  ↔ compositional depletion
#   C: HCN1 expression ↔ genetics
#   E: HCN1 expression ↔ intrinsic sag
# The gene-level driver-plot for Sst_25 (previous Panel C) moves to a
# supplementary figure.
# ===================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(glue)
  library(ggplot2); library(cowplot); library(ggrepel); library(scales)
})

# ────────────────────────────────────────────────────────────────────
# Paths
# ────────────────────────────────────────────────────────────────────
REPO     <- "/Users/shreejoy/Github/scz_celltype_paper/genetics"
DATA_DIR <- file.path(REPO, "results", "figures", "r_panels")
OUT_PNG  <- file.path(REPO, "results", "figures", "scz_sst_hcn1_multipanel_R_v3.png")
OUT_PDF  <- file.path(REPO, "results", "figures", "scz_sst_hcn1_multipanel_R_v3.pdf")
# Former Panel A (enrichment landscape) is now a standalone supplement.
OUT_SUPPA_PNG <- file.path(REPO, "results", "figures", "supp_scz_enrichment_landscape_R.png")
OUT_SUPPA_PDF <- file.path(REPO, "results", "figures", "supp_scz_enrichment_landscape_R.pdf")

read_panel <- function(filename) {
  read_csv(file.path(DATA_DIR, filename), show_col_types = FALSE)
}

# ────────────────────────────────────────────────────────────────────
# Visual constants
# ────────────────────────────────────────────────────────────────────
# Sized for a Nature Neuroscience-style double-column figure: 6.5" wide.
# Text scales proportionally LESS than spatial elements so it stays
# legible at print scale (~7 pt body, ~12 pt panel labels).
BASE_FONT_SIZE   <- 7
PANEL_LABEL_SIZE <- 12
FIG_WIDTH_IN     <- 6.5
FIG_HEIGHT_IN    <- 6.36   # extra 0.75 in over the aspect-preserving 5.61;
                            # the extra height is routed to rows 2/3 (see
                            # rel_heights in the assembly block) so the middle
                            # and bottom panels are less vertically crowded.

HCN1_COLOR <- "#1565C0"   # Used in C (gene driver callout) and D (gene track)

# Depletion encoding: thick black outline = depleted in SCZ post-mortem at
# FDR < 0.20; thin grey outline = not depleted. Used in panels B and E.
DEPLETION_OUTLINE <- c(`Depleted (FDR < 0.20)` = "black",
                       `Not depleted`         = "grey60")
DEPLETION_STROKE  <- c(`Depleted (FDR < 0.20)` = 0.7,
                       `Not depleted`         = 0.2)

# Plain-language axis labels, shared across panels for a consistent vocabulary.
LAB_GWAS <- expression("SCZ GWAS enrichment ("*-log[10]~italic(P)*")")
LAB_DEPL <- expression("Cell depletion in SCZ ("*-beta*")")
LAB_HCN1 <- expression(italic("HCN1")*" expression ("*log[2]*" CP10K+1)")

# ────────────────────────────────────────────────────────────────────
# Shared theme + helpers
# ────────────────────────────────────────────────────────────────────

# Common theme. Print-size text, no bold by default (bold reserved for
# in-panel callouts and panel labels). Drop in via `+ theme_panel()`.
theme_panel <- function(base_size = BASE_FONT_SIZE) {
  theme_cowplot(font_size = base_size) +
    theme(
      plot.title        = element_blank(),
      plot.subtitle     = element_blank(),
      panel.grid.major  = element_line(color = "grey92", linewidth = 0.15),
      panel.grid.minor  = element_blank(),
      axis.title        = element_text(size = base_size + 1),
      axis.text         = element_text(size = base_size - 0.5),
      axis.line         = element_line(linewidth = 0.3, color = "grey20"),
      axis.ticks        = element_line(linewidth = 0.25, color = "grey20"),
      legend.title      = element_text(size = base_size - 0.5),
      legend.text       = element_text(size = base_size - 1),
      legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
      legend.key.size   = unit(0.30, "cm"),
      plot.margin       = margin(3, 4, 3, 4))
}

# Add a `depleted_status` factor (with the canonical level order) to a
# supertype-keyed table, for color/stroke mapping in B and E.
add_depletion_status <- function(df) {
  df |> mutate(
    depleted_status = factor(
      ifelse(depleted_fdr20, "Depleted (FDR < 0.20)", "Not depleted"),
      levels = c("Depleted (FDR < 0.20)", "Not depleted")))
}

# Color and stroke scales for the depletion encoding.
# show_legend = TRUE produces a clean legend with grey-filled keys; FALSE
# suppresses the legend (used in panel E to avoid duplicating B's legend).
depletion_scales <- function(show_legend = TRUE) {
  color_guide <- if (show_legend) {
    guide_legend(override.aes = list(
      shape = 21, fill = "grey85", size = 2.5,
      stroke = unname(DEPLETION_STROKE)))
  } else {
    "none"
  }
  list(
    scale_color_manual(
      values = DEPLETION_OUTLINE,
      breaks = names(DEPLETION_OUTLINE),
      name   = NULL,
      guide  = color_guide),
    scale_discrete_manual("stroke",
      values = DEPLETION_STROKE, guide = "none"))
}

# Inset Spearman ρ — corner annotation used in B, C, E.
# Supports all four corners so positive- and negative-slope scatters can
# each park their annotation in their natural empty corner.
inset_spearman <- function(rho, p,
                            corner = c("top-right", "bottom-right",
                                       "top-left",  "bottom-left"),
                            size = 6) {
  corner <- match.arg(corner)
  xval  <- if (grepl("right",  corner))  Inf else -Inf
  yval  <- if (grepl("top",    corner))  Inf else -Inf
  hjust <- if (grepl("right",  corner))  1.05 else -0.05
  vjust <- if (grepl("top",    corner))  1.5  else -1.0
  annotate("text", x = xval, y = yval, hjust = hjust, vjust = vjust,
           label = sprintf("rho == %.2f * ',' ~ italic(p) == %.3g", rho, p),
           parse = TRUE, size = size, color = "#222")
}

# Pretty Spearman correlation and its asymptotic p-value
spearman_rp <- function(x, y) {
  list(rho = cor(x, y, method = "spearman"),
       p   = cor.test(x, y, method = "spearman")$p.value)
}


# ====================================================================
# A) SCZ enrichment landscape
# ====================================================================
# Neuronal families (used when neurons_only = TRUE in the builder).
NEURONAL_FAMILIES <- c(
  "Sst", "Sst Chodl", "Pvalb", "Chandelier", "Vip", "Sncg",
  "Lamp5", "Pax6", "Lamp5 Lhx6",
  "L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3",
  "L5 ET", "L5/6 NP", "L6 CT", "L6b")

build_enrichment_landscape <- function(neurons_only = TRUE, dense_labels = FALSE) {
  enrich   <- read_panel("panel_A_enrichment.csv")
  families <- read_panel("panel_A_family_groups.csv")

  enrich <- enrich |>
    mutate(neg_log10_p = pmin(-log10(pmax(p_value, 1e-300)), 50))

  # Optionally restrict the panel to neuronal supertypes only. Filtering is
  # done in display space; the FDR/Bonferroni thresholds remain those of the
  # original 503-type correction universe (defined just below).
  if (neurons_only) {
    enrich   <- enrich |> filter(family %in% NEURONAL_FAMILIES) |>
                  arrange(x) |> mutate(x = row_number() - 1L)
    families <- families |> filter(family %in% NEURONAL_FAMILIES)
    # Recompute family min/max in the new x-coordinate space.
    fam_ranges <- enrich |> group_by(family) |>
      summarise(min = min(x), max = max(x), size = n(), .groups = "drop")
    families <- families |> select(-min, -max, -size) |>
      left_join(fam_ranges, by = "family")
  }

  y_max <- max(enrich$neg_log10_p) * 1.30
  fam_y <- -(y_max * 0.04)
  fam_h <-  (y_max * 0.025)

  # Thresholds. FDR line at the smallest -log10(p) among BH-significant
  # types; Bonferroni line at -log10(0.05 / 503) over the correction universe.
  bonf_y <- -log10(0.05 / 503)
  fdr_y  <- if (any(enrich$passes_fdr_05))
              min(enrich$neg_log10_p[enrich$passes_fdr_05]) else NA_real_

  ggplot(enrich, aes(x = x, y = neg_log10_p, fill = color)) +
    geom_col(width = 0.85, color = "black", linewidth = 0.08) +
    scale_fill_identity() +
    geom_hline(yintercept = fdr_y,  color = "#888",   linetype = "dashed",  linewidth = 0.25) +
    geom_hline(yintercept = bonf_y, color = "#4B0082", linetype = "dotdash", linewidth = 0.3) +
    annotate("text", x = max(enrich$x) - 0.5, y = bonf_y + 0.5,
             label = "Bonferroni 0.05", hjust = 1, size = 2.0,
             color = "#4B0082", fontface = "bold") +
    annotate("text", x = max(enrich$x) - 0.5, y = fdr_y + 0.5,
             label = "FDR 0.05", hjust = 1, size = 2.0,
             color = "#555", fontface = "bold") +
    # Label every Bonferroni-significant supertype above its bar (90°-rotated).
    geom_text(data = filter(enrich, passes_bonf_503),
              aes(x = x, y = neg_log10_p + y_max * 0.015, label = supertype),
              angle = 90, hjust = 0, vjust = 0.5,
              size = 1.8, color = "black", inherit.aes = FALSE) +
    # Family colorbar at the bottom.
    geom_rect(data = families,
              aes(xmin = min - 0.4, xmax = max + 0.4,
                  ymin = fam_y, ymax = fam_y + fam_h,
                  fill = color),
              color = "black", linewidth = 0.1, inherit.aes = FALSE) +
    # Family labels. Default: wide groups (≥4 types) horizontal, narrow ones 45°.
    # dense_labels (all-137 panel): angle every group so the many narrow,
    # long-named non-neuronal families stay legible instead of over-printing.
    (if (dense_labels)
       geom_text(data = families,
                 aes(x = (min + max) / 2, y = fam_y - fam_h * 0.7, label = family),
                 inherit.aes = FALSE, size = 2.0, hjust = 1, vjust = 1, angle = 45)
     else list(
       geom_text(data = filter(families, size >= 4),
                 aes(x = (min + max) / 2, y = fam_y - fam_h * 0.6, label = family),
                 inherit.aes = FALSE, size = 2.8, hjust = 0.5, vjust = 1),
       geom_text(data = filter(families, size < 4),
                 aes(x = (min + max) / 2, y = fam_y - fam_h * 1.2, label = family),
                 inherit.aes = FALSE, size = 2.4, hjust = 1, vjust = 1, angle = 45))) +
    scale_x_continuous(expand = expansion(add = c(0.5, 0.5))) +
    coord_cartesian(ylim = c(fam_y - fam_h * (if (dense_labels) 3.4 else 2), y_max),
                    clip = "off") +
    labs(x = NULL, y = LAB_GWAS) +
    theme_panel() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          plot.margin = margin(5, 6, if (dense_labels) 16 else 10, 6))
}


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
    geom_text_repel(aes(label = supertype), size = 2.4,
                    box.padding = 0.25, point.padding = 0.20,
                    max.overlaps = 30, segment.color = "#999",
                    segment.size = 0.2, force = 3, min.segment.length = 0) +
    scale_fill_identity() +
    depletion_scales(show_legend = !compact) +
    (if (!compact) scale_size(range = c(1.2, 3), guide = "none")) +
    inset_spearman(rp$rho, rp$p, corner = "top-left", size = 2.4) +
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
# so C and E share a horizontal axis for cross-reference; same y-axis units
# as Panel A and Panel B so the genetics-axis is consistent across the row.
# Source CSVs are panel_B (provides SCZ enrichment) + panel_E (HCN1 expr);
# join on supertype in R rather than duplicating in the Python exporter.
build_hcn1_vs_scz_plot <- function() {
  panel_b <- read_panel("panel_B_genetics_vs_depletion.csv")
  panel_e <- read_panel("panel_E_hcn1_vs_sag.csv")

  df <- panel_b |>
    select(supertype, scz_neg_log10_p, depleted_fdr20, color, exemplar) |>
    inner_join(select(panel_e, supertype, HCN1_expr),
               by = "supertype") |>
    add_depletion_status() |>
    mutate(HCN1_expr_log2 = HCN1_expr / log(2))

  rp <- spearman_rp(df$HCN1_expr_log2, df$scz_neg_log10_p)

  ggplot(df, aes(x = scz_neg_log10_p, y = HCN1_expr_log2)) +
    geom_smooth(method = "lm", formula = y ~ x, color = "#444",
                fill = "#888", alpha = 0.18, linewidth = 0.3,
                linetype = "dashed", se = TRUE) +
    geom_point(aes(fill = color,
                   color = depleted_status, stroke = depleted_status),
               shape = 21, size = 2.6, alpha = 0.92) +
    geom_text_repel(aes(label = supertype), size = 2.4,
                    box.padding = 0.25, point.padding = 0.20,
                    max.overlaps = 30, segment.color = "#999",
                    segment.size = 0.2, force = 3, min.segment.length = 0) +
    scale_fill_identity() +
    # Depletion legend suppressed — Panel B carries it for the row.
    depletion_scales(show_legend = FALSE) +
    inset_spearman(rp$rho, rp$p, corner = "bottom-right", size = 2.4) +
    labs(x = LAB_GWAS, y = LAB_HCN1) +
    theme_panel()
}


# ====================================================================
# D) HCN1 locus zoom (Manhattan + gene track, hg38)
# ====================================================================
build_hcn1_locus_plot <- function() {
  snps    <- read_panel("panel_D_snps.csv")
  cs      <- read_panel("panel_D_credible_set.csv")
  genes   <- read_panel("panel_D_genes.csv")
  exons   <- read_panel("panel_D_exons.csv")
  meta    <- read_panel("panel_D_meta.csv")

  # --- top: Manhattan with FINEMAP credible set ---
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
      name    = "FINEMAP PIP") +
    scale_size(range = c(1.0, 3.2), guide = "none") +
    coord_cartesian(xlim = c(meta$win_lo_mb, meta$win_hi_mb)) +
    labs(x = NULL, y = expression("SCZ GWAS ("*-log[10]~italic(P)*")")) +
    theme_panel() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = c(0.98, 0.98),
          legend.justification = c("right", "top"),
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

  gene_track <- ggplot() +
    geom_segment(data = genes_with_y,
                 aes(x = start_mb, xend = stop_mb, y = y, yend = y,
                     color = factor(is_hcn1)),
                 linewidth = 0.5, alpha = 0.75) +
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
    geom_text_repel(aes(label = supertype), size = 2.4,
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
    inset_spearman(rp$rho, rp$p, corner = "bottom-right", size = 2.4) +
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
  x_gap   <- 80
  shifts  <- c(-(half_w + x_gap / 2), +(half_w + x_gap / 2))
  names(shifts) <- meta$supertype[order(meta$cell_order)]
  segs <- segs |> mutate(x0_sh = x0 + shifts[cell],
                         x1_sh = x1 + shifts[cell])
  meta <- meta |> mutate(soma_x_um_sh = soma_x_um + shifts[supertype])

  # Per-cell label-y: just above each cell's topmost dendrite/axon.
  cell_tops <- segs |> group_by(cell) |>
    summarise(top_y = min(pmin(y0, y1), na.rm = TRUE)) |>
    rename(supertype = cell)
  meta <- meta |> left_join(cell_tops, by = "supertype") |>
    mutate(label_y = top_y - 90)

  # Vertical zoom: from pia (with a bit of headroom for "Pia" label) to
  # just below the deepest dendrite + scale-bar room.
  max_y    <- max(segs$y0, segs$y1)
  y_bottom <- max_y + 320   # bar at +90, label at +180, plus padding
  y_top    <- -130
  x_extent <- c(shifts[1] - half_w * 1.40, shifts[2] + half_w * 1.15)

  # Layer-boundary lines (Pia is solid; L1/L2, L2/L3, L3/L4 are dashed).
  layer_lines <- tibble(
    y         = c(0,
                   LAYER_FRACTIONS$L1_L2,
                   LAYER_FRACTIONS$L2_L3,
                   LAYER_FRACTIONS$L3_L4) * c(1, rep(CORTEX_THICKNESS_UM, 3)),
    linewidth = c(0.25, 0.18, 0.18, 0.18),
    linetype  = c("solid", "dashed", "dashed", "dashed"))

  # Layer text positions (LEFT side of panel).
  layer_label_x <- x_extent[1] + diff(x_extent) * 0.02
  layer_labels  <- tibble(
    text = c("L1", "L2", "L3", "L4+"),
    y    = c( LAYER_FRACTIONS$L1_L2 / 2,
              (LAYER_FRACTIONS$L1_L2 + LAYER_FRACTIONS$L2_L3) / 2,
              (LAYER_FRACTIONS$L2_L3 + LAYER_FRACTIONS$L3_L4) / 2,
              (LAYER_FRACTIONS$L3_L4 + 1) / 2) * CORTEX_THICKNESS_UM)

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
                  label = paste(supertype, layer, sep = " · "),
                  color = supertype),
              size = 2.8, hjust = 0.5, vjust = 1) +
    # Layer text on the LEFT.
    geom_text(data = layer_labels,
              aes(x = layer_label_x, y = y, label = text),
              color = "#777", size = 2.8, hjust = 0, vjust = 0.5,
              fontface = "italic") +
    annotate("text", x = layer_label_x, y = -50, label = "Pia",
             color = "#555", size = 2.8, hjust = 0) +
    # 100 µm scale bar (label offset clearly below the bar so they don't overlap).
    annotate("segment",
             x    = x_extent[1] + diff(x_extent) * 0.05,
             xend = x_extent[1] + diff(x_extent) * 0.05 + 100,
             y    = max_y + 90, yend = max_y + 90,
             linewidth = 0.8) +
    annotate("text", x = x_extent[1] + diff(x_extent) * 0.05 + 50,
             y = max_y + 180, label = "100 µm",
             size = 2.4, hjust = 0.5, vjust = 0) +
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
          plot.margin = margin(8, 24, 8, 8))
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
  stim_y  <- v_range[1] - v_pad * 0.4

  # Per-trace supertype labels, placed mid-step where the two traces are
  # most separated. Sit ~5 mV above each cell's steady-state V — comfortably
  # above the trace, away from the sag-annotation glyphs at the front of
  # the step (V_peak_t_ms ≈ 50 ms for Sst_25).
  trace_labels <- meta |>
    mutate(label_t_ms = 750,
           label_v_mV = V_steady + 5)

  ggplot(traces, aes(x = t_ms, y = v_mV, color = cell)) +
    geom_line(linewidth = 0.45, alpha = 0.95) +
    geom_text(data = trace_labels,
              aes(x = label_t_ms, y = label_v_mV,
                  label = supertype, color = supertype),
              size = 2.6, hjust = 0.5, vjust = 0,
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
              color = sst25_color, fontface = "bold", size = 2.2,
              hjust = 0, inherit.aes = FALSE) +
    # Stim bar at the bottom.
    annotate("segment", x = 0, xend = 1000, y = stim_y, yend = stim_y,
             linewidth = 0.55) +
    annotate("segment", x = 0,    xend = 0,    y = stim_y, yend = stim_y + v_pad * 0.2,
             linewidth = 0.4) +
    annotate("segment", x = 1000, xend = 1000, y = stim_y, yend = stim_y + v_pad * 0.2,
             linewidth = 0.4) +
    scale_color_manual(values = setNames(meta$color, meta$supertype),
                       guide = "none") +
    coord_cartesian(ylim = c(v_range[1] - v_pad * 0.9, v_range[2] + v_pad * 0.1)) +
    labs(x = "Time (ms, zeroed at step onset)",
         y = "Membrane potential (mV)") +
    theme_panel() +
    theme(legend.position = "none")
}


# ====================================================================
# C (restored)) Sst_25 gene-driver scatter
# ====================================================================
# The original gene-driver panel: which genes drive Sst_25's SCZ enrichment.
# X = gene specificity in Sst_25 (log10); Y = -log10 SCZ MAGMA gene p. Drivers
# (Sst_25-specific AND GWAS-significant) in red, other genes grey, HCN1 called
# out in blue. Guide lines: Sst_25 specificity 90th-pct (vertical), MAGMA
# FDR 0.05 and genome-wide 5e-8 (horizontal). Data: panel_C_gene_drivers /
# _top_labels / _meta (from export_for_R.py). Used by the 7-panel variant.
GENOME_WIDE_NEG_LOG10P <- -log10(5e-8)   # 7.30

build_gene_driver_plot <- function() {
  g    <- read_panel("panel_C_gene_drivers.csv")
  lab  <- read_panel("panel_C_top_labels.csv")
  meta <- read_panel("panel_C_meta.csv")

  g_grey <- g |> filter(!is_driver, !is_hcn1)
  g_red  <- g |> filter(is_driver,  !is_hcn1)
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
    geom_text_repel(data = lab_g, aes(label = symbol), size = 2.1,
                    fontface = "italic", color = "#333",
                    box.padding = 0.30, max.overlaps = 40,
                    segment.color = "#bbb", segment.size = 0.2,
                    min.segment.length = 0, force = 4, seed = 1) +
    geom_text_repel(data = g_hcn1, aes(label = symbol), size = 3.3,
                    fontface = "bold.italic", color = HCN1_COLOR,
                    nudge_y = 3.2, nudge_x = 0.25, box.padding = 0.6,
                    segment.color = HCN1_COLOR, segment.size = 0.4,
                    min.segment.length = 0, seed = 1) +
    annotate("text", x = meta$xlim_hi, y = GENOME_WIDE_NEG_LOG10P,
             label = "'Genome-wide' ~ 5 %*% 10^-8", parse = TRUE,
             hjust = 1, vjust = -0.45, size = 1.9, color = "purple") +
    annotate("text", x = meta$xlim_hi, y = meta$fdr_neg_log10p_cutoff,
             label = "MAGMA FDR 0.05", hjust = 1, vjust = -0.45,
             size = 1.9, color = "#999") +
    scale_x_log10(limits = c(meta$xlim_lo, meta$xlim_hi),
                  breaks = c(0.001, 0.003, 0.01, 0.03)) +
    coord_cartesian(ylim = c(meta$ylim_lo, meta$ylim_hi), clip = "off") +
    labs(x = "Gene specificity in Sst_25 (log scale)",
         y = expression(-log[10]*"(SCZ MAGMA gene "*italic(p)*")")) +
    theme_panel()
}


# ====================================================================
# Assemble & save
# ====================================================================
enrichment_landscape    <- build_enrichment_landscape()
genetics_vs_depletion   <- build_genetics_vs_depletion()
hcn1_vs_scz_plot        <- build_hcn1_vs_scz_plot()
hcn1_locus_plot         <- build_hcn1_locus_plot()
hcn1_expression_vs_sag  <- build_hcn1_expression_vs_sag()
morphology_plot         <- build_morphology_plot()
ephys_traces_plot       <- build_ephys_traces_plot()

panel_label_kwargs <- list(
  label_size = PANEL_LABEL_SIZE, label_fontface = "bold",
  label_x = 0.0, label_y = 1.0,  hjust = -0.3, vjust = 1.3)

# Main figure: the former Panel A (enrichment landscape) is moved to a
# supplement; the remaining six panels relabel A–F across two rows.
row_ABC <- do.call(plot_grid,
  c(list(genetics_vs_depletion, hcn1_locus_plot, hcn1_vs_scz_plot,
         nrow = 1, rel_widths = c(1, 1, 1),
         labels = c("A", "B", "C")),
    panel_label_kwargs))

row_DEF <- do.call(plot_grid,
  c(list(hcn1_expression_vs_sag, morphology_plot, ephys_traces_plot,
         nrow = 1, rel_widths = c(1, 1, 1),
         labels = c("D", "E", "F")),
    panel_label_kwargs))

figure <- plot_grid(row_ABC, row_DEF, ncol = 1, rel_heights = c(1, 1))

MAIN_H <- 4.3   # two panel rows (former Panel A removed to the supplement)
ggsave(OUT_PNG, figure, width = FIG_WIDTH_IN, height = MAIN_H, dpi = 350, bg = "white")
ggsave(OUT_PDF, figure, width = FIG_WIDTH_IN, height = MAIN_H, bg = "white")

# Supplement: SCZ enrichment landscape (former Panel A), on its own.
ggsave(OUT_SUPPA_PNG, enrichment_landscape, width = FIG_WIDTH_IN, height = 2.8,
       dpi = 350, bg = "white")
ggsave(OUT_SUPPA_PDF, enrichment_landscape, width = FIG_WIDTH_IN, height = 2.8,
       bg = "white")
message(glue("\n→ Saved {OUT_PNG}"))
message(glue("→ Saved {OUT_PDF}"))
message(glue("→ Saved {OUT_SUPPA_PNG}  (former Panel A)"))


# ====================================================================
# Variant: 7-panel figure that KEEPS the enrichment landscape (A) and
# restores the Sst_25 gene-driver panel (C). Row 1 = landscape (full width);
# rows 2-3 = B depletion / C gene-driver / D locus  and  E sag / F morph /
# G traces. (The v3 six-panel figure above stays as-is.)
# ====================================================================
gene_driver_plot         <- build_gene_driver_plot()
# Panel A here shows ALL 137 supertypes (incl. non-neuronal), matching the
# original 7-panel layout — unlike the v3 supplement, which is neurons-only.
enrichment_landscape_all <- build_enrichment_landscape(neurons_only = FALSE,
                                                       dense_labels = TRUE)

row_A_wd <- plot_grid(enrichment_landscape_all,
  labels = "A", label_size = PANEL_LABEL_SIZE, label_fontface = "bold",
  label_x = 0.0, label_y = 1.0, hjust = -0.3, vjust = 1.3)

row_BCD_wd <- do.call(plot_grid,
  c(list(genetics_vs_depletion, gene_driver_plot, hcn1_locus_plot,
         nrow = 1, rel_widths = c(1, 1, 1), labels = c("B", "C", "D")),
    panel_label_kwargs))

row_EFG_wd <- do.call(plot_grid,
  c(list(hcn1_expression_vs_sag, morphology_plot, ephys_traces_plot,
         nrow = 1, rel_widths = c(1, 1, 1), labels = c("E", "F", "G")),
    panel_label_kwargs))

figure_wd <- plot_grid(row_A_wd, row_BCD_wd, row_EFG_wd,
                       ncol = 1, rel_heights = c(0.82, 1, 1))

OUT_WD_PNG <- file.path(REPO, "results", "figures", "scz_sst_hcn1_multipanel_withdriver.png")
OUT_WD_PDF <- file.path(REPO, "results", "figures", "scz_sst_hcn1_multipanel_withdriver.pdf")
WD_H <- 2.5 + MAIN_H   # landscape row + the two scatter rows
ggsave(OUT_WD_PNG, figure_wd, width = FIG_WIDTH_IN, height = WD_H, dpi = 350, bg = "white")
ggsave(OUT_WD_PDF, figure_wd, width = FIG_WIDTH_IN, height = WD_H, bg = "white")
message(glue("→ Saved {OUT_WD_PNG}  (7-panel, with gene-driver)"))


# ====================================================================
# Standalone: HCN1 expression vs SCZ Sst depletion (16 Sst supertypes)
# ====================================================================
# Simple Fig-4-style scatter: the more HCN1 a Sst supertype expresses, the
# more it is depleted in SCZ. x = mean HCN1 expression (log2 CP10K+1) from
# panel_E; y = compositional depletion (−β) from panel_B; joined on supertype.
build_hcn1_vs_depletion <- function() {
  panel_b <- read_panel("panel_B_genetics_vs_depletion.csv")
  panel_e <- read_panel("panel_E_hcn1_vs_sag.csv")

  df <- panel_b |>
    select(supertype, comp_beta, depleted_fdr20, color) |>
    inner_join(select(panel_e, supertype, HCN1_expr), by = "supertype") |>
    add_depletion_status() |>
    mutate(HCN1_expr_log2 = HCN1_expr / log(2),
           depletion      = -comp_beta)   # depleted supertypes point up

  rp <- spearman_rp(df$HCN1_expr_log2, df$depletion)

  ggplot(df, aes(x = HCN1_expr_log2, y = depletion)) +
    geom_smooth(method = "lm", formula = y ~ x, color = "#444",
                fill = "#888", alpha = 0.18, linewidth = 0.3,
                linetype = "dashed", se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "#aaa", linewidth = 0.2) +
    geom_point(aes(fill = color, color = depleted_status, stroke = depleted_status),
               shape = 21, size = 2.8, alpha = 0.92) +
    geom_text_repel(aes(label = supertype), size = 2.4,
                    box.padding = 0.25, point.padding = 0.20,
                    max.overlaps = 30, segment.color = "#999",
                    segment.size = 0.2, force = 3, min.segment.length = 0) +
    scale_fill_identity() +
    depletion_scales(show_legend = TRUE) +
    inset_spearman(rp$rho, rp$p, corner = "top-left", size = 2.6) +
    labs(x = LAB_HCN1, y = LAB_DEPL) +
    theme_panel() +
    theme(legend.position = c(0.98, 0.02),
          legend.justification = c("right", "bottom"))
}

hcn1_vs_depletion_plot <- build_hcn1_vs_depletion()
OUT_HD_PNG <- file.path(REPO, "results", "figures", "hcn1_vs_depletion_sst.png")
OUT_HD_PDF <- file.path(REPO, "results", "figures", "hcn1_vs_depletion_sst.pdf")
ggsave(OUT_HD_PNG, hcn1_vs_depletion_plot, width = 3.4, height = 3.1, dpi = 350, bg = "white")
ggsave(OUT_HD_PDF, hcn1_vs_depletion_plot, width = 3.4, height = 3.1, bg = "white")
message(glue("→ Saved {OUT_HD_PNG}  (HCN1 expr vs Sst depletion)"))


# ====================================================================
# Variant: compact Fig 4 — no enrichment-landscape row; the two scatter
# panels (depletion↔genetics, sag↔HCN1) get uniform dot size and no legend;
# the HCN1-expression-vs-depletion scatter is added as the last column of
# the top row. Bottom row (sag / morphology / traces) unchanged.
# ====================================================================
genetics_vs_depletion_c <- build_genetics_vs_depletion(compact = TRUE)
sag_c                   <- build_hcn1_expression_vs_sag(compact = TRUE)

# Top row = the genetic case (cell-type → gene → variant); bottom row = HCN1's
# functional consequences (depletion, sag) grounded in exemplar cells. D and E
# sit adjacent and share the HCN1-expression x-axis.
row_top_c <- do.call(plot_grid,
  c(list(genetics_vs_depletion_c, gene_driver_plot, hcn1_locus_plot,
         nrow = 1, rel_widths = c(1, 1, 1),
         labels = c("A", "B", "C")),
    panel_label_kwargs))

row_bot_c <- do.call(plot_grid,
  c(list(hcn1_vs_depletion_plot, sag_c, morphology_plot, ephys_traces_plot,
         nrow = 1, rel_widths = c(1, 1, 1, 1),
         labels = c("D", "E", "F", "G")),
    panel_label_kwargs))

figure_c <- plot_grid(row_top_c, row_bot_c, ncol = 1, rel_heights = c(1, 1))

OUT_C_PNG <- file.path(REPO, "results", "figures", "scz_sst_hcn1_multipanel_compact.png")
OUT_C_PDF <- file.path(REPO, "results", "figures", "scz_sst_hcn1_multipanel_compact.pdf")
C_W <- 8.0    # slightly wider than the 6.5" default: the top row now has 4 panels
C_H <- 4.5
ggsave(OUT_C_PNG, figure_c, width = C_W, height = C_H, dpi = 350, bg = "white")
ggsave(OUT_C_PDF, figure_c, width = C_W, height = C_H, bg = "white")
message(glue("→ Saved {OUT_C_PNG}  (compact: no landscape, +HCN1-vs-depletion)"))
