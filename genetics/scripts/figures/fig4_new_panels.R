# Builder functions for the Figure 4 marker/AD panels.
#
# Kept separate from the render script so that BOTH the standalone panel
# renderer (render_fig4_new_panels.R) and the multipanel assembly
# (scz_sst_hcn1_story.R) build these panels from one definition.
#
# Panel data: genetics/scripts/figures/export_fig4_new_panels.py

suppressPackageStartupMessages({
  library(ggplot2); library(cowplot); library(ggrepel)
  library(dplyr); library(readr)
})

if (!exists("FIG4_REPO")) FIG4_REPO <- "/Users/shreejoy/Github/scz_celltype_paper"
# Follow the caller's DATA_DIR when there is one. This was previously a
# hardcoded absolute path, which silently kept panels g-i on an old snapshot
# while a-f moved to a new one -- the two halves of the figure disagreed and
# nothing errored.
FIG4_PANELS <- if (exists("DATA_DIR")) DATA_DIR else
  file.path(FIG4_REPO, "genetics", "results", "figures", "r_panels")
source(file.path(REPO, "scripts", "figures", "fig4_style.R"))

# Distinct name: the assembly script has its own read_panel() bound to DATA_DIR.
read_new_panel <- function(f) read_csv(file.path(FIG4_PANELS, f), show_col_types = FALSE)

# =====================================================================
# Volcano - normative marker genes of the SCZ-depleted Sst supertypes
# =====================================================================
# level = "donor" (default) puts the donor-paired P on y -- the design that
# guards against pseudoreplication across only 5 donors.
# level = "cell" shows the raw cell-level Wilcoxon instead: 12,393 cells make
# the p-values astronomically small, 1,091 genes underflow to exactly 0, and
# 63% of genes clear FDR 0.05 -- so the axis saturates and both CALB1 and HCN1
# land on the clipped ceiling. Provided for comparison.
# Effect-size gate for the volcano, in log2 fold change (avg_log2FC computed
# Seurat-style on the un-logged CP10K scale). Loosened from 0.5 to 0.25 so that
# HCN1 (log2FC = +0.434) is called differentially expressed; at 0.5 it fell just
# short despite an FDR of 2e-138. Cost: 183 -> 580 called genes.
LOG2FC_GATE <- 0.25
N_DRIVER_LABELS <- 4

build_marker_volcano <- function(level = c("donor", "cell")) {
  level <- match.arg(level)
  v <- read_new_panel("panel_volcano_vulnerable_vs_notdepleted.csv")

  # Underflowed p-values (exactly 0) are clipped to the smallest observed
  # non-zero p and drawn as triangles, so the ceiling is not mistaken for data.
  cell_floor <- min(v$p_wilcoxon[v$p_wilcoxon > 0], na.rm = TRUE)

  v <- v |>
    mutate(
      clipped = if (level == "cell") p_wilcoxon <= 0 else FALSE,
      donor_logp = if (level == "cell")
                     -log10(pmax(p_wilcoxon, cell_floor))
                   else -log10(pmax(p_donor_paired, 1e-300)),
      # Both axes thresholded, as on a conventional volcano: donor-level BH
      # FDR < 0.10 (the test that respects the 5 donors) AND the effect-size
      # gate |log2FC| > 0.5. FDR alone would colour ~2,900 genes, many sitting
      # at log2FC ~ 0 -- consistent across donors but not markers.
      sig = if (level == "cell") fdr_wilcoxon < 0.05 & abs(log2FC) > LOG2FC_GATE
            else fdr_donor < 0.10 & abs(log2FC) > LOG2FC_GATE,
      score = abs(log2FC) * donor_logp)

  v_ns   <- filter(v, !sig, !gene %in% c("CALB1", "HCN1"))
  v_up   <- filter(v,  sig, log2FC > 0, !gene %in% c("CALB1", "HCN1"))
  v_dn   <- filter(v,  sig, log2FC < 0, !gene %in% c("CALB1", "HCN1"))
  v_key  <- filter(v, gene %in% c("CALB1", "HCN1")) |>
              mutate(key_color = ifelse(gene == "HCN1", HCN1_COLOR, CALB1_COLOR))
  # Was 9 per direction (18 driver labels). Cut to 4 so HCN1 and CALB1 read as
  # the callouts rather than as two names among twenty.
  lab    <- bind_rows(slice_max(v_up, score, n = N_DRIVER_LABELS),
                      slice_max(v_dn, score, n = N_DRIVER_LABELS))

  # All labels in one repel pass (separate passes cannot see each other's boxes).
  # HCN1/CALB1 are additionally nudged off the dense significant cloud into
  # adjacent empty space, and keep a leader line back to their point.
  lab_all <- bind_rows(
    transmute(lab, gene, log2FC, donor_logp, seg = "#bbb",
              lab_color = "#333", lab_face = "italic", lab_size = LBL_GENE,
              nx = 0, ny = 0),
    transmute(v_key, gene, log2FC, donor_logp, seg = key_color,
              lab_color = key_color, lab_face = "bold.italic", lab_size = LBL_CALL,
              nx = ifelse(gene == "HCN1", -0.62,  0.72),
              ny = ifelse(gene == "HCN1",  0.85, -0.80)))

  ggplot(mapping = aes(log2FC, donor_logp)) +
    geom_vline(xintercept = c(-LOG2FC_GATE, LOG2FC_GATE), linetype = "dashed",
               color = "#c4c4c4", linewidth = 0.25) +
    geom_vline(xintercept = 0, linetype = "dotted",
               color = "#e0e0e0", linewidth = 0.2) +
    geom_hline(yintercept = min(v$donor_logp[v$sig]),
               linetype = "dashed", color = "#c4c4c4", linewidth = 0.25) +
    geom_point(data = v_ns, color = COL_NS,      size = 0.35, alpha = 0.5) +
    geom_point(data = v_dn, aes(shape = clipped), color = COL_UP_NOTD,
               size = 0.7, alpha = 0.85) +
    geom_point(data = v_up, aes(shape = clipped), color = COL_UP_VULN,
               size = 0.7, alpha = 0.85) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    geom_point(data = v_key, aes(fill = key_color), color = "white",
               shape = 21, size = 2.8, stroke = 0.3) +
    scale_fill_identity() +
    # One repel pass over ALL labels: separate passes cannot see each other's
    # boxes, which let the top-driver labels collide with the HCN1/CALB1 callouts.
    geom_text_repel(
      data = lab_all,
      aes(label = gene, color = lab_color, fontface = lab_face, size = lab_size,
          segment.colour = seg),
      nudge_x = lab_all$nx, nudge_y = lab_all$ny,
      box.padding = 0.35, max.overlaps = 40,
      segment.size = 0.25, min.segment.length = 0, force = 4, seed = 1,
      # White halo: labels sit over the dense significant-gene cloud, and the
      # HCN1/CALB1 callouts in particular were hard to read against it.
      # bg.r is relative to font size, so the bold callouts get a wider halo.
      bg.color = "white", bg.r = 0.15) +
    scale_color_identity() +
    scale_size_identity() +
    annotate("text", x = -Inf, y = min(v$donor_logp[v$sig]),
             label = if (level == "cell") "FDR == 0.05" else "FDR == 0.10",
             parse = TRUE, hjust = -0.1, vjust = -0.45, size = LBL_STAT,
             color = "#999") +
    labs(x = expression(log[2]~"fold change (depleted / not-depleted Sst)"),
         y = expression(-log[10]~italic(P))) +
    theme_panel()
}

# =====================================================================
# CALB1 violin - fill carries the group, outline carries depletion status
# =====================================================================
build_violin <- function(gene_name = "CALB1") {
  d <- read_new_panel("panel_violin_calb1_percell.csv") |>
    filter(gene == gene_name) |>
    mutate(group = factor(ifelse(group == "Vulnerable", "Depleted", group),
                         levels = c("Not depleted", "Depleted")))
  st <- read_new_panel("panel_violin_stats.csv") |> filter(gene == gene_name)

  # Report the cell-level BH FDR, which is the statistic the adjacent volcano
  # uses to define its gene set. Annotating this panel with the donor-paired
  # test instead left the two panels reporting different statistics for the
  # same comparison. fdr_donor is the fallback for older exports that predate
  # the fdr_cell column.
  q <- if ("fdr_cell" %in% names(st)) st$fdr_cell[1] else st$fdr_donor[1]
  plab <- if (q < 0.001) "FDR < 0.001" else paste0("FDR == ", signif(q, 2))

  ggplot(d, aes(group, expression, fill = group,
                color = group, linewidth = group)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.10, fill = NA, colour = "grey20",
                 outlier.shape = NA, linewidth = 0.35) +
    annotate("text", x = 1.5, y = max(d$expression) * 1.06,
             label = plab, parse = TRUE, size = LBL_STAT, color = "#222") +
    scale_fill_manual(values = c("Not depleted" = COL_UP_NOTD,
                                 "Depleted"     = COL_UP_VULN)) +
    # same outline encoding as the supertype scatters
    scale_color_manual(values = c("Not depleted" = unname(DEPLETION_OUTLINE[2]),
                                  "Depleted"     = unname(DEPLETION_OUTLINE[1]))) +
    # Same outline *encoding* as the scatters, but the point-stroke values
    # (0.7/0.2) read far too heavy at violin scale, so scale them down.
    scale_discrete_manual("linewidth",
                          values = c("Not depleted" = 0.15, "Depleted" = 0.4)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    # Wrapped so the two category labels do not collide in this narrow panel.
    scale_x_discrete(labels = c("Not depleted" = "Not\ndepleted",
                                "Depleted"     = "Depleted")) +
    labs(x = NULL, y = bquote(italic(.(gene_name))*" expression ("*log[2]*" CP10K+1)")) +
    theme_panel() +
    theme(legend.position = "none")
}

# =====================================================================
# SCZ vs AD abundance change across the 16 Sst supertypes
# =====================================================================
build_ad_concordance <- function(show_legend = TRUE) {
  # SEA-AD supertype colours + depletion flags come from panel_B, joined here
  # rather than duplicated in the Python exporter (same pattern as panel C/E).
  pal <- read_new_panel("panel_B_genetics_vs_depletion.csv") |>
    select(supertype, color, depleted_fdr20)
  d <- read_new_panel("panel_ad_concordance_sst.csv") |>
    select(-any_of(c("group"))) |>
    left_join(pal, by = "supertype") |>
    add_depletion_status() |>
    # Plot depletion (−β) on both axes, as panels a and d do, so that depleted
    # types point up/right throughout the figure. Flipping both axes leaves the
    # correlation unchanged.
    mutate(scz_depletion = -scz_beta, ad_depletion = -ad_slope)
  rp <- spearman_rp(d$scz_depletion, d$ad_depletion)

  ggplot(d, aes(scz_depletion, ad_depletion)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "#aaa", linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "#aaa", linewidth = 0.2) +
    trend_line() +
    supertype_points() +
    supertype_labels() +
    scale_fill_identity() +
    depletion_scales(show_legend = show_legend) +
    inset_spearman(rp$rho, rp$p, corner = "bottom-right") +
    # x is the same quantity as the y-axis of panels a and d -- share its range.
    coord_cartesian(xlim = DEPLETION_LIM) +
    # Titles kept short: a long rotated y-title runs the full panel height and
    # collides with the bold panel label in the assembled multipanel.
    labs(x = LAB_DEPL,
         # Region (DLPFC) is stated in the figure legend, not the axis, and this
         # keeps the wording symmetric with the x-axis ("...in SCZ").
         y = expression("Cell depletion in AD ("*-beta*"/SD CPS)")) +
    theme_panel() +
    (if (show_legend)
       theme(legend.position = c(0.02, 0.98),
             legend.justification = c("left", "top"),
             legend.background = element_blank())
     else theme(legend.position = "none"))
}
