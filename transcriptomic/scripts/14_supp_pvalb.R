#!/usr/bin/env Rscript
# ============================================================================
# 14_supp_pvalb.R — SUPPLEMENTARY figure: PVALB mRNA in Pvalb cells (parvalbumin
# interneurons) in SCZ, the three cross-platform views, mirroring how the
# composite (scripts/09) now treats SST and FGFR3:
#   a  Forest    : 7 snRNA-seq cohorts + pooled meta diamond + Xenium replication
#   b  Boxplot   : per-donor library-normalised Xenium expression (CP1K) + edgeR p
#   c  Exemplars : one Control + one SCZ Pvalb cell (boundary, dashed nucleus,
#                  PVALB transcript molecules) at the pooled group-median grain density
#
# PVALB was shown in the composite earlier but moved here so the main figure
# features SST (interneuron) and FGFR3 (astrocyte); this supplement preserves the
# full PVALB story. Builders are copied from scripts/09_composite_figure.R.
#
# Inputs (same as the composite):
#   data/DE_genes_all_cells_scz.csv          meta-analytic snRNA-seq DE
#   data/meta_results_cohorts_subclass.csv   per-cohort snRNA-seq DE (7 cohorts)
#   ../spatial/output/de/de_results_subclass.csv   Xenium spatial DE
#   results/tables/marker_norm_expr*.csv     per-donor CP1K + edgeR p  (scripts/12)
#   results/tables/exemplar_PVALB_*.csv      exemplar cell inputs       (scripts/10)
# Output: results/figures/S_pvalb.{png,pdf}
# ============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(cowplot); library(metafor); library(ggsignif)
})

GENE <- "PVALB"; CELL <- "Pvalb"
BASE <- 9                                   # standalone supplement -> a bit larger than the composite's 7
DOWN_DARK <- "#0072B2"; UP_DARK <- "#D55E00"
COL_POOL <- "black"; COL_XEN <- "#117733"; COL_COHORT <- "grey35"
NEXPR_COL <- c(Control = DOWN_DARK, SCZ = UP_DARK); DOT_COL <- "#B2182B"
TAB <- "results/tables"

ast    <- function(fdr) dplyr::case_when(
  is.na(fdr) ~ "", fdr < 0.01 ~ "***", fdr < 0.05 ~ "**", fdr < 0.10 ~ "*", TRUE ~ "")
is_dot <- function(p, fdr) !is.na(p) & p < 0.05 & (is.na(fdr) | fdr >= 0.10)

# ---- data (same derivation as the composite) ------------------------------
meta_tbl <- read_csv("data/DE_genes_all_cells_scz.csv", show_col_types = FALSE)
ct_map <- c("Astrocyte"="Astro","L2/3 IT"="L2_3 IT","L5/6 NP"="L5_6 NP",
            "Microglia-PVM"="Micro-PVM","Oligodendrocyte"="Oligo","Endothelial"="Endo")
coh <- read_csv("data/meta_results_cohorts_subclass.csv", show_col_types = FALSE) |>
  filter(genes == GENE, !is.na(logFC), !is.na(t), t != 0) |> mutate(SE = logFC / t)
xen_all <- read_csv("../spatial/output/de/de_results_subclass.csv", show_col_types = FALSE) |>
  mutate(cell_type = ifelse(celltype %in% names(ct_map), ct_map[celltype], celltype))
xen_forest <- xen_all |> filter(gene == GENE, !is.na(logFC), !is.na(F), F > 0) |>
  mutate(t = sign(logFC) * sqrt(F), SE = logFC / t) |>
  transmute(genes = gene, logFC, t, P.Value = PValue, adj.P.Val = FDR,
            cell_type, cohort = "Xenium", SE)
coh <- bind_rows(coh |> select(genes, logFC, t, P.Value, adj.P.Val, cell_type, cohort, SE),
                 xen_forest)

# ---- a. forest (copied from scripts/09 build_forest) ----------------------
build_forest <- function(gene_sym, cell, ttl, show_xlab = FALSE) {
  sub    <- coh |> filter(genes == gene_sym, cell_type == cell)
  snrna  <- sub |> filter(cohort != "Xenium")
  xen_in <- sub |> filter(cohort == "Xenium")
  if (nrow(snrna) == 0) return(NULL)
  re <- tryCatch(rma(yi = snrna$logFC, sei = snrna$SE, method = "DL",
                     control = list(maxiter = 500)), error = function(e) NULL)
  cohort_df <- snrna |>
    transmute(cohort, est = logFC, lo = logFC-1.96*SE, hi = logFC+1.96*SE,
              p = P.Value, padj = adj.P.Val,
              sig = ast(adj.P.Val), dot = is_dot(P.Value, adj.P.Val)) |>
    arrange(desc(est)) |>
    mutate(y = rev(seq_along(est)) + 3, role = "cohort", colour = COL_COHORT)
  meta_fdr <- {r <- meta_tbl |> filter(genes == gene_sym, cell_type == cell)
               if (nrow(r) == 1) r$padj else NA_real_}
  pool <- if (!is.null(re)) tibble(
    cohort = "meta", est = as.numeric(re$beta), lo = re$ci.lb, hi = re$ci.ub,
    p = re$pval, padj = meta_fdr, sig = ast(meta_fdr), dot = FALSE,
    y = 3, role = "pool", colour = COL_POOL) else NULL
  has_xen <- nrow(xen_in) > 0
  xen_row <- if (has_xen) { xr <- xen_in[1,]; tibble(
    cohort = "Xenium", est = xr$logFC, lo = xr$logFC-1.96*xr$SE,
    hi = xr$logFC+1.96*xr$SE, p = xr$P.Value, padj = xr$adj.P.Val,
    sig = ast(xr$adj.P.Val), dot = is_dot(xr$P.Value, xr$adj.P.Val),
    y = 1, role = "xenium", colour = COL_XEN) } else NULL
  pd <- bind_rows(cohort_df, pool, xen_row)
  data_lo <- min(pd$lo, na.rm=TRUE); data_hi <- max(pd$hi, na.rm=TRUE)
  rng <- data_hi - data_lo; off <- rng*0.03
  pd <- pd |> mutate(ast_x = ifelse(est>=0, hi+off, lo-off),
                     ast_h = ifelse(est>=0, 0, 1),
                     has_m = sig != "" | dot | (role=="pool" & sig==""))
  rt <- any(pd$has_m & pd$ast_h==0); lf <- any(pd$has_m & pd$ast_h==1)
  xlo <- data_lo - rng*(if (lf) 0.22 else 0.05)
  xhi <- data_hi + rng*(if (rt) 0.22 else 0.05)
  xlo <- min(xlo, -off); xhi <- max(xhi, off)
  ylab <- pd |> select(y, cohort) |> arrange(y)
  if (has_xen) ylab <- bind_rows(ylab, tibble(y = 2, cohort = ""))
  ggplot(pd, aes(est, y, colour = I(colour))) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.25) +
    {if (has_xen) geom_hline(yintercept = 2, colour = "grey85", linetype = "dotted", linewidth = 0.3)} +
    geom_segment(aes(x = lo, xend = hi, yend = y), linewidth = 0.45) +
    geom_point(data = cohort_df, aes(est, y), shape = 15, size = 1.4) +
    {if (!is.null(pool)) geom_point(data = pool, aes(est, y), shape = 23, size = 2.4, fill = COL_POOL, colour = COL_POOL)} +
    {if (has_xen) geom_point(data = xen_row, aes(est, y), shape = 17, size = 2.1, colour = COL_XEN)} +
    geom_text(data = pd, aes(label = sig, x = ast_x, hjust = ast_h),
              size = BASE*0.34, vjust = 0.75, colour = "grey15", fontface = "bold") +
    geom_point(data = dplyr::filter(pd, dot), aes(ast_x, y), shape = 16, size = 0.9, colour = "grey15", inherit.aes = FALSE) +
    geom_text(data = dplyr::filter(pd, role=="pool" & sig==""),
              aes(ast_x, y, hjust = ast_h), label = "n.s.", size = BASE*0.28,
              vjust = 0.5, colour = "grey45", fontface = "italic", inherit.aes = FALSE) +
    scale_y_continuous(breaks = ylab$y, labels = ylab$cohort, expand = expansion(add = 0.6)) +
    scale_x_continuous(limits = c(xlo, xhi), breaks = scales::pretty_breaks(3), expand = expansion(mult = 0.01)) +
    labs(x = if (show_xlab) expression("SCZ log"[2]~"FC") else NULL, y = NULL, subtitle = ttl) +
    theme_cowplot(font_size = BASE) +
    theme(plot.subtitle = element_text(size = BASE, face = "plain", margin = margin(b = 1)),
          axis.text  = element_text(size = BASE - 2),
          axis.title.x = element_text(size = BASE - 0.5, margin = margin(t = 1.5)),
          axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          plot.margin = margin(2, 6, 2, 2))
}

# ---- b. CP1K boxplot (copied from scripts/09 build_normexpr) ---------------
NEXPR <- read_csv(sprintf("%s/marker_norm_expr.csv", TAB), show_col_types = FALSE)
NEXPR$dx <- factor(NEXPR$dx, levels = c("Control", "SCZ"))
NSTAT <- read_csv(sprintf("%s/marker_norm_expr_stats.csv", TAB), show_col_types = FALSE)
build_normexpr <- function(gene, title) {
  df <- NEXPR[NEXPR$gene == gene, ]; p <- NSTAT$p[NSTAT$gene == gene]; yr <- range(df$cp1k)
  ggplot(df, aes(dx, cp1k)) +
    geom_boxplot(aes(fill = dx), width = 0.62, outlier.shape = NA, alpha = 0.55, linewidth = 0.4) +
    geom_jitter(aes(fill = dx), shape = 21, colour = "grey25", size = 1.7, stroke = 0.3,
                width = 0.13, height = 0, alpha = 0.9) +
    geom_signif(comparisons = list(c("Control", "SCZ")),
                annotations = sprintf("italic(p)=='%.3f'", p), parse = TRUE,
                y_position = yr[2] + 0.06 * diff(yr), tip_length = 0.02, textsize = BASE * 0.32, vjust = -0.1) +
    scale_fill_manual(values = NEXPR_COL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22)), n.breaks = 4) +
    labs(x = NULL, y = "counts / 1,000 tx", subtitle = title) +
    theme_cowplot(font_size = BASE) +
    theme(legend.position = "none",
          plot.subtitle = element_text(size = BASE - 1, face = "plain", hjust = 0.5, margin = margin(b = 1)),
          axis.title.y = element_text(size = BASE - 1), axis.text = element_text(size = BASE - 1.5),
          plot.margin = margin(2, 3, 2, 3))
}

# ---- c. exemplar cells (copied from scripts/09 build_exemplar) -------------
build_exemplar <- function(gene, dx, lim, scalebar_lab = FALSE) {
  bd <- read_csv(sprintf("%s/exemplar_%s_%s_boundary.csv", TAB, gene, dx), show_col_types = FALSE)
  nu <- read_csv(sprintf("%s/exemplar_%s_%s_nucleus.csv",  TAB, gene, dx), show_col_types = FALSE)
  dt <- read_csv(sprintf("%s/exemplar_%s_%s_dots.csv",     TAB, gene, dx), show_col_types = FALSE)
  nd <- nrow(dt); sb <- 5
  ggplot() +
    geom_polygon(data = bd, aes(x, y), fill = "grey93", colour = "grey45", linewidth = 0.3) +
    geom_polygon(data = nu, aes(x, y), fill = NA, colour = "grey30", linewidth = 0.25, linetype = "22") +
    {if (nd > 0) geom_point(data = dt, aes(x, y), colour = DOT_COL, size = 0.6, alpha = 0.85)} +
    annotate("text", x = -lim*0.98, y = lim*0.98, label = nd, hjust = 0, vjust = 1,
             size = BASE*0.34, colour = DOT_COL, fontface = "bold") +
    annotate("segment", x = lim*0.96 - sb, xend = lim*0.96, y = -lim*0.94, yend = -lim*0.94, linewidth = 0.7) +
    {if (scalebar_lab) annotate("text", x = lim*0.96 - sb/2, y = -lim*0.80, label = "5~mu*m",
                                parse = TRUE, size = BASE*0.28, vjust = 1)} +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    theme_void(base_size = BASE) + theme(plot.margin = margin(1, 2, 1, 2))
}
ex_hdr <- function(t) ggdraw() + draw_label(t, size = BASE - 0.5)
ex_lim <- 1.15 * max(vapply(c("Control", "SCZ"), function(dx) {
  bd <- read_csv(sprintf("%s/exemplar_%s_%s_boundary.csv", TAB, GENE, dx), show_col_types = FALSE)
  max(abs(c(bd$x, bd$y))) }, numeric(1)))

# ---- assemble -------------------------------------------------------------
cat("Building PVALB supplement...\n")
pForest <- build_forest(GENE, CELL, "PVALB / Pvalb", show_xlab = TRUE)
pBox    <- build_normexpr(GENE, "PVALB mRNA in Pvalb cells")
exC <- build_exemplar(GENE, "Control", ex_lim)
exS <- build_exemplar(GENE, "SCZ", ex_lim, scalebar_lab = TRUE)
ex_panel <- plot_grid(plot_grid(ex_hdr("Control"), ex_hdr("SCZ"), ncol = 2),
                      plot_grid(exC, exS, ncol = 2),
                      ncol = 1, rel_heights = c(0.13, 1))
row   <- plot_grid(pForest, pBox, ex_panel, ncol = 3, rel_widths = c(1.05, 0.72, 1.15),
                   labels = c("a", "b", "c"), label_size = BASE + 2, label_fontface = "bold")
title <- ggdraw() + draw_label(
  "Supplementary | PVALB mRNA in Pvalb cells (parvalbumin interneurons) in schizophrenia",
  size = BASE, fontface = "plain", x = 0.5)
final <- plot_grid(title, row, ncol = 1, rel_heights = c(0.10, 1))

ggsave("results/figures/S_pvalb.png", final, width = 8.2, height = 2.9, dpi = 300, bg = "white")
ggsave("results/figures/S_pvalb.pdf", final, width = 8.2, height = 2.9, bg = "white")
cat("Saved results/figures/S_pvalb.{png,pdf}\n")
