#!/usr/bin/env Rscript
# ============================================================================
# 09_composite_figure.R — Publication composite (7.1 x 6.7 in, 400 dpi)
# ============================================================================
# Assembles, in one figure, the cross-platform SCZ DE story:
#   a  Butterfly: up/down DE counts per cell type, FDR<0.10 (light) vs
#      FDR<0.05 (dark overlay)
#   b-e Volcanoes: Sst (b), L2/3 IT (c), Astro (d), Micro-PVM (e); select genes
#      labelled; each on tight, data-driven axes (limits NOT shared across panels)
#   f-i Forest plots for 4 (gene, cell) pairs (SST/Sst, BDNF/L2_3 IT,
#      FGFR3/Astro, FKBP5/OPC): 7 snRNA-seq cohorts + pooled meta + Xenium repl.
#   j  Concordance scatter: snRNA-seq meta logFC vs Xenium logFC
#   k  Library-normalised expression (CP1K, counts/1,000 tx) per donor, Ctrl vs
#      SCZ, for SST-in-Sst and FGFR3-in-Astro, titled per row, edgeR p (scripts/12)
#   l  Xenium exemplar cells: cell boundary + dashed nucleus + marker molecules
#      for one Control and one SCZ cell each for SST (Sst) and FGFR3 (Astro).
#      Each exemplar is the representative cell at the pooled group-median grain
#      density (size-matched, typical eccentricity); see scripts/10.
#
# DATA PROVENANCE (see notes/figures_crossplatform_validation.md for detail):
#   data/DE_genes_all_cells_scz.csv          meta-analytic snRNA-seq DE
#   data/meta_results_cohorts_subclass.csv   per-cohort snRNA-seq DE (7 cohorts)
#   ../spatial/output/de/de_results_subclass.csv   Xenium spatial DE (symlink to
#      ~/Github/SCZ_Xenium/output/de/; set INPUT_XENIUM if it moves)
#   results/tables/marker_norm_expr.csv      panel k input — scripts/12
#   results/tables/exemplar_*.csv            panel l inputs — produced by
#      scripts/10_xenium_exemplar_cells.py (run that FIRST; it reads the
#      Xenium h5ad + boundary + transcript exports).
#
# Output: results/09_composite.{png,pdf}   (7.1 x 6.7 in)
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(cowplot); library(ggrepel); library(metafor); library(ggsignif)
})

INPUT_META   <- "data/DE_genes_all_cells_scz.csv"
INPUT_COH    <- "data/meta_results_cohorts_subclass.csv"
INPUT_XENIUM <- "../spatial/output/de/de_results_subclass.csv"  # internal cross-ref (was ~/Github/SCZ_Xenium/...)

FIG_W <- 7.1    # max total width (inches)

# ---- palette / groupings ---------------------------------------------------
UP_DARK   <- "#D55E00"; UP_LIGHT   <- "#F2B58C"
DOWN_DARK <- "#0072B2"; DOWN_LIGHT <- "#9FCAE6"
COL_NS    <- "grey75"
COL_POOL  <- "black";   COL_XEN    <- "#117733"; COL_COHORT <- "grey35"
EXC <- c("L2_3 IT","L4 IT","L5 IT","L5 ET","L5_6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")
INH <- c("Lamp5","Pax6","Pvalb","Sncg","Sst","Sst Chodl","Vip","Chandelier")
GLI <- c("Astro","Oligo","OPC","Micro-PVM","Endo","VLMC")
class_of <- function(ct) dplyr::case_when(
  ct %in% EXC ~ "Excitatory", ct %in% INH ~ "Inhibitory",
  ct %in% GLI ~ "Glia",       TRUE        ~ NA_character_)
CLASS_COL <- c(Excitatory = "#117733", Inhibitory = "#882255", Glia = "#DDCC77")

BASE <- 7   # base font size (Nature/NN: all figure text 5–7 pt). Most text is
            # sized relative to this; geom_text multipliers tuned so nothing
            # exceeds 7 pt or drops below 5 pt at the 7.1 × 6.7 in print size.

# significance helpers (shared)
ast    <- function(fdr) dplyr::case_when(
  is.na(fdr) ~ "", fdr < 0.01 ~ "***", fdr < 0.05 ~ "**", fdr < 0.10 ~ "*", TRUE ~ "")
is_dot <- function(p, fdr) !is.na(p) & p < 0.05 & (is.na(fdr) | fdr >= 0.10)

# ============================================================================
# Load data
# ============================================================================
cat("Loading meta-analytic DE...\n")
meta_tbl <- read_csv(INPUT_META, show_col_types = FALSE)

# genes/cells needed for the 6 forest panels
FOREST <- tibble::tribble(
  ~gene,    ~cell,        ~lab,
  "SST",    "Sst",        "SST / Sst",
  "BDNF",   "L2_3 IT",    "BDNF / L2/3 IT",
  "FGFR3",  "Astro",      "FGFR3 / Astro",
  "FKBP5",  "OPC",        "FKBP5 / OPC"
)

cat("Loading per-cohort table (large) and filtering to forest genes...\n")
coh <- read_csv(INPUT_COH, show_col_types = FALSE) |>
  filter(genes %in% FOREST$gene, !is.na(logFC), !is.na(t), t != 0) |>
  mutate(SE = logFC / t)

cat("Loading + harmonizing Xenium DE...\n")
ct_map <- c("Astrocyte"="Astro","L2/3 IT"="L2_3 IT","L5/6 NP"="L5_6 NP",
            "Microglia-PVM"="Micro-PVM","Oligodendrocyte"="Oligo","Endothelial"="Endo")
xen_all <- read_csv(INPUT_XENIUM, show_col_types = FALSE) |>
  mutate(cell_type = ifelse(celltype %in% names(ct_map), ct_map[celltype], celltype))
# Xenium rows for the forest genes, in cohort schema
xen_forest <- xen_all |>
  filter(gene %in% FOREST$gene, !is.na(logFC), !is.na(F), F > 0) |>
  mutate(t = sign(logFC) * sqrt(F), SE = logFC / t) |>
  transmute(genes = gene, logFC, t, P.Value = PValue, adj.P.Val = FDR,
            cell_type, cohort = "Xenium", SE)
coh <- bind_rows(coh |> select(genes, logFC, t, P.Value, adj.P.Val, cell_type,
                               cohort = cohort, SE),
                 xen_forest)

# ============================================================================
# Panel A — butterfly (ggplot)
# ============================================================================
build_butterfly <- function() {
  d <- meta_tbl |> mutate(direction = ifelse(estimate > 0, "up", "down"))
  cnt <- function(thr) d |> filter(padj < thr) |>
    count(cell_type, direction) |>
    pivot_wider(names_from = direction, values_from = n, values_fill = 0)
  c10 <- cnt(0.10); c05 <- cnt(0.05)
  for (cc in c("up","down")) { if (!cc %in% names(c10)) c10[[cc]] <- 0
                               if (!cc %in% names(c05)) c05[[cc]] <- 0 }
  m <- c10 |> transmute(cell_type, up10 = up, down10 = down) |>
    left_join(c05 |> transmute(cell_type, up05 = up, down05 = down), by = "cell_type") |>
    mutate(across(c(up05, down05), ~tidyr::replace_na(.x, 0)),
           total = up10 + down10) |>
    arrange(total)
  m$cell_type <- factor(m$cell_type, levels = m$cell_type)

  # fill mapped to a 4-level key so a legend is generated; light (FDR<0.10)
  # drawn first, dark (FDR<0.05) overlaid
  ggplot(m) +
    geom_col(aes(y = cell_type, x =  up10,   fill = "Up, FDR < 0.10"),   width = 0.78) +
    geom_col(aes(y = cell_type, x = -down10, fill = "Down, FDR < 0.10"), width = 0.78) +
    geom_col(aes(y = cell_type, x =  up05,   fill = "Up, FDR < 0.05"),   width = 0.78) +
    geom_col(aes(y = cell_type, x = -down05, fill = "Down, FDR < 0.05"), width = 0.78) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "black") +
    geom_text(aes(y = cell_type, x =  up10,   label = up10),
              hjust = -0.15, size = BASE*0.26) +
    geom_text(aes(y = cell_type, x = -down10, label = down10),
              hjust = 1.15, size = BASE*0.26) +
    scale_fill_manual(
      values = c("Up, FDR < 0.10" = UP_LIGHT, "Down, FDR < 0.10" = DOWN_LIGHT,
                 "Up, FDR < 0.05" = UP_DARK,  "Down, FDR < 0.05" = DOWN_DARK),
      breaks = c("Up, FDR < 0.05", "Up, FDR < 0.10",
                 "Down, FDR < 0.05", "Down, FDR < 0.10"), name = NULL) +
    scale_x_continuous(labels = function(v) abs(v),
                       expand = expansion(mult = 0.12)) +
    labs(x = "Number of DE genes", y = NULL) +
    theme_cowplot(font_size = BASE) +
    theme(axis.text.y = element_text(size = BASE - 1.5),
          axis.text.x = element_text(size = BASE - 1.5),
          axis.title.x = element_text(size = BASE - 1),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank(),
          legend.position = c(0.70, 0.16),
          legend.text  = element_text(size = BASE - 1.5),
          legend.key.size = unit(9, "pt"),
          legend.spacing.y = unit(0, "pt"),
          plot.margin  = margin(2, 4, 2, 2))
}

# ============================================================================
# Panels B-E — volcano for one cell type (Sst, L2/3 IT, Astro, Micro-PVM)
# ============================================================================
build_volcano <- function(cell, highlight) {
  # colour by the same FDR tiers as the butterfly (panel A): up/down x
  # FDR<0.05 (dark) / 0.05-0.10 (light); NS = grey
  d <- meta_tbl |> filter(cell_type == cell) |>
    mutate(nlp = -log10(padj),
           tier = case_when(
             padj < 0.05 & estimate > 0 ~ "Up, FDR < 0.05",
             padj < 0.10 & estimate > 0 ~ "Up, FDR < 0.10",
             padj < 0.05 & estimate < 0 ~ "Down, FDR < 0.05",
             padj < 0.10 & estimate < 0 ~ "Down, FDR < 0.10",
             TRUE                       ~ "NS") |>
             factor(levels = c("Down, FDR < 0.05","Down, FDR < 0.10","NS",
                               "Up, FDR < 0.10","Up, FDR < 0.05")))
  lab <- d |> filter(genes %in% highlight)
  # Tight, data-driven limits PER PANEL (not symmetric, not shared across the four
  # volcanoes); gene labels are kept inside the panel via ggrepel xlim/ylim +
  # coord_cartesian, so the axes crop close to the data without clipping a label.
  xr  <- range(d$estimate, na.rm = TRUE); xpad <- diff(xr) * 0.07
  xlo <- xr[1] - xpad; xhi <- xr[2] + xpad
  yhi <- max(d$nlp, na.rm = TRUE) * 1.08
  tier_cols <- c("Up, FDR < 0.05" = UP_DARK,  "Up, FDR < 0.10" = UP_LIGHT,
                 "Down, FDR < 0.05" = DOWN_DARK, "Down, FDR < 0.10" = DOWN_LIGHT,
                 "NS" = COL_NS)
  ggplot(d, aes(estimate, nlp)) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.25) +
    geom_hline(yintercept = -log10(0.1), colour = "grey70",
               linetype = "dashed", linewidth = 0.25) +
    geom_point(aes(colour = tier), data = ~filter(.x, tier == "NS"),
               size = 0.35, alpha = 0.3) +
    geom_point(aes(colour = tier), data = ~filter(.x, tier != "NS"),
               size = 0.5, alpha = 0.85) +
    geom_point(data = lab, shape = 21, fill = NA, colour = "black",
               size = 1.4, stroke = 0.4) +
    geom_text_repel(data = lab, aes(label = genes), size = BASE*0.32,
                    fontface = "italic", min.segment.length = 0,
                    segment.size = 0.25, segment.colour = "grey55",
                    box.padding = 0.4, point.padding = 0.3, force = 4,
                    max.overlaps = Inf, seed = 2, colour = "grey10",
                    xlim = c(xlo, xhi), ylim = c(0, yhi)) +
    scale_colour_manual(values = tier_cols, guide = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    coord_cartesian(xlim = c(xlo, xhi), ylim = c(0, yhi), clip = "on") +
    labs(x = expression("SCZ log"[2]~"FC"),
         y = expression(-log[10]~"FDR"),
         subtitle = gsub("_", "/", cell)) +
    theme_cowplot(font_size = BASE) +
    theme(plot.subtitle = element_text(size = BASE, face = "italic", hjust = 0.5,
                                       margin = margin(b = 1)),
          axis.text  = element_text(size = BASE - 1.5),
          axis.title = element_text(size = BASE - 1),
          plot.margin = margin(2, 3, 2, 2))
}

# ============================================================================
# Panels F-H — compact forest
# ============================================================================
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
    {if (has_xen) geom_hline(yintercept = 2, colour = "grey85",
                             linetype = "dotted", linewidth = 0.3)} +
    geom_segment(aes(x = lo, xend = hi, yend = y), linewidth = 0.45) +
    geom_point(data = cohort_df, aes(est, y), shape = 15, size = 1.2) +
    {if (!is.null(pool)) geom_point(data = pool, aes(est, y), shape = 23,
                                    size = 2.1, fill = COL_POOL, colour = COL_POOL)} +
    {if (has_xen) geom_point(data = xen_row, aes(est, y), shape = 17,
                             size = 1.8, colour = COL_XEN)} +
    geom_text(data = pd, aes(label = sig, x = ast_x, hjust = ast_h),
              size = BASE*0.34, vjust = 0.75, colour = "grey15", fontface = "bold") +
    geom_point(data = dplyr::filter(pd, dot), aes(ast_x, y), shape = 16,
               size = 0.8, colour = "grey15", inherit.aes = FALSE) +
    geom_text(data = dplyr::filter(pd, role=="pool" & sig==""),
              aes(ast_x, y, hjust = ast_h), label = "n.s.", size = BASE*0.28,
              vjust = 0.5, colour = "grey45", fontface = "italic", inherit.aes = FALSE) +
    scale_y_continuous(breaks = ylab$y, labels = ylab$cohort,
                       expand = expansion(add = 0.6)) +
    scale_x_continuous(limits = c(xlo, xhi), breaks = scales::pretty_breaks(3),
                       expand = expansion(mult = 0.01)) +
    labs(x = if (show_xlab) expression("SCZ log"[2]~"FC") else NULL,
         y = NULL, subtitle = ttl) +
    theme_cowplot(font_size = BASE) +
    theme(plot.subtitle = element_text(size = BASE, face = "plain",
                                       margin = margin(b = 1)),
          axis.text.y = element_text(size = BASE - 2),
          axis.text.x = element_text(size = BASE - 2),
          axis.title.x = element_text(size = BASE - 0.5, margin = margin(t = 1.5)),
          axis.line   = element_line(linewidth = 0.3),
          axis.ticks  = element_line(linewidth = 0.3),
          plot.margin = margin(2, 4, 2, 2))
}

# ============================================================================
# Panel I — concordance scatter
# ============================================================================
build_scatter <- function() {
  meta <- meta_tbl |> filter(padj < 0.10) |>
    transmute(genes, cell_type, meta_est = estimate, meta_padj = padj)
  xen <- xen_all |> transmute(genes = gene, cell_type, xen_logFC = logFC, xen_p = PValue)
  pr <- inner_join(meta, xen, by = c("genes","cell_type")) |>
    mutate(class = class_of(cell_type),
           concordant = sign(meta_est) == sign(xen_logFC),
           fdr_bin = factor(ifelse(meta_padj < 0.05, "< 0.05", "0.05-0.10"),
                            levels = c("< 0.05","0.05-0.10")))
  rr <- cor(pr$meta_est, pr$xen_logFC)
  pc <- round(100*mean(pr$concordant))
  lab_pairs <- tibble::tribble(~genes,~cell_type,
    "SST","Sst","BDNF","L2_3 IT","FKBP5","OPC","CX3CR1","Micro-PVM",
    "SMAD1","Pvalb","SERPING1","Astro","FGFR3","Astro")
  lab_keys <- paste(lab_pairs$genes, lab_pairs$cell_type)
  lab <- pr |> inner_join(lab_pairs, by = c("genes","cell_type"))
  lim <- max(as.numeric(quantile(abs(c(pr$meta_est, pr$xen_logFC)), 0.97)),
             max(abs(c(lab$meta_est, lab$xen_logFC)))) * 1.08
  # Labels via ggrepel auto-placement (robust to panel resizing); seed fixed for
  # determinism. Constrained to the panel; r / % concordant sit in the sparse
  # bottom-right corner.
  pr <- pr |> mutate(
    tag = sprintf("%s (%s)", genes, gsub("_", "/", cell_type)),
    rep_label = ifelse(paste(genes, cell_type) %in% lab_keys, tag, ""))

  ggplot(pr, aes(meta_est, xen_logFC)) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.25) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.25) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey80", linewidth = 0.3) +
    geom_smooth(method = "lm", se = TRUE, colour = "black", fill = "grey85",
                linewidth = 0.5, formula = y ~ x) +
    geom_point(aes(colour = class, size = fdr_bin), shape = 16, alpha = 0.8) +
    geom_point(data = lab, shape = 21, fill = NA, colour = "black",
               size = 1.8, stroke = 0.5) +
    geom_text_repel(data = pr, aes(label = rep_label), size = BASE*0.28,
                    min.segment.length = 0, segment.size = 0.25,
                    segment.colour = "grey55", box.padding = 0.45,
                    point.padding = 0.3, force = 2.5, force_pull = 0.15,
                    xlim = c(-lim, lim), ylim = c(-lim, lim),
                    max.overlaps = Inf, seed = 7, colour = "grey10") +
    annotate("text", x = lim*0.97, y = -lim*0.80,
             label = sprintf("italic(r)=='%.2f'", rr), parse = TRUE,
             hjust = 1, size = BASE*0.34) +
    annotate("text", x = lim*0.97, y = -lim*0.93,
             label = sprintf("'%d%% concordant'", pc), parse = TRUE,
             hjust = 1, size = BASE*0.30) +
    scale_colour_manual(values = CLASS_COL, name = NULL) +
    scale_size_manual(values = c("< 0.05" = 1.8, "0.05-0.10" = 0.7),
                      name = "meta FDR") +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(x = expression("snRNA-seq meta  log"[2]~"FC"),
         y = expression("Xenium  log"[2]~"FC")) +
    theme_cowplot(font_size = BASE) +
    # cell-class (colour) and meta-FDR (size) legends removed from the panel;
    # the encoding is described in the figure legend instead.
    theme(legend.position = "none",
          axis.text = element_text(size = BASE - 1.5),
          axis.title = element_text(size = BASE - 1),
          aspect.ratio = 1,            # keep the concordance scatter square
          plot.margin = margin(2, 3, 2, 2))
}

# ============================================================================
# Panel K — Xenium exemplar cells (boundary + marker molecules)
# Requires scripts/10_xenium_exemplar_cells.py to have been run.
# ============================================================================
TAB <- "results/tables"
DOT_COL <- "#B2182B"   # marker-molecule dots
# `lim` is SHARED across all four cells so the zoom is identical and the 5 um
# scale bar renders at the same physical length in every panel (consistent
# scale bar). Gene / condition / cell type are carried by the matrix labels in
# the assembly below, so no per-panel title here.
build_exemplar <- function(gene, dx, lim, scalebar_lab = FALSE) {
  bd <- read_csv(sprintf("%s/exemplar_%s_%s_boundary.csv", TAB, gene, dx),
                 show_col_types = FALSE)
  nu <- read_csv(sprintf("%s/exemplar_%s_%s_nucleus.csv", TAB, gene, dx),
                 show_col_types = FALSE)
  dt <- read_csv(sprintf("%s/exemplar_%s_%s_dots.csv", TAB, gene, dx),
                 show_col_types = FALSE)
  nd  <- nrow(dt)
  sb  <- 5  # scale-bar length, microns
  ggplot() +
    geom_polygon(data = bd, aes(x, y), fill = "grey93", colour = "grey45",
                 linewidth = 0.3) +
    geom_polygon(data = nu, aes(x, y), fill = NA, colour = "grey30",
                 linewidth = 0.25, linetype = "22") +     # nucleus = dashed
    {if (nd > 0) geom_point(data = dt, aes(x, y), colour = DOT_COL,
                            size = 0.5, alpha = 0.85)} +
    annotate("text", x = -lim*0.98, y = lim*0.98, label = nd,
             hjust = 0, vjust = 1, size = BASE*0.32, colour = DOT_COL,
             fontface = "bold") +
    annotate("segment", x = lim*0.96 - sb, xend = lim*0.96,
             y = -lim*0.94, yend = -lim*0.94, linewidth = 0.7) +
    {if (scalebar_lab) annotate("text", x = lim*0.96 - sb/2, y = -lim*0.80,
                                label = "5~mu*m", parse = TRUE,
                                size = BASE*0.26, vjust = 1)} +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    theme_void(base_size = BASE) +
    theme(plot.margin = margin(1, 2, 1, 2))
}

# ============================================================================
# Panel J — library-normalised expression (CP1K = counts per 1,000 transcripts)
# with the edgeR DE p-value. Same Xenium edgeR DE as the forests (F–H) and
# scatter (I): per-donor CP1K (TMM) from scripts/12_marker_norm_expr.R, p-value =
# edgeR glmQLFTest (exactly the value driving the forest/scatter significance).
# ============================================================================
NEXPR <- read_csv(sprintf("%s/marker_norm_expr.csv", TAB), show_col_types = FALSE)
NEXPR$dx <- factor(NEXPR$dx, levels = c("Control", "SCZ"))
NSTAT <- read_csv(sprintf("%s/marker_norm_expr_stats.csv", TAB), show_col_types = FALSE)
NEXPR_COL <- c(Control = DOWN_DARK, SCZ = UP_DARK)

build_normexpr <- function(gene, title, show_x = FALSE) {
  df <- NEXPR[NEXPR$gene == gene, ]
  p  <- NSTAT$p[NSTAT$gene == gene]
  yr <- range(df$cp1k)
  ggplot(df, aes(dx, cp1k)) +
    geom_boxplot(aes(fill = dx), width = 0.62, outlier.shape = NA, alpha = 0.55, linewidth = 0.4) +
    geom_jitter(aes(fill = dx), shape = 21, colour = "grey25", size = 1.4, stroke = 0.3,
                width = 0.13, height = 0, alpha = 0.9) +
    geom_signif(comparisons = list(c("Control", "SCZ")),
                annotations = sprintf("italic(p)=='%.3f'", p), parse = TRUE,
                y_position = yr[2] + 0.06 * diff(yr), tip_length = 0.02,
                textsize = BASE * 0.32, vjust = -0.1) +
    scale_fill_manual(values = NEXPR_COL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22)), n.breaks = 4) +
    labs(x = NULL, y = "counts / 1,000 tx", subtitle = title) +
    theme_cowplot(font_size = BASE) +
    theme(legend.position = "none",
          plot.subtitle = element_text(size = BASE - 1, face = "plain", hjust = 0.5,
                                       margin = margin(b = 1)),
          axis.text.x  = if (show_x) element_text(size = BASE - 1.5) else element_blank(),
          axis.ticks.x = if (show_x) element_line() else element_blank(),
          axis.line.x  = if (show_x) element_line() else element_blank(),
          axis.title.y = element_text(size = BASE - 1.5),
          axis.text.y  = element_text(size = BASE - 1.5),
          plot.margin = margin(2, 3, 2, 3))
}

# ============================================================================
# Assemble
# ============================================================================
cat("Building panels...\n")
pA <- build_butterfly()
pB <- build_volcano("Sst",       c("SST","NAT16","SMAD1","AFG3L2"))
pC <- build_volcano("L2_3 IT",   c("BDNF","SMAD1","VWA5B2","ADAMTS9-AS2","ST6GAL2"))
pD <- build_volcano("Astro",     c("SERPING1","CHI3L1","FGFR3","NOTCH1"))
pE <- build_volcano("Micro-PVM", c("C1QA","C1QB","CX3CR1","P2RY12","SORL1"))
# 3 forests in a single row -> all carry the "SCZ log2 FC" x-axis label
forests <- Map(build_forest, FOREST$gene, FOREST$cell, FOREST$lab,
               show_xlab = c(TRUE, TRUE, TRUE, TRUE))
pJ <- build_scatter()
# shared coordinate limit -> consistent zoom + identical 5 um scale bar across all 4 cells
ex_lim <- 1.15 * max(vapply(
  list(c("SST","Control"), c("SST","SCZ"), c("FGFR3","Control"), c("FGFR3","SCZ")),
  function(p) { bd <- read_csv(sprintf("%s/exemplar_%s_%s_boundary.csv", TAB, p[1], p[2]),
                               show_col_types = FALSE); max(abs(c(bd$x, bd$y))) },
  numeric(1)))
ex <- list(build_exemplar("SST","Control", ex_lim), build_exemplar("SST","SCZ", ex_lim),
           build_exemplar("FGFR3","Control", ex_lim),
           build_exemplar("FGFR3","SCZ", ex_lim, scalebar_lab = TRUE))

# Row 1: butterfly (a, ~40% width) | 2x2 volcano grid (~60% width)
#   b Sst       | c L2/3 IT
#   d Astro     | e Micro-PVM
volc_grid <- plot_grid(pB, pC, pD, pE, ncol = 2, labels = c("b","c","d","e"),
                       label_size = 8, label_fontface = "bold")
row1 <- plot_grid(pA, volc_grid, ncol = 2, rel_widths = c(2, 3),
                  labels = c("a",""), label_size = 8, label_fontface = "bold")

# Row 2: 4 forest panels in one row (f = SST/Sst, g = BDNF/L2/3 IT,
#        h = FGFR3/Astro, i = FKBP5/OPC)
row2 <- plot_grid(plotlist = forests, ncol = 4,
                  labels = c("f","g","h","i"),
                  label_size = 8, label_fontface = "bold")

# Row 3: scatter (i) | CP1K boxplots (j, titled per row) | exemplar matrix (k)
# Exemplar matrix (k): columns = condition (Control / SCZ), rows = the marker gene
# shown in its cell type (SST in Sst, FGFR3 in astrocytes). Panel j now carries the
# per-row title (replacing the old rotated shared label), aligned with k's rows.
# Each k cell: solid grey = cell boundary, dashed = nucleus, red dots = marker mRNA.
ex_hdr  <- function(t) ggdraw() + draw_label(t, size = BASE - 0.5)
RH <- c(0.14, 1, 1)   # header | SST row | FGFR3 row (shared across J / K)

# J — library-normalised expression boxplots, each TITLED (the title also names
# the aligned exemplar row in k, so no separate left label is needed)
k_col <- plot_grid(NULL,
                   build_normexpr("SST",   "SST mRNA in Sst cells"),
                   build_normexpr("FGFR3", "FGFR3 mRNA in Astrocytes", show_x = TRUE),
                   ncol = 1, rel_heights = RH)
# K — Xenium exemplar cells (Control / SCZ headers + two cell rows)
l_col <- plot_grid(plot_grid(ex_hdr("Control"), ex_hdr("SCZ"), ncol = 2),
                   plot_grid(ex[[1]], ex[[2]], ncol = 2),
                   plot_grid(ex[[3]], ex[[4]], ncol = 2),
                   ncol = 1, rel_heights = RH)
row3 <- plot_grid(pJ, k_col, l_col, ncol = 3,
                  rel_widths = c(1.0, 0.62, 0.85),
                  labels = c("j", "k", "l"), label_size = 8, label_fontface = "bold")

full <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(1.0, 0.55, 1.05))

FIG_H <- 6.7   # max height (in); width 7.1. Forests in one row (f–i); scatter (j)
               # column sized so the square scatter fills it (minimal dead space)
ggsave("results/09_composite.png", full, width = FIG_W, height = FIG_H,
       dpi = 400, bg = "white")
ggsave("results/09_composite.pdf", full, width = FIG_W, height = FIG_H, bg = "white")
cat(sprintf("Saved results/09_composite.{png,pdf}  (%.1f x %.1f in)\n", FIG_W, FIG_H))
