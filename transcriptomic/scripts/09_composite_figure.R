#!/usr/bin/env Rscript
# ============================================================================
# 09_composite_figure.R — Publication composite (total width = 6.5 inches)
# ============================================================================
# Assembles, in one figure, the cross-platform SCZ DE story:
#   A  Butterfly: up/down DE counts per cell type, FDR<0.10 (light) vs
#      FDR<0.05 (dark overlay)
#   B  Volcano, Sst                (select genes labelled)
#   C  Volcano, L2/3 IT            (select genes labelled)
#   D-I Forest plots for 6 (gene, cell) pairs: 7 snRNA-seq cohorts + pooled
#      meta + Xenium replication
#   J  Concordance scatter: snRNA-seq meta logFC vs Xenium logFC
#   K  Xenium exemplar cells: boundary + marker molecules for one Control and
#      one SCZ cell each for SST (Sst) and RASGRF2 (Pvalb) -- both strongly
#      meta-DE-down genes with enough per-cell expression to show as dots
#
# DATA PROVENANCE (see notes/figures_crossplatform_validation.md for detail):
#   data/DE_genes_all_cells_scz.csv          meta-analytic snRNA-seq DE
#   data/meta_results_cohorts_subclass.csv   per-cohort snRNA-seq DE (7 cohorts)
#   ~/Github/SCZ_Xenium/output/de/de_results_subclass.csv   Xenium spatial DE
#   results/tables/exemplar_*.csv            panel K inputs — produced by
#      scripts/10_xenium_exemplar_cells.py (run that FIRST; it reads the
#      Xenium h5ad + boundary + transcript exports).
#
# Output: results/09_composite.{png,pdf}   (6.5 in wide)
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(cowplot); library(ggrepel); library(metafor)
})

INPUT_META   <- "data/DE_genes_all_cells_scz.csv"
INPUT_COH    <- "data/meta_results_cohorts_subclass.csv"
INPUT_XENIUM <- "../spatial/output/de/de_results_subclass.csv"  # internal cross-ref (was ~/Github/SCZ_Xenium/...)

FIG_W <- 6.5    # required total width (inches)

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

BASE <- 9   # base font size for the whole composite (most text sized relative to this)

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
  "PVALB",  "Pvalb",      "PVALB / Pvalb",
  "BDNF",   "L2_3 IT",    "BDNF / L2/3 IT",
  "FKBP5",  "OPC",        "FKBP5 / OPC",
  "CX3CR1", "Micro-PVM",  "CX3CR1 / Micro-PVM",
  "RASGRF2","Pvalb",      "RASGRF2 / Pvalb"
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
              hjust = -0.15, size = BASE*0.20) +
    geom_text(aes(y = cell_type, x = -down10, label = down10),
              hjust = 1.15, size = BASE*0.20) +
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
          legend.text  = element_text(size = BASE - 2.5),
          legend.key.size = unit(9, "pt"),
          legend.spacing.y = unit(0, "pt"),
          plot.margin  = margin(2, 4, 2, 2))
}

# ============================================================================
# Panels B/C — volcano for one cell type
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
  xm  <- max(abs(d$estimate), na.rm = TRUE)
  lab <- d |> filter(genes %in% highlight)
  ymax <- max(d$nlp, na.rm = TRUE)
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
                    box.padding = 0.5, point.padding = 0.3, force = 5,
                    max.overlaps = Inf, seed = 2, colour = "grey10") +
    scale_colour_manual(values = tier_cols, guide = "none") +
    # roomy limits so edge/top labels are not clipped
    scale_x_continuous(limits = c(-xm*1.30, xm*1.30)) +
    scale_y_continuous(limits = c(0, ymax*1.18)) +
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
# Panels D-I — compact forest
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
              size = BASE*0.42, vjust = 0.75, colour = "grey15", fontface = "bold") +
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
    theme(plot.subtitle = element_text(size = BASE - 0.5, face = "bold",
                                       margin = margin(b = 1)),
          axis.text.y = element_text(size = BASE - 2),
          axis.text.x = element_text(size = BASE - 2),
          axis.title.x = element_text(size = BASE - 0.5, margin = margin(t = 1.5)),
          axis.line   = element_line(linewidth = 0.3),
          axis.ticks  = element_line(linewidth = 0.3),
          plot.margin = margin(2, 4, 2, 2))
}

# ============================================================================
# Panel J — concordance scatter
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
    "SMAD1","Pvalb","SERPING1","Astro","RASGRF2","Pvalb")
  lab_keys <- paste(lab_pairs$genes, lab_pairs$cell_type)
  lab <- pr |> inner_join(lab_pairs, by = c("genes","cell_type"))
  lim <- max(as.numeric(quantile(abs(c(pr$meta_est, pr$xen_logFC)), 0.97)),
             max(abs(c(lab$meta_est, lab$xen_logFC)))) * 1.08
  # Labels are placed manually (low repel force + explicit per-gene nudges) so
  # the crowded scatter is deterministic: nudges below are (desired_label_pos -
  # point) for each gene. SERPING1/FKBP5 (top-right points) label into the freed
  # top-left; CX3CR1/BDNF are stacked vertically in the sparse lower area.
  pr <- pr |> mutate(
    tag = sprintf("%s (%s)", genes, gsub("_", "/", cell_type)),
    rep_label = ifelse(paste(genes, cell_type) %in% lab_keys, tag, ""),
    nx = dplyr::case_when(rep_label == "" ~ 0,
                          genes == "SERPING1" ~ -0.56, genes == "FKBP5" ~ -0.26,
                          genes == "SMAD1" ~  0.27, genes == "SST" ~ -0.16,
                          genes == "CX3CR1" ~  0.24, genes == "BDNF" ~ -0.09,
                          genes == "RASGRF2" ~ -0.10, TRUE ~ 0),
    ny = dplyr::case_when(rep_label == "" ~ 0,
                          genes == "SERPING1" ~  0.01, genes == "FKBP5" ~  0.01,
                          genes == "SMAD1" ~  0.20, genes == "SST" ~  0.33,
                          genes == "CX3CR1" ~ -0.21, genes == "BDNF" ~ -0.20,
                          genes == "RASGRF2" ~  0.36, TRUE ~ 0))

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
    geom_text_repel(data = pr, aes(label = rep_label), size = BASE*0.30,
                    nudge_x = pr$nx, nudge_y = pr$ny,
                    min.segment.length = 0, segment.size = 0.25,
                    segment.colour = "grey55", box.padding = 0.5,
                    point.padding = 0.3, force = 1, force_pull = 0.1,
                    xlim = c(-lim, lim), ylim = c(-lim, lim),
                    max.overlaps = Inf, seed = 7, colour = "grey10") +
    annotate("text", x = lim*0.98, y = -lim*0.56,
             label = sprintf("italic(r)=='%.2f'", rr), parse = TRUE,
             hjust = 1, size = BASE*0.42) +
    annotate("text", x = lim*0.98, y = -lim*0.71,
             label = sprintf("'%d%% concordant'", pc), parse = TRUE,
             hjust = 1, size = BASE*0.34) +
    scale_colour_manual(values = CLASS_COL, name = NULL) +
    scale_size_manual(values = c("< 0.05" = 1.8, "0.05-0.10" = 0.7),
                      name = "meta FDR") +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(x = expression("snRNA-seq meta  log"[2]~"FC"),
         y = expression("Xenium  log"[2]~"FC")) +
    guides(colour = guide_legend(override.aes = list(size = 1.8), order = 1, nrow = 1),
           size = guide_legend(order = 2, nrow = 1)) +
    theme_cowplot(font_size = BASE) +
    # legend moved OUT of the plot (to the bottom): the top-left corner is the
    # only place the wide top-right labels (SERPING1, FKBP5) can extend into.
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.margin = margin(t = 0, b = 0), legend.box.spacing = unit(2, "pt"),
          legend.box.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(1, "pt"),
          legend.text = element_text(size = BASE - 2),
          legend.title = element_text(size = BASE - 1.5),
          legend.key.height = unit(8, "pt"), legend.key.width = unit(8, "pt"),
          legend.spacing.x = unit(3, "pt"),
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
  dt <- read_csv(sprintf("%s/exemplar_%s_%s_dots.csv", TAB, gene, dx),
                 show_col_types = FALSE)
  nd  <- nrow(dt)
  sb  <- 5  # scale-bar length, microns
  ggplot() +
    geom_polygon(data = bd, aes(x, y), fill = "grey93", colour = "grey45",
                 linewidth = 0.3) +
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
# Assemble
# ============================================================================
cat("Building panels...\n")
pA <- build_butterfly()
pB <- build_volcano("Sst",    c("SST","NAT16","SMAD1","AFG3L2"))
pC <- build_volcano("L2_3 IT", c("BDNF","SMAD1","VWA5B2","ADAMTS9-AS2","ST6GAL2"))
# x-axis label ("SCZ log2 FC") only on the bottom row of the 3x2 forest grid (G,H,I)
forests <- Map(build_forest, FOREST$gene, FOREST$cell, FOREST$lab,
               show_xlab = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
pJ <- build_scatter()
# shared coordinate limit -> consistent zoom + identical 5 um scale bar across all 4 cells
ex_lim <- 1.15 * max(vapply(
  list(c("SST","Control"), c("SST","SCZ"), c("RASGRF2","Control"), c("RASGRF2","SCZ")),
  function(p) { bd <- read_csv(sprintf("%s/exemplar_%s_%s_boundary.csv", TAB, p[1], p[2]),
                               show_col_types = FALSE); max(abs(c(bd$x, bd$y))) },
  numeric(1)))
ex <- list(build_exemplar("SST","Control", ex_lim), build_exemplar("SST","SCZ", ex_lim),
           build_exemplar("RASGRF2","Control", ex_lim),
           build_exemplar("RASGRF2","SCZ", ex_lim, scalebar_lab = TRUE))

# Row 1: butterfly (~65% width) | (volcano Sst over volcano L2/3)
volc_col <- plot_grid(pB, pC, ncol = 1, labels = c("B","C"),
                      label_size = 12, label_fontface = "bold")
row1 <- plot_grid(pA, volc_col, ncol = 2, rel_widths = c(1.85, 1),
                  labels = c("A",""), label_size = 12, label_fontface = "bold")

# Row 2: 6 forest panels (3 cols x 2 rows)
row2 <- plot_grid(plotlist = forests, ncol = 3,
                  labels = c("D","E","F","G","H","I"),
                  label_size = 12, label_fontface = "bold")

# Row 3: scatter (J) | exemplar-cell matrix (K)
# K as a labelled matrix: columns = condition (Control / SCZ), rows = the
# marker gene shown in its cell type (SST in Sst, RASGRF2 in Pvalb).
ex_hdr  <- function(t) ggdraw() + draw_label(t, size = BASE - 0.5)
ex_rlab <- function(t) ggdraw() + draw_label(t, size = BASE - 1, angle = 90,
                                             fontface = "italic")
RW <- c(0.13, 1, 1)   # row-label | Control | SCZ  (shared so columns align)
ex_cols <- plot_grid(NULL, ex_hdr("Control"), ex_hdr("SCZ"), ncol = 3, rel_widths = RW)
ex_r1 <- plot_grid(ex_rlab("SST in Sst"),       ex[[1]], ex[[2]], ncol = 3, rel_widths = RW)
ex_r2 <- plot_grid(ex_rlab("RASGRF2 in Pvalb"), ex[[3]], ex[[4]], ncol = 3, rel_widths = RW)
ex_block <- plot_grid(ex_cols, ex_r1, ex_r2, ncol = 1, rel_heights = c(0.14, 1, 1))
row3 <- plot_grid(pJ, ex_block, ncol = 2, rel_widths = c(1, 1),
                  labels = c("J","K"), label_size = 12, label_fontface = "bold")

full <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(1.15, 1.05, 1.0))

FIG_H <- 9.2   # width is locked at 6.5"; grow height to give the larger text room
ggsave("results/09_composite.png", full, width = FIG_W, height = FIG_H,
       dpi = 400, bg = "white")
ggsave("results/09_composite.pdf", full, width = FIG_W, height = FIG_H, bg = "white")
cat(sprintf("Saved results/09_composite.{png,pdf}  (%.1f x %.1f in)\n", FIG_W, FIG_H))
