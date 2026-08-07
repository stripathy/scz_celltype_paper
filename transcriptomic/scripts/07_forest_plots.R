#!/usr/bin/env Rscript
# ============================================================================
# 07_forest_plots.R — Cross-platform forest plots per (gene, cell type)
# ============================================================================
# For one gene in one SEA-AD subclass, shows the SCZ-vs-control effect in each
# snRNA-seq discovery cohort, the random-effects pooled meta-analytic estimate,
# and the independent Xenium spatial replication (held OUT of the pool).
#
# Layout (top -> bottom): individual snRNA-seq cohorts, the pooled "snRNA-seq
# meta" diamond, a dotted separator, then the Xenium replication triangle.
#
# Significance markers (per row):
#   ***  FDR < 0.01     **  FDR < 0.05     *  FDR < 0.10
#   dot  uncorrected p < 0.05 in an individual dataset (cohort or Xenium)
#        that does NOT reach FDR < 0.10
#   n.s. shown on the meta row only, when the pooled estimate is not FDR-sig
#
# Usage:
#   Rscript scripts/07_forest_plots.R                   # 6-panel composite
#   Rscript scripts/07_forest_plots.R GENE "CELL TYPE"  # single pair (e.g. SMAD1 Sst)
#
# Outputs:
#   results/07_forest_composite.{png,pdf}          (composite mode)
#   results/07_forest_<gene>_<celltype>.{png,pdf}  (single-pair mode)
#
# DATA PROVENANCE — update these paths/columns if the DE results change:
#   INPUT_COH    Per-cohort snRNA-seq DE (limma). One row per cohort x cell x gene.
#                Used cols: genes, logFC, t, P.Value, adj.P.Val, cell_type, cohort.
#                7 cohorts: Bat, HBCC, MSSM, Mclean, MtSinai, Multi, OFC.
#   INPUT_META   Meta-analytic DE (metafor). Used ONLY for the canonical meta FDR
#                that sets the diamond's asterisks. Used cols: genes, cell_type, padj.
#   INPUT_XENIUM Xenium spatial pseudobulk DE (edgeR QL F-test), Kwon 2026,
#                12 SCZ vs 12 control DLPFC, 300-gene panel. Used cols:
#                gene, logFC, F, PValue, FDR, celltype.
#
# Per-row SE: logFC / t (limma cohorts); logFC / (sign(logFC)*sqrt(F)) for
# Xenium edgeR. The meta diamond is recomputed here via metafor rma(DL) from
# the 7 cohort logFC +/- SE; its asterisks use INPUT_META's padj (not the
# recomputed p) so they match the project's canonical meta-analysis.
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(cowplot); library(metafor)
})

source("scripts/_figure_inputs.R")     # committed snapshots + staleness guard
INPUT_COH    <- fig_input("meta_results_cohorts_subclass_forest.csv")
INPUT_META   <- fig_input("DE_genes_all_cells_scz.csv")
INPUT_XENIUM <- fig_input("de_results_subclass.csv")
FIG_DIR      <- "results"

COL_POOL   <- "black"      # meta-analysis diamond
COL_XEN    <- "#117733"    # Xenium replication triangle (green = independent)
COL_COHORT <- "grey35"     # individual snRNA-seq cohort squares

# ---- Loaders (read once at script level) ----
cat("Loading cohort-level table (large file, ~310 MB)...\n")
coh <- read_csv(INPUT_COH, show_col_types = FALSE) |>
  rename(row_id = 1) |>
  filter(!is.na(logFC), !is.na(t), t != 0) |>
  mutate(SE = logFC / t)
cat(sprintf("  %d rows, %d unique genes\n", nrow(coh), n_distinct(coh$genes)))

cat("Loading meta-analytic table...\n")
meta_tbl <- read_csv(INPUT_META, show_col_types = FALSE)

# ---- Load + harmonize Xenium (8th cohort, 300-gene panel) ----
if (file.exists(INPUT_XENIUM)) {
  cat("Loading Xenium DE (Kwon 2026 spatial, 300-gene panel)...\n")
  xen_raw <- read_csv(INPUT_XENIUM, show_col_types = FALSE)
  # Subclass-name harmonization
  ct_map <- c("Astrocyte"="Astro", "L2/3 IT"="L2_3 IT", "L5/6 NP"="L5_6 NP",
              "Microglia-PVM"="Micro-PVM", "Oligodendrocyte"="Oligo",
              "Endothelial"="Endo")
  xen <- xen_raw |>
    mutate(cell_type = ifelse(celltype %in% names(ct_map),
                              ct_map[celltype], celltype),
           # For edgeR QL F-test with 1-df numerator, t = sign(logFC) * sqrt(F)
           t  = sign(logFC) * sqrt(F),
           SE = ifelse(t != 0, logFC / t, NA_real_)) |>
    filter(!is.na(SE), is.finite(SE)) |>
    transmute(row_id = paste0(gene, "...xen"),
              genes   = gene,
              logFC, AveExpr = logCPM, t,
              P.Value = PValue, adj.P.Val = FDR, B = F,
              cell_type, cohort = "Xenium", SE)
  cat(sprintf("  Xenium: %d rows, %d unique genes (subset of 300-gene panel)\n",
              nrow(xen), n_distinct(xen$genes)))
  coh <- bind_rows(coh, xen)
}

# ---- Forest-plot builder for ONE (gene, cell_type) pair ------------------
make_forest <- function(gene_sym, cell, label = NULL) {
  sub <- coh |> filter(genes == gene_sym, cell_type == cell)
  if (nrow(sub) == 0) {
    cat(sprintf("  [%s in %s] no cohort data\n", gene_sym, cell))
    return(NULL)
  }

  # ---- SPLIT: snRNA-seq cohorts (for meta) vs Xenium (replication only) -----
  snrna  <- sub |> filter(cohort != "Xenium")
  xen_in <- sub |> filter(cohort == "Xenium")

  # Random-effects meta on snRNA-seq ONLY (Xenium is held out as replication)
  re <- tryCatch(
    rma(yi = snrna$logFC, sei = snrna$SE, method = "DL",
        control = list(maxiter = 500)),
    error = function(e) NULL
  )

  # ---- Build rows (Cochrane convention: cohorts on top, diamond below) -----
  # Significance markers:
  #   ***  FDR < 0.01      **  FDR < 0.05      *  FDR < 0.10
  #   dot  uncorrected p < 0.05 in an individual dataset (cohort / Xenium)
  #        that does NOT reach FDR < 0.10
  ast <- function(fdr) dplyr::case_when(
    is.na(fdr)  ~ "",
    fdr < 0.01  ~ "***",
    fdr < 0.05  ~ "**",
    fdr < 0.10  ~ "*",
    TRUE        ~ ""
  )
  is_dot <- function(p, fdr) !is.na(p) & p < 0.05 & (is.na(fdr) | fdr >= 0.10)

  # snRNA-seq cohort rows — asterisks from per-cohort FDR, dot from nominal p
  cohort_df <- snrna |>
    transmute(cohort, est = logFC, lo = logFC - 1.96 * SE, hi = logFC + 1.96 * SE,
              p = P.Value, padj = adj.P.Val,
              sig = ast(adj.P.Val), dot = is_dot(P.Value, adj.P.Val)) |>
    arrange(desc(est)) |>
    mutate(y = rev(seq_along(est)) + 3,
           role = "snRNA-seq cohort",
           colour = COL_COHORT)

  # Pooled diamond — asterisks from the canonical meta-analysis FDR (meta CSV)
  meta_row <- meta_tbl |> filter(genes == gene_sym, cell_type == cell)
  meta_fdr <- if (nrow(meta_row) == 1) meta_row$padj else NA_real_
  pool <- if (!is.null(re)) {
    tibble(
      cohort = "snRNA-seq meta",
      est = as.numeric(re$beta), lo = re$ci.lb, hi = re$ci.ub,
      p = re$pval, padj = meta_fdr,
      sig = ast(meta_fdr), dot = FALSE,   # meta is a synthesis, never gets a dot
      y = 3, role = "pool", colour = COL_POOL)
  } else NULL

  # Xenium replication — asterisks from Xenium FDR, dot from Xenium nominal p
  has_xen <- nrow(xen_in) > 0
  xen_row <- if (has_xen) {
    xr <- xen_in[1, ]
    tibble(
      cohort = "Xenium",
      est = xr$logFC, lo = xr$logFC - 1.96*xr$SE, hi = xr$logFC + 1.96*xr$SE,
      p = xr$P.Value, padj = xr$adj.P.Val,
      sig = ast(xr$adj.P.Val), dot = is_dot(xr$P.Value, xr$adj.P.Val),
      y = 1, role = "xenium", colour = COL_XEN)
  } else NULL

  # Layout — data-driven, asymmetric x-limits (not forced symmetric around 0,
  # which wastes space). Pad only the side(s) that carry markers, enough to
  # fit the asterisks / dot / n.s. text.
  plot_data <- bind_rows(cohort_df, pool, xen_row)
  data_lo <- min(plot_data$lo, na.rm = TRUE)
  data_hi <- max(plot_data$hi, na.rm = TRUE)
  rng     <- data_hi - data_lo
  off     <- rng * 0.03                       # marker sits just outside the CI
  ast_data <- plot_data |>
    mutate(ast_x = ifelse(est >= 0, hi + off, lo - off),
           ast_h = ifelse(est >= 0, 0, 1),
           has_marker = sig != "" | dot | (role == "pool" & sig == ""))
  right_marker <- any(ast_data$has_marker & ast_data$ast_h == 0)
  left_marker  <- any(ast_data$has_marker & ast_data$ast_h == 1)
  xlo <- data_lo - rng * (if (left_marker)  0.20 else 0.04)
  xhi <- data_hi + rng * (if (right_marker) 0.20 else 0.04)
  xlo <- min(xlo, -off); xhi <- max(xhi, off)  # always keep the x = 0 reference in view

  # y-axis labels (cohort names) — gap row at y=2 for visual separation
  y_labels_df <- plot_data |> select(y, cohort) |> arrange(y)
  if (has_xen) y_labels_df <- bind_rows(y_labels_df, tibble(y = 2, cohort = ""))

  # Identifier subtitle (clean, minimal): gene + cell, italicised gene name
  panel_letter <- if (!is.null(label)) sub("\\..*", "", label) else NA_character_
  id_subtitle  <- bquote(italic(.(gene_sym)) ~ " in " ~ .(cell))

  p <- ggplot(plot_data, aes(x = est, y = y, colour = I(colour))) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey60", linewidth = 0.5) +
    {if (has_xen)
       list(geom_hline(yintercept = 2, colour = "grey80",
                       linewidth = 0.6, linetype = "dotted"))
     } +
    geom_segment(aes(x = lo, xend = hi, yend = y), linewidth = 1.1) +
    geom_point(data = cohort_df, aes(x = est, y = y),
               shape = 15, size = 5) +
    {if (!is.null(pool))
       list(geom_point(data = pool, aes(x = est, y = y),
                       shape = 23, size = 9,
                       fill = COL_POOL, colour = COL_POOL))
     } +
    {if (has_xen)
       list(geom_point(data = xen_row, aes(x = est, y = y),
                       shape = 17, size = 7, colour = COL_XEN))
     } +
    # Significance asterisks (from FDR)
    geom_text(data = ast_data, aes(label = sig, x = ast_x, hjust = ast_h),
              size = 11, vjust = 0.72, colour = "grey15", fontface = "bold",
              show.legend = FALSE) +
    # Dot = nominal p<0.05 (not FDR<0.1) in an individual dataset
    geom_point(data = dplyr::filter(ast_data, dot),
               aes(x = ast_x, y = y), shape = 16, size = 3.4,
               colour = "grey15", inherit.aes = FALSE, show.legend = FALSE) +
    # "n.s." on the meta row when the pooled estimate is not FDR-significant
    geom_text(data = dplyr::filter(ast_data, role == "pool" & sig == ""),
              aes(x = ast_x, y = y, hjust = ast_h), label = "n.s.",
              size = 6.5, vjust = 0.5, colour = "grey40", fontface = "italic",
              inherit.aes = FALSE, show.legend = FALSE) +
    scale_y_continuous(breaks = y_labels_df$y, labels = y_labels_df$cohort,
                       expand = expansion(add = 0.7)) +
    scale_x_continuous(limits = c(xlo, xhi),
                       breaks = scales::pretty_breaks(n = 4),
                       expand = expansion(mult = 0.01)) +
    labs(x = expression("SCZ log"[2]~"FC"), y = NULL,
         subtitle = id_subtitle) +
    theme_cowplot(font_size = 24) +
    theme(plot.subtitle    = element_text(size = 28, face = "plain",
                                          margin = margin(b = 10)),
          axis.text.y      = element_text(size = 22),
          axis.text.x      = element_text(size = 20),
          axis.title.x     = element_text(size = 22),
          axis.line        = element_line(linewidth = 0.8),
          axis.ticks       = element_line(linewidth = 0.8),
          plot.margin      = margin(12, 20, 10, 12))

  # Bold panel letter in the top-left margin (outside the axes), via plot tag
  if (!is.na(panel_letter)) {
    p <- p +
      labs(tag = panel_letter) +
      theme(plot.tag = element_text(face = "bold", size = 30,
                                    hjust = 0, vjust = 1),
            plot.tag.position = c(0.005, 0.995))
  }
  p
}

# ---- Mode dispatch -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2) {
  # Single-pair mode
  gene_sym <- args[[1]]; cell <- args[[2]]
  p <- make_forest(gene_sym, cell)
  if (is.null(p)) quit(status = 1)
  slug <- gsub("[^A-Za-z0-9]+", "_", paste(gene_sym, cell, sep = "_"))
  ggsave(file.path(FIG_DIR, sprintf("07_forest_%s.png", slug)),
         p, width = 9, height = 4.5, dpi = 220, bg = "white")
  ggsave(file.path(FIG_DIR, sprintf("07_forest_%s.pdf", slug)),
         p, width = 9, height = 4.5, bg = "white")
  cat(sprintf("Saved results/07_forest_%s.{png,pdf}\n", slug))
} else {
  # Multi-panel composite (publication-quality). Detail goes to figure legend.
  panels <- tibble::tribble(
    ~gene,      ~cell,        ~label,
    "SST",      "Sst",        "A.",
    "PVALB",    "Pvalb",      "B.",
    "BDNF",     "L2_3 IT",    "C.",
    "FKBP5",    "OPC",        "D.",
    "CX3CR1",   "Micro-PVM",  "E.",
    "SMAD1",    "Pvalb",      "F."
  )

  cat("\n=== Building 6-panel composite forest plot ===\n")
  plots <- vector("list", nrow(panels))
  for (i in seq_len(nrow(panels))) {
    cat(sprintf("[%d/%d] %s in %s ... ", i, nrow(panels),
                panels$gene[i], panels$cell[i]))
    p <- make_forest(panels$gene[i], panels$cell[i], label = panels$label[i])
    if (!is.null(p)) cat("ok\n") else cat("(no data)\n")
    plots[[i]] <- p
  }
  plots <- Filter(Negate(is.null), plots)

  grid <- plot_grid(plotlist = plots, ncol = 2, align = "v")
  ggsave(file.path(FIG_DIR, "07_forest_composite.png"),
         grid, width = 22, height = 20, dpi = 220, bg = "white",
         limitsize = FALSE)
  ggsave(file.path(FIG_DIR, "07_forest_composite.pdf"),
         grid, width = 22, height = 20, bg = "white",
         limitsize = FALSE)
  cat("\nSaved results/07_forest_composite.{png,pdf}\n")
}
