#!/usr/bin/env Rscript
# Volcano plot of SCZ DE genes for a single cell type with ggrepel labels.
#
# Usage:
#   Rscript scripts/01_volcano_per_celltype.R                       # default: Sst
#   Rscript scripts/01_volcano_per_celltype.R Pvalb                 # any cell type
#   Rscript scripts/01_volcano_per_celltype.R Pvalb KCNC1,KCNAB3,AFG3L2  # also force-label these
#
# Following Shreejoy's CLAUDE.md: no in-plot title, data-driven limits,
# Okabe-Ito palette, large legible text, threshold line labels.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(ggrepel); library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
CELL_TYPE     <- if (length(args) >= 1) args[[1]] else "Sst"
FORCE_LABEL   <- if (length(args) >= 2) strsplit(args[[2]], ",")[[1]] else character(0)

source("scripts/_figure_inputs.R")     # committed snapshots + staleness guard
INPUT         <- fig_input("DE_genes_all_cells_scz.csv")
PADJ_SIG      <- 0.1
PADJ_STRICT   <- 0.05
slug          <- gsub("[^A-Za-z0-9]+", "_", CELL_TYPE)
OUT_FIG       <- sprintf("results/01_volcano_%s.png", slug)
OUT_PDF       <- sprintf("results/01_volcano_%s.pdf", slug)

# Okabe-Ito
COL_UP    <- "#D55E00"
COL_DOWN  <- "#0072B2"
COL_NS    <- "grey75"
COL_FORCE <- "#117733"   # green for force-labeled (curated) genes

# ---- Data ------------------------------------------------------------------
df <- read_csv(INPUT, show_col_types = FALSE) |>
  filter(cell_type == CELL_TYPE) |>
  mutate(
    neglog10_padj = -log10(padj),
    direction = case_when(
      padj < PADJ_SIG & estimate > 0 ~ "Up",
      padj < PADJ_SIG & estimate < 0 ~ "Down",
      TRUE                           ~ "NS"
    ) |> factor(levels = c("Down", "NS", "Up"))
  )

n_up   <- sum(df$direction == "Up")
n_down <- sum(df$direction == "Down")
cat(sprintf("%s: %d up, %d down at padj < %g (of %d tested)\n",
            CELL_TYPE, n_up, n_down, PADJ_SIG, nrow(df)))

# ---- Pick label set --------------------------------------------------------
# Top by padj in each direction, plus extreme |estimate|, plus forced labels
top_up_sig   <- df |> filter(direction == "Up")   |> slice_min(padj, n = 10)
top_down_sig <- df |> filter(direction == "Down") |> slice_min(padj, n = 10)
top_effect   <- df |> filter(direction != "NS")   |> slice_max(abs(estimate), n = 8)
forced       <- df |> filter(genes %in% FORCE_LABEL) |>
                  mutate(direction = case_when(estimate > 0 ~ factor("Up", levels = c("Down","NS","Up")),
                                                estimate < 0 ~ factor("Down", levels = c("Down","NS","Up")),
                                                TRUE         ~ direction))

label_df <- bind_rows(top_up_sig, top_down_sig, top_effect, forced) |>
  distinct(genes, .keep_all = TRUE) |>
  mutate(is_forced = genes %in% FORCE_LABEL)

cat(sprintf("Labeling %d genes (%d forced via curated list)\n",
            nrow(label_df), sum(label_df$is_forced)))
if (length(FORCE_LABEL) > 0) {
  cat("  Forced labels:", paste(FORCE_LABEL, collapse = ", "), "\n")
  missing_forced <- setdiff(FORCE_LABEL, df$genes)
  if (length(missing_forced) > 0)
    cat("  (not in DE table):", paste(missing_forced, collapse = ", "), "\n")
}

# ---- Plot ------------------------------------------------------------------
xlim_abs <- max(abs(df$estimate), na.rm = TRUE) * 1.05
ymax     <- max(df$neglog10_padj, na.rm = TRUE) * 1.05

p <- ggplot(df, aes(x = estimate, y = neglog10_padj)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = -log10(PADJ_SIG),
             colour = "grey40", linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(PADJ_STRICT),
             colour = "grey40", linetype = "dotted", linewidth = 0.4) +
  # NS background points
  geom_point(aes(colour = direction),
             data = df |> filter(direction == "NS"),
             alpha = 0.35, size = 1.1) +
  # Sig (un-forced) points
  geom_point(aes(colour = direction),
             data = df |> filter(direction != "NS"),
             alpha = 0.85, size = 1.6) +
  # Force-labeled points highlighted with a black ring
  geom_point(data = df |> filter(genes %in% FORCE_LABEL),
             shape = 21, colour = "black", fill = NA,
             size = 3, stroke = 0.7) +
  # Force-labeled gene names in green
  geom_text_repel(
    data = label_df |> filter(is_forced),
    aes(label = genes),
    colour = COL_FORCE, size = 4.4, fontface = "bold.italic",
    max.overlaps = Inf, box.padding = 0.55, point.padding = 0.3,
    segment.size = 0.4, segment.colour = COL_FORCE,
    min.segment.length = 0, seed = 7
  ) +
  # Auto-labeled gene names coloured by direction
  geom_text_repel(
    data = label_df |> filter(!is_forced),
    aes(label = genes, colour = direction),
    size = 4.0, fontface = "italic",
    max.overlaps = Inf, box.padding = 0.45, point.padding = 0.25,
    segment.size = 0.3, segment.colour = "grey50",
    min.segment.length = 0, seed = 42,
    show.legend = FALSE
  ) +
  annotate("text", x = xlim_abs * 0.98, y = -log10(PADJ_SIG),
           label = sprintf("padj = %.2g", PADJ_SIG),
           hjust = 1, vjust = -0.4, size = 3.6, colour = "grey30") +
  annotate("text", x = xlim_abs * 0.98, y = -log10(PADJ_STRICT),
           label = sprintf("padj = %.2g", PADJ_STRICT),
           hjust = 1, vjust = -0.4, size = 3.6, colour = "grey30") +
  annotate("text", x = -xlim_abs * 0.96, y = ymax * 0.98,
           label = sprintf("Down: %d", n_down),
           hjust = 0, vjust = 1, size = 5, colour = COL_DOWN, fontface = "bold") +
  annotate("text", x =  xlim_abs * 0.96, y = ymax * 0.98,
           label = sprintf("Up: %d", n_up),
           hjust = 1, vjust = 1, size = 5, colour = COL_UP,   fontface = "bold") +
  scale_colour_manual(
    values = c("Down" = COL_DOWN, "NS" = COL_NS, "Up" = COL_UP),
    breaks = c("Down", "Up"),
    labels = c("Down (estimate < 0)", "Up (estimate > 0)"),
    name = NULL
  ) +
  scale_x_continuous(limits = c(-xlim_abs, xlim_abs),
                     expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, ymax),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    x = "Effect size (meta-analytic estimate)",
    y = expression(-log[10]~"(adjusted p-value)")
  ) +
  theme_cowplot(font_size = 15) +
  theme(
    legend.position    = c(0.02, 0.02),
    legend.justification = c(0, 0),
    legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.text        = element_text(size = 12),
    axis.title         = element_text(size = 15),
    axis.text          = element_text(size = 13),
    plot.margin        = margin(10, 14, 8, 10)
  )

ggsave(OUT_FIG, p, width = 8, height = 7, dpi = 250, bg = "white")
ggsave(OUT_PDF, p, width = 8, height = 7, bg = "white")
cat(sprintf("Saved %s and %s\n", OUT_FIG, OUT_PDF))
