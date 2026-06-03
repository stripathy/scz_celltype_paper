#!/usr/bin/env Rscript
# ============================================================================
# 08_meta_vs_xenium_scatter.R — Cross-platform concordance scatter
# ============================================================================
# snRNA-seq meta-analytic logFC (x) vs Xenium spatial logFC (y) for every
# (gene, cell type) pair where the gene is BOTH meta-DE (padj < PADJ_THR) AND
# in the Xenium 300-gene panel AND tested in the matching subclass. Quantifies
# how well snRNA-seq discovery DE replicates on an independent spatial platform
# (Pearson r, % directionally concordant, binomial test for >50%).
#
# Two output styles per run:
#   08_meta_vs_xenium_scatter[_fdrNN]           dot size = meta FDR bin
#                                               (< 0.05 large, 0.05-0.10 small)
#   08b_meta_vs_xenium_scatter_uniform[_fdrNN]  all dots same size
# Both colour by cell class and label a curated gene set (FOREST_PAIRS + extras).
#
# Usage:
#   Rscript scripts/08_meta_vs_xenium_scatter.R         # PADJ_THR = 0.10
#   Rscript scripts/08_meta_vs_xenium_scatter.R 0.05    # stricter cutoff
# Non-default thresholds get an _fdrNN filename suffix so cutoffs coexist.
#
# DATA PROVENANCE — update if the DE results change:
#   INPUT_META   Meta-analytic snRNA-seq DE (metafor). Used cols:
#                genes, cell_type, estimate, pval, padj.
#   INPUT_XENIUM Xenium spatial pseudobulk DE (edgeR QL F-test), Kwon 2026,
#                12 SCZ vs 12 control DLPFC, 300-gene panel. Used cols:
#                gene, logFC, PValue, FDR, celltype.
# Xenium subclass names are harmonized to the SEA-AD convention via ct_map.
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(cowplot); library(ggrepel)
})

INPUT_META   <- "data/DE_genes_all_cells_scz.csv"
INPUT_XENIUM <- "~/Github/SCZ_Xenium/output/de/de_results_subclass.csv"
FIG_DIR      <- "results"

# snRNA-seq meta FDR threshold. Optional CLI arg (default 0.1). Non-default
# thresholds get a filename suffix so outputs at different cutoffs coexist.
.args    <- commandArgs(trailingOnly = TRUE)
PADJ_THR <- if (length(.args) >= 1) as.numeric(.args[[1]]) else 0.1
SUFFIX   <- if (abs(PADJ_THR - 0.1) < 1e-9) "" else
              sprintf("_fdr%02d", round(PADJ_THR * 100))

EXC <- c("L2_3 IT","L4 IT","L5 IT","L5 ET","L5_6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")
INH <- c("Lamp5","Pax6","Pvalb","Sncg","Sst","Sst Chodl","Vip","Chandelier")
GLI <- c("Astro","Oligo","OPC","Micro-PVM","Endo","VLMC")
class_of <- function(ct) dplyr::case_when(
  ct %in% EXC ~ "Excitatory", ct %in% INH ~ "Inhibitory",
  ct %in% GLI ~ "Glia",       TRUE        ~ NA_character_)
CLASS_COL <- c(Excitatory = "#117733", Inhibitory = "#882255", Glia = "#DDCC77")

cat("Loading meta-analysis (filtering padj < ", PADJ_THR, ")...\n", sep="")
meta <- read_csv(INPUT_META, show_col_types = FALSE) |>
  filter(padj < PADJ_THR) |>
  transmute(genes, cell_type, meta_est = estimate, meta_padj = padj)

cat("Loading + harmonizing Xenium DE...\n")
ct_map <- c("Astrocyte"="Astro","L2/3 IT"="L2_3 IT","L5/6 NP"="L5_6 NP",
            "Microglia-PVM"="Micro-PVM","Oligodendrocyte"="Oligo",
            "Endothelial"="Endo")
xen <- read_csv(INPUT_XENIUM, show_col_types = FALSE) |>
  mutate(cell_type = ifelse(celltype %in% names(ct_map),
                            ct_map[celltype], celltype)) |>
  transmute(genes = gene, cell_type, xen_logFC = logFC, xen_p = PValue)

# Merge — inner join gives the testable (gene x cell) pairs
pairs <- inner_join(meta, xen, by = c("genes", "cell_type")) |>
  mutate(class = class_of(cell_type),
         concordant = sign(meta_est) == sign(xen_logFC))

cat(sprintf("Testable pairs (meta padj<%g + gene in Xenium panel + cell-type tested): %d\n",
            PADJ_THR, nrow(pairs)))
cat(sprintf("Directionally concordant: %d (%.1f%%)\n",
            sum(pairs$concordant), 100*mean(pairs$concordant)))

# Statistics
pear_r <- cor(pairs$meta_est, pairs$xen_logFC, method = "pearson")
spear_r <- cor(pairs$meta_est, pairs$xen_logFC, method = "spearman")
bt <- binom.test(sum(pairs$concordant), nrow(pairs), p = 0.5, alternative = "greater")
fit <- lm(xen_logFC ~ meta_est, data = pairs)
slope <- coef(fit)[2]; rsq <- summary(fit)$r.squared

cat(sprintf("Pearson r:  %.3f\nSpearman r: %.3f\nLinear-fit slope: %.2f (R²=%.2f)\nBinomial test for >50%% concordance: p = %.2g\n",
            pear_r, spear_r, slope, rsq, bt$p.value))

# ---- Explicit label set (few, legible) -------------------------------------
# Two-level snRNA-seq meta FDR bin drives dot size (all pairs are < PADJ_THR).
pairs <- pairs |>
  mutate(tag = sprintf("%s (%s)", genes, cell_type),
         meta_fdr_bin = factor(ifelse(meta_padj < 0.05, "< 0.05", "0.05-0.10"),
                               levels = c("< 0.05", "0.05-0.10")))

# The 5 detailed-forest-panel pairs that replicate significantly (Fig. forest)
FOREST_PAIRS <- tibble::tribble(
  ~genes,    ~cell_type,
  "SST",     "Sst",
  "BDNF",    "L2_3 IT",
  "FKBP5",   "OPC",
  "CX3CR1",  "Micro-PVM",
  "SMAD1",   "Pvalb"
)
# A few additional pairs spanning the other validated stories
EXTRA_PAIRS <- tibble::tribble(
  ~genes,     ~cell_type,
  "SERPING1", "Astro",      # astrocyte inflammation UP
  "ABCG2",    "Endo",       # BBB transporter DOWN
  "VGF",      "Lamp5"       # neuropeptide DOWN
)
keep_pairs <- bind_rows(FOREST_PAIRS, EXTRA_PAIRS)

label_df <- pairs |> inner_join(keep_pairs, by = c("genes", "cell_type"))
missing <- anti_join(keep_pairs, pairs, by = c("genes", "cell_type"))
if (nrow(missing) > 0)
  cat("WARNING: requested label pairs not in testable set:\n",
      paste(sprintf("  %s (%s)", missing$genes, missing$cell_type),
            collapse = "\n"), "\n")
cat(sprintf("Labeling %d points (5 forest panels + %d extras)\n",
            nrow(label_df), nrow(EXTRA_PAIRS)))

# Flag which rows get a text label. We pass the FULL dataset to ggrepel
# (with empty labels elsewhere) so it repels text away from ALL 166 points,
# not just the labelled ones.
label_keys <- paste(label_df$genes, label_df$cell_type)
pairs <- pairs |>
  mutate(rep_label = ifelse(paste(genes, cell_type) %in% label_keys, tag, ""),
         # Directional nudges fan the dense lower-left cluster into empty margins
         # and lift the central SMAD1 label; repel then fine-tunes to dodge points.
         nx = dplyr::case_when(
           rep_label == "" ~ 0,
           genes == "CX3CR1" ~ -0.34,   # -> left margin
           genes == "SST"    ~ -0.06,
           genes == "ABCG2"  ~ -0.02,
           genes == "SMAD1"  ~ -0.18,   # up-left into open space
           genes == "BDNF"   ~  0.10,
           TRUE ~ 0),
         ny = dplyr::case_when(
           rep_label == "" ~ 0,
           genes == "SST"    ~  0.22,   # up, out of the cluster
           genes == "CX3CR1" ~  0.04,
           genes == "ABCG2"  ~ -0.22,   # down, out of the cluster
           genes == "SMAD1"  ~  0.26,   # up
           genes == "BDNF"   ~ -0.14,
           TRUE ~ 0))

# ---- Plot (publication quality) --------------------------------------------
# Crop the view to where the data lives (most pairs are near the origin).
# coord_cartesian zooms WITHOUT dropping points from the fit / stats.
lim <- max(
  as.numeric(quantile(abs(c(pairs$meta_est, pairs$xen_logFC)), 0.97)),
  max(abs(c(label_df$meta_est, label_df$xen_logFC)))   # keep labelled points inside
) * 1.08
n_outside <- sum(abs(pairs$meta_est) > lim | abs(pairs$xen_logFC) > lim)
cat(sprintf("View cropped to +/-%.2f; %d/%d points outside the view (still in fit + stats)\n",
            lim, n_outside, nrow(pairs)))

# Quadrant breakdown (reported for the figure caption; not drawn on the plot)
cat(sprintf("Quadrants: both-up=%d  both-down=%d  discordant=%d\n",
            sum(pairs$meta_est > 0 & pairs$xen_logFC > 0),
            sum(pairs$meta_est < 0 & pairs$xen_logFC < 0),
            sum(!pairs$concordant)))

# ---- Plot builder: two variants ---------------------------------------------
# uniform = FALSE : dot size = snRNA-seq meta FDR bin (< 0.05 vs 0.05-0.10),
#                   colour = cell class, all circles
# uniform = TRUE  : every point same shape + size, colour = cell class only
build_scatter <- function(uniform) {

  point_layers <- if (uniform) {
    list(
      geom_point(aes(colour = class), shape = 16, size = 4.3, alpha = 0.8),
      geom_point(data = label_df, shape = 21, fill = NA, colour = "black",
                 size = 5.6, stroke = 1.1, show.legend = FALSE)
    )
  } else {
    list(
      geom_point(aes(colour = class, size = meta_fdr_bin),
                 shape = 16, alpha = 0.8),
      # black outline marks the labelled points (same size as underlying dot)
      geom_point(data = label_df, aes(size = meta_fdr_bin), shape = 21,
                 fill = NA, colour = "black", stroke = 1.1, show.legend = FALSE)
    )
  }

  extra_scales <- if (uniform) {
    list()   # colour scale only
  } else {
    list(
      scale_size_manual(values = c("< 0.05" = 7, "0.05-0.10" = 3.2),
                        name = "snRNA-seq meta FDR")
    )
  }

  guide_spec <- if (uniform) {
    guides(colour = guide_legend(override.aes = list(size = 6)))
  } else {
    guides(colour = guide_legend(override.aes = list(size = 6), order = 1),
           size   = guide_legend(order = 2))
  }

  ggplot(pairs, aes(x = meta_est, y = xen_logFC)) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.5) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.5) +
    geom_abline(slope = 1, intercept = 0,
                colour = "grey75", linetype = "dashed", linewidth = 0.6) +
    geom_smooth(method = "lm", se = TRUE, colour = "black",
                fill = "grey80", linewidth = 1.0, formula = y ~ x) +
    point_layers +
    geom_text_repel(data = pairs, aes(label = rep_label),
                    size = 6.5, fontface = "plain", colour = "grey10",
                    max.overlaps = Inf, seed = 7,
                    nudge_x = pairs$nx, nudge_y = pairs$ny,
                    box.padding = 0.7, point.padding = 0.35,
                    force = 6, force_pull = 0.08,
                    xlim = c(-lim, lim), ylim = c(-lim, lim),
                    max.iter = 1e5, max.time = 3,
                    segment.size = 0.4, segment.colour = "grey45",
                    min.segment.length = 0) +
    annotate("text", x = lim*0.97, y = -lim*0.66,
             label = sprintf("italic(r)=='%.2f'", pear_r), parse = TRUE,
             hjust = 1, vjust = 1, size = 8, colour = "black") +
    annotate("text", x = lim*0.97, y = -lim*0.82,
             label = sprintf("'%d%% concordant'", round(100*mean(pairs$concordant))),
             parse = TRUE, hjust = 1, vjust = 1, size = 6.5, colour = "black") +
    scale_colour_manual(values = CLASS_COL, name = NULL) +
    extra_scales +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    labs(x = "snRNA-seq meta-analysis  log2 fold change",
         y = "Xenium spatial  log2 fold change") +
    guide_spec +
    theme_cowplot(font_size = 24) +
    theme(legend.position    = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.text        = element_text(size = 20),
          legend.spacing.y   = unit(2, "pt"),
          legend.key.height  = unit(22, "pt"),
          axis.title         = element_text(size = 24),
          axis.text          = element_text(size = 21),
          axis.line          = element_line(linewidth = 0.9),
          axis.ticks         = element_line(linewidth = 0.9),
          aspect.ratio       = 1,
          plot.margin        = margin(14, 18, 12, 12))
}

# Original (size = significance, shape = concordance)
p_sized <- build_scatter(uniform = FALSE)
ggsave(file.path(FIG_DIR, sprintf("08_meta_vs_xenium_scatter%s.png", SUFFIX)),
       p_sized, width = 12, height = 12, dpi = 300, bg = "white")
ggsave(file.path(FIG_DIR, sprintf("08_meta_vs_xenium_scatter%s.pdf", SUFFIX)),
       p_sized, width = 12, height = 12, bg = "white")
cat(sprintf("\nSaved results/08_meta_vs_xenium_scatter%s.{png,pdf}\n", SUFFIX))

# Uniform glyph (same shape + size; colour = cell class only)
p_uniform <- build_scatter(uniform = TRUE)
ggsave(file.path(FIG_DIR, sprintf("08b_meta_vs_xenium_scatter_uniform%s.png", SUFFIX)),
       p_uniform, width = 12, height = 12, dpi = 300, bg = "white")
ggsave(file.path(FIG_DIR, sprintf("08b_meta_vs_xenium_scatter_uniform%s.pdf", SUFFIX)),
       p_uniform, width = 12, height = 12, bg = "white")
cat(sprintf("Saved results/08b_meta_vs_xenium_scatter_uniform%s.{png,pdf}\n", SUFFIX))

# Save the table of pairs for downstream use
write_csv(pairs |> arrange(meta_padj),
          file.path(FIG_DIR, sprintf("08_meta_vs_xenium_pairs%s.csv", SUFFIX)))
cat(sprintf("Saved results/08_meta_vs_xenium_pairs%s.csv (%d rows)\n", SUFFIX, nrow(pairs)))
