#!/usr/bin/env Rscript
# Cell-typing validation supplement: Xenium vs SEA-AD MERFISH concordance.
# 2x2 grid -- proportions (top) and median depth (bottom), subclass (left) /
# supertype (right). Reads the flat CSVs built by
# build_celltyping_validation_data.py (one symmetric recipe for all panels).
#
# Conventions: proportions on log-log axes; depth on linear 0-1 axes; identity
# dashed line; points colored by class; one shared legend; plain-text r (Pearson)
# and rho (Spearman) in a corner (no boxed inset, no in-plot title).

suppressPackageStartupMessages({
  library(ggplot2); library(cowplot); library(ggrepel); library(scales)
})
cowplot::set_null_device("agg")  # unicode-capable null device for alignment pass

# ---- constants ----
.args <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- if (length(.args) > 0) .args[1] else "output/celltyping_supplement/data"
# output figure basename: 2nd arg or default
OUT_DIR <- file.path(dirname(DATA_DIR))
FIG_BASE <- if (length(.args) > 1) .args[2] else "celltyping_validation_supplement"

CLASS_COLORS <- c("Glutamatergic" = "#E65100",
                  "GABAergic"     = "#2E7D32",
                  "Non-neuronal"  = "#1565C0")
PT_SIZE <- 2.6
LAB_SIZE <- 2.5
BASE_FS <- 9

theme_panel <- function() {
  theme_classic(base_size = BASE_FS) +
    theme(axis.title = element_text(size = BASE_FS),
          axis.text = element_text(size = BASE_FS - 1, color = "black"),
          plot.margin = margin(6, 8, 4, 6),
          legend.position = "none",
          aspect.ratio = 1)
}

stat_label <- function(r, rho, n, r_tag = "") {
  sprintf("r = %.2f%s\nρ = %.2f\nn = %d", r, r_tag, rho, n)
}

# ---- proportion panel (log-log) ----
prop_panel <- function(df, xlab, ylab, label_all = FALSE, n_label = 10, note = "") {
  df$lr <- abs(log10((df$xenium_prop) / (df$merfish_prop)))
  df$lab <- if (label_all) df$celltype else
            ifelse(rank(-df$lr) <= n_label, df$celltype, "")
  r_log  <- cor(log10(df$merfish_prop), log10(df$xenium_prop))
  rho    <- cor(df$merfish_prop, df$xenium_prop, method = "spearman")
  lims <- range(c(df$merfish_prop, df$xenium_prop))
  p <- ggplot(df, aes(merfish_prop, xenium_prop, fill = klass)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey55", linewidth = 0.4) +
    geom_point(shape = 21, size = PT_SIZE, color = "white", stroke = 0.3, alpha = 0.95) +
    geom_text_repel(aes(label = lab), size = LAB_SIZE, color = "grey20",
                    max.overlaps = 30, min.segment.length = 0.1,
                    segment.size = 0.2, segment.color = "grey70", seed = 1) +
    annotate("text", x = lims[1], y = lims[2], hjust = 0, vjust = 1,
             label = stat_label(r_log, rho, nrow(df), " (log)"),
             size = 2.9, lineheight = 0.95) +
    scale_x_log10(labels = label_percent(accuracy = 0.1)) +
    scale_y_log10(labels = label_percent(accuracy = 0.1)) +
    scale_fill_manual(values = CLASS_COLORS, name = NULL) +
    labs(x = xlab, y = ylab) +
    theme_panel()
  if (nzchar(note)) p <- p + annotate("text", x = lims[2], y = lims[1], hjust = 1, vjust = 0,
                                       label = note, size = 2.7, fontface = "italic", color = "grey35")
  p
}

# ---- depth panel (linear 0-1) ----
depth_panel <- function(df, xlab, ylab, label_all = FALSE, n_label = 8, note = "") {
  df$dev <- abs(df$xenium_depth - df$merfish_depth)
  df$lab <- if (label_all) df$celltype else
            ifelse(rank(-df$dev) <= n_label, df$celltype, "")
  r   <- cor(df$merfish_depth, df$xenium_depth)
  rho <- cor(df$merfish_depth, df$xenium_depth, method = "spearman")
  p <- ggplot(df, aes(merfish_depth, xenium_depth, fill = klass)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey55", linewidth = 0.4) +
    geom_point(shape = 21, size = PT_SIZE, color = "white", stroke = 0.3, alpha = 0.95) +
    geom_text_repel(aes(label = lab), size = LAB_SIZE, color = "grey20",
                    max.overlaps = 30, min.segment.length = 0.1,
                    segment.size = 0.2, segment.color = "grey70", seed = 1) +
    annotate("text", x = 0.97, y = 0.03, hjust = 1, vjust = 0,
             label = stat_label(r, rho, nrow(df)), size = 2.9, lineheight = 0.95) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_fill_manual(values = CLASS_COLORS, name = NULL) +
    labs(x = xlab, y = ylab) +
    theme_panel()
  if (nzchar(note)) p <- p + annotate("text", x = 0.03, y = 0.97, hjust = 0, vjust = 1,
                                       label = note, size = 2.7, fontface = "italic", color = "grey35")
  p
}

# ---- load ----
ps <- read.csv(file.path(DATA_DIR, "prop_subclass.csv"))
pt <- read.csv(file.path(DATA_DIR, "prop_supertype.csv"))
ds <- read.csv(file.path(DATA_DIR, "depth_subclass.csv"))
dt <- read.csv(file.path(DATA_DIR, "depth_supertype.csv"))

pa <- prop_panel(ps, "MERFISH proportion", "Xenium proportion", label_all = TRUE)
# supertype panels (b, d): neurons only -- supertype resolution is used only for
# neuronal analyses in this study; non-neuronal cells are analyzed at subclass level.
pt_neu <- pt[pt$klass != "Non-neuronal", ]
dt_neu <- dt[dt$klass != "Non-neuronal", ]
pb <- prop_panel(pt_neu, "MERFISH proportion", "Xenium proportion", label_all = FALSE,
                 n_label = 10)
pc <- depth_panel(ds, "MERFISH median depth (manual)", "Xenium median depth (predicted)", label_all = TRUE)
pdp <- depth_panel(dt_neu, "MERFISH median depth (manual)", "Xenium median depth (predicted)",
                   label_all = FALSE, n_label = 8)

# shared legend (from a panel with legend turned on)
leg <- get_legend(pa + theme(legend.position = "bottom") +
                  guides(fill = guide_legend(override.aes = list(size = 3.5))))

grid <- plot_grid(pa, pb, pc, pdp, ncol = 2, labels = c("a", "b", "c", "d"),
                  label_size = 13, label_fontface = "bold", align = "hv",
                  hjust = -0.2, vjust = 1.4)
fig <- plot_grid(grid, leg, ncol = 1, rel_heights = c(1, 0.06))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
# png (preview) + pdf (vector, editable text) -- the convention used across the
# paper's supplemental figures.
png_path <- file.path(OUT_DIR, paste0(FIG_BASE, ".png"))
pdf_path <- file.path(OUT_DIR, paste0(FIG_BASE, ".pdf"))
ggsave(png_path, fig, width = 7.2, height = 7.8, dpi = 400, bg = "white",
       device = ragg::agg_png)
ggsave(pdf_path, fig, width = 7.2, height = 7.8, bg = "white")
cat("Saved:", png_path, "\n")
cat("Saved:", pdf_path, "\n")
