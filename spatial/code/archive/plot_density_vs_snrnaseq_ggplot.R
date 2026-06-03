#!/usr/bin/env Rscript
#
# Scatter: snRNAseq meta-analysis betas vs Xenium density logFC
# Uses ggplot2 + cowplot theme + ggrepel labels
#
# Input:
#   output/density_analysis/snrnaseq_vs_density_comparison.csv
#
# Output:
#   output/density_analysis/snrnaseq_vs_density_supertype_ggplot.png

library(ggplot2)
library(cowplot)
library(ggrepel)

base_dir <- path.expand("~/Github/SCZ_Xenium")
in_path <- file.path(base_dir, "output", "density_analysis",
                     "snrnaseq_vs_density_comparison.csv")
out_path <- file.path(base_dir, "output", "density_analysis",
                      "snrnaseq_vs_density_supertype_ggplot.png")

# ── Load data ──
df <- read.csv(in_path, stringsAsFactors = FALSE)
df <- df[complete.cases(df$beta_snrnaseq, df$logFC_density), ]
cat(sprintf("Loaded %d supertypes\n", nrow(df)))

# ── Class colors ──
class_colors <- c(
  "Glut" = "#00ADF8",
  "GABA" = "#F05A28",
  "NN"   = "#808080"
)
df$class <- factor(df$class, levels = c("Glut", "GABA", "NN"))

# ── Correlations ──
r_all <- cor.test(df$beta_snrnaseq, df$logFC_density)
neur <- df[df$class %in% c("Glut", "GABA"), ]
r_neur <- cor.test(neur$beta_snrnaseq, neur$logFC_density)

subtitle_text <- sprintf(
  "All: r = %.3f (p = %.1e)  |  Neuronal: r = %.3f (p = %.1e)  |  n = %d",
  r_all$estimate, r_all$p.value,
  r_neur$estimate, r_neur$p.value,
  nrow(df)
)

# ── Label selection: FDR-sig in either, large effects, or nom-sig in both ──
df$label <- ""
for (i in seq_len(nrow(df))) {
  is_sig_either <- (df$padj_snrnaseq[i] < 0.05) | (!is.na(df$fdr_density[i]) & df$fdr_density[i] < 0.10)
  is_large <- (abs(df$beta_snrnaseq[i]) > 0.15) | (abs(df$logFC_density[i]) > 0.5)
  is_nom_both <- (df$pval_snrnaseq[i] < 0.05) & (df$pval_density[i] < 0.05)
  # Always label specific types of interest
  always_label <- c("Sst_20", "L2/3 IT_8", "L2/3 IT_10")
  if (is_sig_either | (is_large & is_nom_both) | (abs(df$logFC_density[i]) > 0.5) |
      (df$celltype[i] %in% always_label)) {
    df$label[i] <- df$celltype[i]
  }
}

# ── Regression line coefficients ──
fit <- lm(logFC_density ~ beta_snrnaseq, data = df)

# ── Count per class for legend ──
class_n <- table(df$class)
class_labels <- sprintf("%s (n=%d)", names(class_n), as.integer(class_n))
names(class_labels) <- names(class_n)

# ── Plot ──
p <- ggplot(df, aes(x = beta_snrnaseq, y = logFC_density, color = class)) +
  # Error bars (snRNAseq SE horizontal)
  geom_errorbarh(aes(xmin = beta_snrnaseq - se,
                     xmax = beta_snrnaseq + se),
                 height = 0, alpha = 0.15, linewidth = 0.4) +
  # Error bars (density SE vertical)
  geom_errorbar(aes(ymin = logFC_density - se_density,
                    ymax = logFC_density + se_density),
                width = 0, alpha = 0.15, linewidth = 0.4) +
  # Reference lines
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  # Regression line
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2],
              color = "grey30", linewidth = 0.8, linetype = "solid", alpha = 0.4) +
  # Points
  geom_point(size = 2.5, alpha = 0.8) +
  # Labels
  geom_text_repel(
    aes(label = label),
    size = 3, max.overlaps = 30,
    segment.color = "grey50", segment.size = 0.3,
    segment.alpha = 0.5,
    min.segment.length = 0.1,
    box.padding = 0.4, point.padding = 0.2,
    color = "grey20",
    seed = 42
  ) +
  # Colors
  scale_color_manual(
    values = class_colors,
    labels = class_labels
  ) +
  # Axis labels
  labs(
    x = "snRNAseq meta-analysis beta (SCZ effect)",
    y = "Xenium density logFC (SCZ vs Control)\n(cells/mm², cortical)",
    title = "SCZ Supertype Effects: snRNAseq vs Xenium Density",
    subtitle = subtitle_text,
    color = "Class"
  ) +
  # Theme
  theme_cowplot(font_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    legend.position.inside = c(0.82, 0.12),
    legend.position = "inside",
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    legend.key.size = unit(0.5, "cm")
  )

ggsave(out_path, p, width = 9, height = 8, dpi = 150)
cat(sprintf("Saved: %s\n", out_path))
