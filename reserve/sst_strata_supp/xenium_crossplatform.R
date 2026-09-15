#!/usr/bin/env Rscript
# Supplement | Cross-platform check of the stratum-level differential expression.
#
# WHY THIS EXISTS. Two jobs, one supporting and one disclosing.
#   a  The shared component replicates on an independent platform and cohort:
#      Xenium is whole-cell (snRNA-seq is nucleus-restricted), a different gene
#      panel, a different classifier, and different donors (LIBD, 24 DLPFC
#      sections). SST and VGF are reduced in all three strata on both platforms.
#   b  The graded components CANNOT be tested this way, and the figure should say
#      so in its own panel rather than only in the text: a 300-gene cell-typing
#      panel carries 1 of 130 translation and 2 of 114 OxPhos leading-edge genes.
#
# Inputs come from 07_xenium_stratum.R, which documents the stratum aggregation.
source("transcriptomic/scripts/sst_strata/_common.R")
suppressPackageStartupMessages({ library(ggplot2); library(cowplot); library(ggrepel) })

conc <- read_csv(file.path(P$out, "xenium_stratum_concordance.csv"), show_col_types = FALSE)
g    <- read_csv(file.path(P$pb, "gsea_all_signatures.csv"), show_col_types = FALSE)
mods <- modules(g)
S_LAB <- c(depleted = "Depleted", intermediate = "Intermediate", non_depleted = "Non-depleted")

# ---- a | per-stratum agreement across platforms ------------------------------
d <- conc |> mutate(stratum = factor(S_LAB[stratum], unname(S_LAB)))
stats <- d |> group_by(stratum) |>
  summarise(n = n(), r = cor(sn_z, xen_z),
            pct = 100 * mean(sign(sn_z) == sign(xen_z)), .groups = "drop") |>
  mutate(lab = sprintf("r = %.2f\n%.0f%% same sign\nn = %d", r, pct, n))
hl <- d |> filter(gene %in% c("SST", "VGF"))

p_a <- ggplot(d, aes(sn_z, xen_z)) +
  geom_hline(yintercept = 0, colour = "grey88") +
  geom_vline(xintercept = 0, colour = "grey88") +
  geom_abline(linetype = "dashed", colour = "grey65") +
  geom_point(size = 1.5, alpha = 0.45, colour = "grey40") +
  geom_point(data = hl, colour = "black", fill = COL_DOWN, shape = 21, size = 3, stroke = 0.8) +
  geom_text_repel(data = hl, aes(label = gene), size = 4, fontface = "italic",
                  seed = 2, min.segment.length = 0, box.padding = 0.6) +
  geom_text(data = stats, aes(x = -Inf, y = Inf, label = lab), hjust = -0.12, vjust = 1.15,
            size = 3.7, colour = "grey25", lineheight = 0.95, inherit.aes = FALSE) +
  facet_wrap(~stratum, nrow = 1) +
  labs(x = "snRNA-seq meta-analytic z (SCZ vs control)",
       y = "Xenium stratum-level z") +
  theme_cowplot(font_size = FS) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"))

# ---- b | why the graded modules are untestable on this panel -----------------
panel <- unique(conc$gene)
cov <- tibble(module = names(mods),
              on_panel = sapply(mods, function(m) sum(m %in% panel)),
              total = lengths(mods)) |>
  filter(module %in% c("synaptic", "translation", "oxphos")) |>
  mutate(pct = 100 * on_panel / total,
         label = factor(FAM_LABELS[module], unname(FAM_LABELS)),
         txt = sprintf("%d / %d", on_panel, total))

p_b <- ggplot(cov, aes(pct, label, fill = label)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = txt), hjust = -0.15, size = 4.2, colour = "grey20") +
  scale_fill_manual(values = setNames(unname(FAM_COLS), unname(FAM_LABELS)), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(x = "Leading-edge genes on the 300-gene Xenium panel (%)", y = NULL,
       subtitle = "The graded modules are essentially absent from the panel") +
  theme_cowplot(font_size = FS) +
  theme(plot.subtitle = element_text(size = FS - 3, colour = "grey35"))

fig <- plot_grid(p_a, p_b, ncol = 1, rel_heights = c(1, 0.62),
                 labels = c("a", "b"), label_size = 20)
ggsave(file.path(P$supp, "figS_xenium_crossplatform.png"), fig, width = 13, height = 8.4,
       dpi = 200, bg = "white")
say("wrote figS_xenium_crossplatform.png")
print(stats |> select(stratum, n, r, pct))
print(cov |> select(module, on_panel, total, pct))
