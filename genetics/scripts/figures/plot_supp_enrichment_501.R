# Supplementary figure: the SCZ common-variant enrichment landscape behind Fig. 4a,
# for all 501 types of the combined SEA-AD + Siletti taxonomy (Bigdeli 2026 GWAS,
# SEA-AD DLPFC expression reference; MAGMA gene-property, one-sided).
#
#   a  the 125 SEA-AD supertypes (points), grouped by subclass (Fig. 3a order, then
#      non-neuronal), colored with the SEA-AD supertype palette of Figs 1b/3a/4;
#      the five Sst supertypes depleted in SCZ outlined in black
#   b  the 376 Siletti-only clusters (points), grouped by Siletti supercluster
# Both rows share the significance lines (Bonferroni 0.05 over 501 tests; FDR 0.05).
#
# Input : results/tables/scz_enrichment_501_bigdeli_dlpfc.csv
#         (scripts/figures/export_supp_enrichment_501.py)
# Output: results/figures/supp_scz_enrichment_501.{png,pdf}
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(cowplot); library(ggrepel)
})

GEN   <- "/Users/shreejoy/Github/scz_celltype_paper/genetics"
d     <- read_csv(file.path(GEN, "results/tables/scz_enrichment_501_bigdeli_dlpfc.csv"),
                  show_col_types = FALSE)
BASE  <- 7
LAB   <- 1.9          # geom_text size (~5.4 pt)
LAB_A <- 1.55         # point labels in a (~4.4 pt), plain, no repel
PT    <- 1.3          # point size, both rows
N     <- nrow(d)
BONF  <- -log10(0.05 / N)
FDR05 <- -log10(max(d$p[d$fdr < 0.05]))
COL_BONF <- "#6a3d9a"; COL_FDR <- "grey55"
cat(sprintf("N = %d; Bonferroni line %.2f; FDR line %.2f\n", N, BONF, FDR05))

SEAAD_ORDER <- c("Lamp5_Lhx6","Lamp5","Pax6","Sncg","Vip","Sst Chodl","Sst","Pvalb",
                 "Chandelier","L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L5 ET",
                 "L6 CT","L6b","L5/6 NP","Astro","Oligo","OPC","Micro-PVM","Endo","VLMC")

theme_land <- function() theme_cowplot(font_size = BASE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_text(size = BASE - 0.5),
        axis.text.y = element_text(size = BASE - 1),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(4, 6, 2, 4))

# Group strip under the axis: a bar per group plus a rotated group name.
group_strip <- function(df, ymax, y_bar, y_lab, size = LAB, repel = TRUE) {
  g <- df |> group_by(group) |>
    summarise(xmin = min(x) - 0.45, xmax = max(x) + 0.45, xmid = mean(x),
              col = first(strip_col), .groups = "drop")
  list(
    geom_rect(data = g, aes(xmin = xmin, xmax = xmax, ymin = y_bar[1], ymax = y_bar[2], fill = col),
              inherit.aes = FALSE, colour = "white", linewidth = 0.2),
    # vertical group names; repelled sideways only where tiny adjacent groups
    # would otherwise collide (b), plain text in a
    if (repel) geom_text_repel(data = g, aes(x = xmid, y = y_lab, label = group), inherit.aes = FALSE,
                    angle = 90, hjust = 1, vjust = 0.5, size = size, colour = "grey15",
                    direction = "x", ylim = c(-Inf, y_lab), box.padding = 0.08,
                    point.padding = 0, segment.size = 0.15, segment.colour = "grey60",
                    min.segment.length = 0.3, max.overlaps = Inf, seed = 1)
    else geom_text(data = g, aes(x = xmid, y = y_lab, label = group), inherit.aes = FALSE,
                   angle = 90, hjust = 1, vjust = 0.5, size = size, colour = "grey15"),
    scale_fill_identity(),
    coord_cartesian(ylim = c(0, ymax), clip = "off"))
}
thresholds <- list(
  geom_hline(yintercept = BONF, linetype = "dotdash", colour = COL_BONF, linewidth = 0.35),
  geom_hline(yintercept = FDR05, linetype = "dashed", colour = COL_FDR, linewidth = 0.3))

# ---------------------------------------------------------------- a: SEA-AD
a <- d |> filter(source == "SEA-AD") |>
  mutate(group = factor(group, levels = SEAAD_ORDER)) |>
  arrange(group, desc(neglog10p)) |>
  mutate(x = row_number(),
         label = ifelse(bonferroni < 0.05, cell_type, ""),
         strip_col = color)
# subclass strip color = the median-lightness member colour is fine; use the first
a <- a |> group_by(group) |> mutate(strip_col = first(color)) |> ungroup()
ymax_a <- max(a$neglog10p) * 1.28
pa <- ggplot(a, aes(x, neglog10p)) +
  thresholds +
  # points, as in panel b and Fig. 4 (SEA-AD palette fill; thick outline = depleted in SCZ)
  geom_point(aes(fill = color, colour = depleted_scz, stroke = depleted_scz),
             shape = 21, size = PT) +
  # plain vertical labels above the points (no repel); smaller text so
  # neighbouring labels do not collide
  geom_text(data = filter(a, label != ""), aes(y = neglog10p + 0.22, label = label),
            angle = 90, hjust = 0, vjust = 0.5, size = LAB_A, colour = "grey15") +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey30"), guide = "none") +
  scale_discrete_manual("stroke", values = c(`TRUE` = 0.7, `FALSE` = 0.2), guide = "none") +
  annotate("text", x = nrow(a) + 0.5, y = BONF, label = "Bonferroni 0.05", hjust = 1,
           vjust = -0.4, size = LAB, colour = COL_BONF) +
  annotate("text", x = nrow(a) + 0.5, y = FDR05, label = "FDR 0.05", hjust = 1,
           vjust = -0.4, size = LAB, colour = COL_FDR) +
  group_strip(a, ymax_a, y_bar = c(-0.9, -0.3), y_lab = -1.15, size = LAB_A, repel = FALSE) +
  scale_x_continuous(expand = expansion(add = 0.8)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(y = expression("SCZ GWAS enrichment ("*-log[10]~italic(P)*")")) +
  theme_land() + theme(plot.margin = margin(4, 6, 40, 4))

# --------------------------------------------------------------- b: Siletti
SHORT <- c("Committed oligodendrocyte precursor" = "Committed OPC",
           "Oligodendrocyte precursor" = "OPC",
           "Deep-layer corticothalamic and 6b" = "Deep-layer CT and 6b",
           "Upper-layer intratelencephalic" = "Upper-layer IT",
           "Deep-layer intratelencephalic" = "Deep-layer IT",
           "Deep-layer near-projecting" = "Deep-layer NP",
           "Eccentric medium spiny neuron" = "Eccentric MSN",
           "Medium spiny neuron" = "MSN",
           "Hippocampal dentate gyrus" = "Hippocampal DG",
           "Midbrain-derived inhibitory" = "Midbrain inhibitory",
           "LAMP5-LHX6 and Chandelier" = "LAMP5-LHX6/Chandelier")
b <- d |> filter(source == "Siletti") |>
  mutate(group = dplyr::recode(group, !!!SHORT)) |>
  group_by(group) |> mutate(gmax = max(neglog10p)) |> ungroup() |>
  arrange(desc(gmax), group, desc(neglog10p)) |>
  mutate(group = factor(group, levels = unique(group)),
         x = row_number(),
         gi = as.integer(group),
         strip_col = ifelse(gi %% 2 == 1, "grey45", "grey75"),
         label = ifelse(neglog10p >= 8, cell_type, ""))
ymax_b <- max(b$neglog10p) * 1.22
pb <- ggplot(b, aes(x, neglog10p)) +
  thresholds +
  geom_point(aes(fill = strip_col), shape = 21, size = PT, stroke = 0.15, colour = "grey20") +
  geom_text_repel(data = filter(b, label != ""), aes(y = neglog10p, label = label),
                  angle = 90, hjust = 0, vjust = 0.5, size = LAB, colour = "grey15",
                  direction = "both", nudge_y = 0.6, ylim = c(NA, Inf), box.padding = 0.1,
                  point.padding = 0.1, segment.size = 0.15, segment.colour = "grey60",
                  min.segment.length = 0.2, max.overlaps = Inf, seed = 1) +
  group_strip(b, ymax_b, y_bar = c(-1.5, -0.5), y_lab = -2.0) +
  scale_x_continuous(expand = expansion(add = 2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(y = expression("SCZ GWAS enrichment ("*-log[10]~italic(P)*")")) +
  theme_land() + theme(plot.margin = margin(4, 6, 84, 4))

fig <- plot_grid(pa, pb, ncol = 1, rel_heights = c(1, 1.25), labels = c("a", "b"),
                 label_size = 8, label_fontface = "bold")
for (ext in c("png", "pdf"))
  ggsave(file.path(GEN, "results/figures", paste0("supp_scz_enrichment_501.", ext)), fig,
         width = 7.1, height = 6.3, dpi = 400, bg = "white")
cat("wrote results/figures/supp_scz_enrichment_501.{png,pdf}\n")
