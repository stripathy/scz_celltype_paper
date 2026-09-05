# Supplementary figure: the SCZ common-variant enrichment landscape behind Fig. 4a,
# for the 125 supertypes of the SEA-AD DLPFC taxonomy (Bigdeli 2026 GWAS, SEA-AD
# DLPFC expression reference; MAGMA gene-property, one-sided).
#
# Supersedes the retired 501-type landscape: the combined SEA-AD + Siletti taxonomy
# was retired (L. Duncan's advice), so the Siletti panel is gone and the
# significance lines are corrected over the 125 supertypes actually tested.
#
# Points grouped by subclass (Fig. 3a order, then non-neuronal), colored with the
# SEA-AD supertype palette of Figs 1b/3a/4; the five Sst supertypes depleted in
# SCZ outlined in black; supertypes passing Bonferroni labeled.
#
# Input : results/tables/scz_enrichment_seaad125_bigdeli_dlpfc.csv
#         (scripts/figures/export_supp_enrichment_seaad125.py)
# Output: manuscript/figures/supplementary/S09_scz_enrichment_seaad125.{png,pdf}
#         (the submission set -- see that folder's README; the S-number lives
#         in FIGSTEM below)
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(cowplot); library(ggrepel)
})

# Repo-relative: resolve this script's own location, so it runs from any clone.
# Rscript exposes the path via --file=; the fallback assumes the repo root cwd.
.script_dir <- function() {
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else NA_character_
}
.sd <- .script_dir()
GEN   <- normalizePath(if (!is.na(.sd)) file.path(.sd, "..", "..") else "genetics")
d     <- read_csv(file.path(GEN, "results/tables/scz_enrichment_seaad125_bigdeli_dlpfc.csv"),
                  show_col_types = FALSE)
BASE  <- 7
LAB   <- 1.9          # threshold annotations (~5.4 pt)
LAB_A <- 1.55         # point labels (~4.4 pt), plain, no repel
PT    <- 1.5
N     <- nrow(d)
BONF  <- -log10(0.05 / N)
FDR05 <- -log10(max(d$p[d$fdr < 0.05]))
COL_BONF <- "#6a3d9a"; COL_FDR <- "grey55"
cat(sprintf("N = %d; Bonferroni line %.2f; FDR line %.2f\n", N, BONF, FDR05))

SEAAD_ORDER <- c("Lamp5_Lhx6","Lamp5","Pax6","Sncg","Vip","Sst Chodl","Sst","Pvalb",
                 "Chandelier","L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L5 ET",
                 "L6 CT","L6b","L5/6 NP","Astro","Oligo","OPC","Micro-PVM","Endo","VLMC")

a <- d |>
  mutate(group = factor(group, levels = SEAAD_ORDER)) |>
  arrange(group, desc(neglog10p)) |>
  mutate(x = row_number(),
         label = ifelse(bonferroni < 0.05, cell_type, "")) |>
  group_by(group) |> mutate(strip_col = first(color)) |> ungroup()
ymax_a <- max(a$neglog10p) * 1.20

strip <- a |> group_by(group) |>
  summarise(xmin = min(x) - 0.45, xmax = max(x) + 0.45, xmid = mean(x),
            col = first(strip_col), .groups = "drop")

pa <- ggplot(a, aes(x, neglog10p)) +
  geom_hline(yintercept = BONF, linetype = "dotdash", colour = COL_BONF, linewidth = 0.35) +
  geom_hline(yintercept = FDR05, linetype = "dashed", colour = COL_FDR, linewidth = 0.3) +
  # points as in Fig. 4 (SEA-AD palette fill; thick outline = depleted in SCZ)
  geom_point(aes(fill = color, colour = depleted_scz, stroke = depleted_scz),
             shape = 21, size = PT) +
  # vertical labels above the points, repelled sideways only, so labels of
  # adjacent same-subclass points (Sst_2/Sst_20/Sst_23) do not overprint
  geom_text_repel(data = filter(a, label != ""),
                  aes(y = neglog10p, label = label),
                  angle = 90, hjust = 0, vjust = 0.5, size = LAB_A, colour = "grey15",
                  direction = "x", nudge_y = 0.14, box.padding = 0.05,
                  point.padding = 0, segment.size = 0.15, segment.colour = "grey60",
                  min.segment.length = 0.25, max.overlaps = Inf, seed = 1) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey30"), guide = "none") +
  scale_discrete_manual("stroke", values = c(`TRUE` = 0.7, `FALSE` = 0.2), guide = "none") +
  annotate("text", x = nrow(a) + 0.5, y = BONF, label = "Bonferroni 0.05", hjust = 1,
           vjust = -0.4, size = LAB, colour = COL_BONF) +
  annotate("text", x = nrow(a) + 0.5, y = FDR05, label = "FDR 0.05", hjust = 1,
           vjust = 1.4, size = LAB, colour = COL_FDR) +
  # subclass strip under the axis: a bar per subclass plus a rotated name
  geom_rect(data = strip, aes(xmin = xmin, xmax = xmax, ymin = -0.75, ymax = -0.25,
                              fill = col),
            inherit.aes = FALSE, colour = "white", linewidth = 0.2) +
  geom_text(data = strip, aes(x = xmid, y = -0.95, label = group), inherit.aes = FALSE,
            angle = 90, hjust = 1, vjust = 0.5, size = LAB_A, colour = "grey15") +
  scale_fill_identity() +
  coord_cartesian(ylim = c(0, ymax_a), clip = "off") +
  scale_x_continuous(expand = expansion(add = 0.8)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(y = expression("SCZ GWAS enrichment ("*-log[10]~italic(P)*")")) +
  theme_cowplot(font_size = BASE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_text(size = BASE - 0.5),
        axis.text.y = element_text(size = BASE - 1),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(4, 6, 52, 4))

SUPPFIG <- file.path(GEN, "..", "manuscript/figures/supplementary")
FIGSTEM <- "S09_scz_enrichment_seaad125"
for (ext in c("png", "pdf"))
  ggsave(file.path(SUPPFIG, paste0(FIGSTEM, ".", ext)), pa,
         width = 7.1, height = 3.2, dpi = 400, bg = "white")
cat(sprintf("wrote %s/%s.{png,pdf}\n", SUPPFIG, FIGSTEM))
