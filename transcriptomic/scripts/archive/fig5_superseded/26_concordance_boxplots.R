#!/usr/bin/env Rscript
# Concordance v3: exemplar boxplots. Top row, shared genes (SST, VGF): SCZ =
# Jens depleted-stratum donor pseudobulks (CTRL vs SCZ), AD = A9 depleted-
# stratum donor pseudobulks binned by CPS tertile. Bottom row, divergent
# MODULE SCORES (translation, OxPhos): SCZ side uses the Jens AFFECTED stratum
# (depleted+intermediate — the boundary Jens labels draw reliably; single
# divergent genes have no signal in jens-depleted due to label porosity),
# AD side = A9 depleted stratum. Interim panel: swap SCZ side to Nicole's
# 7-cohort pseudobulks (depleted stratum proper) when they land.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(arrow); library(ggplot2); library(cowplot)
})
OUT  <- "transcriptomic/results/sst_strata_gsea"
AOUT <- file.path(OUT, "seaad_a9")
JOUT <- file.path(OUT, "jens")
FS <- 13
SCZ_COLS <- c(CTRL = "#0072B2", SCZ = "#D55E00")
CPS_COLS <- c("#c6c6c6", "#8a8a8a", "#4d4d4d")

g <- read_csv(file.path(OUT, "gsea_all_signatures.csv"), show_col_types = FALSE)
le_union <- function(pws) g |> filter(pathway %in% pws, signature == "depleted") |>
  pull(leadingEdge) |> strsplit("|", fixed = TRUE) |> unlist() |> unique()
MODS <- list(
  "Translation module" = le_union(c("REACTOME_TRANSLATION", "GOCC_RIBOSOMAL_SUBUNIT",
                                    "GOBP_CYTOPLASMIC_TRANSLATION",
                                    "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION")),
  "OxPhos module" = le_union(c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                               "GOBP_OXIDATIVE_PHOSPHORYLATION",
                               "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT")))

load_pb <- function(path_pb, path_meta) {
  pb <- as.data.frame(read_parquet(path_pb))
  idx <- grep("__index_level_0__", colnames(pb), value = TRUE)
  rownames(pb) <- pb[[idx]]; pb[[idx]] <- NULL
  list(pb = as.matrix(pb),
       meta = read_csv(path_meta, show_col_types = FALSE) |> filter(n_cells >= 10))
}
jens <- load_pb(file.path(JOUT, "jens_stratum_pseudobulk.parquet"),
                file.path(JOUT, "jens_stratum_pseudobulk_meta.csv"))
a9 <- load_pb(file.path(AOUT, "a9_stratum_pseudobulk.parquet"),
              file.path(AOUT, "a9_stratum_pseudobulk_meta.csv"))
a9$meta <- a9$meta |>
  mutate(severe = as.character(severe) %in% c("Y", "True", "TRUE", "1")) |>
  filter(!severe, stratum == "depleted") |>
  mutate(bin = cut(CPS, quantile(CPS, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
                   labels = c("low", "mid", "high")))
cat("A9 CPS tertile ranges:\n")
print(a9$meta |> group_by(bin) |> summarise(n = n(), min = min(CPS), max = max(CPS)))

logcpm <- function(pb, keys) {
  x <- pb[, keys, drop = FALSE]
  log1p(t(t(x) / colSums(x)) * 1e6)
}
val <- function(pb, keys, genes) {
  x <- logcpm(pb, keys)
  if (length(genes) == 1) x[genes, ]
  else colMeans(x[intersect(genes, rownames(x)), , drop = FALSE])
}

panel_pair <- function(title, genes, jens_stratum, note = NULL) {
  jm <- jens$meta |> filter(stratum %in% jens_stratum) |>
    group_by(donor_id) |> summarise(diagnosis = first(diagnosis),
                                    keys = list(key), .groups = "drop")
  jv <- map_dbl(jm$keys, function(k) {
    v <- val(jens$pb, k, genes); if (length(k) > 1) mean(v) else v })
  dj <- tibble(grp = factor(jm$diagnosis, c("CTRL", "SCZ")), y = jv)
  da <- tibble(grp = a9$meta$bin,
               y = val(a9$pb, a9$meta$key, genes))
  p_s <- ggplot(dj, aes(grp, y, colour = grp)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.55) +
    geom_jitter(width = 0.13, size = 0.8, alpha = 0.5) +
    scale_colour_manual(values = SCZ_COLS, guide = "none") +
    labs(x = NULL, y = "log1p(CPM)",
         title = title, subtitle = "SCZ (Jens)") +
    theme_cowplot(font_size = FS) +
    theme(plot.title = element_text(face = ifelse(length(genes) == 1, "italic", "bold"),
                                    size = FS + 2),
          plot.subtitle = element_text(size = FS - 1, colour = "grey30"))
  p_a <- ggplot(da, aes(grp, y, colour = grp)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.55) +
    geom_jitter(width = 0.13, size = 0.8, alpha = 0.6) +
    scale_colour_manual(values = setNames(CPS_COLS, levels(da$grp)), guide = "none") +
    labs(x = "CPS tertile", y = NULL, title = "", subtitle = "AD (SEA-AD A9)") +
    theme_cowplot(font_size = FS) +
    theme(plot.subtitle = element_text(size = FS - 1, colour = "grey30"))
  if (!is.null(note))
    p_s <- p_s + labs(caption = note) +
      theme(plot.caption = element_text(size = FS - 3, colour = "grey40"))
  plot_grid(p_s, p_a, nrow = 1, rel_widths = c(0.45, 0.55), align = "h")
}

row1 <- plot_grid(panel_pair("SST", "SST", "depleted"),
                  panel_pair("VGF", "VGF", "depleted"), nrow = 1)
row2 <- plot_grid(
  panel_pair("RPL36", "RPL36", "depleted"),
  panel_pair("NDUFS8", "NDUFS8", "depleted"),
  nrow = 1)
foot <- ggdraw() + draw_label(
  paste0("PLACEHOLDER SCZ side (bottom row): the Jens depleted-stratum labels cannot resolve this signal ",
         "(7-cohort meta z: RPL36 = -3.8, NDUFS8 = -3.3) — rebuild and verify when the 7-cohort pseudobulks land."),
  size = FS - 2, colour = "grey35", x = 0.02, hjust = 0)
hdr1 <- ggdraw() + draw_label("Shared between SCZ and AD (depleted stratum)",
                              fontface = "bold", size = FS + 2, x = 0.02, hjust = 0)
hdr2 <- ggdraw() + draw_label("Divergent: down in SCZ (meta), flat in AD depleted cells",
                              fontface = "bold", size = FS + 2, x = 0.02, hjust = 0)
fig <- plot_grid(hdr1, row1, hdr2, row2, foot, ncol = 1,
                 rel_heights = c(0.07, 1, 0.07, 1, 0.06))
ggsave(file.path(AOUT, "figS_gene_concordance_v3_boxplots.png"), fig,
       width = 12.5, height = 8.6, dpi = 200, bg = "white")
cat("wrote figS_gene_concordance_v3_boxplots.png\n")
