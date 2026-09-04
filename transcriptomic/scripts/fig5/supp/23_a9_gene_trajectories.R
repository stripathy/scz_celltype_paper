#!/usr/bin/env Rscript
# Gabitto-style LOESS trajectories: per-donor stratum-pseudobulk expression
# along CPS for exemplar genes, colored by our depletion strata. Shows the
# inversion directly: ribosomal/OxPhos genes decline along CPS in the spared
# strata but stay flat in the depleted stratum; SST/VGF decline everywhere.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(arrow); library(ggplot2); library(cowplot)
})
OUT  <- "transcriptomic/results/sst_strata_gsea"
AOUT <- file.path(OUT, "seaad_a9")
FS <- 14
GENES <- c("SST", "VGF", "RPL4", "UQCRH")
STRATA_LEVELS <- c("depleted", "intermediate", "non_depleted")
STRATA_LABELS <- c("Depleted", "Intermediate", "Non-depleted")
COLS <- c(Depleted = "#D55E00", Intermediate = "#E69F00", `Non-depleted` = "#0072B2")

pb <- as.data.frame(read_parquet(file.path(AOUT, "a9_stratum_pseudobulk.parquet")))
idx <- grep("__index_level_0__", colnames(pb), value = TRUE)
rownames(pb) <- pb[[idx]]; pb[[idx]] <- NULL; pb <- as.matrix(pb)
meta <- read_csv(file.path(AOUT, "a9_stratum_pseudobulk_meta.csv"), show_col_types = FALSE) |>
  mutate(severe = as.character(severe) %in% c("Y", "True", "TRUE", "1")) |>
  filter(!severe, n_cells >= 10)

cpm <- t(t(pb[, meta$key]) / colSums(pb[, meta$key]) * 1e6)
d <- as_tibble(t(log1p(cpm[GENES, , drop = FALSE])), rownames = "key") |>
  left_join(meta |> select(key, stratum, CPS), by = "key") |>
  pivot_longer(all_of(GENES), names_to = "gene", values_to = "expr") |>
  mutate(stratum = factor(stratum, STRATA_LEVELS, STRATA_LABELS),
         gene = factor(gene, GENES))

p <- ggplot(d, aes(CPS, expr, colour = stratum, fill = stratum)) +
  geom_point(size = 0.9, alpha = 0.4) +
  geom_smooth(method = "loess", span = 1, se = TRUE, alpha = 0.18, linewidth = 1.1) +
  geom_vline(xintercept = 0.6, linetype = "dashed", colour = "grey65") +
  facet_wrap(~gene, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = COLS, name = NULL) +
  scale_fill_manual(values = COLS, name = NULL) +
  labs(x = "Continuous pseudo-progression score (CPS)",
       y = "Stratum pseudobulk expression, log1p(CPM)") +
  theme_cowplot(font_size = FS) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "italic", size = FS + 1))
ggsave(file.path(AOUT, "figS_a9_gene_trajectories.png"), p,
       width = 15, height = 4.8, dpi = 200, bg = "white")
cat("wrote figS_a9_gene_trajectories.png\n")
