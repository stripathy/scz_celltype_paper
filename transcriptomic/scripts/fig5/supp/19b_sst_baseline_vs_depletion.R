#!/usr/bin/env Rscript
# Does baseline SST expression per supertype predict depletion in SCZ?
# Tomoda et al. 2022 (preproSST processing load -> ER stress -> vulnerability)
# predicts: higher baseline SST -> more depleted (negative correlation with the
# crumblr estimate). Three baselines: SEA-AD reference mean expression (the Fig
# 4e source), SEA-AD A9/DLPFC CP10K means, and Xenium control cells (12 donors).
# Depth is entangled with depletion (rho = 0.74), so partial correlations
# controlling depth are reported alongside. n = 16 supertypes -> descriptive.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(cowplot); library(ggrepel)
})
OUT <- "transcriptomic/results/sst_strata_gsea"
FS  <- 14
STRATA_LEVELS <- c("depleted", "intermediate", "non_depleted")
STRATA_LABELS <- c("Depleted (5)", "Intermediate (6)", "Non-depleted (5)")
STRAT_COLS <- c("Depleted (5)" = "#D55E00", "Intermediate (6)" = "#E69F00",
                "Non-depleted (5)" = "#0072B2")

strata <- read_csv(file.path(OUT, "strata_definition.csv"), show_col_types = FALSE) |>
  mutate(stratum = factor(stratum, STRATA_LEVELS, STRATA_LABELS))

sst_row <- function(path) {
  m <- read_csv(path, show_col_types = FALSE)
  m |> filter(.data[[colnames(m)[1]]] == "SST") |>
    select(-1) |> pivot_longer(everything(), names_to = "CellType", values_to = "sst")
}
seaad <- sst_row("genetics/results/intermediates/seaad_sst_supertype_mean_expression.csv") |>
  rename(seaad_ref = sst)
a9    <- sst_row("genetics/results/intermediates/a9_sst_supertype_mean_expression_cp10k.csv") |>
  rename(a9_cp10k = sst)
xen   <- read_csv(file.path(OUT, "xenium_sst_baseline_controls.csv"), show_col_types = FALSE) |>
  select(CellType = supertype, xen_ctrl_cp10k = xen_sst_cp10k)

d <- strata |>
  left_join(seaad, by = "CellType") |>
  left_join(a9, by = "CellType") |>
  left_join(xen, by = "CellType")
write_csv(d, file.path(OUT, "sst_baseline_vs_depletion.csv"))

partial_spearman <- function(x, y, z) {
  rx <- resid(lm(rank(x) ~ rank(z))); ry <- resid(lm(rank(y) ~ rank(z)))
  cor(rx, ry)
}
sources <- c(seaad_ref = "SEA-AD reference (mean expr)",
             a9_cp10k = "SEA-AD DLPFC (CP10K)",
             xen_ctrl_cp10k = "Xenium controls (CP10K)")
stats <- imap_dfr <- NULL
res <- lapply(names(sources), function(v) {
  ct <- cor.test(d[[v]], d$estimate, method = "spearman", exact = TRUE)
  tibble(source = sources[[v]], var = v,
         rho_depletion = unname(ct$estimate), p = ct$p.value,
         rho_depth = cor(d[[v]], d$depth_xenium, method = "spearman",
                         use = "complete.obs"),
         rho_partial_depth = partial_spearman(d[[v]], d$estimate, d$depth_xenium))
}) |> bind_rows()
cat("Baseline SST vs crumblr depletion estimate (n = 16 supertypes):\n")
cat("(Tomoda prediction: rho_depletion < 0)\n")
print(res |> mutate(across(where(is.numeric), ~round(.x, 3))))

cat("\nCross-source agreement of the baselines (Spearman):\n")
print(round(cor(d |> select(seaad_ref, a9_cp10k, xen_ctrl_cp10k),
                method = "spearman"), 2))

panel <- function(v, xlab) {
  s <- res |> filter(var == v)
  ggplot(d, aes(.data[[v]], estimate, colour = stratum)) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_errorbar(aes(ymin = estimate - se, ymax = estimate + se),
                  width = 0, linewidth = 0.5, alpha = 0.6) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = CellType), size = 3.8, seed = 1,
                    show.legend = FALSE, max.overlaps = 20) +
    annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.8, size = 4.6,
             label = sprintf("rho = %.2f, p = %.2g\npartial (| depth) = %.2f",
                             s$rho_depletion, s$p, s$rho_partial_depth)) +
    scale_colour_manual(values = STRAT_COLS, name = NULL) +
    labs(x = xlab, y = "Abundance change in SCZ (crumblr estimate ± SE)") +
    theme_cowplot(font_size = FS) +
    theme(legend.position = "none")
}
fig <- plot_grid(
  panel("seaad_ref", "Baseline SST, SEA-AD reference (mean expr)"),
  panel("a9_cp10k", "Baseline SST, SEA-AD DLPFC (CP10K)") +
    theme(axis.title.y = element_blank()),
  panel("xen_ctrl_cp10k", "Baseline SST, Xenium controls (CP10K)") +
    theme(axis.title.y = element_blank()),
  nrow = 1, rel_widths = c(1.06, 1, 1))
leg <- get_legend(panel("seaad_ref", "") +
                    theme(legend.position = "bottom") +
                    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3))))
fig <- plot_grid(fig, leg, ncol = 1, rel_heights = c(1, 0.08))
ggsave(file.path(OUT, "figS_sst_baseline_vs_depletion.png"), fig,
       width = 15, height = 5.6, dpi = 200, bg = "white")
cat("\nWrote figS_sst_baseline_vs_depletion.png\n")
