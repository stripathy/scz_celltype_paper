#!/usr/bin/env Rscript
# Supplementary figure | The graded dysregulation is not a consequence of the
# depleted groups containing fewer nuclei from cases.
#
# This is the figure cited in the results as (Fig. SXXX) after the sentence on
# downsampling control nuclei. Three panels, in the order the argument runs:
#   a  the confound exists and is itself graded, and matching removes it
#   b  gene-set burden per group, unmatched versus five matched draws
#
# Inputs: subsample_cells.py (matched pseudobulks) and
# sensitivity_cellsubsample.R (their DE and GSEA).
source("transcriptomic/scripts/fig5/_common.R")
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })

SUB  <- file.path(P$supp, "cellsub")
reps <- sort(as.integer(str_remove(basename(Sys.glob(file.path(SUB, "rep*"))), "rep")))
S_LAB <- c(depleted = "Depleted", intermediate = "Intermediate", non_depleted = "Non-depleted")
lab_group <- function(x) factor(S_LAB[x], unname(S_LAB))
# panel a is a horizontal forest, so reverse so Depleted reads at the top
lab_group_rev <- function(x) factor(S_LAB[x], rev(unname(S_LAB)))

# ---- a | the imbalance, before and after matching ----------------------------
# Case-control difference in nuclei per donor, within cohort. Before matching the
# difference is large and graded exactly like the result, which is the reason the
# control is needed at all.
imb <- function(meta, state) {
  map_dfr(STRATA_LEVELS, function(st) {
    d <- meta |> filter(stratum == st) |> mutate(l = log(n_cells))
    co <- summary(lm(l ~ diagnosis + cohort, data = d))$coefficients["diagnosisSCZ", ]
    tibble(stratum = st, state = state, est = co[1], se = co[2], p = co[4])
  })
}
before <- all_donor_meta()
after  <- map_dfr(export_cohorts(), ~ read_csv(
  file.path(SUB, "rep1", paste0(.x, "_donor_stratum_meta.csv")),
  col_types = cols(donor = col_character(), .default = col_guess())))
d_imb <- bind_rows(imb(before, "Observed"), imb(after, "Cell-matched")) |>
  mutate(state = factor(state, c("Observed", "Cell-matched")),
         group = lab_group_rev(stratum),
         pct = 100 * (exp(est) - 1),
         lo = 100 * (exp(est - 1.96 * se) - 1), hi = 100 * (exp(est + 1.96 * se) - 1))
say("panel a, case-control difference in nuclei per donor, percent:")
print(d_imb |> transmute(group, state, pct = round(pct, 1), p = signif(p, 2)))

p_a <- ggplot(d_imb, aes(pct, group, colour = state)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.6 * GEO,
                 position = position_dodge(width = 0.55)) +
  geom_point(size = 2.4 * GEO, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c(Observed = "grey35", `Cell-matched` = "#B2182B"), name = NULL) +
  labs(x = "Nuclei per donor, SCZ vs control (%)", y = NULL) +
  theme_fig5() +
  theme(legend.position = c(0.02, 0.14), legend.justification = c(0, 0),
        legend.background = element_blank())

# ---- b | gene-set burden, unmatched vs matched -------------------------------
gsea_un <- read_csv(file.path(P$pb, "gsea_all_signatures.csv"), show_col_types = FALSE)
sub <- read_csv(file.path(P$supp, "sensitivity_cellsub_gsea.csv"), show_col_types = FALSE)
b_un <- gsea_un |> filter(signature %in% STRATA_LEVELS) |> group_by(stratum = signature) |>
  summarise(n = sum(padj < 0.10), .groups = "drop") |> mutate(group = lab_group(stratum))
b_mt <- sub |> group_by(stratum, rep) |> summarise(n = sum(padj < 0.10), .groups = "drop") |>
  mutate(group = lab_group(stratum))

p_b <- ggplot(b_un, aes(group, n)) +
  geom_col(fill = "grey82", width = 0.62) +
  geom_point(data = b_mt, colour = "#B2182B", size = 1.9 * GEO, alpha = 0.85,
             position = position_jitter(width = 0.10, height = 0, seed = 1)) +
  labs(x = NULL, y = "Gene sets at FDR < 0.10") +
  theme_fig5() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

fig <- plot_grid(p_a, p_b, nrow = 1, rel_widths = c(0.58, 0.42), align = "h", axis = "tb",
                 labels = c("a", "b"), label_size = PANEL_LABEL)
ggsave(file.path(P$supp, "figS_cellcount_control.png"), fig, width = FIG_W, height = 2.6,
       dpi = FIG_DPI, bg = "white")
ggsave(file.path(P$supp, "figS_cellcount_control.pdf"), fig, width = FIG_W, height = 2.6,
       bg = "white")
say("wrote figS_cellcount_control.(png|pdf)")
say("matched draws: %s", paste(reps, collapse = ", "))
