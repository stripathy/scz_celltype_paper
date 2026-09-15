#!/usr/bin/env Rscript
# Supplement | Module scores as single per-donor measurements.
#
# WHY THIS EXISTS. No individual gene in the translation or OxPhos modules
# reaches gene-level FDR within a stratum, so "these programs are suppressed"
# rests on gene-set enrichment. This panel answers the obvious objection without
# any gene-level multiplicity at all: collapse each module to ONE number per
# donor-stratum, then test that number. Two things become visible that a
# heatmap cannot show -- the per-cohort estimates (is one dataset driving it?)
# and the effect size in interpretable SD units.
#
# a  diagnosis effect on each module score, per stratum, with all 7 cohorts shown
# b  the same modules tested WITHIN donor (depleted - non-depleted), which is the
#    paired form of the interaction and cancels donor-level confounds
#
# Deliberately a forest, not a violin: the quantity being tested is the
# meta-analytic effect, and per-cohort consistency is the point.
source("transcriptomic/scripts/sst_strata/_common.R")
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })

tests  <- read_csv(file.path(P$pb, "module_score_tests.csv"), show_col_types = FALSE)
percoh <- read_csv(file.path(P$pb, "module_score_percohort.csv"), show_col_types = FALSE)

MOD_LEV <- c("synaptic", "translation", "oxphos")
MOD_LAB <- c(synaptic = "Synaptic\n(shared)", translation = "Translation\n(graded)",
             oxphos = "OxPhos/mito\n(graded)")
ROW_LEV <- c("depleted", "intermediate", "non_depleted", "interaction")
ROW_LAB <- c(depleted = "Depleted", intermediate = "Intermediate",
             non_depleted = "Non-depleted",
             interaction = "Depleted − non-depleted\n(within donor)")
ROW_COLS <- c(STRAT_COLS, interaction = "grey25")

prep <- function(d, is_test) {
  d |> mutate(stratum = ifelse(is.na(stratum) | stratum == "interaction",
                               "interaction", stratum),
              row = factor(ROW_LAB[stratum], rev(unname(ROW_LAB))),
              module = factor(MOD_LAB[module], unname(MOD_LAB[MOD_LEV])),
              col = ROW_COLS[stratum])
}
tt <- prep(tests |> mutate(stratum = ifelse(test == "within_donor_interaction",
                                           "interaction", stratum)), TRUE) |>
  mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se,
         lab = sprintf("%.2f  P=%s  %d/%d↓", estimate,
                       format(signif(pval, 2), scientific = TRUE), n_neg, k))
pc <- prep(percoh, FALSE)

p <- ggplot(tt, aes(estimate, row)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
  # every cohort's own estimate, so consistency is visible rather than asserted
  geom_point(data = pc, aes(x = est, y = row), colour = "grey55",
             size = 1.5, alpha = 0.75,
             position = position_nudge(y = 0.20)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = col), height = 0, linewidth = 0.9) +
  geom_point(aes(colour = col), size = 3.6, shape = 18) +
  geom_text(aes(label = lab), y = as.numeric(tt$row) - 0.30, size = 3.5,
            colour = "grey20", hjust = 0.5) +
  facet_wrap(~module, nrow = 1) +
  scale_colour_identity() +
  scale_y_discrete(labels = function(x) x) +
  labs(x = "Diagnosis effect on module score (SD units, SCZ − control)", y = NULL,
       caption = paste("Diamonds, random-effects meta-analysis across the seven",
                       "cohorts (whiskers, 95% CI); grey points, per-cohort estimates.",
                       "\nn/7↓ = cohorts with a negative estimate. Module score = mean of",
                       "within-cohort z-scored log-CPM over the module's leading-edge genes.")) +
  theme_cowplot(font_size = FS) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold", size = FS - 1),
        axis.text.y = element_text(size = FS - 2, lineheight = 0.9),
        plot.caption = element_text(size = FS - 4, colour = "grey35", hjust = 0))

ggsave(file.path(P$supp, "figS_module_scores.png"), p, width = 14.5, height = 6.4,
       dpi = 200, bg = "white")
say("wrote figS_module_scores.png")
print(tt |> select(module, row, estimate, se, pval, n_neg, k))
