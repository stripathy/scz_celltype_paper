#!/usr/bin/env Rscript
# Step 6 | Module scores: each module as ONE measurement per donor-stratum.
#
# WHY. No individual module gene reaches gene-level FDR within a stratum, which
# makes "translation is suppressed" uncomfortable to assert from GSEA alone. This
# collapses each module to a single variable per donor-stratum (mean of
# within-cohort z-scored log-CPM over the module's genes) and tests it like a
# biomarker, with no gene-level multiplicity at all:
#   (1) per stratum: score ~ dx + covars per cohort -> REML meta across cohorts
#   (2) interaction: within-donor (depleted - non_depleted) score ~ dx + covars
#       -> REML meta. Paired, so donor confounds cancel.
# Also reports per-cohort direction consistency (how many of 7 are negative),
# which is the honest answer to "is this driven by one dataset?".
#
# Runtime ~4 min.
source("transcriptomic/scripts/sst_strata/_common.R")

gsea <- read_csv(file.path(P$pb, "gsea_all_signatures.csv"), show_col_types = FALSE)
MODS <- modules(gsea)[c("translation", "oxphos", "synaptic")]
say("module sizes: %s", paste(sprintf("%s %d", names(MODS), lengths(MODS)), collapse = ", "))

covars <- load_covars()
scores <- map_dfr(donor_cohorts(), function(co) {
  d <- load_donor_stratum(co)
  m <- d$meta
  keys <- paste(m$donor, m$stratum, sep = "|")
  lc <- log1p(t(t(d$counts[, keys]) / colSums(d$counts[, keys])) * 1e6)
  lc <- lc[rowSums(lc > 0) >= 0.8 * ncol(lc), ]     # expressed genes only
  lcz <- t(scale(t(lc)))                            # z per gene within cohort
  map_dfr(names(MODS), function(mod) {
    gg <- intersect(MODS[[mod]], rownames(lcz))
    tibble(cohort = co, donor = m$donor, stratum = m$stratum, diagnosis = m$diagnosis,
           module = mod, score = colMeans(lcz[gg, , drop = FALSE]), n_genes = length(gg))
  })
}) |> left_join(covars, by = c("cohort", "donor"))
write_csv(scores, file.path(P$pb, "module_scores_donor.csv"))

# per-cohort effect of diagnosis on one score column
cohort_lm <- function(d) {
  d <- prep_model_frame(d)
  f <- if (attr(d, "use_pmi")) y ~ dx + age_s + sex + pmi_s else y ~ dx + age_s + sex
  co <- summary(lm(f, data = d))$coefficients
  tibble(est = co["dxSCZ", 1], se = co["dxSCZ", 2])
}
meta_of <- function(percoh) fit_meta(percoh$est, percoh$se) |>
  mutate(n_neg = sum(percoh$est < 0))

# per-cohort estimates are kept as well: the supplement forest plot shows them,
# and "negative in 7 of 7 cohorts" is the answer to "is one dataset driving this?"
percoh_all <- list()

say("=== (1) module score ~ diagnosis, per stratum (meta over 7 cohorts) ===")
res1 <- scores |> group_by(module, stratum) |>
  group_modify(~ {
    pc <- .x |> rename(y = score) |> group_by(cohort) |>
      group_modify(~ cohort_lm(.x)) |> ungroup()
    percoh_all[[length(percoh_all) + 1]] <<- pc |>
      mutate(module = .y$module, stratum = .y$stratum, test = "per_stratum")
    meta_of(pc)
  }) |> ungroup() |>
  mutate(across(c(estimate, se, zval), ~ round(.x, 3)), pval = signif(pval, 2)) |>
  arrange(module, match(stratum, STRATA_LEVELS))
print(res1, n = 12)

say("=== (2) within-donor interaction: (depleted - non_depleted) score ~ dx ===")
res2 <- scores |>
  select(cohort, donor, stratum, diagnosis, module, score, sex, age, pmi) |>
  pivot_wider(names_from = stratum, values_from = score) |>
  filter(!is.na(depleted), !is.na(non_depleted)) |>
  mutate(y = depleted - non_depleted) |>
  group_by(module) |>
  group_modify(~ {
    pc <- .x |> group_by(cohort) |> group_modify(~ cohort_lm(.x)) |> ungroup()
    percoh_all[[length(percoh_all) + 1]] <<- pc |>
      mutate(module = .y$module, stratum = "interaction", test = "within_donor_interaction")
    meta_of(pc)
  }) |> ungroup() |>
  mutate(across(c(estimate, se, zval), ~ round(.x, 3)), pval = signif(pval, 2))
print(res2, n = 5)
write_csv(bind_rows(res1 |> mutate(test = "per_stratum"),
                    res2 |> mutate(test = "within_donor_interaction", stratum = NA)),
          file.path(P$pb, "module_score_tests.csv"))
write_csv(bind_rows(percoh_all), file.path(P$pb, "module_score_percohort.csv"))
say("n_neg = cohorts (of 7) with a negative estimate")
