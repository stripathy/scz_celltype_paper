#!/usr/bin/env Rscript
# Module-score analysis: turn each module (translation, OxPhos, synaptic) into
# ONE variable per donor-stratum (mean of within-cohort z-scored log-CPM over
# module genes), then test it like a biomarker:
#   (1) per stratum: score ~ dx + covars per cohort -> metafor across cohorts
#       (single test per module x stratum; no gene-level multiplicity)
#   (2) interaction: within-donor difference (depleted - non-depleted score)
#       ~ dx + covars per cohort -> metafor (paired, donor confounds cancel)
# Also reports per-cohort direction consistency (x of 7 cohorts negative).
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(arrow); library(metafor)
})
OUT  <- "transcriptomic/results/sst_strata_gsea"
POUT <- file.path(OUT, "pseudobulk")
DOUT <- file.path(POUT, "donor_stratum")
EXP  <- "transcriptomic/data/stratum_pseudobulks_export"

g <- read_csv(file.path(OUT, "gsea_all_signatures.csv"), show_col_types = FALSE)
le_union <- function(pws) g |> filter(pathway %in% pws, signature == "depleted") |>
  pull(leadingEdge) |> strsplit("|", fixed = TRUE) |> unlist() |> unique()
MODS <- list(
  translation = le_union(c("REACTOME_TRANSLATION", "GOCC_RIBOSOMAL_SUBUNIT",
                           "GOBP_CYTOPLASMIC_TRANSLATION", "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION")),
  oxphos = le_union(c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "GOBP_OXIDATIVE_PHOSPHORYLATION",
                      "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT")),
  synaptic = le_union(c("GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING", "GOCC_SYNAPTIC_MEMBRANE",
                        "GOBP_NEUROTRANSMITTER_SECRETION", "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES")))
cat(sprintf("module sizes: %s\n",
            paste(names(MODS), lengths(MODS), collapse = ", ")))

cohorts <- str_remove(basename(Sys.glob(file.path(DOUT, "*_donor_stratum_meta.csv"))),
                      "_donor_stratum_meta\\.csv$")
covars <- map_dfr(cohorts, function(co)
  read_csv(file.path(EXP, paste0(co, "_groups.csv")),
           col_types = cols(donor = col_character(), .default = col_guess())) |>
    distinct(donor, sex, age, pmi) |>
    mutate(cohort = co,
           age = suppressWarnings(as.numeric(as.character(age))),
           pmi = suppressWarnings(as.numeric(as.character(pmi)))))

scores <- map_dfr(cohorts, function(co) {
  pb <- as.data.frame(read_parquet(file.path(DOUT, paste0(co, "_donor_stratum_counts.parquet"))))
  rn <- pb$gene; pb$gene <- NULL
  M <- vapply(pb, as.numeric, numeric(nrow(pb))); rownames(M) <- rn
  m <- read_csv(file.path(DOUT, paste0(co, "_donor_stratum_meta.csv")),
                col_types = cols(donor = col_character(), .default = col_guess()))
  keys <- paste(m$donor, m$stratum, sep = "|")
  lc <- log1p(t(t(M[, keys]) / colSums(M[, keys])) * 1e6)
  keep <- rowSums(lc > 0) >= 0.8 * ncol(lc)          # expressed genes only
  lcz <- t(scale(t(lc[keep, ])))                     # z per gene within cohort
  map_dfr(names(MODS), function(mod) {
    gg <- intersect(MODS[[mod]], rownames(lcz))
    tibble(cohort = co, donor = m$donor, stratum = m$stratum,
           diagnosis = m$diagnosis, module = mod,
           score = colMeans(lcz[gg, , drop = FALSE]), n_genes = length(gg))
  })
}) |> left_join(covars, by = c("cohort", "donor"))
write_csv(scores, file.path(POUT, "module_scores_donor.csv"))

fit_meta <- function(d) {   # d: one cohort-level estimate per row (est, se)
  fit <- tryCatch(rma(yi = d$est, sei = d$se, method = "REML",
                      control = list(maxiter = 1000)),
                  error = function(e) rma(yi = d$est, sei = d$se, method = "DL"))
  tibble(est = fit$beta[1], se = fit$se, z = fit$zval, p = fit$pval, k = fit$k,
         n_neg = sum(d$est < 0))
}
cohort_lm <- function(d, formula_rhs) {
  d <- d |> mutate(dx = factor(diagnosis, c("Control", "SCZ")), sex = factor(sex),
                   age_s = scale(ifelse(is.na(age), median(age, na.rm = TRUE), age))[, 1])
  use_pmi <- mean(is.na(d$pmi)) < 0.5
  if (use_pmi) d$pmi_s <- scale(ifelse(is.na(d$pmi), median(d$pmi, na.rm = TRUE), d$pmi))[, 1]
  f <- as.formula(paste("y ~ dx + age_s + sex", if (use_pmi) "+ pmi_s" else ""))
  fit <- lm(f, data = d)
  co <- summary(fit)$coefficients
  tibble(est = co["dxSCZ", 1], se = co["dxSCZ", 2])
}

cat("\n=== (1) module score ~ diagnosis, per stratum (meta over 7 cohorts) ===\n")
res1 <- scores |>
  group_by(module, stratum) |>
  group_modify(~ {
    percoh <- .x |> group_by(cohort) |>
      group_modify(~ cohort_lm(.x |> rename(y = score))) |> ungroup()
    fit_meta(percoh)
  }) |> ungroup() |>
  mutate(across(c(est, se, z), ~ round(.x, 3)), p = signif(p, 2))
print(res1 |> arrange(module, stratum), n = 12)

cat("\n=== (2) within-donor interaction: (depleted - non-depleted score) ~ dx ===\n")
res2 <- scores |>
  select(cohort, donor, stratum, diagnosis, module, score, sex, age, pmi) |>
  pivot_wider(names_from = stratum, values_from = score) |>
  filter(!is.na(depleted), !is.na(non_depleted)) |>
  mutate(y = depleted - non_depleted) |>
  group_by(module) |>
  group_modify(~ {
    percoh <- .x |> group_by(cohort) |> group_modify(~ cohort_lm(.x)) |> ungroup()
    fit_meta(percoh)
  }) |> ungroup() |>
  mutate(across(c(est, se, z), ~ round(.x, 3)), p = signif(p, 2))
print(res2, n = 5)
cat("\nn_neg = cohorts (of 7) with negative estimate\n")
