# Composition pooling-strategy sensitivity for the snRNA-seq case/control analysis.
#
# Per-dataset crumblr models (as in the paper) are pooled three ways, and compared
# with "mega" models that put all 469 donors into one crumblr/dream fit:
#   FE meta   fixed-effect meta-analysis of per-dataset estimates   (paper's choice)
#   RE meta   random-effects (REML) meta-analysis
#   mega-fixed   ~ dx + age + sex + PMI + dataset
#   mega-RI      ~ dx + age + sex + PMI + (1 | dataset)
#   mega-RS      ~ dx + age + sex + PMI + (1 + dx | dataset)
# plus leave-one-dataset-out FE meta-analysis and per-dataset heterogeneity (I2, Q).
#
# Input : data/neuron_counts_469donors.csv (Nicole's raw per-donor neuronal supertype
#         counts, 109 supertypes; built by the prep step from 7_cohorts_metadata_names.csv)
# Output: results/composition_pooling_sensitivity.csv (long: supertype x method)
#         results/composition_per_dataset_estimates.csv
#         results/composition_lodo.csv
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(crumblr)
  library(variancePartition); library(metafor)
})
args <- commandArgs(trailingOnly = FALSE)
fp   <- sub("^--file=", "", args[grep("^--file=", args)])
HERE <- if (length(fp)) normalizePath(file.path(dirname(fp), "..")) else normalizePath(getwd())
ROOT <- normalizePath(file.path(HERE, "..", ".."))          # repo root
d   <- read_csv(file.path(HERE, "data/neuron_counts_469donors.csv"), show_col_types = FALSE)
sup <- read_csv(file.path(HERE, "data/neuronal_supertypes_109.csv"), show_col_types = FALSE)$supertype
d[sup][is.na(d[sup])] <- 0
d <- d |> mutate(dx = factor(dx, levels = c("Control", "SCZ")), sex = factor(sex),
                 dataset = factor(dataset))
cat(sprintf("%d donors, %d datasets, %d supertypes\n", nrow(d), nlevels(d$dataset), length(sup)))

extract <- function(fit, coef = "dxSCZ") {
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  tibble(supertype = rownames(tt), beta = tt$logFC, se = tt$logFC / tt$t, p = tt$P.Value)
}

# ---------------------------------------------------------------- per dataset
per <- list()
for (ds in levels(d$dataset)) {
  sub <- d[d$dataset == ds, ]
  X <- as.matrix(sub[, sup]); rownames(X) <- sub$donor
  X <- X[, colSums(X) > 0, drop = FALSE]                  # supertypes absent from this dataset
  info <- data.frame(dx = sub$dx, age_z = as.numeric(scale(sub$age)), sex = droplevels(sub$sex),
                     pmi_z = as.numeric(scale(sub$PMI)), row.names = sub$donor)
  form <- if (all(is.na(info$pmi_z))) ~ dx + age_z + sex else ~ dx + age_z + sex + pmi_z
  cobj <- crumblr(X)
  fit  <- eBayes(dream(cobj, form, info, quiet = TRUE))
  per[[ds]] <- extract(fit) |> mutate(dataset = ds, n_donors = nrow(sub))
  cat(sprintf("  %-9s n = %3d, %d supertypes, PMI %s\n", ds, nrow(sub), ncol(X),
              if (length(all.vars(form)) == 4) "in" else "omitted"))
}
per <- bind_rows(per)
write_csv(per, file.path(HERE, "results/composition_per_dataset_estimates.csv"))

# ---------------------------------------------------------------- meta-analyses
meta_one <- function(df, method) {
  if (nrow(df) < 2) return(tibble(beta = NA, se = NA, p = NA, k = nrow(df), I2 = NA, Q_p = NA))
  m <- tryCatch(rma(yi = beta, sei = se, data = df, method = method), error = function(e) NULL)
  if (is.null(m)) return(tibble(beta = NA, se = NA, p = NA, k = nrow(df), I2 = NA, Q_p = NA))
  tibble(beta = as.numeric(m$b), se = m$se, p = m$pval, k = m$k, I2 = m$I2, Q_p = m$QEp)
}
meta <- bind_rows(
  per |> group_by(supertype) |> group_modify(~ meta_one(.x, "FE"))   |> mutate(method = "FE meta"),
  per |> group_by(supertype) |> group_modify(~ meta_one(.x, "REML")) |> mutate(method = "RE meta")) |>
  ungroup()

# leave-one-dataset-out, FE
lodo <- bind_rows(lapply(levels(d$dataset), function(drop)
  per |> filter(dataset != drop) |> group_by(supertype) |>
    group_modify(~ meta_one(.x, "FE")) |> ungroup() |> mutate(dropped = drop)))
write_csv(lodo, file.path(HERE, "results/composition_lodo.csv"))

# ---------------------------------------------------------------- mega models
X <- as.matrix(d[, sup]); rownames(X) <- d$donor
pmi <- d$PMI; pmi[!is.na(pmi) & pmi == 0] <- NA             # one MSSM 2 donor coded 0 h = missing
pmi_z <- as.numeric(scale(pmi)); pmi_z[is.na(pmi_z)] <- 0     # Multiome (no PMI) and that donor -> mean (0)
# PMI is in hours in every dataset (ranges 2-72 h); dataset means differ (Fröhlich ~34 h, MSSM 2 ~15 h)
# and are absorbed by the dataset term.
info <- data.frame(dx = d$dx, age_z = as.numeric(scale(d$age)), sex = d$sex, pmi_z = pmi_z,
                   dataset = d$dataset, row.names = d$donor)
cobj <- crumblr(X)
mega <- bind_rows(
  extract(eBayes(dream(cobj, ~ dx + age_z + sex + pmi_z + dataset,          info, quiet = TRUE))) |> mutate(method = "mega: dataset fixed"),
  extract(eBayes(dream(cobj, ~ dx + age_z + sex + pmi_z + (1 | dataset),    info, quiet = TRUE))) |> mutate(method = "mega: dataset random intercept"),
  extract(eBayes(dream(cobj, ~ dx + age_z + sex + pmi_z + (1 + dx | dataset), info, quiet = TRUE))) |> mutate(method = "mega: random intercept + SCZ slope"))

# ------------------------------------------------- published FE meta (Nicole's result)
# Optional cross-check against the published meta-analysis, read through the canonical
# seam in shared/snrnaseq_de/ (git-ignored; see its README). When the file is there the
# rebuilt FE estimates are compared against it and it supplies the "published FE meta"
# row; otherwise the rebuilt FE meta-analysis is used for that row.
pub_path <- file.path(ROOT, "shared/snrnaseq_de/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv")
if (file.exists(pub_path)) {
  pub <- read_csv(pub_path, show_col_types = FALSE) |>
    transmute(supertype = CellType, beta = estimate, se = se, p = pval, method = "published FE meta")
  chk <- inner_join(pub |> select(supertype, b_pub = beta), meta |> filter(method == "FE meta") |> select(supertype, b_new = beta), by = "supertype")
  cat(sprintf("published FE file found: %d supertypes, Spearman rho vs rebuilt FE = %.4f, max |diff| = %.2e\n",
              nrow(chk), cor(chk$b_pub, chk$b_new, method = "spearman"), max(abs(chk$b_pub - chk$b_new))))
} else {
  cat("published FE file not found; using the rebuilt fixed-effect meta-analysis for that row\n")
  pub <- meta |> filter(method == "FE meta") |> select(supertype, beta, se, p) |> mutate(method = "published FE meta")
}

res <- bind_rows(pub, meta |> select(supertype, beta, se, p, k, I2, Q_p, method), mega) |>
  group_by(method) |> mutate(fdr = p.adjust(p, method = "BH")) |> ungroup()
write_csv(res, file.path(HERE, "results/composition_pooling_sensitivity.csv"))

# ---------------------------------------------------------------- quick report
w <- res |> select(supertype, method, beta) |> pivot_wider(names_from = method, values_from = beta)
cat("\nSpearman rho of beta vs published FE:\n")
for (m in setdiff(names(w), c("supertype", "published FE meta")))
  cat(sprintf("  %-38s rho = %.3f\n", m, cor(w$`published FE meta`, w[[m]], method = "spearman", use = "complete.obs")))
cat("\nFDR < 0.10 counts by method:\n"); print(res |> group_by(method) |> summarise(n = sum(fdr < 0.10, na.rm = TRUE)))
hl <- c("Sst_2","Sst_3","Sst_20","Sst_22","Sst_25","L6b_1","L6b_4")
cat("\nheadline supertypes:\n")
print(res |> filter(supertype %in% hl) |> select(supertype, method, beta, p, fdr, I2) |>
        mutate(across(c(beta, I2), ~ round(.x, 2)), across(c(p, fdr), ~ signif(.x, 2))) |>
        arrange(factor(supertype, levels = hl), method), n = 60)
