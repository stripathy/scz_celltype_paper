# Subclass-level companion to 01: the same per-dataset crumblr model on supertype counts
# aggregated to subclass (18 neuronal subclasses), fixed-effect meta-analysis and
# leave-one-dataset-out. Gives the composite Sst estimate quoted next to the supertype
# leave-one-out results in the pooling-sensitivity figure.
# Output: results/composition_subclass_estimates.csv, results/composition_subclass_lodo.csv
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(crumblr); library(variancePartition); library(metafor)
})
args <- commandArgs(trailingOnly = FALSE)
fp   <- sub("^--file=", "", args[grep("^--file=", args)])
HERE <- if (length(fp)) normalizePath(file.path(dirname(fp), "..")) else normalizePath(getwd())
ROOT <- normalizePath(file.path(HERE, "..", ".."))          # repo root
d   <- read_csv(file.path(HERE, "data/neuron_counts_469donors.csv"), show_col_types = FALSE)
sup <- read_csv(file.path(HERE, "data/neuronal_supertypes_109.csv"), show_col_types = FALSE)$supertype
d[sup][is.na(d[sup])] <- 0
d <- d |> mutate(dx = factor(dx, levels = c("Control", "SCZ")), sex = factor(sex), dataset = factor(dataset))
subclass <- sub("(_\\d+)+$", "", sup)
S <- sapply(unique(subclass), function(s) rowSums(as.matrix(d[, sup[subclass == s], drop = FALSE])))
rownames(S) <- d$donor
cat(sprintf("%d donors, %d subclasses\n", nrow(S), ncol(S)))

extract <- function(fit, coef = "dxSCZ") {
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  tibble(subclass = rownames(tt), beta = tt$logFC, se = tt$logFC / tt$t, p = tt$P.Value)
}
per <- list()
for (ds in levels(d$dataset)) {
  i <- d$dataset == ds
  X <- S[i, colSums(S[i, , drop = FALSE]) > 0, drop = FALSE]
  info <- data.frame(dx = d$dx[i], age_z = as.numeric(scale(d$age[i])), sex = droplevels(d$sex[i]),
                     pmi_z = as.numeric(scale(d$PMI[i])), row.names = d$donor[i])
  form <- if (all(is.na(info$pmi_z))) ~ dx + age_z + sex else ~ dx + age_z + sex + pmi_z
  per[[ds]] <- extract(eBayes(dream(crumblr(X), form, info, quiet = TRUE))) |>
    mutate(dataset = ds, n_donors = sum(i))
}
per <- bind_rows(per)
fe <- function(df) { m <- rma(yi = beta, sei = se, data = df, method = "FE")
  tibble(beta = as.numeric(m$b), se = m$se, p = m$pval, k = m$k, I2 = m$I2, Q_p = m$QEp) }
full <- per |> group_by(subclass) |> group_modify(~ fe(.x)) |> ungroup() |>
  mutate(fdr = p.adjust(p, "BH"), dropped = "none")
lodo <- bind_rows(lapply(levels(d$dataset), function(dr)
  per |> filter(dataset != dr) |> group_by(subclass) |> group_modify(~ fe(.x)) |> ungroup() |>
    mutate(fdr = p.adjust(p, "BH"), dropped = dr)))
write_csv(per, file.path(HERE, "results/composition_subclass_estimates.csv"))
write_csv(bind_rows(full, lodo), file.path(HERE, "results/composition_subclass_lodo.csv"))
cat("FE meta, subclass level (FDR over 18):\n"); print(as.data.frame(full |> arrange(p) |> select(-dropped)), digits = 3)
cat("\nSst leave-one-dataset-out:\n"); print(as.data.frame(lodo |> filter(subclass == "Sst") |> select(dropped, beta, se, p, fdr)), digits = 3)
