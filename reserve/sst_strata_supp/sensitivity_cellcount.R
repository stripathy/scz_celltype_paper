#!/usr/bin/env Rscript
# Sensitivity | Is the graded dysregulation just a consequence of the depleted
# group containing fewer SCZ nuclei?
#
# THE CONCERN. The groups are DEFINED by SCZ having proportionally fewer cells,
# so within the depleted group SCZ donors contribute ~30% fewer nuclei to their
# pseudobulk than controls (18% intermediate, 4% and n.s. non-depleted). The
# imbalance is therefore graded exactly like the result, and the within-donor
# interaction design does NOT remove it: it cancels donor-level factors, but the
# cell-count difference is a donor-by-group factor.
#
# TEST 1  Refit the differential expression with log(nuclei per donor) as an
#         additional covariate, then repeat GSEA. If the translation and
#         oxidative-phosphorylation signals survive adjustment for the very
#         quantity that defines the groups, the count imbalance is not what
#         produces them. Caveat: nuclei count is collinear with diagnosis by
#         construction, so this adjustment is conservative and will absorb some
#         genuine signal.
# TEST 2  Artifact signature. A library-size or normalisation artifact acts on
#         genes according to how highly they are expressed. Ribosomal and
#         electron-transport genes are among the most highly expressed in any
#         neuron, so if the effect is technical it should appear as a monotonic
#         relationship between a gene's expression level and its diagnosis
#         effect, and the modules should be unremarkable once that is accounted
#         for.
source("transcriptomic/scripts/fig5/_common.R")
suppressPackageStartupMessages({ library(edgeR); library(limma); library(parallel) })
set.seed(42)
N_CORES <- max(1, detectCores() - 2)
OUT <- P$supp

covars <- load_covars()
cohorts <- donor_cohorts()

# per-cohort DE within one group, optionally adjusting for log nuclei count
de_group <- function(co, st, adjust_n) {
  d <- load_donor_stratum(co)
  m <- d$meta |> filter(stratum == st) |>
    left_join(covars |> filter(cohort == co) |> select(-cohort), by = "donor") |>
    mutate(key = paste(donor, stratum, sep = "|"), log_n = log(n_cells)) |>
    prep_model_frame()
  if (nrow(m) < 8 || n_distinct(m$diagnosis) < 2) return(NULL)
  use_pmi <- attr(m, "use_pmi")
  M <- d$counts[, m$key, drop = FALSE]
  y <- calcNormFactors(DGEList(M[rowSums(M >= 1) >= 0.8 * nrow(m), ]))
  rhs <- c("dx", "age_s", "sex", if (use_pmi) "pmi_s", if (adjust_n) "log_n")
  des <- model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), data = m)
  stopifnot(nrow(des) == ncol(y))
  fit <- eBayes(lmFit(voom(y, des), des))
  tt <- topTable(fit, coef = "dxSCZ", number = Inf, sort.by = "none")
  tibble(cohort = co, stratum = st, gene = rownames(tt), logFC = tt$logFC,
         SE = tt$logFC / tt$t)
}
# min_k is relaxed for the matched analysis: matching drops whole cohorts that
# have no overlapping count range, so k >= 5 can leave no testable genes.
meta_of <- function(pc, min_k = 5) {
  genes <- pc |> count(gene) |> filter(n >= min_k) |> pull(gene)
  if (!length(genes)) return(tibble(gene = character(), zval = numeric()))
  res <- mclapply(genes, function(gn) { d <- pc[pc$gene == gn, ]; fit_meta(d$logFC, d$SE) },
                  mc.cores = N_CORES)
  bind_rows(setNames(res, genes), .id = "gene")
}

gs <- msigdb_sets()
say("=== TEST 1: adjusting for log(nuclei per donor) ===")
res <- map_dfr(c(FALSE, TRUE), function(adj) {
  map_dfr(STRATA_LEVELS, function(st) {
    pc <- map_dfr(cohorts, ~ de_group(.x, st, adj))
    mt <- meta_of(pc)
    run_fgsea(mt$gene, mt$zval, st, gs$sets) |>
      mutate(stratum = st, adjusted = adj)
  })
})
write_csv(res |> select(-leadingEdge), file.path(OUT, "sensitivity_cellcount_gsea.csv"))

say("burden at FDR < 0.10, unadjusted vs adjusted for nuclei count:")
print(res |> group_by(stratum, adjusted) |> summarise(n10 = sum(padj < 0.10), .groups = "drop") |>
        pivot_wider(names_from = adjusted, values_from = n10,
                    names_prefix = "adj_") |> arrange(match(stratum, STRATA_LEVELS)))

say("the panel-d blocks, unadjusted vs adjusted:")
print(res |> inner_join(BLOCKS, by = "pathway") |>
        mutate(cell = sprintf("%.2f (%.2g)", NES, padj)) |>
        select(block, pathway, stratum, adjusted, cell) |>
        pivot_wider(names_from = c(stratum, adjusted), values_from = cell) |>
        arrange(match(block, BLOCK_KEYS)), n = 14, width = 230)

# TEST 3  Count-matched donors. Regression adjustment (Test 1) is ambiguous here
#         because the cell loss could be a mediator of the disease effect rather
#         than a confounder, and adjusting for a mediator removes real signal.
#         Matching sidesteps that: within each cohort, pair each SCZ donor with
#         the control donor closest in log nuclei count (caliper 0.25 SD), so the
#         two arms have the same count distribution by design, then rerun.
say("=== TEST 3: count-matched donors ===")
match_group <- function(co, st) {
  d <- load_donor_stratum(co)
  m <- d$meta |> filter(stratum == st) |> mutate(log_n = log(n_cells))
  if (n_distinct(m$diagnosis) < 2) return(NULL)
  cal <- 0.25 * sd(m$log_n)
  scz <- m |> filter(diagnosis == "SCZ"); ctl <- m |> filter(diagnosis == "Control")
  used <- rep(FALSE, nrow(ctl)); keep_s <- c(); keep_c <- c()
  for (i in order(scz$log_n)) {
    dd <- abs(ctl$log_n - scz$log_n[i]); dd[used] <- Inf
    j <- which.min(dd)
    if (is.finite(dd[j]) && dd[j] <= cal) { used[j] <- TRUE; keep_s <- c(keep_s, i); keep_c <- c(keep_c, j) }
  }
  if (length(keep_s) < 4) return(NULL)
  bind_rows(scz[keep_s, ], ctl[keep_c, ])
}
de_matched <- function(co, st) {
  mm <- match_group(co, st); if (is.null(mm)) return(NULL)
  d <- load_donor_stratum(co)
  m <- mm |> left_join(covars |> filter(cohort == co) |> select(-cohort), by = "donor") |>
    mutate(key = paste(donor, stratum, sep = "|")) |> prep_model_frame()
  use_pmi <- attr(m, "use_pmi")
  M <- d$counts[, m$key, drop = FALSE]
  y <- calcNormFactors(DGEList(M[rowSums(M >= 1) >= 0.8 * nrow(m), ]))
  des <- model.matrix(if (use_pmi) ~ dx + age_s + sex + pmi_s else ~ dx + age_s + sex, data = m)
  fit <- eBayes(lmFit(voom(y, des), des))
  tt <- topTable(fit, coef = "dxSCZ", number = Inf, sort.by = "none")
  tibble(cohort = co, gene = rownames(tt), logFC = tt$logFC, SE = tt$logFC / tt$t,
         n_pairs = nrow(m) / 2)
}
mres <- map_dfr(STRATA_LEVELS, function(st) {
  pc <- map_dfr(cohorts, ~ de_matched(.x, st))
  say("  %s: %d matched pairs across %d cohorts", st,
      sum(pc |> distinct(cohort, n_pairs) |> pull(n_pairs)), n_distinct(pc$cohort))
  mt <- meta_of(pc, min_k = 3)
  say("    %d genes with k >= 3", nrow(mt))
  run_fgsea(mt$gene, mt$zval, st, gs$sets) |> mutate(stratum = st)
})
write_csv(mres |> select(-leadingEdge), file.path(OUT, "sensitivity_cellcount_matched_gsea.csv"))
say("matched-donor burden at FDR < 0.10:")
print(mres |> group_by(stratum) |> summarise(n10 = sum(padj < 0.10), .groups = "drop") |>
        arrange(match(stratum, STRATA_LEVELS)))
say("matched-donor results for the panel-d blocks:")
print(mres |> inner_join(BLOCKS, by = "pathway") |>
        mutate(cell = sprintf("%.2f (%.2g)", NES, padj)) |>
        select(block, pathway, stratum, cell) |>
        pivot_wider(names_from = stratum, values_from = cell) |>
        arrange(match(block, BLOCK_KEYS)), n = 14, width = 200)

# TEST 4  Depth-matched by downsampling (Shreejoy's proposal). Rather than
#         adjusting statistically or discarding donors, equalise the thing that
#         differs: binomially thin each CONTROL pseudobulk so that the library-size
#         distribution of controls matches that of cases, within each cohort and
#         group, then rerun unchanged. Targets are assigned by quantile matching
#         so the whole distribution matches, not just the mean.
#
#         What this does and does not simulate. SCZ pseudobulks in the depleted
#         group carry 36% fewer counts and come from 30% fewer nuclei, and library
#         size and nuclei count correlate at r = 0.80. Thinning reproduces the
#         sequencing-depth half of that exactly. It cannot reproduce the loss of
#         cellular sampling, since a pseudobulk cannot be resampled back into
#         cells. A definitive version needs the cell-level matrices.
say("=== TEST 4: controls downsampled to the case library-size distribution ===")
thin_cohort_group <- function(co, st, seed = 1) {
  set.seed(seed)
  d <- load_donor_stratum(co)
  m <- d$meta |> filter(stratum == st) |>
    left_join(covars |> filter(cohort == co) |> select(-cohort), by = "donor") |>
    mutate(key = paste(donor, stratum, sep = "|")) |> prep_model_frame()
  if (nrow(m) < 8 || n_distinct(m$diagnosis) < 2) return(NULL)
  M <- d$counts[, m$key, drop = FALSE]
  lib <- colSums(M)
  is_scz <- m$dx == "SCZ"
  if (sum(is_scz) < 3 || sum(!is_scz) < 3) return(NULL)
  # quantile-match control libraries onto the case distribution
  ctl <- which(!is_scz)
  pr <- (rank(lib[ctl], ties.method = "average")) / (length(ctl) + 1)
  target <- quantile(lib[is_scz], probs = pr, names = FALSE)
  for (i in seq_along(ctl)) {
    j <- ctl[i]; p <- min(1, target[i] / lib[j])
    if (p < 1) M[, j] <- rbinom(nrow(M), size = round(M[, j]), prob = p)
  }
  y <- calcNormFactors(DGEList(M[rowSums(M >= 1) >= 0.8 * nrow(m), ]))
  use_pmi <- attr(m, "use_pmi")
  des <- model.matrix(if (use_pmi) ~ dx + age_s + sex + pmi_s else ~ dx + age_s + sex, data = m)
  fit <- eBayes(lmFit(voom(y, des), des))
  tt <- topTable(fit, coef = "dxSCZ", number = Inf, sort.by = "none")
  attr_lib <- c(scz = median(lib[is_scz]), ctl_before = median(lib[ctl]),
                ctl_after = median(colSums(M)[ctl]))
  list(de = tibble(cohort = co, gene = rownames(tt), logFC = tt$logFC, SE = tt$logFC / tt$t),
       lib = attr_lib)
}
tres <- map_dfr(STRATA_LEVELS, function(st) {
  out <- map(cohorts, ~ thin_cohort_group(.x, st)) |> compact()
  libs <- do.call(rbind, map(out, "lib"))
  say("  %s: median library size, cases %s | controls before %s | controls after %s",
      st, format(round(median(libs[, "scz"])), big.mark = ","),
      format(round(median(libs[, "ctl_before"])), big.mark = ","),
      format(round(median(libs[, "ctl_after"])), big.mark = ","))
  mt <- meta_of(map_dfr(out, "de"))
  run_fgsea(mt$gene, mt$zval, st, gs$sets) |> mutate(stratum = st)
})
write_csv(tres |> select(-leadingEdge), file.path(OUT, "sensitivity_cellcount_thinned_gsea.csv"))
say("depth-matched burden at FDR < 0.10:")
print(tres |> group_by(stratum) |> summarise(n10 = sum(padj < 0.10), .groups = "drop") |>
        arrange(match(stratum, STRATA_LEVELS)))
say("depth-matched results for the panel-d blocks:")
print(tres |> inner_join(BLOCKS, by = "pathway") |>
        mutate(cell = sprintf("%.2f (%.2g)", NES, padj)) |>
        select(block, pathway, stratum, cell) |>
        pivot_wider(names_from = stratum, values_from = cell) |>
        arrange(match(block, BLOCK_KEYS)), n = 14, width = 200)
