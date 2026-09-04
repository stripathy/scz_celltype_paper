#!/usr/bin/env Rscript
# Stratum-level pseudobulk DE + meta-analysis from the 7-cohort export
# (transcriptomic/data/stratum_pseudobulks_export/), replicating the paper's
# snRNA-seq DE framework (Methods, "Cell type-specific differential gene
# expression"):
#   per cohort: donors with >= 10 cells of the group; genes with >= 1 count in
#   >= 80% of retained donors; TMM -> voom -> limma moderated t;
#   ~ Diagnosis + scale(age) + sex + scale(PMI)  (PMI omitted for Multiome);
#   meta: metafor::rma REML over log2FC/SE, k >= 5 cohorts; BH-FDR within
#   cell type (here: within stratum).
# Strata: depleted / intermediate / non_depleted (crumblr FDR < 0.20 definition)
# plus all_sst (pooled 16 supertypes).
# Gene space: symbols; HBCC/MSSM2 (Ensembl-native) use the export's GENCODE v44
# gene_symbol column (flagged in the MANIFEST); duplicate symbols are summed at
# the count level; unmapped genes dropped (counted below).
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(tibble); library(arrow); library(edgeR); library(limma)
  library(metafor); library(parallel)
})
EXP <- "transcriptomic/data/stratum_pseudobulks_export"
OUT <- "transcriptomic/results/sst_strata_gsea/pseudobulk"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
N_CORES <- max(1, detectCores() - 2)

STRATA <- list(
  depleted     = c("Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25"),
  intermediate = c("Sst_9", "Sst_11", "Sst_12", "Sst_13", "Sst_19", "Sst_23"),
  non_depleted = c("Sst_1", "Sst_4", "Sst_5", "Sst_7", "Sst_10"))
st_of <- setNames(rep(names(STRATA), lengths(STRATA)), unlist(STRATA, use.names = FALSE))
GROUPS <- c(names(STRATA), "all_sst")
MIN_CELLS <- 10

t0 <- Sys.time()
say <- function(...) cat(sprintf("[%6.1f min] ", as.numeric(difftime(Sys.time(), t0, units = "mins"))), sprintf(...), "\n", sep = "")

cohorts <- str_remove(basename(Sys.glob(file.path(EXP, "*_groups.csv"))), "_groups\\.csv$")
say("cohorts: %s", paste(cohorts, collapse = ", "))
stopifnot(length(cohorts) == 7)

# ---- per-cohort: aggregate to strata, then DE ---------------------------------
run_cohort <- function(co) {
  gr <- read_csv(file.path(EXP, paste0(co, "_groups.csv")),
                 col_types = cols(donor = col_character(), .default = col_guess())) |>
    filter(supertype %in% names(st_of), n_cells > 0)
  need <- gr$key
  f_pq <- file.path(EXP, paste0(co, "_pseudobulk_counts.parquet"))
  idcols <- intersect(c("gene", "gene_symbol"), names(open_dataset(f_pq)))
  pb <- as.data.frame(read_parquet(f_pq, col_select = all_of(c(idcols, need))))
  sym <- if ("gene_symbol" %in% names(pb)) coalesce(pb$gene_symbol, NA_character_) else pb$gene
  if (is.null(sym)) sym <- pb$gene
  # fall back to native id when it already looks like a symbol and mapping is NA
  looks_sym <- !grepl("^ENSG", pb$gene)
  sym[is.na(sym) & looks_sym] <- pb$gene[is.na(sym) & looks_sym]
  drop_n <- sum(is.na(sym))
  keep <- !is.na(sym)
  M <- vapply(pb[keep, need, drop = FALSE], as.numeric, numeric(sum(keep)))
  rownames(M) <- NULL
  M <- rowsum(M, group = sym[keep])              # sum duplicate symbols
  say("%s: %d genes (symbols; %d unmapped dropped), %d Sst groups",
      co, nrow(M), drop_n, ncol(M))

  # aggregate supertypes -> strata per donor
  gr <- gr |> mutate(stratum = st_of[supertype])
  agg <- function(strata_set, label) {
    gsub <- gr |> filter(stratum %in% strata_set)
    donors <- gsub |> group_by(donor) |>
      summarise(n_cells = sum(n_cells),
                diagnosis = first(diagnosis), sex = first(sex),
                age = suppressWarnings(as.numeric(first(age))),
                pmi = suppressWarnings(as.numeric(first(pmi))), .groups = "drop") |>
      filter(n_cells >= MIN_CELLS)
    if (nrow(donors) < 8 || n_distinct(donors$diagnosis) < 2) return(NULL)
    cols <- sapply(donors$donor, function(d)
      rowSums(M[, gsub$key[gsub$donor == d], drop = FALSE]))
    list(counts = cols, meta = donors, label = label)
  }
  sets <- c(map(names(STRATA), ~ agg(.x, .x)), list(agg(names(STRATA), "all_sst")))
  sets <- compact(sets)

  map_dfr(sets, function(sx) {
    m <- sx$meta |> mutate(dx = factor(diagnosis, c("Control", "SCZ")),
                           sex = factor(sex),
                           age_s = scale(ifelse(is.na(age), median(age, na.rm = TRUE), age))[, 1])
    use_pmi <- mean(is.na(m$pmi)) < 0.5
    if (use_pmi) m$pmi_s <- scale(ifelse(is.na(m$pmi), median(m$pmi, na.rm = TRUE), m$pmi))[, 1]
    y <- DGEList(sx$counts[, m$donor, drop = FALSE])
    keep_g <- rowSums(y$counts >= 1) >= 0.8 * nrow(m)    # >=1 count in >=80% donors
    y <- calcNormFactors(y[keep_g, ])
    f <- if (use_pmi) ~ dx + age_s + sex + pmi_s else ~ dx + age_s + sex
    des <- model.matrix(f, data = m)
    fit <- eBayes(lmFit(voom(y, des), des))
    tt <- topTable(fit, coef = "dxSCZ", number = Inf, sort.by = "none")
    tibble(cohort = co, stratum = sx$label, gene = rownames(tt),
           logFC = tt$logFC, SE = tt$logFC / tt$t, t = tt$t, P.Value = tt$P.Value,
           n_donors = nrow(m), pmi_in_model = use_pmi)
  })
}

percohort <- map_dfr(cohorts, run_cohort)
write_csv(percohort, file.path(OUT, "stratum_percohort_de.csv"))
say("per-cohort DE done: %d rows", nrow(percohort))
print(percohort |> distinct(cohort, stratum, n_donors, pmi_in_model) |>
        pivot_wider(names_from = stratum, values_from = n_donors), n = 30)

# ---- meta-analysis (metafor REML, k >= 5) --------------------------------------
say("meta-analysis ...")
meta_one <- function(d) {
  fit <- tryCatch(rma(yi = d$logFC, sei = d$SE, method = "REML",
                      control = list(maxiter = 1000)),
                  error = function(e) tryCatch(rma(yi = d$logFC, sei = d$SE, method = "DL"),
                                               error = function(e2) NULL))
  if (is.null(fit)) return(NULL)
  tibble(estimate = fit$beta[1], se = fit$se, zval = fit$zval, pval = fit$pval,
         ci.lb = fit$ci.lb, ci.ub = fit$ci.ub, k = fit$k, tau2 = fit$tau2, I2 = fit$I2)
}
meta_all <- map_dfr(GROUPS, function(st) {
  dd <- percohort |> filter(stratum == st)
  genes <- dd |> count(gene) |> filter(n >= 5) |> pull(gene)
  say("  %s: %d genes with k >= 5", st, length(genes))
  res <- mclapply(genes, function(gn) meta_one(dd[dd$gene == gn, ]),
                  mc.cores = N_CORES)
  bind_rows(setNames(res, genes), .id = "gene") |> mutate(stratum = st)
})
meta_all <- meta_all |> group_by(stratum) |>
  mutate(padj = p.adjust(pval, "BH")) |> ungroup()
write_csv(meta_all, file.path(OUT, "stratum_meta_de.csv"))
say("meta done: %s",
    paste(capture.output(print(count(meta_all, stratum)))[-(1:3)], collapse = "; "))
for (st in GROUPS) {
  s <- meta_all |> filter(stratum == st)
  say("  %s: %d genes, %d at FDR<0.05 (%d down)", st, nrow(s),
      sum(s$padj < 0.05), sum(s$padj < 0.05 & s$estimate < 0))
}
say("SST/VGF/RPL36/NDUFS8 check (zval by stratum):")
print(meta_all |> filter(gene %in% c("SST", "VGF", "RPL36", "NDUFS8", "CALB1")) |>
        select(gene, stratum, zval) |> mutate(zval = round(zval, 2)) |>
        pivot_wider(names_from = stratum, values_from = zval))
