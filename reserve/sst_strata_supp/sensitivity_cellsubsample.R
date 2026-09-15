#!/usr/bin/env Rscript
# Sensitivity | DE and GSEA on CELL-MATCHED pseudobulks.
#
# Input comes from subsample_cells.py, which subsamples control donors' nuclei
# down to the case distribution within each cohort and stratum and rebuilds the
# pseudobulks. This is the definitive form of the cell-count control: the
# imbalance that defines the groups is removed by construction rather than by
# modelling, so there is no confounder-versus-mediator ambiguity and no donor is
# discarded. Everything downstream is the main pipeline, unchanged.
#
# Reported across independent random draws, since each draw is one realisation.
source("transcriptomic/scripts/sst_strata/_common.R")
suppressPackageStartupMessages({ library(arrow); library(edgeR); library(limma); library(parallel) })
N_CORES <- max(1, detectCores() - 2)
SUB  <- file.path(P$supp, "cellsub")
reps <- sort(as.integer(str_remove(basename(Sys.glob(file.path(SUB, "rep*"))), "rep")))
covars <- load_covars(); gs <- msigdb_sets()
say("replicates found: %s", paste(reps, collapse = ", "))

load_rep <- function(rep, co) {
  f <- file.path(SUB, paste0("rep", rep), paste0(co, "_donor_stratum_counts.parquet"))
  pb <- as.data.frame(read_parquet(f)); rn <- pb$gene; pb$gene <- NULL
  M <- vapply(pb, as.numeric, numeric(nrow(pb))); rownames(M) <- rn
  m <- read_csv(file.path(SUB, paste0("rep", rep), paste0(co, "_donor_stratum_meta.csv")),
                col_types = cols(donor = col_character(), .default = col_guess()))
  list(counts = M, meta = m)
}
de_one <- function(rep, co, st) {
  d <- load_rep(rep, co)
  m <- d$meta |> filter(stratum == st) |>
    left_join(covars |> filter(cohort == co) |> select(-cohort), by = "donor") |>
    mutate(key = paste(donor, stratum, sep = "|")) |> prep_model_frame()
  if (nrow(m) < 8 || n_distinct(m$diagnosis) < 2) return(NULL)
  M <- d$counts[, m$key, drop = FALSE]
  y <- calcNormFactors(DGEList(M[rowSums(M >= 1) >= 0.8 * nrow(m), ]))
  des <- model.matrix(if (attr(m, "use_pmi")) ~ dx + age_s + sex + pmi_s else ~ dx + age_s + sex,
                      data = m)
  stopifnot(nrow(des) == ncol(y))
  tt <- topTable(eBayes(lmFit(voom(y, des), des)), coef = "dxSCZ", number = Inf, sort.by = "none")
  tibble(cohort = co, gene = rownames(tt), logFC = tt$logFC, SE = tt$logFC / tt$t)
}

res <- map_dfr(reps, function(rep) {
  map_dfr(STRATA_LEVELS, function(st) {
    pc <- map_dfr(export_cohorts(), ~ de_one(rep, .x, st))
    genes <- pc |> count(gene) |> filter(n >= 5) |> pull(gene)
    mt <- bind_rows(setNames(mclapply(genes, function(g) {
      x <- pc[pc$gene == g, ]; fit_meta(x$logFC, x$SE) }, mc.cores = N_CORES), genes), .id = "gene")
    suppressMessages(run_fgsea(mt$gene, mt$zval, st, gs$sets)) |>
      mutate(stratum = st, rep = rep)
  })
})
write_csv(res |> select(-leadingEdge), file.path(P$supp, "sensitivity_cellsub_gsea.csv"))

say("=== gene sets at FDR < 0.10, per replicate (cell-matched) ===")
print(res |> group_by(stratum, rep) |> summarise(n = sum(padj < 0.10), .groups = "drop") |>
        pivot_wider(names_from = rep, values_from = n, names_prefix = "rep") |>
        arrange(match(stratum, STRATA_LEVELS)))

say("=== panel-d blocks: NES (FDR) per replicate ===")
print(res |> inner_join(BLOCKS, by = "pathway") |>
        mutate(cell = sprintf("%.2f (%.1g)", NES, padj)) |>
        select(block, pathway, stratum, rep, cell) |>
        pivot_wider(names_from = rep, values_from = cell, names_prefix = "rep") |>
        arrange(match(block, BLOCK_KEYS), match(stratum, STRATA_LEVELS)), n = 40, width = 250)

say("=== summary across replicates for the two graded blocks ===")
print(res |> inner_join(BLOCKS, by = "pathway") |> filter(block %in% c("translation", "oxphos")) |>
        group_by(block, stratum) |>
        summarise(median_NES = round(median(NES), 2),
                  n_sig_of_total = sprintf("%d/%d", sum(padj < 0.10), n()), .groups = "drop") |>
        arrange(match(block, BLOCK_KEYS), match(stratum, STRATA_LEVELS)))
