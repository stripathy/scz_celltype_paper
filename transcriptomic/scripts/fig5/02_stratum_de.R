#!/usr/bin/env Rscript
# Step 2 | Stratum-level pseudobulk DE + across-cohort meta-analysis.
#
# Replicates the paper's snRNA-seq DE framework (Methods, "Cell type-specific
# differential gene expression") one aggregation level up -- strata instead of
# supertypes -- so Figure 5 is directly comparable with Figure 2:
#   per cohort: donors with >= MIN_CELLS nuclei of the group; genes with >= 1
#   count in >= 80% of retained donors; TMM -> voom -> limma moderated t;
#   ~ dx + scale(age) + sex + scale(PMI)   (PMI omitted where >50% missing)
#   meta: metafor REML over log2FC/SE, k >= 5 cohorts; BH-FDR within stratum.
#
# Gene space: symbols. HBCC/MSSM2 are Ensembl-native and use the export's
# GENCODE v44 gene_symbol column (flagged in the export MANIFEST); duplicate
# symbols are summed at the count level and unmapped genes dropped.
#
# Runtime ~10 min on 8 cores. Outputs stratum_percohort_de.csv (large, kept for
# leave-one-cohort-out) and stratum_meta_de.csv (the canonical gene-level table).
source("transcriptomic/scripts/fig5/_common.R")
suppressPackageStartupMessages({
  library(arrow); library(edgeR); library(limma); library(parallel)
})
N_CORES <- max(1, detectCores() - 2)

st_of  <- stratum_of()
GROUPS <- c(STRATA_LEVELS, "all_sst")
cohorts <- export_cohorts()
say("cohorts: %s", paste(cohorts, collapse = ", "))
stopifnot(length(cohorts) == 7)

run_cohort <- function(co) {
  gr <- read_csv(file.path(P$exp, paste0(co, "_groups.csv")),
                 col_types = cols(donor = col_character(), .default = col_guess())) |>
    filter(supertype %in% names(st_of), n_cells > 0) |>
    mutate(stratum = st_of[supertype])
  f_pq <- file.path(P$exp, paste0(co, "_pseudobulk_counts.parquet"))
  idcols <- intersect(c("gene", "gene_symbol"), names(open_dataset(f_pq)))
  pb <- as.data.frame(read_parquet(f_pq, col_select = all_of(c(idcols, gr$key))))

  # resolve to symbols: prefer the mapped column, fall back to a native id that
  # already looks like a symbol
  sym <- if ("gene_symbol" %in% names(pb)) pb$gene_symbol else pb$gene
  looks_sym <- !grepl("^ENSG", pb$gene)
  sym[is.na(sym) & looks_sym] <- pb$gene[is.na(sym) & looks_sym]
  keep <- !is.na(sym)
  M <- vapply(pb[keep, gr$key, drop = FALSE], as.numeric, numeric(sum(keep)))
  rownames(M) <- NULL
  M <- rowsum(M, group = sym[keep])
  say("%s: %d genes (%d unmapped dropped), %d Sst groups", co, nrow(M), sum(!keep), ncol(M))

  # pool a stratum's supertypes into one pseudobulk per donor
  agg <- function(strata_set, label) {
    gsub <- gr |> filter(stratum %in% strata_set)
    donors <- gsub |> group_by(donor) |>
      summarise(n_cells = sum(n_cells), diagnosis = first(diagnosis), sex = first(sex),
                age = suppressWarnings(as.numeric(first(age))),
                pmi = suppressWarnings(as.numeric(first(pmi))), .groups = "drop") |>
      filter(n_cells >= MIN_CELLS)
    if (nrow(donors) < 8 || n_distinct(donors$diagnosis) < 2) return(NULL)
    cols <- sapply(donors$donor, function(d)
      rowSums(M[, gsub$key[gsub$donor == d], drop = FALSE]))
    list(counts = cols, meta = donors, label = label)
  }
  sets <- compact(c(map(STRATA_LEVELS, ~ agg(.x, .x)), list(agg(STRATA_LEVELS, "all_sst"))))

  map_dfr(sets, function(sx) {
    m <- prep_model_frame(sx$meta)
    use_pmi <- attr(m, "use_pmi")
    y <- DGEList(sx$counts[, m$donor, drop = FALSE])
    y <- calcNormFactors(y[rowSums(y$counts >= 1) >= 0.8 * nrow(m), ])
    des <- model.matrix(if (use_pmi) ~ dx + age_s + sex + pmi_s else ~ dx + age_s + sex,
                        data = m)
    stopifnot(nrow(des) == ncol(y))          # guard: NA covariates drop design rows
    fit <- eBayes(lmFit(voom(y, des), des))
    tt <- topTable(fit, coef = "dxSCZ", number = Inf, sort.by = "none")
    tibble(cohort = co, stratum = sx$label, gene = rownames(tt), logFC = tt$logFC,
           SE = tt$logFC / tt$t, t = tt$t, P.Value = tt$P.Value,
           n_donors = nrow(m), pmi_in_model = use_pmi)
  })
}

percohort <- map_dfr(cohorts, run_cohort)
write_csv(percohort, file.path(P$pb, "stratum_percohort_de.csv"))
say("per-cohort DE: %d rows", nrow(percohort))
print(percohort |> distinct(cohort, stratum, n_donors, pmi_in_model) |>
        pivot_wider(names_from = stratum, values_from = n_donors), n = 30)

say("meta-analysis (REML, k >= 5) ...")
meta_all <- map_dfr(GROUPS, function(st) {
  dd <- percohort |> filter(stratum == st)
  genes <- dd |> count(gene) |> filter(n >= 5) |> pull(gene)
  say("  %s: %d genes with k >= 5", st, length(genes))
  res <- mclapply(genes, function(gn) {
    d <- dd[dd$gene == gn, ]; fit_meta(d$logFC, d$SE, extra = TRUE)
  }, mc.cores = N_CORES)
  bind_rows(setNames(res, genes), .id = "gene") |> mutate(stratum = st)
}) |> group_by(stratum) |> mutate(padj = p.adjust(pval, "BH")) |> ungroup()
write_csv(meta_all, file.path(P$pb, "stratum_meta_de.csv"))

for (st in GROUPS) {
  s <- meta_all |> filter(stratum == st)
  say("  %s: %d genes | FDR<0.05 %d | FDR<0.10 %d", st, nrow(s),
      sum(s$padj < 0.05), sum(s$padj < 0.10))
}
print(meta_all |> filter(gene %in% c("SST", "VGF", "RPL36", "NDUFS8", "COX5B")) |>
        transmute(gene, stratum, z = round(zval, 2)) |>
        pivot_wider(names_from = stratum, values_from = z))
