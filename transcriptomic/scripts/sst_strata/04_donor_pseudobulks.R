#!/usr/bin/env Rscript
# Step 4 | Persist per-donor x per-stratum pseudobulks.
#
# Infrastructure only, no statistics. Steps 2 and 3 collapse each stratum to one
# pseudobulk per donor internally; the interaction model (step 5) and the module
# scores (step 6) instead need every donor to contribute one pseudobulk PER
# stratum, so that the diagnosis x stratum contrast is taken within donor. Those
# matrices are expensive to rebuild from the 7-cohort export, so they are written
# once here and reused.
#
# Runtime ~3 min. Writes <cohort>_donor_stratum_{counts.parquet,meta.csv}.
source("transcriptomic/scripts/sst_strata/_common.R")
suppressPackageStartupMessages(library(arrow))

st_of <- stratum_of()
for (co in export_cohorts()) {
  gr <- read_csv(file.path(P$exp, paste0(co, "_groups.csv")),
                 col_types = cols(donor = col_character(), .default = col_guess())) |>
    filter(supertype %in% names(st_of), n_cells > 0) |>
    mutate(stratum = st_of[supertype])
  f_pq <- file.path(P$exp, paste0(co, "_pseudobulk_counts.parquet"))
  idcols <- intersect(c("gene", "gene_symbol"), names(open_dataset(f_pq)))
  pb <- as.data.frame(read_parquet(f_pq, col_select = all_of(c(idcols, gr$key))))

  sym <- if ("gene_symbol" %in% names(pb)) pb$gene_symbol else pb$gene
  looks_sym <- !grepl("^ENSG", pb$gene)
  sym[is.na(sym) & looks_sym] <- pb$gene[is.na(sym) & looks_sym]
  keep <- !is.na(sym)
  M <- vapply(pb[keep, gr$key, drop = FALSE], as.numeric, numeric(sum(keep)))
  rownames(M) <- NULL
  M <- rowsum(M, group = sym[keep])

  ds <- gr |> group_by(donor, stratum) |>
    summarise(n_cells = sum(n_cells), diagnosis = first(diagnosis), .groups = "drop") |>
    filter(n_cells >= MIN_CELLS)
  agg <- sapply(seq_len(nrow(ds)), function(i)
    rowSums(M[, gr$key[gr$donor == ds$donor[i] & gr$stratum == ds$stratum[i]], drop = FALSE]))
  colnames(agg) <- paste(ds$donor, ds$stratum, sep = "|")

  write_parquet(as.data.frame(agg) |> mutate(gene = rownames(agg), .before = 1),
                file.path(P$donor, paste0(co, "_donor_stratum_counts.parquet")))
  write_csv(ds |> mutate(cohort = co), file.path(P$donor, paste0(co, "_donor_stratum_meta.csv")))
  say("%s: %d genes x %d donor-stratum samples (%d donors)",
      co, nrow(agg), ncol(agg), n_distinct(ds$donor))
}

# coverage summary -- the donor and nuclei counts quoted in the results text
ds_all <- all_donor_meta()
say("=== stratum coverage (donors with >= %d nuclei) ===", MIN_CELLS)
print(ds_all |> group_by(stratum) |>
        summarise(donors = n(), nuclei = sum(n_cells), median_per_donor = median(n_cells)) |>
        arrange(match(stratum, STRATA_LEVELS)))
print(ds_all |> distinct(cohort, donor, diagnosis) |> count(diagnosis))
