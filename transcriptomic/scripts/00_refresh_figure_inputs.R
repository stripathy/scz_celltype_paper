#!/usr/bin/env Rscript
# ============================================================================
# 00_refresh_figure_inputs.R — resync the committed Figure 2 input snapshots
# ============================================================================
# Copies each canonical analysis output into data/figure_inputs/ and records
# its checksum in MANIFEST.tsv. Every Figure 2 script reads through
# fig_input() (see _figure_inputs.R), which refuses to run when a snapshot no
# longer matches the source it came from — so this script is the single step
# that has to be run after any upstream DE, crumblr, or meta-analysis rerun.
#
# One entry is not a plain copy: the per-cohort DE table is ~2.2M rows, so only
# the SST + PVALB rows the forest panels need are stored.
#
# Usage:  Rscript scripts/00_refresh_figure_inputs.R        (from transcriptomic/)
# ============================================================================
suppressPackageStartupMessages({ library(readr); library(dplyr) })

OUT <- "data/figure_inputs"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

XEN <- path.expand("~/Github/SCZ_Xenium/output")

# file in figure_inputs        canonical source                                        how
SPEC <- list(
  # The snRNA-seq pipeline has lived in ../snrnaseq/ since 2026-09-09, so its
  # outputs are the sources now; the shared/snrnaseq_de/ copies predate the
  # 9 Sep rerun and must not be used, or a refresh would silently revert the
  # snapshots to the old numbers.
  list(file = "DE_genes_all_cells_scz.csv",
       source = "../snrnaseq/snRNAseq_DE/Files/DE_genes_all_cells_scz.csv",
       how = "copy",
       what = "meta-analytic snRNA-seq DE (7 datasets), full gene x subclass table"),
  list(file = "plotdata.csv",
       source = "../snrnaseq/Final_figures/Data/plotdata.csv",
       how = "copy",
       what = "composition meta-analysis plot table; only the per-dataset donor n for Sst_25 is used (forest labels, panels b, f)"),
  list(file = "xen_Sst_proportions.csv",
       source = "../snrnaseq/Final_figures/Data/xen_Sst_proportions.csv",
       how = "copy",
       what = "Xenium per-donor Sst supertype proportions; only the Sst_25 donor count is used (forest labels)"),
  list(file = "snrnaseq_subclass_mean_prop.csv",
       source = "../snrnaseq/Compositional_analysis/Files/7_cohorts_metadata_names.csv",
       how = "subclass_prop",
       what = "mean per-donor share of nuclei per SEA-AD subclass across the 7 datasets, 469 donors (panel i inset)"),
  list(file = "meta_results_cohorts_subclass_forest.csv",
       source = "../snrnaseq/snRNAseq_DE/Files/DE_genes_all_cohorts_subclass.csv",   # git-lfs; 2.2M rows
       how = "subset_forest",
       what = "per-cohort snRNA-seq DE, SST + PVALB rows only (forest panels b, f)"),
  list(file = "de_results_subclass.csv",
       source = file.path(XEN, "de/de_results_subclass.csv"),
       how = "copy",
       what = "Xenium pseudobulk DE, subclass level (panels b, f, j)"),
  list(file = "crumblr_input_subclass_corr.csv",
       source = file.path(XEN, "crumblr/crumblr_input_subclass_corr.csv"),
       how = "copy",
       what = "Xenium per-donor subclass composition (not read by Figure 2 since 2026-09-14; kept for the standalone scripts and KNOWN_ISSUES #18)"),
  list(file = "pseudobulk_subclass.csv",
       source = file.path(XEN, "de/pseudobulk_subclass.csv"),
       how = "copy",
       what = "Xenium pseudobulk counts (normalised marker expression, panels c, g)"),
  list(file = "pseudobulk_subclass_samples.csv",
       source = file.path(XEN, "de/pseudobulk_subclass_samples.csv"),
       how = "copy",
       what = "Xenium pseudobulk sample metadata (donor, diagnosis, sex, age, pmi)")
)

FOREST_GENES <- c("SST", "PVALB")

rows <- list()
for (s in SPEC) {
  src <- path.expand(s$source)
  dst <- file.path(OUT, s$file)
  if (!file.exists(src)) {
    cat(sprintf("  SKIP  %-42s source not found: %s\n", s$file, s$source))
    next
  }
  if (s$how == "subset_forest") {
    d <- read_csv(src, show_col_types = FALSE)
    gcol <- if ("genes" %in% names(d)) "genes" else "gene"
    d <- d[d[[gcol]] %in% FOREST_GENES, ]
    write_csv(d, dst)
    note <- sprintf("subset to %s (%d rows)", paste(FOREST_GENES, collapse = "+"), nrow(d))
  } else if (s$how == "subclass_prop") {
    # Per-donor supertype nucleus counts -> per-donor share of nuclei per subclass
    # -> mean over donors. Supertype collapses to subclass by the pipeline's own
    # rule (snrnaseq/snRNAseq_DE/*/1_Pseudobulk.r: gsub("_[0-9].*$", "", id)),
    # which also folds the "-SEAAD" extended supertypes; then the DE table's
    # spellings, so the inset joins on cell_type without a lookup.
    d <- read_csv(src, show_col_types = FALSE)
    cnt <- d[, !(names(d) %in% c("Cohort", "Donor", "Age", "Sex", "Diagnosis", "PMI")) &
                !grepl("^\\.\\.\\.", names(d))]
    sub <- gsub("_[0-9].*$", "", names(cnt))
    sub <- dplyr::recode(sub, "L2/3 IT" = "L2_3 IT", "L5/6 NP" = "L5_6 NP", "Lamp5 Lhx6" = "Lamp5_Lhx6")
    m <- as.matrix(cnt); m[is.na(m)] <- 0
    by_sub <- sapply(split(seq_len(ncol(m)), sub), function(j) rowSums(m[, j, drop = FALSE]))
    tot <- rowSums(by_sub); keep <- tot > 0
    prop <- colMeans(by_sub[keep, , drop = FALSE] / tot[keep])
    out <- tibble::tibble(CellType = names(prop), mean_prop = unname(prop), n_donors = sum(keep)) |>
      arrange(desc(mean_prop))
    write_csv(out, dst)
    note <- sprintf("%d subclasses, mean over %d donors", nrow(out), sum(keep))
  } else {
    file.copy(src, dst, overwrite = TRUE)
    note <- sprintf("%.1f MB", file.size(dst) / 1e6)
  }
  rows[[length(rows) + 1]] <- data.frame(
    file = s$file,
    source = s$source,
    source_md5 = unname(tools::md5sum(src)),
    snapshot_md5 = unname(tools::md5sum(dst)),
    source_mtime = format(file.mtime(src), "%Y-%m-%d %H:%M:%S"),
    refreshed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    description = s$what,
    stringsAsFactors = FALSE)
  cat(sprintf("  ok    %-42s %s\n", s$file, note))
}

man <- bind_rows(rows)
write.table(man, file.path(OUT, "MANIFEST.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nWrote %s/MANIFEST.tsv (%d entries)\n", OUT, nrow(man)))
cat("Now re-render the figures: scripts 01, 07, 08, 09, 12, 16\n")
