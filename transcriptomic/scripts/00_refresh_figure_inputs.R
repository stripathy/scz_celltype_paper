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
  list(file = "DE_genes_all_cells_scz.csv",
       source = "../shared/snrnaseq_de/DE_genes_all_cells_scz.csv",
       how = "copy",
       what = "meta-analytic snRNA-seq DE (7 cohorts), full gene x subclass table"),
  list(file = "meta_results_cohorts_subclass_forest.csv",
       source = "../shared/snrnaseq_de/meta_results_cohorts_subclass.csv",
       how = "subset_forest",
       what = "per-cohort snRNA-seq DE, SST + PVALB rows only (forest panels b, f)"),
  list(file = "de_results_subclass.csv",
       source = file.path(XEN, "de/de_results_subclass.csv"),
       how = "copy",
       what = "Xenium pseudobulk DE, subclass level (panels b, f, j)"),
  list(file = "crumblr_input_subclass_corr.csv",
       source = file.path(XEN, "crumblr/crumblr_input_subclass_corr.csv"),
       how = "copy",
       what = "Xenium per-donor subclass composition (panel i inset)"),
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
