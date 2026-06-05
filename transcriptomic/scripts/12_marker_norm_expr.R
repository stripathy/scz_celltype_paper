#!/usr/bin/env Rscript
# ============================================================================
# Per-donor library-normalised expression (CP1K = counts per 1,000 transcripts,
# TMM-normalised) for the composite panel-K marker pairs — SST in Sst, PVALB in
# Pvalb — together with the edgeR DE p-value. Computed from the same edgeR object
# the rest of the figure uses: DGEList -> filterByExpr(design ~ dx+sex+age) ->
# calcNormFactors(TMM) -> cpm() / 1000.
#
# Input : ~/Github/SCZ_Xenium/output/de/pseudobulk_subclass{,_samples}.csv  (build_de_input.py)
#         ../spatial/output/de/de_results_subclass.csv                       (run_de.R)
# Output: results/tables/marker_norm_expr.csv        gene, celltype, donor, dx, cp1k
#         results/tables/marker_norm_expr_stats.csv  gene, celltype, logFC, p, fdr
# ============================================================================
suppressPackageStartupMessages({library(edgeR); library(dplyr); library(tidyr); library(readr)})

PB  <- path.expand("~/Github/SCZ_Xenium/output/de")
DE  <- "../spatial/output/de/de_results_subclass.csv"
OUT <- "results/tables"
PAIRS <- list(c("SST", "Sst"), c("FGFR3", "Astrocyte"), c("PVALB", "Pvalb"))

long <- read_csv(file.path(PB, "pseudobulk_subclass.csv"), show_col_types = FALSE)
samp <- read_csv(file.path(PB, "pseudobulk_subclass_samples.csv"), show_col_types = FALSE)
de   <- read_csv(DE, show_col_types = FALSE)

expr_rows <- list(); stat_rows <- list()
for (pr in PAIRS) {
  gene <- pr[1]; ct <- pr[2]
  w <- long[long$celltype == ct, c("gene", "donor", "count")] |>
       pivot_wider(names_from = donor, values_from = count)
  mat <- as.matrix(w[, -1]); rownames(mat) <- w$gene
  sm <- samp[samp$celltype == ct, ]; sm <- sm[match(colnames(mat), sm$donor), ]
  dx <- factor(sm$diagnosis, levels = c("Control", "SCZ"))
  sex <- factor(sm$sex); age_c <- sm$age - mean(sm$age)
  design <- model.matrix(~ dx + sex + age_c)
  d <- DGEList(counts = mat)
  keep <- filterByExpr(d, design); d <- d[keep, , keep.lib.sizes = FALSE]
  d <- calcNormFactors(d, method = "TMM")
  cp1k <- cpm(d, log = FALSE)[gene, ] / 1000        # counts per 1,000 transcripts (TMM-normalised)
  expr_rows[[paste(gene, ct)]] <- tibble(gene = gene, celltype = ct, donor = colnames(mat),
                                         dx = as.character(dx), cp1k = as.numeric(cp1k))
  dr <- de[de$gene == gene & de$celltype == ct, ]
  stat_rows[[paste(gene, ct)]] <- tibble(gene = gene, celltype = ct,
                                         logFC = dr$logFC[1], p = dr$PValue[1], fdr = dr$FDR[1])
}
write_csv(bind_rows(expr_rows), file.path(OUT, "marker_norm_expr.csv"))
write_csv(bind_rows(stat_rows), file.path(OUT, "marker_norm_expr_stats.csv"))
cat("Saved marker_norm_expr.csv + marker_norm_expr_stats.csv\n")
print(as.data.frame(bind_rows(stat_rows)), digits = 3)
