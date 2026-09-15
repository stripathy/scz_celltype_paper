#!/usr/bin/env Rscript
# Step 3 | Preranked GSEA on the stratum meta-DE.
#
# Ranks are z = estimate/SE from step 2; house settings (Hallmark + GO BP/CC/MF +
# Reactome, size 10-500) match the subclass pipeline so burdens are comparable
# with Figure 2. The pooled 16-supertype analysis is emitted under the label
# "Sst_subclass" and serves as the reference column in panels b, d and e.
#
# Runtime ~8 min. Writes the two canonical tables the figure reads.
source("transcriptomic/scripts/sst_strata/_common.R")
set.seed(42)

meta <- read_csv(file.path(P$pb, "stratum_meta_de.csv"), show_col_types = FALSE)

# gene-level signature table (long, one row per gene x stratum)
sig <- meta |> filter(stratum %in% STRATA_LEVELS) |>
  transmute(stratum, gene, k, estimate, se, z = zval, padj)
write_csv(sig, file.path(P$pb, "stratum_gene_signatures.csv"))

gs <- msigdb_sets()
gsea <- map_dfr(c(STRATA_LEVELS, "all_sst"), function(st) {
  d <- meta |> filter(stratum == st)
  run_fgsea(d$gene, d$zval, if (st == "all_sst") "Sst_subclass" else st, gs$sets)
}) |> left_join(gs$source, by = c(pathway = "gs_name"))
write_csv(gsea |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = "|")),
          file.path(P$pb, "gsea_all_signatures.csv"))

say("=== gene-set burden ===")
print(gsea |> group_by(signature) |>
        summarise(n05 = sum(padj < 0.05), n05_dn = sum(padj < 0.05 & NES < 0),
                  n10 = sum(padj < 0.10), n10_dn = sum(padj < 0.10 & NES < 0),
                  .groups = "drop") |>
        arrange(match(signature, SIG_LEVELS)))

say("=== the blocks shown in panel d ===")
print(gsea |> inner_join(BLOCKS, by = "pathway") |>
        mutate(cell = sprintf("%.2f%s", NES, stars_of(padj))) |>
        select(block, pathway, signature, cell) |>
        pivot_wider(names_from = signature, values_from = cell) |>
        arrange(match(block, BLOCK_KEYS)), n = 15, width = 200)

mods <- modules(gsea)
say("module sizes (leading-edge unions): %s",
    paste(sprintf("%s %d", names(mods), lengths(mods)), collapse = ", "))
