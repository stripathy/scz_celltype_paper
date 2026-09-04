#!/usr/bin/env Rscript
# Step 5b | Gene-set enrichment on the interaction statistics.
#
# Split from 05 because fitting the interaction models takes ~35 min while the
# enrichment takes ~1 min; coupling them meant any change to the gene-set
# collections forced a full refit.
#
# Ranks are the meta-analytic interaction z from 05. A negative NES means the
# diagnosis effect is more negative in the depleted group than in the
# non-depleted one. Only the depleted-vs-non-depleted and linear-trend
# coefficients are ranked: the intermediate-vs-non-depleted coefficient has no
# gene at FDR < 0.10, so its ranking carries no signal worth enriching.
source("transcriptomic/scripts/fig5/_common.R")
set.seed(42)

meta <- read_csv(file.path(P$pb, "interaction_meta.csv"), show_col_types = FALSE)
gs <- msigdb_sets()
gsea <- map_dfr(c("dxSCZ:stratumdepleted", "dxSCZ:score"), function(cf) {
  d <- meta |> filter(coef == cf)
  run_fgsea(d$gene, d$zval, cf, gs$sets) |> rename(coef = signature)
})
write_csv(gsea |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = "|")),
          file.path(P$pb, "interaction_gsea.csv"))

say("=== interaction GSEA on the panel-d blocks ===")
say("    (negative NES = diagnosis effect more negative in depleted)")
print(gsea |> inner_join(BLOCKS, by = "pathway") |>
        mutate(cell = sprintf("%.2f%s", NES, stars_of(padj))) |>
        select(block, pathway, coef, cell) |>
        pivot_wider(names_from = coef, values_from = cell) |>
        arrange(match(block, BLOCK_KEYS)), n = 15, width = 200)
say("=== top 10 interaction gene sets (depleted vs non-depleted) ===")
print(gsea |> filter(coef == "dxSCZ:stratumdepleted") |> arrange(padj) |>
        transmute(pathway = str_trunc(pathway, 58), NES = round(NES, 2),
                  padj = signif(padj, 2)) |> head(10), n = 10)
