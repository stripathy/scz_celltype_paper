#!/usr/bin/env Rscript
# Visualize how OxPhos comes up as enriched in Sst even though GOT2
# (and most OxPhos genes) don't pass individual padj < 0.05.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(cowplot)
  library(fgsea); library(msigdbr)
})

df <- read_csv("data/DE_genes_all_cells_scz.csv", show_col_types = FALSE) |>
  filter(!is.na(estimate), !is.na(se), se > 0) |>
  mutate(z = estimate / se)

oxphos <- msigdbr(species = "Homo sapiens", collection = "H") |>
  filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") |>
  pull(gene_symbol) |> unique()

# Sst ranks
sst <- df |> filter(cell_type == "Sst") |>
  arrange(desc(abs(z))) |> distinct(genes, .keep_all = TRUE)
ranks <- setNames(sst$z, sst$genes)
ranks <- sort(ranks, decreasing = TRUE)

# fgsea — get the leading edge for OxPhos
set.seed(42)
res <- fgsea(pathways = list(OXPHOS = oxphos),
             stats = ranks, minSize = 10, maxSize = 500,
             nPermSimple = 10000)
le_genes <- res$leadingEdge[[1]]
cat(sprintf("Sst OXPHOS GSEA: NES = %.2f, padj = %.2g, leading-edge size = %d\n",
            res$NES, res$padj, length(le_genes)))

# Build enrichment plot manually so we can annotate GOT2
n <- length(ranks)
hit <- names(ranks) %in% oxphos
N_R <- sum(abs(ranks[hit]))
inc <- ifelse(hit, abs(ranks)/N_R, -1/(n - sum(hit)))
runES <- cumsum(inc)
peak_idx <- which.min(runES)   # negative enrichment

es_df <- tibble(idx = seq_along(ranks), z = ranks, hit = hit, runES = runES)
hit_df <- es_df |> filter(hit)
got2_idx <- which(names(ranks) == "GOT2")

p_top <- ggplot(es_df, aes(idx, runES)) +
  geom_line(colour = "#0072B2", linewidth = 0.7) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = peak_idx, colour = "red", linetype = "dashed",
             linewidth = 0.4) +
  annotate("point", x = peak_idx, y = min(runES), colour = "red", size = 2) +
  annotate("text", x = peak_idx, y = min(runES) - 0.04,
           label = sprintf("Peak ES = %.2f (at rank %d)\nLeading-edge = ranks 1 - %d (%d OxPhos genes)",
                           min(runES), peak_idx, peak_idx, length(le_genes)),
           hjust = -0.05, vjust = 1, size = 3.5, colour = "red") +
  labs(x = NULL, y = "Running enrichment score",
       title = "Sst: HALLMARK_OXIDATIVE_PHOSPHORYLATION leading-edge analysis",
       subtitle = sprintf("NES = %.2f, padj = %.2g  -  164 OxPhos genes among %d ranked",
                          res$NES, res$padj, n)) +
  theme_cowplot(font_size = 12) +
  theme(plot.subtitle = element_text(size = 11, colour = "grey20"))

# tick mark plot
p_ticks <- ggplot(hit_df, aes(x = idx)) +
  geom_segment(aes(xend = idx, y = 0, yend = 1),
               colour = "#0072B2", linewidth = 0.3, alpha = 0.7) +
  geom_segment(data = tibble(x = got2_idx),
               aes(x = x, xend = x, y = -0.2, yend = 1.2),
               colour = "darkred", linewidth = 0.8) +
  annotate("text", x = got2_idx, y = 1.5, label = "GOT2",
           colour = "darkred", size = 3.6, hjust = 0.2, fontface = "bold") +
  scale_y_continuous(limits = c(-0.3, 1.8), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_void()

p_z <- ggplot(es_df, aes(idx, z)) +
  geom_col(fill = "grey70", width = 1) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
  geom_vline(xintercept = peak_idx, colour = "red", linetype = "dashed",
             linewidth = 0.3) +
  geom_vline(xintercept = got2_idx, colour = "darkred", linewidth = 0.4) +
  geom_hline(yintercept = c(-1.96, 1.96), colour = "grey40",
             linetype = "dotted", linewidth = 0.3) +
  annotate("text", x = n*0.98, y = 2.4, label = "|z|>1.96 (p<0.05 nominal)",
           hjust = 1, size = 3, colour = "grey40") +
  labs(x = "Gene rank (high z = up in SCZ, low z = down)",
       y = "DE z-score") +
  theme_cowplot(font_size = 11)

grid <- plot_grid(p_top, p_ticks, p_z, ncol = 1, align = "v",
                  rel_heights = c(2.5, 0.6, 1.5))
ggsave("results/06_explainer_Sst_OXPHOS.png", grid,
       width = 12, height = 8, dpi = 220, bg = "white")
ggsave("results/06_explainer_Sst_OXPHOS.pdf", grid,
       width = 12, height = 8, bg = "white")
cat("Wrote gsea_explainer_Sst_OXPHOS.{png,pdf}\n")

# Also: confirm GOT2 is in the leading edge by definition
cat(sprintf("\nGOT2 in leading-edge? %s\n", "GOT2" %in% le_genes))
cat(sprintf("GOT2 rank in Sst: %d / %d (z = %.3f, padj = %.3f)\n",
            got2_idx, n, ranks["GOT2"], df |>
              filter(cell_type == "Sst", genes == "GOT2") |>
              pull(padj)))
cat(sprintf("Peak ES rank: %d - any OxPhos gene with rank <= %d is leading-edge\n",
            peak_idx, peak_idx))

# Show the spread of padj values among leading-edge OXPHOS genes
le_stats <- df |> filter(cell_type == "Sst", genes %in% le_genes) |>
  select(genes, z, pval, padj) |> arrange(z)
cat(sprintf("\n=== %d OxPhos genes in leading edge: distribution of individual sig ===\n",
            nrow(le_stats)))
cat(sprintf("  padj < 0.05:  %d (%.0f%%)\n",
            sum(le_stats$padj < 0.05), 100*mean(le_stats$padj < 0.05)))
cat(sprintf("  padj < 0.10:  %d (%.0f%%)\n",
            sum(le_stats$padj < 0.10), 100*mean(le_stats$padj < 0.10)))
cat(sprintf("  padj < 0.25:  %d (%.0f%%)\n",
            sum(le_stats$padj < 0.25), 100*mean(le_stats$padj < 0.25)))
cat(sprintf("  padj >= 0.25: %d (%.0f%%) - these would all be filtered out by ORA!\n",
            sum(le_stats$padj >= 0.25), 100*mean(le_stats$padj >= 0.25)))
