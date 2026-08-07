#!/usr/bin/env Rscript
# Two analyses in one script:
#   (1) GSEA of PGC3 SCZ GWAS gene set across cell types (de-novo fgsea run)
#   (2) Leading-edge analysis of major themes from prior GSEA (reuses cache)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(cowplot); library(fgsea); library(msigdbr)
})

INPUT     <- "data/DE_genes_all_cells_scz.csv"
GWAS_FILE <- "../genetics/data/gwas/scz_gwas_gene_set_no_mhc.csv"  # internal: genetics/ owns the GWAS set
GWAS_MHC  <- "../genetics/data/gwas/scz_gwas_gene_set.csv"
OUT_DIR   <- "results"
MIN_SIZE  <- 5
MAX_SIZE  <- 1000

EXC <- c("L2_3 IT","L4 IT","L5 IT","L5 ET","L5_6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")
INH <- c("Lamp5","Pax6","Pvalb","Sncg","Sst","Sst Chodl","Vip","Chandelier")
GLI <- c("Astro","Oligo","OPC","Micro-PVM","Endo","VLMC")
order_ct <- function(cts) {
  ord <- c(intersect(EXC, cts), intersect(INH, cts), intersect(GLI, cts))
  factor(cts, levels = ord)
}

# ============================================================================
# PART 1 — SCZ GWAS gene set GSEA
# ============================================================================
cat("=== Part 1: SCZ GWAS gene set GSEA ===\n")

gwas_nomhc <- read_csv(GWAS_FILE, show_col_types = FALSE)$gene
gwas_with  <- read_csv(GWAS_MHC,  show_col_types = FALSE)$gene
cat(sprintf("PGC3 SCZ GWAS genes (no MHC): %d\n", length(gwas_nomhc)))
cat(sprintf("PGC3 SCZ GWAS genes (with MHC): %d\n", length(gwas_with)))

# Also build a couple of negative-control reference sets
hallmark_translation <- msigdbr(species = "Homo sapiens", collection = "C2",
                                subcollection = "CP:REACTOME") |>
  filter(gs_name == "REACTOME_TRANSLATION") |> pull(gene_symbol) |> unique()
hallmark_oxphos <- msigdbr(species = "Homo sapiens", collection = "H") |>
  filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") |> pull(gene_symbol) |> unique()

custom_sets <- list(
  PGC3_SCZ_GWAS_no_MHC = gwas_nomhc,
  PGC3_SCZ_GWAS_with_MHC = gwas_with,
  REF_REACTOME_TRANSLATION = hallmark_translation,
  REF_HALLMARK_OXPHOS      = hallmark_oxphos
)

df <- read_csv(INPUT, show_col_types = FALSE) |>
  filter(!is.na(estimate), !is.na(se), se > 0) |>
  mutate(z = estimate / se)

# How many GWAS genes appear in each cell type's tested gene universe?
cat("\nGWAS gene coverage in tested DE gene universe (per cell type):\n")
coverage <- df |> group_by(cell_type) |>
  summarise(n_tested = n_distinct(genes),
            gwas_in_tested = sum(unique(genes) %in% gwas_nomhc),
            .groups = "drop") |>
  arrange(desc(gwas_in_tested))
print(coverage, n = 30)

# Run fgsea per cell type with custom sets
cts <- sort(unique(df$cell_type))
all_gwas <- list()
for (i in seq_along(cts)) {
  ct  <- cts[i]
  sub <- df |> filter(cell_type == ct) |>
    arrange(desc(abs(z))) |> distinct(genes, .keep_all = TRUE)
  ranks <- setNames(sub$z, sub$genes)
  ranks <- sort(ranks, decreasing = TRUE)
  res <- fgsea(pathways = custom_sets, stats = ranks,
               minSize = MIN_SIZE, maxSize = MAX_SIZE,
               nPermSimple = 50000)  # more perms for precise small-set p-values
  res$cell_type <- ct
  all_gwas[[ct]] <- as_tibble(res)
  cat(sprintf("[%d/%d] %s\n", i, length(cts), ct))
}
gwas_results <- bind_rows(all_gwas) |>
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = "|"))
write_csv(gwas_results, file.path(OUT_DIR, "gwas_gsea_results.csv"))

cat("\n=== PGC3 SCZ GWAS (no MHC) — GSEA per cell type ===\n")
print(gwas_results |> filter(pathway == "PGC3_SCZ_GWAS_no_MHC") |>
        arrange(padj) |>
        select(cell_type, size, NES, pval, padj), n = 30)

# Plot GWAS enrichment heatmap
plot_df <- gwas_results |> filter(pathway %in% c("PGC3_SCZ_GWAS_no_MHC",
                                                  "PGC3_SCZ_GWAS_with_MHC")) |>
  mutate(cell_type = order_ct(cell_type),
         pathway   = factor(pathway,
                            levels = c("PGC3_SCZ_GWAS_with_MHC",
                                       "PGC3_SCZ_GWAS_no_MHC"),
                            labels = c("PGC3 (with MHC)", "PGC3 (no MHC)")),
         sig = case_when(padj < 0.001 ~ "***",
                          padj < 0.01  ~ "**",
                          padj < 0.05  ~ "*",
                          padj < 0.1   ~ "+",
                          TRUE         ~ ""))

p_gwas <- ggplot(plot_df, aes(x = cell_type, y = pathway, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig), size = 4.2, vjust = 0.6) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, name = "NES",
                       limits = c(-2.5, 2.5), oob = scales::squish) +
  labs(x = NULL, y = NULL,
       title = "Are SCZ DE genes enriched for PGC3 SCZ GWAS risk genes?",
       subtitle = "GSEA of 484 PGC3 (Trubetskoy 2022) genes against per-cell-type DE z-scores",
       caption = "+ padj<0.1   * padj<0.05   ** padj<0.01   *** padj<0.001    orange=up-shifted, blue=down-shifted") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y   = element_text(size = 11),
        plot.title    = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, colour = "grey20"),
        plot.caption  = element_text(size = 10, colour = "grey40"),
        panel.grid    = element_blank(), axis.line = element_blank(),
        axis.ticks    = element_blank())

ggsave(file.path(OUT_DIR, "GWAS_enrichment_heatmap.png"),
       p_gwas, width = 13, height = 4.5, dpi = 220, bg = "white")
ggsave(file.path(OUT_DIR, "GWAS_enrichment_heatmap.pdf"),
       p_gwas, width = 13, height = 4.5, bg = "white")
cat("Wrote GWAS_enrichment_heatmap.{png,pdf}\n")

# ============================================================================
# PART 2 — Leading-edge analysis
# ============================================================================
cat("\n=== Part 2: Leading-edge analysis ===\n")

cache <- file.path(OUT_DIR, "gsea_cache.rds")
if (!file.exists(cache)) stop("No gsea_cache.rds — run gsea_all_celltypes.R first.")
all_gsea <- readRDS(cache)
combined <- bind_rows(all_gsea)

# Major themes to inspect
themes <- list(
  mitochondrial = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  synaptic      = "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES",
  cholesterol   = "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
  bmp           = "REACTOME_SIGNALING_BY_BMP",
  axon_guidance = "REACTOME_SIGNALING_BY_ROBO_RECEPTORS",
  growth_stress = "HALLMARK_MYC_TARGETS_V1",
  translation   = "REACTOME_TRANSLATION"
)

# Pull leading-edge gene sets per (cell type, theme), with padj for filtering
le_long <- combined |> filter(pathway %in% unlist(themes)) |>
  filter(padj < 0.1)  # only cell types where the theme is at least marginally sig
le_long <- le_long |> rowwise() |>
  mutate(genes = list(leadingEdge)) |>
  ungroup() |>
  select(cell_type, pathway, NES, padj, genes)

# Save expanded long-form table: one row per gene per cell-type per theme
le_genes <- le_long |>
  unnest(genes) |>
  rename(gene = genes)
write_csv(le_genes, file.path(OUT_DIR, "leading_edge_genes_long.csv"))
cat(sprintf("Saved %s (%d rows)\n",
            file.path(OUT_DIR, "leading_edge_genes_long.csv"), nrow(le_genes)))

# For each theme, find genes that appear in leading edge of >= N cell types
rank_recurrent <- function(theme_name, gs_name, n_min = 3) {
  d <- le_genes |> filter(pathway == gs_name) |>
    group_by(gene) |>
    summarise(n_cells = n(),
              cells   = paste(sort(unique(cell_type)), collapse = ","),
              direction = if_else(mean(NES) < 0, "down", "up"),
              .groups = "drop") |>
    arrange(desc(n_cells)) |>
    filter(n_cells >= n_min)
  d$theme <- theme_name
  d$pathway <- gs_name
  d
}

recurrent <- imap_dfr(themes, ~ rank_recurrent(.y, .x, n_min = 3))
write_csv(recurrent, file.path(OUT_DIR, "leading_edge_recurrent_genes.csv"))

cat("\n=== Top recurrent leading-edge genes per theme (in >=3 cell types) ===\n")
for (t in names(themes)) {
  d <- recurrent |> filter(theme == t) |> head(15)
  if (nrow(d) == 0) { cat(sprintf("\n[%s] no genes meeting threshold\n", t)); next }
  cat(sprintf("\n[%s] direction = %s (top by recurrence)\n", t, d$direction[1]))
  print(d |> select(gene, n_cells, cells), n = 20)
}

# ============================================================================
# PART 3 — Do leading-edge genes overlap PGC3 SCZ GWAS?
# ============================================================================
cat("\n=== Part 3: Overlap of leading-edge themes with PGC3 GWAS genes ===\n")

gwas_set <- gwas_nomhc

overlap_summary <- le_genes |>
  group_by(pathway, cell_type) |>
  summarise(n_le        = n(),
            n_le_gwas   = sum(gene %in% gwas_set),
            le_gwas_genes = paste(sort(unique(gene[gene %in% gwas_set])), collapse = ","),
            .groups = "drop") |>
  filter(n_le_gwas > 0) |>
  arrange(desc(n_le_gwas))

write_csv(overlap_summary,
          file.path(OUT_DIR, "leading_edge_x_gwas_overlap.csv"))
cat(sprintf("Saved %s\n",
            file.path(OUT_DIR, "leading_edge_x_gwas_overlap.csv")))

cat("\nTop 25 (theme x cell type) leading-edge sets with PGC3 GWAS overlap:\n")
print(overlap_summary |> head(25), n = 25)

# Per theme: which GWAS genes recur in leading edges?
cat("\n=== Recurrent PGC3 GWAS genes within leading-edge of each theme ===\n")
gwas_recurrent <- le_genes |>
  filter(gene %in% gwas_set) |>
  group_by(pathway, gene) |>
  summarise(n_cells = n_distinct(cell_type),
            cells   = paste(sort(unique(cell_type)), collapse = ","),
            .groups = "drop") |>
  arrange(pathway, desc(n_cells))
for (t in names(themes)) {
  gs <- themes[[t]]
  d <- gwas_recurrent |> filter(pathway == gs) |> head(15)
  if (nrow(d) == 0) { cat(sprintf("\n[%s] no GWAS overlap\n", t)); next }
  cat(sprintf("\n[%s = %s]\n", t, gs))
  print(d |> select(gene, n_cells, cells), n = 15)
}
write_csv(gwas_recurrent, file.path(OUT_DIR, "gwas_in_leading_edge_recurrent.csv"))

cat("\nDone.\n")
