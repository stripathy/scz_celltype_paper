#!/usr/bin/env Rscript
# Preranked GSEA across all cell types using fgsea.
# Ranks all 12k+ genes per cell type by z = estimate/se, tests MSigDB
# Hallmarks + GO BP/MF/CC + Reactome. Builds direct comparison to ORA.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(cowplot); library(fgsea); library(msigdbr)
})

INPUT     <- "data/DE_genes_all_cells_scz.csv"
OUT_DIR   <- "results"
dir.create(OUT_DIR, showWarnings = FALSE)
MIN_SIZE  <- 10
MAX_SIZE  <- 500
N_PERM    <- "auto"   # fgsea adaptive
RANK_NAME <- "z = estimate / se"

# Cell-class grouping (for heatmap ordering)
EXC <- c("L2_3 IT","L4 IT","L5 IT","L5 ET","L5_6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")
INH <- c("Lamp5","Pax6","Pvalb","Sncg","Sst","Sst Chodl","Vip","Chandelier")
GLI <- c("Astro","Oligo","OPC","Micro-PVM","Endo","VLMC")

# ---- Load DE data ----------------------------------------------------------
df <- read_csv(INPUT, show_col_types = FALSE) |>
  filter(!is.na(estimate), !is.na(se), se > 0) |>
  mutate(z = estimate / se)
cat(sprintf("Loaded %d rows; %d cell types\n",
            nrow(df), length(unique(df$cell_type))))

# ---- Build gene set list ---------------------------------------------------
cat("Pulling MSigDB Hallmarks + GO + Reactome ...\n")
collect <- function(cat, subcat = NULL) {
  args <- list(species = "Homo sapiens", collection = cat)
  if (!is.null(subcat)) args$subcollection <- subcat
  do.call(msigdbr, args) |>
    transmute(gs_name, gene = gene_symbol, source = paste(cat, subcat, sep = ":"))
}
gs_h    <- collect("H")
gs_gobp <- collect("C5", "GO:BP")
gs_gocc <- collect("C5", "GO:CC")
gs_gomf <- collect("C5", "GO:MF")
gs_reac <- collect("C2", "CP:REACTOME")

all_gs <- bind_rows(gs_h, gs_gobp, gs_gocc, gs_gomf, gs_reac) |>
  distinct(gs_name, gene, .keep_all = TRUE)

# Save a mapping of gs_name -> source for later annotation
gs_source <- all_gs |> distinct(gs_name, source)

# Convert to named list for fgsea
gene_sets <- split(all_gs$gene, all_gs$gs_name)
cat(sprintf("Built %d gene sets total\n", length(gene_sets)))

# ---- Run fgsea per cell type ------------------------------------------------
cache <- file.path(OUT_DIR, "gsea_cache.rds")
if (file.exists(cache)) {
  cat("Loading cached fgsea results\n")
  all_gsea <- readRDS(cache)
} else {
  all_gsea <- list()
  cts <- sort(unique(df$cell_type))
  for (i in seq_along(cts)) {
    ct <- cts[i]
    sub <- df |> filter(cell_type == ct)
    if (nrow(sub) < 500) {
      cat(sprintf("[%d/%d] %s — skip (only %d genes tested)\n",
                  i, length(cts), ct, nrow(sub)))
      next
    }
    # Handle duplicate gene names (take max |z|)
    sub <- sub |> arrange(desc(abs(z))) |> distinct(genes, .keep_all = TRUE)
    ranks <- setNames(sub$z, sub$genes)
    ranks <- sort(ranks, decreasing = TRUE)
    cat(sprintf("[%d/%d] %s — %d genes ranked, running fgsea ...\n",
                i, length(cts), ct, length(ranks)))
    res <- fgsea(pathways = gene_sets, stats = ranks,
                 minSize = MIN_SIZE, maxSize = MAX_SIZE,
                 nPermSimple = 10000)
    res$cell_type <- ct
    all_gsea[[ct]] <- as_tibble(res)
  }
  saveRDS(all_gsea, cache)
}

combined <- bind_rows(all_gsea) |>
  left_join(gs_source, by = c("pathway" = "gs_name"))

# leadingEdge is a list column — flatten for CSV
combined_csv <- combined |>
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = "|"))
write_csv(combined_csv, file.path(OUT_DIR, "gsea_all_celltypes.csv"))
cat(sprintf("\nWrote gsea_all_celltypes.csv (%d rows)\n", nrow(combined_csv)))

# ---- Quick summary: significant pathways per cell type --------------------
sig_summary <- combined |>
  filter(padj < 0.05) |>
  group_by(cell_type) |>
  summarise(n_sig_up   = sum(NES > 0),
            n_sig_down = sum(NES < 0),
            .groups = "drop") |>
  arrange(desc(n_sig_up + n_sig_down))
cat("\n=== Significant GSEA hits per cell type (padj < 0.05) ===\n")
print(sig_summary, n = 30)
write_csv(sig_summary, file.path(OUT_DIR, "gsea_sig_summary.csv"))

# ---- Helper: order cell types ----------------------------------------------
order_ct <- function(cts) {
  ord <- c(intersect(EXC, cts), intersect(INH, cts), intersect(GLI, cts))
  factor(cts, levels = ord)
}

# ---- HEATMAP 1: focused on Sst-derived themes ------------------------------
# Map Sst ORA terms to MSigDB equivalents (gs_name is uppercased, underscored)
focus_terms <- list(
  presynaptic_vesicle = c(
    "GOCC_SYNAPSE","GOCC_PRESYNAPSE","GOCC_SYNAPTIC_VESICLE",
    "GOCC_EXOCYTIC_VESICLE","GOCC_SECRETORY_VESICLE",
    "GOCC_HIPPOCAMPAL_MOSSY_FIBER","GOCC_TERMINAL_BOUTON",
    "GOCC_TRANSPORT_VESICLE","GOCC_SECRETORY_GRANULE_MEMBRANE",
    "GOCC_SYNAPTIC_VESICLE_MEMBRANE","GOCC_POSTSYNAPSE",
    "GOCC_POSTSYNAPTIC_DENSITY","GOCC_DENDRITIC_SPINE",
    "GOBP_SYNAPTIC_VESICLE_CYCLE","GOBP_NEUROTRANSMITTER_SECRETION",
    "GOBP_CHEMICAL_SYNAPTIC_TRANSMISSION","GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING",
    "GOBP_SYNAPSE_ORGANIZATION","GOBP_SYNAPSE_ASSEMBLY",
    "REACTOME_NEURONAL_SYSTEM","REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES"
  ),
  bmp_signaling = c(
    "GOMF_I_SMAD_BINDING","GOBP_BMP_SIGNALING_PATHWAY",
    "GOBP_RESPONSE_TO_BMP","HALLMARK_TGF_BETA_SIGNALING",
    "REACTOME_SIGNALING_BY_BMP","GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY"
  ),
  cholesterol_metabolism = c(
    "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
    "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
    "REACTOME_METABOLISM_OF_STEROIDS",
    "GOBP_CHOLESTEROL_BIOSYNTHETIC_PROCESS",
    "GOBP_CHOLESTEROL_METABOLIC_PROCESS",
    "GOBP_STEROL_BIOSYNTHETIC_PROCESS",
    "WP_CHOLESTEROL_BIOSYNTHESIS_PATHWAY"
  ),
  mitochondrial = c(
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
    "GOCC_MITOCHONDRION","GOCC_MITOCHONDRIAL_INNER_MEMBRANE",
    "GOBP_MITOCHONDRIAL_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",
    "GOBP_OXIDATIVE_PHOSPHORYLATION"
  )
)

plot_focused <- function(term_set, panel_title) {
  d <- combined |> filter(pathway %in% term_set) |>
    select(cell_type, pathway, NES, padj)
  if (nrow(d) == 0) return(NULL)
  d$cell_type <- order_ct(d$cell_type)
  d$pathway   <- factor(d$pathway, levels = rev(term_set))
  ggplot(d, aes(x = cell_type, y = pathway, fill = NES)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(padj < 0.01, "**",
                          ifelse(padj < 0.05, "*",
                          ifelse(padj < 0.1,  "+", "")))),
              size = 4, vjust = 0.7) +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                         midpoint = 0, name = "NES",
                         limits = c(-3.5, 3.5), oob = scales::squish) +
    labs(x = NULL, y = NULL, title = panel_title) +
    theme_cowplot(font_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 9),
          plot.title  = element_text(size = 12, face = "bold"),
          panel.grid  = element_blank(),
          axis.line   = element_blank(),
          axis.ticks  = element_blank())
}

p_syn  <- plot_focused(focus_terms$presynaptic_vesicle,
                       "Presynaptic / synaptic gene sets — GSEA NES")
p_bmp  <- plot_focused(focus_terms$bmp_signaling,
                       "BMP / TGF-beta signaling gene sets — GSEA NES")
p_chol <- plot_focused(focus_terms$cholesterol_metabolism,
                       "Cholesterol metabolism gene sets — GSEA NES")
p_mito <- plot_focused(focus_terms$mitochondrial,
                       "Mitochondrial / OxPhos gene sets — GSEA NES")

grid_focus <- plot_grid(p_syn, p_chol, p_bmp, p_mito,
                        ncol = 1, align = "v",
                        rel_heights = c(1.3, 0.8, 0.7, 0.8))
ggsave(file.path(OUT_DIR, "heatmap_focused_GSEA.png"),
       grid_focus, width = 13, height = 18, dpi = 200, bg = "white")
ggsave(file.path(OUT_DIR, "heatmap_focused_GSEA.pdf"),
       grid_focus, width = 13, height = 18, bg = "white")

# ---- HEATMAP 2: discovery - top consistent hits across cell types ----------
top_consistent <- combined |>
  group_by(pathway, source) |>
  summarise(n_sig = sum(padj < 0.05),
            min_padj = min(padj),
            mean_NES = mean(NES, na.rm = TRUE),
            size = first(size),
            .groups = "drop") |>
  filter(n_sig >= 3, size >= 15, size <= 500) |>
  arrange(desc(n_sig), min_padj)

top_up   <- top_consistent |> filter(mean_NES > 0) |> head(25)
top_down <- top_consistent |> filter(mean_NES < 0) |> head(25)

plot_disc <- function(top_terms, panel_title) {
  d <- combined |> filter(pathway %in% top_terms$pathway) |>
    mutate(label = sprintf("%s | %s",
                           gsub("C[0-9]:|CP:", "", source),
                           pathway))
  d$cell_type <- order_ct(d$cell_type)
  d$label <- factor(d$label,
    levels = rev(sprintf("%s | %s",
                         gsub("C[0-9]:|CP:", "", top_terms$source),
                         top_terms$pathway)))
  ggplot(d, aes(x = cell_type, y = label, fill = NES)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(padj < 0.01, "**",
                          ifelse(padj < 0.05, "*",
                          ifelse(padj < 0.1,  "+", "")))),
              size = 3.5, vjust = 0.7) +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                         midpoint = 0, name = "NES",
                         limits = c(-3.5, 3.5), oob = scales::squish) +
    labs(x = NULL, y = NULL, title = panel_title) +
    theme_cowplot(font_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 8),
          plot.title  = element_text(size = 12, face = "bold"),
          panel.grid  = element_blank(),
          axis.line   = element_blank(),
          axis.ticks  = element_blank())
}

p_top_down <- plot_disc(top_down, "Top 25 consistently DOWN-regulated gene sets (GSEA, sig in >= 3 cell types)")
p_top_up   <- plot_disc(top_up,   "Top 25 consistently UP-regulated gene sets (GSEA, sig in >= 3 cell types)")

ggsave(file.path(OUT_DIR, "heatmap_discovery_down_GSEA.png"),
       p_top_down, width = 14, height = 11, dpi = 200, bg = "white")
ggsave(file.path(OUT_DIR, "heatmap_discovery_up_GSEA.png"),
       p_top_up,   width = 14, height = 11, dpi = 200, bg = "white")

# ---- Brief comparison to ORA -----------------------------------------------
cat("\n=== Side-by-side: Sst — GSEA vs ORA on key terms ===\n")
sst_gsea <- combined |> filter(cell_type == "Sst",
                               pathway %in% unlist(focus_terms)) |>
  arrange(padj) |> head(15) |>
  select(pathway, size, NES, pval, padj)
print(sst_gsea, n = 20)

cat("\nDone. Outputs in", OUT_DIR, "\n")
