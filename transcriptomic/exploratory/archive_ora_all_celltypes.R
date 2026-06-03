#!/usr/bin/env Rscript
# Enrichment across ALL cell types: does the Sst presynaptic/vesicle signature
# generalize? Per-cell-type background, FDR correction.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(cowplot); library(gprofiler2)
})

INPUT     <- "data/DE_genes_all_cells_scz.csv"
PADJ_THR  <- 0.1
MIN_GENES <- 20                     # skip lists with fewer than this
SOURCES   <- c("GO:BP","GO:MF","GO:CC","REAC","KEGG","WP")
OUT_DIR   <- "results/exploratory/archive_ora_all_celltypes"
dir.create(OUT_DIR, showWarnings = FALSE)

# Cell-class grouping (for heatmap ordering)
EXC <- c("L2_3 IT","L4 IT","L5 IT","L5 ET","L5_6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")
INH <- c("Lamp5","Pax6","Pvalb","Sncg","Sst","Sst Chodl","Vip","Chandelier")
GLI <- c("Astro","Oligo","OPC","Micro-PVM","Endo","VLMC")

df <- read_csv(INPUT, show_col_types = FALSE)

celltype_lists <- df |>
  filter(padj < PADJ_THR) |>
  mutate(dir = if_else(estimate > 0, "up", "down")) |>
  group_by(cell_type, dir) |>
  summarise(genes = list(genes), n = n(), .groups = "drop") |>
  filter(n >= MIN_GENES) |>
  mutate(key = sprintf("%s|%s", cell_type, dir))

cat(sprintf("Will run %d enrichments (cell_type x direction with n>=%d)\n",
            nrow(celltype_lists), MIN_GENES))
print(celltype_lists |> select(cell_type, dir, n) |>
        arrange(dir, desc(n)), n = 50)

# Cache background gene lists per cell type
bg_by_ct <- df |> group_by(cell_type) |>
  summarise(bg = list(unique(genes)), .groups = "drop") |>
  { \(x) setNames(x$bg, x$cell_type) }()

run_one <- function(genes, bg, label) {
  res <- tryCatch(
    gost(query = genes, organism = "hsapiens",
         custom_bg = bg, sources = SOURCES,
         correction_method = "fdr",
         user_threshold = 1, significant = FALSE,
         evcodes = FALSE),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(res) || is.null(res$result) || nrow(res$result) == 0) return(NULL)
  res$result |> as_tibble() |>
    select(source, term_name, term_id, p_value, term_size, intersection_size) |>
    mutate(query_set = label)
}

cache_file <- file.path(OUT_DIR, "enrich_cache.rds")
if (file.exists(cache_file)) {
  cat("Loading cached results from", cache_file, "\n")
  all_res <- readRDS(cache_file)
} else {
  all_res <- list()
  for (i in seq_len(nrow(celltype_lists))) {
    row <- celltype_lists[i, ]
    cat(sprintf("[%d/%d] %s  (n=%d)\n", i, nrow(celltype_lists), row$key, row$n))
    bg <- bg_by_ct[[row$cell_type]]
    r  <- run_one(row$genes[[1]], bg, row$key)
    if (!is.null(r)) all_res[[row$key]] <- r
    Sys.sleep(0.5)  # be gentle on the web API
  }
  saveRDS(all_res, cache_file)
}

combined <- bind_rows(all_res) |>
  separate(query_set, c("cell_type","dir"), sep = "\\|", remove = FALSE)
write_csv(combined, file.path(OUT_DIR, "enrich_all_celltypes_fdr.csv"))
cat(sprintf("\nWrote %s (%d rows)\n",
            file.path(OUT_DIR, "enrich_all_celltypes_fdr.csv"), nrow(combined)))

# ---------- HEATMAP 1: Sst-derived terms across all cell types ---------------
sst_terms <- list(
  down = c(
    "secretory vesicle", "synapse", "exocytic vesicle", "presynapse",
    "synaptic vesicle", "secretory granule membrane", "hippocampal mossy fiber",
    "terminal bouton", "transport vesicle", "endocytic vesicle",
    "synaptic vesicle membrane", "cell junction"
  ),
  up = c(
    "I-SMAD binding", "neurexin family protein binding",
    "P-type calcium transporter activity",
    "calcium ion transmembrane transporter activity",
    "metal ion transmembrane transporter activity",
    "L-glutamine transmembrane transporter activity",
    "lipid binding"
  )
)

build_heatmap_data <- function(direction, term_set) {
  combined |>
    filter(dir == direction, term_name %in% term_set) |>
    select(cell_type, term_name, p_value, intersection_size, term_size) |>
    mutate(neglog10p = -log10(p_value))
}

order_ct <- function(cts) {
  ord <- c(intersect(EXC, cts), intersect(INH, cts), intersect(GLI, cts))
  factor(cts, levels = ord)
}

plot_focused_heatmap <- function(direction, term_set, fill_col, title) {
  d <- build_heatmap_data(direction, term_set)
  if (nrow(d) == 0) return(NULL)
  d$cell_type <- order_ct(d$cell_type)
  d$term_name <- factor(d$term_name, levels = rev(term_set))
  ggplot(d, aes(x = cell_type, y = term_name, fill = neglog10p)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(p_value < 0.05, "*",
                          ifelse(p_value < 0.1,  "+", ""))),
              size = 5, vjust = 0.7) +
    scale_fill_gradient(low = "grey95", high = fill_col,
                        name = expression(-log[10]~"FDR"),
                        limits = c(0, NA)) +
    labs(x = NULL, y = NULL, title = title) +
    theme_cowplot(font_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11),
          plot.title  = element_text(size = 13, face = "bold"),
          panel.grid  = element_blank(),
          axis.line   = element_blank(),
          axis.ticks  = element_blank(),
          legend.position = "right")
}

p_focus_down <- plot_focused_heatmap(
  "down", sst_terms$down, "#0072B2",
  "Sst-derived 'down' terms: synaptic/vesicle signature across cell types")
p_focus_up <- plot_focused_heatmap(
  "up", sst_terms$up, "#D55E00",
  "Sst-derived 'up' terms: BMP / transport / synaptic adhesion across cell types")

focus_grid <- plot_grid(p_focus_down, p_focus_up, ncol = 1, align = "v",
                        rel_heights = c(1.05, 0.85))
ggsave(file.path(OUT_DIR, "heatmap_sst_terms_across_celltypes.png"),
       focus_grid, width = 13, height = 11, dpi = 220, bg = "white")
ggsave(file.path(OUT_DIR, "heatmap_sst_terms_across_celltypes.pdf"),
       focus_grid, width = 13, height = 11, bg = "white")

# ---------- HEATMAP 2: discovery — top terms across all cell types -----------
make_discovery_heatmap <- function(direction, fill_col, n_terms = 25) {
  # rank terms by best p across cell types, but require sig in >=2 cell types
  sub <- combined |> filter(dir == direction)
  term_stats <- sub |> group_by(term_name, source) |>
    summarise(min_p   = min(p_value),
              n_sig05 = sum(p_value < 0.05),
              n_sig10 = sum(p_value < 0.10),
              tsize   = first(term_size),
              .groups = "drop") |>
    filter(n_sig10 >= 2, tsize >= 10, tsize <= 800) |>
    arrange(min_p) |> head(n_terms)
  if (nrow(term_stats) == 0) return(NULL)
  d <- sub |> filter(term_name %in% term_stats$term_name) |>
    mutate(neglog10p = -log10(p_value))
  d$cell_type <- order_ct(d$cell_type)
  d$term_name <- factor(d$term_name, levels = rev(term_stats$term_name))
  d$label_y   <- factor(sprintf("%s | %s", d$source, d$term_name),
                        levels = rev(sprintf("%s | %s", term_stats$source, term_stats$term_name)))
  ggplot(d, aes(x = cell_type, y = label_y, fill = neglog10p)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(p_value < 0.05, "*",
                          ifelse(p_value < 0.1,  "+", ""))),
              size = 4, vjust = 0.7) +
    scale_fill_gradient(low = "grey95", high = fill_col,
                        name = expression(-log[10]~"FDR"),
                        limits = c(0, NA)) +
    labs(x = NULL, y = NULL,
         title = sprintf("Top %d cross-cell-type %s-regulated terms (sig in >= 2 cell types)",
                          nrow(term_stats), direction)) +
    theme_cowplot(font_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 9),
          plot.title  = element_text(size = 12, face = "bold"),
          panel.grid  = element_blank(),
          axis.line   = element_blank(),
          axis.ticks  = element_blank())
}

p_disc_down <- make_discovery_heatmap("down", "#0072B2", 25)
p_disc_up   <- make_discovery_heatmap("up",   "#D55E00", 25)

ggsave(file.path(OUT_DIR, "heatmap_discovery_down.png"),
       p_disc_down, width = 13, height = 10, dpi = 220, bg = "white")
ggsave(file.path(OUT_DIR, "heatmap_discovery_up.png"),
       p_disc_up,   width = 13, height = 10, dpi = 220, bg = "white")

# ---------- Print a quick text summary --------------------------------------
cat("\n===== Sst presynaptic/vesicle terms across cell types =====\n")
foc_d <- build_heatmap_data("down", sst_terms$down) |>
  mutate(sig = case_when(p_value < 0.05 ~ "**", p_value < 0.1 ~ "*", TRUE ~ "")) |>
  filter(p_value < 0.1) |>
  arrange(term_name, p_value)
print(foc_d, n = 200)

cat("\n===== Sst 'up' BMP/transport/adhesion terms across cell types =====\n")
foc_u <- build_heatmap_data("up", sst_terms$up) |>
  mutate(sig = case_when(p_value < 0.05 ~ "**", p_value < 0.1 ~ "*", TRUE ~ "")) |>
  filter(p_value < 0.2) |>
  arrange(term_name, p_value)
print(foc_u, n = 200)

cat("\nDone. Outputs in", OUT_DIR, "\n")
