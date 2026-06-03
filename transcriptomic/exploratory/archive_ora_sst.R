#!/usr/bin/env Rscript
# gprofiler2 enrichment for Sst SCZ DE genes
# - 4 gene sets: all-up, all-down, Sst-exclusive-up, Sst-exclusive-down
# - Custom background = all 12,492 genes tested in Sst (proper conditioning)
# - Multi-source: GO:BP, GO:MF, GO:CC, Reactome, KEGG, WikiPathways, HPA, TF, HP

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(cowplot); library(gprofiler2)
})

INPUT     <- "data/DE_genes_all_cells_scz.csv"
CELL_TYPE <- "Sst"
PADJ_THR  <- 0.1
SOURCES   <- c("GO:BP","GO:MF","GO:CC","REAC","KEGG","WP","HPA","TF","HP")
OUT_DIR   <- "results/exploratory/archive_ora_sst"
dir.create(OUT_DIR, showWarnings = FALSE)

COL_UP   <- "#D55E00"
COL_DOWN <- "#0072B2"

# ---- Build gene lists ------------------------------------------------------
df  <- read_csv(INPUT, show_col_types = FALSE)
sst <- df |> filter(cell_type == CELL_TYPE)
bg  <- unique(sst$genes)
cat(sprintf("Background (all Sst-tested genes): %d\n", length(bg)))

sst_sig <- sst |> filter(padj < PADJ_THR) |>
  mutate(dir = if_else(estimate > 0, "up", "down"))

# Find which Sst-DE genes are NOT DE in any other cell type in the same direction
sig_all <- df |> filter(padj < PADJ_THR) |>
  mutate(dir = if_else(estimate > 0, "up", "down"))

sst_sig <- sst_sig |> rowwise() |>
  mutate(n_other = sum(sig_all$genes == genes &
                       sig_all$cell_type != CELL_TYPE &
                       sig_all$dir == dir)) |>
  ungroup()

lists <- list(
  up         = sst_sig |> filter(dir == "up")   |> pull(genes),
  down       = sst_sig |> filter(dir == "down") |> pull(genes),
  up_only    = sst_sig |> filter(dir == "up",   n_other == 0) |> pull(genes),
  down_only  = sst_sig |> filter(dir == "down", n_other == 0) |> pull(genes)
)
for (nm in names(lists)) cat(sprintf("  %-10s n = %d\n", nm, length(lists[[nm]])))

# ---- Run g:Profiler --------------------------------------------------------
run_one <- function(genes, label) {
  cat(sprintf("\n>>> g:Profiler [%s]  (%d genes)\n", label, length(genes)))
  res <- gost(
    query           = genes,
    organism        = "hsapiens",
    custom_bg       = bg,
    sources         = SOURCES,
    correction_method = "g_SCS",
    user_threshold  = 0.05,
    evcodes         = TRUE,
    significant     = TRUE
  )
  if (is.null(res) || nrow(res$result) == 0) {
    cat("  no significant terms\n"); return(NULL)
  }
  out <- res$result |>
    select(source, term_name, term_id, p_value,
           term_size, query_size, intersection_size, intersection) |>
    arrange(p_value) |>
    mutate(query_set = label)
  cat(sprintf("  %d sig terms (top: %s, p=%.2g)\n",
              nrow(out), out$term_name[1], out$p_value[1]))
  out
}

results <- imap(lists, ~ run_one(.x, .y))
results <- compact(results)

# Write per-list CSVs + combined
for (nm in names(results)) {
  write_csv(results[[nm]], file.path(OUT_DIR, sprintf("enrich_%s.csv", nm)))
}
all_res <- bind_rows(results)
write_csv(all_res, file.path(OUT_DIR, "enrich_all.csv"))

# ---- Summary plot: top 10 per direction ------------------------------------
plot_top <- function(res, label, fill) {
  top <- res |> arrange(p_value) |> head(15) |>
    mutate(term_name = ifelse(nchar(term_name) > 55,
                              paste0(substr(term_name, 1, 52), "..."),
                              term_name),
           term_label = sprintf("%s | %s", source, term_name),
           term_label = factor(term_label, levels = rev(term_label)),
           neglog10p  = -log10(p_value))
  ggplot(top, aes(x = neglog10p, y = term_label)) +
    geom_col(fill = fill, alpha = 0.85, width = 0.75) +
    geom_text(aes(label = sprintf("%d/%d", intersection_size, term_size)),
              hjust = -0.1, size = 3.4) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = expression(-log[10]~"(adj p)"), y = NULL,
         subtitle = sprintf("%s   (n = %d genes)", label,
                            length(lists[[label]]))) +
    theme_cowplot(font_size = 12) +
    theme(plot.subtitle = element_text(size = 13, face = "bold"),
          axis.text.y   = element_text(size = 10))
}

panels <- list()
if (!is.null(results$up))        panels$up        <- plot_top(results$up,        "up",        COL_UP)
if (!is.null(results$down))      panels$down      <- plot_top(results$down,      "down",      COL_DOWN)
if (!is.null(results$up_only))   panels$up_only   <- plot_top(results$up_only,   "up_only",   COL_UP)
if (!is.null(results$down_only)) panels$down_only <- plot_top(results$down_only, "down_only", COL_DOWN)

grid <- plot_grid(plotlist = panels, ncol = 2, align = "v")
ggsave(file.path(OUT_DIR, "enrichment_summary.png"),
       grid, width = 16, height = 12, dpi = 200, bg = "white")
ggsave(file.path(OUT_DIR, "enrichment_summary.pdf"),
       grid, width = 16, height = 12, bg = "white")
cat(sprintf("\nWrote summary figure and %d CSVs into %s/\n",
            length(results), OUT_DIR))

# ---- Source-level summary (how many sig terms per source per list) --------
src_summary <- all_res |> count(query_set, source) |>
  pivot_wider(names_from = source, values_from = n, values_fill = 0)
cat("\n=== Significant terms per source per list ===\n")
print(src_summary)
write_csv(src_summary, file.path(OUT_DIR, "source_summary.csv"))
