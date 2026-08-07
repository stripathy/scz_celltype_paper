#!/usr/bin/env Rscript
# Final enrichment figure: FDR-corrected, term-size-filtered for interpretability.
# Two panels: Sst-down and Sst-up, top terms per source.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(cowplot); library(gprofiler2)
})

INPUT     <- "data/DE_genes_all_cells_scz.csv"
CELL_TYPE <- "Sst"
PADJ_THR  <- 0.1
SOURCES   <- c("GO:BP","GO:MF","GO:CC","REAC","KEGG","WP")
OUT_DIR   <- "results/exploratory/archive_ora_sst"

COL_UP   <- "#D55E00"
COL_DOWN <- "#0072B2"

df  <- read_csv(INPUT, show_col_types = FALSE)
sst <- df |> filter(cell_type == CELL_TYPE)
bg  <- unique(sst$genes)
sig <- sst |> filter(padj < PADJ_THR) |>
  mutate(dir = if_else(estimate > 0, "up", "down"))

lists <- list(
  up   = sig |> filter(dir == "up")   |> pull(genes),
  down = sig |> filter(dir == "down") |> pull(genes)
)

run_fdr <- function(genes) {
  res <- gost(query = genes, organism = "hsapiens",
              custom_bg = bg, sources = SOURCES,
              correction_method = "fdr",
              user_threshold = 1, significant = FALSE,
              evcodes = TRUE)
  if (is.null(res)) return(tibble())
  res$result |> as_tibble() |>
    select(source, term_name, p_value, term_size, intersection_size, intersection) |>
    arrange(p_value)
}

cat("Querying gprofiler2 (FDR)...\n")
enrich <- lapply(lists, run_fdr)

write_csv(bind_rows(enrich, .id = "query_set"),
          file.path(OUT_DIR, "enrich_fdr_full.csv"))

# Filter to interpretable terms:
#   - term_size between 5 and 800 (avoid huge nonspecific bins + tiny noise)
#   - p_value <= 0.25 (show suggestive + significant, marked)
filter_terms <- function(x, n = 15) {
  x |> filter(term_size >= 5, term_size <= 800, p_value <= 0.25) |>
    arrange(p_value) |> head(n)
}

down_show <- filter_terms(enrich$down, 15)
up_show   <- filter_terms(enrich$up,   15)

cat("\n--- DOWN: terms shown ---\n");  print(down_show |> select(source,term_name,p_value,term_size,intersection_size))
cat("\n--- UP: terms shown ---\n");    print(up_show   |> select(source,term_name,p_value,term_size,intersection_size))

# ---- Plot ------------------------------------------------------------------
make_panel <- function(d, label, fill_col) {
  d2 <- d |> mutate(
    term_short = ifelse(nchar(term_name) > 48,
                        paste0(substr(term_name, 1, 45), "..."), term_name),
    label_y    = sprintf("%s | %s", source, term_short),
    label_y    = factor(label_y, levels = rev(label_y)),
    neglog10p  = -log10(p_value),
    sig_mark   = case_when(p_value < 0.05 ~ "*", p_value < 0.1 ~ "+", TRUE ~ ""),
    label_text = sprintf("%d / %d  %s", intersection_size, term_size, sig_mark)
  )
  ggplot(d2, aes(x = neglog10p, y = label_y)) +
    geom_col(aes(alpha = -log10(p_value)),
             fill = fill_col, width = 0.72, show.legend = FALSE) +
    geom_text(aes(label = label_text), hjust = -0.08, size = 3.6) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    annotate("text", x = -log10(0.05), y = 0.6,
             label = "FDR 0.05", hjust = -0.05, vjust = 1,
             size = 3.4, colour = "grey30") +
    scale_alpha(range = c(0.45, 1)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
    labs(x = expression(-log[10]~"(FDR)"), y = NULL,
         subtitle = sprintf("Sst %s-regulated in SCZ  (n = %d genes)",
                            label, length(lists[[label]]))) +
    theme_cowplot(font_size = 13) +
    theme(plot.subtitle = element_text(size = 14, face = "bold",
                                       colour = fill_col),
          axis.text.y   = element_text(size = 11),
          plot.margin   = margin(8, 14, 8, 6))
}

p_down <- make_panel(down_show, "down", COL_DOWN)
p_up   <- make_panel(up_show,   "up",   COL_UP)

grid <- plot_grid(p_down, p_up, ncol = 1, align = "v",
                  rel_heights = c(1, 1))

ggsave(file.path(OUT_DIR, "enrichment_fdr_clean.png"),
       grid, width = 11, height = 11, dpi = 220, bg = "white")
ggsave(file.path(OUT_DIR, "enrichment_fdr_clean.pdf"),
       grid, width = 11, height = 11, bg = "white")
cat("\nSaved enrichment_fdr_clean.{png,pdf}\n")

# ---- Print actual gene-level overlaps for the top biologically meaningful terms ----
cat("\n=== Genes driving the top synaptic terms (DOWN) ===\n")
for (term in c("synapse","secretory vesicle","presynapse","synaptic vesicle","exocytic vesicle")) {
  row <- enrich$down |> filter(term_name == term) |> slice(1)
  if (nrow(row) > 0) cat(sprintf("\n[%s]  %d/%d  p=%.3g\n  %s\n",
                                  term, row$intersection_size, row$term_size,
                                  row$p_value, row$intersection))
}
cat("\n=== Genes driving the top BMP/adhesion terms (UP) ===\n")
for (term in c("I-SMAD binding","neurexin family protein binding","P-type calcium transporter activity")) {
  row <- enrich$up |> filter(term_name == term) |> slice(1)
  if (nrow(row) > 0) cat(sprintf("\n[%s]  %d/%d  p=%.3g\n  %s\n",
                                  term, row$intersection_size, row$term_size,
                                  row$p_value, row$intersection))
}
