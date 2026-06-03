#!/usr/bin/env Rscript
# Follow-up: look at top-ranked terms regardless of g_SCS significance,
# and try FDR (less conservative) correction.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr)
  library(gprofiler2)
})

INPUT     <- "data/DE_genes_all_cells_scz.csv"
CELL_TYPE <- "Sst"
PADJ_THR  <- 0.1
SOURCES   <- c("GO:BP","GO:MF","GO:CC","REAC","KEGG","WP")
OUT_DIR   <- "results/exploratory/archive_ora_sst"

df  <- read_csv(INPUT, show_col_types = FALSE)
sst <- df |> filter(cell_type == CELL_TYPE)
bg  <- unique(sst$genes)
sig <- sst |> filter(padj < PADJ_THR) |>
  mutate(dir = if_else(estimate > 0, "up", "down"))

lists <- list(
  up   = sig |> filter(dir == "up")   |> pull(genes),
  down = sig |> filter(dir == "down") |> pull(genes)
)

run <- function(genes, label, method) {
  res <- gost(query = genes, organism = "hsapiens",
              custom_bg = bg, sources = SOURCES,
              correction_method = method,
              user_threshold = 1,
              significant = FALSE,
              evcodes = TRUE)
  if (is.null(res)) return(tibble())
  res$result |> as_tibble() |>
    select(source, term_name, p_value, term_size, intersection_size, intersection) |>
    arrange(p_value) |> head(20) |>
    mutate(query_set = label, correction = method)
}

cat("=== Top 20 terms per direction, g_SCS correction (gprofiler default) ===\n")
for (nm in names(lists)) {
  r <- run(lists[[nm]], nm, "g_SCS")
  cat(sprintf("\n--- %s (n=%d) ---\n", nm, length(lists[[nm]])))
  print(r |> select(source, term_name, p_value, term_size, intersection_size), n = 20)
}

cat("\n\n=== Top 20 terms per direction, FDR correction ===\n")
for (nm in names(lists)) {
  r <- run(lists[[nm]], nm, "fdr")
  cat(sprintf("\n--- %s (n=%d) ---\n", nm, length(lists[[nm]])))
  print(r |> select(source, term_name, p_value, term_size, intersection_size), n = 20)
  write_csv(r, file.path(OUT_DIR, sprintf("top20_fdr_%s.csv", nm)))
}
