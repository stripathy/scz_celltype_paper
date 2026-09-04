#!/usr/bin/env Rscript
# Violin exemplars from the REAL 7-cohort per-donor stratum pseudobulks:
# control vs SCZ, stratified by depletion group, for shared (SST, VGF) and
# divergent (RPL36, NDUFS8) genes. Values are log1p(CPM) centered within
# cohort x stratum (removes cohort baselines; what remains is the within-cohort
# case-control contrast the meta-analysis tests). Also saves the per-donor
# stratum pseudobulk matrices for reuse (sensitivities S2/S9/S10, future panels).
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(arrow); library(ggplot2); library(cowplot)
})
EXP  <- "transcriptomic/data/stratum_pseudobulks_export"
POUT <- "transcriptomic/results/sst_strata_gsea/pseudobulk"
DOUT <- file.path(POUT, "donor_stratum")
dir.create(DOUT, showWarnings = FALSE, recursive = TRUE)
FS <- 14
MIN_CELLS <- 10
GENES <- c("SST", "VGF", "RPL36", "NDUFS8")
STRATA <- list(
  depleted     = c("Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25"),
  intermediate = c("Sst_9", "Sst_11", "Sst_12", "Sst_13", "Sst_19", "Sst_23"),
  non_depleted = c("Sst_1", "Sst_4", "Sst_5", "Sst_7", "Sst_10"))
st_of <- setNames(rep(names(STRATA), lengths(STRATA)), unlist(STRATA, use.names = FALSE))

cohorts <- str_remove(basename(Sys.glob(file.path(EXP, "*_groups.csv"))), "_groups\\.csv$")
vals <- map_dfr(cohorts, function(co) {
  gr <- read_csv(file.path(EXP, paste0(co, "_groups.csv")),
                 col_types = cols(donor = col_character(), .default = col_guess())) |>
    filter(supertype %in% names(st_of), n_cells > 0) |>
    mutate(stratum = st_of[supertype])
  f_pq <- file.path(EXP, paste0(co, "_pseudobulk_counts.parquet"))
  idcols <- intersect(c("gene", "gene_symbol"), names(open_dataset(f_pq)))
  pb <- as.data.frame(read_parquet(f_pq, col_select = all_of(c(idcols, gr$key))))
  sym <- if ("gene_symbol" %in% names(pb)) coalesce(pb$gene_symbol, NA_character_) else pb$gene
  looks_sym <- !grepl("^ENSG", pb$gene)
  sym[is.na(sym) & looks_sym] <- pb$gene[is.na(sym) & looks_sym]
  keep <- !is.na(sym)
  M <- vapply(pb[keep, gr$key, drop = FALSE], as.numeric, numeric(sum(keep)))
  M <- rowsum(M, group = sym[keep])
  ds <- gr |> group_by(donor, stratum) |>
    summarise(n_cells = sum(n_cells), diagnosis = first(diagnosis), .groups = "drop") |>
    filter(n_cells >= MIN_CELLS)
  agg <- sapply(seq_len(nrow(ds)), function(i)
    rowSums(M[, gr$key[gr$donor == ds$donor[i] & gr$stratum == ds$stratum[i]],
              drop = FALSE]))
  colnames(agg) <- paste(ds$donor, ds$stratum, sep = "|")
  # persist for reuse
  write_parquet(as.data.frame(agg) |> mutate(gene = rownames(agg), .before = 1),
                file.path(DOUT, paste0(co, "_donor_stratum_counts.parquet")))
  write_csv(ds |> mutate(cohort = co), file.path(DOUT, paste0(co, "_donor_stratum_meta.csv")))
  cpm <- t(t(agg) / colSums(agg)) * 1e6
  gg <- intersect(GENES, rownames(cpm))
  as_tibble(t(log1p(cpm[gg, , drop = FALSE])), rownames = "key") |>
    separate(key, c("donor", "stratum"), sep = "\\|") |>
    left_join(ds, by = c("donor", "stratum")) |>
    pivot_longer(all_of(gg), names_to = "gene", values_to = "logcpm") |>
    mutate(cohort = co)
})
cat(sprintf("values: %d donor-stratum-gene rows, %d cohorts\n",
            nrow(vals), n_distinct(vals$cohort)))

d <- vals |>
  group_by(cohort, stratum, gene) |>
  mutate(centered = logcpm - mean(logcpm)) |> ungroup() |>
  group_by(gene) |>
  mutate(centered = pmax(pmin(centered, quantile(centered, 0.995)),
                         quantile(centered, 0.005))) |> ungroup() |>
  mutate(stratum = factor(stratum, names(STRATA),
                          c("Depleted", "Intermediate", "Non-depleted")),
         diagnosis = factor(diagnosis, c("Control", "SCZ")),
         gene = factor(gene, GENES))

meta <- read_csv(file.path(POUT, "stratum_meta_de.csv"), show_col_types = FALSE) |>
  filter(gene %in% GENES, stratum != "all_sst") |>
  mutate(stratum = factor(stratum, names(STRATA),
                          c("Depleted", "Intermediate", "Non-depleted")),
         gene = factor(gene, GENES),
         lab = sprintf("z = %.1f", zval))

p <- ggplot(d, aes(stratum, centered, fill = diagnosis, colour = diagnosis)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_violin(position = position_dodge(width = 0.75), width = 0.7,
              alpha = 0.25, linewidth = 0.5, scale = "width") +
  geom_point(position = position_jitterdodge(jitter.width = 0.10,
                                             dodge.width = 0.75),
             size = 0.55, alpha = 0.4, show.legend = FALSE) +
  stat_summary(fun = median, geom = "crossbar", width = 0.55, linewidth = 0.55,
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  geom_text(data = meta, aes(x = stratum, y = Inf, label = lab),
            inherit.aes = FALSE, vjust = 1.4, size = 3.9, colour = "grey25") +
  facet_wrap(~gene, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#0072B2", SCZ = "#D55E00"), name = NULL) +
  scale_colour_manual(values = c(Control = "#0072B2", SCZ = "#D55E00"), name = NULL) +
  labs(x = NULL, y = "Donor pseudobulk log1p(CPM),\ncentered within cohort × stratum",
       caption = paste("7-cohort per-donor stratum pseudobulks (469 donors);",
                       "z = stratum meta-analytic z (SCZ vs control);",
                       "display clipped at the 0.5-99.5 percentile per gene")) +
  theme_cowplot(font_size = FS) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 25, hjust = 1),
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "italic", size = FS + 1))
ggsave(file.path(POUT, "figS_stratum_violin_exemplars.png"), p,
       width = 14.5, height = 5.4, dpi = 200, bg = "white")
cat("wrote figS_stratum_violin_exemplars.png\n")
