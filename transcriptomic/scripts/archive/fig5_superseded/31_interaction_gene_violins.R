#!/usr/bin/env Rscript
# Violin plots for the FDR-significant diagnosis x stratum interaction genes
# (S11): control vs SCZ by stratum, from the saved per-donor stratum
# pseudobulks. For an interaction gene the case-control gap should CHANGE
# across strata. Annotations: per-stratum meta z (top) and the depleted-vs-
# non-depleted interaction z (facet strip).
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(arrow); library(ggplot2); library(cowplot)
})
POUT <- "transcriptomic/results/sst_strata_gsea/pseudobulk"
DOUT <- file.path(POUT, "donor_stratum")
FS <- 13
GENES <- c("CELF2", "GRIK4", "L3MBTL3", "HMGXB4", "QKI", "GRIN3A", "FHOD3", "FAM151B")
S_LEV <- c("depleted", "intermediate", "non_depleted")
S_LAB <- c("Depl.", "Interm.", "Non-depl.")

cohorts <- str_remove(basename(Sys.glob(file.path(DOUT, "*_donor_stratum_meta.csv"))),
                      "_donor_stratum_meta\\.csv$")
vals <- map_dfr(cohorts, function(co) {
  pb <- as.data.frame(read_parquet(file.path(DOUT, paste0(co, "_donor_stratum_counts.parquet"))))
  rn <- pb$gene; pb$gene <- NULL
  M <- vapply(pb, as.numeric, numeric(nrow(pb))); rownames(M) <- rn
  m <- read_csv(file.path(DOUT, paste0(co, "_donor_stratum_meta.csv")),
                col_types = cols(donor = col_character(), .default = col_guess()))
  cpm <- t(t(M[, paste(m$donor, m$stratum, sep = "|")]) /
             colSums(M[, paste(m$donor, m$stratum, sep = "|")])) * 1e6
  gg <- intersect(GENES, rownames(cpm))
  as_tibble(t(log1p(cpm[gg, , drop = FALSE]))) |>
    mutate(donor = m$donor, stratum = m$stratum, diagnosis = m$diagnosis, cohort = co) |>
    pivot_longer(all_of(gg), names_to = "gene", values_to = "logcpm")
})
d <- vals |>
  group_by(cohort, stratum, gene) |>
  mutate(centered = logcpm - mean(logcpm)) |> ungroup() |>
  group_by(gene) |>
  mutate(centered = pmax(pmin(centered, quantile(centered, 0.995)),
                         quantile(centered, 0.005))) |> ungroup() |>
  mutate(stratum = factor(stratum, S_LEV, S_LAB),
         diagnosis = factor(diagnosis, c("Control", "SCZ")))

im <- read_csv(file.path(POUT, "interaction_meta.csv"), show_col_types = FALSE) |>
  filter(coef == "dxSCZ:stratumdepleted", gene %in% GENES)
strip_lab <- setNames(sprintf("%s  (int z = %.1f, padj = %.2g)",
                              im$gene, im$zval, im$padj), im$gene)
sm <- read_csv(file.path(POUT, "stratum_meta_de.csv"), show_col_types = FALSE) |>
  filter(gene %in% GENES, stratum != "all_sst") |>
  mutate(stratum = factor(stratum, S_LEV, S_LAB), lab = sprintf("z=%.1f", zval))
d$gene <- factor(d$gene, GENES); sm$gene <- factor(sm$gene, GENES)

p <- ggplot(d, aes(stratum, centered, fill = diagnosis, colour = diagnosis)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_violin(position = position_dodge(width = 0.75), width = 0.7,
              alpha = 0.25, linewidth = 0.45, scale = "width") +
  geom_point(position = position_jitterdodge(jitter.width = 0.09,
                                             dodge.width = 0.75),
             size = 0.45, alpha = 0.35, show.legend = FALSE) +
  stat_summary(fun = median, geom = "crossbar", width = 0.55, linewidth = 0.5,
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  geom_text(data = sm, aes(x = stratum, y = Inf, label = lab),
            inherit.aes = FALSE, vjust = 1.35, size = 3.4, colour = "grey25") +
  facet_wrap(~gene, nrow = 2, scales = "free_y",
             labeller = labeller(gene = strip_lab)) +
  scale_fill_manual(values = c(Control = "#0072B2", SCZ = "#D55E00"), name = NULL) +
  scale_colour_manual(values = c(Control = "#0072B2", SCZ = "#D55E00"), name = NULL) +
  labs(x = NULL, y = "Donor pseudobulk log1p(CPM), centered within cohort × stratum",
       caption = paste("FDR-significant diagnosis × stratum interaction genes (S11);",
                       "strip: depleted-vs-non-depleted interaction z and padj;",
                       "in-panel z: per-stratum meta-analytic diagnosis effect;",
                       "display clipped at 0.5-99.5 percentile per gene")) +
  theme_cowplot(font_size = FS) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 25, hjust = 1),
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "italic", size = FS - 1))
ggsave(file.path(POUT, "figS_interaction_gene_violins.png"), p,
       width = 15.5, height = 8.6, dpi = 200, bg = "white")
cat("wrote figS_interaction_gene_violins.png\n")
