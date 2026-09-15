#!/usr/bin/env Rscript
# Supplement | Volcanoes for the stratum DE and the interaction.
#
# WHY THIS EXISTS. The main figure asserts a graded burden from gene-set
# statistics. This shows the gene level honestly: significant genes are few, they
# do NOT track the gene-set gradient, and the module genes sit as a coherent
# low-amplitude cloud rather than as a handful of large effects. That is the
# evidence for "coordinated shifts across many genes", and it is also the
# disclosure a reviewer will ask for.
#
# a  SCZ vs control within each stratum
# b  the interaction contrasts (depleted vs non-depleted, intermediate vs
#    non-depleted, linear trend across strata)
# Points coloured by module membership, using the same palette as Figure S8.
source("transcriptomic/scripts/sst_strata/_common.R")
suppressPackageStartupMessages({ library(ggplot2); library(cowplot); library(ggrepel) })

g  <- read_csv(file.path(P$pb, "gsea_all_signatures.csv"), show_col_types = FALSE)
sm <- read_csv(file.path(P$pb, "stratum_meta_de.csv"),     show_col_types = FALSE)
im <- read_csv(file.path(P$pb, "interaction_meta.csv"),    show_col_types = FALSE)
mods <- modules(g)[MODULE_PRIORITY]

volcano <- function(d, title, xlab, fdr = 0.05) {
  d <- d |> mutate(module = factor(FAM_LABELS[module_of(gene, mods)], unname(FAM_LABELS)),
                   nlp = -log10(pval),
                   est = pmax(pmin(estimate, quantile(estimate, 0.999)),
                              quantile(estimate, 0.001)))
  lab <- d |> filter(padj < fdr)
  n_lab <- nrow(lab)
  if (n_lab > 25) lab <- lab |> slice_min(pval, n = 25)
  thr <- if (n_lab > 0) max(d$pval[d$padj < fdr]) else NA_real_
  p <- ggplot(d, aes(est, nlp, colour = module)) +
    geom_point(data = ~ filter(.x, module == FAM_LABELS[["other"]]), size = 0.5, alpha = 0.3) +
    geom_point(data = ~ filter(.x, module != FAM_LABELS[["other"]]), size = 1.0, alpha = 0.75) +
    geom_vline(xintercept = 0, colour = "grey85") +
    scale_colour_manual(values = setNames(unname(FAM_COLS), unname(FAM_LABELS)),
                        name = NULL, drop = FALSE) +
    labs(x = xlab, y = expression(-log[10] ~ "p"), title = title,
         subtitle = sprintf("%d genes at FDR < %.2f", n_lab, fdr)) +
    theme_cowplot(font_size = FS - 1) +
    theme(plot.title = element_text(size = FS), legend.position = "none",
          plot.subtitle = element_text(size = FS - 3, colour = "grey30"))
  if (is.finite(thr))
    p <- p + geom_hline(yintercept = -log10(thr), linetype = "dotted", colour = "grey50") +
      geom_point(data = lab, aes(fill = module), colour = "black", shape = 21,
                 size = 1.8, stroke = 0.6, show.legend = FALSE) +
      scale_fill_manual(values = setNames(unname(FAM_COLS), unname(FAM_LABELS)), guide = "none") +
      geom_text_repel(data = lab, aes(label = gene), colour = "black", size = 3.2,
                      fontface = "italic", max.overlaps = 30, seed = 3,
                      min.segment.length = 0.1)
  p
}

row1 <- plot_grid(plotlist = lapply(STRATA_LEVELS, function(st)
  volcano(sm |> filter(stratum == st),
          sprintf("%s stratum", c(depleted = "Depleted", intermediate = "Intermediate",
                                  non_depleted = "Non-depleted")[[st]]),
          "log2FC (SCZ vs control)")), nrow = 1)
IC <- c("dxSCZ:stratumdepleted", "dxSCZ:stratumintermediate", "dxSCZ:score")
IT <- c("Interaction: depleted vs non-depleted", "Interaction: intermediate vs non-depleted",
        "Interaction: linear trend across strata")
row2 <- plot_grid(plotlist = map2(IC, IT, ~ volcano(im |> filter(coef == .x), .y,
                                                    "interaction log2FC")), nrow = 1)
leg <- get_legend(
  ggplot(tibble(m = factor(unname(FAM_LABELS), unname(FAM_LABELS)), x = 1),
         aes(x, x, colour = m)) + geom_point(size = 3) +
    scale_colour_manual(values = setNames(unname(FAM_COLS), unname(FAM_LABELS)), name = NULL) +
    theme_cowplot(font_size = FS) + theme(legend.position = "bottom"))
fig <- plot_grid(row1, row2, leg, ncol = 1, rel_heights = c(1, 1, 0.08),
                 labels = c("a", "b", ""), label_size = 20)
ggsave(file.path(P$supp, "figS_volcanoes.png"), fig, width = 15, height = 10.2,
       dpi = 200, bg = "white")
say("wrote figS_volcanoes.png")

say("FDR<0.05 genes per contrast:")
print(bind_rows(sm |> transmute(contrast = paste0("stratum_", stratum), gene, padj),
                im |> transmute(contrast = coef, gene, padj)) |>
        filter(padj < 0.05) |> group_by(contrast) |>
        summarise(n = n(), genes = paste(head(gene[order(padj)], 10), collapse = ", ")),
      width = 200)
