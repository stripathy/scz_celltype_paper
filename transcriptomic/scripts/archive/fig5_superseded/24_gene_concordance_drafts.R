#!/usr/bin/env Rscript
# DRAFTS for a prospective Fig-5 last-row element: genes whose depleted-stratum
# DE is SHARED between SCZ and AD versus disorder-specific. Two layouts to
# iterate on before anything enters the figure:
#   draft 1: class-grouped mini-heatmap (rows = curated exemplars, cols = SCZ z / AD t)
#   draft 2: annotated scatter (all genes grey; exemplars colored by class)
# Classes (thresholds stated, applied to depleted-stratum effects):
#   shared down:      z_scz < -2.5 & t_ad < -1.5
#   SCZ-specific down z_scz < -3   & t_ad > -0.3
#   shared up:        z_scz >  2   & t_ad >  1.5
#   SCZ up / AD down: z_scz >  2   & t_ad < -1.5
# AD side = our A9 recalculation (linear CPS, severely affected excluded).
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(cowplot); library(ggrepel)
})
OUT  <- "transcriptomic/results/sst_strata_gsea"
AOUT <- file.path(OUT, "seaad_a9")
FS <- 14

sig <- read_csv(file.path(OUT, "stratum_gene_signatures.csv"), show_col_types = FALSE)
a9  <- read_csv(file.path(AOUT, "de_a9_strata.csv"), show_col_types = FALSE) |>
  filter(signature == "a9_depleted") |> select(gene, t_ad = t)
m <- sig |> filter(stratum == "depleted") |> select(gene, z_scz = z) |>
  inner_join(a9, by = "gene") |>
  mutate(class = case_when(
    z_scz < -2.5 & t_ad < -1.5 ~ "Shared down",
    z_scz < -3.0 & t_ad > -0.3 ~ "SCZ-specific down",
    z_scz >  2.0 & t_ad >  1.5 ~ "Shared up",
    z_scz >  2.0 & t_ad < -1.5 ~ "SCZ up / AD down",
    TRUE ~ "Other"))
write_csv(m |> filter(class != "Other") |> arrange(class, z_scz),
          file.path(AOUT, "gene_concordance_classes.csv"))
cat("class sizes:\n"); print(count(m, class))

CLS_LEVELS <- c("Shared down", "SCZ-specific down", "Shared up", "SCZ up / AD down")
CLS_COLS <- c("Shared down" = "#009E73", "SCZ-specific down" = "#CC79A7",
              "Shared up" = "#56B4E9", "SCZ up / AD down" = "#E69F00",
              Other = "grey85")
# curated exemplars: top by |z| within class, favoring annotated genes
EX <- list(
  "Shared down"       = c("SST", "VGF", "CBLN4", "UNC13A", "CIRBP", "RASGRF2"),
  "SCZ-specific down" = c("RPL36", "RPS21", "NDUFS8", "HDAC4", "MEG3", "NTM"),
  "Shared up"         = c("INSYN2B", "GLCCI1", "PTPRD", "CELF2"),
  "SCZ up / AD down"  = c("TCF7L2", "ACAT1", "SYTL2"))

# ---- draft 1: class-grouped mini-heatmap --------------------------------------
d1 <- imap_dfr <- NULL
d1 <- purrr::imap(EX, function(gg, cl) m |> filter(gene %in% gg) |> mutate(class = cl)) |>
  bind_rows() |>
  pivot_longer(c(z_scz, t_ad), names_to = "disorder", values_to = "eff") |>
  mutate(disorder = factor(disorder, c("z_scz", "t_ad"),
                           c("SCZ\n(z)", "AD\n(t, CPS)")),
         class = factor(class, CLS_LEVELS))
ord <- d1 |> filter(disorder == "SCZ\n(z)") |> arrange(class, eff) |> pull(gene)
d1$gene <- factor(d1$gene, levels = rev(unique(ord)))
p1 <- ggplot(d1, aes(disorder, gene, fill = eff)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", eff)), size = 3.9) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, limits = c(-4, 4), oob = scales::squish,
                       name = "Signed effect\n(z or t)") +
  facet_grid(class ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(class = label_wrap_gen(12))) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_cowplot(font_size = FS) +
  theme(axis.text.y = element_text(face = "italic"),
        strip.background = element_rect(fill = "grey92"),
        strip.text.y = element_text(angle = 270, size = FS - 2, face = "bold"),
        axis.ticks = element_blank(), axis.line = element_blank())

# ---- draft 2: annotated scatter -------------------------------------------------
hl <- m |> filter(gene %in% unlist(EX)) |> mutate(class = factor(class, CLS_LEVELS))
p2 <- ggplot(m, aes(z_scz, t_ad)) +
  geom_point(colour = "grey85", size = 0.5, alpha = 0.3) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(data = hl, aes(colour = class), size = 2.4) +
  geom_text_repel(data = hl, aes(label = gene, colour = class), size = 4,
                  fontface = "italic", seed = 7, max.overlaps = 30,
                  show.legend = FALSE) +
  scale_colour_manual(values = CLS_COLS, name = NULL) +
  labs(x = "Gene z, SCZ depleted stratum (case vs control)",
       y = "Gene t, AD depleted stratum (along CPS)") +
  theme_cowplot(font_size = FS) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 3)))

fig <- plot_grid(p1, p2, nrow = 1, rel_widths = c(0.34, 0.66),
                 labels = c("draft 1", "draft 2"), label_size = 15)
ggsave(file.path(AOUT, "figS_gene_concordance_DRAFTS.png"), fig,
       width = 14, height = 6.5, dpi = 200, bg = "white")
cat("wrote figS_gene_concordance_DRAFTS.png\n")

# ---- draft 3: per-donor exemplar genes from the Jens cohort --------------------
# (real pseudobulks; stand-in for the 7-cohort pseudobulks until Nicole's land)
suppressPackageStartupMessages(library(arrow))
JOUT <- file.path(OUT, "jens")
pbj <- as.data.frame(read_parquet(file.path(JOUT, "jens_stratum_pseudobulk.parquet")))
idx <- grep("__index_level_0__", colnames(pbj), value = TRUE)
rownames(pbj) <- pbj[[idx]]; pbj[[idx]] <- NULL; pbj <- as.matrix(pbj)
mj <- read_csv(file.path(JOUT, "jens_stratum_pseudobulk_meta.csv"),
               show_col_types = FALSE) |> filter(n_cells >= 10)
GX <- c("SST", "VGF", "RPL36", "INSYN2B")
cpm <- t(t(pbj[, mj$key]) / colSums(pbj[, mj$key]) * 1e6)
d3 <- as_tibble(t(log1p(cpm[GX, , drop = FALSE])), rownames = "key") |>
  left_join(mj |> select(key, stratum, diagnosis), by = "key") |>
  pivot_longer(all_of(GX), names_to = "gene", values_to = "expr") |>
  mutate(gene = factor(gene, GX),
         stratum = factor(stratum, c("depleted", "intermediate", "non_depleted"),
                          c("Depl.", "Interm.", "Non-depl.")),
         diagnosis = factor(diagnosis, c("CTRL", "SCZ")))
p3 <- ggplot(d3, aes(stratum, expr, colour = diagnosis)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12,
                                             dodge.width = 0.65),
             size = 0.9, alpha = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.5,
               position = position_dodge(width = 0.65)) +
  facet_wrap(~gene, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c(CTRL = "#0072B2", SCZ = "#D55E00"), name = NULL) +
  labs(x = NULL, y = "Donor pseudobulk, log1p(CPM)",
       caption = "Jens cohort (83 donors; independent of the 7-cohort meta-analysis)") +
  theme_cowplot(font_size = FS) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "italic", size = FS + 1))
ggsave(file.path(AOUT, "figS_gene_exemplars_jens_DRAFT.png"), p3,
       width = 13, height = 4.6, dpi = 200, bg = "white")
cat("wrote figS_gene_exemplars_jens_DRAFT.png\n")
