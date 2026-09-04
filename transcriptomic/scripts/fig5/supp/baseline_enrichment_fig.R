#!/usr/bin/env Rscript
# Two-panel mockup for the baseline-signature enrichment (Fig. 4g companion):
# a) top ORA terms per direction; b) baseline log2FC of the Fig. 5 disease-state
# modules. Direction colors reuse the Fig. 5a stratum semantics (purple =
# depleted side, green = spared side); module colors reuse Fig. 5 FAM_COLS.
source("transcriptomic/scripts/fig5/_common.R")
suppressPackageStartupMessages({ library(forcats); library(ggplot2); library(cowplot) })
OUT <- file.path(P$supp, "baseline_marker_enrichment")
# depleted/spared reuse the Figure 5a stratum colours; modules reuse FAM_COLS
DIR_COLS <- c(`Higher in depleted` = unname(STRAT_COLS[["depleted"]]),
              `Higher in spared`   = unname(STRAT_COLS[["non_depleted"]]))
MOD_COLS <- c(Synaptic = unname(FAM_COLS[["synaptic"]]),
              Translation = unname(FAM_COLS[["translation"]]),
              `OxPhos/mito` = unname(FAM_COLS[["oxphos"]]))

ora <- read_csv(file.path(OUT, "ora_580_baseline_signature.csv"), show_col_types = FALSE) |>
  mutate(direction = ifelse(direction == "higher_in_depleted",
                            "Higher in depleted", "Higher in spared"))
pick <- ora |> group_by(direction) |> slice_min(padj, n = 8, with_ties = FALSE) |> ungroup() |>
  mutate(lab = pathway |> str_remove("^GOBP_|^GOCC_|^GOMF_|^REACTOME_|^HALLMARK_") |>
           str_replace_all("_", " ") |> str_to_sentence() |>
           str_replace("^Srp dependent.*", "SRP-dependent cotransl. targeting") |>
           str_replace("^Response of eif2ak4.*", "GCN2 response (ribosomal-gene set)") |>
           str_replace("Gpcr|G protein coupled receptor", "GPCR") |> str_wrap(34),
         enr = overlap / expected,
         stars = case_when(padj < 1e-15 ~ "***", padj < 1e-9 ~ "**", padj < 0.05 ~ "*", TRUE ~ ""))
p_a <- ggplot(pick, aes(enr, fct_reorder(lab, enr), fill = direction)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf(" %d/%d %s", overlap, size, stars)), hjust = 0, size = 3.4, colour = "grey20") +
  facet_wrap(~direction, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = DIR_COLS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(x = "Fold enrichment over expressed-gene background", y = NULL) +
  theme_cowplot(font_size = FS) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = FS - 4, lineheight = 0.85))

v <- read_csv("genetics/results/figures/r_panels/panel_volcano_vulnerable_vs_notdepleted.csv",
              show_col_types = FALSE) |>
  filter(pct_vulnerable >= 0.10 | pct_not_depleted >= 0.10)
g5 <- read_csv(file.path(P$pb, "gsea_all_signatures.csv"), show_col_types = FALSE)
le_union <- function(ps) g5 |> filter(pathway %in% ps, signature == "depleted") |>
  pull(leadingEdge) |> strsplit("|", fixed = TRUE) |> unlist() |> unique()
MODS <- list(
  Translation = le_union(c("REACTOME_TRANSLATION","GOCC_RIBOSOMAL_SUBUNIT",
                           "GOBP_CYTOPLASMIC_TRANSLATION","REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION")),
  `OxPhos/mito` = le_union(c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","GOBP_OXIDATIVE_PHOSPHORYLATION",
                             "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT")),
  Synaptic = le_union(c("GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING","GOCC_SYNAPTIC_MEMBRANE",
                        "GOBP_NEUROTRANSMITTER_SECRETION","REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES")))
md <- bind_rows(lapply(names(MODS), function(m)
  v |> filter(gene %in% MODS[[m]]) |> transmute(module = m, log2FC)))
stats <- md |> group_by(module) |> summarise(
  med = median(log2FC),
  p = wilcox.test(log2FC, v$log2FC[!v$gene %in% MODS[[cur_group()$module]]])$p.value)
md$module <- factor(md$module, names(MOD_COLS)); stats$module <- factor(stats$module, names(MOD_COLS))
p_b <- ggplot(md, aes(module, log2FC, fill = module)) +
  geom_hline(yintercept = median(v$log2FC), colour = "grey55", linetype = "dashed") +
  geom_violin(width = 0.8, alpha = 0.55, colour = NA, scale = "width") +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.5) +
  geom_text(data = stats, aes(y = 0.62, label = sprintf("med %.2f\np %s", med,
            format(signif(p, 1), scientific = TRUE))), size = 3.6, colour = "grey20", vjust = 0) +
  scale_fill_manual(values = MOD_COLS, guide = "none") +
  coord_cartesian(ylim = c(-0.55, 0.78)) +
  labs(x = NULL, y = "Baseline log2FC, depleted vs spared\n(neurotypical SEA-AD DLPFC)") +
  theme_cowplot(font_size = FS) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
annot <- ggdraw(p_b) + draw_label("dashed: background median", x = 0.98, y = 0.97,
                                  hjust = 1, size = FS - 4, colour = "grey40")
fig <- plot_grid(p_a, annot, nrow = 1, rel_widths = c(1.45, 1), labels = c("a", "b"), label_size = 20)
ggsave(file.path(OUT, "figS_baseline_marker_enrichment.png"), fig, width = 13.5, height = 7.6,
       dpi = 200, bg = "white")
cat("wrote figS_baseline_marker_enrichment.png\n")
