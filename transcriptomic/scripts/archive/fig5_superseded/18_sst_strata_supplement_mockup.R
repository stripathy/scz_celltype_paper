#!/usr/bin/env Rscript
# Mock-up of the Sst-strata pathway supplement (Etienne-review response, C2).
#
# Layout (2 rows x 3 columns, landscape):
#   a  strata definition: per-supertype crumblr depletion vs median depth,
#      colored by stratum, layer gutter, FDR<0.20 grouping boundary
#   b  pathway burden per group (padj<0.05), pooled subclass as reference
#   c  pairwise GSEA NES scatters (depleted vs non-depleted / vs intermediate),
#      colored by module family
#   d  pathway NES heatmap, shared vs graded blocks, subclass reference column
#   e  leading-edge exemplar genes per block (+ SST pinned for reference)
#
# Reads cached outputs of 17_sst_strata_gsea.R — no fgsea rerun.
# NOTE: gene ranks come from the IVW aggregation (script 17 caveat); swap in
# Nicole's stratum pseudobulk meta-DE before publication.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(cowplot); library(ggrepel)
})
OUT      <- "transcriptomic/results/sst_strata_gsea"
SUBCLASS <- "shared/snrnaseq_de/DE_genes_all_cells_scz.csv"
FIG4_PANELS <- "genetics/results/figures/r_panels"   # Fig 4 supertype colors
FS  <- 14
# Figure 2 conventions (transcriptomic/scripts/01_volcano_per_celltype.R)
COL_UP   <- "#D55E00"; COL_DOWN <- "#0072B2"
# Figure 4 conventions (genetics/scripts/figures/fig4_style.R): depletion is
# encoded by point OUTLINE + STROKE over per-supertype SEA-AD fills
DEPLETION_OUTLINE <- c(`Depleted (FDR < 0.20)` = "black", `Not depleted` = "grey60")
DEPLETION_STROKE  <- c(`Depleted (FDR < 0.20)` = 0.7, `Not depleted` = 0.2)

g      <- read_csv(file.path(OUT, "gsea_all_signatures.csv"),   show_col_types = FALSE)
sig    <- read_csv(file.path(OUT, "stratum_gene_signatures.csv"), show_col_types = FALSE)
strata <- read_csv(file.path(OUT, "strata_definition.csv"),     show_col_types = FALSE)

STRATA_LEVELS <- c("depleted", "intermediate", "non_depleted")
STRATA_LABELS <- c("Depleted (5)", "Intermediate (6)", "Non-depleted (5)")
SIG_LEVELS <- c("Sst_subclass", STRATA_LEVELS)          # a & d carry the pooled
SIG_LABELS <- c("All Sst (subclass)", STRATA_LABELS)    # subclass as reference
SIG_SHORT  <- c("All Sst", "Depleted", "Intermediate", "Non-depleted")  # b, d, e
# strata palette: purple/green, deliberately outside the blue-orange (down/up)
# and magenta/teal (module) families used in b-e
STRAT_COLS <- c("Depleted (5)" = "#762A83", "Intermediate (6)" = "#C2A5CF",
                "Non-depleted (5)" = "#1B7837")
# module palette: maximize contrast between Translation and OxPhos/mito
FAM_COLS   <- c(Synaptic = "#7FA6CC", Translation = "#AA3377",
                `OxPhos/mito` = "#EE7733", Deubiquitination = "#009988",
                Other = "grey85")
# Figure 2 two-tier significance palette (09_composite_figure.R)
UP_DARK <- "#D55E00"; UP_LIGHT <- "#F2B58C"
DOWN_DARK <- "#0072B2"; DOWN_LIGHT <- "#9FCAE6"
nes_fill <- function(lims) scale_fill_gradient2(
  low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
  limits = lims, oob = scales::squish, name = "NES (SCZ vs control)")

# ---- a | strata definition scatter -------------------------------------------
# Depth vs crumblr abundance change for the 16 Sst supertypes, colored by
# stratum, with the layer gutter and the FDR<0.20 grouping boundary. (Fig 3f
# plots the same two variables; this panel adds the strata and the layers.)
p_def <- local({
  d <- strata |> mutate(stratum = factor(stratum, STRATA_LEVELS, STRATA_LABELS))
  # boundary between depleted (5 most negative, crumblr FDR<0.20) and the rest
  cut_fdr <- mean(c(max(d$estimate[d$stratum == STRATA_LABELS[1]]),
                    min(d$estimate[d$stratum == STRATA_LABELS[2]])))
  # cortical layer boundaries manually proposed from the Xenium data
  lb <- read_csv("spatial/output/depth_proportions/proposed_layer_boundaries.csv",
                 show_col_types = FALSE)$xenium
  layers <- tibble(name = c("L1", "L2/3", "L4", "L5", "L6"),
                   mid = (c(0, lb[-length(lb)]) + lb) / 2)
  XR <- 0.155                       # right edge of the annotation gutter
  ggplot(d, aes(estimate, depth_xenium, colour = stratum)) +
    annotate("segment", x = XR - 0.012, xend = XR,
             y = c(0, lb), yend = c(0, lb), colour = "grey55", linewidth = 0.5) +
    annotate("text", x = XR, y = layers$mid, label = layers$name,
             hjust = 1, size = 3.9, colour = "grey30") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = cut_fdr, linetype = "dashed", colour = "grey60") +
    annotate("text", x = cut_fdr - 0.005, y = 0.985, hjust = 1, size = 4.2,
             colour = "grey30", label = "crumblr FDR < 0.20") +
    # strata named in-panel (they partition the x axis), replacing the legend
    annotate("text", x = c(-0.235, -0.092, 0.045), y = 0.035,
             label = STRATA_LABELS, colour = unname(STRAT_COLS),
             fontface = "bold", size = 4.3) +
    geom_errorbarh(aes(xmin = estimate - se, xmax = estimate + se),
                   height = 0, linewidth = 0.5, alpha = 0.6) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = CellType), size = 4, seed = 1,
                    show.legend = FALSE, max.overlaps = 20) +
    scale_y_reverse() +
    scale_colour_manual(values = STRAT_COLS, name = NULL) +
    labs(x = "Abundance change in SCZ (crumblr estimate ± SE)",
         y = "Median cortical depth (0 = pia)") +
    theme_cowplot(font_size = FS) +
    theme(legend.position = "none")
})

# ---- b | pathway burden, Figure 2 two-tier significance ----------------------
p_burden <- local({
  d <- g |>
    filter(signature %in% SIG_LEVELS) |>
    mutate(dir = ifelse(NES > 0, "Up", "Down")) |>
    group_by(signature, dir) |>
    summarise(n05 = sum(padj < 0.05), n10 = sum(padj < 0.10), .groups = "drop") |>
    complete(signature = SIG_LEVELS, dir = c("Down", "Up"),
             fill = list(n05 = 0, n10 = 0)) |>
    mutate(group = factor(signature, SIG_LEVELS, SIG_SHORT),
           dir = factor(dir, c("Down", "Up")))
  cat("Panel b two-tier burden:\n"); print(d |> select(group, dir, n05, n10))
  dodge <- position_dodge(width = 0.78)
  ggplot(d, aes(group, group = dir)) +
    geom_col(aes(y = n10, fill = paste0(dir, ", FDR < 0.10")),
             position = dodge, width = 0.7) +
    geom_col(aes(y = n05, fill = paste0(dir, ", FDR < 0.05")),
             position = dodge, width = 0.7) +
    geom_text(aes(y = n10, label = n10), position = dodge,
              vjust = -0.35, size = 4.4, colour = "grey25") +
    scale_fill_manual(
      values = c("Up, FDR < 0.05" = UP_DARK, "Up, FDR < 0.10" = UP_LIGHT,
                 "Down, FDR < 0.05" = DOWN_DARK, "Down, FDR < 0.10" = DOWN_LIGHT),
      breaks = c("Down, FDR < 0.05", "Down, FDR < 0.10",
                 "Up, FDR < 0.05", "Up, FDR < 0.10"), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = NULL, y = "Significant gene sets (count)") +
    guides(fill = guide_legend(ncol = 1)) +
    theme_cowplot(font_size = FS) +
    theme(legend.position = c(1, 1), legend.justification = c(1, 1),
          legend.text = element_text(size = FS - 3),
          legend.key.size = unit(0.9, "lines"),
          axis.text.x = element_text(angle = 30, hjust = 1))
})

# ---- module gene lists (for c colors and e exemplars) -------------------------
blocks <- tribble(
  ~pathway, ~block, ~label,
  "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING", "Shared: synaptic program", "Regulation of trans-synaptic signaling (GO:BP)",
  "GOCC_SYNAPTIC_MEMBRANE",                      "Shared: synaptic program", "Synaptic membrane (GO:CC)",
  "GOBP_NEUROTRANSMITTER_SECRETION",             "Shared: synaptic program", "Neurotransmitter secretion (GO:BP)",
  "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES","Shared: synaptic program","Transmission across chemical synapses (Reactome)",
  "GOBP_PROTEIN_K63_LINKED_DEUBIQUITINATION",    "Shared: ubiquitin", "Protein K63-linked deubiquitination (GO:BP)",
  "GOMF_UBIQUITIN_CONJUGATING_ENZYME_BINDING",   "Shared: ubiquitin", "Ubiquitin-conjugating enzyme binding (GO:MF)",
  "REACTOME_TRANSLATION",                        "Graded: cytosolic translation", "Translation (Reactome)",
  "GOCC_RIBOSOMAL_SUBUNIT",                      "Graded: cytosolic translation", "Ribosomal subunit (GO:CC)",
  "GOBP_CYTOPLASMIC_TRANSLATION",                "Graded: cytosolic translation", "Cytoplasmic translation (GO:BP)",
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",  "Graded: cytosolic translation", "Eukaryotic translation elongation (Reactome)",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",          "Graded: oxidative phosphorylation", "Oxidative phosphorylation (Hallmark)",
  "GOBP_OXIDATIVE_PHOSPHORYLATION",              "Graded: oxidative phosphorylation", "Oxidative phosphorylation (GO:BP)",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",     "Graded: oxidative phosphorylation", "Respiratory electron transport (Reactome)")
BLOCK_LEVELS <- c("Shared: synaptic program", "Shared: ubiquitin",
                  "Graded: cytosolic translation", "Graded: oxidative phosphorylation")

le_union <- function(pws) {
  g |> filter(pathway %in% pws, signature == "depleted") |>
    pull(leadingEdge) |> strsplit("|", fixed = TRUE) |> unlist() |> unique()
}
mods <- lapply(setNames(BLOCK_LEVELS, BLOCK_LEVELS),
               function(b) le_union(blocks$pathway[blocks$block == b]))
# name the modules explicitly: positional indexing silently mislabels panels
# whenever a block is added or reordered
M_SYN <- mods[["Shared: synaptic program"]]
M_DUB <- mods[["Shared: ubiquitin"]]
M_TRA <- mods[["Graded: cytosolic translation"]]
M_OXP <- mods[["Graded: oxidative phosphorylation"]]

# ---- shared exemplar selection (panels c and e) -------------------------------
sw_all <- sig |> select(gene, stratum, z) |>
  pivot_wider(names_from = stratum, values_from = z) |> drop_na()
pick_dn <- function(genes, n = 8, shared = FALSE) {
  d <- sw_all |> filter(gene %in% genes)
  if (shared) d <- d |> filter(depleted < 0, intermediate < 0, non_depleted < 0) |>
      arrange(pmax(depleted, intermediate, non_depleted))
  else d <- d |> arrange(depleted)
  head(d$gene, n)
}
# graded blocks: core members of the GO group that are also most DISCORDANT
# between strata (down in depleted, at/above zero in non-depleted)
pick_graded <- function(genes, core_rx, n = 8, max_dep = -1.2) {
  sw_all |> filter(gene %in% genes, str_detect(gene, core_rx), depleted < max_dep) |>
    arrange(desc(non_depleted - depleted)) |> head(n) |> pull(gene)
}
pick_up <- function(genes, n = 6) {
  sw_all |> filter(gene %in% genes, depleted > 0, non_depleted > 0) |>
    arrange(desc(pmin(depleted, non_depleted))) |> head(n) |> pull(gene)
}
SEL <- list(unique(c("SST", pick_dn(M_SYN, 8, shared = TRUE))),
            pick_up(M_DUB, 6),
            pick_graded(M_TRA, "^RP[LS][0-9]|^EIF3", 8),
            pick_graded(M_OXP, "^NDUF|^UQCR|^COX[0-9]|^CYC[S1]|^ATP5|^SDH", 8)) |>
  setNames(BLOCK_LEVELS)
# genes labelled in c are the top two of each block in e (bolded there)
HL <- unlist(lapply(SEL, head, 2), use.names = FALSE)
cat("panel c labels / bolded in e:", paste(HL, collapse = ", "), "\n")

# ---- c | gene-level depleted vs non-depleted scatter --------------------------
# Synaptic (the largest module) is drawn small and translucent so the two graded
# modules read clearly.
p_gsc <- local({
  sw <- sw_all |>
    mutate(module = case_when(gene %in% M_DUB ~ "Deubiquitination",
                              gene %in% M_TRA ~ "Translation",
                              gene %in% M_OXP ~ "OxPhos/mito",
                              gene %in% M_SYN ~ "Synaptic",
                              TRUE ~ "Other"),
           module = factor(module, names(FAM_COLS)))
  r <- cor(sw$depleted, sw$non_depleted)
  hl <- sw |> filter(gene %in% HL)
  # manual label positions: the four graded-module genes sit inside the dense
  # cloud, so their labels are parked in the (near-empty) upper-left corner with
  # leader lines. hj = 1 -> text runs leftward from the anchor.
  labpos <- tribble(
    ~gene,       ~lx,    ~ly, ~hj,
    "EIF3G",   -3.45,  4.05,   1,
    "COX5B",   -3.45,  3.35,   1,
    "RPL36",   -3.45,  2.65,   1,
    "NDUFS8",  -3.45,  1.95,   1,
    "SST",     -3.70, -1.30,   1,
    "VGF",     -3.45, -3.75,   1,
    "STAMBPL1", 2.55,  4.25,   0,
    "USP8",     3.85,  1.10,   0) |>
    left_join(hl |> select(gene, depleted, non_depleted), by = "gene") |>
    mutate(sx = lx + ifelse(hj == 1, 0.10, -0.10))
  ggplot(sw, aes(depleted, non_depleted)) +
    geom_hline(yintercept = 0, colour = "grey88") +
    geom_vline(xintercept = 0, colour = "grey88") +
    geom_abline(linetype = "dashed", colour = "grey65") +
    geom_point(data = ~ filter(.x, module == "Other"), aes(colour = module),
               size = 0.4, alpha = 0.22) +
    geom_point(data = ~ filter(.x, module == "Synaptic"), aes(colour = module),
               size = 0.9, alpha = 0.40) +
    geom_point(data = ~ filter(.x, module %in% c("Deubiquitination", "Translation",
                                                 "OxPhos/mito")),
               aes(colour = module), size = 1.7, alpha = 0.92) +
    geom_point(data = hl, aes(fill = module), colour = "black", shape = 21,
               size = 2.4, stroke = 0.8, show.legend = FALSE) +
    geom_segment(data = labpos, aes(x = sx, y = ly, xend = depleted,
                                    yend = non_depleted),
                 inherit.aes = FALSE, colour = "grey50", linewidth = 0.28) +
    geom_text(data = labpos, aes(x = lx, y = ly, label = gene, hjust = hj),
              inherit.aes = FALSE, colour = "black", size = 4,
              fontface = "italic") +
    scale_colour_manual(values = FAM_COLS, name = NULL,
                        breaks = c("Synaptic", "Translation", "OxPhos/mito",
                                   "Deubiquitination", "Other")) +
    scale_fill_manual(values = FAM_COLS, guide = "none") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1.5, size = 4.4,
             label = sprintf("r = %.2f", r)) +
    coord_cartesian(xlim = c(-5.4, 5.4), ylim = c(-5.4, 4.9)) +
    labs(x = "Gene z, depleted stratum", y = "Gene z, non-depleted stratum") +
    theme_cowplot(font_size = FS) +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 1,
                                 override.aes = list(size = 3, alpha = 1)))
})

# ---- d | pathway NES heatmap (blocks x groups) --------------------------------
D_SHOW <- c("GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING", "GOCC_SYNAPTIC_MEMBRANE",
            "GOBP_PROTEIN_K63_LINKED_DEUBIQUITINATION",
            "GOMF_UBIQUITIN_CONJUGATING_ENZYME_BINDING",
            "GOBP_CYTOPLASMIC_TRANSLATION", "GOCC_RIBOSOMAL_SUBUNIT",
            "GOBP_OXIDATIVE_PHOSPHORYLATION", "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT")
p_path <- local({
  d <- g |>
    filter(signature %in% SIG_LEVELS, pathway %in% D_SHOW) |>
    inner_join(blocks, by = "pathway") |>
    mutate(label = str_wrap(label, 26),
           group = factor(signature, SIG_LEVELS, SIG_LABELS),
           block = factor(block, BLOCK_LEVELS),
           stars = ifelse(padj < 0.01, "**", ifelse(padj < 0.05, "*",
                   ifelse(padj < 0.1, "+", ""))))
  ord <- d |> filter(signature == "depleted") |> arrange(block, NES) |> pull(label)
  d$label <- factor(d$label, levels = rev(ord))
  levels(d$group) <- SIG_SHORT
  ggplot(d, aes(group, label, fill = NES)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = stars), size = 6, vjust = 0.75) +
    nes_fill(c(-3, 3)) +
    facet_grid(block ~ ., scales = "free_y", space = "free_y",
               labeller = labeller(block = label_wrap_gen(12))) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_cowplot(font_size = FS) +
    theme(axis.text.x.top = element_text(size = FS, angle = 45, hjust = 0),
          strip.background = element_rect(fill = "grey92"),
          strip.text.y = element_text(angle = 270, size = FS - 3, face = "bold"),
          axis.ticks = element_blank(), axis.line = element_blank(),
          legend.position = "bottom")
})

# ---- e | leading-edge exemplar genes (c-labelled genes in bold) ---------------
p_gene <- local({
  sub_z <- read_csv(SUBCLASS, show_col_types = FALSE) |>
    filter(cell_type == "Sst", !is.na(estimate), !is.na(se), se > 0) |>
    mutate(z = estimate / se) |>
    arrange(desc(abs(z))) |> distinct(genes, .keep_all = TRUE) |>
    transmute(gene = genes, stratum = "Sst_subclass", z)
  d <- imap(SEL, function(gg, b) {
    bind_rows(sig |> filter(gene %in% gg) |> select(gene, stratum, z),
              sub_z |> filter(gene %in% gg)) |> mutate(block = b)
  }) |>
    bind_rows() |>
    mutate(stratum = factor(stratum, SIG_LEVELS, SIG_SHORT),
           block = factor(block, BLOCK_LEVELS))
  ord <- d |> filter(stratum == "Depleted") |> arrange(block, z) |> pull(gene)
  d$gene <- factor(d$gene, levels = rev(unique(ord)))
  lab <- d |> distinct(gene, block) |> mutate(hl = as.character(gene) %in% HL)
  ggplot(d, aes(stratum, gene, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f", z)), size = 3.9) +
    geom_text(data = lab, aes(x = 0.44, y = gene, label = gene),
              inherit.aes = FALSE, hjust = 1, size = 3.9, colour = "grey12",
              fontface = ifelse(lab$hl, "bold.italic", "italic")) +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                         midpoint = 0, limits = c(-4, 4), oob = scales::squish,
                         name = "Meta-DE z (SCZ vs control)") +
    facet_grid(block ~ ., scales = "free_y", space = "free_y",
               labeller = labeller(block = label_wrap_gen(12))) +
    scale_x_discrete(position = "top", expand = expansion(add = c(1.32, 0.55))) +
    labs(x = NULL, y = NULL) +
    theme_cowplot(font_size = FS) +
    theme(axis.text.x.top = element_text(size = FS, angle = 45, hjust = 0),
          axis.text.y = element_blank(), axis.ticks = element_blank(),
          strip.background = element_rect(fill = "grey92"),
          strip.text.y = element_text(angle = 270, size = FS - 3, face = "bold"),
          axis.line = element_blank(), legend.position = "bottom")
})

# ---- assembly ------------------------------------------------------------------
top <- plot_grid(p_def, p_burden, p_gsc, nrow = 1,
                 rel_widths = c(0.37, 0.24, 0.39),
                 labels = c("a", "b", "c"), label_size = 24)
bot <- plot_grid(p_path, p_gene, nrow = 1, rel_widths = c(0.52, 0.48),
                 labels = c("d", "e"), label_size = 24)
fig <- plot_grid(top, bot, ncol = 1, rel_heights = c(0.42, 0.58))
ggsave(file.path(OUT, "figS_sst_strata_mockup.png"), fig,
       width = 17.5, height = 13, dpi = 200, bg = "white")
ggsave(file.path(OUT, "figS_sst_strata_mockup.pdf"), fig,
       width = 17.5, height = 13, bg = "white")
cat("\nWrote figS_sst_strata_mockup.(png|pdf)\n")
