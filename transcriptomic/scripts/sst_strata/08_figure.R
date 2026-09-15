#!/usr/bin/env Rscript
# Step 8 | Supplementary figure S8. Laid out at main-figure size, so it can be
# promoted to a main figure by changing FIGSTEM and the output directory at the
# bottom.
#
#   a  stratum definition: crumblr abundance change vs median cortical depth for
#      the 16 Sst supertypes, with the layer gutter and the FDR < 0.20 boundary
#   b  gene-set burden per group, two FDR tiers, pooled subclass as reference
#   c  gene-level z, depleted vs non-depleted, coloured by module
#   d  NES heatmap: shared blocks (synaptic, ubiquitin) vs graded blocks
#      (translation, oxidative phosphorylation)
#   e  exemplar leading-edge genes per block; the genes labelled in c are bold
#
# Purely a rendering step -- no statistics are computed here. Panels a and b use
# the strata palette and the Figure 2 direction palette respectively; c-e use the
# module palette. All three come from _common.R so the supplements match.
source("transcriptomic/scripts/sst_strata/_common.R")
suppressPackageStartupMessages({
  library(ggplot2); library(cowplot); library(ggrepel)
})

g      <- read_csv(file.path(P$pb, "gsea_all_signatures.csv"), show_col_types = FALSE)
sig    <- read_csv(file.path(P$pb, "stratum_gene_signatures.csv"), show_col_types = FALSE)
strata <- load_strata()

S_LAB    <- strata_labels()                       # "Depleted (5)" etc, counted from the table
STRAT_BY_LABEL <- setNames(unname(STRAT_COLS[STRATA_LEVELS]), S_LAB[STRATA_LEVELS])
mods <- modules(g)[MODULE_PRIORITY]               # priority order matters for c
say("module sizes: %s", paste(sprintf("%s %d", names(mods), lengths(mods)), collapse = ", "))

# ---- a | stratum definition ---------------------------------------------------
p_def <- local({
  d <- strata |> mutate(stratum = factor(stratum, STRATA_LEVELS, S_LAB[STRATA_LEVELS]))
  # midpoint between the least-depleted "depleted" supertype and the most
  # negative "intermediate" one -- i.e. where crumblr FDR crosses 0.20
  cut_fdr <- mean(c(max(d$estimate[d$stratum == S_LAB[["depleted"]]]),
                    min(d$estimate[d$stratum == S_LAB[["intermediate"]]])))
  lb <- read_csv(P$layers, show_col_types = FALSE)$xenium
  layers <- tibble(name = c("L1", "L2/3", "L4", "L5", "L6"),
                   mid = (c(0, lb[-length(lb)]) + lb) / 2)
  XR <- 0.155                                     # right edge of the layer gutter
  ggplot(d, aes(estimate, depth_xenium, colour = stratum)) +
    annotate("segment", x = XR - 0.012, xend = XR, y = c(0, lb), yend = c(0, lb),
             colour = "grey55", linewidth = 0.5 * GEO) +
    annotate("text", x = XR, y = layers$mid, label = layers$name,
             hjust = 1, size = LBL_SMALL, colour = "grey30") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = cut_fdr, linetype = "dashed", colour = "grey60") +
    # "crumblr" dropped from the annotation: the x axis title already names the
    # estimate, and the full string overruns the panel at 7.1 in
    annotate("text", x = cut_fdr - 0.006, y = 0.985, hjust = 1, size = LBL_SMALL,
             colour = "grey30", label = "FDR < 0.20") +
    # strata named in-panel: they partition the x axis, so a legend is redundant
    # counts dropped and the middle label staggered: at 7.1 in the three full
    # labels are wider than the x range and collide
    annotate("text", x = c(-0.255, -0.088, 0.052), y = c(0.035, 0.105, 0.035),
             label = c("Depleted", "Intermediate", "Non-depleted"),
             colour = unname(STRAT_COLS[STRATA_LEVELS]),
             fontface = "bold", size = LBL_CALL) +
    geom_errorbarh(aes(xmin = estimate - se, xmax = estimate + se),
                   height = 0, linewidth = 0.5 * GEO, alpha = 0.6) +
    geom_point(size = 3 * GEO) +
    geom_text_repel(aes(label = CellType), size = LBL_SMALL, seed = 1,
                    show.legend = FALSE, max.overlaps = 40,
                    box.padding = 0.32, point.padding = 0.12, force = 3,
                    min.segment.length = 0.2, segment.size = 0.18,
                    segment.colour = "grey55") +
    scale_y_reverse() +
    scale_colour_manual(values = STRAT_BY_LABEL, name = NULL) +
    labs(x = "Abundance change in SCZ (crumblr estimate ± SE)",
         y = "Median cortical depth (0 = pia)") +
    theme_strata() +
    theme(legend.position = "none")
})

# ---- b | gene-set burden, two FDR tiers --------------------------------------
p_burden <- local({
  d <- g |> filter(signature %in% SIG_LEVELS) |>
    mutate(dir = ifelse(NES > 0, "Up", "Down")) |>
    group_by(signature, dir) |>
    summarise(n05 = sum(padj < 0.05), n10 = sum(padj < 0.10), .groups = "drop") |>
    complete(signature = SIG_LEVELS, dir = c("Down", "Up"),
             fill = list(n05 = 0, n10 = 0)) |>
    mutate(group = factor(signature, SIG_LEVELS, SIG_SHORT),
           dir = factor(dir, c("Down", "Up")))
  say("panel b burden (labels are the FDR<0.10 tier):")
  print(d |> select(group, dir, n05, n10))
  dodge <- position_dodge(width = 0.78)
  ggplot(d, aes(group, group = dir)) +
    geom_col(aes(y = n10, fill = paste0(dir, ", FDR < 0.10")), position = dodge, width = 0.7) +
    geom_col(aes(y = n05, fill = paste0(dir, ", FDR < 0.05")), position = dodge, width = 0.7) +
    geom_text(aes(y = n10, label = n10), position = dodge, vjust = -0.35,
              size = LBL_CALL, colour = "grey25") +
    scale_fill_manual(
      values = c("Up, FDR < 0.05" = COL_UP, "Up, FDR < 0.10" = UP_LIGHT,
                 "Down, FDR < 0.05" = COL_DOWN, "Down, FDR < 0.10" = DOWN_LIGHT),
      breaks = c("Down, FDR < 0.05", "Down, FDR < 0.10",
                 "Up, FDR < 0.05", "Up, FDR < 0.10"), name = NULL) +
    # extra headroom so the inset legend clears the tallest bar and its label
    scale_y_continuous(expand = expansion(mult = c(0, 0.24))) +
    labs(x = NULL, y = "Significant gene sets (count)") +
    guides(fill = guide_legend(ncol = 1)) +
    theme_strata() +
    # inset top-right: the Non-depleted column is empty and Intermediate reaches
    # only 59, so this corner is free once the y axis has headroom
    theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
          legend.text = element_text(size = LEGEND_TEXT),
          legend.key.size = unit(7 * FIG_SCALE, "pt"),
          legend.spacing.y = unit(0, "pt"),
          legend.background = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
})

# ---- exemplar gene selection (shared by c and e) ------------------------------
sw_all <- sig |> select(gene, stratum, z) |>
  pivot_wider(names_from = stratum, values_from = z) |> drop_na()

# shared blocks: genes moving the same way in all three strata
pick_shared <- function(genes, n = 8) {
  sw_all |> filter(gene %in% genes, depleted < 0, intermediate < 0, non_depleted < 0) |>
    arrange(pmax(depleted, intermediate, non_depleted)) |> head(n) |> pull(gene)
}
pick_up <- function(genes, n = 6) {
  sw_all |> filter(gene %in% genes, depleted > 0, non_depleted > 0) |>
    arrange(desc(pmin(depleted, non_depleted))) |> head(n) |> pull(gene)
}
# graded blocks: core members of the GO group (regex) that are also the most
# DISCORDANT between strata -- down in depleted, at or above zero in
# non-depleted. Without the regex the picks drift to peripheral members
# (mitoribosomal genes for "cytosolic translation", TCA enzymes for OxPhos).
pick_graded <- function(genes, core_rx, n = 8, max_dep = -1.2) {
  sw_all |> filter(gene %in% genes, str_detect(gene, core_rx), depleted < max_dep) |>
    arrange(desc(non_depleted - depleted)) |> head(n) |> pull(gene)
}
SEL <- list(
  synaptic    = unique(c("SST", pick_shared(mods$synaptic, 8))),
  ubiquitin   = pick_up(mods$ubiquitin, 6),
  translation = pick_graded(mods$translation, "^RP[LS][0-9]|^EIF3", 8),
  oxphos      = pick_graded(mods$oxphos, "^NDUF|^UQCR|^COX[0-9]|^CYC[S1]|^ATP5|^SDH", 8)
)[BLOCK_KEYS]
HL <- unlist(lapply(SEL, head, 2), use.names = FALSE)   # labelled in c, bold in e
say("panel c labels / bold in e: %s", paste(HL, collapse = ", "))

# ---- c | gene z, depleted vs non-depleted ------------------------------------
p_gsc <- local({
  sw <- sw_all |>
    mutate(module = factor(FAM_LABELS[module_of(gene, mods)], unname(FAM_LABELS)))
  say("panel c module counts: %s",
      paste(capture.output(print(count(sw, module))), collapse = " | "))
  # Spearman, matching the correlation reported in Figs 2, 3 and 4
  rho <- cor(sw$depleted, sw$non_depleted, method = "spearman")
  hl <- sw |> filter(gene %in% HL)
  # The four graded-module exemplars sit inside the dense cloud, so their labels
  # are placed by hand in the near-empty upper-left corner with leader lines,
  # ordered so the leaders fan out without crossing. Coordinates are data units;
  # revisit if the exemplar genes change.
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
    inner_join(hl |> select(gene, depleted, non_depleted), by = "gene") |>
    mutate(sx = lx + ifelse(hj == 1, 0.10, -0.10))
  if (nrow(labpos) < length(HL))
    say("NOTE: %d exemplar(s) have no manual label position: %s",
        length(HL) - nrow(labpos), paste(setdiff(HL, labpos$gene), collapse = ", "))
  big <- unname(FAM_LABELS[c("ubiquitin", "translation", "oxphos")])
  ggplot(sw, aes(depleted, non_depleted)) +
    geom_hline(yintercept = 0, colour = "grey88") +
    geom_vline(xintercept = 0, colour = "grey88") +
    geom_abline(linetype = "dashed", colour = "grey65") +
    geom_point(data = ~ filter(.x, module == FAM_LABELS[["other"]]),
               aes(colour = module), size = 0.4 * GEO, alpha = 0.22) +
    geom_point(data = ~ filter(.x, module == FAM_LABELS[["synaptic"]]),
               aes(colour = module), size = 0.9 * GEO, alpha = 0.40) +
    geom_point(data = ~ filter(.x, module %in% big), aes(colour = module),
               size = 1.7 * GEO, alpha = 0.92) +
    geom_point(data = hl, aes(fill = module), colour = "black", shape = 21,
               size = 2.4 * GEO, stroke = 0.8 * GEO, show.legend = FALSE) +
    geom_segment(data = labpos, aes(x = sx, y = ly, xend = depleted, yend = non_depleted),
                 inherit.aes = FALSE, colour = "grey50", linewidth = 0.28 * GEO) +
    geom_text(data = labpos, aes(x = lx, y = ly, label = gene, hjust = hj),
              inherit.aes = FALSE, colour = "black", size = LBL_GENE, fontface = "italic") +
    scale_colour_manual(values = setNames(unname(FAM_COLS), unname(FAM_LABELS)),
                        name = NULL, breaks = unname(FAM_LABELS)) +
    scale_fill_manual(values = setNames(unname(FAM_COLS), unname(FAM_LABELS)), guide = "none") +
    # plotmath rather than a literal rho: the base pdf device cannot encode U+03C1
    annotate("text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1.5, size = LBL_CALL,
             label = sprintf("rho == %.2f", rho), parse = TRUE) +
    coord_cartesian(xlim = c(-5.4, 5.4), ylim = c(-5.4, 4.9)) +
    labs(x = "Gene z, depleted stratum", y = "Gene z, non-depleted stratum") +
    theme_strata() +
    # inset bottom-right: verified empty, 0 of 9,631 genes fall in x > 2.2, y < -2.2
    theme(legend.position = c(1.0, 0.0), legend.justification = c(1, 0),
          legend.text = element_text(size = LEGEND_TEXT),
          legend.key.size = unit(7 * FIG_SCALE, "pt"),
          legend.spacing.y = unit(0, "pt"),
          legend.background = element_blank()) +
    guides(colour = guide_legend(ncol = 1,
                                 override.aes = list(size = 1.6, alpha = 1)))
})

# ---- d | NES heatmap ---------------------------------------------------------
# Two gene sets per block: the clearest GO/Reactome representative of each.
D_SHOW <- c("GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING", "GOCC_SYNAPTIC_MEMBRANE",
            "GOBP_PROTEIN_K63_LINKED_DEUBIQUITINATION",
            "GOMF_UBIQUITIN_CONJUGATING_ENZYME_BINDING",
            "GOBP_CYTOPLASMIC_TRANSLATION", "GOCC_RIBOSOMAL_SUBUNIT",
            "GOBP_OXIDATIVE_PHOSPHORYLATION", "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT")
p_path <- local({
  d <- g |> filter(signature %in% SIG_LEVELS, pathway %in% D_SHOW) |>
    inner_join(BLOCKS, by = "pathway") |>
    mutate(label = str_wrap(label, 30),
           group = factor(signature, SIG_LEVELS, SIG_SHORT),
           block = factor(BLOCK_LABELS[block], unname(BLOCK_LABELS)),
           stars = stars_of(padj))
  ord <- d |> filter(signature == "depleted") |> arrange(block, NES) |> pull(label)
  d$label <- factor(d$label, levels = rev(ord))
  ggplot(d, aes(group, label, fill = NES)) +
    geom_tile(colour = "white", linewidth = LW) +
    geom_text(aes(label = stars), size = LBL_STAR, vjust = 0.75) +
    nes_fill(c(-3, 3)) +
    facet_grid(block ~ ., scales = "free_y", space = "free_y",
               labeller = labeller(block = label_wrap_gen(12))) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_strata() +
    theme(axis.text.x.top = element_text(size = AXIS_TEXT, angle = 45, hjust = 0),
          strip.background = element_rect(fill = "grey92"),
          strip.text.y = element_text(angle = 270, size = LEGEND_TEXT, face = "bold"),
          axis.ticks = element_blank(), axis.line = element_blank(),
          legend.position = "bottom")
})

# ---- e | exemplar gene z heatmap --------------------------------------------
p_gene <- local({
  # pooled-subclass column comes from the Fig 2 gene-level table
  sub_z <- read_csv(P$subclass, show_col_types = FALSE) |>
    filter(cell_type == "Sst", !is.na(estimate), !is.na(se), se > 0) |>
    mutate(z = estimate / se) |>
    arrange(desc(abs(z))) |> distinct(genes, .keep_all = TRUE) |>
    transmute(gene = genes, stratum = "Sst_subclass", z)
  d <- imap(SEL, function(gg, b) {
    bind_rows(sig |> filter(gene %in% gg) |> select(gene, stratum, z),
              sub_z |> filter(gene %in% gg)) |> mutate(block = b)
  }) |> bind_rows() |>
    mutate(stratum = factor(stratum, SIG_LEVELS, SIG_SHORT),
           block = factor(BLOCK_LABELS[block], unname(BLOCK_LABELS)))
  ord <- d |> filter(stratum == "Depleted") |> arrange(block, z) |> pull(gene)
  d$gene <- factor(d$gene, levels = rev(unique(ord)))
  # gene names drawn as data (not axis text) so the exemplars shared with panel c
  # can be bolded individually
  lab <- d |> distinct(gene, block) |> mutate(hl = as.character(gene) %in% HL)
  ggplot(d, aes(stratum, gene, fill = z)) +
    geom_tile(colour = "white", linewidth = LW) +
    geom_text(aes(label = sprintf("%.1f", z)), size = LBL_SMALL) +
    geom_text(data = lab, aes(x = 0.44, y = gene, label = gene), inherit.aes = FALSE,
              hjust = 1, size = LBL_GENE, colour = "grey12",
              fontface = ifelse(lab$hl, "bold.italic", "italic")) +
    nes_fill(c(-4, 4), name = "Meta-DE z (SCZ vs control)") +
    facet_grid(block ~ ., scales = "free_y", space = "free_y",
               labeller = labeller(block = label_wrap_gen(12))) +
    scale_x_discrete(position = "top", expand = expansion(add = c(1.75, 0.55))) +
    labs(x = NULL, y = NULL) +
    theme_strata() +
    theme(axis.text.x.top = element_text(size = AXIS_TEXT, angle = 45, hjust = 0),
          axis.text.y = element_blank(), axis.ticks = element_blank(),
          strip.background = element_rect(fill = "grey92"),
          strip.text.y = element_text(angle = 270, size = LEGEND_TEXT, face = "bold"),
          axis.line = element_blank(), legend.position = "bottom")
})

# ---- assembly ----------------------------------------------------------------
# align = "h" + axis = "tb" puts the three plot panels on a common top and
# bottom edge, so the x axes line up despite a having no legend and b and c
# carrying two-row legends underneath
top <- plot_grid(p_def, p_burden, p_gsc, nrow = 1, rel_widths = c(0.35, 0.27, 0.38),
                 align = "h", axis = "tb",
                 labels = c("a", "b", "c"), label_size = PANEL_LABEL)
bot <- plot_grid(p_path, p_gene, nrow = 1, rel_widths = c(0.52, 0.48),
                 align = "h", axis = "tb",
                 labels = c("d", "e"), label_size = PANEL_LABEL)
FIG_H <- 9.8   # at 8.0 in wide this prints as 7.1 x 8.7 in; panel e carries 31 gene
               # rows, and aligning a-c costs panel a the height b and c give to legends
fig <- plot_grid(top, bot, ncol = 1, rel_heights = c(0.41, 0.59))
# S-number lives here, per manuscript/figures/supplementary/README.md; the figure
# is written directly into the submission folder like the other supplements
FIGSTEM <- "S08_sst_strata"
for (ext in c("png", "pdf"))
  ggsave(file.path(P$suppfig, paste0(FIGSTEM, ".", ext)), fig,
         width = FIG_W, height = FIG_H, dpi = FIG_DPI, bg = "white")
say("wrote %s/%s.(png|pdf)", P$suppfig, FIGSTEM)
