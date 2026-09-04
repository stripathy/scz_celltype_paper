#!/usr/bin/env Rscript
# Exploratory / decision-support: stratified DE + GSEA across Sst supertype groups.
#
# Question (Etienne's review, 2026-08): are dysregulated pathways shared or
# distinct between depleted, intermediate, and non-depleted Sst supertypes?
# Shared signal cannot be a survivor-selection artifact (the non-depleted group
# is not losing cells); depleted-specific signal is the candidate vulnerability
# signature (with the selection caveat).
#
# Inputs: per-supertype 7-cohort meta-analytic DE (sst_meta_de/SST_*_meta.csv),
# SCZ supertype crumblr (final_results_crumblr_7_cohorts.csv), supertype depth
# medians (spatial depth_platform), and Nicole's pre-existing aggregated-group
# DE runs (Unaffected / Vulnerable_AD / Vulnerable_SCZAD) as a sensitivity.
#
# Stratification DECISIONS (documented, not re-litigated downstream):
#   depleted     = crumblr estimate < 0 & padj < 0.20
#                  -> Sst_2, Sst_3, Sst_20, Sst_22, Sst_25
#                  (FDR<0.20 per Shreejoy 2026-08-29, so borderline Sst_20 —
#                   nominal p = 0.019, padj = 0.19 — counts as depleted)
#   intermediate = estimate < 0 & padj >= 0.20
#   non_depleted = estimate >= 0
# Stratum gene signature = inverse-variance-weighted meta of per-supertype
# estimates (gene kept if present in >= 2 supertypes of the stratum). Supertypes
# share donors, so IVW SEs are anticonservative; the z is used as a RANKING for
# preranked GSEA (house convention), not for gene-level inference.
#
# Run from repo root:  Rscript transcriptomic/scripts/17_sst_strata_gsea.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(tibble); library(ggplot2); library(cowplot); library(ggrepel)
  library(fgsea); library(msigdbr)
})

set.seed(42)

DE_DIR    <- "shared/snrnaseq_de/nicole_scz_snrnaseq_betas/sst_meta_de"
CRUMBLR   <- "shared/snrnaseq_de/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv"
DEPTH     <- "spatial/output/depth_platform/supertype_depth_platform_summary.csv"
SUBCLASS  <- "shared/snrnaseq_de/DE_genes_all_cells_scz.csv"
OUT_DIR   <- "transcriptomic/results/sst_strata_gsea"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_SIZE <- 10; MAX_SIZE <- 500          # house GSEA convention (script 02)
FS_LAB <- 16; FS_TITLE <- 18             # figure text standards

STRATUM_COLS <- c(depleted = "#D55E00", intermediate = "#E69F00",
                  non_depleted = "#0072B2", Sst_subclass = "grey40")

t0 <- Sys.time()
say <- function(...) cat(sprintf("[%5.1f min] ", as.numeric(difftime(Sys.time(), t0, units = "mins"))), sprintf(...), "\n", sep = "")

# ---- 1. Strata from SCZ crumblr ---------------------------------------------
crum <- read_csv(CRUMBLR, show_col_types = FALSE) |>
  filter(str_starts(CellType, "Sst"), !str_detect(CellType, "Chodl")) |>
  mutate(stratum = case_when(
    estimate < 0 & padj < 0.20 ~ "depleted",
    estimate < 0               ~ "intermediate",
    TRUE                       ~ "non_depleted"))

depth <- read_csv(DEPTH, show_col_types = FALSE) |>
  filter(subclass == "Sst") |>
  select(supertype, depth_xenium = median_Xenium, depth_merfish = median_MERFISH)

strata <- crum |>
  left_join(depth, by = c(CellType = "supertype")) |>
  arrange(estimate)
write_csv(strata, file.path(OUT_DIR, "strata_definition.csv"))
say("Strata (n = %d supertypes):", nrow(strata))
print(strata |> select(CellType, estimate, padj, stratum, depth_xenium), n = 20)
say("Spearman(crumblr estimate, Xenium depth) = %.2f",
    cor(strata$estimate, strata$depth_xenium, method = "spearman", use = "complete.obs"))
say("Mean Xenium depth by stratum:")
print(strata |> group_by(stratum) |> summarise(n = n(), mean_depth = mean(depth_xenium, na.rm = TRUE)))

# ---- 2. Per-supertype DE -> long table --------------------------------------
files <- list.files(DE_DIR, pattern = "^SST_\\d+_meta\\.csv$", full.names = TRUE)
de <- map_dfr(files, function(f) {
  st <- str_replace(basename(f), "^SST_(\\d+)_meta\\.csv$", "Sst_\\1")
  read_csv(f, show_col_types = FALSE) |>
    transmute(supertype = st, gene = Gene, estimate, se, zval)
}) |>
  filter(!is.na(estimate), !is.na(se), se > 0)
de <- de |> inner_join(strata |> select(supertype = CellType, stratum), by = "supertype")
say("Loaded %d gene x supertype rows across %d supertypes",
    nrow(de), n_distinct(de$supertype))

# sanity: zval should equal estimate/se
stopifnot(max(abs(de$zval - de$estimate / de$se), na.rm = TRUE) < 1e-6)

# ---- 3. Supertype x supertype DE concordance heatmap ------------------------
zw <- de |> select(gene, supertype, zval) |>
  pivot_wider(names_from = supertype, values_from = zval)
sts <- strata |> arrange(factor(stratum, levels = names(STRATUM_COLS)), depth_xenium) |> pull(CellType)
sts <- intersect(sts, colnames(zw))
cm <- cor(as.matrix(zw[, sts]), use = "pairwise.complete.obs", method = "spearman")
cm_long <- as_tibble(cm, rownames = "a") |>
  pivot_longer(-a, names_to = "b", values_to = "rho") |>
  mutate(a = factor(a, levels = sts), b = factor(b, levels = rev(sts)))
lab_col <- STRATUM_COLS[strata$stratum[match(sts, strata$CellType)]]
p_cor <- ggplot(cm_long, aes(a, b, fill = rho)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3.2) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho (z)") +
  labs(x = NULL, y = NULL,
       title = "Cross-supertype concordance of SCZ DE (gene-level z)",
       subtitle = "Ordered by stratum (depleted, intermediate, non-depleted), then depth") +
  theme_cowplot(font_size = FS_LAB) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = lab_col),
        axis.text.y = element_text(colour = rev(lab_col)),
        plot.title = element_text(size = FS_TITLE))
ggsave(file.path(OUT_DIR, "fig1_supertype_z_correlation.png"), p_cor,
       width = 11.5, height = 10, dpi = 200, bg = "white")
say("Median within-stratum rho vs between-stratum rho:")
same <- outer(strata$stratum[match(sts, strata$CellType)],
              strata$stratum[match(sts, strata$CellType)], "==")
diag(same) <- NA
say("  within = %.3f | between = %.3f",
    median(cm[which(same)], na.rm = TRUE), median(cm[which(!same)], na.rm = TRUE))

# ---- 4. Stratum signatures (IVW across supertypes; ranking only) ------------
sig <- de |>
  group_by(stratum, gene) |>
  summarise(n_st = n(),
            est_ivw = sum(estimate / se^2) / sum(1 / se^2),
            se_ivw  = sqrt(1 / sum(1 / se^2)),
            .groups = "drop") |>
  filter(n_st >= 2) |>
  mutate(z = est_ivw / se_ivw)
write_csv(sig, file.path(OUT_DIR, "stratum_gene_signatures.csv"))
say("Stratum signatures: %s",
    paste(capture.output(print(count(sig, stratum)))[-(1:3)], collapse = "; "))

# gene-level cross-stratum scatters
sw <- sig |> select(gene, stratum, z) |>
  pivot_wider(names_from = stratum, values_from = z) |> drop_na()
pair_plot <- function(xv, yv) {
  r  <- cor(sw[[xv]], sw[[yv]])
  rs <- cor(sw[[xv]], sw[[yv]], method = "spearman")
  hl <- sw |> filter(gene %in% c("SST", "HCN1", "CALB1", "GRIN2A", "BDNF"))
  ggplot(sw, aes(.data[[xv]], .data[[yv]])) +
    geom_point(alpha = 0.15, size = 0.7) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_vline(xintercept = 0, colour = "grey80") +
    geom_point(data = hl, colour = "#D55E00", size = 2.2) +
    ggrepel::geom_text_repel(data = hl, aes(label = gene), colour = "#D55E00",
                             size = 5, fontface = "italic") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 5.5,
             label = sprintf("r = %.2f, rho = %.2f, n = %s genes",
                             r, rs, format(nrow(sw), big.mark = ","))) +
    labs(x = paste("z,", xv), y = paste("z,", yv)) +
    theme_cowplot(font_size = FS_LAB)
}
p_sc <- plot_grid(pair_plot("depleted", "non_depleted"),
                  pair_plot("depleted", "intermediate"),
                  pair_plot("intermediate", "non_depleted"), nrow = 1)
ggsave(file.path(OUT_DIR, "fig2_stratum_gene_scatter.png"), p_sc,
       width = 18, height = 6.2, dpi = 200, bg = "white")
say("Gene-level z correlations: dep~non = %.2f, dep~int = %.2f, int~non = %.2f",
    cor(sw$depleted, sw$non_depleted), cor(sw$depleted, sw$intermediate),
    cor(sw$intermediate, sw$non_depleted))

# ---- 5. Gene sets (house collections, script 02 convention) -----------------
say("Pulling MSigDB Hallmark + GO(BP/CC/MF) + Reactome ...")
collect <- function(cat, subcat = NULL) {
  args <- list(species = "Homo sapiens", collection = cat)
  if (!is.null(subcat)) args$subcollection <- subcat
  do.call(msigdbr, args) |>
    transmute(gs_name, gene = gene_symbol, source = paste(cat, subcat, sep = ":"))
}
all_gs <- bind_rows(collect("H"), collect("C5", "GO:BP"), collect("C5", "GO:CC"),
                    collect("C5", "GO:MF"), collect("C2", "CP:REACTOME")) |>
  distinct(gs_name, gene, .keep_all = TRUE)
gs_source <- all_gs |> distinct(gs_name, source)
gene_sets <- split(all_gs$gene, all_gs$gs_name)
say("Built %d gene sets", length(gene_sets))

run_gsea <- function(genes, zs, label) {
  d <- tibble(gene = genes, z = zs) |>
    filter(!is.na(z)) |>
    arrange(desc(abs(z))) |> distinct(gene, .keep_all = TRUE)
  ranks <- sort(setNames(d$z, d$gene), decreasing = TRUE)
  say("fgsea: %s (%d genes)", label, length(ranks))
  fgsea(gene_sets, ranks, minSize = MIN_SIZE, maxSize = MAX_SIZE,
        nPermSimple = 10000) |>
    as_tibble() |> mutate(signature = label)
}

# ---- 6. GSEA: three strata + Sst subclass anchor + Nicole's groups ----------
gsea_strata <- sig |>
  group_by(stratum) |>
  group_map(~ run_gsea(.x$gene, .x$z, .y$stratum)) |>
  bind_rows()

sst_sub <- read_csv(SUBCLASS, show_col_types = FALSE) |>
  filter(cell_type == "Sst", !is.na(estimate), !is.na(se), se > 0)
gsea_sub <- run_gsea(sst_sub$genes, sst_sub$estimate / sst_sub$se, "Sst_subclass")

nicole_files <- c(Unaffected = "Unaffected_DE.csv",
                  Vulnerable_AD = "Vulnerable_DE_AD.csv",
                  Vulnerable_SCZAD = "Vulnerable_DE_SCZAD.csv")
nicole <- imap(nicole_files, function(f, nm) {
  read_csv(file.path(DE_DIR, f), show_col_types = FALSE) |>
    filter(!is.na(zval)) |>
    transmute(gene = Gene, z = zval, group = paste0("nicole_", nm))
})
gsea_nicole <- map(nicole, ~ run_gsea(.x$gene, .x$z, .x$group[1])) |> bind_rows()

# forensic: which supertypes does each of Nicole's groups resemble?
say("Correlation of Nicole's group z with each supertype z (top 6 each):")
for (nm in names(nicole)) {
  m <- nicole[[nm]] |> inner_join(de |> select(gene, supertype, zval), by = "gene") |>
    group_by(supertype) |> summarise(r = cor(z, zval, use = "complete.obs")) |>
    arrange(desc(r)) |> head(6)
  say("  %s: %s", nm, paste(sprintf("%s %.2f", m$supertype, m$r), collapse = ", "))
}

gsea_all <- bind_rows(gsea_strata, gsea_sub, gsea_nicole) |>
  left_join(gs_source, by = c(pathway = "gs_name"))
write_csv(gsea_all |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = "|")),
          file.path(OUT_DIR, "gsea_all_signatures.csv"))

# ---- 7. Shared vs stratum-specific pathways ---------------------------------
wide <- gsea_all |>
  filter(signature %in% c("depleted", "intermediate", "non_depleted")) |>
  select(pathway, source, signature, NES, padj) |>
  pivot_wider(names_from = signature, values_from = c(NES, padj)) |>
  drop_na(NES_depleted, NES_non_depleted) |>
  mutate(category = case_when(
    padj_depleted < 0.05 & padj_non_depleted < 0.05 &
      sign(NES_depleted) == sign(NES_non_depleted) ~ "shared",
    padj_depleted < 0.05 & padj_non_depleted > 0.20 &
      abs(NES_depleted - NES_non_depleted) > 1     ~ "depleted_specific",
    padj_non_depleted < 0.05 & padj_depleted > 0.20 &
      abs(NES_depleted - NES_non_depleted) > 1     ~ "non_depleted_specific",
    TRUE                                           ~ "other"))
write_csv(wide, file.path(OUT_DIR, "gsea_stratum_comparison.csv"))

say("NES correlation depleted vs non-depleted (all %d pathways): r = %.2f",
    nrow(wide), cor(wide$NES_depleted, wide$NES_non_depleted))
say("Category counts:")
print(count(wide, category))
for (cc in c("shared", "depleted_specific", "non_depleted_specific")) {
  say("--- top %s (by depleted padj) ---", cc)
  print(wide |> filter(category == cc) |>
          arrange(pmin(padj_depleted, padj_non_depleted)) |>
          select(pathway, NES_depleted, NES_intermediate, NES_non_depleted,
                 padj_depleted, padj_non_depleted) |> head(15), width = 200)
}

lab_set <- wide |>
  filter(category != "other") |>
  group_by(category) |>
  slice_min(pmin(padj_depleted, padj_non_depleted), n = 6) |> ungroup()
p_nes <- ggplot(wide, aes(NES_depleted, NES_non_depleted)) +
  geom_point(data = filter(wide, category == "other"), colour = "grey75",
             alpha = 0.4, size = 0.9) +
  geom_point(data = filter(wide, category != "other"), aes(colour = category), size = 1.8) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  ggrepel::geom_text_repel(data = lab_set,
                           aes(label = str_trunc(pathway, 42), colour = category),
                           size = 3.6, max.overlaps = 25, show.legend = FALSE) +
  scale_colour_manual(values = c(shared = "#009E73", depleted_specific = "#D55E00",
                                 non_depleted_specific = "#0072B2")) +
  labs(x = "NES, depleted Sst supertypes", y = "NES, non-depleted Sst supertypes",
       colour = NULL,
       title = "Pathway dysregulation: depleted vs non-depleted Sst strata") +
  theme_cowplot(font_size = FS_LAB) +
  theme(legend.position = "bottom", plot.title = element_text(size = FS_TITLE))
ggsave(file.path(OUT_DIR, "fig3_nes_depleted_vs_nondepleted.png"), p_nes,
       width = 11, height = 10, dpi = 200, bg = "white")

# ---- 8. Focused-theme heatmap (paper themes, script 02) ---------------------
focus_terms <- list(
  synaptic = c("GOCC_SYNAPTIC_VESICLE", "GOBP_SYNAPTIC_VESICLE_CYCLE",
               "GOBP_NEUROTRANSMITTER_SECRETION", "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING",
               "GOBP_CHEMICAL_SYNAPTIC_TRANSMISSION", "GOCC_PRESYNAPSE", "GOCC_POSTSYNAPSE",
               "REACTOME_NEURONAL_SYSTEM", "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES"),
  mitochondrial = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
                    "GOCC_MITOCHONDRIAL_INNER_MEMBRANE", "GOBP_OXIDATIVE_PHOSPHORYLATION",
                    "GOBP_MITOCHONDRIAL_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT"),
  cholesterol = c("HALLMARK_CHOLESTEROL_HOMEOSTASIS", "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
                  "GOBP_CHOLESTEROL_BIOSYNTHETIC_PROCESS", "GOBP_STEROL_BIOSYNTHETIC_PROCESS"),
  bmp_tgfb = c("GOBP_BMP_SIGNALING_PATHWAY", "HALLMARK_TGF_BETA_SIGNALING",
               "REACTOME_SIGNALING_BY_BMP"))
sig_order <- c("Sst_subclass", "depleted", "intermediate", "non_depleted")
d_focus <- gsea_all |>
  filter(pathway %in% unlist(focus_terms), signature %in% sig_order) |>
  mutate(theme = names(unlist(map(focus_terms, ~ .x)))[match(pathway, unlist(focus_terms))],
         theme = str_remove(theme, "\\d+$"),
         signature = factor(signature, levels = sig_order),
         pathway = factor(pathway, levels = rev(unlist(focus_terms))))
p_heat <- ggplot(d_focus, aes(signature, pathway, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(padj < 0.01, "**", ifelse(padj < 0.05, "*",
                        ifelse(padj < 0.1, "+", "")))), size = 5, vjust = 0.75) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
                       limits = c(-3.5, 3.5), oob = scales::squish, name = "NES") +
  labs(x = NULL, y = NULL, title = "Paper pathway themes across Sst strata",
       subtitle = "** padj<0.01, * padj<0.05, + padj<0.1") +
  theme_cowplot(font_size = FS_LAB) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(size = FS_TITLE))
ggsave(file.path(OUT_DIR, "fig4_theme_heatmap.png"), p_heat,
       width = 12, height = 10.5, dpi = 200, bg = "white")

say("Done. Outputs in %s", OUT_DIR)
