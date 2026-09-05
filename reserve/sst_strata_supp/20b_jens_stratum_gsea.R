#!/usr/bin/env Rscript
# Independent-cohort replication of the Sst-strata GSEA (Jens dataset; not in
# the 7-cohort meta-analysis). Pseudobulk limma-voom per stratum (SCZ vs CTRL,
# + sex, age, PMI), preranked fgsea with the house collections, then compare
# the top-level conclusions against the meta-analytic (IVW) strata results:
# (1) burden gradient depleted > intermediate > non-depleted, down-dominated;
# (2) shared synaptic/SST/VGF component; (3) graded translation + OxPhos.
# Sensitivity: supertype_confidence >= 0.5 cells only (label-quality concern:
# 4/5 depleted supertypes flagged CONFIDENTLY_WRONG in this cohort's mapping).
# Run from repo root: Rscript reserve/sst_strata_supp/20b_jens_stratum_gsea.R
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(tibble); library(arrow); library(edgeR); library(limma)
  library(fgsea); library(msigdbr); library(ggplot2); library(cowplot)
})
set.seed(42)
OUT  <- "transcriptomic/results/sst_strata_gsea"
JOUT <- file.path(OUT, "jens")
FS <- 14
MIN_CELLS <- 10
STRATA <- c("depleted", "intermediate", "non_depleted")

t0 <- Sys.time()
say <- function(...) cat(sprintf("[%5.1f min] ", as.numeric(difftime(Sys.time(), t0, units="mins"))), sprintf(...), "\n", sep = "")

load_pb <- function(tag) {
  pb <- as.data.frame(read_parquet(file.path(JOUT, sprintf("jens_stratum_pseudobulk%s.parquet", tag))))
  idx <- grep("__index_level_0__", colnames(pb), value = TRUE)
  rownames(pb) <- pb[[idx]]; pb[[idx]] <- NULL
  meta <- read_csv(file.path(JOUT, sprintf("jens_stratum_pseudobulk_meta%s.csv", tag)),
                   show_col_types = FALSE)
  list(pb = as.matrix(pb), meta = meta)
}

de_stratum <- function(pb, meta, keys) {
  m <- meta |> filter(key %in% keys, n_cells >= MIN_CELLS) |>
    mutate(dx = factor(diagnosis, c("CTRL", "SCZ")), sex = factor(sex),
           age_s = scale(ifelse(is.na(age), median(age, na.rm = TRUE), age))[, 1],
           pmi_s = scale(ifelse(is.na(pmi), median(pmi, na.rm = TRUE), pmi))[, 1]) |>
    filter(!is.na(dx), !is.na(sex))
  say("  n groups used: %d (%s)", nrow(m),
      paste(capture.output(table(m$dx)), collapse = " "))
  y <- DGEList(pb[, m$key])
  keep <- filterByExpr(y, group = m$diagnosis)
  y <- calcNormFactors(y[keep, ])
  des <- model.matrix(~ dx + sex + age_s + pmi_s, data = m)
  stopifnot(nrow(des) == ncol(y))
  v <- voom(y, des)
  fit <- eBayes(lmFit(v, des))
  topTable(fit, coef = "dxSCZ", number = Inf, sort.by = "none") |>
    rownames_to_column("gene") |> as_tibble()
}

# ---- gene sets (house collections) -------------------------------------------
say("building gene sets ...")
collect <- function(cat, subcat = NULL) {
  args <- list(species = "Homo sapiens", collection = cat)
  if (!is.null(subcat)) args$subcollection <- subcat
  do.call(msigdbr, args) |> transmute(gs_name, gene = gene_symbol)
}
all_gs <- bind_rows(collect("H"), collect("C5", "GO:BP"), collect("C5", "GO:CC"),
                    collect("C5", "GO:MF"), collect("C2", "CP:REACTOME")) |>
  distinct(gs_name, gene)
gene_sets <- split(all_gs$gene, all_gs$gs_name)

run_gsea <- function(de, label) {
  d <- de |> filter(!is.na(t)) |> arrange(desc(abs(t))) |> distinct(gene, .keep_all = TRUE)
  ranks <- sort(setNames(d$t, d$gene), decreasing = TRUE)
  say("fgsea: %s (%d genes)", label, length(ranks))
  fgsea(gene_sets, ranks, minSize = 10, maxSize = 500, nPermSimple = 10000) |>
    as_tibble() |> mutate(signature = label)
}

# ---- run both variants ---------------------------------------------------------
res_all <- list(); de_all <- list()
for (tag in c("", "_conf05")) {
  d <- load_pb(tag)
  for (st in STRATA) {
    keys <- d$meta |> filter(stratum == st) |> pull(key)
    de <- de_stratum(d$pb, d$meta, keys)
    lab <- paste0("jens_", st, tag)
    de_all[[lab]] <- de |> mutate(signature = lab)
    res_all[[lab]] <- run_gsea(de, lab)
  }
  if (tag == "") {   # pooled all-Sst for the dilution comparison
    agg <- d$meta |> mutate(donor = donor_id)
    pb_don <- sapply(split(agg$key, agg$donor), function(k)
      rowSums(d$pb[, k, drop = FALSE]))
    meta_don <- agg |> group_by(donor_id) |>
      summarise(n_cells = sum(n_cells), diagnosis = first(diagnosis),
                sex = first(sex), age = first(age), pmi = first(pmi)) |>
      mutate(key = donor_id, stratum = "all")
    de <- de_stratum(pb_don, meta_don, meta_don$key)
    de_all[["jens_all_sst"]] <- de |> mutate(signature = "jens_all_sst")
    res_all[["jens_all_sst"]] <- run_gsea(de, "jens_all_sst")
  }
}
gsea_j <- bind_rows(res_all)
write_csv(gsea_j |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = "|")),
          file.path(JOUT, "gsea_jens.csv"))
write_csv(bind_rows(de_all), file.path(JOUT, "de_jens_strata.csv"))

# ---- (1) burden ---------------------------------------------------------------
say("=== burden (padj<0.05) ===")
burden <- gsea_j |> filter(padj < 0.05) |>
  count(signature, dir = ifelse(NES > 0, "up", "down")) |>
  pivot_wider(names_from = dir, values_from = n, values_fill = 0)
print(burden)

# ---- (2)+(3): key blocks vs the meta results -----------------------------------
g_meta <- read_csv(file.path(OUT, "gsea_all_signatures.csv"), show_col_types = FALSE) |>
  filter(signature %in% STRATA) |>
  transmute(pathway, signature = paste0("meta_", signature), NES, padj)
key_paths <- c("GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING", "GOCC_PRESYNAPSE",
               "GOBP_NEUROTRANSMITTER_SECRETION", "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES",
               "REACTOME_TRANSLATION", "GOCC_RIBOSOMAL_SUBUNIT", "GOBP_CYTOPLASMIC_TRANSLATION",
               "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
               "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "GOBP_OXIDATIVE_PHOSPHORYLATION",
               "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT")
comb <- bind_rows(g_meta,
                  gsea_j |> filter(signature %in% paste0("jens_", STRATA)) |>
                    transmute(pathway, signature, NES, padj))
say("=== key pathway blocks, meta vs Jens (NES [stars]) ===")
tabkey <- comb |> filter(pathway %in% key_paths) |>
  mutate(cell = sprintf("%.1f%s", NES, ifelse(padj < 0.01, "**",
                 ifelse(padj < 0.05, "*", ifelse(padj < 0.1, "+", ""))))) |>
  select(pathway, signature, cell) |>
  pivot_wider(names_from = signature, values_from = cell) |>
  select(pathway, meta_depleted, jens_depleted, meta_intermediate,
         jens_intermediate, meta_non_depleted, jens_non_depleted)
print(tabkey, n = 15, width = 250)
write_csv(tabkey, file.path(JOUT, "key_blocks_meta_vs_jens.csv"))

# NES correlations across all common pathways, per stratum
say("=== cross-cohort NES correlations (all pathways | meta padj<0.05 pathways) ===")
for (st in STRATA) {
  m <- g_meta |> filter(signature == paste0("meta_", st)) |>
    inner_join(gsea_j |> filter(signature == paste0("jens_", st)) |>
                 select(pathway, NES_j = NES), by = "pathway")
  say("%s: r = %.2f (n=%d) | r = %.2f among meta-sig (n=%d)", st,
      cor(m$NES, m$NES_j), nrow(m),
      cor(m$NES[m$padj < 0.05], m$NES_j[m$padj < 0.05]), sum(m$padj < 0.05))
}

# ---- gene-level checks ----------------------------------------------------------
say("=== key genes, Jens t by stratum (primary) ===")
print(bind_rows(de_all) |>
        filter(gene %in% c("SST", "VGF", "CALB1", "RASGRF2", "NTRK2", "GAD2"),
               !str_detect(signature, "conf05|all_sst")) |>
        select(gene, signature, t) |>
        mutate(t = round(t, 2)) |>
        pivot_wider(names_from = signature, values_from = t))

# ---- comparison heatmap ----------------------------------------------------------
lab_map <- read_csv(file.path(OUT, "gsea_all_signatures.csv"), show_col_types = FALSE) |>
  distinct(pathway) # placeholder no-op to keep pathway naming consistent
d <- comb |> filter(pathway %in% key_paths) |>
  mutate(block = case_when(
    str_detect(pathway, "SYNAP|PRESYN|NEUROTRANS|TRANSMISSION") ~ "Shared: synaptic",
    str_detect(pathway, "TRANSLAT|RIBOSOM") ~ "Graded: translation",
    TRUE ~ "Graded: OxPhos"),
    cohort = ifelse(str_starts(signature, "meta_"), "Meta (7 cohorts, IVW)", "Jens (pseudobulk)"),
    stratum = str_remove(signature, "^(meta|jens)_"),
    stratum = factor(stratum, STRATA, c("Depleted", "Intermediate", "Non-depleted")),
    stars = ifelse(padj < 0.01, "**", ifelse(padj < 0.05, "*", ifelse(padj < 0.1, "+", ""))))
d$pathway <- factor(d$pathway, levels = rev(key_paths))
p <- ggplot(d, aes(stratum, pathway, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = stars), size = 5.5, vjust = 0.75) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
                       limits = c(-3, 3), oob = scales::squish, name = "NES") +
  facet_grid(block ~ cohort, scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL,
       title = "Sst-strata pathway replication: 7-cohort meta vs independent Jens cohort") +
  theme_cowplot(font_size = FS) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "grey92"),
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
        axis.ticks = element_blank(), axis.line = element_blank(),
        plot.title = element_text(size = FS + 2))
ggsave(file.path(JOUT, "fig_jens_replication.png"), p, width = 13.5, height = 7,
       dpi = 200, bg = "white")
say("wrote %s", file.path(JOUT, "fig_jens_replication.png"))
