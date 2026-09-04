#!/usr/bin/env Rscript
# Sensitivity for the Jens replication: the mapping audit shows depleted <->
# intermediate label porosity (runner-up crosses that boundary ~30-40% of the
# time) while depleted <-> non-depleted confusion is ~0. So test the boundary
# the labels CAN draw: merged "affected" (depleted + intermediate) vs
# non-depleted. Also print the key-block NES for the conf05 variant.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(tibble)
  library(arrow); library(edgeR); library(limma); library(fgsea); library(msigdbr)
})
set.seed(42)
OUT  <- "transcriptomic/results/sst_strata_gsea"
JOUT <- file.path(OUT, "jens")
MIN_CELLS <- 10
KEY <- c("GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING", "GOBP_NEUROTRANSMITTER_SECRETION",
         "REACTOME_TRANSLATION", "GOCC_RIBOSOMAL_SUBUNIT", "GOBP_CYTOPLASMIC_TRANSLATION",
         "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
         "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "GOBP_OXIDATIVE_PHOSPHORYLATION",
         "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT")

pb <- as.data.frame(read_parquet(file.path(JOUT, "jens_stratum_pseudobulk.parquet")))
idx <- grep("__index_level_0__", colnames(pb), value = TRUE)
rownames(pb) <- pb[[idx]]; pb[[idx]] <- NULL; pb <- as.matrix(pb)
meta <- read_csv(file.path(JOUT, "jens_stratum_pseudobulk_meta.csv"), show_col_types = FALSE)

# merge depleted + intermediate pseudobulks per donor
mk_group <- function(strata, label) {
  m <- meta |> filter(stratum %in% strata) |> group_by(donor_id) |>
    summarise(n_cells = sum(n_cells), diagnosis = first(diagnosis),
              sex = first(sex), age = first(age), pmi = first(pmi), .groups = "drop") |>
    mutate(key = paste0(donor_id, "|", label))
  cols <- sapply(m$donor_id, function(d) {
    k <- meta |> filter(stratum %in% strata, donor_id == d) |> pull(key)
    rowSums(pb[, k, drop = FALSE])
  })
  colnames(cols) <- m$key
  list(pb = cols, meta = m)
}
aff <- mk_group(c("depleted", "intermediate"), "affected")
non <- mk_group("non_depleted", "non_depleted")

de_run <- function(g) {
  m <- g$meta |> filter(n_cells >= MIN_CELLS) |>
    mutate(dx = factor(diagnosis, c("CTRL", "SCZ")), sex = factor(sex),
           age_s = scale(ifelse(is.na(age), median(age, na.rm = TRUE), age))[, 1],
           pmi_s = scale(ifelse(is.na(pmi), median(pmi, na.rm = TRUE), pmi))[, 1])
  y <- DGEList(g$pb[, m$key])
  y <- calcNormFactors(y[filterByExpr(y, group = m$diagnosis), ])
  des <- model.matrix(~ dx + sex + age_s + pmi_s, data = m)
  fit <- eBayes(lmFit(voom(y, des), des))
  topTable(fit, coef = "dxSCZ", number = Inf, sort.by = "none") |>
    rownames_to_column("gene") |> as_tibble()
}
cat("DE: affected (dep+int) and non-depleted ...\n")
de_aff <- de_run(aff); de_non <- de_run(non)

collect <- function(cat, subcat = NULL) {
  args <- list(species = "Homo sapiens", collection = cat)
  if (!is.null(subcat)) args$subcollection <- subcat
  do.call(msigdbr, args) |> transmute(gs_name, gene = gene_symbol)
}
gs <- bind_rows(collect("H"), collect("C5", "GO:BP"), collect("C5", "GO:CC"),
                collect("C5", "GO:MF"), collect("C2", "CP:REACTOME")) |> distinct()
gsl <- split(gs$gene, gs$gs_name)
run_g <- function(de, lab) {
  d <- de |> filter(!is.na(t)) |> arrange(desc(abs(t))) |> distinct(gene, .keep_all = TRUE)
  fgsea(gsl, sort(setNames(d$t, d$gene), decreasing = TRUE),
        minSize = 10, maxSize = 500, nPermSimple = 10000) |>
    as_tibble() |> mutate(signature = lab)
}
g <- bind_rows(run_g(de_aff, "jens_affected"), run_g(de_non, "jens_non_depleted_v2"))
write_csv(g |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = "|")),
          file.path(JOUT, "gsea_jens_affected.csv"))

cat("\nBurden (padj<0.05):\n")
print(g |> filter(padj < 0.05) |> count(signature, dir = ifelse(NES > 0, "up", "down")) |>
        pivot_wider(names_from = dir, values_from = n, values_fill = 0))

fmt <- function(x) sprintf("%.1f%s", x$NES, ifelse(x$padj < 0.01, "**",
                    ifelse(x$padj < 0.05, "*", ifelse(x$padj < 0.1, "+", ""))))
cat("\nKey blocks, affected vs non-depleted:\n")
tab <- g |> filter(pathway %in% KEY) |>
  mutate(cell = fmt(pick(NES, padj))) |>
  select(pathway, signature, cell) |> pivot_wider(names_from = signature, values_from = cell)
print(tab, n = 12)

# conf05 key blocks from the main run
gj <- read_csv(file.path(JOUT, "gsea_jens.csv"), show_col_types = FALSE)
cat("\nKey blocks, conf05 strata (from main run):\n")
print(gj |> filter(pathway %in% KEY, str_detect(signature, "conf05")) |>
        mutate(cell = fmt(pick(NES, padj))) |>
        select(pathway, signature, cell) |>
        pivot_wider(names_from = signature, values_from = cell), n = 12)

# gene-level checks in the merged contrast
cat("\nKey genes (t), affected vs non-depleted:\n")
print(bind_rows(de_aff |> mutate(sig = "affected"), de_non |> mutate(sig = "non_depleted")) |>
        filter(gene %in% c("SST", "VGF", "CALB1", "NTRK2", "RASGRF2")) |>
        select(gene, sig, t) |> mutate(t = round(t, 2)) |>
        pivot_wider(names_from = sig, values_from = t))
