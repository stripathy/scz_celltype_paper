#!/usr/bin/env Rscript
# Shared constants and helpers for the Figure 5 pipeline (Sst depletion strata).
#
# WHY this file exists. Before consolidation, the module definitions (the
# leading-edge unions behind "translation", "OxPhos/mito", "synaptic",
# "deubiquitination") were re-derived independently in nine scripts, the strata
# membership was hard-coded in three, and the MSigDB collections were rebuilt in
# nine. Any edit to a block's gene sets silently desynchronised the figure from
# the statistics quoted in the text -- the bug class that produced the mislabeled
# panel c. Everything downstream now imports from here.
#
# All paths assume the repo root as the working directory.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(purrr); library(tibble)
})

# ---- paths -------------------------------------------------------------------
P <- list(
  out   = "transcriptomic/results/sst_strata_gsea",              # canonical results
  pb    = "transcriptomic/results/sst_strata_gsea/pseudobulk",    # pseudobulk DE/GSEA
  donor = "transcriptomic/results/sst_strata_gsea/pseudobulk/donor_stratum",
  supp  = "transcriptomic/results/sst_strata_gsea/supp",
  exp   = "transcriptomic/data/stratum_pseudobulks_export",       # 7-cohort export
  arch  = "transcriptomic/results/sst_strata_gsea/archive_ivw",   # superseded IVW run
  cache = "transcriptomic/results/sst_strata_gsea/.cache",
  # Submission locations, both tracked in git. Main figures live in the module's
  # results/figures (as Fig 2 does: transcriptomic/results/figures/09_composite);
  # supplements are written DIRECTLY into the submission folder as originals, per
  # manuscript/figures/supplementary/README.md -- no copy step, nothing to drift.
  fig     = "transcriptomic/results/figures",
  suppfig = "manuscript/figures/supplementary",
  # external inputs. NOTE: crumblr here is the SCZ 7-cohort meta (Fig 3a). The
  # similarly-named crossdisorder/results/crumblr_results_supertype_neurons.csv
  # is the AD/CPS model -- not interchangeable.
  crumblr  = "shared/snrnaseq_de/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv",
  depth    = "spatial/output/depth_platform/supertype_depth_platform_summary.csv",
  layers   = "spatial/output/depth_proportions/proposed_layer_boundaries.csv",
  # The canonical copy lives in the git-ignored seam; the committed snapshot under
  # transcriptomic/data/figure_inputs/ was captured from that exact path (see its
  # MANIFEST.tsv, md5 719da3d7519cbc5a3aec1b8c42d438da), so falling back to it lets
  # S8 render from a clean clone without changing which numbers are drawn.
  subclass = local({
    seam <- "shared/snrnaseq_de/DE_genes_all_cells_scz.csv"
    snap <- "transcriptomic/data/figure_inputs/DE_genes_all_cells_scz.csv"
    if (file.exists(seam)) seam else snap
  }),
  spatial_de = "spatial/output/de/de_results_supertype.csv")

for (d in c(P$pb, P$donor, P$supp, P$cache, P$fig, P$suppfig)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- progress ----------------------------------------------------------------
.t0 <- Sys.time()
say <- function(...) cat(sprintf("[%6.1f min] ", as.numeric(difftime(Sys.time(), .t0, units = "mins"))),
                         sprintf(...), "\n", sep = "")

# ---- strata ------------------------------------------------------------------
STRATA_LEVELS <- c("depleted", "intermediate", "non_depleted")
SIG_LEVELS    <- c("Sst_subclass", STRATA_LEVELS)
SIG_SHORT     <- c("All Sst", "Depleted", "Intermediate", "Non-depleted")
MIN_CELLS     <- 10          # donors need >= this many nuclei of a stratum

# Read the stratum assignment produced by 01_strata.R. Membership lives in that
# one CSV so the definition cannot drift between the DE, the interaction model
# and the figure.
load_strata <- function() {
  read_csv(file.path(P$out, "strata_definition.csv"), show_col_types = FALSE)
}
# supertype -> stratum, as a named character vector
stratum_of <- function(strata = load_strata()) {
  setNames(strata$stratum, strata$CellType)
}
# "Depleted (5)" style labels, counts derived from the table rather than typed in
strata_labels <- function(strata = load_strata()) {
  n <- table(factor(strata$stratum, STRATA_LEVELS))
  setNames(sprintf("%s (%d)", c("Depleted", "Intermediate", "Non-depleted"), as.integer(n)),
           STRATA_LEVELS)
}

# ---- gene-set blocks (the canonical module definition) -----------------------
# One row per gene set shown in the figure: which block it belongs to and the
# display label. modules() derives the gene lists from this table, so adding a
# row to a block automatically updates every panel and every module statistic.
BLOCKS <- tribble(
  ~pathway,                                         ~block,   ~label,
  "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING",    "synaptic",    "Regulation of trans-synaptic signaling (GO:BP)",
  "GOCC_SYNAPTIC_MEMBRANE",                         "synaptic",    "Synaptic membrane (GO:CC)",
  "GOBP_NEUROTRANSMITTER_SECRETION",                "synaptic",    "Neurotransmitter secretion (GO:BP)",
  "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES", "synaptic",    "Transmission across chemical synapses (Reactome)",
  "GOBP_PROTEIN_K63_LINKED_DEUBIQUITINATION",       "ubiquitin",   "Protein K63-linked deubiquitination (GO:BP)",
  "GOMF_UBIQUITIN_CONJUGATING_ENZYME_BINDING",      "ubiquitin",   "Ubiquitin-conjugating enzyme binding (GO:MF)",
  "REACTOME_TRANSLATION",                           "translation", "Translation (Reactome)",
  "GOCC_RIBOSOMAL_SUBUNIT",                         "translation", "Ribosomal subunit (GO:CC)",
  "GOBP_CYTOPLASMIC_TRANSLATION",                   "translation", "Cytoplasmic translation (GO:BP)",
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",     "translation", "Eukaryotic translation elongation (Reactome)",
  "GOBP_OXIDATIVE_PHOSPHORYLATION",                 "oxphos",      "Oxidative phosphorylation (GO:BP)",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",        "oxphos",      "Respiratory electron transport (Reactome)")

BLOCK_KEYS   <- c("synaptic", "ubiquitin", "translation", "oxphos")
BLOCK_LABELS <- c(synaptic    = "Shared: synaptic program",
                  ubiquitin   = "Shared: ubiquitin",
                  translation = "Graded: cytosolic translation",
                  oxphos      = "Graded: oxidative phosphorylation")

# Leading-edge union per block, taken from the depleted-stratum GSEA (the
# stratum that defines the modules). Returns a named list keyed by BLOCK_KEYS.
# Accepts leadingEdge either as fgsea's in-memory list-column or as the
# "|"-joined string written to CSV, so it works before and after the write.
modules <- function(gsea, signature = "depleted", keys = BLOCK_KEYS) {
  le_union <- function(pws) {
    le <- gsea |> filter(pathway %in% pws, signature == !!signature) |> pull(leadingEdge)
    if (is.list(le)) unique(unlist(le)) else unique(unlist(strsplit(le, "|", fixed = TRUE)))
  }
  lapply(setNames(keys, keys), function(k) le_union(BLOCKS$pathway[BLOCKS$block == k]))
}
# gene -> block key, for colouring. CONTRACT: a gene in several modules takes the
# EARLIEST module in `mods`, so pass `mods` in priority order (the assignment loop
# runs in reverse so earlier names overwrite later ones). Panel c relies on this:
# ubiquitin > translation > oxphos > synaptic.
module_of <- function(genes, mods) {
  out <- rep("other", length(genes))
  for (k in rev(names(mods))) out[genes %in% mods[[k]]] <- k
  out
}
MODULE_PRIORITY <- c("ubiquitin", "translation", "oxphos", "synaptic")

# ---- MSigDB gene sets (built once, cached) -----------------------------------
# msigdbr + the collection bind takes ~30 s; it was previously repeated in every
# script that ran fgsea.
msigdb_sets <- function(refresh = FALSE) {
  f <- file.path(P$cache, "msigdb_sets.rds")
  if (!refresh && file.exists(f)) return(readRDS(f))
  suppressPackageStartupMessages(library(msigdbr))
  collect <- function(cat, subcat = NULL) {
    args <- list(species = "Homo sapiens", collection = cat)
    if (!is.null(subcat)) args$subcollection <- subcat
    do.call(msigdbr, args) |>
      transmute(gs_name, gene = gene_symbol, source = paste(cat, subcat, sep = ":"))
  }
  # GO (BP/CC/MF) + Reactome. Hallmark is deliberately excluded: its sets are
  # coarse meta-signatures assembled across pathways, which makes them hard to
  # interpret next to the named GO and Reactome pathways the figure reports.
  all_gs <- bind_rows(collect("C5", "GO:BP"), collect("C5", "GO:CC"),
                      collect("C5", "GO:MF"), collect("C2", "CP:REACTOME")) |>
    distinct(gs_name, gene, .keep_all = TRUE)
  res <- list(sets = split(all_gs$gene, all_gs$gs_name),
              source = all_gs |> distinct(gs_name, source))
  saveRDS(res, f)
  res
}

# House GSEA settings, in one place (matches the subclass pipeline).
GSEA_MIN <- 10; GSEA_MAX <- 500; GSEA_NPERM <- 10000
run_fgsea <- function(genes, z, label, sets = NULL) {
  suppressPackageStartupMessages(library(fgsea))
  if (is.null(sets)) sets <- msigdb_sets()$sets
  d <- tibble(gene = genes, z = z) |> filter(!is.na(z)) |>
    arrange(desc(abs(z))) |> distinct(gene, .keep_all = TRUE)
  ranks <- sort(setNames(d$z, d$gene), decreasing = TRUE)
  say("fgsea: %s (%d genes)", label, length(ranks))
  fgsea(sets, ranks, minSize = GSEA_MIN, maxSize = GSEA_MAX,
        nPermSimple = GSEA_NPERM) |> as_tibble() |> mutate(signature = label)
}

# ---- meta-analysis -----------------------------------------------------------
# Random-effects (REML) inverse-variance meta-analysis, DL fallback, as in Fig 2.
fit_meta <- function(yi, sei, extra = FALSE) {
  suppressPackageStartupMessages(library(metafor))
  fit <- tryCatch(metafor::rma(yi = yi, sei = sei, method = "REML",
                               control = list(maxiter = 1000)),
                  error = function(e) tryCatch(metafor::rma(yi = yi, sei = sei, method = "DL"),
                                               error = function(e2) NULL))
  if (is.null(fit)) return(NULL)
  out <- tibble(estimate = fit$beta[1], se = fit$se, zval = fit$zval,
                pval = fit$pval, k = fit$k)
  if (extra) out <- bind_cols(out, tibble(ci.lb = fit$ci.lb, ci.ub = fit$ci.ub,
                                          tau2 = fit$tau2, I2 = fit$I2))
  out
}

# ---- pseudobulk loading ------------------------------------------------------
# Per-donor x stratum counts written by 04_donor_pseudobulks.R. Arrow returns
# int64 columns that break rowsum()/DGEList(); coerce to numeric on load.
load_donor_stratum <- function(cohort) {
  suppressPackageStartupMessages(library(arrow))
  pb <- as.data.frame(arrow::read_parquet(
    file.path(P$donor, paste0(cohort, "_donor_stratum_counts.parquet"))))
  rn <- pb$gene; pb$gene <- NULL
  M <- vapply(pb, as.numeric, numeric(nrow(pb))); rownames(M) <- rn
  meta <- read_csv(file.path(P$donor, paste0(cohort, "_donor_stratum_meta.csv")),
                   col_types = cols(donor = col_character(), .default = col_guess()))
  list(counts = M, meta = meta)
}
donor_cohorts <- function() {
  str_remove(basename(Sys.glob(file.path(P$donor, "*_donor_stratum_meta.csv"))),
             "_donor_stratum_meta\\.csv$")
}
# Donor ids must be read as character: Multiome's are numeric-looking, and a
# bare read_csv types them as double, which breaks any bind_rows across cohorts.
read_donor_meta <- function(cohort) {
  read_csv(file.path(P$donor, paste0(cohort, "_donor_stratum_meta.csv")),
           col_types = cols(donor = col_character(), .default = col_guess()))
}
all_donor_meta <- function() map_dfr(donor_cohorts(), read_donor_meta)
export_cohorts <- function() {
  str_remove(basename(Sys.glob(file.path(P$exp, "*_groups.csv"))), "_groups\\.csv$")
}
# donor-level covariates from the export groups files
load_covars <- function(cohorts = export_cohorts()) {
  map_dfr(cohorts, function(co)
    read_csv(file.path(P$exp, paste0(co, "_groups.csv")),
             col_types = cols(donor = col_character(), .default = col_guess())) |>
      distinct(donor, sex, age, pmi) |>
      mutate(cohort = co,
             age = suppressWarnings(as.numeric(as.character(age))),
             pmi = suppressWarnings(as.numeric(as.character(pmi)))))
}
# Median-impute age/PMI and build the model frame. model.matrix() silently drops
# rows with NA covariates, which desynchronises the design from the count matrix.
prep_model_frame <- function(m) {
  m <- m |> mutate(dx = factor(diagnosis, c("Control", "SCZ")), sex = factor(sex),
                   age_s = scale(ifelse(is.na(age), median(age, na.rm = TRUE), age))[, 1])
  attr(m, "use_pmi") <- mean(is.na(m$pmi)) < 0.5
  if (isTRUE(attr(m, "use_pmi")))
    m$pmi_s <- scale(ifelse(is.na(m$pmi), median(m$pmi, na.rm = TRUE), m$pmi))[, 1]
  m
}

# ---- figure conventions ------------------------------------------------------
# MATCHES FIGURE 4 (genetics/scripts/figures/fig4_style.R) exactly. Figure 4 is
# built at 8.0 in wide and placed in a 7.1 in column, so every text size is its
# nominal value x 7.1/8.0 in print: base 8.0 pt, axis titles 7.5, tick labels 7.0.
# That is deliberately ~1 pt larger than Figure 2, which builds at 7.1 in with a
# 7 pt base. Figure 5 follows Figure 4 rather than Figure 2 so the two figures the
# reader meets last look the same.
FIG_SCALE <- 8.0 / 7.1
FIG_W     <- 8.0     # built here, printed at 7.1 in
FIG_DPI   <- 600     # Nature accepts 300-600; Figure 4 uses 400, this is finer
FS          <- 8 * FIG_SCALE          # 9.01 pt nominal -> 8.0 pt printed
AXIS_TITLE  <- FS - 0.5 * FIG_SCALE   # 8.45 -> 7.5 printed
AXIS_TEXT   <- FS - 1.0 * FIG_SCALE   # 7.89 -> 7.0 printed
LEGEND_TEXT <- FS - 1.5 * FIG_SCALE   # 7.32 -> 6.5 printed
PANEL_LABEL <- 9 * FIG_SCALE          # 10.14 -> 9.0 printed
# In-panel text as multiples of the base, the same ratios fig4_style.R uses.
# ggplot's `size` is millimetres, so pt = size * 2.845.
LBL_SMALL <- FS * 0.28        # 2.52 mm = 7.2 pt nominal
LBL_GENE  <- FS * 0.32        # 2.88 mm = 8.2 pt
LBL_CALL  <- FS * 0.34        # 3.06 mm = 8.7 pt
LBL_STAR  <- FS * 0.43        # significance asterisks
LW        <- 0.25 * FIG_SCALE # axis and tile linewidths, as in fig4_style.R
GEO       <- 0.64             # point-size multiplier, scaled with the base font

# Figure 2 palette (01_volcano_per_celltype.R / 09_composite_figure.R):
# direction = blue/orange, with a light tier for the FDR 0.05-0.10 band
COL_DOWN <- "#0072B2"; COL_UP <- "#D55E00"
DOWN_LIGHT <- "#9FCAE6"; UP_LIGHT <- "#F2B58C"
# Strata: purple -> green, deliberately outside the blue-orange (direction) and
# magenta/teal (module) families so panel a cannot be misread as a direction.
STRAT_COLS <- c(depleted = "#762A83", intermediate = "#C2A5CF", non_depleted = "#1B7837")
# Modules: maximum separation between translation and OxPhos, the two graded blocks
FAM_COLS <- c(synaptic = "#7FA6CC", translation = "#AA3377", oxphos = "#EE7733",
              ubiquitin = "#009988", other = "grey85")
FAM_LABELS <- c(synaptic = "Synaptic", translation = "Translation",
                oxphos = "OxPhos/mito", ubiquitin = "Deubiquitination", other = "Other")

nes_fill <- function(lims = c(-3, 3), name = "NES (SCZ vs control)") {
  suppressPackageStartupMessages(library(ggplot2))
  scale_fill_gradient2(low = COL_DOWN, mid = "white", high = COL_UP, midpoint = 0,
                       limits = lims, oob = scales::squish, name = name)
}

theme_fig5 <- function(base_size = FS) {
  suppressPackageStartupMessages(library(cowplot))
  theme_cowplot(font_size = base_size) +
    theme(axis.title = element_text(size = AXIS_TITLE),
          axis.text  = element_text(size = AXIS_TEXT),
          axis.line  = element_line(linewidth = LW, colour = "grey20"),
          axis.ticks = element_line(linewidth = LW, colour = "grey20"),
          legend.title = element_text(size = AXIS_TEXT),
          legend.text  = element_text(size = LEGEND_TEXT),
          legend.key.size = unit(9 * FIG_SCALE, "pt"),
          plot.margin = margin(3, 4, 3, 4))
}

# House significance coding, identical to Figures 2 and 3 (see their legends in
# the manuscript). Figure 5 previously used a shifted scale in which ** meant
# FDR < 0.01, so the same glyph carried different meanings in different figures.
stars_of <- function(padj) ifelse(padj < 0.01, "***", ifelse(padj < 0.05, "**",
                          ifelse(padj < 0.10, "*", "")))
# Figure 3's extended form minus its nominal-P mark. The trend tier is worth
# having for gene sets; the nominal-P dot is not, because over the ~6,000-set
# universe here P < 0.05 is expected ~313 times by chance (503 observed), so it
# would mark noise -- and it landed only on positive-side synaptic sets, against
# what the text says. Figure 3 keeps the dot legitimately: 109 crumblr tests.
stars_ext <- function(padj) ifelse(padj < 0.01, "***", ifelse(padj < 0.05, "**",
                       ifelse(padj < 0.10, "*", ifelse(padj < 0.20, "+", ""))))
