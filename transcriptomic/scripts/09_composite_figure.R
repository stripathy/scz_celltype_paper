#!/usr/bin/env Rscript
# ============================================================================
# 09_composite_figure.R — Publication composite (7.1 x 6.625 in, 400 dpi)
# ============================================================================
# Cross-platform SCZ DE story for the two canonical interneuron markers, one
# marker PER ROW (volcano -> forest -> per-donor boxplot -> exemplar cells):
#   Row 1  SST in Sst cells:    a volcano (Sst) | b forest (SST/Sst) |
#                               c CP1K boxplot  | d exemplar cells (Control/SCZ)
#   Row 2  PVALB in Pvalb cells: e volcano (Pvalb) | f forest (PVALB/Pvalb) |
#                               g CP1K boxplot  | h exemplar cells (Control/SCZ)
#   Row 3  i butterfly (up/down DE-gene counts per cell type) |
#          j concordance scatter (snRNA-seq meta logFC vs Xenium logFC)
# Volcanoes: tight, data-driven axes. Forests: 7 snRNA cohorts + pooled meta +
# Xenium. Boxplots: per-donor CP1K (counts/1,000 tx) + edgeR p. Exemplars: the
# representative cell at the pooled group-median grain density (scripts/10).
#
# INPUTS — every one a committed snapshot under data/figure_inputs/ (MANIFEST.tsv
# records the source and its checksum; scripts/00_refresh_figure_inputs.R rebuilds
# the set, and fig_input() below refuses to draw from a snapshot whose source has
# moved on). Sources are the snRNA-seq pipeline's own outputs in ../snrnaseq/ and
# the Xenium repo's DE/pseudobulk tables:
#   DE_genes_all_cells_scz.csv                meta-analytic snRNA-seq DE (a, e, i, j)
#   meta_results_cohorts_subclass_forest.csv  per-dataset DE, SST+PVALB rows (b, f)
#   de_results_subclass.csv                   Xenium DE (b, f, j)
#   plotdata.csv, xen_Sst_proportions.csv     per-dataset / Xenium donor n (b, f labels)
#   snrnaseq_subclass_mean_prop.csv           share of nuclei per subclass (i inset)
#   results/tables/marker_norm_expr.csv       boxplot input (c, g) — scripts/12
#   results/tables/exemplar_*.csv             exemplar cells (d, h) — scripts/10
#
# Output: ../manuscript/figures/main/Fig2_cross_platform_de.{png,pdf}, in place.
# ============================================================================
# Run from anywhere: resolve the component root from this script's own location
# (same convention as genetics/scripts/figures/*.R). Falls back to the cwd when
# sourced interactively, where the user is expected to already be in transcriptomic/.
.script_dir <- function() {
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else NA_character_
}
.sd <- .script_dir()
if (!is.na(.sd)) setwd(normalizePath(file.path(.sd, "..")))
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(cowplot); library(ggrepel); library(metafor); library(ggsignif)
})

# Self-contained figure inputs, committed under data/figure_inputs/ (a subdir, so
# it escapes the single-level data/*.csv ignore — see REPRODUCE.md). The cohorts
# table is the SST+PVALB-only subset the forests need (the full 313 MB per-cohort
# table stays external); the meta table is the full per-(cell type × gene) table.
source("scripts/_figure_inputs.R")     # committed snapshots + staleness guard
INPUT_META    <- fig_input("DE_genes_all_cells_scz.csv")                   # meta-analytic snRNA-seq DE (full)
INPUT_COH     <- fig_input("meta_results_cohorts_subclass_forest.csv")     # per-cohort DE, SST+PVALB rows only
INPUT_XENIUM  <- fig_input("de_results_subclass.csv")                      # Xenium spatial DE

# These three used to be read from the analysis working directory on the cluster,
# which was the one reason Figure 2 did not render from a clean clone. They are
# now snapshots under data/figure_inputs/ like everything else (sources: the two
# tables Nicole committed under snrnaseq/Final_figures/Data/ on 9 Sep, and the
# per-donor nucleus counts behind Fig. 3a). The names are kept so the read sites
# below are unchanged; fig_input() stops with a clear message if one is missing.
EXTERNAL <- c(
  sample_n_cohorts = fig_input("plotdata.csv"),                    # per-dataset donor n (forest labels)
  sample_n_xenium  = fig_input("xen_Sst_proportions.csv"),         # Xenium Sst_25 donor n
  subclass_prop    = fig_input("snrnaseq_subclass_mean_prop.csv")  # mean share of nuclei per subclass (panel i inset)
)

FIG_W <- 7.1    # max total width (inches)

# ---- palette / groupings ---------------------------------------------------
UP_DARK   <- "#D55E00"; UP_LIGHT   <- "#F2B58C"
DOWN_DARK <- "#0072B2"; DOWN_LIGHT <- "#9FCAE6"
COL_NS    <- "grey75"
COL_POOL  <- "black";   COL_XEN <- "#1b9e77"; COL_COHORT <- "grey35"
EXC <- c("L2_3 IT","L4 IT","L5 IT","L5 ET","L5_6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")
INH <- c("Lamp5","Lamp5_Lhx6","Pax6","Pvalb","Sncg","Sst","Sst Chodl","Vip","Chandelier")
GLI <- c("Astro","Oligo","OPC","Micro-PVM","Endo","VLMC")
class_of <- function(ct) dplyr::case_when(
  ct %in% EXC ~ "Excitatory", ct %in% INH ~ "Inhibitory",
  ct %in% GLI ~ "Non-neuronal",       TRUE        ~ NA_character_)
CLASS_COL <- c(Excitatory = "#117733", Inhibitory = "#882255", `Non-neuronal` = "#DDCC77")  # Endo and VLMC are in this class, so not "Glia"

BASE <- 7   # base font size (Nature/NN: all figure text 5–7 pt). Most text is
            # sized relative to this; geom_text multipliers tuned so nothing
            # exceeds 7 pt or drops below 5 pt at the 7.1 × 6.625 in print size.
# Shared axis-text sizes — applied uniformly across panels for consistency
# (hierarchy: subtitle BASE=7 > axis title 6.5 > tick labels 6). The inset keeps
# its own smaller sizes by design.
AXIS_TITLE <- BASE - 0.5   # 6.5 pt — axis titles (the size used in the forest, reads well)
AXIS_TEXT  <- BASE - 1     # 6.0 pt — axis tick labels

# significance helpers (shared)
ast    <- function(fdr) dplyr::case_when(
  is.na(fdr) ~ "", fdr < 0.01 ~ "***", fdr < 0.05 ~ "**", fdr < 0.10 ~ "*", TRUE ~ "")
is_dot <- function(p, fdr) !is.na(p) & p < 0.05 & (is.na(fdr) | fdr >= 0.10)

# ============================================================================
# Load data
# ============================================================================
cohort_map <- c("MSSM"="MSSM 2","HBCC"="HBCC","OFC"="Fröhlich","Bat"="Batiuk",
                "Mclean"="McLean","MtSinai"="MSSM 1","Multi"="Multiome","meta"="Meta-analysis", "MSSM1" = "MSSM 1", "MSSM2" ="MSSM 2", "Frohlich" = "Fröhlich")

cat("Loading meta-analytic DE...\n")
meta_tbl <- read_csv(INPUT_META, show_col_types = FALSE)

FOREST <- tibble::tribble(
  ~gene,   ~cell,   ~lab,
  "SST",   "Sst",   "SST / Sst",
  "PVALB", "Pvalb", "PVALB / Pvalb"
)

cat("Loading per-cohort table (large) and filtering to forest genes...\n")
coh <- read_csv(INPUT_COH, show_col_types = FALSE) |>
  filter(genes %in% FOREST$gene, !is.na(logFC), !is.na(t), t != 0) |>
  mutate(cohort = recode(cohort, !!!cohort_map), SE = logFC / t)

cat("Loading + harmonizing Xenium DE...\n")
ct_map <- c("Astrocyte"="Astro","L2/3 IT"="L2_3 IT","L5/6 NP"="L5_6 NP",
            "Microglia-PVM"="Micro-PVM","Oligodendrocyte"="Oligo","Endothelial"="Endo")

xen_all <- read_csv(INPUT_XENIUM, show_col_types = FALSE) |>
  mutate(cell_type = ifelse(celltype %in% names(ct_map), ct_map[celltype], celltype))

xen_forest <- xen_all |>
  filter(gene %in% FOREST$gene, !is.na(logFC), !is.na(F), F > 0) |>
  mutate(t = sign(logFC) * sqrt(F), SE = logFC / t) |>
  transmute(genes = gene, logFC, t, P.Value = PValue, adj.P.Val = FDR,
            cell_type, cohort = "Xenium", SE)

coh <- bind_rows(
  coh |> select(genes, logFC, t, P.Value, adj.P.Val, cell_type, cohort, SE),
  xen_forest
)

sample_n <- read.csv(EXTERNAL[["sample_n_cohorts"]]) |>
  filter(CellType == "Sst_25", Cohort != "Meta-analysis") |>
  transmute(cohort = recode(Cohort, "MtSinai"="MSSM 1","MSSM"="MSSM 2",
                            "OFC"="Fröhlich","Bat"="Batiuk","Mclean"="McLean",
                            "Multi"="Multiome", "MSSM1" = "MSSM 1", "MSSM2" ="MSSM 2", "Frohlich" = "Fröhlich"), n) |>
  distinct(cohort, n)

xen_n <- read.csv(EXTERNAL[["sample_n_xenium"]]) |>
  filter(subtype == "Sst_25") |> nrow()

sample_n <- bind_rows(sample_n, tibble(cohort="Xenium", n=xen_n))
coh <- coh |> left_join(sample_n, by="cohort")

# ============================================================================
# Panel A — butterfly (ggplot)
# ============================================================================
build_butterfly <- function() {
  d <- meta_tbl |> mutate(direction = ifelse(estimate > 0, "up", "down"))
  cnt <- function(thr) d |> filter(padj < thr) |>
    count(cell_type, direction) |>
    pivot_wider(names_from = direction, values_from = n, values_fill = 0)
  c10 <- cnt(0.10); c05 <- cnt(0.05)
  for (cc in c("up","down")) { if (!cc %in% names(c10)) c10[[cc]] <- 0
                               if (!cc %in% names(c05)) c05[[cc]] <- 0 }
  m <- c10 |> transmute(cell_type, up10 = up, down10 = down) |>
    left_join(c05 |> transmute(cell_type, up05 = up, down05 = down), by = "cell_type") |>
    mutate(across(c(up05, down05), ~tidyr::replace_na(.x, 0)),
           total = up10 + down10) |>
    arrange(total)
  m$cell_type <- factor(m$cell_type, levels = m$cell_type)

  # fill mapped to a 4-level key so a legend is generated; light (FDR<0.10)
  # drawn first, dark (FDR<0.05) overlaid
  ggplot(m) +
    geom_col(aes(y = cell_type, x =  up10,   fill = "Up, FDR < 0.10"),   width = 0.78) +
    geom_col(aes(y = cell_type, x = -down10, fill = "Down, FDR < 0.10"), width = 0.78) +
    geom_col(aes(y = cell_type, x =  up05,   fill = "Up, FDR < 0.05"),   width = 0.78) +
    geom_col(aes(y = cell_type, x = -down05, fill = "Down, FDR < 0.05"), width = 0.78) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "black") +
    geom_text(aes(y = cell_type, x =  up10,   label = up10),
              hjust = -0.15, size = BASE*0.26) +
    geom_text(aes(y = cell_type, x = -down10, label = down10),
              hjust = 1.15, size = BASE*0.26) +
    scale_fill_manual(
      values = c("Up, FDR < 0.10" = UP_LIGHT, "Down, FDR < 0.10" = DOWN_LIGHT,
                 "Up, FDR < 0.05" = UP_DARK,  "Down, FDR < 0.05" = DOWN_DARK),
      breaks = c("Up, FDR < 0.05", "Up, FDR < 0.10",
                 "Down, FDR < 0.05", "Down, FDR < 0.10"), name = NULL) +
    scale_x_continuous(labels = function(v) abs(v),
                       expand = expansion(mult = 0.12)) +
    labs(x = "Number of DE genes", y = NULL) +
    theme_cowplot(font_size = BASE) +
    theme(axis.text.y = element_text(size = AXIS_TEXT),
          axis.text.x = element_text(size = AXIS_TEXT),
          axis.title.x = element_text(size = AXIS_TITLE),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank(),
          legend.position = c(0.23, 0.16),
          legend.justification = c(0.5, 0.5),
          legend.text  = element_text(size = BASE - 1.5),
          legend.key.size = unit(9, "pt"),
          legend.spacing.y = unit(0, "pt"),
          plot.margin = margin(2, 4, 2, 6))
}

# ----------------------------------------------------------------------------
# Panel I inset — # meta DE genes (FDR<0.10) vs cell-type proportion (subclass).
# Proportion = mean per-donor share of nuclei per subclass across the seven
# snRNA-seq datasets (469 donors) -- the same nuclei the DE was run on, so the
# panel reads as what it is: DE-gene count scaling with the number of nuclei
# tested, largely a power effect. Until 2026-09-14 the axis used Xenium cell
# proportions as a proxy; those disagree with the nuclei by up to 25-fold for
# subclasses the 300-gene panel resolves poorly (Lamp5_Lhx6, L5 ET), which put
# both far off the trend for reasons that had nothing to do with DE power.
# Minimal styling:
# no point labels, two log10 % stops, "DE genes (#)" / "Cell proportion (%)".
# Standalone version with cell-type labels: scripts/16_de_vs_proportion.R.
# ----------------------------------------------------------------------------

build_de_prop_inset <- function() {
  nde <- meta_tbl |> group_by(cell_type) |> summarise(n_de = sum(padj < 0.10, na.rm = TRUE), .groups = "drop")
  prop <- read_csv(EXTERNAL[["subclass_prop"]], show_col_types = FALSE) |> transmute(cell_type = CellType, prop = mean_prop)
  d <- inner_join(nde, prop, by = "cell_type") |> filter(n_de >= 1) |> mutate(class = class_of(cell_type))
  rs <- cor(d$prop, d$n_de, method = "spearman")
  lab <- d |> filter(cell_type %in% c("Astro", "L5 IT", "Vip", "L6b")) |> mutate(lbl = gsub("_", "/", cell_type))

  ggplot(d, aes(prop, n_de)) +
    geom_smooth(method = "lm", se = FALSE, colour = scales::alpha("grey25", 0.3), linewidth = 0.5, formula = y ~ x) +
    geom_point(aes(colour = class), size = 0.7, alpha = 0.9) +
    geom_text_repel(data = lab, aes(label = lbl), size = BASE * 0.28,
                    min.segment.length = 0, segment.size = 0.2, segment.colour = "grey55",
                    box.padding = 0.28, point.padding = 0.2, force = 2,
                    max.overlaps = Inf, seed = 3, colour = "grey15") +
    annotate("text", x = max(d$prop), y = 0, label = sprintf("rho == %.2f", rs),
             parse = TRUE, hjust = 1, vjust = 0, size = BASE * 0.26, colour = "grey25") +
    scale_colour_manual(values = CLASS_COL, guide = "none") +
    scale_x_log10(breaks = c(0.01, 0.1), labels = c("1%", "10%")) +
    scale_y_continuous(breaks = c(0, 400, 800), expand = expansion(mult = c(0.05, 0.16))) +
    labs(x = "Cell proportion (%)", y = "DE genes (#)") +
    theme_cowplot(font_size = BASE - 1) +
    theme(legend.position = "none",
          axis.title.x = element_text(size = BASE - 1.5, margin = margin(t = 1)),
          axis.title.y = element_text(size = BASE - 1.5, margin = margin(r = 1)),
          axis.text = element_text(size = BASE - 2.5),
          axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25),
          plot.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.3),
          plot.margin = margin(2, 3, 1, 1))
}

# ============================================================================
# Panels B-E — volcano for one cell type (Sst, L2/3 IT, Astro, Micro-PVM)
# ============================================================================
build_volcano <- function(cell, highlight) {
  # colour by the same FDR tiers as the butterfly (panel A): up/down x
  # FDR<0.05 (dark) / 0.05-0.10 (light); NS = grey
  d <- meta_tbl |> filter(cell_type == cell) |>
    mutate(nlp = -log10(pval),          # y = raw significance; colour still encodes FDR (below)
           tier = case_when(
             padj < 0.05 & estimate > 0 ~ "Up, FDR < 0.05",
             padj < 0.10 & estimate > 0 ~ "Up, FDR < 0.10",
             padj < 0.05 & estimate < 0 ~ "Down, FDR < 0.05",
             padj < 0.10 & estimate < 0 ~ "Down, FDR < 0.10",
             TRUE                       ~ "NS") |>
             factor(levels = c("Down, FDR < 0.05","Down, FDR < 0.10","NS",
                               "Up, FDR < 0.10","Up, FDR < 0.05")))
  lab <- d |> filter(genes %in% highlight) |>
    mutate(face = ifelse(genes == toupper(cell), "bold.italic", "italic"))  # bold the cell-type marker (SST / PVALB)
  # Tight, data-driven limits PER PANEL (not symmetric, not shared across the four
  # volcanoes); gene labels are kept inside the panel via ggrepel xlim/ylim +
  # coord_cartesian, so the axes crop close to the data without clipping a label.
  xr  <- range(d$estimate, na.rm = TRUE); xpad <- diff(xr) * 0.07
  xlo <- xr[1] - xpad; xhi <- xr[2] + xpad
  yhi <- max(d$nlp, na.rm = TRUE) * 1.08
  # BH-FDR is monotone in p, so {FDR<0.10} == {p <= p_thr}; the dashed line at this
  # p-boundary keeps "above the line = coloured (FDR<0.10), below = grey" exact.
  p_thr <- if (any(d$padj < 0.10, na.rm = TRUE)) max(d$pval[d$padj < 0.10], na.rm = TRUE) else NA_real_
  tier_cols <- c("Up, FDR < 0.05" = UP_DARK,  "Up, FDR < 0.10" = UP_LIGHT,
                 "Down, FDR < 0.05" = DOWN_DARK, "Down, FDR < 0.10" = DOWN_LIGHT,
                 "NS" = COL_NS)
  ggplot(d, aes(estimate, nlp)) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.25) +
    {if (!is.na(p_thr)) geom_hline(yintercept = -log10(p_thr), colour = "grey70",
               linetype = "dashed", linewidth = 0.25)} +
    geom_point(aes(colour = tier), data = ~filter(.x, tier == "NS"),
               size = 0.35, alpha = 0.3) +
    geom_point(aes(colour = tier), data = ~filter(.x, tier != "NS"),
               size = 0.5, alpha = 0.85) +
    geom_point(data = lab, shape = 21, fill = NA, colour = "black",
               size = 1.4, stroke = 0.4) +
    geom_text_repel(data = lab, aes(label = genes), size = BASE*0.32,
                    fontface = lab$face, min.segment.length = 0,
                    segment.size = 0.25, segment.colour = "grey55",
                    box.padding = 0.4, point.padding = 0.3, force = 4,
                    max.overlaps = Inf, seed = 2, colour = "grey10",
                    xlim = c(xlo, xhi), ylim = c(0, yhi)) +
    scale_colour_manual(values = tier_cols, guide = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    coord_cartesian(xlim = c(xlo, xhi), ylim = c(0, yhi), clip = "on") +
    labs(x = expression("SCZ log"[2]~"FC"),
         y = expression(-log[10]~italic(P)),
         subtitle = paste0(gsub("_", "/", cell), " cells")) +
    theme_cowplot(font_size = BASE) +
    theme(plot.subtitle = element_text(size = BASE, face = "italic", hjust = 0.5,
                                       margin = margin(b = 1)),
          axis.text  = element_text(size = AXIS_TEXT),
          axis.title = element_text(size = AXIS_TITLE),
          plot.margin = margin(2, 3, 2, 2))
}

# ============================================================================
# Panels F-H — compact forest
# ============================================================================
build_forest <- function(gene_sym,cell,ttl,show_xlab=FALSE){
  sub <- coh |> filter(genes==gene_sym,cell_type==cell); snrna <- sub |> filter(cohort!="Xenium"); xen_in <- sub |> filter(cohort=="Xenium")
  if(nrow(snrna)==0) return(NULL)
  re <- tryCatch(rma(yi=snrna$logFC,sei=snrna$SE,method="DL",control=list(maxiter=500)),error=function(e) NULL)

  cohort_df <- snrna |>
    transmute(cohort,est=logFC,lo=logFC-1.96*SE,hi=logFC+1.96*SE,n,p=P.Value,padj=adj.P.Val,
              sig=ast(adj.P.Val),dot=is_dot(P.Value,adj.P.Val)) |>
    arrange(desc(est)) |>
    mutate(y=rev(seq_along(est))+3,role="cohort",colour=COL_COHORT)

  meta_fdr <- {r <- meta_tbl |> filter(genes==gene_sym,cell_type==cell); if(nrow(r)==1) r$padj else NA_real_}
  pool <- if(!is.null(re)) tibble(cohort="Meta-analysis",est=as.numeric(re$beta),lo=re$ci.lb,hi=re$ci.ub,
                                  n=sum(snrna$n,na.rm=TRUE),p=re$pval,padj=meta_fdr,sig=ast(meta_fdr),
                                  dot=FALSE,y=3,role="pool",colour=COL_POOL) else NULL

  has_xen <- nrow(xen_in)>0
  xen_row <- if(has_xen){xr <- xen_in[1,]; tibble(cohort="Xenium",est=xr$logFC,lo=xr$logFC-1.96*xr$SE,
                                                  hi=xr$logFC+1.96*xr$SE,n=xr$n,p=xr$P.Value,padj=xr$adj.P.Val,
                                                  sig=ast(xr$adj.P.Val),dot=is_dot(xr$P.Value,xr$adj.P.Val),
                                                  y=1,role="xenium",colour=COL_XEN)} else NULL

  pd <- bind_rows(cohort_df,pool,xen_row)
  data_lo <- min(pd$lo,na.rm=TRUE); data_hi <- max(pd$hi,na.rm=TRUE); rng <- data_hi-data_lo; off <- rng*.03
  pd <- pd |> mutate(ast_x=ifelse(est>=0,hi+off,lo-off),ast_h=ifelse(est>=0,0,1),
                     has_m=sig!=""|dot|(role=="pool"&sig==""))
  rt <- any(pd$has_m&pd$ast_h==0); lf <- any(pd$has_m&pd$ast_h==1)
  xlo <- data_lo-rng*(if(lf).22 else .05); xhi <- data_hi+rng*(if(rt).22 else .05)
  xlo <- min(xlo,-off); xhi <- max(xhi,off)
  ylab <- pd |> select(y,cohort) |> arrange(y)
  if(has_xen) ylab <- bind_rows(ylab,tibble(y=2,cohort=""))

  ggplot(pd,aes(est,y,colour=I(colour))) +
    geom_vline(xintercept=0,linetype="dashed",colour="grey65",linewidth=.25) +
    {if(has_xen) geom_hline(yintercept=2,colour="grey85",linetype="dotted",linewidth=.3)} +
    geom_segment(aes(x=lo,xend=hi,yend=y),linewidth=.45) +
    geom_point(data=cohort_df,aes(est,y,size=n),shape=21,fill=COL_COHORT,colour=COL_COHORT,stroke=.3) +
    {if(!is.null(pool)) geom_point(data=pool,aes(est,y,size=n),shape=23,fill=COL_POOL,colour=COL_POOL)} +
    {if(has_xen) geom_point(data=xen_row,aes(est,y,size=n),shape=17,colour=COL_XEN)} +
    geom_text(data=pd,aes(label=sig,x=ast_x,hjust=ast_h),size=BASE*.34,vjust=.75,
              colour="grey15",fontface="bold") +
    geom_point(data=filter(pd,dot),aes(ast_x,y),shape=16,size=.8,colour="grey15",inherit.aes=FALSE) +
    geom_text(data=filter(pd,role=="pool"&sig==""),aes(ast_x,y,hjust=ast_h),label="n.s.",
              size=BASE*.28,vjust=.5,colour="grey45",fontface="italic",inherit.aes=FALSE) +
    scale_size_continuous(range=c(1,3.3),breaks=c(50,200,400),name="n") +
    scale_y_continuous(breaks=ylab$y,labels=ylab$cohort,expand=expansion(add=.6)) +
    scale_x_continuous(limits=c(xlo,xhi),breaks=scales::pretty_breaks(3),expand=expansion(mult=.01)) +
    labs(x=if(show_xlab) expression("SCZ log"[2]~"FC (95% CI)") else NULL,y=NULL,subtitle=ttl) +
    theme_cowplot(font_size=BASE) +
    theme(plot.subtitle=element_text(size=BASE,face="plain",margin=margin(b=1)),
          axis.text.y=element_text(size=AXIS_TEXT),axis.text.x=element_text(size=AXIS_TEXT),
          axis.title.x=element_text(size=AXIS_TITLE,margin=margin(t=1.5)),
          axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),
          legend.position="none",plot.margin=margin(2,4,2,2))
}

# ============================================================================
# Panel I — concordance scatter
# ============================================================================
build_scatter <- function() {
  meta <- meta_tbl |> filter(padj < 0.10) |>
    transmute(genes, cell_type, meta_est = estimate, meta_padj = padj)
  xen <- xen_all |> transmute(genes = gene, cell_type, xen_logFC = logFC, xen_p = PValue)
  pr <- inner_join(meta, xen, by = c("genes","cell_type")) |>
    mutate(class = class_of(cell_type),
           concordant = sign(meta_est) == sign(xen_logFC),
           fdr_bin = factor(ifelse(meta_padj < 0.05, "< 0.05", "0.05-0.10"),
                            levels = c("< 0.05","0.05-0.10")))
ct <- cor.test(pr$meta_est, pr$xen_logFC, method="spearman", exact=FALSE)
rho <- unname(ct$estimate)
pval <- ct$p.value
p_lab <- if (pval < 0.001) "p<0.001" else sprintf("p==%.3f",pval)

pc <- round(100*mean(pr$concordant))
  lab_pairs <- tibble::tribble(~genes,~cell_type,
    "SST","Sst","BDNF","L2_3 IT","FKBP5","OPC","CX3CR1","Micro-PVM",
    "SMAD1","Pvalb","SERPING1","Astro","FGFR3","Astro",
    "VGF","Chandelier","CALB1","L6 IT","ATP2B4","Sst")
  lab_keys <- paste(lab_pairs$genes, lab_pairs$cell_type)
  lab <- pr |> inner_join(lab_pairs, by = c("genes","cell_type"))
  lim <- max(as.numeric(quantile(abs(c(pr$meta_est, pr$xen_logFC)), 0.97)),
             max(abs(c(lab$meta_est, lab$xen_logFC)))) * 1.08
  # Labels via ggrepel auto-placement (robust to panel resizing); seed fixed for
  # determinism. Constrained to the panel; r / % concordant sit in the sparse
  # bottom-right corner.
  pr <- pr |> mutate(
    tag = sprintf("%s (%s)", genes, gsub("_", "/", cell_type)),
    rep_label = ifelse(paste(genes, cell_type) %in% lab_keys, tag, ""),
    # ATP2B4/Sst and CALB1/L6 IT are both up-regulated and land close in the
    # up-right cluster; bias ATP2B4's label up-left so the two don't collide.
    nudge_x = ifelse(genes == "ATP2B4" & cell_type == "Sst", -0.16, 0),
    nudge_y = ifelse(genes == "ATP2B4" & cell_type == "Sst",  0.18, 0))

  ggplot(pr, aes(meta_est, xen_logFC)) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.25) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.25) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey80", linewidth = 0.3) +
    geom_point(aes(colour = class, size = fdr_bin), shape = 16, alpha = 0.8) +
    geom_point(data = lab, shape = 21, fill = NA, colour = "black",
               size = 1.8, stroke = 0.5) +
    geom_text_repel(data = pr, aes(label = rep_label), size = BASE*0.28,
                    min.segment.length = 0, segment.size = 0.25,
                    segment.colour = "grey55", box.padding = 0.5,
                    point.padding = 0.3, force = 5, force_pull = 0.12,
                    nudge_x = pr$nudge_x, nudge_y = pr$nudge_y,
                    xlim = c(-lim, lim), ylim = c(-lim, lim),
                    max.overlaps = Inf, seed = 7, colour = "grey10") +
annotate("text", x = lim*0.97, y = -lim*0.76,
         label = paste0("rho=='", sprintf("%.2f", rho), "'*','~~", p_lab),
         parse = TRUE, hjust = 1, size = BASE*0.34) +
annotate("text", x = lim*0.97, y = -lim*0.93,
         label = sprintf("'%d%% concordant'", pc), parse = TRUE,
         hjust = 1, size = BASE*0.30) +
    scale_colour_manual(values = CLASS_COL, name = NULL) +
    scale_size_manual(values = c("< 0.05" = 1.8, "0.05-0.10" = 0.7),
                      name = "meta FDR") +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(x = expression("snRNA-seq meta-analysis  log"[2]~"FC"),
         y = expression("Xenium  log"[2]~"FC")) +
    theme_cowplot(font_size = BASE) +

    guides(
      colour = guide_legend(title=NULL),
      size = "none"
    ) +

    theme(
      legend.position = c(0.98,0.25),
      legend.justification = c(1,0),
      legend.background = element_rect(fill="white",colour=NA),
      axis.text = element_text(size=AXIS_TEXT),
      axis.title = element_text(size=AXIS_TITLE),
      aspect.ratio = 1,
      plot.margin = margin(2,3,2,2)
    )
}
# ============================================================================
# Panels d, h — Xenium exemplar cells (boundary, nucleus, marker transcripts).
# Requires scripts/10_xenium_exemplar_cells.py to have been run. The builder
# lives in scripts/_exemplar_panels.R so a preview draws exactly what the figure
# draws; see that file for what each label is for and why.
# ============================================================================
TAB <- "results/tables"
source("scripts/_exemplar_panels.R")
# How the two outlines are named for the reader (G#123 review item):
#   "oncell"  "cell" / "nucleus" with short leaders on the first cell (d, Control)
#   "key"     cells left clean; a one-line key drawn under panel d instead
EX_STYLE <- "oncell"

# ============================================================================
# Panel J — library-normalised expression (CP1K = counts per 1,000 transcripts)
# with the edgeR DE p-value. Same Xenium edgeR DE as the forests (F–H) and
# scatter (I): per-donor CP1K (TMM) from scripts/12_marker_norm_expr.R, p-value =
# edgeR glmQLFTest (exactly the value driving the forest/scatter significance).
# ============================================================================
NEXPR <- read_csv(sprintf("%s/marker_norm_expr.csv", TAB), show_col_types = FALSE)
NEXPR$dx <- factor(NEXPR$dx, levels = c("Control", "SCZ"))
NSTAT <- read_csv(sprintf("%s/marker_norm_expr_stats.csv", TAB), show_col_types = FALSE)
NEXPR_COL <- c(Control = DOWN_DARK, SCZ = UP_DARK)

build_normexpr <- function(gene, title = NULL, show_x = FALSE) {
  df <- NEXPR[NEXPR$gene == gene, ]
  p  <- NSTAT$p[NSTAT$gene == gene]
  yr <- range(df$cp1k)
  ggplot(df, aes(dx, cp1k)) +
    geom_boxplot(aes(fill = dx), width = 0.62, outlier.shape = NA, alpha = 0.55, linewidth = 0.4) +
    geom_jitter(aes(fill = dx), shape = 21, colour = "grey25", size = 1.4, stroke = 0.3,
                width = 0.13, height = 0, alpha = 0.9) +
    geom_signif(comparisons = list(c("Control", "SCZ")),
                annotations = sprintf("p=='%.3f'", p), parse = TRUE,
                y_position = yr[2] + 0.06 * diff(yr), tip_length = 0.02,
                textsize = BASE * 0.32, vjust = -0.1) +
    scale_fill_manual(values = NEXPR_COL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22)), n.breaks = 4) +
    labs(x = NULL, y = sprintf("%s expr (CP1K)", gene), subtitle = title) +
    theme_cowplot(font_size = BASE) +
    theme(legend.position = "none",
          plot.subtitle = element_text(size = BASE - 1, face = "plain", hjust = 0.5,
                                       margin = margin(b = 1)),
          axis.text.x  = if (show_x) element_text(size = AXIS_TEXT) else element_blank(),
          axis.ticks.x = if (show_x) element_line() else element_blank(),
          axis.line.x  = if (show_x) element_line() else element_blank(),
          axis.title.y = element_text(size = AXIS_TITLE),
          axis.text.y  = element_text(size = AXIS_TEXT),
          plot.margin = margin(2, 3, 2, 3))
}

# ============================================================================
# Assemble
# ============================================================================
cat("Building panels...\n")
# --- per-marker panels (one marker per row) ---
vSst   <- build_volcano("Sst",   c("SST","NAT16","SMAD1","AFG3L2","STAC","KCTD4","SLC9A9","DRD3"))
vPvalb <- build_volcano("Pvalb", c("PVALB","SMAD1","ANXA2","SCN3A","NAT16","VGF","CIRBP","TCAF2","FGF10"))
fSst   <- build_forest("SST",   "Sst",   "SST / Sst",     show_xlab = TRUE)
fPvalb <- build_forest("PVALB", "Pvalb", "PVALB / Pvalb", show_xlab = TRUE)
bSst   <- build_normexpr("SST",   show_x = TRUE)
bPvalb <- build_normexpr("PVALB", show_x = TRUE)
# Butterfly (i) carries the DE-genes-vs-proportion inset in its empty bottom-left.
pA <- ggdraw(build_butterfly()) +
  draw_plot(build_de_prop_inset(), x=0.75, y=0.1, width=0.36, height=0.378)
pJ     <- build_scatter()

# Exemplar cells: one shared coordinate limit -> identical zoom and an identical
# 5 um bar in all four cells. The bar is labelled once, on d/Control, where a
# reader meets the panel first. The outline labels need extra room around the
# largest cell; the key does not.
EX_PAIRS <- list(c("SST","Control"), c("SST","SCZ"), c("PVALB","Control"), c("PVALB","SCZ"))
ex_lim   <- exemplar_lim(EX_PAIRS, TAB, expand = if (EX_STYLE == "oncell") 1.22 else 1.08)
eSst   <- ex_pair("SST",   ex_lim, TAB, BASE, show_hdr = TRUE,  scalebar_lab_on = "Control",
                  label_outlines_on = if (EX_STYLE == "oncell") "Control" else "none",
                  key = if (EX_STYLE == "key") "draw" else "none")
ePvalb <- ex_pair("PVALB", ex_lim, TAB, BASE, show_hdr = FALSE, scalebar_lab_on = "none",
                  key = if (EX_STYLE == "key") "spacer" else "none")

# Row 1 (SST) / Row 2 (PVALB): [volcano | forest | CP1K boxplot] are aligned with
# cowplot align="h"/axis="tb" so the boxplot (c,g) x-axis lines up with the volcano
# and forest (a,b); the exemplar pair is appended as the 4th column.
RW3 <- c(1.05, 0.95, 0.70)
mkrow <- function(volc, forest, box, expair, labs)
  plot_grid(
    plot_grid(volc, forest, box, ncol = 3, rel_widths = RW3, align = "h", axis = "tb",
              labels = labs[1:3], label_size = 8, label_fontface = "bold"),
    expair, ncol = 2, rel_widths = c(sum(RW3), 1.25),
    labels = c("", labs[4]), label_size = 8, label_fontface = "bold")
row1 <- mkrow(vSst,   fSst,   bSst,   eSst,   c("a","b","c","d"))
row2 <- mkrow(vPvalb, fPvalb, bPvalb, ePvalb, c("e","f","g","h"))
# Row 3: butterfly (i) | concordance scatter (j) — equal (50/50) width
row3 <- plot_grid(pA, pJ, ncol = 2, rel_widths = c(1, 1),
                  labels = c("i","j"), label_size = 8, label_fontface = "bold")

# rel_heights given in INCHES (they sum to FIG_H): rows 1-2 (SST, PVALB) = 2.00 in
# each; row 3 (butterfly + scatter) = 2.625 in.
full <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(2, 2, 2.625))

FIG_H <- 6.625 # 7.1 x 6.625 in. Rows 1-2 (SST a-d, PVALB e-h) = 2.00 in each;
               # butterfly + scatter (i,j) = 2.625 in (rel_heights are inches).
# Written in place into the submission folder, like every other renderer in the
# repo (manuscript/figures/main/README.md): one original, no copy step to drift.
FIG_OUT <- "../manuscript/figures/main/Fig2_cross_platform_de"
ggsave(paste0(FIG_OUT, ".png"), full, width = FIG_W, height = FIG_H, dpi = 400, bg = "white")
ggsave(paste0(FIG_OUT, ".pdf"), full, width = FIG_W, height = FIG_H, bg = "white")
cat(sprintf("Saved %s.{png,pdf}  (%.1f x %.1f in)\n", FIG_OUT, FIG_W, FIG_H))
