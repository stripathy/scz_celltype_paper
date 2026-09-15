#!/usr/bin/env Rscript
# RESERVE (not in the current manuscript, 2026-09-01) | Controls and formal tests
# behind the Sst depletion-group analysis (S8). Dropped with the analysis's move
# to a supplement; kept intact so it can be revived if a reviewer asks. Writes to
# results/, not the submission folder.
#
# Combines what were two separate supplements, since both are cited from the same
# stretch of the results and neither fills a page alone:
#   a  the cell-count confound and its removal by matching
#   b  gene-set burden per group, unmatched versus five cell-matched draws
#   c  diagnosis x depletion-group interaction for the Fig. S8d gene sets, plus
#      the two sets named in the text
#
# a and b support the sentence on downsampling control nuclei; c supports the
# sentence on formal interaction models, which otherwise cites no figure.
#
# Inputs: subsample_cells.py, sensitivity_cellsubsample.R, 05b_interaction_gsea.R.
source("transcriptomic/scripts/sst_strata/_common.R")
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })

FIGSTEM <- "reserve_strata_controls"

SUB   <- file.path(P$supp, "cellsub")
reps  <- sort(as.integer(str_remove(basename(Sys.glob(file.path(SUB, "rep*"))), "rep")))
S_LAB <- c(depleted = "Depleted", intermediate = "Intermediate", non_depleted = "Non-depleted")
lab_group     <- function(x) factor(S_LAB[x], unname(S_LAB))
lab_group_rev <- function(x) factor(S_LAB[x], rev(unname(S_LAB)))   # horizontal forest

# ---- a | the imbalance, before and after cell matching -----------------------
imb <- function(meta, state) {
  map_dfr(STRATA_LEVELS, function(st) {
    d <- meta |> filter(stratum == st) |> mutate(l = log(n_cells))
    co <- summary(lm(l ~ diagnosis + cohort, data = d))$coefficients["diagnosisSCZ", ]
    tibble(stratum = st, state = state, est = co[1], se = co[2], p = co[4])
  })
}
after <- map_dfr(export_cohorts(), ~ read_csv(
  file.path(SUB, "rep1", paste0(.x, "_donor_stratum_meta.csv")),
  col_types = cols(donor = col_character(), .default = col_guess())))
d_imb <- bind_rows(imb(all_donor_meta(), "Observed"), imb(after, "Cell-matched")) |>
  mutate(state = factor(state, c("Observed", "Cell-matched")),
         group = lab_group_rev(stratum),
         pct = 100 * (exp(est) - 1),
         lo = 100 * (exp(est - 1.96 * se) - 1), hi = 100 * (exp(est + 1.96 * se) - 1))
say("panel a, case-control difference in nuclei per donor, percent:")
print(d_imb |> transmute(group, state, pct = round(pct, 1), p = signif(p, 2)))

p_a <- ggplot(d_imb, aes(pct, group, colour = state)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.6 * GEO,
                 position = position_dodge(width = 0.55)) +
  geom_point(size = 2.4 * GEO, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c(Observed = "grey35", `Cell-matched` = "#B2182B"), name = NULL) +
  labs(x = "Nuclei per donor, SCZ vs control (%)", y = NULL) +
  theme_strata() +
  theme(legend.position = c(0.02, 0.14), legend.justification = c(0, 0),
        legend.background = element_blank())

# ---- b | burden, unmatched versus cell-matched -------------------------------
gsea_un <- read_csv(file.path(P$pb, "gsea_all_signatures.csv"), show_col_types = FALSE)
sub <- read_csv(file.path(P$supp, "sensitivity_cellsub_gsea.csv"), show_col_types = FALSE)
b_un <- gsea_un |> filter(signature %in% STRATA_LEVELS) |> group_by(stratum = signature) |>
  summarise(n = sum(padj < 0.10), .groups = "drop") |> mutate(group = lab_group(stratum))
b_mt <- sub |> group_by(stratum, rep) |> summarise(n = sum(padj < 0.10), .groups = "drop") |>
  mutate(group = lab_group(stratum))

p_b <- ggplot(b_un, aes(group, n)) +
  geom_col(fill = "grey82", width = 0.62) +
  geom_point(data = b_mt, colour = "#B2182B", size = 1.9 * GEO, alpha = 0.85,
             position = position_jitter(width = 0.10, height = 0, seed = 1)) +
  labs(x = NULL, y = "Gene sets at FDR < 0.10") +
  theme_strata() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ---- c | the interaction -----------------------------------------------------
ig <- read_csv(file.path(P$pb, "interaction_gsea.csv"), show_col_types = FALSE) |>
  filter(coef == "dxSCZ:stratumdepleted")
CITED <- c(REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE = "translation",
           GOBP_ELECTRON_TRANSPORT_CHAIN = "oxphos")
BL <- c(synaptic = "Synaptic", ubiquitin = "Ubiquitin",
        translation = "Cytosolic translation", oxphos = "Oxidative phosphorylation")
pretty_set <- function(x) {
  y <- x |> str_remove("^GOBP_|^GOCC_|^GOMF_|^REACTOME_") |>
    str_replace_all("_", " ") |> str_to_sentence()
  y <- str_replace(y, "^Srp dependent cotranslational protein targeting to membrane",
                   "SRP-dependent cotranslational protein targeting")
  str_replace(y, "k63", "K63")
}
d_c <- bind_rows(BLOCKS |> select(pathway, block) |> mutate(cited = FALSE),
                 tibble(pathway = names(CITED), block = unname(CITED), cited = TRUE)) |>
  inner_join(ig, by = "pathway") |>
  mutate(block = factor(BL[block], unname(BL)), lab = pretty_set(pathway),
         stars = stars_ext(padj))
say("panel c, interaction NES:")
print(d_c |> transmute(block, lab = str_trunc(lab, 44), NES = round(NES, 2),
                       padj = signif(padj, 2), cited))

p_c <- ggplot(d_c, aes(NES, reorder(lab, NES), colour = block)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_segment(aes(x = 0, xend = NES, yend = reorder(lab, NES)), linewidth = 0.7 * GEO) +
  geom_point(aes(shape = cited), size = 2.4 * GEO) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  geom_text(aes(label = stars), hjust = ifelse(d_c$NES < 0, 1.4, -0.4),
            size = LBL_SMALL, colour = "grey20") +
  scale_colour_manual(values = setNames(unname(FAM_COLS[names(BL)]), unname(BL)), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.14, 0.14))) +
  labs(x = "Interaction NES (depleted vs non-depleted)", y = NULL) +
  theme_strata() +
  theme(legend.position = c(0.01, 0.97), legend.justification = c(0, 1),
        legend.background = element_blank(), legend.key.size = unit(7, "pt"),
        axis.text.y = element_text(size = LBL_SMALL * 2.845))

top <- plot_grid(p_a, p_b, nrow = 1, rel_widths = c(0.58, 0.42), align = "h", axis = "tb",
                 labels = c("a", "b"), label_size = PANEL_LABEL)
fig <- plot_grid(top, p_c, ncol = 1, rel_heights = c(0.40, 0.60),
                 labels = c("", "c"), label_size = PANEL_LABEL)
for (ext in c("png", "pdf"))
  ggsave(file.path(P$supp, paste0(FIGSTEM, ".", ext)), fig,
         width = FIG_W, height = 5.6, dpi = FIG_DPI, bg = "white")
say("wrote %s/%s.(png|pdf); matched draws %s", P$supp, FIGSTEM, paste(reps, collapse = ", "))
