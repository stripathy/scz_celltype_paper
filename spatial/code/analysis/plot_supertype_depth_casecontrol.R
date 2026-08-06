#!/usr/bin/env Rscript
# Control vs SCZ cortical depth per NEURONAL supertype, with a per-supertype linear
# mixed-effects model (depth ~ diagnosis + (1|donor)), BH-FDR across supertypes.
# Single composite figure: (a) Glutamatergic row, (b) GABAergic row. Supertypes
# grouped by subclass and ordered pia->WM by median depth. Mirrors the existing
# supertype_depth_violins layout (depth on y, pia at top, layer bands).
#
# Violins show the INNER 95% of cells per supertype x diagnosis (2.5-97.5 pctile) so
# long outlier tails don't dominate; the mixed model uses ALL cells. Significance
# brackets/stars are drawn ONLY for supertypes with FDR < 0.10 (none here).
# (Brackets drawn manually; ggsignif manual mode is broken under ggplot2 4.0.)
#
# Data: build_supertype_depth_casecontrol_data.py -> output/depth_casecontrol/
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(lmerTest); library(cowplot)
})
cowplot::set_null_device("agg")

DD <- "output/depth_casecontrol"                        # per-cell data + results CSV
FIGDIR <- "../manuscript/figures/supplementary"                  # the single home for submission figures
FIGSTEM <- "S10_supertype_depth_by_diagnosis"
DX_COL <- c(Control = "#4C9BD4", SCZ = "#E2625A"); SIG_COL <- "#B11226"
SIG_FDR <- 0.10
LAYER_B <- c(.10, .40, .55, .70, .90)
LAYER_T <- c(.05, .25, .475, .625, .80, .95); LAYER_N <- c("L1", "L2/3", "L4", "L5", "L6", "WM")
GLUT_ORDER <- c("L2/3 IT","L4 IT","L5 IT","L5 ET","L5/6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")
GABA_ORDER <- c("Lamp5","Lamp5 Lhx6","Sncg","Vip","Pax6","Chandelier","Pvalb","Sst","Sst Chodl")
SUBORDER <- list(Glutamatergic = GLUT_ORDER, GABAergic = GABA_ORDER)
BASE <- 13

df <- read_csv(file.path(DD, "supertype_depth_casecontrol.csv.gz"), show_col_types = FALSE)
df$dx <- factor(df$dx, levels = c("Control", "SCZ"))

cnt <- df |> group_by(supertype, dx) |>
  summarise(n = n(), nd = n_distinct(donor), .groups = "drop") |>
  pivot_wider(names_from = dx, values_from = c(n, nd), values_fill = 0)
keep <- cnt |> filter(n_Control >= 50, n_SCZ >= 50, nd_Control >= 3, nd_SCZ >= 3) |> pull(supertype)
df <- df |> filter(supertype %in% keep)

# ---- per-supertype mixed model on ALL cells: depth ~ dx + (1|donor) ----
fitres <- lapply(keep, function(s) {
  d <- df[df$supertype == s, ]
  m <- tryCatch(lmerTest::lmer(depth ~ dx + (1 | donor), data = d),
                error = function(e) NULL, warning = function(w) NULL)
  if (is.null(m)) return(data.frame(supertype = s, beta = NA, p = NA))
  co <- summary(m)$coefficients
  if (!"dxSCZ" %in% rownames(co)) return(data.frame(supertype = s, beta = NA, p = NA))
  data.frame(supertype = s, beta = co["dxSCZ", "Estimate"], p = co["dxSCZ", "Pr(>|t|)"])
}) |> bind_rows()

res <- df |> group_by(class, subclass, supertype) |>
  summarise(med_depth = median(depth), n_ctrl = sum(dx == "Control"),
            n_scz = sum(dx == "SCZ"), .groups = "drop") |>
  left_join(fitres, by = "supertype")
res$fdr <- p.adjust(res$p, method = "BH")
# Breaks must cover SIG_FDR: brackets are drawn for fdr < SIG_FDR (0.10), so a tier
# for [0.05, 0.10) is required or a marginal hit gets a bracket labelled "ns".
# That tier is labelled with its FDR rather than a symbol -- it needs no legend key,
# and the base pdf() device cannot encode a dagger.
res$sig <- ifelse(is.na(res$fdr), "",
                  as.character(cut(res$fdr, c(-Inf, .001, .01, .05, SIG_FDR, Inf),
                                   labels = c("***", "**", "*", "marginal", "ns"))))
marg <- !is.na(res$fdr) & res$sig == "marginal"
res$sig[marg] <- sprintf("FDR=%.2f", res$fdr[marg])
write_csv(res |> arrange(class, subclass, med_depth),
          file.path(DD, "supertype_depth_casecontrol_lmm.csv"))
cat(sprintf("tested %d supertypes | min FDR = %.3f | FDR<0.10: %d | FDR<0.05: %d | nominal p<0.05: %d\n",
            sum(!is.na(res$p)), min(res$fdr, na.rm = TRUE),
            sum(res$fdr < SIG_FDR, na.rm = TRUE), sum(res$fdr < .05, na.rm = TRUE),
            sum(res$p < .05, na.rm = TRUE)))

shorten <- function(sup, sub) sub(paste0("^", sub, "[_ ]"), "", sup)

plot_class <- function(cls) {
  ro <- res |> filter(class == cls, !is.na(p)) |>
    mutate(subclass = factor(subclass, levels = SUBORDER[[cls]])) |>
    arrange(subclass, med_depth) |>
    mutate(xpos = row_number(), short = mapply(shorten, supertype, as.character(subclass)))
  d <- df |> filter(supertype %in% ro$supertype) |>
    left_join(ro |> select(supertype, xpos), by = "supertype")
  # inner-95% clip per supertype x dx (display only; model used all cells)
  dc <- d |> group_by(supertype, dx) |>
    filter(depth >= quantile(depth, .025), depth <= quantile(depth, .975)) |> ungroup()
  bnds <- ro |> group_by(subclass) |> summarise(mx = max(xpos), mid = mean(xpos), .groups = "drop")
  vlines <- head(sort(bnds$mx), -1) + 0.5
  sig_ann <- ro |> filter(!is.na(fdr) & fdr < SIG_FDR) |>
    transmute(xmin = xpos - 0.2, xmax = xpos + 0.2, xmid = xpos, yb = -0.05, lab = sig)
  lab_y <- -0.105; ylim_top <- -0.15
  g <- ggplot(dc, aes(x = xpos, y = depth, fill = dx, group = interaction(xpos, dx))) +
    geom_violin(position = position_dodge(0.8), width = 0.78, alpha = 0.85,
                linewidth = 0.15, color = "grey30", scale = "width") +
    stat_summary(fun = median, geom = "point", position = position_dodge(0.8),
                 size = 0.7, color = "white") +
    geom_hline(yintercept = LAYER_B, linetype = "dashed", color = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = vlines, color = "grey55", linewidth = 0.3)
  if (nrow(sig_ann) > 0) {                                   # brackets only when FDR<0.10
    g <- g +
      geom_segment(data = sig_ann, aes(x = xmin, xend = xmax, y = yb, yend = yb), inherit.aes = FALSE, linewidth = 0.3, color = SIG_COL) +
      geom_segment(data = sig_ann, aes(x = xmin, xend = xmin, y = yb, yend = yb + 0.018), inherit.aes = FALSE, linewidth = 0.3, color = SIG_COL) +
      geom_segment(data = sig_ann, aes(x = xmax, xend = xmax, y = yb, yend = yb + 0.018), inherit.aes = FALSE, linewidth = 0.3, color = SIG_COL) +
      geom_text(data = sig_ann, aes(x = xmid, y = yb - 0.012, label = lab), inherit.aes = FALSE, size = 3.0, vjust = 0, color = SIG_COL)
  }
  g +
    annotate("text", x = bnds$mid, y = lab_y, label = as.character(bnds$subclass),
             fontface = "bold", size = 3.5, vjust = 0) +
    scale_y_reverse(breaks = LAYER_T, labels = LAYER_N, expand = expansion(mult = c(0.02, 0))) +
    scale_x_continuous(breaks = ro$xpos, labels = ro$short, expand = expansion(add = 0.7)) +
    scale_fill_manual(values = DX_COL, name = NULL) +
    coord_cartesian(ylim = c(1.02, ylim_top)) +
    labs(x = NULL, y = "cortical depth (pia at top, WM at bottom)", title = cls) +
    theme_cowplot(font_size = BASE) +
    theme(axis.text.x = element_text(size = BASE - 4, angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = BASE - 1), axis.line = element_line(linewidth = 0.3),
          plot.title = element_text(size = BASE + 1, face = "bold"),
          legend.position = "none", plot.margin = margin(6, 6, 4, 6))
}

pa <- plot_class("Glutamatergic"); pb <- plot_class("GABAergic")
leg <- get_legend(pa + theme(legend.position = "top", legend.justification = "center"))
grid <- plot_grid(pa, pb, ncol = 1, labels = c("a", "b"), label_size = 17, align = "v", axis = "lr")
# Computed, never hardcoded — the previous literal caption ("103 supertypes; 0/103")
# silently went stale when the data were regenerated.
n_tested <- sum(!is.na(res$p)); n_sig <- sum(res$fdr < SIG_FDR, na.rm = TRUE)
note <- sprintf(paste("Linear mixed model depth ~ diagnosis + (1|donor) per supertype, BH-FDR over %d neuronal",
                      "supertypes; %d/%d reached FDR < %.2f. Violins show inner 95%% of cells; white point = median."),
                n_tested, n_sig, n_tested, SIG_FDR)
cap <- ggdraw() + draw_label(note, size = 8.5, colour = "grey30", x = 0.5, hjust = 0.5)
fig <- plot_grid(leg, grid, cap, ncol = 1, rel_heights = c(0.035, 1, 0.028))

ggsave(file.path(FIGDIR, paste0(FIGSTEM, ".png")), fig, width = 16.8, height = 12.4,
       dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIGDIR, paste0(FIGSTEM, ".pdf")), fig, width = 16.8, height = 12.4,
       bg = "white", device = pdf)
cat("saved", file.path(FIGDIR, paste0(FIGSTEM, ".png/.pdf")), "(composite)\n")
