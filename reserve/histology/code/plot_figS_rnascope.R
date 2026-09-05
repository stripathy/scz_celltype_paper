# ============================================================================
# Supplementary figure: RNAscope SST-positive cell density in sgACC
# (re-analysis of Arbabi et al. 2025), drawn in the main-figure conventions
# (theme_cowplot 7 pt, Control blue / SCZ orange as in Fig. 2, lowercase bold
# panel letters, 7.1-in width).
#
#   a  per-subject SST+ density, Control vs SCZ, all frames and by layer
#      (P from the mixed model, read from results/sst_mixed_models.csv)
#   b  representative lipofuscin-suppressed composites, 2 x 2 (Control | SCZ
#      columns; L2/3 | L5/6 rows) with the manually annotated SST+ / VIP+ cells
#   c  mixed-model diagnosis coefficients (cells/mm2 vs Control, 95% CI) for
#      MDD, Bipolar and SCZ, all layers / L2/3 / L5/6
# Layout: a over c in the left column, b in the right column.
#
# Inputs (all produced by code/fit_sst_mixed_models.py and
# code/extract_micrograph_panels.py; nothing is hard-coded here):
#   results/sst_mixed_models.csv, results/sst_subject_density.csv,
#   results/microscopy/panels/panel*.png, results/microscopy/marker_coordinates.csv
# Output: results/figS_rnascope_sst.{png,pdf}
# Run from histology/:  Rscript code/plot_figS_rnascope.R
# ============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(cowplot); library(ggsignif); library(png); library(grid)
})

# ---- constants (match transcriptomic/scripts/09_composite_figure.R) --------
BASE       <- 7
AXIS_TITLE <- BASE - 0.5
AXIS_TEXT  <- BASE - 1
FIG_W      <- 7.1
DX_COL     <- c(Control = "#0072B2", SCZ = "#D55E00")   # Fig. 2 palette
COL_GREY   <- "grey35"
FRAME_PX   <- 800                       # counting frame, pixels
FRAME_UM   <- 800 / 1024 * 333          # ~260 um per frame side
SCALEBAR_UM <- 100

theme_panel <- function() {
  theme_cowplot(font_size = BASE) +
    theme(axis.title = element_text(size = AXIS_TITLE),
          axis.text  = element_text(size = AXIS_TEXT),
          axis.line  = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          strip.background = element_blank(),
          strip.text = element_text(size = BASE, margin = margin(b = 2)),
          legend.position = "none",
          plot.margin = margin(3, 4, 2, 3))
}

# ---- data --------------------------------------------------------------------
models <- read_csv("results/sst_mixed_models.csv", show_col_types = FALSE)
subj   <- read_csv("results/sst_subject_density.csv", show_col_types = FALSE) %>%
  filter(diagnosis %in% c("Control", "SCZ")) %>%
  mutate(diagnosis = factor(diagnosis, levels = c("Control", "SCZ")),
         layer = factor(layer, levels = c("All layers", "L2/3", "L5/6")))
n_ctrl <- n_distinct(subj$subject[subj$diagnosis == "Control"])
n_scz  <- n_distinct(subj$subject[subj$diagnosis == "SCZ"])
cat(sprintf("Control n = %d, SCZ n = %d\n", n_ctrl, n_scz))

# ---- a: per-subject density, three facets sharing one y-axis ---------------
ann <- models %>% filter(diagnosis == "SCZ") %>%
  transmute(layer = factor(model, levels = levels(subj$layer)),
            annotations = sprintf("italic(P)=='%.3f'", p),
            xmin = "Control", xmax = "SCZ")
ymax <- max(subj$density)
ann$y_position <- ymax * 1.04

p_box <- ggplot(subj, aes(diagnosis, density)) +
  geom_boxplot(aes(fill = diagnosis), width = 0.62, outlier.shape = NA,
               alpha = 0.55, linewidth = 0.4) +
  geom_jitter(aes(fill = diagnosis), shape = 21, colour = "grey25", size = 1.4,
              stroke = 0.3, width = 0.13, height = 0, alpha = 0.9) +
  geom_signif(data = ann, aes(xmin = xmin, xmax = xmax, annotations = annotations,
                              y_position = y_position),
              manual = TRUE, parse = TRUE, tip_length = 0.02,
              textsize = BASE * 0.32, vjust = -0.1, inherit.aes = FALSE) +
  facet_wrap(~ layer, nrow = 1) +
  scale_fill_manual(values = DX_COL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.16)), n.breaks = 5) +
  labs(x = NULL, y = expression("SST"^"+"~"cell density (cells/mm"^2*")")) +
  theme_panel()

# ---- b: diagnosis coefficients (cells/mm2 vs Control), by layer -------------
coef <- models %>%
  mutate(diagnosis = factor(diagnosis, levels = c("MDD", "Bipolar", "SCZ")),
         model = factor(model, levels = c("L5/6", "L2/3", "All layers")),
         emph = if_else(diagnosis == "SCZ", "SCZ", "other"))
dodge <- position_dodge(width = 0.6)
p_lab <- coef %>% filter(model == "All layers") %>%
  mutate(label = sprintf("italic(P)=='%.3f'", p))
xr <- range(c(coef$ci_lo_density, coef$ci_hi_density))

p_coef <- ggplot(coef, aes(x = beta_density, y = diagnosis, group = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.3) +
  geom_errorbar(aes(xmin = ci_lo_density, xmax = ci_hi_density, colour = emph),
                width = 0, linewidth = 0.4, position = dodge, orientation = "y") +
  geom_point(aes(shape = model, colour = emph, fill = emph), size = 1.7,
             stroke = 0.4, position = dodge) +
  geom_text(data = p_lab, aes(x = xr[2] + 0.06 * diff(xr), label = label),
            parse = TRUE, hjust = 0, size = BASE * 0.30, colour = "grey30") +
  scale_shape_manual(values = c(`All layers` = 21, `L2/3` = 24, `L5/6` = 25),
                     breaks = c("All layers", "L2/3", "L5/6"), name = NULL) +
  scale_colour_manual(values = c(SCZ = DX_COL[["SCZ"]], other = COL_GREY), guide = "none") +
  scale_fill_manual(values = c(SCZ = DX_COL[["SCZ"]], other = COL_GREY), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.28))) +
  labs(x = expression("Difference vs Control (cells/mm"^2*")"), y = NULL) +
  theme_panel() +
  theme(legend.position = "top", legend.justification = "right", legend.direction = "horizontal",
        legend.text = element_text(size = AXIS_TEXT), legend.key.height = unit(7, "pt"),
        legend.key.width = unit(9, "pt"), legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, -4, 0), legend.background = element_blank())

# ---- c: representative composites with manual markers -----------------------
# 2 x 2 grid: columns Control | SCZ, rows L2/3 (top) | L5/6 (bottom).
markers <- read_csv("results/microscopy/marker_coordinates.csv", show_col_types = FALSE)
panel_meta <- tibble(
  panel = 0:3,
  file  = c("panel0_Control_L23.png", "panel1_SCZ_L23.png",
            "panel2_Control_L56.png", "panel3_SCZ_L56.png"))
MK_COL <- c(sst = "#00FF66", vip = "#FF3355")

CH_COL <- c(DAPI = "#6D8CFF", SST = "#00FF66", VIP = "#FF3355")   # on-image channel key

build_micrograph <- function(i, scalebar = FALSE, key = FALSE) {
  m   <- panel_meta[panel_meta$panel == i, ]
  img <- readPNG(file.path("results/microscopy/panels", m$file))
  mk  <- markers %>% filter(panel == i)
  p <- ggplot() +
    annotation_raster(img, xmin = 0, xmax = FRAME_PX, ymin = 0, ymax = FRAME_PX) +
    geom_point(data = mk, aes(pixel_x, FRAME_PX - pixel_y, fill = type),
               shape = 21, colour = "white", size = 1.25, stroke = 0.35) +
    scale_fill_manual(values = MK_COL, guide = "none") +
    coord_fixed(xlim = c(0, FRAME_PX), ylim = c(0, FRAME_PX), expand = FALSE) +
    theme_void(base_size = BASE) +
    theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5))
  if (scalebar) {
    bar_px <- SCALEBAR_UM / FRAME_UM * FRAME_PX
    p <- p +
      annotate("segment", x = 30, xend = 30 + bar_px, y = 40, yend = 40,
               colour = "white", linewidth = 0.9) +
      annotate("text", x = 30 + bar_px / 2, y = 78, label = sprintf("%d \u00b5m", SCALEBAR_UM),
               colour = "white", size = BASE * 0.30)
  }
  if (key) {
    # channel key, top-right, each channel name in its own colour
    ch <- tibble(lab = c("DAPI", "SST mRNA", "VIP mRNA"),
                 col = CH_COL, y = FRAME_PX - c(35, 80, 125))
    p <- p + annotate("text", x = FRAME_PX - 22, y = ch$y, label = ch$lab,
                      colour = ch$col, hjust = 1, size = BASE * 0.30, fontface = "bold")
    # marker key, top-left: the two marker glyphs with white labels
    mkk <- tibble(type = c("sst", "vip"), y = FRAME_PX - c(35, 80),
                  lab = c("SST+ cell", "VIP+ cell"))
    p <- p +
      geom_point(data = mkk, aes(x = 30, y = y, fill = type), shape = 21,
                 colour = "white", size = 1.25, stroke = 0.35) +
      annotate("text", x = 55, y = mkk$y, label = mkk$lab, colour = "white",
               hjust = 0, size = BASE * 0.30)
  }
  p
}
hdr <- function(txt, angle = 0) ggdraw() +
  draw_label(txt, size = BASE, angle = angle, x = 0.5, y = 0.5)
LBL_W <- 0.07   # width of the row-label strip relative to one image
row_hdr <- plot_grid(NULL, hdr("Control"), hdr("SCZ"), nrow = 1, rel_widths = c(LBL_W, 1, 1))
row_l23 <- plot_grid(hdr("L2/3", 90), build_micrograph(0, scalebar = TRUE, key = TRUE), build_micrograph(1),
                     nrow = 1, rel_widths = c(LBL_W, 1, 1))
row_l56 <- plot_grid(hdr("L5/6", 90), build_micrograph(2), build_micrograph(3),
                     nrow = 1, rel_widths = c(LBL_W, 1, 1))
p_micro <- plot_grid(row_hdr, row_l23, row_l56, ncol = 1, rel_heights = c(0.09, 1, 1))

# ---- assemble ----------------------------------------------------------------
# Left column: a over b. Right column: c (2 x 2 composites).
# Panel letters: a (density), b (micrographs, right column), c (coefficients).
left  <- plot_grid(p_box, p_coef, ncol = 1, rel_heights = c(1, 1),
                   labels = c("a", "c"), label_size = 8, label_fontface = "bold")
right <- plot_grid(p_micro, labels = "b", label_size = 8, label_fontface = "bold")
FIG_H <- 3.95
full <- plot_grid(left, right, ncol = 2, rel_widths = c(1, 1.02))
dir.create("results", showWarnings = FALSE)
ggsave("results/figS_rnascope_sst.png", full, width = FIG_W, height = FIG_H, dpi = 400, bg = "white")
ggsave("results/figS_rnascope_sst.pdf", full, width = FIG_W, height = FIG_H, bg = "white")
cat("wrote results/figS_rnascope_sst.{png,pdf}\n")
