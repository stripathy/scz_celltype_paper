#!/usr/bin/env Rscript
# Render the new Figure 4 elements as STANDALONE panels (layout-independent).
# Builders live in fig4_new_panels.R; the multipanel assembly in
# scz_sst_hcn1_story.R sources the same file.

FIG4_REPO <- "/Users/shreejoy/Github/scz_celltype_paper"
source(file.path(FIG4_REPO, "genetics/scripts/figures/fig4_new_panels.R"))
OUT <- file.path(FIG4_REPO, "genetics/results/figures")

# PNG (preview) + PDF (submission) + SVG (editable vector), as in the assembly script.
save_panel <- function(p, name, w, h) {
  stem <- file.path(OUT, name)
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 400, bg = "white")
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h, bg = "white")
  ggsave(paste0(stem, ".svg"), p, width = w, height = h, bg = "white",
         device = svglite::svglite)
  cat("wrote", name, "(png/pdf/svg)\n")
}

save_panel(build_marker_volcano(),  "fig4_volcano_vulnerable_markers", 3.6, 3.0)
save_panel(build_violin("CALB1"),   "fig4_violin_calb1",               1.7, 2.2)
save_panel(build_violin("SST"),     "fig4_violin_sst",                 1.7, 2.2)
save_panel(build_ad_concordance(),  "fig4_scz_vs_ad_concordance_sst",  3.3, 3.0)
cat("\nDone.\n")

# Comparison variant: cell-level Wilcoxon on the y-axis (see build_marker_volcano).
save_panel(build_marker_volcano("cell"),
           "fig4_volcano_vulnerable_markers_celllevel", 3.6, 3.0)
