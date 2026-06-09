#!/usr/bin/env Rscript
# ============================================================================
# 16_de_vs_proportion.R — Does the number of snRNA-seq META DE genes (FDR<0.10)
# per subclass scale with that cell type's PROPORTION?
#
# Proportion is taken from the XENIUM data for now (a proxy / stand-in); it will
# be swapped for the snRNA-seq cohort proportion when that per-donor data is in
# hand (edit the "proportion" block).
#
# CAVEAT (read before interpreting): the DE-gene COUNT is power-dependent — more
# cells -> more statistical power -> more DE detected — so a POSITIVE relationship
# with proportion is EXPECTED and is partly a detection artifact, not purely
# biology. The interesting signal is the DEVIATION from the trend (cell types with
# many more / fewer DE genes than their abundance predicts, e.g. excitatory
# neurons in Ruzicka 2024).
#
# Inputs : data/DE_genes_all_cells_scz.csv                       (meta DE: cell_type, padj)
#          ../spatial/output/crumblr/crumblr_input_subclass_corr.csv  (Xenium per-donor counts)
# Outputs: results/figures/de_vs_proportion_subclass.{png,pdf}
#          results/tables/de_vs_proportion_subclass.csv
# ============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(ggrepel)
  library(cowplot); library(scales)
})

# --- y: meta DE-gene count (FDR<0.10) per subclass ---
m <- read_csv("data/DE_genes_all_cells_scz.csv", show_col_types = FALSE)
nde <- m |> group_by(cell_type) |>
  summarise(n_de = sum(padj < 0.10, na.rm = TRUE), n_tested = dplyr::n(), .groups = "drop")

# --- x: cell-type proportion (XENIUM, mean per-donor) — swap this block for snRNA later ---
ct_map <- c("Astrocyte"="Astro","L2/3 IT"="L2_3 IT","L5/6 NP"="L5_6 NP",
            "Microglia-PVM"="Micro-PVM","Oligodendrocyte"="Oligo","Endothelial"="Endo")
cr <- read_csv("../spatial/output/crumblr/crumblr_input_subclass_corr.csv", show_col_types = FALSE)
prop <- cr |>
  mutate(p = count / total,
         cell_type = ifelse(celltype %in% names(ct_map), ct_map[celltype], celltype)) |>
  group_by(cell_type) |> summarise(prop_xen = mean(p), .groups = "drop")

# --- join + cell class ---
EXC <- c("L2_3 IT","L4 IT","L5 IT","L5 ET","L5_6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")
INH <- c("Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Sncg","Sst","Sst Chodl","Vip","Chandelier")
GLI <- c("Astro","Oligo","OPC","Micro-PVM","Endo","VLMC")
CLASS_COL <- c(Excitatory = "#117733", Inhibitory = "#882255", Glia = "#DDCC77")
d <- inner_join(nde, prop, by = "cell_type") |>
  mutate(class = dplyr::case_when(cell_type %in% EXC ~ "Excitatory",
                                  cell_type %in% INH ~ "Inhibitory",
                                  cell_type %in% GLI ~ "Glia", TRUE ~ "Other"))
unmatched <- setdiff(union(nde$cell_type, prop$cell_type), d$cell_type)
dropped_zero <- d$cell_type[d$n_de < 1]            # unpowered types with no detected DE
d <- d |> filter(n_de >= 1)                        # analysis set: cell types with >=1 meta DE gene
write_csv(d |> arrange(desc(prop_xen)), "results/tables/de_vs_proportion_subclass.csv")

# --- correlations (Spearman is scale-free; Pearson on log10 matches the fit line) ---
rs  <- cor(d$prop_xen, d$n_de, method = "spearman")
rlp <- cor(log10(d$prop_xen), d$n_de)
cat(sprintf("n=%d subclasses (>=1 DE)%s%s\nSpearman(prop, n_DE)=%.2f | Pearson(log10 prop, n_DE)=%.2f\n",
            nrow(d),
            if (length(dropped_zero)) paste0("; dropped 0-DE: ", paste(dropped_zero, collapse=", ")) else "",
            if (length(unmatched)) paste0("; unmatched: ", paste(unmatched, collapse=", ")) else "",
            rs, rlp))
print(as.data.frame(d |> arrange(desc(prop_xen)) |>
        transmute(cell_type, class, prop_pct = round(100*prop_xen,2), n_de, n_tested)), row.names = FALSE)

# --- plot: x = proportion (log10), y = n DE genes ---
p <- ggplot(d, aes(prop_xen, n_de)) +
  geom_smooth(method = "lm", se = FALSE, colour = alpha("grey25", 0.3),
              linewidth = 0.8, formula = y ~ x) +
  geom_point(aes(colour = class), size = 2.6, alpha = 0.9) +
  geom_text_repel(aes(label = gsub("_", "/", cell_type)), size = 3.0, seed = 1,
                  max.overlaps = Inf, min.segment.length = 0, segment.size = 0.25,
                  box.padding = 0.4, colour = "grey20") +
  annotate("text", x = max(d$prop_xen), y = min(d$n_de),
           label = sprintf("Spearman~rho == %.2f", rs), parse = TRUE, hjust = 1, vjust = 0, size = 3.6) +
  scale_colour_manual(values = CLASS_COL, name = NULL) +
  scale_x_log10(breaks = c(0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3),
                labels = c("0.03%", "0.1%", "0.3%", "1%", "3%", "10%", "30%")) +
  annotation_logticks(sides = "b", colour = "grey55", linewidth = 0.3,
                      short = unit(0.04, "cm"), mid = unit(0.07, "cm"), long = unit(0.12, "cm")) +
  labs(x = "Cell-type proportion (Xenium, mean per donor)",
       y = "snRNA-seq meta DE genes (FDR < 0.10)",
       title = "DE-gene count vs cell-type proportion (subclass)") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = c(0.02, 0.97), legend.justification = c(0, 1),
        legend.text = element_text(size = 10))
ggsave("results/figures/de_vs_proportion_subclass.png", p, width = 6.6, height = 5.0, dpi = 200, bg = "white")
ggsave("results/figures/de_vs_proportion_subclass.pdf", p, width = 6.6, height = 5.0, bg = "white")
cat("Saved results/figures/de_vs_proportion_subclass.{png,pdf} + results/tables/de_vs_proportion_subclass.csv\n")
