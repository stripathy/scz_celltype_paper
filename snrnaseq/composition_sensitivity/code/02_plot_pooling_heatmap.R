# Supplementary figure: robustness of the composition results to pooling strategy.
#   a  heatmap of SCZ beta for the 109 neuronal supertypes (Fig. 3a order) under five
#      pooling strategies, with FDR tier marks, an I2 row, a subclass strip and every
#      supertype named under its column (Fig. 3a styling: red bold FDR < 0.10 and red
#      italic FDR < 0.20 in the paper's fixed-effect meta-analysis)
#   b  beta under each alternative vs the paper's fixed-effect meta-analysis
#   c  leave-one-dataset-out for the supertypes at FDR < 0.20 in the paper
# Input : results/composition_pooling_sensitivity.csv, results/composition_lodo.csv,
#         data/neuronal_supertypes_109.csv
# Output: manuscript/figures/supplementary/S06_composition_pooling.{png,pdf}
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(cowplot); library(jsonlite)
})
args <- commandArgs(trailingOnly = FALSE)
fp   <- sub("^--file=", "", args[grep("^--file=", args)])
HERE <- if (length(fp)) normalizePath(file.path(dirname(fp), "..")) else normalizePath(getwd())
ROOT <- normalizePath(file.path(HERE, "..", ".."))          # repo root
BASE <- 7; LAB <- 1.9
UP <- "#D55E00"; DOWN <- "#0072B2"; RED <- "#D7191C"
THR <- 0.20                                  # robustness threshold used throughout

res  <- read_csv(file.path(HERE, "results/composition_pooling_sensitivity.csv"), show_col_types = FALSE)
lodo <- read_csv(file.path(HERE, "results/composition_lodo.csv"), show_col_types = FALSE)
order_st <- read_csv(file.path(HERE, "data/neuronal_supertypes_109.csv"), show_col_types = FALSE)$supertype
pal <- fromJSON(file.path(HERE, "data/seaad_supertype_colors.json"))   # copy of genetics/data/ (untracked there)

METHODS <- c("published FE meta" = "FE meta-analysis (paper)",
             "RE meta" = "RE meta-analysis",
             "mega: dataset random intercept" = "Mega, dataset random intercept",
             "mega: random intercept + SCZ slope" = "Mega, random intercept + SCZ slope")
tier <- function(f) case_when(is.na(f) ~ "", f < 0.01 ~ "***", f < 0.05 ~ "**", f < 0.10 ~ "*", f < THR ~ "+", TRUE ~ "")

xs <- tibble(supertype = order_st, x = seq_along(order_st), subclass = sub("(_\\d+)+$", "", order_st))
res <- res |> filter(method %in% names(METHODS)) |>
  mutate(method = factor(METHODS[method], levels = rev(METHODS))) |>
  left_join(xs, by = "supertype") |>
  mutate(stars = tier(fdr), beta_c = pmax(pmin(beta, 0.5), -0.5))
paper <- res |> filter(method == METHODS["published FE meta"]) |> select(supertype, x, fdr_paper = fdr, beta_paper = beta)
xs <- xs |> left_join(paper |> select(supertype, fdr_paper), by = "supertype") |>
  mutate(face = case_when(fdr_paper < 0.10 ~ "bold", fdr_paper < THR ~ "italic", TRUE ~ "plain"),
         col  = ifelse(fdr_paper < THR, RED, "grey20"))
HL <- xs$supertype[xs$fdr_paper < THR]                       # the paper's FDR < 0.20 set
n_sig <- res |> group_by(method) |> summarise(n = sum(fdr < THR, na.rm = TRUE), .groups = "drop")
NX <- length(order_st); NM <- length(METHODS)

# subclass blocks: contiguous runs in the taxonomy order (a subclass may appear in more than one run)
runs <- xs |> mutate(run = cumsum(subclass != lag(subclass, default = first(subclass)))) |>
  group_by(run, subclass) |> summarise(xmin = min(x) - 0.5, xmax = max(x) + 0.5, xmid = mean(x), n = n(), .groups = "drop") |>
  mutate(col = sapply(subclass, function(s) pal[[order_st[xs$subclass == s][1]]] %||% "grey70"),
         label = ifelse(n >= 4, subclass, ""))
I2 <- res |> filter(method == METHODS["RE meta"]) |> select(x, I2)

# ---------------------------------------------------------------- a: heatmap
# (The "mega: dataset fixed" model is fit in 01 but not shown: it is numerically the same as the
#  random-intercept model, beta r = 0.9996, identical FDR sets at every threshold.)
# Rows, top to bottom: the seven per-dataset crumblr estimates and the Xenium estimate
# (dataset-independent effects; dot = nominal P < 0.05), the alternative pooling
# strategies, the paper's fixed-effect meta-analysis, and between-dataset I2.
per <- read_csv(file.path(HERE, "results/composition_per_dataset_estimates.csv"), show_col_types = FALSE)
xen <- read_csv(file.path(ROOT, "spatial/output/crumblr/crumblr_results_supertype_neuronal.csv"),
                show_col_types = FALSE) |> transmute(supertype = celltype, beta = logFC, p = P.Value, dataset = "Xenium")
DS <- c("MSSM 2", "HBCC", "Fröhlich", "MSSM 1", "McLean", "Batiuk", "Multiome", "Xenium")
# row labels carry the control / SCZ donor counts of each dataset
cnt <- read_csv(file.path(HERE, "data/neuron_counts_469donors.csv"), show_col_types = FALSE) |>
  count(dataset, dx) |> pivot_wider(names_from = dx, values_from = n)
xen_cnt <- read_csv(file.path(ROOT, "spatial/output/crumblr/crumblr_input_supertype_neuronal.csv"),
                    show_col_types = FALSE) |> distinct(donor, diagnosis) |> count(diagnosis) |>
  pivot_wider(names_from = diagnosis, values_from = n) |> mutate(dataset = "Xenium")
cnt <- bind_rows(cnt, xen_cnt)
row_lab <- setNames(sprintf("%s (%d control, %d SCZ)", cnt$dataset, cnt$Control, cnt$SCZ), cnt$dataset)
stopifnot(all(DS %in% names(row_lab)))
single <- bind_rows(per |> select(supertype, beta, p, dataset), xen) |>
  filter(supertype %in% order_st) |>
  mutate(row = row_lab[dataset], block = "dataset",
         mark = ifelse(!is.na(p) & p < 0.05, "\u2022", ""))
ALT <- c("RE meta-analysis", "Mega, dataset random intercept", "Mega, random intercept + SCZ slope")
pooled <- res |> mutate(row = as.character(method), mark = stars,
                        block = ifelse(row == METHODS["published FE meta"], "paper", "alt")) |>
  select(supertype, beta, p, fdr, row, block, mark, x)
single <- single |> left_join(xs |> select(supertype, x), by = "supertype")
# numeric y with gaps between blocks (top = first)
rows_top <- c(row_lab[DS], ALT, METHODS["published FE meta"])
ypos <- setNames(rev(seq_along(rows_top)), rows_top)          # 1 = bottom
ypos[row_lab["Xenium"]] <- ypos[row_lab["Xenium"]] - 0.35     # Xenium set off from the snRNA-seq datasets
ypos[ALT] <- ypos[ALT] - 0.95 - 0.35                           # heading gap below dataset block
ypos[METHODS["published FE meta"]] <- ypos[METHODS["published FE meta"]] - 0.95 - 0.35 - 0.45  # FE row set off
hm_df <- bind_rows(single |> select(supertype, x, beta, row, mark, block),
                   pooled |> select(supertype, x, beta, row, mark, block)) |>
  mutate(y = ypos[row], beta_c = pmax(pmin(beta, 0.5), -0.5))
Y_I2 <- min(ypos) - 1.05; Y_STRIP <- max(ypos) + 0.95; Y_LAB <- Y_I2 - 0.55
I2 <- res |> filter(method == METHODS["RE meta"]) |> select(x, I2)
n_sig_hm <- n_sig |> mutate(y = ypos[as.character(method)])
NX <- length(order_st)
hm <- ggplot(hm_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = beta_c), colour = "white", linewidth = 0.15, height = 0.92) +
  geom_text(aes(label = mark), size = 1.45, vjust = 0.72, colour = "black") +
  geom_vline(xintercept = runs$xmin[-1], colour = "white", linewidth = 0.7) +
  # block headings, right-aligned in the row-label column
  annotate("text", x = -0.6, y = Y_STRIP, label = "Per dataset", hjust = 1, size = LAB, fontface = "italic", colour = "grey15") +
  annotate("text", x = -0.6, y = max(ypos[ALT]) + 0.95, label = "Pooled", hjust = 1, size = LAB, fontface = "italic", colour = "grey15") +
  # I2 row
  geom_tile(data = I2, aes(x = x, y = Y_I2, alpha = I2), fill = "grey20", height = 0.55,
            colour = "white", linewidth = 0.15, inherit.aes = FALSE) +
  # subclass strip above, names where wide enough
  geom_rect(data = runs, aes(xmin = xmin, xmax = xmax, ymin = Y_STRIP - 0.22, ymax = Y_STRIP + 0.22),
            fill = runs$col, colour = "white", linewidth = 0.3, inherit.aes = FALSE) +
  geom_text(data = runs, aes(x = xmid, y = Y_STRIP, label = label), size = 1.7, colour = "white",
            fontface = "bold", inherit.aes = FALSE) +
  # supertype names under their columns (Fig. 3a styling)
  geom_text(data = xs, aes(x = x, y = Y_LAB, label = supertype), angle = 90, hjust = 1, vjust = 0.5,
            size = 1.32, colour = xs$col, fontface = xs$face, inherit.aes = FALSE) +
  # n at FDR < THR for the pooled rows
  geom_text(data = n_sig_hm, aes(x = NX + 1.6, y = y, label = n), hjust = 0, size = LAB, colour = "grey15", inherit.aes = FALSE) +
  annotate("text", x = NX + 1.6, y = max(ypos[ALT]) + 0.95,
           label = sprintf("n FDR < %.2f", THR), hjust = 0, size = LAB - 0.2, colour = "grey15") +
  scale_fill_gradient2(low = DOWN, mid = "white", high = UP, midpoint = 0, limits = c(-0.5, 0.5),
                       breaks = c(-0.5, -0.25, 0, 0.25, 0.5), labels = c("-0.5", "-0.25", "0", "0.25", "0.5"),
                       name = expression("SCZ coefficient "*beta*" (clipped at ±0.5)"), na.value = "grey92") +
  scale_alpha_continuous(range = c(0.03, 1), limits = c(0, 100), name = expression(I^2~"(%)")) +
  scale_x_continuous(expand = expansion(add = c(1.2, 8))) +
  scale_y_continuous(breaks = c(ypos, Y_I2), labels = c(names(ypos), "Heterogeneity (I²)"),
                     expand = expansion(add = 0)) +
  coord_cartesian(ylim = c(Y_LAB - 0.2, Y_STRIP + 0.35), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_cowplot(font_size = BASE) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        axis.text.y = element_text(size = BASE - 1),
        legend.position = "top", legend.justification = "left", legend.box = "horizontal",
        legend.title = element_text(size = BASE - 1), legend.text = element_text(size = BASE - 1.5),
        legend.key.height = unit(6, "pt"), legend.key.width = unit(12, "pt"),
        legend.box.margin = margin(0, 0, 2, 0), legend.margin = margin(0, 14, 0, 0),
        plot.margin = margin(2, 6, 30, 6))

# ---------------------------------------------------------------- b: agreement scatters
WRAP <- c("RE meta-analysis" = "RE meta-analysis",
          "Mega, dataset random intercept" = "Mega, dataset\nrandom intercept",
          "Mega, random intercept + SCZ slope" = "Mega, random\nintercept + SCZ slope")
wide <- res |> select(supertype, method, beta) |> pivot_wider(names_from = method, values_from = beta)
ref  <- METHODS["published FE meta"]
sc_df <- bind_rows(lapply(setdiff(METHODS, ref), function(m)
  tibble(supertype = wide$supertype, x = wide[[ref]], y = wide[[m]], method = WRAP[[m]]))) |>
  mutate(method = factor(method, levels = WRAP[setdiff(METHODS, ref)]), hl = supertype %in% HL)
sc_stats <- sc_df |> group_by(method) |>
  summarise(rho = cor(x, y, method = "spearman", use = "complete.obs"), .groups = "drop") |>
  mutate(lab = sprintf("rho == %.2f", rho))
sc <- ggplot(sc_df, aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.25) + geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.25) +
  geom_point(data = filter(sc_df, !hl), colour = "grey65", size = 0.9, alpha = 0.9) +
  geom_point(data = filter(sc_df, hl), colour = RED, size = 1.2) +
  ggrepel::geom_text_repel(data = filter(sc_df, hl), aes(label = supertype), size = 1.6, colour = "grey15",
                           box.padding = 0.12, point.padding = 0.1, segment.size = 0.15, segment.colour = "grey60",
                           min.segment.length = 0, max.overlaps = Inf, seed = 2) +
  geom_text(data = sc_stats, aes(x = -Inf, y = Inf, label = lab), parse = TRUE, hjust = -0.15, vjust = 1.4,
            size = 2.3, inherit.aes = FALSE) +
  facet_wrap(~ method, nrow = 1, scales = "free_y") +
  labs(x = expression(beta*", fixed-effect meta-analysis (paper)"), y = expression(beta*", alternative strategy")) +
  theme_cowplot(font_size = BASE) +
  theme(strip.background = element_blank(), strip.text = element_text(size = BASE - 1, lineheight = 0.9),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        panel.spacing = unit(8, "pt"), plot.margin = margin(4, 6, 2, 6))

# ---------------------------------------------------------------- c: leave-one-dataset-out
# Cell type along x (paper order), beta on y. Black diamond + bar = full seven-dataset FE
# estimate with its 95% CI; grey circles = the seven estimates with one dataset left out
# (no dataset identity encoded; which omission attenuates most is stated in the legend).
full_c <- res |> filter(method == METHODS["published FE meta"], supertype %in% HL) |>
  transmute(ct = factor(supertype, levels = HL), beta, lo95 = beta - 1.96 * se, hi95 = beta + 1.96 * se, fdr)
lo_c <- lodo |> filter(supertype %in% HL) |> transmute(ct = factor(supertype, levels = HL), beta, dropped)
# every leave-one-out estimate inside the full estimate's 95% CI?
chk <- lo_c |> left_join(full_c |> select(ct, lo95, hi95), by = "ct") |> mutate(inside = beta > lo95 & beta < hi95)
cat(sprintf("leave-one-out estimates inside the full 95%% CI: %d of %d; signs preserved: %s\n",
            sum(chk$inside), nrow(chk), all(sign(chk$beta) == sign((chk$lo95 + chk$hi95) / 2))))
lodo_p <- ggplot(full_c, aes(x = ct, y = beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.3) +
  geom_linerange(aes(ymin = lo95, ymax = hi95), colour = "black", linewidth = 0.35) +
  geom_point(data = lo_c, aes(x = ct, y = beta), shape = 21, size = 1.5, colour = "grey25", fill = "white", stroke = 0.4,
             position = position_jitter(width = 0.12, height = 0, seed = 3), inherit.aes = FALSE) +
  geom_point(shape = 18, size = 2.6, colour = "black") +
  annotate("text", x = 0.55, y = max(full_c$hi95), hjust = 0, vjust = 1, size = LAB, colour = "grey15", lineheight = 0.95,
           label = "diamond, bar: all seven datasets (95% CI)\ncircles: one dataset left out") +
  labs(x = NULL, y = expression("SCZ coefficient "*beta)) +
  theme_cowplot(font_size = BASE) +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        axis.text.x = element_text(size = BASE - 1, colour = xs$col[match(HL, xs$supertype)],
                                   face = xs$face[match(HL, xs$supertype)]),
        plot.margin = margin(4, 8, 2, 6))

fig <- plot_grid(hm, sc, lodo_p, ncol = 1, rel_heights = c(1.78, 1, 0.78),
                 labels = c("a", "b", "c"), label_size = 8, label_fontface = "bold")
# Written directly into the supplementary submission folder as an original,
# per manuscript/figures/supplementary/README.md; the S-number lives in FIGSTEM.
SUPPFIG <- file.path(ROOT, "manuscript/figures/supplementary")
FIGSTEM <- "S06_composition_pooling"
ggsave(file.path(SUPPFIG, paste0(FIGSTEM, ".png")), fig, width = 7.1, height = 8.6, dpi = 400, bg = "white")
ggsave(file.path(SUPPFIG, paste0(FIGSTEM, ".pdf")), fig, width = 7.1, height = 8.6, bg = "white")
cat(sprintf("wrote %s/%s.{png,pdf}\n", SUPPFIG, FIGSTEM))
cat("paper FDR <", THR, "set:", HL, "\n")
