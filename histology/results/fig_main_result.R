library(tidyverse)
library(cowplot)
library(ggsignif)
library(patchwork)

theme_set(theme_cowplot(font_size = 18))

diag_colors <- c(Control = "#4C72B0", SCHIZ = "#C44E52")

# ── Load data ──
# Single consolidated dataset — filter on vip_pass for primary analysis
df_all <- read_csv("../data/sst_analysis_data.csv", show_col_types = FALSE) %>%
  mutate(
    layer = factor(layer, levels = c("L2/3", "L5/6")),
    diagnosis = factor(diagnosis, levels = c("Control", "MDD", "Bipolar", "SCHIZ"))
  )

# Convert counts to density (cells/mm²)
# FOV = 333 x 333 µm = 0.110889 mm²
fov_area_mm2 <- 0.110889
df_all <- df_all %>% mutate(SST_density = SST / fov_area_mm2)

df_filt <- df_all %>% filter(vip_pass)
n_subj <- n_distinct(df_filt$subject)

# ── Control vs SCHIZ only ──
df_cs <- df_filt %>% filter(diagnosis %in% c("Control", "SCHIZ")) %>% droplevels()
n_ctrl <- n_distinct(df_cs$subject[df_cs$diagnosis == "Control"])
n_schiz <- n_distinct(df_cs$subject[df_cs$diagnosis == "SCHIZ"])

# Per-subject means
subj_means <- df_cs %>%
  group_by(subject, diagnosis) %>%
  summarise(SST_density = mean(SST_density, na.rm = TRUE), .groups = "drop")

subj_layer_means <- df_cs %>%
  group_by(subject, diagnosis, layer) %>%
  summarise(SST_density = mean(SST_density, na.rm = TRUE), .groups = "drop")

density_label <- bquote("SST density (cells/"*mm^2*")")

# ── Panel A: Aggregate ──
p_agg <- ggplot(subj_means, aes(x = diagnosis, y = SST_density, fill = diagnosis)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(color = diagnosis), width = 0.12, size = 2.5, alpha = 0.6) +
  scale_fill_manual(values = diag_colors, guide = "none") +
  scale_color_manual(values = diag_colors, guide = "none") +
  geom_signif(
    comparisons = list(c("Control", "SCHIZ")),
    annotations = "p = 0.071",
    y_position = 43,
    tip_length = 0.02,
    textsize = 5,
    color = "grey40"
  ) +
  labs(
    x = NULL,
    y = density_label,
    title = "Aggregate"
  ) +
  theme(
    axis.text.x = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

# ── Panel B: L2/3 only ──
p_l23 <- ggplot(subj_layer_means %>% filter(layer == "L2/3"),
                aes(x = diagnosis, y = SST_density, fill = diagnosis)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(color = diagnosis), width = 0.12, size = 2.5, alpha = 0.6) +
  scale_fill_manual(values = diag_colors, guide = "none") +
  scale_color_manual(values = diag_colors, guide = "none") +
  geom_signif(
    comparisons = list(c("Control", "SCHIZ")),
    annotations = "p = 0.126",
    y_position = 58,
    tip_length = 0.02,
    textsize = 5,
    color = "grey40"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "L2/3"
  ) +
  theme(
    axis.text.x = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

# ── Panel C: L5/6 only ──
p_l56 <- ggplot(subj_layer_means %>% filter(layer == "L5/6"),
                aes(x = diagnosis, y = SST_density, fill = diagnosis)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(color = diagnosis), width = 0.12, size = 2.5, alpha = 0.6) +
  scale_fill_manual(values = diag_colors, guide = "none") +
  scale_color_manual(values = diag_colors, guide = "none") +
  geom_signif(
    comparisons = list(c("Control", "SCHIZ")),
    annotations = "p = 0.180",
    y_position = 38,
    tip_length = 0.02,
    textsize = 5,
    color = "grey40"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "L5/6"
  ) +
  theme(
    axis.text.x = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

# ── Panel D: Beta coefficients for all diagnoses (aggregate model) ──
# Betas are in count units; convert to density units
count_to_density <- 1 / fov_area_mm2  # ~9.018

coef_df <- tibble(
  diagnosis = factor(c("MDD", "Bipolar", "SCHIZ"), levels = c("MDD", "Bipolar", "SCHIZ")),
  beta = c(0.075, -0.379, -0.597) * count_to_density,
  p = c(0.813, 0.214, 0.071),
  se = abs(c(0.075, -0.379, -0.597) * count_to_density) / abs(qnorm(c(0.813, 0.214, 0.071) / 2))
)

coef_colors <- c(MDD = "#55A868", Bipolar = "#DD8452", SCHIZ = "#C44E52")

p_coef <- ggplot(coef_df, aes(x = diagnosis, y = beta, fill = diagnosis)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_col(alpha = 0.6, width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se),
                width = 0.2, linewidth = 0.8) +
  geom_text(aes(label = sprintf("B = %.1f\np = %.3f", beta, p)),
            vjust = ifelse(coef_df$beta < 0, 1.3, -0.3),
            size = 4.5, lineheight = 0.85) +
  scale_fill_manual(values = coef_colors, guide = "none") +
  labs(
    x = "Diagnosis (vs Control)",
    y = bquote("SST density coefficient (cells/"*mm^2*")"),
    title = "Diagnosis effects"
  ) +
  theme(
    axis.text.x = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

# ── Combine ──
p_combined <- (p_agg | p_l23 | p_l56 | p_coef) +
  plot_layout(widths = c(1, 1, 1, 1)) +
  plot_annotation(
    title = "SST interneuron density in sgACC by diagnosis and cortical layer",
    subtitle = sprintf("VIP-filtered (N = %d; %d Control, %d SCZ) | Mixed model: SST ~ diagnosis + covariates + (1|subject)",
                       n_subj, n_ctrl, n_schiz),
    theme = theme(
      plot.title = element_text(size = 20, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "grey40")
    )
  )

ggsave("fig_main_result.png", p_combined, width = 20, height = 6, dpi = 150)
cat("Saved fig_main_result.png\n")

ggsave("fig_main_result.pdf", p_combined, width = 20, height = 6)
cat("Saved fig_main_result.pdf\n")
