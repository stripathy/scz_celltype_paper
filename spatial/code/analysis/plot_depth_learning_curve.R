#!/usr/bin/env Rscript
# Depth-model donor learning curve (low-CPS cohort): test R2 vs # training donors.
# Reads output/depth_validation/lowcps/learning_curve.csv. House style.
suppressPackageStartupMessages({library(readr);library(dplyr);library(ggplot2);library(cowplot)})
cowplot::set_null_device("agg"); BASE <- 7
DD <- "output/depth_validation/lowcps"
lc <- read_csv(file.path(DD, "learning_curve.csv"), show_col_types = FALSE)
ref <- 0.928  # full low-CPS held-out CV R2
s <- lc %>% group_by(k_train) %>% summarise(m = mean(r2), sd = sd(r2), .groups = "drop") %>%
  mutate(lo = m - sd, hi = m + sd)

g <- ggplot() +
  geom_hline(yintercept = ref, linetype = "dashed", color = "#D55E00", linewidth = 0.3) +
  annotate("text", x = 9.6, y = ref + 0.004, label = "all-13 CV (0.928)", hjust = 1,
           size = 2.3, color = "#D55E00") +
  geom_ribbon(data = s, aes(k_train, ymin = lo, ymax = hi), fill = "#0072B2", alpha = 0.15) +
  geom_jitter(data = lc, aes(k_train, r2), width = 0.12, height = 0, size = 0.7,
              color = "#0072B2", alpha = 0.4) +
  geom_line(data = s, aes(k_train, m), color = "#0072B2", linewidth = 0.6) +
  geom_point(data = s, aes(k_train, m), color = "#0072B2", size = 1.6) +
  annotate("rect", xmin = 8.5, xmax = 10.5, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.3) +
  annotate("text", x = 9.5, y = min(lc$r2) + 0.01, label = "proposed\n(train 9-10)", size = 2.2, color = "grey30") +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "# training donors", y = expression("Held-out test"~italic(R)^2),
       subtitle = "Depth-model learning curve (low-CPS SEA-AD donors)") +
  theme_cowplot(font_size = BASE) +
  theme(plot.subtitle = element_text(size = BASE, hjust = 0.5),
        axis.title = element_text(size = BASE - 0.5), axis.text = element_text(size = BASE - 1, color = "black"),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.15))
ggsave(file.path(DD, "depth_learning_curve.png"), g, width = 4.6, height = 3.0, dpi = 400,
       bg = "white", device = ragg::agg_png)
cat("saved depth_learning_curve.png\n")
print(s)
