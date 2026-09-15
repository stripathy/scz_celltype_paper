
library(ggplot2)
library(dplyr)
library(ggtext) 
setwd("P1_SCZ_paper")


BASE <- 16
AXIS_TITLE <- BASE - 0.5
AXIS_TEXT  <- BASE - 1
GEOM_TEXT  <- BASE * 0.30

final_results <- read.csv("Compositional_analysis/Files/Meta_nonneurons.csv")

colours <- read.csv("Compositional_analysis/Files/cluster_order_and_colors.csv")

final_results$FDR <- final_results$padj

plot_df <- final_results %>%
  left_join(colours, by = c("CellType" = "cluster_label")) %>%
  mutate(signif_label = case_when(FDR < 0.01 ~ "***", FDR < 0.05 ~ "**", FDR < 0.1 ~ "*", FDR >= 0.1 & FDR < 0.2 ~ "+", TRUE ~ ""),
         CellType = factor(CellType, levels = colours$cluster_label))

celltype_labels <- setNames(case_when(
  plot_df$FDR < 0.1 ~ paste0("<span style='color:red;font-size:16pt'><b>", plot_df$CellType, "</b></span>"),
  plot_df$FDR < 0.2 ~ paste0("<span style='color:red'><i>", plot_df$CellType, "</i></span>"),
  TRUE ~ paste0("<span style='color:black'>", plot_df$CellType, "</span>")
), as.character(plot_df$CellType))

p3a <- ggplot(plot_df, aes(x = CellType, y = estimate, fill = cluster_color)) +
  geom_col(color = "black", width = 0.8) +
  geom_errorbar(aes(ymin = estimate - se, ymax = estimate + se), width = 0.5) +
  geom_text(aes(label = signif_label, y = estimate + sign(estimate) * (se + 0.02)), fontface = "bold",
            vjust = ifelse(plot_df$estimate >= 0, 0, 1), size = GEOM_TEXT, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_identity() +
  scale_x_discrete(labels = celltype_labels) +
  annotate("text", x = -Inf, y = 0.58, label = "Increased abundance in SCZ", hjust = -0.1, color = "black", size = GEOM_TEXT) +
  annotate("text", x = -Inf, y = -0.36, label = "Decreased abundance in SCZ", hjust = -0.1, color = "black", size = GEOM_TEXT) +
  theme_classic(base_size = BASE) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1, vjust = 1, size = AXIS_TEXT),
        axis.text.y = element_text(size = AXIS_TEXT), axis.title.y = element_text(size = AXIS_TITLE),
         axis.ticks.y = element_blank(), panel.grid = element_blank(),
        axis.title.x = element_blank(), plot.caption = element_blank()) +
  labs(y = "SCZ abundance change (β ± SE)")

ggsave("Final_figures/Supplemental/Figures/FigureS5_nonneurons.png", p3a, width = 14, height = 6,  dpi = 600)
