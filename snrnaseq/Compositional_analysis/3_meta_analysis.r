#conda activate r_env_meta_analysis
setwd("P1_Compositional_analysis")

library(dplyr)
library(metafor)
library(ggplot2)


# Read combined crumblr results
all_data <- read.csv("Files/crumblr_results_final_one_doc.csv")


# Make SE from t-stat
all_data <- all_data %>%
  mutate(
    SE = abs(logFC / t)
  ) %>%
  filter(
    !is.na(SE),
    SE != 0,
    !is.na(logFC),
    !is.na(CellType)
  )

cell_types <- unique(all_data$CellType)

meta_results <- list()

for (ct in cell_types) {
  
  df <- all_data %>%
    filter(CellType == ct, !is.na(SE), SE != 0)
  
  if (nrow(df) >= 2) {
    
    res <- tryCatch({
      
      model <- rma(
        yi = logFC,
        sei = SE,
        data = df,
        method = "FE"
      )
      
      data.frame(
        CellType = ct,
        estimate = as.numeric(model$b),
        se = model$se,
        zval = model$zval,
        pval = model$pval,
        I2 = model$I2,
        ci.lb = model$ci.lb,
        ci.ub = model$ci.ub,
        k = model$k
      )
      
    }, error = function(e) NULL)
    
    if (!is.null(res)) {
      meta_results[[ct]] <- res
    }
  }
}

final_results <- bind_rows(meta_results)

final_results$padj <- p.adjust(final_results$pval, method = "fdr")

write.csv(final_results, "/scratch/nendresz/P1_Compositional_analysis/Files/crumblr_results_final_meta.csv")