
#conda activate r_env_meta_analysis
setwd("P1_SCZ_paper/Compositional_analysis")

library(dplyr)
library(metafor)
library(ggplot2)

# Read combined crumblr results
all_data <- read.csv("Files/crumblr_results_nonneurons.csv")


# Make SE from t-stat
all_data <- all_data %>%
  mutate(SE = abs(logFC / t)) %>%
  filter(!is.na(SE), SE != 0,
    !is.na(logFC),
    !is.na(CellType))

meta_results <- lapply(unique(all_data$CellType),function(ct){
  df <- all_data %>% filter(CellType==ct,!is.na(SE),SE!=0)
  if(nrow(df)<2) return(NULL)
  tryCatch({
    m <- rma(yi=logFC,sei=SE,data=df,method="FE")
    data.frame(CellType=ct,estimate=as.numeric(m$b),se=m$se,zval=m$zval,pval=m$pval,I2=m$I2,ci.lb=m$ci.lb,ci.ub=m$ci.ub,k=m$k)
  },error=function(e) NULL)
})

final_results <- bind_rows(meta_results) %>% mutate(padj=p.adjust(pval,method="fdr"))

write.csv(final_results,"Files/Meta_nonneurons.csv",row.names=FALSE)

