#conda activate r_env_meta_analysis
library(dplyr)
library(metafor)
library(ggplot2)
library(readr)

setwd("scz_celltype_paper/snrnaseq")

all_data <- read.csv("Compositional_sensitivity_analysis/Files/crumblr_results_noSSTgenes.csv") %>%
  mutate(SE=abs(logFC/t)) %>%
  filter(!is.na(SE),SE!=0,!is.na(logFC),!is.na(CellType))

meta_results <- list()

for(ct in unique(all_data$CellType)){
  df <- all_data %>% filter(CellType==ct,!is.na(SE),SE!=0)
  if(nrow(df)>=2){
    res <- tryCatch({
      model <- rma(yi=logFC,sei=SE,data=df,method="FE")
      data.frame(CellType=ct,estimate=as.numeric(model$b),se=model$se,zval=model$zval,pval=model$pval,
                 I2=model$I2,ci.lb=model$ci.lb,ci.ub=model$ci.ub,k=model$k)
    },error=function(e) NULL)
    if(!is.null(res)) meta_results[[ct]] <- res
  }
}

final_results <- bind_rows(meta_results) %>% mutate(padj=p.adjust(pval,method="fdr"))

write_csv(final_results,"Compositional_sensitivity_analysis/Files/meta_noSST_genes.csv")