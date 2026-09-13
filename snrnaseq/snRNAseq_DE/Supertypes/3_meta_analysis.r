# conda activate r_env_meta_analysis
setwd("P1_SCZ_DE_fresh")
library(dplyr)
library(ggplot2)
library(tibble)
library(ggrepel)
library(metafor)
library(writexl)

files <- list.files("/scratch/nendresz/P1_SCZ_DE_fresh/Files",pattern="^meta_results_.*\\.csv$",full.names=TRUE)
cell_types <- basename(files) %>% sub("^meta_results_","",x=.) %>% sub("\\.csv$","",x=.)

subclasses <- c("Astro","Chandelier","Endo","L2_3 IT","L4 IT","L5 ET","L5 IT","L5_6 NP","L6 CT",
"L6 IT","L6 IT Car3","L6b","Lamp5","Lamp5Lhx6","Micro-PVM","Oligo","OPC","Pax6","Pvalb","Sncg",
"Sst","Sst Chodl","Vip","VLMC","log2cells_Sst")
cell_types_supertype <- setdiff(cell_types,subclasses)
cohorts <- c("Bat","OFC","MtSinai","Mclean","MSSM","HBCC","Multi")

de_list <- list()
for(ct in cell_types_supertype){
  de_list[[ct]] <- list()
  for(cohort in cohorts){
    f <- paste0("Files/DE_results_",cohort,"_",ct,".rds")
    if(!file.exists(f)){message("Skipping missing: ",f); next}
    de_list[[ct]][[cohort]] <- readRDS(f)
  }
}

meta_all <- list()
for(ct in cell_types_supertype){
  if(is.null(de_list[[ct]]) || length(de_list[[ct]])==0) next
  message("Meta-analysis for ",ct)

  all_data <- bind_rows(lapply(names(de_list[[ct]]),function(cohort)
    de_list[[ct]][[cohort]] %>% mutate(cohort=cohort,SE=abs(logFC/t))))

  meta_results <- list()
  for(g in unique(all_data$genes)){
    df <- all_data %>% filter(genes==g,!is.na(SE),SE!=0,!is.na(logFC))
    if(nrow(df)>4){
      res <- tryCatch({
        m <- rma(yi=logFC,sei=SE,data=df,method="REML")
        data.frame(cell_type=ct,genes=g,estimate=as.numeric(m$b),se=m$se,pval=m$pval,
                   ci.lb=m$ci.lb,ci.ub=m$ci.ub,k=m$k,tau2=m$tau2,I2=m$I2)
      },error=function(e) NULL)
      if(!is.null(res)) meta_results[[g]] <- res
    }
  }

  final_results <- bind_rows(meta_results)
  if(nrow(final_results)>0) final_results$padj <- p.adjust(final_results$pval,method="fdr")
  meta_all[[ct]] <- final_results
  write.csv(final_results,paste0("Files/meta_results_",ct,".csv"),row.names=FALSE)
}

meta_supertype <- bind_rows(meta_all)

sup_cohorts <- bind_rows(lapply(names(de_list),function(ct)
  bind_rows(lapply(names(de_list[[ct]]),function(cohort)
    de_list[[ct]][[cohort]] %>% mutate(cell_type=ct,cohort=cohort)))))

write.csv(sup_cohorts,"Files/DE_results_cohorts_supertype.csv",row.names=FALSE)
write.csv(meta_supertype,"Files/DE_meta_results_supertype.csv",row.names=FALSE)

sup_bind <- sup_cohorts %>%
  transmute(gene=genes,logFC,PValue=P.Value,FDR=adj.P.Val,cell_type,cohort) %>%
  mutate(cohort=recode(cohort,"Bat"="Batiuk","OFC"="Fröhlich","MtSinai"="MSSM 1",
                       "Mclean"="McLean","MSSM"="MSSM 2","Multi"="Multiome"))

meta_sup_bind <- meta_supertype %>%
  transmute(gene=genes,logFC=estimate,PValue=pval,FDR=padj,cell_type,cohort="Meta-analysis")

de_all_supertype <- bind_rows(sup_bind,meta_sup_bind)

write.csv(de_all_supertype,"Files/DE_all_cohorts_meta_supertype.csv",row.names=FALSE)


de_all_supertype <- read.csv("Files/DE_all_cohorts_meta_supertype.csv")

library(dplyr)
library(writexl)

library(dplyr)
library(writexl)

de_all_supertype <- de_all_supertype |>
  mutate(subclass = sub("-[0-9]+$", "", cell_type))

table(de_all_supertype$subclass)


de_sheets <- split(de_all_supertype, interaction(de_all_supertype$subclass, de_all_supertype$cohort, sep="_", drop=TRUE))
names(de_sheets) <- make.unique(substr(names(de_sheets),1,31))

write_xlsx(de_sheets, "Files/DE_all_cohorts_meta_supertype_by_subclass_cohort.xlsx")
dim(de_all_supertype)
table(de_all_supertype$cohort)
length(unique(de_all_supertype$cell_type))
