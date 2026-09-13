setwd("scz_celltype_paper/snrnaseq/snRNAseq_DE/Supertypes")
library(dplyr)

cohorts <- c("Batiuk","Frohlich","MSSM1","McLean","MSSM2","HBCC","Multiome")

# Get all supertypes
files <- list.files("Files",pattern="^DE_results_.*\\.rds$")
cell_types <- basename(files)
for(cohort in cohorts) cell_types <- sub(paste0("^DE_results_",cohort,"_"),"",cell_types)
cell_types <- unique(sub("\\.rds$","",cell_types))

# Combine cohort results
cohort_results <- list()
for(ct in cell_types){
  for(cohort in cohorts){
    f <- paste0("Files/DE_results_",cohort,"_",ct,".rds")
    if(!file.exists(f)) next
    cohort_results[[paste(cohort,ct,sep="_")]] <- readRDS(f) %>% mutate(cell_type=ct,cohort=cohort)
  }
}

sup_cohorts <- bind_rows(cohort_results)
write.csv(sup_cohorts,"Files/DE_results_cohorts_supertype.csv",row.names=FALSE)

# Combine with meta-analysis
meta_supertype <- read.csv("Files/DE_meta_results_supertype.csv",check.names=FALSE)

sup_bind <- sup_cohorts %>%
  transmute(gene=genes,logFC,PValue=P.Value,FDR=adj.P.Val,cell_type,cohort) %>%
  mutate(cohort=recode(cohort,"Frohlich"="Fröhlich","MSSM1"="MSSM 1","MSSM2"="MSSM 2"))

meta_sup_bind <- meta_supertype %>%
  transmute(gene=genes,logFC=estimate,PValue=pval,FDR=padj,cell_type,cohort="Meta-analysis")

de_all_supertype <- bind_rows(sup_bind,meta_sup_bind)
write.csv(de_all_supertype,"Files/DE_all_cohorts_meta_supertype.csv",row.names=FALSE)