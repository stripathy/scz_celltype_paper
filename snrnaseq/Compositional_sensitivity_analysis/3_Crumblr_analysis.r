# conda activate crumblr_env
# Load libraries
library(SingleCellExperiment)
library(zellkonverter)
library(ggplot2)
library(scattermore)
library(ggtree)
library(crumblr)
library(aplot)
library(tidyverse)
library(dreamlet)
library(kableExtra)
library(ggcorrplot)
library(RColorBrewer)
library(DelayedArray)
library(Seurat)
library(Matrix)
library(dplyr)
library(SummarizedExperiment)
library(reshape2)
library(ggplot2)
library(ape)
library(readr)
library(tibble)


all_meta_counts <- read.csv("P1_Controls/Files/7_cohorts_metadata_NoSST_DEgenes_names.csv", check.names=FALSE)
types <- read.csv("/project/rrg-shreejoy/nendresz/cluster_order_and_colors.csv")
types_neurons <- types %>% filter(class_label %in% c("Neuronal: GABAergic","Neuronal: Glutamatergic"))

run_crumblr <- function(cohort_name,all_meta_counts){
  message("Running ",cohort_name)
  df <- all_meta_counts %>% filter(Cohort==cohort_name) %>% mutate(Sex=tolower(Sex))
  meta_cols <- c("Cohort","Donor","Age","Sex","Diagnosis","PMI")
  counts <- df %>% select(-all_of(meta_cols)) %>% as.data.frame(); rownames(counts) <- df$Donor
  counts <- counts[,colSums(!is.na(counts))>0,drop=FALSE]
  counts <- counts[,colnames(counts) %in% types_neurons$cluster_label,drop=FALSE]
  meta <- df %>% select(all_of(meta_cols)) %>% as.data.frame(); rownames(meta) <- meta$Donor
  meta$Age <- as.numeric(meta$Age); meta$PMI <- as.numeric(meta$PMI)

  dream_formula <- if(cohort_name=="Multiome") ~scale(Age)+Diagnosis+(1|Sex) else ~scale(Age)+scale(PMI)+Diagnosis+(1|Sex)

  cobj <- crumblr(counts)
  fit <- eBayes(dream(cobj,dream_formula,meta))

  topTable(fit,coef="DiagnosisSchizophrenia",number=Inf) %>%
    dplyr::select(logFC,AveExpr,t,P.Value,adj.P.Val) %>%
    tibble::rownames_to_column("CellType") %>%
    mutate(Cohort=cohort_name)
}

cohorts_use <- c("HBCC","MSSM","OFC","Batiuk","McLean","MtSinai","Multiome")
all_results <- bind_rows(lapply(cohorts_use,\(x) run_crumblr(x,all_meta_counts)))
write_csv(all_results,"P1_Controls/Files/crumblr_results_noSSTgenes.csv")



