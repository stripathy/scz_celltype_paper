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

setwd("scz_celltype_paper/snrnaseq")


all_meta_counts <- read.csv("Compositional_sensitivity_analysis/Files/7_cohorts_metadata_NoSST_DEgenes_names.csv", check.names=FALSE)
types <- read.csv("Compositional_analysis/Files/cluster_order_and_colors.csv")
types_neurons <- types %>% filter(class_label %in% c("Neuronal: GABAergic","Neuronal: Glutamatergic"))

run_crumblr <- function(cohort_name,all_meta_counts){
  message("Running ",cohort_name)
  df <- all_meta_counts %>% filter(Cohort==cohort_name) %>% mutate(Sex=tolower(Sex))
  meta_cols <- c("Cohort","Donor","Age","Sex","Diagnosis","PMI")
  counts <- df %>% select(-all_of(meta_cols)) %>% as.data.frame(); rownames(counts) <- df$Donor
  # Drop cell types absent from this cohort. The committed matrix marks them NA,
  # but a matrix pivoted from the per-cell tables marks them 0; either way they are
  # structurally empty and must not enter crumblr's CLR. Verified a no-op on
  # 7_cohorts_metadata_names.csv (no all-zero column there) -- 2026-09-13.
  counts <- counts[,colSums(!is.na(counts))>0 & colSums(counts,na.rm=TRUE)>0,drop=FALSE]
  counts <- counts[,colnames(counts) %in% types_neurons$cluster_label,drop=FALSE]
  meta <- df %>% select(all_of(meta_cols)) %>% as.data.frame(); rownames(meta) <- meta$Donor
  meta$Age <- as.numeric(meta$Age); meta$PMI <- as.numeric(meta$PMI)

  dream_formula <- if(cohort_name=="Multiome") ~scale(Age)+Diagnosis+(Sex) else ~scale(Age)+scale(PMI)+Diagnosis+(Sex)

  cobj <- crumblr(counts)
  fit <- eBayes(dream(cobj,dream_formula,meta))

  # The Diagnosis contrast is named after the non-reference factor level, so it is
  # "DiagnosisSchizophrenia" for this metadata but "DiagnosisSCZ" for a table using
  # the repo's Control/SCZ vocabulary. Read it off the fit instead of hard-coding.
  dx_coef <- grep("^Diagnosis", colnames(fit$coefficients), value = TRUE)
  # Exactly one Diagnosis column means a clean two-level contrast. More than one
  # means Diagnosis has >2 levels in this fit — which happens if cohorts are pooled,
  # because Multiome spells its controls "control" and the rest spell them "Control".
  # Stop rather than silently contrast against the wrong reference.
  stopifnot(length(dx_coef) == 1)
  dx_coef <- dx_coef[1]
  topTable(fit,coef=dx_coef,number=Inf) %>%
    dplyr::select(logFC,AveExpr,t,P.Value,adj.P.Val) %>%
    tibble::rownames_to_column("CellType") %>%
    mutate(Cohort=cohort_name)
}

cohorts_use <- c("HBCC","MSSM 2","Fröhlich","Batiuk","McLean","MSSM 1","Multiome")
all_results <- bind_rows(lapply(cohorts_use,\(x) run_crumblr(x,all_meta_counts)))
write_csv(all_results,"Compositional_sensitivity_analysis/Files/crumblr_results_noSSTgenes.csv")



