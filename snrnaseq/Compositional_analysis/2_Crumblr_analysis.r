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
setwd("P1_Compositional_analysis")

all_meta_counts <- read.csv("Files/7_cohorts_metadata_names.csv", row.names = 1, check.names = FALSE)

types <- read.csv("/project/rrg-shreejoy/nendresz/cluster_order_and_colors.csv")
unique(types$class_label)
types_neurons <- filter(types, class_label %in% c("Neuronal: GABAergic", "Neuronal: Glutamatergic"))

run_crumblr <- function(cohort_name, all_meta_counts) {
  
  message("Running ", cohort_name)

  df <- all_meta_counts %>%
  filter(Cohort == cohort_name) %>%
  mutate(Sex = tolower(Sex)) 
  
  meta_cols <- c("Cohort", "Donor", "Age", "Sex", "Diagnosis", "PMI")
  
  counts <- df %>%
    select(-all_of(meta_cols)) %>%
    as.data.frame()
  
  rownames(counts) <- df$Donor

  counts <- counts[, colSums(!is.na(counts)) > 0, drop = FALSE]

neuron_labels <- types_neurons$cluster_label

counts <- counts[, colnames(counts) %in% neuron_labels, drop = FALSE]
  
  meta <- df %>%
    select(all_of(meta_cols)) %>%
    as.data.frame()
  
  rownames(meta) <- meta$Donor
  
  meta$Age <- as.numeric(meta$Age)
  meta$PMI <- as.numeric(meta$PMI)


    if (cohort_name == "Multiome") {
    dream_formula <- ~ scale(Age) + Diagnosis + (Sex) 
  } else {
    dream_formula <- ~  scale(Age) + scale(PMI) + Diagnosis + (Sex) 
  }
  
  cobj <- crumblr(counts)
  
  fit <- dream(cobj, dream_formula, meta)
  fit <- eBayes(fit)
  
  res <- topTable(
    fit,
    coef = "DiagnosisSchizophrenia",
    number = Inf
  ) %>%
    dplyr::select(logFC, AveExpr, t, P.Value, adj.P.Val) %>%
    tibble::rownames_to_column("CellType") %>%
    mutate(Cohort = cohort_name)
  
  return(res)}

cohorts_use <- c("HBCC", "MSSM", "OFC", "Batiuk", "McLean", "MtSinai", "Multiome")

all_results <- bind_rows(
  lapply(cohorts_use, function(x) {
    run_crumblr(x, all_meta_counts)}))


write_csv(all_results,"Files/crumblr_results_final_one_doc.csv")

