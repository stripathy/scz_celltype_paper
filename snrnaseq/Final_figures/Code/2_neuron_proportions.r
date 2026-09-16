#Load library
library(dplyr)

#Set working directory 
setwd("scz_celltype_paper/snrnaseq/")

#Load data 
all_meta_counts <- read.csv("Compositional_analysis/Files/7_cohorts_metadata_names.csv", row.names=1, check.names=FALSE)
types <- read.csv("Compositional_analysis/Files/cluster_order_and_colors.csv")

#Get neuron types
types_neurons <- filter(types, class_label %in% c("Neuronal: GABAergic","Neuronal: Glutamatergic"))

neuron_labels <- types_neurons$cluster_label

#Calculate within neuron proportions
neuron_props <- all_meta_counts |>
  select(Cohort, Donor, Age, Sex, Diagnosis, PMI, all_of(neuron_labels)) |>
  mutate(total_neurons=rowSums(across(all_of(neuron_labels)), na.rm=TRUE)) |>
  mutate(across(all_of(neuron_labels), ~ .x/total_neurons)) |>
  select(-total_neurons)

#Save
  write.csv(neuron_props, "Final_figures/Data/neuron_props_7_cohorts.csv", row.names = FALSE)