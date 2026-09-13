# Load packages
library(Seurat)
library(dplyr)

setwd("OFC_cohort")

# Load original OFC
OFC <- readRDS("OFC_seurat_obj_final.rds")

# Combine transferred metadata
query_meta <- do.call(
  rbind,
  lapply(1:5, function(i)
    readRDS(paste0("query_subset_annotated_", i, ".rds"))@meta.data
  )
)

# Add transferred metadata back to original object
OFC <- AddMetaData(OFC, metadata=query_meta)

# Keep SCZ + controls
OFC <- subset(OFC, subset=Classification %in% c("Schizophrenia","Control"))

# Standardize metadata
OFC@meta.data <- OFC@meta.data %>%
  mutate(
    Diagnosis=Classification
  )

# Save
saveRDS(OFC, "/project/rrg-shreejoy/nendresz/Supertypes/OFC_updated.rds")