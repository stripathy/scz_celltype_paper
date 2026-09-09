# Load packages
library(Seurat)
library(Matrix)

setwd("OFC_cohort")

# Load reference
counts_hodge <- Matrix(as.matrix(readRDS("/project/s/shreejoy/nendresz/raw_counts_ref.rds")), sparse=TRUE)
meta_hodge <- readRDS("/project/s/shreejoy/nendresz/Neurotypical_ref_metadata.rds")

ref <- CreateSeuratObject(counts=t(counts_hodge), meta.data=meta_hodge)
ref <- NormalizeData(ref)
ref <- FindVariableFeatures(ref)
ref <- ScaleData(ref)
ref <- RunPCA(ref, npcs=20)

# Load and split OFC
OFC <- readRDS("OFC_seurat_obj_final.rds")

num_parts <- 5
splits <- split(colnames(OFC), cut(seq_along(colnames(OFC)), num_parts, labels=FALSE))

for(i in 1:num_parts){
  
  message("Processing subset ", i)
  
  query <- subset(OFC, cells=splits[[i]])
  query <- NormalizeData(query)
  query <- FindVariableFeatures(query)
  query <- ScaleData(query)
  query <- RunPCA(query, npcs=20)
  
  # Match genes
  common_genes <- intersect(rownames(ref), rownames(query))
  ref_sub <- subset(ref, features=common_genes)
  query_sub <- subset(query, features=common_genes)
  
  # Label transfer
  anchors <- FindTransferAnchors(reference=ref_sub, query=query_sub, dims=1:20)
  predictions <- TransferData(anchorset=anchors, refdata=ref_sub$Supertype, dims=1:20)
  
  query <- AddMetaData(query, metadata=predictions)
  
  # Save annotated subset
  saveRDS(query, paste0("query_subset_annotated_", i, ".rds"))
  
  rm(query, query_sub, ref_sub, anchors, predictions)
  gc()
}