
#Load packages 
library(readxl)
library(readr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(org.Hs.eg.db)

setwd("PsychAD")

# -------------------------
# Load in Reference
# -------------------------
Counts_ref <- readRDS("/project/rrg-shreejoy/nendresz/raw_counts_ref.rds")

counts_hodge_matrix <- as.matrix(Counts_ref)
counts_hodge <- Matrix(counts_hodge_matrix, sparse=TRUE)

meta_hodge <- readRDS("/project/rrg-shreejoy/nendresz/Neurotypical_ref_metadata.rds")

all(rownames(counts_hodge) %in% rownames(meta_hodge))  # Ensure they align

# -------------------------
# Loop through subsets
# -------------------------
for (part in 1:20) {
  message("Processing subset: ", part)

  # Load the testing count matrix
  counts_sn <- readRDS(paste0("/scratch/nendresz/PsychAD/Data/Intermediate_dfs/psychAD_matrix_part_", part, ".rds"))
  counts_sn <- t(counts_sn)
  meta_sn <- read.csv("/scratch/nendresz/PsychAD/Data/psychAD_metadata.csv")
  meta_sn <- as.data.frame(meta_sn)
  rownames(meta_sn) <- meta_sn$barcodekey

  # Subset metadata to match this matrix’s cells
  meta_sn <- meta_sn[rownames(counts_sn), , drop = FALSE]

  # Quick check
  all(rownames(counts_sn) == rownames(meta_sn))

  # Map Ensembl → gene symbol
  ensembl_ids <- colnames(counts_sn)
  symbols <- mapIds(
    org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  colnames(counts_sn) <- symbols

  # Filter counts matrices
  common_genes <- intersect(colnames(counts_hodge), colnames(counts_sn))
  counts_hodge_f <- counts_hodge[, common_genes]
  counts_sn_f <- counts_sn[, common_genes]

  # SEURAT INTEGRATION
  Seu_hodge_for_int <- CreateSeuratObject(counts = t(counts_hodge_f), 
                                          meta.data = meta_hodge) 
  Seu_sn_for_int <- CreateSeuratObject(counts = t(counts_sn_f), meta.data = meta_sn) 

  Seu.list <- c(Seu_hodge_for_int, Seu_sn_for_int)

  rm(meta_sn, common_genes, counts_hodge_f, counts_sn_f)

  Seu.list <- lapply(X = Seu.list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1000000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  })

  # Extract reference and query datasets
  Seu_hodge <- Seu.list[[1]]  
  Seu_sn <- Seu.list[[2]]     

  # Scale, PCA, clustering for reference
  Seu_hodge <- ScaleData(Seu_hodge)
  Seu_hodge <- RunPCA(Seu_hodge)
  Seu_hodge <- FindNeighbors(Seu_hodge, dims = 1:30)
  Seu_hodge <- FindClusters(Seu_hodge)

  # Transfer
  anchors <- FindTransferAnchors(reference = Seu_hodge, query = Seu_sn, dims = 1:30, 
                                 reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors, refdata = Seu_hodge$Supertype, dims = 1:30)
  Seu_sn <- AddMetaData(Seu_sn, metadata = predictions)

# Extract new metadata
new_meta <- Seu_sn@meta.data

# Re-load original Ensembl-count matrix (unique rownames)
counts_sn_raw <- readRDS(paste0("/scratch/nendresz/PsychAD/Data/Intermediate_dfs/psychAD_matrix_part_", part, ".rds"))
counts_sn_raw <- t(counts_sn_raw)

New_seu <- CreateSeuratObject(counts = t(counts_sn_raw), meta.data = new_meta)

  saveRDS(New_seu, paste0("Data/Seu_objs/LT_psychAD_", part, ".rds"))

  rm(New_seu, Seu_sn, Seu_hodge, Seu.list, predictions, anchors)
}
