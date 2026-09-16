# Reference taxonomy: SEA-AD (Gabitto et al. 2024) neurotypical snRNA-seq, whose
# Supertype labels are what `refdata` transfers. Variables are named *_seaad for
# that reason; they were called *_hodge until 2026-09-16, which named the wrong
# atlas (Hodge et al. 2019 is a different taxonomy of a different region).

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

counts_seaad_matrix <- as.matrix(Counts_ref)
counts_seaad <- Matrix(counts_seaad_matrix, sparse=TRUE)

meta_seaad <- readRDS("/project/rrg-shreejoy/nendresz/Neurotypical_ref_metadata.rds")

all(rownames(counts_seaad) %in% rownames(meta_seaad))  # Ensure they align

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
  common_genes <- intersect(colnames(counts_seaad), colnames(counts_sn))
  counts_seaad_f <- counts_seaad[, common_genes]
  counts_sn_f <- counts_sn[, common_genes]

  # SEURAT INTEGRATION
  Seu_seaad_for_int <- CreateSeuratObject(counts = t(counts_seaad_f), 
                                          meta.data = meta_seaad) 
  Seu_sn_for_int <- CreateSeuratObject(counts = t(counts_sn_f), meta.data = meta_sn) 

  Seu.list <- c(Seu_seaad_for_int, Seu_sn_for_int)

  rm(meta_sn, common_genes, counts_seaad_f, counts_sn_f)

  Seu.list <- lapply(X = Seu.list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1000000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  })

  # Extract reference and query datasets
  Seu_seaad <- Seu.list[[1]]  
  Seu_sn <- Seu.list[[2]]     

  # Scale, PCA, clustering for reference
  Seu_seaad <- ScaleData(Seu_seaad)
  Seu_seaad <- RunPCA(Seu_seaad)
  Seu_seaad <- FindNeighbors(Seu_seaad, dims = 1:30)
  Seu_seaad <- FindClusters(Seu_seaad)

  # Transfer
  anchors <- FindTransferAnchors(reference = Seu_seaad, query = Seu_sn, dims = 1:30, 
                                 reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors, refdata = Seu_seaad$Supertype, dims = 1:30)
  Seu_sn <- AddMetaData(Seu_sn, metadata = predictions)

# Extract new metadata
new_meta <- Seu_sn@meta.data

# Re-load original Ensembl-count matrix (unique rownames)
counts_sn_raw <- readRDS(paste0("/scratch/nendresz/PsychAD/Data/Intermediate_dfs/psychAD_matrix_part_", part, ".rds"))
counts_sn_raw <- t(counts_sn_raw)

New_seu <- CreateSeuratObject(counts = t(counts_sn_raw), meta.data = new_meta)

  saveRDS(New_seu, paste0("Data/Seu_objs/LT_psychAD_", part, ".rds"))

  rm(New_seu, Seu_sn, Seu_seaad, Seu.list, predictions, anchors)
}
