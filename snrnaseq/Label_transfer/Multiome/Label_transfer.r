setwd("P1_Brain_scope")

library(readr)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

# Load reference
Counts_ref <- readRDS("/project/rrg-shreejoy/nendresz/raw_counts_ref.rds")
meta_hodge <- readRDS("/project/rrg-shreejoy/nendresz/Neurotypical_ref_metadata.rds")

counts_hodge <- Matrix(as.matrix(Counts_ref), sparse = TRUE)
rownames(counts_hodge) <- rownames(Counts_ref)
colnames(counts_hodge) <- colnames(Counts_ref)

# Load Multiome counts
Counts_sn <- readRDS("/scratch/nendresz/P1_Brain_scope/Files/Multiome_matrix.rds")

rownames(Counts_sn) <- Counts_sn$featurekey
Counts_sn$featurekey <- NULL

counts_sn <- Matrix(as.matrix(Counts_sn), sparse = TRUE)
rownames(counts_sn) <- rownames(Counts_sn)
colnames(counts_sn) <- colnames(Counts_sn)

# Match genes
common_genes <- intersect(colnames(counts_hodge), rownames(counts_sn))

counts_hodge_sub <- counts_hodge[, common_genes]
counts_sn_sub <- counts_sn[common_genes, ]

# Create Seurat objects
Seu_hodge <- CreateSeuratObject(
  counts = t(counts_hodge_sub),
  meta.data = meta_hodge
)

Seu_sn <- CreateSeuratObject(
  counts = counts_sn_sub
)

# Normalize + variable features
Seu_hodge <- NormalizeData(Seu_hodge, scale.factor = 1e6)
Seu_hodge <- FindVariableFeatures(Seu_hodge, nfeatures = 3000)

Seu_sn <- NormalizeData(Seu_sn, scale.factor = 1e6)
Seu_sn <- FindVariableFeatures(Seu_sn, nfeatures = 3000)

# Reference PCA
Seu_hodge <- ScaleData(Seu_hodge)
Seu_hodge <- RunPCA(Seu_hodge)
Seu_hodge <- FindNeighbors(Seu_hodge, dims = 1:30)
Seu_hodge <- FindClusters(Seu_hodge)

# Label transfer
anchors <- FindTransferAnchors(
  reference = Seu_hodge,
  query = Seu_sn,
  dims = 1:30,
  reference.reduction = "pca"
)

predictions <- TransferData(
    anchorset = anchors,
  refdata = Seu_hodge$Supertype,
  dims = 1:30
)

Seu_sn <- AddMetaData(Seu_sn, predictions)


# Recreate object with original counts 
Multi <- CreateSeuratObject(counts = counts_sn, meta.data = Seu_sn@meta.data)


# Add metadata

m <- read_tsv(
  "/project/rrg-shreejoy/nendresz/Brain_scope/PEC2_sample_metadata.txt"
)

Multi$Donor <- sub("_[^_]+$", "", Cells(Multi))

meta <- Multi@meta.data
meta$cell_id <- rownames(meta)

meta <- meta %>%
  left_join(m, by = c("Donor" = "Individual_ID"))

rownames(meta) <- meta$cell_id
Multi@meta.data <- meta

# Save final Multiome object
saveRDS(Multi, "/scratch/nendresz/P1_Brain_scope/Files/Multiome.rds")