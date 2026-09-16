# Reference taxonomy: SEA-AD (Gabitto et al. 2024) neurotypical snRNA-seq, whose
# Supertype labels are what `refdata` transfers. Variables are named *_seaad for
# that reason; they were called *_hodge until 2026-09-16, which named the wrong
# atlas (Hodge et al. 2019 is a different taxonomy of a different region).
setwd("P1_Brain_scope")

library(readr)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

# Load reference
Counts_ref <- readRDS("/project/rrg-shreejoy/nendresz/raw_counts_ref.rds")
meta_seaad <- readRDS("/project/rrg-shreejoy/nendresz/Neurotypical_ref_metadata.rds")

counts_seaad <- Matrix(as.matrix(Counts_ref), sparse = TRUE)
rownames(counts_seaad) <- rownames(Counts_ref)
colnames(counts_seaad) <- colnames(Counts_ref)

# Load Multiome counts
Counts_sn <- readRDS("/scratch/nendresz/P1_Brain_scope/Files/Multiome_matrix.rds")

rownames(Counts_sn) <- Counts_sn$featurekey
Counts_sn$featurekey <- NULL

counts_sn <- Matrix(as.matrix(Counts_sn), sparse = TRUE)
rownames(counts_sn) <- rownames(Counts_sn)
colnames(counts_sn) <- colnames(Counts_sn)

# Match genes
common_genes <- intersect(colnames(counts_seaad), rownames(counts_sn))

counts_seaad_sub <- counts_seaad[, common_genes]
counts_sn_sub <- counts_sn[common_genes, ]

# Create Seurat objects
Seu_seaad <- CreateSeuratObject(
  counts = t(counts_seaad_sub),
  meta.data = meta_seaad
)

Seu_sn <- CreateSeuratObject(
  counts = counts_sn_sub
)

# Normalize + variable features
Seu_seaad <- NormalizeData(Seu_seaad, scale.factor = 1e6)
Seu_seaad <- FindVariableFeatures(Seu_seaad, nfeatures = 3000)

Seu_sn <- NormalizeData(Seu_sn, scale.factor = 1e6)
Seu_sn <- FindVariableFeatures(Seu_sn, nfeatures = 3000)

# Reference PCA
Seu_seaad <- ScaleData(Seu_seaad)
Seu_seaad <- RunPCA(Seu_seaad)
Seu_seaad <- FindNeighbors(Seu_seaad, dims = 1:30)
Seu_seaad <- FindClusters(Seu_seaad)

# Label transfer
anchors <- FindTransferAnchors(
  reference = Seu_seaad,
  query = Seu_sn,
  dims = 1:30,
  reference.reduction = "pca"
)

predictions <- TransferData(
    anchorset = anchors,
  refdata = Seu_seaad$Supertype,
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