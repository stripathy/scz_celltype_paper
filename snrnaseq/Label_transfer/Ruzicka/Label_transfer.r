# Reference taxonomy: SEA-AD (Gabitto et al. 2024) neurotypical snRNA-seq, whose
# Supertype labels are what `refdata` transfers. Variables are named *_seaad for
# that reason; they were called *_hodge until 2026-09-16, which named the wrong
# atlas (Hodge et al. 2019 is a different taxonomy of a different region).

# Load packages
library(Seurat)
library(dplyr)
library(Matrix)

setwd("PsychENCODE_cohort")

# Load reference
counts_seaad <- Matrix(as.matrix(readRDS("raw_counts_ref.rds")), sparse=TRUE)
ref_anno <- readRDS("Neurotypical_ref_gene_symbols.rds")
meta_seaad <- ref_anno@meta.data
meta_seaad <- meta_seaad[match(rownames(counts_seaad), rownames(meta_seaad)),]

# Load unfiltered PsychENCODE counts + metadata
old_counts <- readRDS("data/counts_sn.RData")
actionnet_summary <- readRDS("data/ACTIONet_summary.rds")

counts_sn <- old_counts
meta_sn <- actionnet_summary$metadata
meta_sn <- meta_sn[match(rownames(counts_sn), rownames(meta_sn)),]

# Keep common genes for label transfer
common_genes <- intersect(colnames(counts_seaad), colnames(counts_sn))
counts_seaad <- counts_seaad[,common_genes]
counts_sn <- counts_sn[,common_genes]

# Create reference + query
ref <- CreateSeuratObject(counts=t(counts_seaad), meta.data=meta_seaad)
query <- CreateSeuratObject(counts=t(counts_sn), meta.data=meta_sn)

# Normalize
ref <- NormalizeData(ref, scale.factor=1e6)
ref <- FindVariableFeatures(ref, nfeatures=3000)

query <- NormalizeData(query, scale.factor=1e6)
query <- FindVariableFeatures(query, nfeatures=3000)

# Reference PCA
ref <- ScaleData(ref)
ref <- RunPCA(ref)

# Transfer labels
anchors <- FindTransferAnchors(
  reference=ref,
  query=query,
  reference.reduction="pca",
  dims=1:30
)

predictions <- TransferData(
  anchorset=anchors,
  refdata=ref$Supertype,
  dims=1:30
)

query <- AddMetaData(query, predictions)

# Put transferred metadata back onto full counts
new_metadata <- query@meta.data
new_metadata <- new_metadata[match(rownames(old_counts), rownames(new_metadata)),]

Ruz <- CreateSeuratObject(
  counts=t(old_counts),
  meta.data=new_metadata
)

# Standardize metadata
Ruz@meta.data <- Ruz@meta.data %>%
  mutate(
    Donor=ID,
    Diagnosis=case_when(
      Phenotype=="SZ" ~ "Schizophrenia",
      Phenotype=="CON" ~ "Control",
      TRUE ~ Phenotype
    ),
    Sex=Gender,
    cohort=Cohort
  )

# Save final object
saveRDS(Ruz, "/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds")