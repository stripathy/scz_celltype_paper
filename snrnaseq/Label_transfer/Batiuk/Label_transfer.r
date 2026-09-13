# Load packages
library(Seurat)
library(Matrix)

setwd("Supplemental_Analysis_Dan_paper_011725")

# Load reference
counts_hodge <- Matrix(as.matrix(readRDS("/project/s/shreejoy/nendresz/raw_counts_ref.rds")), sparse=TRUE)
meta_hodge <- readRDS("/project/s/shreejoy/nendresz/Neurotypical_ref_metadata.rds")

# Load unfiltered Batiuk counts + metadata
old_counts <- readRDS("counts_combined_sparse.rds")
counts_sn <- t(old_counts)
meta_sn <- readRDS("Merged_Metadata_Batiuk.rds")
rownames(meta_sn) <- meta_sn$CellBarcode

# Keep common genes for label transfer
common_genes <- intersect(colnames(counts_hodge), colnames(counts_sn))
counts_hodge <- counts_hodge[,common_genes]
counts_sn <- counts_sn[,common_genes]

# Create reference + query
ref <- CreateSeuratObject(counts=t(counts_hodge), meta.data=meta_hodge)
query <- CreateSeuratObject(counts=t(counts_sn), meta.data=meta_sn)

# Normalize + reference PCA
ref <- NormalizeData(ref, scale.factor=1e6)
ref <- FindVariableFeatures(ref, nfeatures=3000)
ref <- ScaleData(ref)
ref <- RunPCA(ref)

query <- NormalizeData(query, scale.factor=1e6)
query <- FindVariableFeatures(query, nfeatures=3000)

# Transfer labels
anchors <- FindTransferAnchors(reference=ref, query=query, reference.reduction="pca", dims=1:30)
predictions <- TransferData(anchorset=anchors, refdata=ref$Supertype, dims=1:30)
query <- AddMetaData(query, predictions)

# Put transferred metadata back onto unfiltered counts
new_metadata <- query@meta.data
new_metadata <- new_metadata[colnames(old_counts),,drop=FALSE]

new_seu <- CreateSeuratObject(
  counts=old_counts,
  meta.data=new_metadata
)

#Update metadata labels for consistency 

Bat <- new_seu 

Bat <- subset(Bat, subset = DonorLabel != "MB8")

Bat@meta.data <- Bat@meta.data %>%
  mutate(
    Donor = DonorLabel,
    Diagnosis = case_when(
      Diagnosis == "Scz" ~ "Schizophrenia",
      Diagnosis == "Ctr" ~ "Control",
      TRUE ~ Diagnosis
    ),
    Sex = case_when(
      Gender == "F" ~ "Female",
      Gender == "M" ~ "Male",
      TRUE ~ Gender
    ),
    AgeOnset = First_schizophrenia_symptoms_age,
    PMI = PMI_hrs,
    cohort = "Batiuk"
  )



saveRDS(Bat, "/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds")