library(Seurat)
library(dplyr)
library(Matrix)

# SST DE GENES 

meta <- read.csv("P1_SCZ_DE_fresh/Files/DE_genes_all_cells_scz.csv")
all_genes <- meta %>% filter(cell_type=="Sst",padj<0.1) %>% pull(genes) %>% unique()
length(all_genes)

# LOAD COHORTS 

HBCC1 <- readRDS("/scratch/nendresz/PsychAD/Data/HBCC_1_symbols.rds")
HBCC2 <- readRDS("/scratch/nendresz/PsychAD/Data/HBCC_2_symbols.rds")

HBCC1$Age <- as.numeric(HBCC1$Age)
HBCC2$Age <- as.numeric(HBCC2$Age)

HBCC1 <- subset(HBCC1,Age>20 & Age<70)
HBCC2 <- subset(HBCC2,Age>20 & Age<70)

MSSM1 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_1_symbols.rds") %>% subset(Age<70)
MSSM2 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_2_symbols.rds") %>% subset(Age<70)
MSSM3 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_3_symbols.rds") %>% subset(Age<70)
MSSM4 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_4_symbols.rds") %>% subset(Age<70)

OFC <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/OFC_updated.rds") %>% subset(Age<70)
Bat <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds") %>% subset(Age<70)

Ruz <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") %>% subset(Age<70)
Ruz_McLean <- subset(Ruz,Cohort=="McLean")
Ruz_MtSinai <- subset(Ruz,Cohort=="MtSinai")

Multi <- readRDS("/scratch/nendresz/P1_Brain_scope/Files/Multiome.rds") %>%
  subset(Age_death<70 & Disorder %in% c("control","Schizophrenia"))

rm(Ruz); gc()

cohorts <- list(
  Batiuk=Bat, Frohlich=OFC,
  HBCC1=HBCC1, HBCC2=HBCC2,
  McLean=Ruz_McLean, MSSM1=Ruz_MtSinai,
  MSSM2_part1=MSSM1, MSSM2_part2=MSSM2,
  MSSM2_part3=MSSM3, MSSM2_part4=MSSM4,
  Multiome=Multi
)

rm(Bat,OFC,HBCC1,HBCC2,Ruz_McLean,Ruz_MtSinai,MSSM1,MSSM2,MSSM3,MSSM4,Multi)
gc()

# REFERENCE 

Counts_ref <- readRDS("/project/rrg-shreejoy/nendresz/raw_counts_ref.rds")
meta_hodge <- readRDS("/project/rrg-shreejoy/nendresz/Neurotypical_ref_metadata.rds")
meta_hodge <- meta_hodge[rownames(Counts_ref),,drop=FALSE]

options(future.globals.maxSize=Inf)

# LABEL TRANSFER 

for(nm in names(cohorts)){

  cat("\nPROCESSING:",nm,"\n")

  obj <- cohorts[[nm]]
  counts_sn <- t(GetAssayData(obj,assay="RNA",layer="counts"))
  meta_sn <- as.data.frame(obj@meta.data)

  stopifnot(identical(rownames(counts_sn),rownames(meta_sn)))

  old_predictions <- if("predicted.id" %in% colnames(meta_sn))
    meta_sn$predicted.id else rep(NA,nrow(meta_sn))

  # Remove SST DE genes
  common_genes <- intersect(colnames(Counts_ref),colnames(counts_sn))
  genes_use <- setdiff(common_genes,all_genes)

  cat(
    "Cells:",nrow(counts_sn),
    "| Common genes:",length(common_genes),
    "| Removed:",sum(common_genes %in% all_genes),
    "| Used:",length(genes_use),"\n"
  )

  counts_sn <- as(counts_sn[,genes_use,drop=FALSE],"dgCMatrix")
  counts_hodge <- as(Counts_ref[,genes_use,drop=FALSE],"dgCMatrix")

  # Seurat objects
  Seu_hodge <- CreateSeuratObject(t(counts_hodge),meta.data=meta_hodge)
  Seu_sn <- CreateSeuratObject(t(counts_sn),meta.data=meta_sn)

  rm(counts_sn,counts_hodge); gc()

  # Normalize
  Seu_hodge <- NormalizeData(Seu_hodge,scale.factor=1e6,verbose=FALSE) %>%
    FindVariableFeatures(nfeatures=min(3000,nrow(Seu_hodge)),verbose=FALSE) %>%
    ScaleData(verbose=FALSE)

  Seu_sn <- NormalizeData(Seu_sn,scale.factor=1e6,verbose=FALSE) %>%
    FindVariableFeatures(nfeatures=min(3000,nrow(Seu_sn)),verbose=FALSE)

  # Reference PCA
  npcs_use <- min(30,length(genes_use)-1)

  Seu_hodge <- RunPCA(Seu_hodge,npcs=npcs_use,verbose=FALSE) %>%
    FindNeighbors(dims=1:npcs_use,verbose=FALSE) %>%
    FindClusters(verbose=FALSE)

  # Transfer labels
  anchors <- FindTransferAnchors(
    reference=Seu_hodge,
    query=Seu_sn,
    dims=1:npcs_use,
    reference.reduction="pca"
  )

  predictions <- TransferData(
    anchorset=anchors,
    refdata=Seu_hodge$Supertype,
    dims=1:npcs_use
  )

  Seu_sn <- AddMetaData(Seu_sn,predictions)
  Seu_sn$old_predictions <- old_predictions

  # Save
  saveRDS(Seu_sn,paste0("P1_Controls/Files/",nm,"_NoSST_DEgenes.rds"))
  write.csv(Seu_sn@meta.data,paste0("P1_Controls/Files/",nm,"_NoSST_DEgenes.csv"))

  cat("SAVED:",nm,"\n")

  rm(obj,Seu_hodge,Seu_sn,anchors,predictions,old_predictions,
     genes_use,common_genes)
  gc()
}

# COMBINE SPLIT COHORTS 

HBCC_meta <- bind_rows(
  read.csv("P1_Controls/Files/HBCC1_NoSST_DEgenes.csv",row.names=1,check.names=FALSE),
  read.csv("P1_Controls/Files/HBCC2_NoSST_DEgenes.csv",row.names=1,check.names=FALSE))
write.csv(HBCC_meta,"P1_Controls/Files/HBCC_NoSST_DEgenes.csv")

MSSM2_meta <- bind_rows(
  read.csv("P1_Controls/Files/MSSM2_part1_NoSST_DEgenes.csv",row.names=1,check.names=FALSE),
  read.csv("P1_Controls/Files/MSSM2_part2_NoSST_DEgenes.csv",row.names=1,check.names=FALSE),
  read.csv("P1_Controls/Files/MSSM2_part3_NoSST_DEgenes.csv",row.names=1,check.names=FALSE),
  read.csv("P1_Controls/Files/MSSM2_part4_NoSST_DEgenes.csv",row.names=1,check.names=FALSE))
write.csv(MSSM2_meta,"P1_Controls/Files/MSSM2_NoSST_DEgenes.csv")

