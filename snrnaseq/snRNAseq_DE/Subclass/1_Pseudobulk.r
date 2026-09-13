library(Seurat) 
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)   

setwd("scz_celltype_paper/snrnaseq/snRNAseq_DE")

################################################################################
# McLean
################################################################################

Ruz <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") %>%
  subset(Age < 70)

# Split Ruz into McLean
Ruz_McLean <- subset(Ruz, subset = Cohort == "McLean")

# Standardize cohort metadata
Ruz_McLean@meta.data$Cohort <- "McLean"

Ruz_McLean@meta.data$Subclass <- gsub("_[0-9].*$", "", Ruz_McLean@meta.data$predicted.id)
Ruz_McLean@meta.data$Subclass <- gsub("^Lamp5_Lhx6$", "Lamp5Lhx6", Ruz_McLean@meta.data$Subclass)

unique(Ruz_McLean@meta.data$Subclass)

bulk <- AggregateExpression(
  Ruz_McLean,
  group.by = c("Subclass", "Donor"),
  return.seurat = TRUE
)

df <- GetAssayData(bulk, assay = "RNA", layer = "counts")

saveRDS(df, "Files/McLean_pseudobulk_SCZ_subclass.rds")

counts_Ruz_McLean <- as.data.frame.matrix(
  table(Ruz_McLean$Donor, Ruz_McLean$Subclass)
)

counts_Ruz_McLean$Donor <- rownames(counts_Ruz_McLean)

meta_Ruz_McLean <- Ruz_McLean@meta.data %>%
  group_by(Donor) %>%
  summarise(
    Cohort = dplyr::first(Cohort),
    Age = dplyr::first(Age),
    Sex = dplyr::first(Sex),
    Diagnosis = dplyr::first(Diagnosis),
    PMI = dplyr::first(PMI),
    .groups = "drop"
  )

meta <- left_join(
  counts_Ruz_McLean,
  meta_Ruz_McLean,
  by = "Donor"
)

write.csv(meta, "Files/Pseudobulk_metadata_subclass_McLean.csv", row.names = FALSE)


################################################################################
# MSSM 1
################################################################################

Ruz <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") %>%
  subset(Age < 70)

# Split Ruz into MtSinai
Ruz_MtSinai <- subset(Ruz, subset = Cohort == "MtSinai")

# Standardize cohort metadata
Ruz_MtSinai@meta.data$Cohort <- "MSSM 1"

Ruz_MtSinai@meta.data$Subclass <- gsub("_[0-9].*$", "", Ruz_MtSinai@meta.data$predicted.id)
Ruz_MtSinai@meta.data$Subclass <- gsub("^Lamp5_Lhx6$", "Lamp5Lhx6", Ruz_MtSinai@meta.data$Subclass)

unique(Ruz_MtSinai@meta.data$Subclass)

bulk <- AggregateExpression(
  Ruz_MtSinai,
  group.by = c("Subclass", "Donor"),
  return.seurat = TRUE
)

df <- GetAssayData(bulk, assay = "RNA", layer = "counts")

saveRDS(df, "Files/MSSM1_pseudobulk_SCZ_subclass.rds")

counts_Ruz_MtSinai <- as.data.frame.matrix(
  table(Ruz_MtSinai$Donor, Ruz_MtSinai$Subclass)
)

counts_Ruz_MtSinai$Donor <- rownames(counts_Ruz_MtSinai)

meta_Ruz_MtSinai <- Ruz_MtSinai@meta.data %>%
  group_by(Donor) %>%
  summarise(
    Cohort = dplyr::first(Cohort),
    Age = dplyr::first(Age),
    Sex = dplyr::first(Sex),
    Diagnosis = dplyr::first(Diagnosis),
    PMI = dplyr::first(PMI),
    .groups = "drop"
  )

meta <- left_join(
  counts_Ruz_MtSinai,
  meta_Ruz_MtSinai,
  by = "Donor"
)

write.csv(meta, "Files/Pseudobulk_metadata_subclass_MSSM1.csv", row.names = FALSE)


################################################################################
# Fröhlich
################################################################################

OFC <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/OFC_updated.rds") %>%
  subset(Age < 70)

OFC$Age <- factor(OFC$Age, levels = sort(unique(OFC$Age)))
OFC$Age <- as.character(OFC$Age)
OFC <- subset(OFC, subset = Age < 70)
OFC$Donor <- droplevels(OFC$Donor)

# Standardize cohort metadata
OFC@meta.data$Cohort <- "Fröhlich"

OFC@meta.data$Subclass <- gsub("_[0-9].*$", "", OFC@meta.data$predicted.id)
OFC@meta.data$Subclass <- gsub("^Lamp5_Lhx6$", "Lamp5Lhx6", OFC@meta.data$Subclass)

unique(OFC@meta.data$Subclass)

bulk <- AggregateExpression(
  OFC,
  group.by = c("Subclass", "Donor"),
  return.seurat = TRUE
)

df <- GetAssayData(bulk, assay = "RNA", layer = "counts")

saveRDS(df, "Files/Frohlich_pseudobulk_SCZ_subclass.rds")

counts_OFC <- as.data.frame.matrix(
  table(OFC$Donor, OFC$Subclass)
)

counts_OFC$Donor <- rownames(counts_OFC)

meta_OFC <- OFC@meta.data %>%
  group_by(Donor) %>%
  summarise(
    Cohort = dplyr::first(Cohort),
    Age = dplyr::first(Age),
    Sex = dplyr::first(Sex),
    Diagnosis = dplyr::first(Diagnosis),
    PMI = dplyr::first(PMI),
    .groups = "drop"
  )

meta <- left_join(
  counts_OFC,
  meta_OFC,
  by = "Donor"
)

write.csv(meta, "Files/Pseudobulk_metadata_subclass_Frohlich.csv", row.names = FALSE)


################################################################################
# Batiuk
################################################################################

Bat <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds") %>%
  subset(Age < 70)

# Standardize cohort metadata
Bat@meta.data$Cohort <- "Batiuk"

Bat@meta.data$Subclass <- gsub("_[0-9].*$", "", Bat@meta.data$predicted.id)
Bat@meta.data$Subclass <- gsub("^Lamp5_Lhx6$", "Lamp5Lhx6", Bat@meta.data$Subclass)

unique(Bat@meta.data$Subclass)

bulk <- AggregateExpression(
  Bat,
  group.by = c("Subclass", "Donor"),
  return.seurat = TRUE
)

df <- GetAssayData(bulk, assay = "RNA", layer = "counts")

saveRDS(df, "Files/Batiuk_pseudobulk_SCZ_subclass.rds")

counts_Bat <- as.data.frame.matrix(
  table(Bat$Donor, Bat$Subclass)
)

counts_Bat$Donor <- rownames(counts_Bat)

meta_Bat <- Bat@meta.data %>%
  group_by(Donor) %>%
  summarise(
    Cohort = dplyr::first(Cohort),
    Age = dplyr::first(Age),
    Sex = dplyr::first(Sex),
    Diagnosis = dplyr::first(Diagnosis),
    PMI = dplyr::first(PMI),
    .groups = "drop"
  )

meta <- left_join(
  counts_Bat,
  meta_Bat,
  by = "Donor"
)

write.csv(meta,"Files/Pseudobulk_metadata_subclass_Batiuk.csv",row.names = FALSE)


################################################################################
# Multiome
################################################################################

Multi <- readRDS("/scratch/nendresz/P1_Brain_scope/Files/Multiome.rds") %>%
  subset(Age_death < 70)

Multi <- subset(
  Multi,
  subset = Disorder %in% c("control", "Schizophrenia")
)

Multi$Age <- Multi$Age_death
Multi$Sex <- Multi$Biological_Sex
Multi$Diagnosis <- Multi$Disorder
Multi$cohort <- "Multiome"

# Standardize cohort metadata
Multi$Cohort <- "Multiome"

Multi$Diagnosis[Multi$Diagnosis == "control"] <- "Control"
Multi$Sex[Multi$Sex == "male"] <- "Male"
Multi$Sex[Multi$Sex == "female"] <- "Female"

Multi@meta.data$Subclass <- gsub("_[0-9].*$", "", Multi@meta.data$predicted.id)
Multi@meta.data$Subclass <- gsub("^Lamp5_Lhx6$", "Lamp5Lhx6", Multi@meta.data$Subclass)

bulk <- AggregateExpression(
  Multi,
  group.by = c("Subclass", "Donor"),
  return.seurat = TRUE
)

df <- GetAssayData(bulk, assay = "RNA", layer = "counts")

saveRDS(df, "Files/Multiome_pseudobulk_SCZ_subclass.rds")

counts_Multi <- as.data.frame.matrix(
  table(Multi$Donor, Multi$Subclass)
)

counts_Multi$Donor <- rownames(counts_Multi)

meta_Multi <- Multi@meta.data %>%
  group_by(Donor) %>%
  summarise(
    Cohort = dplyr::first(Cohort),
    Age = dplyr::first(Age),
    Sex = dplyr::first(Sex),
    Diagnosis = dplyr::first(Diagnosis),
    PMI = dplyr::first(PMI),
    .groups = "drop"
  )

meta <- left_join(
  counts_Multi,
  meta_Multi,
  by = "Donor"
)

write.csv(meta,"Files/Pseudobulk_metadata_subclass_Multiome.csv",row.names = FALSE)


################################################################################
# MSSM 2
################################################################################

MSSM1 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_1.rds") %>%
  subset(Age < 70)

MSSM2 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_2.rds") %>%
  subset(Age < 70)

MSSM3 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_3.rds") %>%
  subset(Age < 70)

MSSM4 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_4.rds") %>%
  subset(Age < 70)

# Standardize cohort metadata
MSSM1@meta.data$Cohort <- "MSSM 2"
MSSM2@meta.data$Cohort <- "MSSM 2"
MSSM3@meta.data$Cohort <- "MSSM 2"
MSSM4@meta.data$Cohort <- "MSSM 2"

MSSM1@meta.data$Subclass <- gsub("_[0-9].*$", "", MSSM1@meta.data$predicted.id)
MSSM2@meta.data$Subclass <- gsub("_[0-9].*$", "", MSSM2@meta.data$predicted.id)
MSSM3@meta.data$Subclass <- gsub("_[0-9].*$", "", MSSM3@meta.data$predicted.id)
MSSM4@meta.data$Subclass <- gsub("_[0-9].*$", "", MSSM4@meta.data$predicted.id)

MSSM1@meta.data$Subclass <- gsub("^Lamp5_Lhx6$", "Lamp5Lhx6", MSSM1@meta.data$Subclass)
MSSM2@meta.data$Subclass <- gsub("^Lamp5_Lhx6$", "Lamp5Lhx6", MSSM2@meta.data$Subclass)
MSSM3@meta.data$Subclass <- gsub("^Lamp5_Lhx6$", "Lamp5Lhx6", MSSM3@meta.data$Subclass)
MSSM4@meta.data$Subclass <- gsub("^Lamp5_Lhx6$", "Lamp5Lhx6", MSSM4@meta.data$Subclass)

bulk1 <- AggregateExpression(
  MSSM1,
  group.by = c("Subclass", "Donor"),
  return.seurat = TRUE
)

bulk2 <- AggregateExpression(
  MSSM2,
  group.by = c("Subclass", "Donor"),
  return.seurat = TRUE
)

bulk3 <- AggregateExpression(
  MSSM3,
  group.by = c("Subclass", "Donor"),
  return.seurat = TRUE
)

bulk4 <- AggregateExpression(
  MSSM4,
  group.by = c("Subclass", "Donor"),
  return.seurat = TRUE
)

# Extract raw counts from each bulked dataset
df1 <- GetAssayData(bulk1, assay = "RNA", layer = "counts")
df2 <- GetAssayData(bulk2, assay = "RNA", layer = "counts")
df3 <- GetAssayData(bulk3, assay = "RNA", layer = "counts")
df4 <- GetAssayData(bulk4, assay = "RNA", layer = "counts")

# Check dupes
intersect(colnames(df1), colnames(df2))
intersect(colnames(df1), colnames(df3))
intersect(colnames(df1), colnames(df4))
intersect(colnames(df2), colnames(df3))
intersect(colnames(df2), colnames(df4))
intersect(colnames(df3), colnames(df4))

# No dupes
df_combined <- cbind(df1, df2, df3, df4)

# Genes from Ensembl
genes <- rownames(df_combined)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = genes,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Replace NAs with original Ensembl IDs
na_genes <- is.na(gene_symbols)

if (any(na_genes)) {
  gene_symbols[na_genes] <- genes[na_genes]
}

# Update the object
rownames(df_combined) <- gene_symbols[rownames(df_combined)]
df_combined <- df_combined[!duplicated(rownames(df_combined)), ]

saveRDS(
  df_combined,
  "Files/MSSM2_pseudobulk_SCZ_subclass.rds"
)

MSSM <- rbind(
  MSSM1@meta.data,
  MSSM2@meta.data,
  MSSM3@meta.data,
  MSSM4@meta.data
)

counts_MSSM <- as.data.frame.matrix(
  table(MSSM$Donor, MSSM$Subclass)
)

counts_MSSM$Donor <- rownames(counts_MSSM)

meta_MSSM <- MSSM %>%
  group_by(Donor) %>%
  summarise(
    Cohort = dplyr::first(Cohort),
    Age = dplyr::first(Age),
    Sex = dplyr::first(Sex),
    Diagnosis = dplyr::first(Diagnosis),
    PMI = dplyr::first(PMI),
    .groups = "drop"
  )

meta <- left_join(
  counts_MSSM,
  meta_MSSM,
  by = "Donor"
)

write.csv(meta,"Files/Pseudobulk_metadata_subclass_MSSM2.csv",row.names = FALSE)


################################################################################
# HBCC
################################################################################

HBCC1 <- readRDS("/scratch/nendresz/PsychAD/Data/HBCC_1.rds")
HBCC2 <- readRDS("/scratch/nendresz/PsychAD/Data/HBCC_2.rds")

HBCC1$Age <- as.numeric(HBCC1$Age)
HBCC2$Age <- as.numeric(HBCC2$Age)

HBCC1 <- subset(HBCC1, subset = Age > 20 & Age < 70)
HBCC2 <- subset(HBCC2, subset = Age > 20 & Age < 70)

# Balance HBCC: take out below 20, youngest with SCZ is 21

# Standardize cohort metadata
HBCC1@meta.data$Cohort <- "HBCC"
HBCC2@meta.data$Cohort <- "HBCC"

HBCC1@meta.data$Subclass <- gsub("_[0-9].*$", "", HBCC1@meta.data$predicted.id)
HBCC2@meta.data$Subclass <- gsub("_[0-9].*$", "", HBCC2@meta.data$predicted.id)

HBCC1@meta.data$Subclass <- gsub("^Lamp5_Lhx6$", "Lamp5Lhx6", HBCC1@meta.data$Subclass)
HBCC2@meta.data$Subclass <- gsub("^Lamp5_Lhx6$", "Lamp5Lhx6", HBCC2@meta.data$Subclass)

bulk1 <- AggregateExpression(
  HBCC1,
  group.by = c("Subclass", "Donor"),
  return.seurat = TRUE
)

bulk2 <- AggregateExpression(
  HBCC2,
  group.by = c("Subclass", "Donor"),
  return.seurat = TRUE
)

# Extract raw counts from each bulked dataset
df1 <- GetAssayData(bulk1, assay = "RNA", layer = "counts")
df2 <- GetAssayData(bulk2, assay = "RNA", layer = "counts")

# Check dupes
intersect(colnames(df1), colnames(df2))

# No dupes
df_combined <- cbind(df1, df2)

# Genes from Ensembl
genes <- rownames(df_combined)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = genes,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Replace NAs with original Ensembl IDs
na_genes <- is.na(gene_symbols)

if (any(na_genes)) {
  gene_symbols[na_genes] <- genes[na_genes]
}

# Update the object
rownames(df_combined) <- gene_symbols[rownames(df_combined)]
df_combined <- df_combined[!duplicated(rownames(df_combined)), ]

saveRDS(df_combined,"Files/HBCC_pseudobulk_SCZ_subclass.rds")

HBCC <- rbind(
  HBCC1@meta.data,
  HBCC2@meta.data
)

counts_HBCC <- as.data.frame.matrix(
  table(HBCC$Donor, HBCC$Subclass)
)

counts_HBCC$Donor <- rownames(counts_HBCC)

meta_HBCC <- HBCC %>%
  group_by(Donor) %>%
  summarise(
    Cohort = dplyr::first(Cohort),
    Age = dplyr::first(Age),
    Sex = dplyr::first(Sex),
    Diagnosis = dplyr::first(Diagnosis),
    PMI = dplyr::first(PMI),
    .groups = "drop"
  )

meta <- left_join(
  counts_HBCC,
  meta_HBCC,
  by = "Donor"
)

write.csv(meta,"Files/Pseudobulk_metadata_subclass_HBCC.csv",row.names = FALSE)










