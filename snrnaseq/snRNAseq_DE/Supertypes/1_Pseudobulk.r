library(Seurat) 
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)   


setwd("scz_celltype_paper/snrnaseq/snRNAseq_DE")

# ==================== McLean ====================

Ruz <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") %>%
  subset(Age < 70)

Ruz_McLean <- subset(Ruz, subset = Cohort == "McLean")
Ruz_McLean@meta.data$Cohort <- "McLean"

bulk <- AggregateExpression(Ruz_McLean, group.by = c("predicted.id", "Donor"), return.seurat = TRUE)
df <- GetAssayData(bulk, assay = "RNA", layer = "counts")

saveRDS(df, "Files/McLean_pseudobulk_SCZ_supertype.rds")

counts_Ruz_McLean <- as.data.frame.matrix(table(Ruz_McLean$Donor, Ruz_McLean$predicted.id))
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

meta <- left_join(counts_Ruz_McLean, meta_Ruz_McLean, by = "Donor")

write.csv(meta, "Files/Pseudobulk_metadata_supertype_McLean.csv", row.names = FALSE)


# ==================== MSSM 1 ====================

Ruz <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") %>%
  subset(Age < 70)

Ruz_MtSinai <- subset(Ruz, subset = Cohort == "MtSinai")
Ruz_MtSinai@meta.data$Cohort <- "MSSM 1"

bulk <- AggregateExpression(Ruz_MtSinai, group.by = c("predicted.id", "Donor"), return.seurat = TRUE)
df <- GetAssayData(bulk, assay = "RNA", layer = "counts")

saveRDS(df, "Files/MSSM1_pseudobulk_SCZ_supertype.rds")

counts_Ruz_MtSinai <- as.data.frame.matrix(table(Ruz_MtSinai$Donor, Ruz_MtSinai$predicted.id))
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

meta <- left_join(counts_Ruz_MtSinai, meta_Ruz_MtSinai, by = "Donor")

write.csv(meta, "Files/Pseudobulk_metadata_supertype_MSSM1.csv", row.names = FALSE)


# ==================== Fröhlich ====================

OFC <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/OFC_updated.rds") %>%
  subset(Age < 70)

OFC$Age <- factor(OFC$Age, levels = sort(unique(OFC$Age)))
OFC$Age <- as.character(OFC$Age)
OFC <- subset(OFC, subset = Age < 70)
OFC$Donor <- droplevels(OFC$Donor)
OFC@meta.data$Cohort <- "Fröhlich"

bulk <- AggregateExpression(OFC, group.by = c("predicted.id", "Donor"), return.seurat = TRUE)
df <- GetAssayData(bulk, assay = "RNA", layer = "counts")

saveRDS(df, "Files/Frohlich_pseudobulk_SCZ_supertype.rds")

counts_OFC <- as.data.frame.matrix(table(OFC$Donor, OFC$predicted.id))
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

meta <- left_join(counts_OFC, meta_OFC, by = "Donor")

write.csv(meta, "Files/Pseudobulk_metadata_supertype_Frohlich.csv", row.names = FALSE)


# ==================== Batiuk ====================

Bat <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds") %>%
  subset(Age < 70)

Bat@meta.data$Cohort <- "Batiuk"

bulk <- AggregateExpression(Bat, group.by = c("predicted.id", "Donor"), return.seurat = TRUE)
df <- GetAssayData(bulk, assay = "RNA", layer = "counts")

saveRDS(df, "Files/Batiuk_pseudobulk_SCZ_supertype.rds")

counts_Bat <- as.data.frame.matrix(table(Bat$Donor, Bat$predicted.id))
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

meta <- left_join(counts_Bat, meta_Bat, by = "Donor")

write.csv(meta, "Files/Pseudobulk_metadata_supertype_Batiuk.csv", row.names = FALSE)


# ==================== Multiome ====================

Multi <- readRDS("/scratch/nendresz/P1_Brain_scope/Files/Multiome.rds") %>%
  subset(Age_death < 70)

Multi <- subset(Multi, subset = Disorder %in% c("control", "Schizophrenia"))

Multi$Age <- Multi$Age_death
Multi$Sex <- Multi$Biological_Sex
Multi$Diagnosis <- Multi$Disorder
Multi$Cohort <- "Multiome"

Multi$Diagnosis[Multi$Diagnosis == "control"] <- "Control"
Multi$Sex[Multi$Sex == "male"] <- "Male"
Multi$Sex[Multi$Sex == "female"] <- "Female"

bulk <- AggregateExpression(Multi, group.by = c("predicted.id", "Donor"), return.seurat = TRUE)
df <- GetAssayData(bulk, assay = "RNA", layer = "counts")

saveRDS(df, "Files/Multiome_pseudobulk_SCZ_supertype.rds")

counts_Multi <- as.data.frame.matrix(table(Multi$Donor, Multi$predicted.id))
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

meta <- left_join(counts_Multi, meta_Multi, by = "Donor")

write.csv(meta, "Files/Pseudobulk_metadata_supertype_Multiome.csv", row.names = FALSE)


# ==================== MSSM 2 ====================

MSSM1 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_1.rds") %>% subset(Age < 70)
MSSM2 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_2.rds") %>% subset(Age < 70)
MSSM3 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_3.rds") %>% subset(Age < 70)
MSSM4 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_4.rds") %>% subset(Age < 70)

MSSM1@meta.data$Cohort <- "MSSM 2"
MSSM2@meta.data$Cohort <- "MSSM 2"
MSSM3@meta.data$Cohort <- "MSSM 2"
MSSM4@meta.data$Cohort <- "MSSM 2"

bulk1 <- AggregateExpression(MSSM1, group.by = c("predicted.id", "Donor"), return.seurat = TRUE)
bulk2 <- AggregateExpression(MSSM2, group.by = c("predicted.id", "Donor"), return.seurat = TRUE)
bulk3 <- AggregateExpression(MSSM3, group.by = c("predicted.id", "Donor"), return.seurat = TRUE)
bulk4 <- AggregateExpression(MSSM4, group.by = c("predicted.id", "Donor"), return.seurat = TRUE)

df1 <- GetAssayData(bulk1, assay = "RNA", layer = "counts")
df2 <- GetAssayData(bulk2, assay = "RNA", layer = "counts")
df3 <- GetAssayData(bulk3, assay = "RNA", layer = "counts")
df4 <- GetAssayData(bulk4, assay = "RNA", layer = "counts")

# Check duplicate columns
intersect(colnames(df1), colnames(df2))
intersect(colnames(df1), colnames(df3))
intersect(colnames(df1), colnames(df4))
intersect(colnames(df2), colnames(df3))
intersect(colnames(df2), colnames(df4))
intersect(colnames(df3), colnames(df4))

df_combined <- cbind(df1, df2, df3, df4)

# Convert Ensembl IDs to gene symbols
genes <- rownames(df_combined)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = genes,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

na_genes <- is.na(gene_symbols)

if (any(na_genes)) {
  gene_symbols[na_genes] <- genes[na_genes]
}

rownames(df_combined) <- gene_symbols[rownames(df_combined)]
df_combined <- df_combined[!duplicated(rownames(df_combined)), ]

saveRDS(df_combined, "Files/MSSM2_pseudobulk_SCZ_supertype.rds")

MSSM <- rbind(
  MSSM1@meta.data,
  MSSM2@meta.data,
  MSSM3@meta.data,
  MSSM4@meta.data
)

counts_MSSM <- as.data.frame.matrix(table(MSSM$Donor, MSSM$predicted.id))
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

meta <- left_join(counts_MSSM, meta_MSSM, by = "Donor")

write.csv(meta, "Files/Pseudobulk_metadata_supertype_MSSM2.csv", row.names = FALSE)


# ==================== HBCC ====================

HBCC1 <- readRDS("/scratch/nendresz/PsychAD/Data/HBCC_1.rds")
HBCC2 <- readRDS("/scratch/nendresz/PsychAD/Data/HBCC_2.rds")

HBCC1$Age <- as.numeric(HBCC1$Age)
HBCC2$Age <- as.numeric(HBCC2$Age)

HBCC1 <- subset(HBCC1, subset = Age > 20 & Age < 70)
HBCC2 <- subset(HBCC2, subset = Age > 20 & Age < 70)

HBCC1@meta.data$Cohort <- "HBCC"
HBCC2@meta.data$Cohort <- "HBCC"

bulk1 <- AggregateExpression(HBCC1, group.by = c("predicted.id", "Donor"), return.seurat = TRUE)
bulk2 <- AggregateExpression(HBCC2, group.by = c("predicted.id", "Donor"), return.seurat = TRUE)

df1 <- GetAssayData(bulk1, assay = "RNA", layer = "counts")
df2 <- GetAssayData(bulk2, assay = "RNA", layer = "counts")

# Check duplicate columns
intersect(colnames(df1), colnames(df2))

df_combined <- cbind(df1, df2)

# Convert Ensembl IDs to gene symbols
genes <- rownames(df_combined)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = genes,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

na_genes <- is.na(gene_symbols)

if (any(na_genes)) {
  gene_symbols[na_genes] <- genes[na_genes]
}

rownames(df_combined) <- gene_symbols[rownames(df_combined)]
df_combined <- df_combined[!duplicated(rownames(df_combined)), ]

saveRDS(df_combined, "Files/HBCC_pseudobulk_SCZ_supertype.rds")

HBCC <- rbind(
  HBCC1@meta.data,
  HBCC2@meta.data
)

counts_HBCC <- as.data.frame.matrix(table(HBCC$Donor, HBCC$predicted.id))
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

meta <- left_join(counts_HBCC, meta_HBCC, by = "Donor")

write.csv(meta, "Files/Pseudobulk_metadata_supertype_HBCC.csv", row.names = FALSE)