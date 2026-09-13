# conda activate crumblr_env

# Packages
library(Seurat)
library(dplyr)
library(stringr)
library(tibble)
setwd("scz_celltype_paper/snrnaseq/Compositional_analysis")

# Load + age filter
HBCC1 <- readRDS("/scratch/nendresz/PsychAD/Data/HBCC_1.rds")
HBCC2 <- readRDS("/scratch/nendresz/PsychAD/Data/HBCC_2.rds")

HBCC1$Age <- as.numeric(HBCC1$Age)
HBCC2$Age <- as.numeric(HBCC2$Age)

HBCC1 <- subset(HBCC1, Age > 20 & Age < 70)
HBCC2 <- subset(HBCC2, Age > 20 & Age < 70)

MSSM1 <- subset(readRDS("/scratch/nendresz/PsychAD/Data/MSSM_1.rds"), Age < 70)
MSSM2 <- subset(readRDS("/scratch/nendresz/PsychAD/Data/MSSM_2.rds"), Age < 70)
MSSM3 <- subset(readRDS("/scratch/nendresz/PsychAD/Data/MSSM_3.rds"), Age < 70)
MSSM4 <- subset(readRDS("/scratch/nendresz/PsychAD/Data/MSSM_4.rds"), Age < 70)

OFC <- subset(readRDS("/project/rrg-shreejoy/nendresz/Supertypes/OFC_updated.rds"), Age < 70)
Bat <- subset(readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds"), Age < 70)
Ruz <- subset(readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds"), Age < 70)

Multi <- subset(
  readRDS("/scratch/nendresz/P1_Brain_scope/Files/Multiome.rds"),
  Age_death < 70 & Disorder %in% c("control", "Schizophrenia")
)

# Extract metadata
Ruz_McLean <- subset(Ruz, Cohort == "McLean")@meta.data
Ruz_MSSM1  <- subset(Ruz, Cohort == "MtSinai")@meta.data

hbcc <- rbind(HBCC1@meta.data, HBCC2@meta.data)
mssm <- rbind(MSSM1@meta.data, MSSM2@meta.data, MSSM3@meta.data, MSSM4@meta.data)
ofc <- OFC@meta.data
bat <- Bat@meta.data
Multi <- Multi@meta.data

# Donor metadata
collapse_meta <- function(df) {
  df %>%
    group_by(Donor) %>%
    summarise(
      Age = first(Age),
      Sex = first(Sex),
      Diagnosis = first(Diagnosis),
      PMI = first(PMI),
      .groups = "drop"
    )
}

hbcc_meta <- collapse_meta(hbcc)
mssm_meta <- collapse_meta(mssm)
ofc_meta <- collapse_meta(ofc)
bat_meta <- collapse_meta(bat)
Ruz_McLean_meta <- collapse_meta(Ruz_McLean)
Ruz_MSSM1_meta <- collapse_meta(Ruz_MSSM1)

hbcc_meta$Age <- as.numeric(hbcc_meta$Age)
mssm_meta$Age <- as.numeric(mssm_meta$Age)

# Multiome metadata
Multi_meta <- Multi %>%
  group_by(Donor) %>%
  summarise(
    Age = as.numeric(first(Age_death)),
    Sex = first(Biological_Sex),
    Diagnosis = first(Disorder),
    .groups = "drop"
  )

# Combine donors
all_meta <- bind_rows(
  HBCC = hbcc_meta,
  `MSSM 2` = mssm_meta,
  Fröhlich = ofc_meta,
  Batiuk = bat_meta,
  McLean = Ruz_McLean_meta,
  `MSSM 1` = Ruz_MSSM1_meta,
  Multiome = Multi_meta,
  .id = "Cohort"
)

# Cell counts
make_counts <- function(df, cohort) {
  as.data.frame.matrix(table(df$Donor, df$predicted.id)) %>%
    rownames_to_column("Donor") %>%
    mutate(Cohort = cohort)
}

counts_all <- bind_rows(
  make_counts(hbcc, "HBCC"),
  make_counts(mssm, "MSSM 2"),
  make_counts(ofc, "Fröhlich"),
  make_counts(bat, "Batiuk"),
  make_counts(Ruz_McLean, "McLean"),
  make_counts(Ruz_MSSM1, "MSSM 1"),
  make_counts(Multi, "Multiome")
)

# Add counts
all_meta_counts <- left_join(
  all_meta,
  counts_all,
  by = c("Donor", "Cohort")
)

# Standardize cell-type names
standardize_names <- function(x) {
  x %>%
    str_replace_all("\\.SEAAD$", "-SEAAD") %>%
    str_replace_all("Micro\\.PVM", "Micro-PVM") %>%
    str_replace_all("Sst\\.Chodl", "Sst Chodl") %>%
    str_replace_all("L2\\.3", "L2/3") %>%
    str_replace_all("L5\\.6", "L5/6") %>%
    str_replace_all("\\.", " ")
}

colnames(all_meta_counts) <- standardize_names(colnames(all_meta_counts))

# Save
write.csv(all_meta_counts, "Files/7_cohorts_metadata_names.csv", row.names = FALSE)