library(Seurat)
library(dplyr)



path <- "/scratch/nendresz/PsychAD/Data/Seu_objs"
files <- list.files(path, pattern = "LT_psychAD_.*\\.rds", full.names = TRUE)


#Check which objects have which cohorts
for (f in files) {
  x <- readRDS(f)
  cat(basename(f), "→", unique(x$Source), "\n")
  rm(x); gc()
}

#MSSM
#1-14

#HBCC
#14-18

#RADC
#18-20



s1 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_1.rds")
s2  <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_2.rds")
s3  <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_3.rds")
s4  <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_4.rds")
s5  <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_5.rds")
s6  <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_6.rds")
s7  <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_7.rds")
s8  <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_8.rds")
s9  <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_9.rds")
s10 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_10.rds")
s11 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_11.rds")
s12 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_12.rds")
s13 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_13.rds")
s14 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_14.rds")
s15 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_15.rds")
s16 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_16.rds")
s17 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_17.rds")
s18 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_18.rds")
s19 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_19.rds")
s20 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_20.rds")


meta_subset <- read.csv("/project/rrg-shreejoy/PsychAD_NPS/meta/NPS-AD_individual_metadata.csv")


# ---------------- s1 ----------------
s1 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_1.rds")

s1@meta.data <- s1@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s1@meta.data) <- s1@meta.data$barcodekey
s1 <- subset(s1, subset = Source == "MSSM" & Diagnosis %in% c("Schizophrenia", "Control"))


# ---------------- s2 ----------------
s2 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_2.rds")

s2@meta.data <- s2@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s2@meta.data) <- s2@meta.data$barcodekey
s2 <- subset(s2, subset = Source == "MSSM" & Diagnosis %in% c("Schizophrenia", "Control"))


# ---------------- s3 ----------------
s3 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_3.rds")

s3@meta.data <- s3@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s3@meta.data) <- s3@meta.data$barcodekey
s3 <- subset(s3, subset = Source == "MSSM" & Diagnosis %in% c("Schizophrenia", "Control"))


# ---------------- s4 ----------------
s4 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_4.rds")

s4@meta.data <- s4@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s4@meta.data) <- s4@meta.data$barcodekey
s4 <- subset(s4, subset = Source == "MSSM" & Diagnosis %in% c("Schizophrenia", "Control"))


# ---------------- s5 ----------------
s5 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_5.rds")

s5@meta.data <- s5@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s5@meta.data) <- s5@meta.data$barcodekey
s5 <- subset(s5, subset = Source == "MSSM" & Diagnosis %in% c("Schizophrenia", "Control"))


# ---------------- s6 ----------------
s6 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_6.rds")

s6@meta.data <- s6@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>% select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s6@meta.data) <- s6@meta.data$barcodekey
s6 <- subset(s6, subset = Source == "MSSM" & Diagnosis %in% c("Schizophrenia", "Control"))


# ---------------- s7 ----------------
s7 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_7.rds")

s7@meta.data <- s7@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>% select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s7@meta.data) <- s7@meta.data$barcodekey
s7 <- subset(s7, subset = Source == "MSSM" & Diagnosis %in% c("Schizophrenia", "Control"))


# ---------------- s8 ----------------
s8 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_8.rds")

s8@meta.data <- s8@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>% select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s8@meta.data) <- s8@meta.data$barcodekey
s8 <- subset(s8, subset = Source == "MSSM" & Diagnosis %in% c("Schizophrenia", "Control"))


# ---------------- s9 ----------------
s9 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_9.rds")

s9@meta.data <- s9@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>% select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s9@meta.data) <- s9@meta.data$barcodekey
s9 <- subset(s9, subset = Source == "MSSM" & Diagnosis %in% c("Schizophrenia", "Control"))


# ---------------- s10 ----------------
s10 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_10.rds")

s10@meta.data <- s10@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>% select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s10@meta.data) <- s10@meta.data$barcodekey
s10 <- subset(s10, subset = Source == "MSSM" & Diagnosis %in% c("Schizophrenia", "Control"))


# ---------------- s11 ----------------
s11 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_11.rds")

s11@meta.data <- s11@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>% select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s11@meta.data) <- s11@meta.data$barcodekey
s11 <- subset(s11, subset = Source == "MSSM" & Diagnosis %in% c("Schizophrenia", "Control"))


# ---------------- s12 ----------------
s12 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_12.rds")

s12@meta.data <- s12@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>% select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s12@meta.data) <- s12@meta.data$barcodekey
s12 <- subset(s12, subset = Source == "MSSM" & Diagnosis %in% c("Schizophrenia", "Control"))


# ---------------- s13 ----------------
s13 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_13.rds")

s13@meta.data <- s13@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>% select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s13@meta.data) <- s13@meta.data$barcodekey
s13 <- subset(s13, subset = Source == "MSSM" & Diagnosis %in% c("Schizophrenia", "Control"))


# ---------------- s14 ----------------
s14 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_14.rds")

s14@meta.data <- s14@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>% select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s14@meta.data) <- s14@meta.data$barcodekey
s14 <- subset(s14, subset = Source == "MSSM" & Diagnosis %in% c("Schizophrenia", "Control"))



lapply(list(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14), dim)


merged_MSSM_1 <- merge(s1, y = list(s2, s3))

merged_MSSM_1 <- JoinLayers(merged_MSSM_1)

merged_MSSM_2 <- merge(s4, y = list(s5, s6, s7))

merged_MSSM_2 <- JoinLayers(merged_MSSM_2)

merged_MSSM_3 <- merge(s8, y = list(s9, s10))

merged_MSSM_3 <- JoinLayers(merged_MSSM_3)

merged_MSSM_4 <- merge(s11, y = list(s12, s13, s14))

merged_MSSM_4 <- JoinLayers(merged_MSSM_4)


Layers(merged_MSSM_1)
Layers(merged_MSSM_2)
Layers(merged_MSSM_3)
Layers(merged_MSSM_4)


lapply(list(merged_MSSM_1, merged_MSSM_2, merged_MSSM_3, merged_MSSM_4), dim)

saveRDS(merged_MSSM_1, "/scratch/nendresz/PsychAD/Data/MSSM_1.rds") 
saveRDS(merged_MSSM_2, "/scratch/nendresz/PsychAD/Data/MSSM_2.rds") 
saveRDS(merged_MSSM_3, "/scratch/nendresz/PsychAD/Data/MSSM_3.rds") 
saveRDS(merged_MSSM_4, "/scratch/nendresz/PsychAD/Data/MSSM_4.rds") 




##################################
##################################
##################################
###############HBCC###############
##################################
##################################
##################################


# --- Update metadata ---
s14@meta.data <- s14@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s14@meta.data) <- s14@meta.data$barcodekey

s14 <- subset(s14,subset = Source == "HBCC" & Diagnosis %in% c("Schizophrenia", "Control"))


# --- Update metadata ---
s15@meta.data <- s15@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s15@meta.data) <- s15@meta.data$barcodekey

s15 <- subset(s15,subset = Source == "HBCC" & Diagnosis %in% c("Schizophrenia", "Control"))


  # --- Update metadata ---
s16@meta.data <- s16@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s16@meta.data) <- s16@meta.data$barcodekey

s16 <- subset(s16,subset = Source == "HBCC" & Diagnosis %in% c("Schizophrenia", "Control"))

  # --- Update metadata ---
s17@meta.data <- s17@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s17@meta.data) <- s17@meta.data$barcodekey

s17 <- subset(s17,subset = Source == "HBCC" & Diagnosis %in% c("Schizophrenia", "Control"))




# --- Update metadata ---
s18@meta.data <- s18@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s18@meta.data) <- s18@meta.data$barcodekey

s18 <- subset(s18,subset = Source == "HBCC" & Diagnosis %in% c("Schizophrenia", "Control"))


merged_HBCC_1 <- merge(x = s14, y = list(s15))

merged_HBCC_1 <- JoinLayers(merged_HBCC_1)

merged_HBCC_2 <- merge(x = s16, y = list(s17, s18))

merged_HBCC_2 <- JoinLayers(merged_HBCC_2)

Layers(merged_HBCC_1)
Layers(merged_HBCC_2)

saveRDS(merged_HBCC_1, "/scratch/nendresz/PsychAD/Data/HBCC_1.rds") 
saveRDS(merged_HBCC_2, "/scratch/nendresz/PsychAD/Data/HBCC_2.rds") 




##################################################
##################################################
##################################################
##################################################
##################################################
# --- Update metadata ---
s18@meta.data <- s18@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s18@meta.data) <- s18@meta.data$barcodekey

s18 <- subset(s18,subset = Source == "RADC" & Diagnosis %in% c("AD", "Control"))


# --- Update metadata ---
s19@meta.data <- s19@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s19@meta.data) <- s19@meta.data$barcodekey

s19 <- subset(s19,subset = Source == "RADC" & Diagnosis %in% c("AD", "Control"))


  # --- Update metadata ---
s20@meta.data <- s20@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s20@meta.data) <- s20@meta.data$barcodekey

s20 <- subset(s20,subset = Source == "RADC" & Diagnosis %in% c("AD", "Control"))


merged_RADC <- merge(x = s18, y = list(s19, s20))

merged_RADC <- JoinLayers(merged_RADC)

Layers(merged_RADC)

saveRDS(merged_RADC, "/scratch/nendresz/PsychAD/Data/RADC.rds") 



#########################################
#########################################
#########################################
################AD######################
#########################################
#########################################

# --- Load metadata ---
meta_subset <- read.csv("/project/rrg-shreejoy/PsychAD_NPS/meta/NPS-AD_individual_metadata.csv")

# ---------------- s1 ----------------
s1 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_1.rds")

s1@meta.data <- s1@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>% select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s1@meta.data) <- s1@meta.data$barcodekey
s1 <- subset(s1, subset = Source == "MSSM" & Diagnosis %in% c("AD", "Control"))


# ---------------- s2 ----------------
s2 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_2.rds")

s2@meta.data <- s2@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>% select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s2@meta.data) <- s2@meta.data$barcodekey
s2 <- subset(s2, subset = Source == "MSSM" & Diagnosis %in% c("AD", "Control"))


# ---------------- s3 ----------------
s3 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_3.rds")

s3@meta.data <- s3@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>% select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s3@meta.data) <- s3@meta.data$barcodekey
s3 <- subset(s3, subset = Source == "MSSM" & Diagnosis %in% c("AD", "Control"))


# ---------------- s4 ----------------
s4 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_4.rds")

s4@meta.data <- s4@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>% select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s4@meta.data) <- s4@meta.data$barcodekey
s4 <- subset(s4, subset = Source == "MSSM" & Diagnosis %in% c("AD", "Control"))


# ---------------- s5 ----------------
s5 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_5.rds")

s5@meta.data <- s5@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>% select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s5@meta.data) <- s5@meta.data$barcodekey
s5 <- subset(s5, subset = Source == "MSSM" & Diagnosis %in% c("AD", "Control"))

# ---------------- s6 ----------------
s6 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_6.rds")

s6@meta.data <- s6@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s6@meta.data) <- s6@meta.data$barcodekey
s6 <- subset(s6, subset = Source == "MSSM" & Diagnosis %in% c("AD", "Control"))


# ---------------- s7 ----------------
s7 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_7.rds")

s7@meta.data <- s7@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s7@meta.data) <- s7@meta.data$barcodekey
s7 <- subset(s7, subset = Source == "MSSM" & Diagnosis %in% c("AD", "Control"))


# ---------------- s8 ----------------
s8 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_8.rds")

s8@meta.data <- s8@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s8@meta.data) <- s8@meta.data$barcodekey
s8 <- subset(s8, subset = Source == "MSSM" & Diagnosis %in% c("AD", "Control"))


# ---------------- s9 ----------------
s9 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_9.rds")

s9@meta.data <- s9@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s9@meta.data) <- s9@meta.data$barcodekey
s9 <- subset(s9, subset = Source == "MSSM" & Diagnosis %in% c("AD", "Control"))


# ---------------- s10 ----------------
s10 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_10.rds")

s10@meta.data <- s10@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s10@meta.data) <- s10@meta.data$barcodekey
s10 <- subset(s10, subset = Source == "MSSM" & Diagnosis %in% c("AD", "Control"))


# ---------------- s11 ----------------
s11 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_11.rds")

s11@meta.data <- s11@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s11@meta.data) <- s11@meta.data$barcodekey
s11 <- subset(s11, subset = Source == "MSSM" & Diagnosis %in% c("AD", "Control"))


# ---------------- s12 ----------------
s12 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_12.rds")

s12@meta.data <- s12@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s12@meta.data) <- s12@meta.data$barcodekey
s12 <- subset(s12, subset = Source == "MSSM" & Diagnosis %in% c("AD", "Control"))


# ---------------- s13 ----------------
s13 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_13.rds")

s13@meta.data <- s13@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s13@meta.data) <- s13@meta.data$barcodekey
s13 <- subset(s13, subset = Source == "MSSM" & Diagnosis %in% c("AD", "Control"))


# ---------------- s14 ----------------
s14 <- readRDS("/scratch/nendresz/PsychAD/Data/Seu_objs/LT_psychAD_14.rds")

s14@meta.data <- s14@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s14@meta.data) <- s14@meta.data$barcodekey
s14 <- subset(s14, subset = Source == "MSSM" & Diagnosis %in% c("AD", "Control"))

lapply(list(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14), dim)


meta_all <- do.call(rbind, lapply(list(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14),
                                    function(x) x@meta.data))
dim(meta_all)

# create split IDs
split_ids <- rep(1:ceiling(nrow(meta_all) / 550000), each = 550000)[1:nrow(meta_all)]
meta_split <- split(meta_all, split_ids)

# Check how many cells per group
sapply(meta_split, nrow)


cell_splits <- lapply(meta_split, rownames)

objs <- list(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14)
names(objs) <- paste0("s", 1:14)

# Create a list to store merged subsets
subset_list <- list()

for (i in seq_along(cell_splits)) {
  message("Processing Subset ", i, " ...")
  cells_i <- cell_splits[[i]]
  
  # Find which objects contain those cells
  objs_with_cells <- lapply(objs, function(x) intersect(cells_i, colnames(x)))
  objs_with_cells <- objs_with_cells[sapply(objs_with_cells, length) > 0]
  
  # Subset each relevant Seurat object
  parts <- mapply(function(x, cids) subset(x, cells = cids),
                  objs[names(objs_with_cells)],
                  objs_with_cells,
                  SIMPLIFY = FALSE)
  
  # Merge all pieces for that subset
  merged_subset <- Reduce(merge, parts)
  message("Subset ", i, " → ", ncol(merged_subset), " cells.")
  
  subset_list[[i]] <- merged_subset
}

# Check results
sapply(subset_list, ncol)

subset_list <- lapply(seq_along(subset_list), function(i) {
  message("Joining layers for subset ", i, " ...")
  JoinLayers(subset_list[[i]])
})

sapply(subset_list, dim)
sapply(subset_list, Layers)


for (i in seq_along(subset_list)) {
  saveRDS(subset_list[[i]], paste0("/scratch/nendresz/PsychAD/Data/MSSM_AD_part", i, ".rds"))
}






##################################
##################################
##################################
###########HBCC BPD###############
##################################
##################################
##################################


# --- Update metadata ---
s14@meta.data <- s14@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s14@meta.data) <- s14@meta.data$barcodekey

s14 <- subset(s14,subset = Source == "HBCC" & Diagnosis %in% c("Bipolar", "Control"))


# --- Update metadata ---
s15@meta.data <- s15@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s15@meta.data) <- s15@meta.data$barcodekey

s15 <- subset(s15,subset = Source == "HBCC" & Diagnosis %in% c("Bipolar", "Control"))


  # --- Update metadata ---
s16@meta.data <- s16@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s16@meta.data) <- s16@meta.data$barcodekey

s16 <- subset(s16,subset = Source == "HBCC" & Diagnosis %in% c("Bipolar", "Control"))

  # --- Update metadata ---
s17@meta.data <- s17@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s17@meta.data) <- s17@meta.data$barcodekey

s17 <- subset(s17,subset = Source == "HBCC" & Diagnosis %in% c("Bipolar", "Control"))




# --- Update metadata ---
s18@meta.data <- s18@meta.data %>%
  mutate(
    Diagnosis = case_when(
      AD_status == "Yes" ~ "AD",
      Schizophrenia == "Yes" ~ "Schizophrenia",
      Bipolar_Disorder == "Yes" ~ "Bipolar",
      Parkinson_disease == "Yes" ~ "Parkinson",
      Tardive_dyskinesia == "Yes" ~ "Tardive_dyskinesia",
      DLBD_status == "Yes" ~ "DLBD",
      FTD_status == "Yes" ~ "FTD",
      Tauopathy_status == "Yes" ~ "Tauopathy",
      Vascular_status == "Yes" ~ "Vascular",
      ASHCVD_status == "Yes" ~ "ASHCVD",
      TRUE ~ "Control"
    )
  ) %>%
  left_join(
    meta_subset %>%
      select(individualID, sex, ageDeath, apoeGenotype, diagnosis, PMI),
    by = c("SubID_export_synapse" = "individualID")
  ) %>%
  rename(Donor = SubID_export_synapse, Sex = sex, Age = ageDeath) %>%
  mutate(PMI = PMI / 60)

rownames(s18@meta.data) <- s18@meta.data$barcodekey

s18 <- subset(s18,subset = Source == "HBCC" & Diagnosis %in% c("Bipolar", "Control"))


merged_HBCC_1 <- merge(x = s14, y = list(s15))

merged_HBCC_1 <- JoinLayers(merged_HBCC_1)

#s16 on it's own 

merged_HBCC_2 <- merge(x = s17, y = list(s18))

merged_HBCC_2 <- JoinLayers(merged_HBCC_2)

Layers(merged_HBCC_1)
Layers(merged_HBCC_2)

saveRDS(merged_HBCC_1, "/scratch/nendresz/PsychAD/Data/HBCC_BPD_1.rds") 
saveRDS(s16, "/scratch/nendresz/PsychAD/Data/HBCC_BPD_2.rds") 
saveRDS(merged_HBCC_2, "/scratch/nendresz/PsychAD/Data/HBCC_BPD_3.rds") 

