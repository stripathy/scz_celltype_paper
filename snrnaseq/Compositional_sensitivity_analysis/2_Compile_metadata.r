library(dplyr)
library(tibble)
library(stringr)

# Load new metadata
hbcc <- read.csv("P1_Controls/Files/HBCC_SST_DEgenes.csv",row.names=1,check.names=FALSE)
mssm <- read.csv("P1_Controls/Files/MSSM2_SST_DEgenes.csv",row.names=1,check.names=FALSE)
ofc <- read.csv("P1_Controls/Files/Frohlich_NoSST_DEgenes.csv",row.names=1,check.names=FALSE)
bat <- read.csv("P1_Controls/Files/Batiuk_NoSST_DEgenes.csv",row.names=1,check.names=FALSE)
Ruz_McLean <- read.csv("P1_Controls/Files/McLean_NoSST_DEgenes.csv",row.names=1,check.names=FALSE)
Ruz_MtSinai <- read.csv("P1_Controls/Files/MSSM1_NoSST_DEgenes.csv",row.names=1,check.names=FALSE)
Multi <- read.csv("P1_Controls/Files/Multiome_NoSST_DEgenes.csv",row.names=1,check.names=FALSE)

# Donor metadata
collapse_meta <- function(df) df %>% group_by(Donor) %>% summarise(Age=first(Age),Sex=first(Sex),Diagnosis=first(Diagnosis),PMI=first(PMI),.groups="drop")

hbcc_meta <- collapse_meta(hbcc); mssm_meta <- collapse_meta(mssm); ofc_meta <- collapse_meta(ofc); bat_meta <- collapse_meta(bat)
Ruz_McLean_meta <- collapse_meta(Ruz_McLean); Ruz_MtSinai_meta <- collapse_meta(Ruz_MtSinai)
hbcc_meta$Age <- as.numeric(hbcc_meta$Age); mssm_meta$Age <- as.numeric(mssm_meta$Age)

Multi_meta <- Multi %>% group_by(Donor) %>% summarise(Age=first(Age_death),Sex=first(Biological_Sex),Diagnosis=first(Disorder),.groups="drop")
Multi_meta$Age <- as.numeric(Multi_meta$Age); Multi_meta$Donor <- as.character(Multi_meta$Donor)

all_meta <- bind_rows(list(HBCC=hbcc_meta,MSSM=mssm_meta,OFC=ofc_meta,Batiuk=bat_meta,McLean=Ruz_McLean_meta,MtSinai=Ruz_MtSinai_meta,Multiome=Multi_meta),.id="Cohort")

# Cell counts from new predicted.id
make_counts <- function(df,cohort) as.data.frame.matrix(table(df$Donor,df$predicted.id)) %>% rownames_to_column("Donor") %>% mutate(Cohort=cohort)

counts_all <- bind_rows(make_counts(hbcc,"HBCC"),make_counts(mssm,"MSSM"),make_counts(ofc,"OFC"),make_counts(bat,"Batiuk"),
                        make_counts(Ruz_McLean,"McLean"),make_counts(Ruz_MtSinai,"MtSinai"),make_counts(Multi,"Multiome"))

all_meta$Donor <- as.character(all_meta$Donor); counts_all$Donor <- as.character(counts_all$Donor)
all_meta_counts <- left_join(all_meta,counts_all,by=c("Donor","Cohort"))

# Standardize cell-type names
standardize_names <- function(x) x %>% str_replace_all("\\.SEAAD$","-SEAAD") %>% str_replace_all("Micro\\.PVM","Micro-PVM") %>%
  str_replace_all("Sst\\.Chodl","Sst Chodl") %>% str_replace_all("L2\\.3","L2/3") %>% str_replace_all("L5\\.6","L5/6") %>% str_replace_all("\\."," ")

colnames(all_meta_counts) <- standardize_names(colnames(all_meta_counts))

write.csv(all_meta_counts,"P1_Controls/Files/7_cohorts_metadata_NoSST_DEgenes_names.csv",row.names=FALSE)