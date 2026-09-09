# run on conda activate de_env
setwd("P1_SCZ_DE_fresh")
library(ggrepel)
library(cowplot)
library(limma)
library(dplyr)
library(edgeR)
library(tidyr)
library(dplyr)
library(EnhancedVolcano)
library(ggplot2)

bulk_files <- c("Files/Ruz_Mclean_pseudobulk_SCZ_supertype.rds", "Files/Ruz_MtSinai_pseudobulk_SCZ_supertype.rds", "Files/OFC_pseudobulk_SCZ_supertype.rds",
"Files/Bat_pseudobulk_SCZ_supertype.rds", "Files/MSSM_pseudobulk_SCZ_supertype.rds", "Files/HBCC_pseudobulk_SCZ_supertype.rds")

meta_files <- c("Files/Pseudobulk_metadata_supertype_Mclean.csv", "Files/Pseudobulk_metadata_supertype_MtSinai.csv", "Files/Pseudobulk_metadata_supertype_OFC.csv",
"Files/Pseudobulk_metadata_supertype_Bat.csv","Files/Pseudobulk_metadata_supertype_MSSM.csv","Files/Pseudobulk_metadata_supertype_HBCC.csv")


# Loop over each cell type

# loop over pseudobulk files
for (i in seq_along(bulk_files)) {
  
  bulk <- readRDS(bulk_files[i])
  meta <- read.csv(meta_files[i])
  
 cohort <- sub("Pseudobulk_metadata_supertype_(.*)\\.csv", "\\1", basename(meta_files[i]))

# make dots spaces

colnames(meta)  <- gsub("\\.", " ", colnames(meta))
colnames(meta)  <- gsub("\\_", "-", colnames(meta))
colnames(meta)  <- gsub("^L2 3", "L2/3", colnames(meta))
colnames(meta)  <- gsub("^L5 6", "L5/6", colnames(meta))
colnames(meta) <- gsub("^Micro PVM", "Micro-PVM", colnames(meta))

#get cell types
meta_cols <- c("Donor", "Age", "Sex", "Diagnosis", "PMI")
cell_types <- setdiff(colnames(meta), meta_cols)
print(cell_types)


for (type in cell_types) {
  
  message("Processing supertype: ", type, " from ", cohort)

#For HBCC and MSSM 
meta <- meta %>%
  mutate(Donor = gsub("_", "-", Donor))

  # Keep donors with 10 of that cell type
donors_keep <- meta %>% filter(.data[[type]] >= 10) %>% pull(Donor)

bulk_subset <- bulk[, grep(paste0("^", type, "_"), colnames(bulk)), drop = FALSE]

if (ncol(bulk_subset) == 0) next
if (Matrix::nnzero(bulk_subset) == 0) next
  
  pseudobulk <- bulk_subset %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID_Celltype") %>%
    separate(ID_Celltype, into = c("supertype", "Donor"), sep = "_", extra = "merge", fill = "right") %>%
    relocate(supertype, Donor) %>%
    mutate(Donor = gsub("_", "-", Donor)) %>%
    filter(Donor %in% donors_keep)
  
  # Create matrix for DGEList
  mat <- pseudobulk %>%
    tibble::column_to_rownames("Donor") %>%
    dplyr::select(-supertype) %>%
    as.matrix() %>%
    t()
  
  # Reorder metadata to match columns
  meta_subset <- meta[match(colnames(mat), meta$Donor), ]
  
     # Skip if Diagnosis has <2 levels
  if (length(unique(meta_subset$Diagnosis)) < 2) {
    message("Skipping ", type, " because Diagnosis has <2 levels")
    next
  }

  if (length(unique(meta_subset$Donor)) < 2) {
    message("Skipping ", type, " because Donor has <2 levels")
    next
  }

    if (length(unique(meta_subset$Sex)) < 2) {
    message("Skipping ", type, " because Sex has <2 levels")
    next
  }

 if(nrow(meta_subset) < 5) {
  message(type, ": skipping (too few samples)")
  next
}

  # DGEList
  dge <- DGEList(counts = mat, genes = rownames(mat))
  
  # Keep genes with ≥1 count in ≥80% of samples
  min_samples <- ncol(mat) * 0.8
  dge <- dge[rowSums(dge$counts >= 1) >= min_samples, ]
  
  # Normalization
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Design matrix
  design <- model.matrix(~ scale(Age) + Sex + Diagnosis + scale(PMI), data = meta_subset)
  
if (nrow(design) <= ncol(design)) {
  message("Skipping ", type, " (design not estimable: no residual df)")
  next
}

  # voom + lmFit + eBayes
  vm <- voom(dge, design, plot = FALSE)
  fit <- lmFit(vm, design)
  fit <- eBayes(fit)
  
  # Extract DE results
  DE <- topTable(
    fit,
    coef = "DiagnosisSchizophrenia",
    n = Inf,
    adjust.method = "BH"
  )
  
  safe_type <- gsub("/", "_", type)
  # Save results per supertype
  saveRDS(DE, paste0("Files/DE_results_", cohort, "_", safe_type, ".rds"))
  
}}









#####Multi has no pmi 


bulk_files <- c("Files/Multiome_pseudobulk_SCZ_supertype.rds")

meta_files <- c("Files/Pseudobulk_metadata_supertype_Multi.csv")

for (i in seq_along(bulk_files)) {
  
  bulk <- readRDS(bulk_files[i])
  meta <- read.csv(meta_files[i])
  
 cohort <- sub("Pseudobulk_metadata_supertype_(.*)\\.csv", "\\1", basename(meta_files[i]))

# make dots spaces

colnames(meta)  <- gsub("\\.", " ", colnames(meta))
colnames(meta)  <- gsub("\\_", "-", colnames(meta))
colnames(meta)  <- gsub("^L2 3", "L2/3", colnames(meta))
colnames(meta)  <- gsub("^L5 6", "L5/6", colnames(meta))
colnames(meta) <- gsub("^Micro PVM", "Micro-PVM", colnames(meta))

#get cell types
meta_cols <- c("Donor", "Age", "Sex", "Diagnosis", "PMI")
cell_types <- setdiff(colnames(meta), meta_cols)
print(cell_types)


for (type in cell_types) {
  
  message("Processing supertype: ", type, " from ", cohort)

#For HBCC and MSSM 
meta <- meta %>%
  mutate(Donor = gsub("_", "-", Donor))

  # Keep donors with 10 of that cell type
donors_keep <- meta %>% filter(.data[[type]] >= 10) %>% pull(Donor)

bulk_subset <- bulk[, grep(paste0("^", type, "_"), colnames(bulk)), drop = FALSE]

if (ncol(bulk_subset) == 0) next
if (Matrix::nnzero(bulk_subset) == 0) next
  
  pseudobulk <- bulk_subset %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID_Celltype") %>%
    separate(ID_Celltype, into = c("supertype", "Donor"), sep = "_", extra = "merge", fill = "right") %>%
    relocate(supertype, Donor) %>%
    mutate(Donor = gsub("_", "-", Donor)) %>%
    filter(Donor %in% donors_keep)
  
  # Create matrix for DGEList
  mat <- pseudobulk %>%
    tibble::column_to_rownames("Donor") %>%
    dplyr::select(-supertype) %>%
    as.matrix() %>%
    t()
  
  # Reorder metadata to match columns
  meta_subset <- meta[match(colnames(mat), meta$Donor), ]
  
     # Skip if Diagnosis has <2 levels
  if (length(unique(meta_subset$Diagnosis)) < 2) {
    message("Skipping ", type, " because Diagnosis has <2 levels")
    next
  }

  if (length(unique(meta_subset$Donor)) < 2) {
    message("Skipping ", type, " because Donor has <2 levels")
    next
  }

    if (length(unique(meta_subset$Sex)) < 2) {
    message("Skipping ", type, " because Sex has <2 levels")
    next
  }

 if(nrow(meta_subset) < 5) {
  message(type, ": skipping (too few samples)")
  next
}

  # DGEList
  dge <- DGEList(counts = mat, genes = rownames(mat))
  
  # Keep genes with ≥1 count in ≥80% of samples
  min_samples <- ncol(mat) * 0.8
  dge <- dge[rowSums(dge$counts >= 1) >= min_samples, ]
  
  # Normalization
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Design matrix
  design <- model.matrix(~ scale(Age) + Sex + Diagnosis, data = meta_subset)
  
if (nrow(design) <= ncol(design)) {
  message("Skipping ", type, " (design not estimable: no residual df)")
  next
}

  # voom + lmFit + eBayes
  vm <- voom(dge, design, plot = FALSE)
  fit <- lmFit(vm, design)
  fit <- eBayes(fit)
  
  # Extract DE results
  DE <- topTable(
    fit,
    coef = "DiagnosisSchizophrenia",
    n = Inf,
    adjust.method = "BH"
  )
  
  safe_type <- gsub("/", "_", type)
  # Save results per supertype
  saveRDS(DE, paste0("Files/DE_results_", cohort, "_", safe_type, ".rds"))
  
}}



