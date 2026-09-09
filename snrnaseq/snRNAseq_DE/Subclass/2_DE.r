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

bulk_files <- c("Files/Ruz_Mclean_pseudobulk_SCZ_subclass.rds", "Files/Ruz_MtSinai_pseudobulk_SCZ_subclass.rds", "Files/OFC_pseudobulk_SCZ_subclass.rds",
"Files/Bat_pseudobulk_SCZ_subclass.rds", "Files/MSSM_pseudobulk_SCZ_subclass.rds", "Files/HBCC_pseudobulk_SCZ_subclass.rds")

meta_files <- c("Files/Pseudobulk_metadata_subclass_Mclean.csv", "Files/Pseudobulk_metadata_subclass_MtSinai.csv", "Files/Pseudobulk_metadata_subclass_OFC.csv",
"Files/Pseudobulk_metadata_subclass_Bat.csv","Files/Pseudobulk_metadata_subclass_MSSM.csv","Files/Pseudobulk_metadata_subclass_HBCC.csv")


# Loop over each cell type

# loop over pseudobulk files
for (i in seq_along(bulk_files)) {
  
  bulk <- readRDS(bulk_files[i])
  meta <- read.csv(meta_files[i])
  
 cohort <- sub("Pseudobulk_metadata_subclass_(.*)\\.csv", "\\1", basename(meta_files[i]))

# make dots spaces

colnames(meta) <- gsub("\\.", " ", colnames(meta))
colnames(meta) <- gsub("^L2 3", "L2/3", colnames(meta))
colnames(meta) <- gsub("^L5 6", "L5/6", colnames(meta))
colnames(meta) <- gsub("^Micro PVM", "Micro-PVM", colnames(meta))


# get cell types 
cell_types <- colnames(meta)[1:24]
print(cell_types)

# total cells per donor
meta$total_cells <- rowSums(meta[, cell_types])
cell_types <- colnames(meta)[1:24]

meta$total_cells <- rowSums(meta [1:24])

for (type in cell_types) {
  
  message("Processing subclass: ", type, " from ", cohort)

#For HBCC and MSSM 
meta <- meta %>%
  mutate(Donor = gsub("_", "-", Donor))

  # Keep donors with 500 total cells
 donors_keep <- meta %>% filter(total_cells >= 500) %>% pull(Donor)

  # Subset bulk for cell type
  bulk_subset <- bulk[, grep(paste0("^", type, "_"), colnames(bulk))]
  
  if (ncol(bulk_subset) == 0) next
  
  pseudobulk <- bulk_subset %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID_Celltype") %>%
    separate(ID_Celltype, into = c("subclass", "Donor"), sep = "_", extra = "merge", fill = "right") %>%
    relocate(subclass, Donor) %>%
    mutate(Donor = gsub("_", "-", Donor)) %>%
    filter(Donor %in% donors_keep)
  
  # Create matrix for DGEList
  mat <- pseudobulk %>%
    tibble::column_to_rownames("Donor") %>%
    dplyr::select(-subclass) %>%
    as.matrix() %>%
    t()
  
  # Reorder metadata to match columns
  meta_subset <- meta[match(colnames(mat), meta$Donor), ]
  
  #add log cells
  meta_subset$log2_cells <- scale(log2(meta_subset[[type]]))

     # Skip if Diagnosis has <2 levels
  if (length(unique(meta_subset$Diagnosis)) < 2) {
    message("Skipping ", type, " because Diagnosis has <2 levels")
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
  # Save results per subclass
  saveRDS(DE, paste0("Files/DE_results_", cohort, "_", safe_type, ".rds"))
  
}}













#####Multi has no pmi 


bulk_files <- c("Files/Multiome_pseudobulk_SCZ_subclass.rds")

meta_files <- c("Files/Pseudobulk_metadata_subclass_Multi.csv")

# Loop over each cell type

# loop over pseudobulk files
for (i in seq_along(bulk_files)) {
  
  bulk <- readRDS(bulk_files[i])
  meta <- read.csv(meta_files[i])
  
 cohort <- sub("Pseudobulk_metadata_subclass_(.*)\\.csv", "\\1", basename(meta_files[i]))
colnames(meta) <- gsub("\\.", " ", colnames(meta))
colnames(meta) <- gsub("^L2 3", "L2/3", colnames(meta))
colnames(meta) <- gsub("^L5 6", "L5/6", colnames(meta))
colnames(meta) <- gsub("^Micro PVM", "Micro-PVM", colnames(meta))

# get cell types 
cell_types <- colnames(meta)[1:23]

print(cell_types)

meta$total_cells <- rowSums(meta [1:23])

for (type in cell_types) {
  
  message("Processing subclass: ", type, " from ", cohort)

#For HBCC and MSSM 
meta <- meta %>%
  mutate(Donor = gsub("_", "-", Donor))

  # Keep donors with 500 total cells
 donors_keep <- meta %>% filter(total_cells >= 500) %>% pull(Donor)

  # Subset bulk for cell type
  bulk_subset <- bulk[, grep(paste0("^", type, "_"), colnames(bulk))]
  
  if (ncol(bulk_subset) == 0) next
  
  pseudobulk <- bulk_subset %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID_Celltype") %>%
    separate(ID_Celltype, into = c("subclass", "Donor"), sep = "_", extra = "merge", fill = "right") %>%
    relocate(subclass, Donor) %>%
    mutate(Donor = gsub("_", "-", Donor)) %>%
    filter(Donor %in% donors_keep)
  
  # Create matrix for DGEList
  mat <- pseudobulk %>%
    tibble::column_to_rownames("Donor") %>%
    dplyr::select(-subclass) %>%
    as.matrix() %>%
    t()
  
  # Reorder metadata to match columns
  meta_subset <- meta[match(colnames(mat), meta$Donor), ]

  #add log cells
  meta_subset$log2_cells <- scale(log2(meta_subset[[type]]))
  
     # Skip if Diagnosis has <2 levels
  if (length(unique(meta_subset$Diagnosis)) < 2) {
    message("Skipping ", type, " because Diagnosis has <2 levels")
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
  # Save results per subclass
  saveRDS(DE, paste0("Files/DE_results_", cohort, "_", safe_type, ".rds"))
  
}}

