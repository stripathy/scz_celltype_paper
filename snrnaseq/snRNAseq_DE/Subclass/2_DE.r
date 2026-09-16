# run on conda activate de_env
setwd("scz_celltype_paper/snrnaseq/snRNAseq_DE")
library(ggrepel)
library(cowplot)
library(limma)
library(dplyr)
library(edgeR)
library(tidyr)
library(dplyr)
library(EnhancedVolcano)
library(ggplot2)

# Cell types are every column that is not donor metadata. 1_Pseudobulk.r builds this
# table as left_join(table(Donor, Subclass), donor metadata), so the cell-type columns
# come first and the five metadata columns last — but selecting by position ([1:24])
# silently takes "Donor" as a cell type for any dataset with fewer subclasses.
META_COLS <- c("Donor", "Age", "Sex", "Diagnosis", "PMI")

bulk_files <- c("Files/McLean_pseudobulk_SCZ_subclass.rds", "Files/MSSM1_pseudobulk_SCZ_subclass.rds", "Files/Frohlich_pseudobulk_SCZ_subclass.rds",
                "Files/Batiuk_pseudobulk_SCZ_subclass.rds", "Files/MSSM2_pseudobulk_SCZ_subclass.rds", "Files/HBCC_pseudobulk_SCZ_subclass.rds")

meta_files <- c("Files/Pseudobulk_metadata_subclass_McLean.csv", "Files/Pseudobulk_metadata_subclass_MSSM1.csv", "Files/Pseudobulk_metadata_subclass_Frohlich.csv",
                "Files/Pseudobulk_metadata_subclass_Batiuk.csv",  "Files/Pseudobulk_metadata_subclass_MSSM2.csv", "Files/Pseudobulk_metadata_subclass_HBCC.csv")
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
cell_types <- setdiff(colnames(meta), META_COLS)
print(cell_types)

# total cells per donor
meta$total_cells <- rowSums(meta[, cell_types])

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

  # Diagnosis contrast named after the non-reference level — "DiagnosisSchizophrenia"
  # here, "DiagnosisSCZ" under the repo's Control/SCZ vocabulary. Read it off the design.
  dx_coef <- grep("^Diagnosis", colnames(design), value = TRUE)
  # Exactly one Diagnosis column means a clean two-level contrast. More than one
  # means Diagnosis has >2 levels in this fit — which happens if cohorts are pooled,
  # because Multiome spells its controls "control" and the rest spell them "Control".
  # Stop rather than silently contrast against the wrong reference.
  stopifnot(length(dx_coef) == 1)
  dx_coef <- dx_coef[1]
  
  # Extract DE results
  DE <- topTable(
    fit,
    coef = dx_coef,
    n = Inf,
    adjust.method = "BH"
  )
  
  safe_type <- gsub("/", "_", type)
  # Save results per subclass
  saveRDS(DE, paste0("Files/DE_results_", cohort, "_", safe_type, ".rds"))
  
}}


#####Multi has no pmi 

bulk_files <- c("Files/Multiome_pseudobulk_SCZ_subclass.rds")

meta_files <- c("Files/Pseudobulk_metadata_subclass_Multiome.csv")

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
cell_types <- setdiff(colnames(meta), META_COLS)

print(cell_types)

meta$total_cells <- rowSums(meta[, cell_types])

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

  # Diagnosis contrast named after the non-reference level — "DiagnosisSchizophrenia"
  # here, "DiagnosisSCZ" under the repo's Control/SCZ vocabulary. Read it off the design.
  dx_coef <- grep("^Diagnosis", colnames(design), value = TRUE)
  # Exactly one Diagnosis column means a clean two-level contrast. More than one
  # means Diagnosis has >2 levels in this fit — which happens if cohorts are pooled,
  # because Multiome spells its controls "control" and the rest spell them "Control".
  # Stop rather than silently contrast against the wrong reference.
  stopifnot(length(dx_coef) == 1)
  dx_coef <- dx_coef[1]
  
  # Extract DE results
  DE <- topTable(
    fit,
    coef = dx_coef,
    n = Inf,
    adjust.method = "BH"
  )
  
  
  safe_type <- gsub("/", "_", type)
  # Save results per subclass
  saveRDS(DE, paste0("Files/DE_results_", cohort, "_", safe_type, ".rds"))
  
}}

