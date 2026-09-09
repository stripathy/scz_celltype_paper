# Load required libraries
library(Matrix)

# Define the subset file parts
subset_parts <- as.character(1:20)

# Loop through each subset
for (part in subset_parts) {
    
    message("Processing subset: ", part)

    # File paths (adjust if inputs are already in Intermediate_dfs)
    matrix_file <- paste0("/scratch/nendresz/PsychAD/Data/Intermediate_dfs/psychAD_counts_part_", part, ".mtx")
    cell_file   <- paste0("/scratch/nendresz/PsychAD/Data/Intermediate_dfs/psychAD_cells_part_", part, ".csv")
    gene_file   <- paste0("/scratch/nendresz/PsychAD/Data/Intermediate_dfs/psychAD_genes_part_", part, ".csv")

    # Load sparse matrix
    sparse_matrix <- readMM(matrix_file)

    sparse_matrix <- t(sparse_matrix)

    # Load cell/gene metadata
    cell_metadata <- read.csv(cell_file, header = FALSE, stringsAsFactors = FALSE)
    gene_metadata <- read.csv(gene_file, header = FALSE, stringsAsFactors = FALSE)

    # Fix row/col names (genes as rows, cells as cols)
    rownames(sparse_matrix) <- gene_metadata[,1]   # genes
    colnames(sparse_matrix) <- cell_metadata[,1]   # cells

    # Convert to correct sparse format
    sparse_matrix <- as(sparse_matrix, "dgCMatrix")

    # Save into Intermediate_dfs
    saveRDS(
      sparse_matrix, 
      paste0("/scratch/nendresz/PsychAD/Data/Intermediate_dfs/psychAD_matrix_part_", part, ".rds")
    )

    message("Finished processing subset: ", part)
}

message("✅ All subsets processed and saved successfully in Intermediate_dfs!")



