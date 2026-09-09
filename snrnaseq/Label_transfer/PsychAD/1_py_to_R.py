#conda activate py_anndata_env 

import anndata as ad
import os
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

# Define the output directory
output_dir = "/scratch/nendresz/P1_SCZ_paper/Label_transfer/PsychAD/data"
os.makedirs(output_dir, exist_ok=True)

#PsychAD dataset directory
os.chdir("/project/rrg-shreejoy/nendresz/PsychAD_Data")

# Load dataset in backed mode
adata = ad.read_h5ad("merged_final_clean.h5ad", backed="r")

# Parameters
num_cells = adata.n_obs
num_splits = 20  # Keeping 20 subsets for memory efficiency
split_size = num_cells // num_splits  # ~316,000 cells per part

for i in range(num_splits):
    start = i * split_size
    end = start + split_size if i < num_splits - 1 else num_cells
    
    print(f"Processing subset {i+1}/{num_splits} (cells {start} to {end})...")

    # Load subset into memory first (to speed up processing)
    adata_subset = adata[start:end, :].to_memory()

    # Convert to sparse matrix
    sparse_matrix = csr_matrix(adata_subset.layers["counts"])


    # Save as Matrix Market format
    mtx_filename = os.path.join(output_dir, f"psychAD_counts_part_{i+1}.mtx")
    mmwrite(mtx_filename, sparse_matrix)

    # Save gene names as CSV
    genes_filename = os.path.join(output_dir, f"psychAD_genes_part_{i+1}.csv")
    pd.DataFrame(adata_subset.var.index).to_csv(genes_filename, index=False, header=False)

    # Save cell names as CSV
    cells_filename = os.path.join(output_dir, f"psychAD_cells_part_{i+1}.csv")
    pd.DataFrame(adata_subset.obs.index).to_csv(cells_filename, index=False, header=False)
    
    # Free up memory

    del adata_subset, sparse_matrix

    print(f"Saved: {mtx_filename}, {genes_filename}, {cells_filename}")

# after the for loop
adata.file.close()
print("✅ All subsets processed and H5AD file closed cleanly")