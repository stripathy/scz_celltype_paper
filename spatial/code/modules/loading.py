"""
Data loading utilities for Xenium spatial transcriptomics data.

Handles loading 10x Xenium .h5 files, cell boundary CSVs, and
constructing AnnData objects with spatial coordinates.
"""

import os
import re
import glob
import numpy as np
import pandas as pd
import anndata as ad
import h5py
from scipy.sparse import csc_matrix


def load_xenium_h5(h5_path):
    """
    Load a 10x Xenium cell_feature_matrix.h5 into an AnnData object.

    Reads the sparse count matrix, gene names, and cell barcodes.
    Filters to 'Gene Expression' features only (excludes controls).

    Parameters
    ----------
    h5_path : str
        Path to the .h5 file.

    Returns
    -------
    anndata.AnnData
        AnnData with raw counts in .X, gene names in .var_names,
        cell barcodes in .obs_names.
    """
    with h5py.File(h5_path, "r") as f:
        matrix = f["matrix"]
        barcodes = [b.decode() for b in matrix["barcodes"][:]]
        features = matrix["features"]
        gene_names = [g.decode() for g in features["name"][:]]
        feature_type = [t.decode() for t in features["feature_type"][:]]

        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        shape = matrix["shape"][:]

        X = csc_matrix((data, indices, indptr), shape=shape).T

    # Filter to Gene Expression features only
    gene_mask = np.array([t == "Gene Expression" for t in feature_type])
    gene_names_filtered = [g for g, m in zip(gene_names, gene_mask) if m]
    X_filtered = X[:, gene_mask]

    adata = ad.AnnData(
        X=X_filtered,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=gene_names_filtered),
    )
    return adata


def load_cell_boundaries(csv_path):
    """
    Load cell boundary polygon CSV and compute cell centroids.

    Parameters
    ----------
    csv_path : str
        Path to cell_boundaries.csv.gz file.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['cell_id', 'x', 'y'] where x and y
        are centroid coordinates.
    """
    df = pd.read_csv(csv_path)
    centroids = df.groupby("cell_id").agg(
        x=("vertex_x", "mean"),
        y=("vertex_y", "mean"),
    ).reset_index()
    return centroids


def load_xenium_sample(h5_path, boundaries_path):
    """
    Load a single Xenium sample with spatial coordinates.

    Combines the expression matrix from the .h5 file with spatial
    centroids computed from the cell boundary CSV.

    Parameters
    ----------
    h5_path : str
        Path to the cell_feature_matrix.h5 file.
    boundaries_path : str
        Path to the cell_boundaries.csv.gz file.

    Returns
    -------
    anndata.AnnData
        AnnData with counts in .X, spatial coords in .obsm['spatial'],
        and 'sample_id' in .obs.
    """
    adata = load_xenium_h5(h5_path)
    centroids = load_cell_boundaries(boundaries_path)

    # Match centroids to barcodes by position (1-indexed cell_id)
    centroids = centroids.sort_values("cell_id").reset_index(drop=True)

    # Align: cell_id in boundaries is 1-indexed, matching barcode order
    n_cells = min(len(adata), len(centroids))
    adata = adata[:n_cells].copy()
    centroids = centroids.iloc[:n_cells]

    adata.obsm["spatial"] = centroids[["x", "y"]].values

    # Extract sample ID from filename
    sample_id = _extract_sample_id(h5_path)
    adata.obs["sample_id"] = sample_id

    return adata


def discover_samples(data_dir):
    """
    Discover all Xenium samples in a directory.

    Looks for GSM*-cell_feature_matrix.h5 files and their corresponding
    cell_boundaries.csv.gz files.

    Parameters
    ----------
    data_dir : str
        Directory containing the Xenium data files.

    Returns
    -------
    list of dict
        Each dict has keys: 'sample_id', 'h5_path', 'boundaries_path'.
    """
    h5_files = sorted(glob.glob(os.path.join(data_dir, "GSM*-cell_feature_matrix.h5")))
    samples = []
    for h5_path in h5_files:
        sample_id = _extract_sample_id(h5_path)
        # Find matching boundaries file
        prefix = os.path.basename(h5_path).split("-cell_feature_matrix")[0]
        bounds_path = os.path.join(data_dir, f"{prefix}-cell_boundaries.csv.gz")
        if os.path.exists(bounds_path):
            samples.append({
                "sample_id": sample_id,
                "h5_path": h5_path,
                "boundaries_path": bounds_path,
            })
        else:
            print(f"Warning: no boundaries file for {sample_id}, skipping")
    return samples


def load_all_samples(data_dir):
    """
    Load all Xenium samples and merge into a single AnnData.

    Parameters
    ----------
    data_dir : str
        Directory containing all Xenium data files.

    Returns
    -------
    anndata.AnnData
        Merged AnnData with all samples. Raw counts in .layers['counts'],
        spatial coords in .obsm['spatial'], sample_id in .obs.
    """
    samples = discover_samples(data_dir)
    print(f"Found {len(samples)} samples")

    adatas = []
    for info in samples:
        print(f"  Loading {info['sample_id']}...")
        adata = load_xenium_sample(info["h5_path"], info["boundaries_path"])
        # Make barcodes unique across samples
        adata.obs_names = [
            f"{info['sample_id']}_{bc}" for bc in adata.obs_names
        ]
        adatas.append(adata)

    combined = ad.concat(adatas, join="inner")
    combined.layers["counts"] = combined.X.copy()

    # Reconstruct spatial from individual samples
    spatial = np.vstack([a.obsm["spatial"] for a in adatas])
    combined.obsm["spatial"] = spatial

    print(f"Combined: {combined.shape[0]:,} cells x {combined.shape[1]} genes")
    return combined


def _extract_sample_id(h5_path):
    """Extract sample ID (e.g. 'Br8667') from a Xenium h5 filename."""
    basename = os.path.basename(h5_path)
    match = re.search(r"(Br\d+)", basename)
    if match:
        return match.group(1)
    return basename.split("-")[0]
