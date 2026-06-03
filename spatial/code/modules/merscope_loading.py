"""
Data loading utilities for MERSCOPE spatial transcriptomics data.

Handles loading Fang et al. MERSCOPE datasets (genes.csv / features.csv /
matrix.csv triplet format) and constructing AnnData objects with spatial
coordinates and pre-existing cluster labels.

Data source: Fang et al. (Dryad: doi:10.5061/dryad.x3ffbg7mw)
Format: sparse triplet CSVs (1-indexed row=gene, col=cell, val=count)
"""

import os
import re
import glob
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse


def load_merscope_sample(prefix):
    """
    Load a single MERSCOPE sample from triplet CSV files into AnnData.

    Reads three files:
      - {prefix}.genes.csv: gene index -> gene name mapping
      - {prefix}.features.csv: cell metadata (x/y, cluster labels)
      - {prefix}.matrix.csv: sparse triplet (row=gene, col=cell, val=count), 1-indexed

    Parameters
    ----------
    prefix : str
        Path prefix (without .genes.csv / .features.csv / .matrix.csv).

    Returns
    -------
    anndata.AnnData
        AnnData with raw counts in .X (sparse CSR), spatial coords in
        .obsm['spatial'], cluster labels in .obs, and sample_id in .obs.
    """
    genes_path = f"{prefix}.genes.csv"
    features_path = f"{prefix}.features.csv"
    matrix_path = f"{prefix}.matrix.csv"

    for p in [genes_path, features_path, matrix_path]:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing: {p}")

    # Gene names
    genes_df = pd.read_csv(genes_path)
    gene_names = list(genes_df["name"])

    # Cell metadata (x, y, cluster labels)
    features_df = pd.read_csv(features_path, index_col=0)
    n_cells = len(features_df)
    n_genes = len(gene_names)

    # Sparse count matrix (triplet format, 1-indexed)
    matrix_df = pd.read_csv(matrix_path)
    rows = matrix_df["row"].values - 1  # gene index (0-based)
    cols = matrix_df["col"].values - 1  # cell index (0-based)
    vals = matrix_df["val"].values

    X = sparse.csr_matrix((vals, (cols, rows)), shape=(n_cells, n_genes))

    # Build AnnData
    adata = ad.AnnData(X=X, obs=features_df.copy())
    adata.var_names = gene_names
    adata.var_names_make_unique()

    # Spatial coordinates
    if "global.x" in features_df.columns and "global.y" in features_df.columns:
        adata.obsm["spatial"] = features_df[["global.x", "global.y"]].values
    elif "adjusted.x" in features_df.columns and "adjusted.y" in features_df.columns:
        adata.obsm["spatial"] = features_df[["adjusted.x", "adjusted.y"]].values

    # Extract sample ID from filename
    sample_id = _extract_sample_id(prefix)
    adata.obs["sample_id"] = sample_id

    # Make obs_names unique across samples
    adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]

    return adata


def discover_merscope_samples(data_dir, panel_size=None):
    """
    Discover all MERSCOPE samples in a directory.

    Looks for *.genes.csv files and extracts sample metadata from filenames.
    Filename format: {donor}.{region}.{panel_size}.{expand_type}.{rep}.genes.csv
    Example: H18.06.006.MTG.4000.expand.rep1.genes.csv

    Parameters
    ----------
    data_dir : str
        Directory containing MERSCOPE data files.
    panel_size : int or None
        If specified, filter to only this panel size (250 or 4000).

    Returns
    -------
    list of dict
        Each dict has: 'prefix', 'sample_id', 'donor', 'region',
        'panel_size', 'expand_type', 'replicate'.
    """
    gene_files = sorted(glob.glob(os.path.join(data_dir, "*.genes.csv")))
    samples = []

    for gf in gene_files:
        prefix = gf.replace(".genes.csv", "")
        basename = os.path.basename(prefix)
        parts = basename.split(".")

        # Parse: H18.06.006.MTG.4000.expand.rep1
        if len(parts) < 7:
            print(f"  Warning: unexpected filename format: {basename}")
            continue

        donor = ".".join(parts[:3])
        region = parts[3]
        psize = int(parts[4])
        expand_type = parts[5]
        replicate = parts[6]

        if panel_size is not None and psize != panel_size:
            continue

        sample_id = f"{donor}_{region}_{psize}_{replicate}"

        samples.append({
            "prefix": prefix,
            "sample_id": sample_id,
            "donor": donor,
            "region": region,
            "panel_size": psize,
            "expand_type": expand_type,
            "replicate": replicate,
        })

    return samples


def load_all_merscope_samples(data_dir, panel_size=4000):
    """
    Load all MERSCOPE samples with a given panel size and merge.

    Parameters
    ----------
    data_dir : str
        Directory containing MERSCOPE data files.
    panel_size : int
        Panel size to load (default: 4000).

    Returns
    -------
    anndata.AnnData
        Merged AnnData with raw counts in .X, spatial coords in
        .obsm['spatial'], sample_id in .obs.
    """
    samples = discover_merscope_samples(data_dir, panel_size=panel_size)
    print(f"Found {len(samples)} MERSCOPE {panel_size}-gene samples")

    adatas = []
    for info in samples:
        print(f"  Loading {info['sample_id']}...", end=" ", flush=True)
        adata = load_merscope_sample(info["prefix"])
        print(f"{adata.n_obs:,} cells x {adata.n_vars} genes")
        adatas.append(adata)

    if not adatas:
        raise ValueError(f"No MERSCOPE samples found with panel_size={panel_size}")

    combined = ad.concat(adatas, join="inner")
    combined.layers["counts"] = combined.X.copy()

    # Reconstruct spatial from individual samples
    spatial = np.vstack([a.obsm["spatial"] for a in adatas])
    combined.obsm["spatial"] = spatial

    print(f"Combined: {combined.n_obs:,} cells x {combined.n_vars} genes")
    return combined


def _extract_sample_id(prefix):
    """Extract a short sample ID from a MERSCOPE prefix path."""
    basename = os.path.basename(prefix)
    parts = basename.split(".")
    if len(parts) >= 7:
        donor = ".".join(parts[:3])
        region = parts[3]
        psize = parts[4]
        rep = parts[6]
        return f"{donor}_{region}_{psize}_{rep}"
    return basename
