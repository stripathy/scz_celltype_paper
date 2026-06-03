"""
spatial.py — Layer-stratified analysis of SST subtype SCZ enrichment.

Uses spatial transcriptomics depth data (MERFISH + Xenium) to classify SST
supertypes by their cortical position. Normalized depth runs from 0 (pial surface)
to 1 (white matter boundary).

Upper-layer SST interneurons (~L2/3) are of particular interest because:
  - They show strong SCZ GWAS enrichment
  - They have distinct electrophysiology (high sag, HCN channel activity)
  - They express specific SCZ risk genes (CACNA1I, TPH2, BCL11B, HCN1)

This module identifies which SCZ risk genes are preferentially expressed
in upper-layer vs deep-layer SST interneurons.
"""
import numpy as np
import pandas as pd

from ..config import UPPER_LAYER_DEPTH_THRESHOLD, MIDDLE_LAYER_DEPTH_THRESHOLD


def load_depth_data(depth_csv_path):
    """
    Load spatial transcriptomics depth data per supertype.

    The file contains median cortical depth from two platforms (MERFISH, Xenium)
    plus cell counts per platform.

    Parameters
    ----------
    depth_csv_path : str or Path
        Path to median_depth_supertype.csv.

    Returns
    -------
    pd.DataFrame
        Columns: supertype, median_depth_merfish, n_cells_merfish,
                 median_depth_xenium, n_cells_xenium, avg_depth
    """
    df = pd.read_csv(str(depth_csv_path))

    # Compute weighted average depth across platforms
    df["avg_depth"] = (
        df["median_depth_merfish"] * df["n_cells_merfish"]
        + df["median_depth_xenium"] * df["n_cells_xenium"]
    ) / (df["n_cells_merfish"] + df["n_cells_xenium"])

    print(f"  Loaded depth data for {len(df)} supertypes")
    return df


def get_sst_supertypes(depth_df):
    """
    Extract SST-related supertypes from the depth data.

    Parameters
    ----------
    depth_df : pd.DataFrame
        Depth data with 'supertype' column.

    Returns
    -------
    pd.DataFrame
        Subset of depth_df containing only SST supertypes.
    """
    sst_mask = depth_df["supertype"].str.startswith("Sst")
    return depth_df[sst_mask].copy()


def classify_sst_layers(
    depth_df,
    upper_threshold=UPPER_LAYER_DEPTH_THRESHOLD,
    middle_threshold=MIDDLE_LAYER_DEPTH_THRESHOLD,
):
    """
    Classify SST supertypes into upper/middle/deep based on average depth.

    Upper:  avg_depth < upper_threshold (default: 0.35)
    Middle: upper_threshold <= avg_depth < middle_threshold (default: 0.55)
    Deep:   avg_depth >= middle_threshold

    Parameters
    ----------
    depth_df : pd.DataFrame
        Depth data for SST supertypes (from get_sst_supertypes).
    upper_threshold : float
        Depth threshold for upper layer.
    middle_threshold : float
        Depth threshold for middle layer.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with added 'layer_group' column.
    """
    df = depth_df.copy()

    df["layer_group"] = "Deep"
    df.loc[df["avg_depth"] < middle_threshold, "layer_group"] = "Middle"
    df.loc[df["avg_depth"] < upper_threshold, "layer_group"] = "Upper"

    # Print classification summary
    for group in ["Upper", "Middle", "Deep"]:
        types = df[df["layer_group"] == group]["supertype"].tolist()
        if types:
            depths = df[df["layer_group"] == group]["avg_depth"]
            print(f"  {group:6s}: {', '.join(types)}")
            print(f"          depth range: [{depths.min():.3f}, {depths.max():.3f}]")

    return df


def compute_upper_layer_gene_scores(
    specificity, gwas_df, upper_sst_types, deep_sst_types, all_sst_types,
    scz_gene_set=None,
):
    """
    Compute gene-level scores for upper-layer SST specificity × GWAS signal.

    For each gene, computes:
      - upper_sst_spec: sum of specificity across upper SST types
      - deep_sst_spec: sum of specificity across deep SST types
      - total_sst_spec: sum specificity across ALL SST types
      - upper_fraction: upper_sst_spec / total_sst_spec (0–1)
      - upper_gwas_score: upper_sst_spec × gwas_z

    Parameters
    ----------
    specificity : pd.DataFrame
        Full specificity matrix (genes × supertypes).
    gwas_df : pd.DataFrame
        GWAS results with SYMBOL, ZSTAT, P columns.
    upper_sst_types : list of str
        Upper-layer SST supertype names.
    deep_sst_types : list of str
        Deep-layer SST supertype names.
    all_sst_types : list of str
        All SST supertype names.
    scz_gene_set : set, optional
        SCZ GWAS gene set for annotation.

    Returns
    -------
    pd.DataFrame
        Sorted by upper_gwas_score (descending), with columns:
        gene, gwas_zstat, gwas_p, upper_sst_spec, deep_sst_spec,
        total_sst_spec, upper_fraction, is_scz_gwas, upper_gwas_score
    """
    gwas_indexed = gwas_df.set_index("SYMBOL")
    genes = specificity.index

    # Filter to types that exist in the specificity matrix
    upper_types = [t for t in upper_sst_types if t in specificity.columns]
    deep_types = [t for t in deep_sst_types if t in specificity.columns]
    all_types = [t for t in all_sst_types if t in specificity.columns]

    upper_spec = specificity[upper_types].sum(axis=1)
    deep_spec = specificity[deep_types].sum(axis=1)
    total_spec = specificity[all_types].sum(axis=1)

    # Upper fraction (avoid division by zero)
    upper_fraction = np.where(
        total_spec > 0,
        upper_spec / total_spec,
        0
    )

    zstat = gwas_indexed.loc[genes, "ZSTAT"].values
    pvals = gwas_indexed.loc[genes, "P"].values

    result = pd.DataFrame({
        "gene": genes,
        "gwas_zstat": zstat,
        "gwas_p": pvals,
        "upper_sst_spec": upper_spec.values,
        "deep_sst_spec": deep_spec.values,
        "total_sst_spec": total_spec.values,
        "upper_fraction": upper_fraction,
        "upper_gwas_score": upper_spec.values * zstat,
    })

    if scz_gene_set is not None:
        result["is_scz_gwas"] = result["gene"].isin(scz_gene_set)

    result = result.sort_values("upper_gwas_score", ascending=False).reset_index(drop=True)

    # Print top genes
    print(f"\n  Top genes by upper-layer SST × GWAS score:")
    for _, row in result.head(15).iterrows():
        gwas_flag = " *SCZ*" if scz_gene_set and row.get("is_scz_gwas", False) else ""
        print(f"    {row['gene']:15s}  Z={row['gwas_zstat']:6.2f}  "
              f"upper_frac={row['upper_fraction']:.0%}  "
              f"score={row['upper_gwas_score']:.4f}{gwas_flag}")

    return result
