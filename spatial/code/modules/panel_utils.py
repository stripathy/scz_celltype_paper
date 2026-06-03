"""
Xenium panel gene list loading and detection efficiency lookup utilities.

Used by: hierarchical_probe_selection.py, supertype_markers_panel_overlap.py,
         merscope_panel_assessment.py, nsforest_supertype_markers.py
"""

import os
import pandas as pd
import numpy as np


# Default panel paths (~/Downloads — these are 10x metadata CSVs)
_DEFAULT_PANEL_5K = os.path.expanduser(
    "~/Downloads/XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv")
_DEFAULT_PANEL_V1 = os.path.expanduser(
    "~/Downloads/Xenium_hBrain_v1_metadata.csv")


def load_xenium_panels(panel_5k_path=None, panel_v1_path=None):
    """Load Xenium 5K Prime and v1 Brain gene sets from metadata CSVs.

    Parameters
    ----------
    panel_5k_path : str, optional
        Path to Xenium 5K Prime panel metadata CSV.
    panel_v1_path : str, optional
        Path to Xenium v1 Brain panel metadata CSV.

    Returns
    -------
    dict : {'xenium_5k': set of gene names, 'xenium_v1': set of gene names}
    """
    path_5k = panel_5k_path or _DEFAULT_PANEL_5K
    path_v1 = panel_v1_path or _DEFAULT_PANEL_V1

    panel_5k = pd.read_csv(path_5k)
    panel_5k_genes = set(panel_5k["gene_name"].str.strip())

    panel_v1 = pd.read_csv(path_v1)
    panel_v1_genes = set(panel_v1["Genes"].str.strip())

    print(f"  Loaded panels: 5K Prime = {len(panel_5k_genes)} genes, "
          f"v1 Brain = {len(panel_v1_genes)} genes")

    return {"xenium_5k": panel_5k_genes, "xenium_v1": panel_v1_genes}


def load_detection_efficiency(path):
    """Load per-gene snRNAseq vs spatial detection efficiency table.

    Expects a CSV with columns: gene, detection_efficiency (or frac_snrna, frac_spatial).

    Parameters
    ----------
    path : str
        Path to detection efficiency CSV (e.g., snrnaseq_vs_merscope4k_detection.csv).

    Returns
    -------
    dict : gene_name -> detection_efficiency (float)
    float : median detection efficiency (for fallback)
    """
    df = pd.read_csv(path)

    # Handle different column name conventions
    if "detection_efficiency" in df.columns:
        eff_col = "detection_efficiency"
    elif "frac_spatial" in df.columns and "frac_snrna" in df.columns:
        df["detection_efficiency"] = df["frac_spatial"] / (df["frac_snrna"] + 1e-10)
        eff_col = "detection_efficiency"
    else:
        raise ValueError(f"Cannot find detection_efficiency column in {path}. "
                         f"Available columns: {list(df.columns)}")

    gene_col = "gene" if "gene" in df.columns else df.columns[0]
    lookup = dict(zip(df[gene_col], df[eff_col]))
    median_eff = float(np.median(list(lookup.values())))

    print(f"  Loaded detection efficiency for {len(lookup)} genes "
          f"(median = {median_eff:.3f})")

    return lookup, median_eff


def load_spatial_validation(xenium_corr_path=None, merfish_corr_path=None,
                            r_threshold=0.7):
    """Load cross-platform gene correlations for spatial validation.

    Genes with high cross-platform Pearson r (> threshold) on Xenium or MERFISH
    are considered spatially validated.

    Parameters
    ----------
    xenium_corr_path : str, optional
        Path to snrnaseq_vs_xenium_gene_corr.csv.
    merfish_corr_path : str, optional
        Path to snrnaseq_vs_merfish_gene_corr.csv.
    r_threshold : float
        Minimum Pearson r to consider a gene validated (default 0.7).

    Returns
    -------
    set : gene names that pass validation on at least one platform
    dict : gene_name -> best Pearson r across platforms
    """
    validated_genes = set()
    best_r = {}

    for path, platform in [(xenium_corr_path, "Xenium"),
                           (merfish_corr_path, "MERFISH")]:
        if path is None or not os.path.exists(path):
            continue
        df = pd.read_csv(path, index_col=0)
        for gene in df.index:
            r = df.loc[gene, "pearson_r"]
            if r >= r_threshold:
                validated_genes.add(gene)
            if gene not in best_r or r > best_r[gene]:
                best_r[gene] = r

    if validated_genes:
        print(f"  Spatially validated genes (r >= {r_threshold}): {len(validated_genes)}")

    return validated_genes, best_r


def load_gene_quality(path):
    """Load gene properties vs cross-platform correlation table.

    Parameters
    ----------
    path : str
        Path to gene_properties_vs_correlation.csv.

    Returns
    -------
    DataFrame with gene quality metrics (pearson_r, biotype, CV, entropy, etc.)
    """
    df = pd.read_csv(path, index_col=0)
    print(f"  Loaded gene quality data for {len(df)} genes")
    return df
