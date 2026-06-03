"""
clusters.py — Siletti cluster mapping and metadata utilities.

Provides mappings between:
  - Conti-style numeric cluster IDs (Cluster0–Cluster460) and named clusters
  - Cluster names and supercluster assignments
  - Cluster name prefixes and broad cell-class labels
"""
import numpy as np

from ..utils import Timer


def build_cluster_number_map(cluster_stats_path):
    """
    Build bidirectional mappings between cluster names and Conti-style numbers.

    The Siletti loom labels clusters CS202210140_1 through CS202210140_461.
    Conti's Cluster0 maps to the cluster with the lowest label number (label 1),
    Cluster1 → label 2, ..., Cluster460 → label 461.

    Parameters
    ----------
    cluster_stats_path : str or Path
        Path to siletti_cluster_level_stats.npz.

    Returns
    -------
    name_to_num : dict
        Cluster name → integer index (0-based).
    num_to_name : dict
        Integer index → cluster name.
    """
    d = np.load(str(cluster_stats_path), allow_pickle=True)
    names = d["cluster_names"]
    labels = d["cluster_labels"]

    label_nums = np.array([int(l.split("_")[-1]) for l in labels])
    sorted_idx = np.argsort(label_nums)
    sorted_names = names[sorted_idx]

    name_to_num = {sorted_names[i]: i for i in range(len(sorted_names))}
    num_to_name = {i: sorted_names[i] for i in range(len(sorted_names))}

    return name_to_num, num_to_name


def load_supercluster_info(cluster_stats_path):
    """
    Load supercluster assignments for all Siletti clusters.

    Parameters
    ----------
    cluster_stats_path : str or Path
        Path to siletti_cluster_level_stats.npz.

    Returns
    -------
    dict
        Cluster name → supercluster name.
    """
    d = np.load(str(cluster_stats_path), allow_pickle=True)
    return dict(zip(d["cluster_names"], d["cluster_superclusters"]))


# Mapping from Siletti cluster name prefix to broad cell class.
# Used for coloring and grouping in figures.
SILETTI_CLASS_MAP = {
    # GABAergic
    "MGE": "GABAergic (MGE)",
    "LLC": "GABAergic (MGE)",
    "Thex": "GABAergic (MGE)",
    "CGE": "GABAergic (CGE)",
    "Splat": "GABAergic (mixed/other)",
    "Midi": "GABAergic (mixed/other)",
    "CBI": "GABAergic (mixed/other)",
    "LRL": "GABAergic (mixed/other)",
    "Mmb": "GABAergic (mixed/other)",
    # Glutamatergic
    "Amex": "Glutamatergic",
    "Misc": "Glutamatergic",
    "DLIT": "Glutamatergic",
    "DLCT6b": "Glutamatergic",
    "DLNP": "Glutamatergic",
    "L5ET": "Glutamatergic",
    "ULIT": "Glutamatergic",
    "DG": "Glutamatergic",
    "CA13": "Glutamatergic",
    "CA4": "Glutamatergic",
    "MSN": "Glutamatergic",
    "EMSN": "Glutamatergic",
    "URL": "Glutamatergic",
    # Non-neuronal
    "Astro": "Non-neuronal",
    "Oligo": "Non-neuronal",
    "OPC": "Non-neuronal",
    "COP": "Non-neuronal",
    "Mgl": "Non-neuronal",
    "Epen": "Non-neuronal",
    "Fbl": "Non-neuronal",
    "Vend": "Non-neuronal",
    "Vsmc": "Non-neuronal",
    "Per": "Non-neuronal",
    "Bgl": "Non-neuronal",
    "Chrp": "Non-neuronal",
    "Bcell": "Non-neuronal",
    "Tcell": "Non-neuronal",
    "Mono": "Non-neuronal",
    "Nkcell": "Non-neuronal",
}


def classify_siletti_cluster(cluster_name):
    """
    Classify a Siletti cluster into a broad cell class based on its name prefix.

    Parameters
    ----------
    cluster_name : str
        Cluster name (e.g., "MGE_240", "Astro_52").

    Returns
    -------
    str
        Broad class label (e.g., "GABAergic (MGE)", "Glutamatergic", "Non-neuronal").
    """
    prefix = cluster_name.split("_")[0]
    return SILETTI_CLASS_MAP.get(prefix, "Other")
