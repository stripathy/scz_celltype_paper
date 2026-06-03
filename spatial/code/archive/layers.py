"""
Cortical layer segmentation from spatial cell type assignments.

Uses 2D density maps of layer-specific excitatory neuron types to
estimate dominant cortical layers and extract boundary contours.
Handles curved/folded tissue sections.
"""

import numpy as np
from scipy.ndimage import gaussian_filter
from skimage.measure import find_contours


# Default mapping from cortical layers to excitatory neuron subclasses
LAYER_SUBCLASS_MAP = {
    "L2/3": ["L2/3 IT"],
    "L4": ["L4 IT"],
    "L5": ["L5 IT", "L5 ET", "L5/6 NP"],
    "L6": ["L6 IT", "L6 IT Car3", "L6 CT"],
    "L6b": ["L6b"],
    "WM": ["Oligodendrocyte"],
}

LAYER_NAMES = ["L2/3", "L4", "L5", "L6", "L6b", "WM"]

LAYER_COLORS = {
    "L2/3": "#00FFFF",
    "L4": "#00CC00",
    "L5": "#FFFF00",
    "L6": "#FF8800",
    "L6b": "#FF4400",
    "WM": "#888888",
}

LAYER_CONTOUR_COLORS = {
    "L2/3": "#00DDDD",
    "L4": "#00AA00",
    "L5": "#CCCC00",
    "L6": "#DD7700",
    "L6b": "#DD3300",
    "WM": "#666666",
}


def compute_layer_densities(x, y, subclass_labels, bin_size=50,
                            sigma=3, layer_map=None):
    """
    Compute 2D density fraction maps for each cortical layer.

    For each layer, computes a smoothed histogram of the corresponding
    excitatory cell types, then normalizes by total cell density to
    get the local fraction.

    Parameters
    ----------
    x, y : np.ndarray
        Spatial coordinates of all cells.
    subclass_labels : np.ndarray
        Subclass label for each cell.
    bin_size : float
        Histogram bin size in spatial units (µm).
    sigma : float
        Gaussian smoothing sigma in bins.
    layer_map : dict, optional
        Mapping from layer names to subclass lists.
        Defaults to LAYER_SUBCLASS_MAP.

    Returns
    -------
    dict
        Keys: 'fraction_maps' (dict of layer -> 2D array),
              'dominant' (2D int array, -1 = no cells),
              'x_centers', 'y_centers' (1D arrays of bin centers),
              'x_edges', 'y_edges' (1D arrays of bin edges),
              'total_density' (2D smoothed total density).
    """
    if layer_map is None:
        layer_map = LAYER_SUBCLASS_MAP

    x_edges = np.arange(x.min() - bin_size, x.max() + 2 * bin_size, bin_size)
    y_edges = np.arange(y.min() - bin_size, y.max() + 2 * bin_size, bin_size)
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2

    # Per-layer density
    density_maps = {}
    for layer_name, types in layer_map.items():
        mask = np.isin(subclass_labels, types)
        H, _, _ = np.histogram2d(x[mask], y[mask], bins=[x_edges, y_edges])
        density_maps[layer_name] = gaussian_filter(H.T, sigma=sigma)

    # Total density
    H_total, _, _ = np.histogram2d(x, y, bins=[x_edges, y_edges])
    total_smooth = gaussian_filter(H_total.T, sigma=sigma)
    total_smooth[total_smooth < 0.5] = np.nan

    # Fraction maps
    fraction_maps = {}
    for ln in LAYER_NAMES:
        fraction_maps[ln] = density_maps[ln] / total_smooth

    # Dominant layer
    ny, nx = fraction_maps[LAYER_NAMES[0]].shape
    dominant = np.full((ny, nx), -1, dtype=int)
    max_frac = np.zeros((ny, nx))
    for i, ln in enumerate(LAYER_NAMES):
        fm = np.nan_to_num(fraction_maps[ln], nan=0)
        better = fm > max_frac
        dominant[better] = i
        max_frac[better] = fm[better]
    dominant[np.isnan(total_smooth)] = -1

    return {
        "fraction_maps": fraction_maps,
        "dominant": dominant,
        "x_centers": x_centers,
        "y_centers": y_centers,
        "x_edges": x_edges,
        "y_edges": y_edges,
        "total_density": total_smooth,
    }


def extract_layer_contours(dominant, x_centers, y_centers,
                           smooth_sigma=1.5, min_length=20):
    """
    Extract boundary contours from the dominant layer zone map.

    Creates a binary mask for each layer zone, smooths it slightly,
    and finds the 0.5 contour.

    Parameters
    ----------
    dominant : np.ndarray
        2D integer array of dominant layer indices (-1 = no cells).
    x_centers, y_centers : np.ndarray
        Bin center coordinates.
    smooth_sigma : float
        Gaussian sigma for smoothing binary masks before contouring.
    min_length : int
        Minimum number of contour vertices to keep.

    Returns
    -------
    dict
        {layer_name: list of (cx, cy) contour arrays in spatial coords}.
    """
    contours_by_layer = {}
    for i, ln in enumerate(LAYER_NAMES):
        binary_mask = (dominant == i).astype(float)
        binary_smooth = gaussian_filter(binary_mask, smooth_sigma)
        raw_contours = find_contours(binary_smooth, 0.5)

        layer_contours = []
        for contour in raw_contours:
            if len(contour) >= min_length:
                cy = np.interp(
                    contour[:, 0], np.arange(len(y_centers)), y_centers
                )
                cx = np.interp(
                    contour[:, 1], np.arange(len(x_centers)), x_centers
                )
                layer_contours.append((cx, cy))
        contours_by_layer[ln] = layer_contours

    return contours_by_layer


def segment_layers(adata, bin_size=50, sigma=3, contour_sigma=1.5):
    """
    Full layer segmentation pipeline for a single annotated sample.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated AnnData with 'subclass_label' in .obs and
        spatial coords in .obsm['spatial'].
    bin_size : float
        Histogram bin size in µm.
    sigma : float
        Smoothing sigma for density maps.
    contour_sigma : float
        Smoothing sigma for contour extraction.

    Returns
    -------
    dict
        Contains 'densities' (from compute_layer_densities),
        'contours' (from extract_layer_contours).
    """
    coords = adata.obsm["spatial"]
    x, y = coords[:, 0], coords[:, 1]
    subclass = adata.obs["subclass_label"].values

    densities = compute_layer_densities(x, y, subclass,
                                        bin_size=bin_size, sigma=sigma)
    contours = extract_layer_contours(
        densities["dominant"], densities["x_centers"], densities["y_centers"],
        smooth_sigma=contour_sigma,
    )

    return {"densities": densities, "contours": contours}
