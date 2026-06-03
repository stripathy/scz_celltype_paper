#!/usr/bin/env python3
"""
Curved Cortex Strip Identification.

Uses L1 border cells (from BANKSY) to trace the pia surface as a smooth curve,
then defines cortical strips perpendicular to the local pia tangent at each point.
This handles:
  - Curved cortex (pia surface follows gyral/sulcal curvature)
  - Sulcus-gyrus-sulcus geometry (multiple cortical regions in one section)
  - Variable cortical thickness

Algorithm:
  1. Extract L1 cells → clean with DBSCAN → detect cortical segments
  2. Fit smooth pia curve per segment (arc-length parameterized spline)
  3. Compute local tangent/normal directions
  4. Assign all cortical cells to strips via projection onto pia curve
  5. Score strips for layer completeness, ordering, and purity

Usage:
    python3 -u code/analysis/banksy_05_curved_strips.py --sample Br2421
    python3 -u code/analysis/banksy_05_curved_strips.py --sample Br6437 --strip-width 600
"""

import os
import sys
import time
import argparse
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.collections import LineCollection
from scipy.interpolate import UnivariateSpline
from scipy.spatial import cKDTree
from scipy.stats import kendalltau
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from sklearn.neighbors import kneighbors_graph
from scipy.sparse.csgraph import shortest_path

# Project imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import H5AD_DIR, SAMPLE_TO_DX, LAYER_COLORS, LAYER_ORDER

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "modules"))
from depth_model import LAYER_BINS

# ── Constants ─────────────────────────────────────────────────────────

OUT_DIR = os.path.expanduser("~/Github/SCZ_Xenium/output/banksy")
os.makedirs(OUT_DIR, exist_ok=True)

# L1 point cloud cleaning
DBSCAN_EPS = 200          # μm — neighborhood radius for L1 cell clustering
DBSCAN_MIN_SAMPLES = 5    # minimum L1 neighbors within eps
MIN_SEGMENT_L1_CELLS = 20 # minimum L1 cells for a valid cortical segment (low for discontinuous pia)

# Pia curve fitting
PIA_BIN_SPACING = 125      # μm — window spacing for centroid computation
PIA_SPLINE_SMOOTH = None   # auto-determined from data
PIA_SAMPLE_SPACING = 10    # μm — spacing for dense spline sampling

# Strip parameters
STRIP_WIDTH_ALONG_PIA = 600  # μm — strip width measured along pia curve
STRIP_MAX_DEPTH_UM = 3000    # μm — maximum perpendicular depth from pia
MIN_CELLS_PER_LAYER = 15     # minimum cells per layer in a strip

# Fold detection (for splitting pia at gyral crowns / sulcal bottoms)
FOLD_ANGLE_THRESH_DEG = 150   # cumulative turning angle to trigger a split
FOLD_MIN_BANK_LENGTH = 2000   # μm — minimum pia arc-length for a bank to be valid

# Layer definitions
REQUIRED_LAYERS = ["L1", "L2/3", "L4", "L5", "L6"]
LAYER_DEPTH_ORDER = {"L1": 0, "L2/3": 1, "L4": 2, "L5": 3, "L6": 4}


# ── Step 1: Extract and Clean L1 Point Cloud ─────────────────────────

# Direction-aware segmentation parameters
L1_OUTWARD_SAMPLE_DIST = 300   # μm — sample point for outward density
L1_OUTWARD_DENSITY_RADIUS = 200  # μm — radius for density measurement
L1_LOW_DENSITY_THRESH = 30     # cells — below this = "facing nothing"
L1_L23_SEARCH_RADIUS = 300     # μm — radius for finding nearby L2/3
L1_DIRECTION_SCALE = 500       # μm — scale for directional clustering


def extract_l1_points(adata):
    """Extract L1 cell coordinates and detect distinct pial surfaces.

    Uses a direction-aware approach: each L1 cell's "outward direction"
    (away from L2/3, toward empty space) determines which pia surface it
    belongs to. This correctly splits L1 bands that serve multiple cortices
    (e.g., dual cortex with L2/3 → L1 → gap → L1 → L2/3).

    Algorithm:
    1. DBSCAN for noise removal
    2. For each L1 cell, compute direction toward nearest L2/3 (inward)
    3. Verify outward direction faces low-density region (true pia border)
    4. Cluster in (x, y, outward_x*scale, outward_y*scale) space
       → L1 cells facing the same direction = same pia surface
       → L1 cells at a fold facing opposite directions = different surfaces

    Returns
    -------
    l1_clean : ndarray (M, 2) — cleaned L1 cell coordinates
    l1_segment_labels : ndarray (M,) — segment ID for each L1 cell
    l1_clean_idx : ndarray (M,) — indices into original adata for cleaned L1 cells
    """
    coords = adata.obsm["spatial"]

    # Get L1 mask
    if "banksy_is_l1" in adata.obs.columns:
        is_l1 = adata.obs["banksy_is_l1"].values.astype(bool)
    else:
        is_l1 = adata.obs["banksy_domain"].values == "Meningeal"

    if is_l1.sum() == 0:
        # Fallback: use depth-based L1
        pred_depth = pd.to_numeric(adata.obs["predicted_norm_depth"], errors="coerce").values
        is_l1 = pred_depth < 0.10
        print(f"  WARNING: No banksy_is_l1 cells, using depth<0.10 fallback ({is_l1.sum()} cells)")

    # Also identify L2/3 cells for adjacency detection
    layers = adata.obs["layer"].values.astype(str)
    is_l23 = layers == "L2/3"

    l1_coords = coords[is_l1]
    l1_idx = np.where(is_l1)[0]
    n_l1 = len(l1_coords)
    print(f"  L1 cells: {n_l1}, L2/3 cells: {is_l23.sum()}")

    if n_l1 < MIN_SEGMENT_L1_CELLS:
        print(f"  WARNING: Only {n_l1} L1 cells — too few for pia curve fitting")
        return np.empty((0, 2)), np.array([]), np.array([])

    # Step A: Find L1 cells adjacent to L2/3
    l23_coords = coords[is_l23]
    if len(l23_coords) > 0:
        l23_tree = cKDTree(l23_coords)
        dist_to_l23, _ = l23_tree.query(l1_coords)
        l1_touches_l23 = dist_to_l23 < 200
        n_touching = l1_touches_l23.sum()
        print(f"  L1 cells adjacent to L2/3 (within 200μm): {n_touching} / {n_l1}")
    else:
        l1_touches_l23 = np.ones(n_l1, dtype=bool)
        print(f"  WARNING: No L2/3 cells found, using all L1 cells")

    # Step B: DBSCAN for noise removal
    db = DBSCAN(eps=DBSCAN_EPS, min_samples=DBSCAN_MIN_SAMPLES).fit(l1_coords)
    core_mask = db.labels_ >= 0
    n_noise = (~core_mask).sum()
    print(f"  DBSCAN noise removal: {core_mask.sum()} core L1 cells, {n_noise} noise")

    # Step C: Compute outward direction for each L1 cell
    # "Outward" = away from L2/3 = toward empty space / pia surface
    valid_for_dir = core_mask & l1_touches_l23
    outward_dirs = np.zeros((n_l1, 2))
    has_direction = np.zeros(n_l1, dtype=bool)

    if len(l23_coords) > 0 and valid_for_dir.sum() > 0:
        for i in range(n_l1):
            if not valid_for_dir[i]:
                continue
            nearby = l23_tree.query_ball_point(l1_coords[i], L1_L23_SEARCH_RADIUS)
            if len(nearby) >= 3:
                centroid = l23_coords[nearby].mean(axis=0)
                inward = centroid - l1_coords[i]
                norm = np.linalg.norm(inward)
                if norm > 1:
                    outward_dirs[i] = -inward / norm  # outward = away from L2/3
                    has_direction[i] = True

    n_with_dir = has_direction.sum()
    print(f"  L1 cells with outward direction: {n_with_dir} / {n_l1}")

    # Step D: Direction-aware clustering
    # Cluster in (x, y, outward_x * scale, outward_y * scale) space
    # This splits L1 bands where outward direction changes sharply (dual cortex)
    # while keeping curved single-cortex bands together (gradual rotation)
    usable = core_mask & has_direction
    if usable.sum() >= MIN_SEGMENT_L1_CELLS:
        usable_coords = l1_coords[usable]
        usable_dirs = outward_dirs[usable]

        # Build feature matrix: (x, y, outward_x * scale, outward_y * scale)
        features = np.column_stack([
            usable_coords,
            usable_dirs * L1_DIRECTION_SCALE
        ])

        # DBSCAN on combined spatial+directional features
        # Using generous eps to keep curved cortex L1 bands together.
        # Spatial gaps are handled by post-processing (split_segments_at_gaps).
        db_dir = DBSCAN(eps=DBSCAN_EPS * 1.5, min_samples=3).fit(features)

        full_labels = np.full(n_l1, -1, dtype=int)
        full_labels[usable] = db_dir.labels_

        # For core L1 cells without direction, assign to nearest usable cell's segment
        core_no_dir = core_mask & ~has_direction
        if core_no_dir.sum() > 0 and usable.sum() > 0:
            usable_tree = cKDTree(usable_coords)
            _, nearest_usable = usable_tree.query(l1_coords[core_no_dir])
            full_labels[core_no_dir] = db_dir.labels_[nearest_usable]
    else:
        # Fall back to spatial-only DBSCAN
        full_labels = db.labels_.copy()
        print(f"  WARNING: Few L1 cells with direction, using spatial-only segmentation")

    # Filter to valid (non-noise) cells
    valid = full_labels >= 0
    l1_clean = l1_coords[valid]
    l1_segment_labels = full_labels[valid]
    l1_clean_idx = l1_idx[valid]

    # Step E: Split segments at spatial gaps
    # The generous DBSCAN eps can merge L1 clusters across real tissue gaps.
    # Check each segment for internal spatial gaps (disconnected sub-clusters)
    # and split when both sides of a gap are confirmed pia borders
    # (each side has L2/3 nearby = cortex on that side).
    l1_segment_labels = _split_segments_at_gaps(
        l1_clean, l1_segment_labels, l23_coords
    )

    # Step F: Split "outer" vs "inner" L1 cells within each segment
    # Detects the pattern: nothing → L1 → L2/3 → L1 → L2/3 (stacked cortices)
    # "Outer" L1 directly faces empty space (low density in outward direction)
    # "Inner" L1 has L2/3 between it and the tissue edge (another cortex above)
    # These are different pia surfaces even though they're spatially adjacent.
    l1_segment_labels = _split_inner_outer_l1(
        l1_clean, l1_segment_labels, outward_dirs[valid], has_direction[valid],
        coords, l23_coords
    )

    # Renumber segments to be contiguous (0, 1, 2, ...)
    unique_segs = np.unique(l1_segment_labels)
    seg_map = {old: new for new, old in enumerate(unique_segs)}
    l1_segment_labels = np.array([seg_map[s] for s in l1_segment_labels])

    # Report segments
    segments = np.unique(l1_segment_labels)
    print(f"  Detected {len(segments)} pial surface(s) (direction-aware + gap-split):")
    for seg_id in segments:
        seg_mask = l1_segment_labels == seg_id
        n_seg = seg_mask.sum()
        seg_dirs = outward_dirs[valid][seg_mask]
        seg_has_dir = has_direction[valid][seg_mask]
        if seg_has_dir.sum() > 0:
            mean_dir = seg_dirs[seg_has_dir].mean(axis=0)
            angle = np.degrees(np.arctan2(mean_dir[1], mean_dir[0]))
            print(f"    Surface {seg_id}: {n_seg} L1 cells, "
                  f"mean outward angle: {angle:.0f}°")
        else:
            print(f"    Surface {seg_id}: {n_seg} L1 cells (no direction info)")

    return l1_clean, l1_segment_labels, l1_clean_idx


def _split_segments_at_gaps(l1_clean, l1_segment_labels, l23_coords,
                            gap_threshold=200, min_fragment_size=20):
    """Split L1 segments at spatial gaps (tissue boundaries).

    For each segment, uses spatial DBSCAN (eps=gap_threshold) to find
    spatially disconnected sub-clusters. If a sub-cluster has enough cells
    AND has L2/3 nearby (confirming it's a real pia border, not noise),
    it becomes its own segment.

    This handles the "nothing → L1 → L2/3" pattern: two L1 regions separated
    by a gap each form a distinct pia surface if they both border cortex.

    Parameters
    ----------
    l1_clean : ndarray (M, 2) — cleaned L1 coordinates
    l1_segment_labels : ndarray (M,) — segment IDs
    l23_coords : ndarray (K, 2) — L2/3 cell coordinates
    gap_threshold : float — spatial distance to consider as a gap (μm)
    min_fragment_size : int — minimum L1 cells to form a valid segment
    """
    if len(l23_coords) == 0:
        return l1_segment_labels

    l23_tree = cKDTree(l23_coords)
    new_labels = l1_segment_labels.copy()
    next_label = l1_segment_labels.max() + 1

    segments = np.unique(l1_segment_labels)
    n_splits = 0

    for seg_id in segments:
        seg_mask = l1_segment_labels == seg_id
        n_seg = seg_mask.sum()
        if n_seg < 2 * min_fragment_size:
            continue  # too small to split

        seg_coords = l1_clean[seg_mask]
        seg_indices = np.where(seg_mask)[0]

        # DBSCAN with strict spatial eps to find disconnected sub-clusters
        db_gap = DBSCAN(eps=gap_threshold, min_samples=3).fit(seg_coords)

        unique_clusters = set(db_gap.labels_)
        unique_clusters.discard(-1)  # remove noise

        if len(unique_clusters) <= 1:
            continue  # no spatial gaps → nothing to split

        # Verify each sub-cluster is a real pia border (has L2/3 nearby)
        valid_clusters = []
        for cl in sorted(unique_clusters):
            cl_mask = db_gap.labels_ == cl
            cl_coords = seg_coords[cl_mask]
            n_cl = cl_mask.sum()

            if n_cl < min_fragment_size:
                continue

            # Check if L2/3 exists nearby (within 300μm of any L1 cell in this cluster)
            dists, _ = l23_tree.query(cl_coords)
            n_touching_l23 = (dists < 300).sum()
            frac_touching = n_touching_l23 / n_cl

            if frac_touching > 0.5:
                valid_clusters.append((cl, n_cl, frac_touching))

        if len(valid_clusters) <= 1:
            continue  # only one valid sub-cluster → no real split needed

        # Split: keep the largest cluster with the original segment ID,
        # assign others new segment IDs
        valid_clusters.sort(key=lambda x: -x[1])  # largest first

        for i, (cl, n_cl, frac) in enumerate(valid_clusters):
            cl_mask = db_gap.labels_ == cl
            cl_indices = seg_indices[cl_mask]

            if i == 0:
                # Keep original label for largest
                pass
            else:
                # Assign new label
                new_labels[cl_indices] = next_label
                next_label += 1
                n_splits += 1

    if n_splits > 0:
        print(f"  Gap splitting: {n_splits} new segment(s) created from spatial gaps")

    return new_labels


def _split_inner_outer_l1(l1_clean, l1_segment_labels, outward_dirs, has_direction,
                           all_coords, l23_coords, min_fragment_size=30,
                           thickness_threshold=400):
    """Split thick L1 bands that contain stacked cortices.

    Distinguishes TRUE stacked cortex from curved cortex using bimodality:
    - Stacked cortex: L1 PC2 distribution is BIMODAL (two bands with a gap)
      → L1_band_A | L2/3 gap | L1_band_B
    - Curved cortex: L1 PC2 distribution is continuous (cells fill the bend)
      → appears thick in PCA but L1 density never drops to zero

    Algorithm:
    1. Compute PCA thickness — skip if < threshold (normal thin L1)
    2. Check L1 PC2 distribution for BIMODALITY (gap in L1 density)
    3. If bimodal, verify L2/3 fills the gap between L1 peaks
    4. Split at the gap center

    Parameters
    ----------
    thickness_threshold : float — minimum L1 band thickness (μm) to trigger check
    """
    if len(l23_coords) == 0:
        return l1_segment_labels

    l23_tree = cKDTree(l23_coords)
    new_labels = l1_segment_labels.copy()
    next_label = l1_segment_labels.max() + 1
    n_splits = 0

    segments = np.unique(l1_segment_labels)

    for seg_id in segments:
        seg_mask = l1_segment_labels == seg_id
        n_seg = seg_mask.sum()
        if n_seg < 2 * min_fragment_size:
            continue

        seg_coords = l1_clean[seg_mask]
        seg_indices = np.where(seg_mask)[0]

        # Step 1: Compute L1 band thickness via PCA
        pca = PCA(n_components=2).fit(seg_coords)
        pc2_dir = pca.components_[1]  # perpendicular to band
        pc1_dir = pca.components_[0]  # along band
        center = pca.mean_

        # Project L1 cells onto PC2 (perpendicular axis)
        pc2_proj = (seg_coords - center) @ pc2_dir

        # Band thickness = range of PC2 values (5th to 95th percentile)
        p5, p95 = np.percentile(pc2_proj, [5, 95])
        thickness = p95 - p5

        if thickness < thickness_threshold:
            continue  # Normal thin L1 band

        # Step 2: Check for BIMODALITY in L1 PC2 distribution
        # True stacked cortex has a gap (valley) in L1 density along PC2.
        # Curved cortex has continuous L1 density even though PCA spread is large.
        n_bins = max(30, int(thickness / 50))  # ~50μm per bin
        bin_edges = np.linspace(p5, p95, n_bins + 1)
        bin_counts, _ = np.histogram(pc2_proj, bins=bin_edges)

        # Find the gap: consecutive bins where L1 count is very low
        # "Low" = less than 15% of the mean count of non-empty bins
        nonempty_mean = bin_counts[bin_counts > 0].mean() if (bin_counts > 0).any() else 1
        low_threshold = max(2, nonempty_mean * 0.15)
        is_low = bin_counts < low_threshold

        # Find the longest run of consecutive low bins
        # This represents the gap between the two L1 bands
        best_gap_start = -1
        best_gap_len = 0
        current_start = -1
        current_len = 0
        for i in range(len(is_low)):
            if is_low[i]:
                if current_start < 0:
                    current_start = i
                current_len += 1
            else:
                if current_len > best_gap_len:
                    best_gap_start = current_start
                    best_gap_len = current_len
                current_start = -1
                current_len = 0
        if current_len > best_gap_len:
            best_gap_start = current_start
            best_gap_len = current_len

        # Gap must be at least 3 bins (~150μm) and not at the edges
        gap_width_um = best_gap_len * (p95 - p5) / n_bins
        gap_is_internal = (best_gap_start > 1 and
                           best_gap_start + best_gap_len < n_bins - 1)

        if best_gap_len < 3 or gap_width_um < 100 or not gap_is_internal:
            # No clear gap → continuous L1 band → curved cortex, not stacked
            if thickness > 500:
                print(f"  Segment {seg_id}: thick ({thickness:.0f}μm) but no L1 gap "
                      f"(best gap: {best_gap_len} bins, {gap_width_um:.0f}μm) → curved cortex")
            continue

        # Gap center in PC2 coordinates
        gap_center_bin = best_gap_start + best_gap_len / 2
        gap_center_pc2 = p5 + gap_center_bin * (p95 - p5) / n_bins

        # Step 3: Check for L2/3 cells filling the gap
        seg_pc1 = (seg_coords - center) @ pc1_dir
        pc1_min, pc1_max = seg_pc1.min() - 200, seg_pc1.max() + 200
        seg_span = seg_pc1.max() - seg_pc1.min()

        # Find L2/3 cells near the L1 band
        seg_center = seg_coords.mean(axis=0)
        l23_near = l23_tree.query_ball_point(seg_center,
                                              max(1000, thickness * 2))

        if len(l23_near) < 10:
            continue

        l23_near_coords = l23_coords[l23_near]
        l23_pc2 = (l23_near_coords - center) @ pc2_dir
        l23_pc1 = (l23_near_coords - center) @ pc1_dir

        # Filter L2/3 to those within the band's PC1 extent
        in_band_pc1 = (l23_pc1 >= pc1_min) & (l23_pc1 <= pc1_max)

        # Count L2/3 cells in the gap region
        gap_lo = bin_edges[best_gap_start]
        gap_hi = bin_edges[best_gap_start + best_gap_len]
        l23_in_gap = in_band_pc1 & (l23_pc2 >= gap_lo) & (l23_pc2 <= gap_hi)
        n_l23_gap = l23_in_gap.sum()

        if n_l23_gap < 20:
            print(f"  Segment {seg_id}: L1 gap found ({gap_width_um:.0f}μm) but "
                  f"only {n_l23_gap} L2/3 cells in gap → not stacked")
            continue

        # Verify L2/3 in gap spans enough of the band length
        l23_gap_pc1 = l23_pc1[l23_in_gap]
        l23_gap_span = l23_gap_pc1.max() - l23_gap_pc1.min() if n_l23_gap > 1 else 0

        if l23_gap_span < seg_span * 0.3:
            print(f"  Segment {seg_id}: L1 gap found but L2/3 span too small "
                  f"({l23_gap_span:.0f}/{seg_span:.0f}μm) → not stacked")
            continue

        # Count L1 on each side of the gap
        group_a = pc2_proj < gap_center_pc2
        group_b = pc2_proj >= gap_center_pc2
        n_a, n_b = group_a.sum(), group_b.sum()

        if n_a < min_fragment_size or n_b < min_fragment_size:
            continue

        print(f"  Stacked cortex in segment {seg_id}: "
              f"thickness={thickness:.0f}μm, "
              f"L1 gap={gap_width_um:.0f}μm wide ({best_gap_len} bins), "
              f"{n_l23_gap} L2/3 cells in gap, "
              f"L2/3 gap span={l23_gap_span:.0f}/{seg_span:.0f}μm, "
              f"L1 groups: {n_a} vs {n_b}")

        # Step 4: Split at the gap center
        if n_a < n_b:
            for i in range(n_seg):
                if group_a[i]:
                    new_labels[seg_indices[i]] = next_label
        else:
            for i in range(n_seg):
                if group_b[i]:
                    new_labels[seg_indices[i]] = next_label

        next_label += 1
        n_splits += 1

    if n_splits > 0:
        print(f"  Stacked cortex splitting: {n_splits} segment(s) split")

    return new_labels


# ── Step 2: Fit Pia Curve Per Segment ────────────────────────────────

def fit_pia_curve(seg_coords, segment_id=0):
    """Fit a smooth pia curve to one cortical segment's L1 cells.

    Uses PCA projection for elongated segments or arc-length parameterization
    for curved (U-shaped) segments.

    Returns
    -------
    spline_x, spline_y : UnivariateSpline — parameterized pia curve
    t_range : tuple (t_min, t_max) — parameter domain
    method : str — "pca" or "arc_length"
    """
    n = len(seg_coords)
    if n < 10:
        return None, None, (0, 0), "too_few"

    # Decide PCA vs arc-length based on variance ratio
    pca = PCA(n_components=2)
    pca.fit(seg_coords)
    var_ratio = pca.explained_variance_ratio_[1] / max(pca.explained_variance_ratio_[0], 1e-10)

    print(f"    Segment {segment_id}: PC2/PC1 variance ratio = {var_ratio:.3f}")

    if var_ratio > 0.3:
        # Curved segment — use arc-length parameterization
        method = "arc_length"
        spline_x, spline_y, t_range = _fit_arc_length_spline(seg_coords)
    else:
        # Elongated segment — PCA projection works
        method = "pca"
        spline_x, spline_y, t_range = _fit_pca_spline(seg_coords, pca)

    print(f"    Method: {method}, curve length: {t_range[1] - t_range[0]:.0f} μm")
    return spline_x, spline_y, t_range, method


def _fit_pca_spline(seg_coords, pca):
    """Fit smoothing spline using PCA projection along dominant axis."""
    pc1_dir = pca.components_[0]
    center = pca.mean_

    # Project onto PC1
    proj = (seg_coords - center) @ pc1_dir

    # Sort and bin into windows
    n_bins = max(10, int((proj.max() - proj.min()) / PIA_BIN_SPACING))
    bin_edges = np.linspace(proj.min(), proj.max(), n_bins + 1)

    centroids_t = []
    centroids_x = []
    centroids_y = []

    for i in range(n_bins):
        in_bin = (proj >= bin_edges[i]) & (proj < bin_edges[i + 1])
        if in_bin.sum() >= 3:
            centroids_t.append((bin_edges[i] + bin_edges[i + 1]) / 2)
            centroids_x.append(np.mean(seg_coords[in_bin, 0]))
            centroids_y.append(np.mean(seg_coords[in_bin, 1]))

    if len(centroids_t) < 4:
        return None, None, (0, 0)

    t = np.array(centroids_t)
    cx = np.array(centroids_x)
    cy = np.array(centroids_y)

    # Convert to arc-length parameter for more natural parameterization
    dists = np.sqrt(np.diff(cx)**2 + np.diff(cy)**2)
    arc = np.concatenate([[0], np.cumsum(dists)])

    # Smoothing: s = number of points × variance gives moderate smoothing
    s_val = len(t) * 500
    spline_x = UnivariateSpline(arc, cx, s=s_val, k=min(3, len(t) - 1))
    spline_y = UnivariateSpline(arc, cy, s=s_val, k=min(3, len(t) - 1))

    return spline_x, spline_y, (arc[0], arc[-1])


def _fit_arc_length_spline(seg_coords):
    """Fit smoothing spline using arc-length parameterization via KNN graph.

    For U-shaped or C-shaped L1 point clouds where PCA fails.
    """
    n = len(seg_coords)
    k = min(10, n - 1)

    # Build KNN graph weighted by Euclidean distance
    G = kneighbors_graph(seg_coords, n_neighbors=k, mode='distance',
                         include_self=False)
    G = (G + G.T) / 2  # symmetrize

    # Find endpoints: the two cells with maximum geodesic distance
    # Use shortest_path from a random starting node first, then from the farthest node
    dist_from_0 = shortest_path(G, directed=False, indices=0)
    i_start = np.argmax(dist_from_0)
    dist_from_start = shortest_path(G, directed=False, indices=i_start)
    i_end = np.argmax(dist_from_start)

    # Get shortest path from start to end
    dist_mat, predecessors = shortest_path(G, directed=False, indices=i_start,
                                            return_predecessors=True)

    # Reconstruct path
    path = [i_end]
    visited = {i_end}
    while path[-1] != i_start:
        pred = predecessors[path[-1]]
        if pred < 0 or pred in visited:
            break  # disconnected graph or cycle
        path.append(pred)
        visited.add(pred)
    path = path[::-1]

    if len(path) < 4:
        # Fallback to PCA if path reconstruction failed
        pca = PCA(n_components=2).fit(seg_coords)
        return _fit_pca_spline(seg_coords, pca)

    # Compute arc-length along path
    path_coords = seg_coords[path]
    diffs = np.diff(path_coords, axis=0)
    segment_lengths = np.sqrt((diffs**2).sum(axis=1))
    arc_length = np.concatenate([[0], np.cumsum(segment_lengths)])

    # Subsample at regular intervals for spline fitting
    target_spacing = PIA_BIN_SPACING
    n_sample = max(10, int(arc_length[-1] / target_spacing))
    t_regular = np.linspace(0, arc_length[-1], n_sample)
    x_interp = np.interp(t_regular, arc_length, path_coords[:, 0])
    y_interp = np.interp(t_regular, arc_length, path_coords[:, 1])

    # Smooth spline
    s_val = n_sample * 500
    spline_x = UnivariateSpline(t_regular, x_interp, s=s_val, k=min(3, n_sample - 1))
    spline_y = UnivariateSpline(t_regular, y_interp, s=s_val, k=min(3, n_sample - 1))

    return spline_x, spline_y, (0, arc_length[-1])


# ── Step 3: Compute Tangents and Inward Normals ─────────────────────

def compute_normals(spline_x, spline_y, t_range, cortical_coords, cortical_depths,
                    n_sample=200):
    """Compute inward-pointing normals along the pia curve.

    Uses a depth-gradient approach: for each pia point, find nearby cortical
    cells, compute the local depth gradient, and orient the normal to align
    with increasing depth (pia→WM direction).

    Returns
    -------
    sample_t : ndarray (S,) — parameter values at sampled points
    sample_xy : ndarray (S, 2) — pia curve coordinates at sampled points
    normals : ndarray (S, 2) — unit inward-pointing normals
    tangents : ndarray (S, 2) — unit tangent vectors
    """
    t_min, t_max = t_range
    sample_t = np.linspace(t_min, t_max, n_sample)

    # Evaluate spline and derivatives
    x = spline_x(sample_t)
    y = spline_y(sample_t)
    dx = spline_x.derivative()(sample_t)
    dy = spline_y.derivative()(sample_t)

    sample_xy = np.column_stack([x, y])

    # Tangent vectors (unit)
    tangent_lengths = np.sqrt(dx**2 + dy**2)
    tangent_lengths = np.where(tangent_lengths > 0, tangent_lengths, 1.0)
    tangents = np.column_stack([dx / tangent_lengths, dy / tangent_lengths])

    # Two candidate normals at each point (perpendicular to tangent)
    normals_a = np.column_stack([-tangents[:, 1], tangents[:, 0]])   # CCW 90deg
    normals_b = np.column_stack([tangents[:, 1], -tangents[:, 0]])    # CW 90deg

    # Determine inward direction using local depth gradient
    # For each pia point, find cortical cells within a wide radius and compute
    # the centroid-weighted-by-depth direction (= local depth gradient)
    cortical_tree = cKDTree(cortical_coords)

    # Test at multiple sample points with depth-gradient approach
    test_indices = np.linspace(0, n_sample - 1, min(30, n_sample), dtype=int)
    votes_a = 0
    votes_b = 0

    for idx in test_indices:
        pia_pt = sample_xy[idx]
        # Find ALL cortical cells within 1500μm (covers most of cortex)
        nearby_idx = cortical_tree.query_ball_point(pia_pt, 1500)
        if len(nearby_idx) < 20:
            continue

        nearby_coords = cortical_coords[nearby_idx]
        nearby_depths = cortical_depths[nearby_idx]

        # Compute depth-weighted centroid relative to pia point
        # Deeper cells pull the centroid in the "inward" direction
        valid = ~np.isnan(nearby_depths)
        if valid.sum() < 20:
            continue

        weights = nearby_depths[valid]  # depth as weight (0 at pia, 1 at WM)
        # Skip if all weights are ~same (no depth gradient)
        if weights.std() < 0.05:
            continue

        vec_to_cells = nearby_coords[valid] - pia_pt
        # Depth-weighted mean direction = direction toward deeper cortex
        grad_dir = np.average(vec_to_cells, axis=0, weights=weights)
        grad_dir = grad_dir / max(np.linalg.norm(grad_dir), 1e-10)

        # Which candidate normal aligns better with the depth gradient?
        dot_a = np.dot(normals_a[idx], grad_dir)
        dot_b = np.dot(normals_b[idx], grad_dir)

        if dot_a > dot_b:
            votes_a += 1
        elif dot_b > dot_a:
            votes_b += 1

    # Select the normal direction that aligns with depth gradient
    if votes_a >= votes_b:
        normals = normals_a.copy()
        print(f"    Normal direction: CCW (votes {votes_a} vs {votes_b})")
    else:
        normals = normals_b.copy()
        print(f"    Normal direction: CW (votes {votes_b} vs {votes_a})")

    # Enforce consistency along the curve (prevent flipping)
    for i in range(1, len(normals)):
        if np.dot(normals[i], normals[i - 1]) < 0:
            normals[i] = -normals[i]

    return sample_t, sample_xy, normals, tangents


# ── Step 3b: Split Pia at Folds (Gyral Crowns / Sulcal Bottoms) ───

def split_pia_at_folds(pia_seg, fold_thresh_deg=FOLD_ANGLE_THRESH_DEG,
                       min_bank_length=FOLD_MIN_BANK_LENGTH):
    """Split a single pia segment into banks at fold points.

    When a pia curve wraps around a gyrus, the cortex on either side of the
    gyral crown represents separate cortical columns. This function detects
    points where the tangent direction has rotated significantly (cumulative
    turning angle exceeds threshold) and splits the pia into sub-segments
    ("banks"), each representing a distinct cortical column.

    Parameters
    ----------
    pia_seg : dict — a pia segment with 'sample_t', 'sample_xy', 'normals',
        'tangents', 'spline_x', 'spline_y', 't_range'
    fold_thresh_deg : float — cumulative turning angle to trigger a split
    min_bank_length : float — minimum arc-length (μm) for a bank to be valid

    Returns
    -------
    banks : list of dict — sub-segments, each with same keys as input plus
        'bank_id' and 'parent_segment_id'
    """
    tangents = pia_seg["tangents"]
    sample_t = pia_seg["sample_t"]
    sample_xy = pia_seg["sample_xy"]
    normals = pia_seg["normals"]
    spline_x = pia_seg["spline_x"]
    spline_y = pia_seg["spline_y"]
    t_range = pia_seg["t_range"]
    seg_id = pia_seg.get("segment_id", 0)

    n = len(tangents)
    if n < 4:
        return [pia_seg]

    # Compute incremental turning angle between consecutive tangent vectors
    # Using atan2 of cross product / dot product for signed angle
    angles = np.zeros(n - 1)
    for i in range(n - 1):
        cross = tangents[i, 0] * tangents[i + 1, 1] - tangents[i, 1] * tangents[i + 1, 0]
        dot = tangents[i, 0] * tangents[i + 1, 0] + tangents[i, 1] * tangents[i + 1, 1]
        angles[i] = np.degrees(np.arctan2(cross, dot))

    # Cumulative turning angle (reset at each split point)
    cumulative = np.cumsum(angles)

    # Find split points where cumulative angle crosses threshold
    fold_thresh = fold_thresh_deg
    split_indices = [0]  # always start with index 0

    last_split_angle = 0.0
    for i in range(len(cumulative)):
        if abs(cumulative[i] - last_split_angle) >= fold_thresh:
            split_indices.append(i + 1)  # +1 because angles are between consecutive points
            last_split_angle = cumulative[i]

    split_indices.append(n)  # always end with last index

    # Remove duplicate / too-close splits
    split_indices = sorted(set(split_indices))

    if len(split_indices) <= 2:
        # No folds detected — return as single bank
        pia_seg_with_bank = dict(pia_seg)
        pia_seg_with_bank["bank_id"] = 0
        pia_seg_with_bank["parent_segment_id"] = seg_id
        return [pia_seg_with_bank]

    # Create sub-segments (banks)
    banks = []
    bank_id = 0
    for bi in range(len(split_indices) - 1):
        i_start = split_indices[bi]
        i_end = split_indices[bi + 1]

        if i_end - i_start < 4:
            continue

        sub_t = sample_t[i_start:i_end]
        sub_xy = sample_xy[i_start:i_end]
        sub_normals = normals[i_start:i_end]
        sub_tangents = tangents[i_start:i_end]

        # Check minimum bank length
        bank_length = sub_t[-1] - sub_t[0]
        if bank_length < min_bank_length:
            continue

        # Create a sub-segment with its own t_range
        sub_t_range = (sub_t[0], sub_t[-1])

        banks.append({
            "segment_id": seg_id,
            "bank_id": bank_id,
            "parent_segment_id": seg_id,
            "spline_x": spline_x,  # same spline, different t_range
            "spline_y": spline_y,
            "t_range": sub_t_range,
            "sample_t": sub_t,
            "sample_xy": sub_xy,
            "normals": sub_normals,
            "tangents": sub_tangents,
            "method": pia_seg.get("method", "unknown"),
        })
        bank_id += 1

    if not banks:
        # All sub-segments too short — return original
        pia_seg_with_bank = dict(pia_seg)
        pia_seg_with_bank["bank_id"] = 0
        pia_seg_with_bank["parent_segment_id"] = seg_id
        return [pia_seg_with_bank]

    return banks


# ── Step 4: Assign Cells to Strips ──────────────────────────────────

def assign_cells_to_strips(coords, domains, pia_segments, strip_width=STRIP_WIDTH_ALONG_PIA,
                           max_depth=STRIP_MAX_DEPTH_UM):
    """Assign cortical cells to curved strips via projection onto pia curves.

    Parameters
    ----------
    coords : ndarray (N, 2) — all cell coordinates
    domains : ndarray (N,) str — domain labels
    pia_segments : list of dict, each with keys:
        'sample_xy', 'sample_t', 'normals', 'tangents', 'spline_x', 'spline_y',
        't_range', and optionally 'bank_id', 'parent_segment_id'

    Returns
    -------
    strip_ids : ndarray (N,) int — strip ID per cell (-1 if not in strip)
    strip_segment : ndarray (N,) int — which pia segment/bank the cell belongs to (-1)
    strip_bank : ndarray (N,) int — bank ID per cell (-1 if not in strip)
    perp_depth_um : ndarray (N,) float — perpendicular distance from pia in μm (NaN)
    pia_arc_pos : ndarray (N,) float — arc-length position along pia curve (NaN)
    """
    n = len(coords)
    strip_ids = np.full(n, -1, dtype=int)
    strip_segment = np.full(n, -1, dtype=int)
    strip_bank = np.full(n, -1, dtype=int)
    perp_depth_um = np.full(n, np.nan)
    pia_arc_pos = np.full(n, np.nan)

    cortical_mask = domains == "Cortical"
    cortical_idx = np.where(cortical_mask)[0]
    cortical_coords = coords[cortical_mask]

    if len(cortical_idx) == 0:
        return strip_ids, strip_segment, perp_depth_um, pia_arc_pos

    global_strip_offset = 0

    for seg_i, seg in enumerate(pia_segments):
        sample_xy = seg["sample_xy"]
        sample_t = seg["sample_t"]
        normals = seg["normals"]
        t_range = seg["t_range"]

        if sample_xy is None or len(sample_xy) < 2:
            continue

        # Densely sample the pia spline for KDTree projection
        spline_x = seg["spline_x"]
        spline_y = seg["spline_y"]
        t_min, t_max = t_range
        n_dense = max(100, int((t_max - t_min) / PIA_SAMPLE_SPACING))
        t_dense = np.linspace(t_min, t_max, n_dense)
        xy_dense = np.column_stack([spline_x(t_dense), spline_y(t_dense)])

        # Precompute normals at dense pia samples by interpolation
        normal_t_samples = seg["sample_t"]
        normals_dense = np.zeros((n_dense, 2))
        for di in range(n_dense):
            ni = np.searchsorted(normal_t_samples, t_dense[di])
            ni = np.clip(ni, 0, len(normals) - 1)
            normals_dense[di] = normals[ni]

        # Build KDTree on dense pia samples
        pia_tree = cKDTree(xy_dense)

        n_strips_seg = max(1, int((t_max - t_min) / strip_width))

        # Normal-aware assignment: for each cortical cell, find the pia point
        # whose normal ray passes closest to the cell. This makes strip
        # boundaries curve with the cortex rather than cutting straight across.
        #
        # At sharp curves (e.g., bottom of a sulcus), the nearest pia point
        # by Euclidean distance may be on the wrong "side" of the curve.
        # Instead, we find K candidates and pick the one where the cell is
        # most aligned with the pia normal direction.
        K_CANDIDATES = 20
        k_actual = min(K_CANDIDATES, len(xy_dense))
        dists_k, idx_k = pia_tree.query(cortical_coords, k=k_actual)

        # For single-candidate case (very short pia), fall back to nearest
        if k_actual == 1:
            nearest_idx = idx_k.ravel()
            distances = dists_k.ravel()
        else:
            # Vectorized normal-aware selection across all cells × K candidates
            N_c = len(cortical_coords)
            # Gather candidate pia points and normals: (N, K, 2)
            cand_xy = xy_dense[idx_k]          # (N, K, 2)
            cand_normals = normals_dense[idx_k] # (N, K, 2)

            # Vector from each candidate pia point to cell: (N, K, 2)
            cell_expanded = cortical_coords[:, np.newaxis, :]  # (N, 1, 2)
            vecs = cell_expanded - cand_xy                      # (N, K, 2)

            # Cross product magnitude = perpendicular distance from normal ray
            cross_mag = np.abs(vecs[:, :, 0] * cand_normals[:, :, 1] -
                               vecs[:, :, 1] * cand_normals[:, :, 0])  # (N, K)

            # Dot product = projection along normal (positive = inward)
            dot = (vecs[:, :, 0] * cand_normals[:, :, 0] +
                   vecs[:, :, 1] * cand_normals[:, :, 1])  # (N, K)

            # Score: cross_mag, with large penalty if cell is on the wrong side
            score = cross_mag.copy()
            score[dot < 0] += 1e6

            # Best candidate per cell
            best_k = np.argmin(score, axis=1)  # (N,)
            nearest_idx = idx_k[np.arange(N_c), best_k]
            distances = dists_k[np.arange(N_c), best_k]

        nearest_t = t_dense[nearest_idx]

        # Compute perpendicular depth (projection along normal)
        vecs = cortical_coords - xy_dense[nearest_idx]  # (N, 2)
        dots = np.sum(vecs * normals_dense[nearest_idx], axis=1)  # (N,)

        # Assign cells based on normal-aware pia projection
        arc_pos = nearest_t - t_min
        s_ids = np.clip((arc_pos / strip_width).astype(int), 0, n_strips_seg - 1)

        # Get bank ID for this segment (defaults to seg_i if no bank splitting)
        bank_id = seg.get("bank_id", seg_i)

        for ci in range(len(cortical_idx)):
            global_idx = cortical_idx[ci]
            # For multi-segment: keep assignment from closest pia segment
            if strip_ids[global_idx] >= 0:
                if distances[ci] >= abs(perp_depth_um[global_idx]):
                    continue
            strip_ids[global_idx] = s_ids[ci] + global_strip_offset
            strip_segment[global_idx] = seg_i
            strip_bank[global_idx] = bank_id
            perp_depth_um[global_idx] = dots[ci]
            pia_arc_pos[global_idx] = arc_pos[ci]

        global_strip_offset += n_strips_seg

    n_assigned = (strip_ids >= 0).sum()
    n_cortical = cortical_mask.sum()
    n_banks = len(set(strip_bank[strip_bank >= 0]))
    print(f"  Cell assignment: {n_assigned:,} / {n_cortical:,} cortical cells "
          f"({n_assigned/max(1,n_cortical)*100:.1f}%) assigned to strips "
          f"across {n_banks} bank(s)")

    return strip_ids, strip_segment, strip_bank, perp_depth_um, pia_arc_pos


# ── Step 5: Score Strips ─────────────────────────────────────────────

def score_strips(strip_ids, pred_depth, domains, min_cells_per_layer=MIN_CELLS_PER_LAYER):
    """Score each strip for layer completeness, ordering, and purity.

    Returns
    -------
    strip_scores : list of dict
    complete_ids : set — strips with all L1-L6 layers
    partial_ids : set — strips with 4/5 layers
    """
    # Compute depth-bin layers for all cells
    depth_layers = np.full(len(pred_depth), "Unknown", dtype=object)
    for layer_name, (lo, hi) in LAYER_BINS.items():
        mask = (pred_depth >= lo) & (pred_depth < hi)
        depth_layers[mask] = layer_name

    cortical_mask = domains == "Cortical"
    unique_strips = np.unique(strip_ids[strip_ids >= 0])

    strip_scores = []
    complete_ids = set()
    partial_ids = set()

    for s in unique_strips:
        in_strip = (strip_ids == s) & cortical_mask
        n_total = in_strip.sum()

        if n_total < min_cells_per_layer * 3:
            continue

        strip_layers = depth_layers[in_strip]
        strip_depths = pred_depth[in_strip]

        # Layer counts
        layer_counts = {l: (strip_layers == l).sum() for l in REQUIRED_LAYERS}
        layers_present = [l for l in REQUIRED_LAYERS
                         if layer_counts.get(l, 0) >= min_cells_per_layer]
        completeness = len(layers_present) / len(REQUIRED_LAYERS)

        # Depth coverage
        valid_d = strip_depths[~np.isnan(strip_depths)]
        depth_range = (np.percentile(valid_d, 95) - np.percentile(valid_d, 5)) if len(valid_d) > 10 else 0

        # Laminar order (Kendall tau)
        order_score = 0.0
        if len(layers_present) >= 3:
            expected_order = sorted(layers_present, key=lambda l: LAYER_DEPTH_ORDER[l])
            actual_medians = {l: np.nanmedian(strip_depths[strip_layers == l])
                             for l in layers_present
                             if (strip_layers == l).sum() >= min_cells_per_layer}
            if len(actual_medians) >= 3:
                actual_order = sorted(actual_medians.keys(),
                                     key=lambda l: actual_medians[l])
                expected_ranks = [expected_order.index(l) for l in actual_order]
                tau, _ = kendalltau(range(len(expected_ranks)), expected_ranks)
                order_score = max(0, tau)

        # Purity
        all_in_strip = strip_ids == s
        purity = cortical_mask[all_in_strip].sum() / max(1, all_in_strip.sum())

        score = {
            "strip_id": s,
            "n_cells": n_total,
            "completeness": completeness,
            "depth_range": depth_range,
            "order_score": order_score,
            "purity": purity,
            "layers_present": layers_present,
            "layer_counts": layer_counts,
        }
        strip_scores.append(score)

        # Tier classification
        if completeness == 1.0 and order_score > 0.8 and purity > 0.75:
            complete_ids.add(s)
        elif completeness >= 0.80 and order_score > 0.6 and purity > 0.60:
            partial_ids.add(s)

    # Print summary
    if strip_scores:
        print(f"\n  Strip Quality Summary ({len(strip_scores)} scored strips):")
        print(f"  {'S':>4} | {'N':>6} | {'Cmpl':>5} | {'Order':>5} | {'Purity':>6} "
              f"| {'DepRng':>6} | {'Layers':<25} | Tier")
        print(f"  {'-'*90}")
        for score in strip_scores:
            layers_str = ",".join(score.get("layers_present", []))
            tier = ""
            if score["strip_id"] in complete_ids:
                tier = "COMPLETE"
            elif score["strip_id"] in partial_ids:
                tier = "PARTIAL"
            print(f"  {score['strip_id']:>4} | {score['n_cells']:>6,} | "
                  f"{score['completeness']:>5.2f} | {score['order_score']:>5.2f} | "
                  f"{score['purity']*100:>5.1f}% | {score['depth_range']:>6.2f} | "
                  f"{layers_str:<25} | {tier}")

    n_complete = len(complete_ids)
    n_partial = len(partial_ids)
    n_in_complete = sum((strip_ids == s).sum() for s in complete_ids)
    n_in_partial = sum((strip_ids == s).sum() for s in partial_ids)
    cortical_n = (domains == "Cortical").sum()
    total_in_selected = n_in_complete + n_in_partial
    print(f"\n  Selected: {n_complete} complete + {n_partial} partial strips")
    print(f"  Cells in selected strips: {total_in_selected:,} / {cortical_n:,} "
          f"({total_in_selected/max(1,cortical_n)*100:.1f}% of cortical)")

    return strip_scores, complete_ids, partial_ids


# ── Step 6: Diagnostic Figure ────────────────────────────────────────

def plot_curved_strips(adata, sample_id, pia_segments, strip_ids,
                       strip_scores, complete_ids, partial_ids,
                       perp_depth_um):
    """Generate a 2×3 diagnostic figure for curved cortex strips."""
    coords = adata.obsm["spatial"]
    x, y = coords[:, 0], coords[:, 1]
    pred_depth = pd.to_numeric(adata.obs["predicted_norm_depth"], errors="coerce").values
    domains = adata.obs["banksy_domain"].values.astype(str)

    # Get L1 mask
    if "banksy_is_l1" in adata.obs.columns:
        is_l1 = adata.obs["banksy_is_l1"].values.astype(bool)
    else:
        is_l1 = domains == "Meningeal"

    fig = plt.figure(figsize=(24, 16))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.25)
    dx_label = SAMPLE_TO_DX.get(sample_id, "?")

    def setup_spatial_ax(ax, title):
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_facecolor("#0a0a0a")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title, fontsize=14, fontweight="bold", color="white",
                     pad=8, backgroundcolor="#333333")

    # Panel 1: L1 cells + pia curves
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.scatter(x[~is_l1], y[~is_l1], c="#222222", s=0.1, alpha=0.2, rasterized=True)
    ax1.scatter(x[is_l1], y[is_l1], c="#FF6B00", s=2, alpha=0.8, rasterized=True,
                label=f"L1 ({is_l1.sum():,})")

    # Overlay pia curves
    colors = ["#00FF88", "#00AAFF", "#FF00AA", "#FFFF00"]
    for seg_i, seg in enumerate(pia_segments):
        if seg["sample_xy"] is not None:
            sx = seg["sample_xy"][:, 0]
            sy = seg["sample_xy"][:, 1]
            color = colors[seg_i % len(colors)]
            ax1.plot(sx, sy, color=color, linewidth=2.5, label=f"Pia curve {seg_i}",
                     zorder=5)

    setup_spatial_ax(ax1, "L1 Cells + Pia Curves")
    ax1.legend(loc="upper right", fontsize=10, markerscale=5, framealpha=0.8,
               facecolor="#333333", labelcolor="white")

    # Panel 2: Pia curves + normals
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.scatter(x, y, c="#161616", s=0.1, alpha=0.1, rasterized=True)

    for seg_i, seg in enumerate(pia_segments):
        if seg["sample_xy"] is not None and seg["normals"] is not None:
            sx = seg["sample_xy"][:, 0]
            sy = seg["sample_xy"][:, 1]
            nx = seg["normals"][:, 0]
            ny = seg["normals"][:, 1]
            color = colors[seg_i % len(colors)]
            ax2.plot(sx, sy, color=color, linewidth=2, zorder=5)

            # Draw normals as arrows (every 10th point)
            step = max(1, len(sx) // 20)
            arrow_len = 300  # μm
            for i in range(0, len(sx), step):
                ax2.annotate("", xy=(sx[i] + nx[i]*arrow_len, sy[i] + ny[i]*arrow_len),
                            xytext=(sx[i], sy[i]),
                            arrowprops=dict(arrowstyle="->", color=color, lw=1.5),
                            zorder=6)

    setup_spatial_ax(ax2, "Pia Curves + Inward Normals")

    # Panel 3: All cortical cells colored by depth
    ax3 = fig.add_subplot(gs[0, 2])
    cortical_mask = domains == "Cortical"
    non_cortical = ~cortical_mask
    ax3.scatter(x[non_cortical], y[non_cortical], c="#222222", s=0.1, alpha=0.1,
                rasterized=True)
    if cortical_mask.sum() > 0:
        depth_colors = plt.cm.viridis(np.clip(pred_depth[cortical_mask], 0, 1))
        ax3.scatter(x[cortical_mask], y[cortical_mask], c=depth_colors, s=0.3,
                    alpha=0.5, rasterized=True)
    setup_spatial_ax(ax3, "Cortical Cells (depth)")

    # Panel 4: Cells colored by strip assignment
    ax4 = fig.add_subplot(gs[1, 0])
    selected = complete_ids | partial_ids
    in_selected = np.array([s in selected for s in strip_ids])
    not_in_strip = strip_ids < 0
    in_strip_not_selected = (strip_ids >= 0) & ~in_selected

    if not_in_strip.sum() > 0:
        ax4.scatter(x[not_in_strip], y[not_in_strip], c="#222222", s=0.1, alpha=0.2,
                    rasterized=True)
    if in_strip_not_selected.sum() > 0:
        ax4.scatter(x[in_strip_not_selected], y[in_strip_not_selected], c="#444444",
                    s=0.3, alpha=0.3, rasterized=True)
    if in_selected.sum() > 0:
        # Color by strip ID
        unique_selected = sorted(selected)
        strip_cmap = plt.cm.tab20(np.linspace(0, 1, max(1, len(unique_selected))))
        strip_color_map = {s: strip_cmap[i % len(strip_cmap)] for i, s in enumerate(unique_selected)}
        for s_id in unique_selected:
            mask = strip_ids == s_id
            tier_label = "C" if s_id in complete_ids else "P"
            ax4.scatter(x[mask], y[mask], c=[strip_color_map[s_id]], s=0.5, alpha=0.6,
                        rasterized=True, label=f"S{s_id}({tier_label})")

    # Overlay pia curves
    for seg_i, seg in enumerate(pia_segments):
        if seg["sample_xy"] is not None:
            sx = seg["sample_xy"][:, 0]
            sy = seg["sample_xy"][:, 1]
            ax4.plot(sx, sy, color="white", linewidth=1.5, zorder=5)

    n_sel = len(selected)
    setup_spatial_ax(ax4, f"Selected Strips ({n_sel} strips)")
    if n_sel <= 20:
        ax4.legend(loc="upper right", fontsize=8, markerscale=8, framealpha=0.8,
                   facecolor="#333333", labelcolor="white", ncol=2)

    # Panel 5: Cells in selected strips colored by depth-bin layer
    ax5 = fig.add_subplot(gs[1, 1])
    if not_in_strip.sum() > 0:
        ax5.scatter(x[not_in_strip | in_strip_not_selected],
                    y[not_in_strip | in_strip_not_selected],
                    c="#222222", s=0.1, alpha=0.1, rasterized=True)

    if in_selected.sum() > 0:
        layers = adata.obs["layer"].values.astype(str)
        for layer_name in LAYER_ORDER:
            if layer_name not in LAYER_COLORS:
                continue
            mask = in_selected & (layers == layer_name)
            if mask.sum() > 0:
                ax5.scatter(x[mask], y[mask], c=[LAYER_COLORS[layer_name]], s=0.5,
                            alpha=0.6, rasterized=True, label=f"{layer_name} ({mask.sum():,})")

    for seg_i, seg in enumerate(pia_segments):
        if seg["sample_xy"] is not None:
            sx = seg["sample_xy"][:, 0]
            sy = seg["sample_xy"][:, 1]
            ax5.plot(sx, sy, color="white", linewidth=1.5, zorder=5)

    setup_spatial_ax(ax5, "Selected Strips — Layers")
    ax5.legend(loc="upper right", fontsize=9, markerscale=8, framealpha=0.8,
               facecolor="#333333", labelcolor="white")

    # Panel 6: Strip quality bar chart
    ax6 = fig.add_subplot(gs[1, 2])
    if strip_scores:
        scored = [s for s in strip_scores if s["n_cells"] >= MIN_CELLS_PER_LAYER * 3]
        if scored:
            bar_x = np.arange(len(scored))
            width = 0.25
            completeness_vals = [s["completeness"] for s in scored]
            order_vals = [s["order_score"] for s in scored]
            purity_vals = [s["purity"] for s in scored]

            ax6.bar(bar_x - width, completeness_vals, width, label="Completeness",
                    color="#4CAF50", alpha=0.8)
            ax6.bar(bar_x, order_vals, width, label="Order",
                    color="#2196F3", alpha=0.8)
            ax6.bar(bar_x + width, purity_vals, width, label="Purity",
                    color="#FF9800", alpha=0.8)

            # Mark selected strips
            for i, s in enumerate(scored):
                if s["strip_id"] in complete_ids:
                    ax6.annotate("C", (i, 1.02), ha="center", fontsize=8,
                                fontweight="bold", color="#4CAF50")
                elif s["strip_id"] in partial_ids:
                    ax6.annotate("P", (i, 1.02), ha="center", fontsize=8,
                                fontweight="bold", color="#FF9800")

            ax6.set_xlabel("Strip Index", fontsize=13)
            ax6.set_ylabel("Score", fontsize=13)
            ax6.axhline(0.8, color="red", linestyle="--", linewidth=1, alpha=0.5)
            ax6.set_ylim(0, 1.1)
            ax6.legend(fontsize=11)

    ax6.set_title("Strip Quality Scores", fontsize=14, fontweight="bold")
    ax6.spines["top"].set_visible(False)
    ax6.spines["right"].set_visible(False)

    fig.suptitle(f"{sample_id} ({dx_label}) — Curved Cortex Strips",
                 fontsize=20, fontweight="bold")

    return fig


# ── Main Pipeline ────────────────────────────────────────────────────

def process_sample(sample_id, strip_width=STRIP_WIDTH_ALONG_PIA, save_h5ad=False):
    """Full curved cortex strip pipeline for one sample.

    Parameters
    ----------
    sample_id : str
    strip_width : float — strip width along pia in μm
    save_h5ad : bool — if True, save strip columns back to h5ad file
    """
    t0 = time.time()

    print(f"\n{'='*70}")
    print(f"CURVED CORTEX STRIPS: {sample_id} ({SAMPLE_TO_DX.get(sample_id, '?')})")
    print(f"{'='*70}")

    # Load sample
    path = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    adata = ad.read_h5ad(path)

    # QC filter
    if "hybrid_qc_pass" in adata.obs.columns:
        qc_mask = adata.obs["hybrid_qc_pass"].values.astype(bool)
    elif "qc_pass" in adata.obs.columns:
        qc_mask = adata.obs["qc_pass"].values.astype(bool)
    else:
        qc_mask = np.ones(adata.n_obs, dtype=bool)

    # Keep track of full-size qc_mask for h5ad saving
    qc_indices = np.where(qc_mask)[0]
    adata = adata[qc_mask].copy()
    print(f"  {adata.n_obs:,} QC-pass cells")

    # Step 1: Extract and clean L1 points
    print("\n  Step 1: Extracting L1 boundary...")
    l1_clean, l1_segment_labels, l1_clean_idx = extract_l1_points(adata)

    if len(l1_clean) == 0:
        print("  ERROR: No valid L1 points found. Cannot fit pia curve.")
        return None

    # Step 2+3: Fit pia curve and compute normals per segment
    print("\n  Step 2+3: Fitting pia curves and computing normals...")
    domains = adata.obs["banksy_domain"].values.astype(str)
    pred_depth = pd.to_numeric(adata.obs["predicted_norm_depth"], errors="coerce").values

    cortical_mask = domains == "Cortical"
    cortical_coords = adata.obsm["spatial"][cortical_mask]
    cortical_depths = pred_depth[cortical_mask]

    segments = np.unique(l1_segment_labels)
    pia_segments = []
    global_bank_counter = 0

    for seg_id in segments:
        seg_mask = l1_segment_labels == seg_id
        n_seg = seg_mask.sum()

        if n_seg < MIN_SEGMENT_L1_CELLS:
            print(f"    Skipping segment {seg_id}: only {n_seg} L1 cells")
            continue

        seg_coords = l1_clean[seg_mask]
        print(f"\n    Segment {seg_id}: {n_seg} L1 cells")

        # Fit pia curve
        spline_x, spline_y, t_range, method = fit_pia_curve(seg_coords, seg_id)

        if spline_x is None:
            print(f"    WARNING: Failed to fit pia curve for segment {seg_id}")
            continue

        # Compute normals
        n_sample = max(50, int((t_range[1] - t_range[0]) / 50))
        sample_t, sample_xy, normals, tangents = compute_normals(
            spline_x, spline_y, t_range, cortical_coords, cortical_depths,
            n_sample=n_sample
        )

        fitted_seg = {
            "segment_id": seg_id,
            "spline_x": spline_x,
            "spline_y": spline_y,
            "t_range": t_range,
            "sample_t": sample_t,
            "sample_xy": sample_xy,
            "normals": normals,
            "tangents": tangents,
            "method": method,
        }

        # Split at folds (gyral crowns / sulcal bottoms)
        banks = split_pia_at_folds(fitted_seg)
        # Assign globally unique bank IDs
        for b in banks:
            b["bank_id"] = global_bank_counter
            global_bank_counter += 1
        if len(banks) > 1:
            print(f"    Split into {len(banks)} banks at folds")
            for b in banks:
                bl = b["t_range"][1] - b["t_range"][0]
                print(f"      Bank {b['bank_id']}: arc-length {bl:.0f} μm, "
                      f"{len(b['sample_t'])} sample points")
        else:
            bl = banks[0]["t_range"][1] - banks[0]["t_range"][0]
            print(f"    Single bank (bank {banks[0]['bank_id']}): arc-length {bl:.0f} μm")
        pia_segments.extend(banks)

    if not pia_segments:
        print("  ERROR: No valid pia curves fitted. Cannot define strips.")
        return None

    # Step 4: Assign cells to strips
    print(f"\n  Step 4: Assigning cells to strips (width={strip_width}μm)...")
    strip_ids, strip_segment, strip_bank, perp_depth_um, pia_arc_pos = assign_cells_to_strips(
        adata.obsm["spatial"], domains, pia_segments, strip_width=strip_width
    )

    # Step 5: Score strips
    print("\n  Step 5: Scoring strips...")
    strip_scores, complete_ids, partial_ids = score_strips(strip_ids, pred_depth, domains)

    # Step 6: Diagnostic figure
    print("\n  Step 6: Generating diagnostic figure...")
    fig = plot_curved_strips(adata, sample_id, pia_segments, strip_ids,
                             strip_scores, complete_ids, partial_ids, perp_depth_um)
    fig_path = os.path.join(OUT_DIR, f"curved_strips_{sample_id}.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {fig_path}")

    elapsed = time.time() - t0
    n_cortical = (domains == "Cortical").sum()
    n_selected = sum((strip_ids == s).sum() for s in (complete_ids | partial_ids))
    coverage = n_selected / max(1, n_cortical) * 100

    n_banks = len(set(strip_bank[strip_bank >= 0]))
    print(f"\n  SUMMARY:")
    print(f"    Banks (distinct cortical columns): {n_banks}")
    print(f"    Pia segments (before fold-splitting): {len(set(s.get('parent_segment_id', s.get('segment_id', i)) for i, s in enumerate(pia_segments)))}")
    print(f"    Complete strips: {len(complete_ids)}")
    print(f"    Partial strips: {len(partial_ids)}")
    print(f"    Cells in selected strips: {n_selected:,} ({coverage:.1f}% of cortical)")
    print(f"    Time: {elapsed:.0f}s")

    # Save strip columns back to h5ad if requested
    if save_h5ad:
        print(f"\n  Saving strip columns to h5ad...")
        adata_full = ad.read_h5ad(path)

        # Initialize columns
        full_strip_id = np.full(adata_full.n_obs, -1, dtype=int)
        full_in_strip = np.zeros(adata_full.n_obs, dtype=bool)
        full_strip_tier = np.full(adata_full.n_obs, "", dtype=object)
        full_strip_bank = np.full(adata_full.n_obs, -1, dtype=int)

        # Map QC-pass indices
        for qi, fi in enumerate(qc_indices):
            sid = strip_ids[qi]
            full_strip_id[fi] = sid
            full_strip_bank[fi] = strip_bank[qi]
            if sid in complete_ids:
                full_in_strip[fi] = True
                full_strip_tier[fi] = "complete"
            elif sid in partial_ids:
                full_in_strip[fi] = True
                full_strip_tier[fi] = "partial"

        adata_full.obs["curved_strip_id"] = full_strip_id
        adata_full.obs["in_curved_strip"] = full_in_strip
        adata_full.obs["curved_strip_tier"] = full_strip_tier
        adata_full.obs["curved_strip_bank"] = full_strip_bank
        adata_full.write_h5ad(path)
        print(f"  Saved curved_strip_id, in_curved_strip, curved_strip_tier, "
              f"curved_strip_bank to {path}")

    return {
        "sample_id": sample_id,
        "n_banks": n_banks,
        "n_pia_segments": len(pia_segments),
        "n_complete_strips": len(complete_ids),
        "n_partial_strips": len(partial_ids),
        "n_cells_selected": n_selected,
        "coverage_pct": coverage,
        "time_sec": elapsed,
    }


# ── CLI ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Curved cortex strip identification using L1/pia boundary")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sample", type=str,
                       help="Sample ID (e.g., Br2421)")
    group.add_argument("--all", action="store_true",
                       help="Process all samples sequentially")
    parser.add_argument("--strip-width", type=float, default=STRIP_WIDTH_ALONG_PIA,
                        help=f"Strip width along pia in μm (default: {STRIP_WIDTH_ALONG_PIA})")
    parser.add_argument("--save", action="store_true",
                        help="Save strip columns (curved_strip_id, in_curved_strip, "
                             "curved_strip_tier) back to h5ad files")
    args = parser.parse_args()

    if args.all:
        # Process all samples
        all_samples = sorted(SAMPLE_TO_DX.keys())
        print(f"Processing {len(all_samples)} samples...")
        t_total = time.time()

        results = []
        for sample_id in all_samples:
            try:
                result = process_sample(sample_id, strip_width=args.strip_width,
                                        save_h5ad=args.save)
                if result:
                    results.append(result)
            except Exception as e:
                print(f"\n  ERROR processing {sample_id}: {e}")
                import traceback
                traceback.print_exc()
                results.append({
                    "sample_id": sample_id,
                    "n_banks": 0,
                    "n_pia_segments": 0,
                    "n_complete_strips": 0,
                    "n_partial_strips": 0,
                    "n_cells_selected": 0,
                    "coverage_pct": 0,
                    "time_sec": 0,
                })

        # Save batch summary
        if results:
            df = pd.DataFrame(results)
            df["diagnosis"] = df["sample_id"].map(SAMPLE_TO_DX)
            cols = ["sample_id", "diagnosis", "n_banks", "n_pia_segments",
                    "n_complete_strips", "n_partial_strips",
                    "n_cells_selected", "coverage_pct", "time_sec"]
            df = df[cols]
            summary_path = os.path.join(OUT_DIR, "curved_strips_batch_summary.csv")
            df.to_csv(summary_path, index=False)
            print(f"\n{'='*70}")
            print(f"BATCH SUMMARY ({len(results)} samples, {time.time()-t_total:.0f}s total)")
            print(f"{'='*70}")
            print(df.to_string(index=False))
            print(f"\nMean coverage: {df.coverage_pct.mean():.1f}%")
            print(f"Median coverage: {df.coverage_pct.median():.1f}%")
            print(f"Min coverage: {df.coverage_pct.min():.1f}% ({df.loc[df.coverage_pct.idxmin(), 'sample_id']})")
            print(f"Max coverage: {df.coverage_pct.max():.1f}% ({df.loc[df.coverage_pct.idxmax(), 'sample_id']})")
            print(f"\nSaved: {summary_path}")
    else:
        result = process_sample(args.sample, strip_width=args.strip_width,
                                save_h5ad=args.save)
        if result:
            print(f"\nDone! Check {OUT_DIR}/curved_strips_{args.sample}.png")


if __name__ == "__main__":
    main()
