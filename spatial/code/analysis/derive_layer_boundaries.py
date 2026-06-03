#!/usr/bin/env python3
"""
Derive cortical layer boundaries from MERFISH manual annotations and
validate against Xenium predicted depths.

Method:
  For each pair of adjacent cortical layers (L2/3↔L4, L4↔L5, L5↔L6),
  compute the depth at which the pairwise proportion of excitatory marker
  neurons A/(A+B) crosses 0.5. This is the depth where the upper layer's
  neurons become less common than the lower layer's neurons.

  L1/L2-3 boundary: depth where excitatory neuron density first exceeds
  25% of its peak (going superficial → deep), marking the onset of L2/3.

  L6/WM boundary: depth where L6 excitatory neuron density drops below
  10% of its peak (going deep → WM).

  Boundaries are computed independently in both MERFISH (manual depth
  annotations = ground truth) and Xenium (predicted depths). The Xenium
  boundaries are used for the pipeline since they reflect the depth model's
  actual output space.

Layer marker types:
  L2/3: L2/3 IT
  L4:   L4 IT
  L5:   L5 IT, L5 ET, L5/6 NP
  L6:   L6 CT, L6 IT, L6 IT Car3, L6b

Output:
  output/depth_proportions/proposed_layer_boundaries.csv
  output/depth_proportions/layer_boundary_comparison.png
  Prints proposed LAYER_BINS dict for copy-paste into depth_model.py

Usage:
    python3 -u derive_layer_boundaries.py
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR, load_merfish_cortical

DEPTH_DIR = os.path.join(BASE_DIR, "output", "depth_proportions")
CELL_LEVEL_PATH = os.path.join(DEPTH_DIR, "cell_level_data.csv")

# Excitatory marker types used to define layer boundaries
LAYER_MARKERS = {
    "L2/3": ["L2/3 IT"],
    "L4": ["L4 IT"],
    "L5": ["L5 IT", "L5 ET", "L5/6 NP"],
    "L6": ["L6 CT", "L6 IT", "L6 IT Car3", "L6b"],
}
ALL_EXC_TYPES = [t for types in LAYER_MARKERS.values() for t in types]

# Current boundaries for comparison
CURRENT_BOUNDS = {
    "L1/L2-3": 0.1225,
    "L2-3/L4": 0.4696,
    "L4/L5": 0.5443,
    "L5/L6": 0.7079,
    "L6/WM": 0.9275,
}

# Fine histogram bins for KDE-like smoothing
N_HIST_BINS = 200
SMOOTH_SIGMA = 3  # Gaussian filter sigma (in bin units)


def pairwise_crossover(depths_a, depths_b, bin_edges, bin_mids):
    """Find the depth where A/(A+B) crosses 0.5 from above.

    Parameters
    ----------
    depths_a, depths_b : array-like
        Depth values for the upper and lower layer marker types.

    Returns
    -------
    float or np.nan
        Crossover depth, or NaN if no crossing found.
    """
    ha, _ = np.histogram(depths_a, bins=bin_edges)
    hb, _ = np.histogram(depths_b, bins=bin_edges)
    ha_s = gaussian_filter1d(ha.astype(float), sigma=SMOOTH_SIGMA)
    hb_s = gaussian_filter1d(hb.astype(float), sigma=SMOOTH_SIGMA)
    total = ha_s + hb_s

    # Pairwise fraction: A / (A + B)
    frac_a = np.where(total > 5, ha_s / total, np.nan)

    # Find first downward crossing of 0.5
    for i in range(len(frac_a) - 1):
        if not np.isnan(frac_a[i]) and not np.isnan(frac_a[i + 1]):
            if frac_a[i] > 0.5 and frac_a[i + 1] <= 0.5:
                # Linear interpolation for sub-bin precision
                f = (frac_a[i] - 0.5) / (frac_a[i] - frac_a[i + 1])
                return bin_mids[i] + f * (bin_mids[i + 1] - bin_mids[i])
    return np.nan


def compute_all_boundaries(data, depth_col, type_col):
    """Compute all layer boundaries from one dataset.

    Parameters
    ----------
    data : DataFrame
        Must contain depth_col and type_col columns.
    depth_col : str
        Column with normalized depth values (0=pia, 1=WM).
    type_col : str
        Column with subclass labels.

    Returns
    -------
    dict
        Boundary name -> depth value.
    """
    bin_edges = np.linspace(0, 1, N_HIST_BINS + 1)
    bin_mids = (bin_edges[:-1] + bin_edges[1:]) / 2

    exc = data[data[type_col].isin(ALL_EXC_TYPES)].copy()
    exc_depths = exc[depth_col].dropna().clip(0, 1).values

    boundaries = {}

    # ── L1/L2-3: where excitatory density first exceeds 25% of peak ──
    h_all, _ = np.histogram(exc_depths, bins=bin_edges)
    h_all_s = gaussian_filter1d(h_all.astype(float), sigma=SMOOTH_SIGMA)
    peak = h_all_s.max()
    l1_boundary = np.nan
    for i in range(len(h_all_s)):
        if h_all_s[i] > peak * 0.25:
            l1_boundary = bin_mids[i]
            break
    boundaries["L1/L2-3"] = l1_boundary

    # ── Pairwise crossovers between adjacent layers ──
    pairs = [
        ("L2-3/L4", "L2/3", "L4"),
        ("L4/L5", "L4", "L5"),
        ("L5/L6", "L5", "L6"),
    ]
    for name, upper, lower in pairs:
        a_types = LAYER_MARKERS[upper]
        b_types = LAYER_MARKERS[lower]
        a_depths = data[data[type_col].isin(a_types)][depth_col].dropna().clip(0, 1).values
        b_depths = data[data[type_col].isin(b_types)][depth_col].dropna().clip(0, 1).values
        boundaries[name] = pairwise_crossover(a_depths, b_depths,
                                               bin_edges, bin_mids)

    # ── L6/WM: where L6 marker density drops below 10% of peak ──
    l6_depths = data[data[type_col].isin(LAYER_MARKERS["L6"])][depth_col].dropna().clip(0, 1).values
    h6, _ = np.histogram(l6_depths, bins=bin_edges)
    h6_s = gaussian_filter1d(h6.astype(float), sigma=SMOOTH_SIGMA)
    peak6 = h6_s.max()
    peak6_idx = np.argmax(h6_s)
    l6_wm = np.nan
    for i in range(peak6_idx, len(h6_s)):
        if h6_s[i] < peak6 * 0.10:
            l6_wm = bin_mids[i]
            break
    boundaries["L6/WM"] = l6_wm

    return boundaries


def plot_boundary_derivation(merfish_raw, cells_raw, mer_bounds, xen_bounds):
    """Three-panel figure: MERFISH crossovers, Xenium crossovers, comparison."""
    bin_edges = np.linspace(0, 1, N_HIST_BINS + 1)
    bin_mids = (bin_edges[:-1] + bin_edges[1:]) / 2

    layer_colors = {"L2/3": "#56B4E9", "L4": "#009E73",
                    "L5": "#CC79A7", "L6": "#D55E00"}

    fig, axes = plt.subplots(1, 3, figsize=(22, 8))

    # Panels 1-2: pairwise fractions for MERFISH and Xenium
    datasets = [
        (axes[0], merfish_raw, "depth", "subclass",
         f"MERFISH (n={merfish_raw['donor'].nunique()})", mer_bounds),
        (axes[1], cells_raw, "predicted_norm_depth", "subclass_label",
         "Xenium (Controls only)", xen_bounds),
    ]

    pairs = [
        ("L2/3", LAYER_MARKERS["L2/3"], "L4", LAYER_MARKERS["L4"]),
        ("L4", LAYER_MARKERS["L4"], "L5", LAYER_MARKERS["L5"]),
        ("L5", LAYER_MARKERS["L5"], "L6", LAYER_MARKERS["L6"]),
    ]

    for ax, data, dcol, tcol, title, bounds in datasets:
        for l1_name, l1_types, l2_name, l2_types in pairs:
            a = data[data[tcol].isin(l1_types)][dcol].dropna().clip(0, 1).values
            b = data[data[tcol].isin(l2_types)][dcol].dropna().clip(0, 1).values
            ha, _ = np.histogram(a, bins=bin_edges)
            hb, _ = np.histogram(b, bins=bin_edges)
            ha_s = gaussian_filter1d(ha.astype(float), sigma=SMOOTH_SIGMA)
            hb_s = gaussian_filter1d(hb.astype(float), sigma=SMOOTH_SIGMA)
            total = ha_s + hb_s
            frac_a = np.where(total > 5, ha_s / total, np.nan)

            ax.plot(frac_a, bin_mids, color=layer_colors[l1_name], linewidth=2.0,
                    label=f"{l1_name} / ({l1_name}+{l2_name})")

        ax.axvline(0.5, color="#999999", linewidth=1.0, linestyle="--", alpha=0.5)

        # Mark computed boundaries
        bnames = ["L2-3/L4", "L4/L5", "L5/L6"]
        bcolors = ["#56B4E9", "#009E73", "#CC79A7"]
        for bname, bcolor in zip(bnames, bcolors):
            val = bounds.get(bname, np.nan)
            if not np.isnan(val):
                ax.axhline(val, color=bcolor, linewidth=2.0, alpha=0.7)
                ax.text(0.03, val + 0.008, f"{val:.3f}", fontsize=11,
                       color=bcolor, fontweight="bold", va="bottom")

        # Current boundaries as dotted
        for cval in CURRENT_BOUNDS.values():
            ax.axhline(cval, color="#333333", linewidth=0.8,
                      linestyle=":", alpha=0.3)

        ax.set_ylim(1.0, 0.0)
        ax.set_xlim(0, 1)
        ax.set_ylabel("Cortical Depth (pia → WM)", fontsize=14)
        ax.set_xlabel("Pairwise fraction", fontsize=14)
        ax.set_title(title, fontsize=15, fontweight="bold")
        ax.legend(fontsize=11, loc="lower left")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=12)

    # Panel 3: summary comparison table
    ax = axes[2]
    ax.set_xlim(0, 10)
    ax.set_ylim(1.0, 0.0)

    layer_fill = {"L1": "#FFF3E0", "L2/3": "#E3F2FD", "L4": "#E8F5E9",
                  "L5": "#FFF9C4", "L6": "#FCE4EC"}
    layer_names = ["L1", "L2/3", "L4", "L5", "L6"]
    bound_names = ["L1/L2-3", "L2-3/L4", "L4/L5", "L5/L6", "L6/WM"]

    # Current (left column)
    curr_vals = [0.0] + [CURRENT_BOUNDS[b] for b in bound_names]
    for i, layer in enumerate(layer_names):
        lo, hi = curr_vals[i], curr_vals[i + 1]
        ax.axhspan(lo, hi, xmin=0.02, xmax=0.30, alpha=0.4,
                   color=layer_fill[layer])
        ax.text(1.6, (lo + hi) / 2, layer, fontsize=13, ha="center",
               va="center", fontweight="bold")
        ax.axhline(lo, xmin=0.02, xmax=0.30, color="#333333",
                  linewidth=1.0, zorder=5)
    ax.axhline(curr_vals[-1], xmin=0.02, xmax=0.30, color="#333333",
              linewidth=1.0, zorder=5)
    ax.text(1.6, 0.96, "Current", fontsize=12, ha="center",
           fontweight="bold", color="#555555")

    # Proposed (right column)
    prop_vals = [0.0] + [xen_bounds.get(b, np.nan) for b in bound_names]
    for i, layer in enumerate(layer_names):
        lo, hi = prop_vals[i], prop_vals[i + 1]
        if np.isnan(lo) or np.isnan(hi):
            continue
        ax.axhspan(lo, hi, xmin=0.70, xmax=0.98, alpha=0.4,
                   color=layer_fill[layer])
        ax.text(8.4, (lo + hi) / 2, layer, fontsize=13, ha="center",
               va="center", fontweight="bold")
        ax.axhline(lo, xmin=0.70, xmax=0.98, color="#333333",
                  linewidth=1.0, zorder=5)
    if not np.isnan(prop_vals[-1]):
        ax.axhline(prop_vals[-1], xmin=0.70, xmax=0.98, color="#333333",
                  linewidth=1.0, zorder=5)
    ax.text(8.4, 0.96, "Proposed", fontsize=12, ha="center",
           fontweight="bold", color="#555555")

    # Arrows
    for i in range(1, len(curr_vals)):
        c = curr_vals[i]
        p = prop_vals[i]
        if np.isnan(p) or abs(p - c) < 0.003:
            continue
        mid_y = (c + p) / 2
        ax.annotate("", xy=(7.0, p), xytext=(3.0, c),
                    arrowprops=dict(arrowstyle="->", color="#CC3311",
                                   linewidth=1.5,
                                   connectionstyle="arc3,rad=0.05"))
        ax.text(5.0, mid_y, f"{p - c:+.3f}", fontsize=11, ha="center",
               va="center", color="#CC3311", fontweight="bold",
               bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                        edgecolor="#CC3311", alpha=0.8))

    ax.set_xticks([])
    ax.set_title("Boundary Comparison", fontsize=15, fontweight="bold")
    ax.set_ylabel("Cortical Depth (pia → WM)", fontsize=14)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(labelsize=12)

    fig.suptitle("Layer Boundary Derivation from Excitatory Marker Crossovers\n"
                 "(dotted = current boundaries)",
                 fontsize=18, fontweight="bold", y=1.02)
    fig.tight_layout()

    outpath = os.path.join(DEPTH_DIR, "layer_boundary_comparison.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {outpath}")
    plt.close(fig)


def main():
    os.makedirs(DEPTH_DIR, exist_ok=True)

    # ── Load MERFISH reference (manual depth annotations = ground truth) ──
    print("Loading MERFISH reference...")
    merfish = load_merfish_cortical()
    merfish["depth"] = merfish["depth"].clip(0, 1)
    print(f"  {len(merfish):,} cells, {merfish['donor'].nunique()} donors")
    print(f"  Excitatory markers: {sum(merfish['subclass'].isin(ALL_EXC_TYPES)):,} cells")

    # ── Load Xenium Controls ──
    print("\nLoading Xenium Controls...")
    cells = pd.read_csv(CELL_LEVEL_PATH)
    cells_ctrl = cells[cells["diagnosis"] == "Control"].copy()
    print(f"  {len(cells_ctrl):,} cells, {cells_ctrl['sample_id'].nunique()} donors")

    # ── Compute boundaries ──
    print("\nComputing MERFISH boundaries (ground truth)...")
    mer_bounds = compute_all_boundaries(merfish, "depth", "subclass")
    print("Computing Xenium boundaries (predicted depths)...")
    xen_bounds = compute_all_boundaries(cells_ctrl, "predicted_norm_depth",
                                         "subclass_label")

    # ── Print comparison table ──
    bound_names = ["L1/L2-3", "L2-3/L4", "L4/L5", "L5/L6", "L6/WM"]
    print(f"\n{'Boundary':12s}  {'Current':>8s}  {'MERFISH':>8s}  {'Xenium':>8s}  "
          f"{'Xen-Curr':>10s}")
    print("-" * 55)
    for name in bound_names:
        c = CURRENT_BOUNDS[name]
        m = mer_bounds[name]
        x = xen_bounds[name]
        shift = x - c
        print(f"{name:12s}  {c:8.3f}  {m:8.4f}  {x:8.4f}  {shift:+10.4f}")

    # ── Print proposed LAYER_BINS for depth_model.py ──
    print("\n" + "=" * 60)
    print("PROPOSED LAYER_BINS (copy to code/modules/depth_model.py):")
    print("=" * 60)
    b = xen_bounds
    # Round to 4 decimal places for precision
    print(f"""
LAYER_BINS = {{
    'L1':   (-np.inf, {b['L1/L2-3']:.4f}),
    'L2/3': ({b['L1/L2-3']:.4f}, {b['L2-3/L4']:.4f}),
    'L4':   ({b['L2-3/L4']:.4f}, {b['L4/L5']:.4f}),
    'L5':   ({b['L4/L5']:.4f}, {b['L5/L6']:.4f}),
    'L6':   ({b['L5/L6']:.4f}, {b['L6/WM']:.4f}),
    'WM':   ({b['L6/WM']:.4f}, np.inf),
}}

DEPTH_STRATA = {{
    'L2/3': ({b['L1/L2-3']:.4f}, {b['L2-3/L4']:.4f}),
    'L4':   ({b['L2-3/L4']:.4f}, {b['L4/L5']:.4f}),
    'L5':   ({b['L4/L5']:.4f}, {b['L5/L6']:.4f}),
    'L6':   ({b['L5/L6']:.4f}, {b['L6/WM']:.4f}),
}}""")

    # ── Print proposed layer ranges ──
    print(f"\nProposed layer ranges:")
    layers = ["L1", "L2/3", "L4", "L5", "L6", "WM"]
    bound_vals = [0.0, b["L1/L2-3"], b["L2-3/L4"], b["L4/L5"],
                  b["L5/L6"], b["L6/WM"], 1.0]
    for i, layer in enumerate(layers):
        lo, hi = bound_vals[i], bound_vals[i + 1]
        width = hi - lo
        print(f"  {layer:5s}: {lo:.4f} - {hi:.4f}  (width={width:.4f})")

    # ── Save CSV ──
    rows = []
    for name in bound_names:
        rows.append({
            "boundary": name,
            "current": CURRENT_BOUNDS[name],
            "merfish_crossover": mer_bounds[name],
            "xenium_crossover": xen_bounds[name],
            "shift_from_current": xen_bounds[name] - CURRENT_BOUNDS[name],
        })
    df = pd.DataFrame(rows)
    outpath = os.path.join(DEPTH_DIR, "proposed_layer_boundaries.csv")
    df.to_csv(outpath, index=False)
    print(f"\nSaved: {outpath}")

    # ── Plot ──
    print("\nPlotting boundary derivation figure...")
    plot_boundary_derivation(merfish, cells_ctrl, mer_bounds, xen_bounds)

    print("\nDone!")


if __name__ == "__main__":
    main()
