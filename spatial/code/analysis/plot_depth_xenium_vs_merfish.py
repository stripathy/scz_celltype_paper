#!/usr/bin/env python3
"""
Compare Xenium vs SEA-AD MERFISH depth profiles and derive optimal
layer boundaries.

Three figures:
  1. Neuronal depth profiles: Xenium Control vs MERFISH overlay
  2. Non-neuronal depth profiles: same
  3. Layer boundary derivation: shows pairwise excitatory crossovers
     in both datasets, current vs proposed boundaries

Boundary method: for each pair of adjacent layers, compute the depth
where the smoothed pairwise proportion A/(A+B) = 0.5, using excitatory
marker neuron types.

Output:
  output/depth_proportions/xenium_vs_merfish_neuronal.png
  output/depth_proportions/xenium_vs_merfish_nonneuronal.png
  output/depth_proportions/layer_boundary_comparison.png
  output/depth_proportions/proposed_layer_boundaries.csv
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR, infer_class, load_merfish_cortical

CRUMBLR_DIR = os.path.join(BASE_DIR, "output", "crumblr")
DEPTH_DIR = os.path.join(BASE_DIR, "output", "depth_proportions")

XENIUM_COLOR = "#CC3311"
MERFISH_COLOR = "#33BBEE"

LAYER_BOUNDS_CURRENT = {
    "L1": (0.00, 0.1225),
    "L2/3": (0.1225, 0.4696),
    "L4": (0.4696, 0.5443),
    "L5": (0.5443, 0.7079),
    "L6": (0.7079, 0.9275),
}

LAYER_TYPES = {
    "L2/3 IT": "L2/3", "L4 IT": "L4",
    "L5 IT": "L5", "L5 ET": "L5", "L5/6 NP": "L5",
    "L6 CT": "L6", "L6 IT": "L6", "L6 IT Car3": "L6", "L6b": "L6",
}


def compute_mean_profile(df, celltype, class_filter=None):
    """Compute mean proportion (%) per depth bin across donors."""
    data = df.copy()
    if class_filter is not None:
        class_cts = [ct for ct in data["celltype"].unique()
                     if _matches_class(ct, class_filter)]
        data = data[data["celltype"].isin(class_cts)]
        class_totals = (data.groupby(["donor", "depth_bin"])["count"]
                        .sum().reset_index(name="class_total"))
        data = data.merge(class_totals, on=["donor", "depth_bin"])
        data["pct"] = 100.0 * data["count"] / data["class_total"]
    else:
        data["pct"] = 100.0 * data["count"] / data["total"]

    ct_data = data[data["celltype"] == celltype]
    if len(ct_data) == 0:
        return None
    summary = (ct_data.groupby("depth_midpoint")["pct"]
               .agg(["mean", "sem"]).reset_index().sort_values("depth_midpoint"))
    return summary


def _matches_class(celltype, class_filter):
    cls = infer_class(celltype)
    if class_filter == "neuronal":
        return cls in ("Glut", "GABA")
    elif class_filter == "nonneuronal":
        return cls in ("NN", "Other")
    return True


def compute_boundaries(data, depth_col, type_col):
    """Compute layer boundaries from excitatory marker crossovers.

    Returns dict with boundary names and depths.
    """
    bin_edges = np.linspace(0, 1, 201)
    bin_mids = (bin_edges[:-1] + bin_edges[1:]) / 2

    def pairwise_cross(a_types, b_types):
        a = data[data[type_col].isin(a_types)][depth_col].dropna().clip(0, 1).values
        b = data[data[type_col].isin(b_types)][depth_col].dropna().clip(0, 1).values
        ha, _ = np.histogram(a, bins=bin_edges)
        hb, _ = np.histogram(b, bins=bin_edges)
        ha_s = gaussian_filter1d(ha.astype(float), sigma=3)
        hb_s = gaussian_filter1d(hb.astype(float), sigma=3)
        total = ha_s + hb_s
        frac_a = np.where(total > 5, ha_s / total, 0.5)
        for i in range(len(frac_a) - 1):
            if frac_a[i] > 0.5 and frac_a[i + 1] <= 0.5:
                f = (frac_a[i] - 0.5) / (frac_a[i] - frac_a[i + 1])
                return bin_mids[i] + f * (bin_mids[i + 1] - bin_mids[i])
        return np.nan

    # L1/L2-3: where excitatory density first exceeds 25% of peak
    exc_types = list(LAYER_TYPES.keys())
    exc_depths = data[data[type_col].isin(exc_types)][depth_col].dropna().clip(0, 1).values
    h, _ = np.histogram(exc_depths, bins=bin_edges)
    h_s = gaussian_filter1d(h.astype(float), sigma=3)
    peak = h_s.max()
    l1_l23 = np.nan
    for i in range(len(h_s)):
        if h_s[i] > peak * 0.25:
            l1_l23 = bin_mids[i]
            break

    # Pairwise crossovers
    l23_l4 = pairwise_cross(["L2/3 IT"], ["L4 IT"])
    l4_l5 = pairwise_cross(["L4 IT"], ["L5 IT", "L5 ET", "L5/6 NP"])
    l5_l6 = pairwise_cross(["L5 IT", "L5 ET", "L5/6 NP"],
                            ["L6 CT", "L6 IT", "L6 IT Car3", "L6b"])

    # L6/WM: where L6 marker density drops to 10% of peak
    l6_types = ["L6 CT", "L6 IT", "L6 IT Car3", "L6b"]
    l6_depths = data[data[type_col].isin(l6_types)][depth_col].dropna().clip(0, 1).values
    h6, _ = np.histogram(l6_depths, bins=bin_edges)
    h6_s = gaussian_filter1d(h6.astype(float), sigma=3)
    peak6 = h6_s.max()
    peak6_idx = np.argmax(h6_s)
    l6_wm = np.nan
    for i in range(peak6_idx, len(h6_s)):
        if h6_s[i] < peak6 * 0.10:
            l6_wm = bin_mids[i]
            break

    return {
        "L1/L2-3": l1_l23,
        "L2-3/L4": l23_l4,
        "L4/L5": l4_l5,
        "L5/L6": l5_l6,
        "L6/WM": l6_wm,
    }


def plot_comparison_profiles(xenium_df, merfish_df, celltypes, class_label,
                              class_filter, filename, proposed_bounds,
                              ncols=6):
    """Depth profiles overlaying Xenium Control and MERFISH."""
    n_ct = len(celltypes)
    nrows = int(np.ceil(n_ct / ncols))
    panel_w = 3.0
    panel_h = 4.2 if nrows == 1 else 3.8
    fig_w = panel_w * ncols + 2.2
    fig_h = panel_h * nrows + 2.5

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(nrows, ncols, figure=fig,
                           wspace=0.06, hspace=0.28,
                           left=0.13, right=0.97,
                           top=0.90, bottom=0.08)

    for idx, ct in enumerate(celltypes):
        row, col = divmod(idx, ncols)
        ax = fig.add_subplot(gs[row, col])
        is_left = (col == 0)
        is_bottom = (row == nrows - 1) or (idx + ncols >= n_ct)

        # Layer shading using PROPOSED boundaries
        shaded = False
        layer_names = ["L1", "L2/3", "L4", "L5", "L6"]
        bounds_list = list(proposed_bounds.values())
        layer_ranges = list(zip([0.0] + bounds_list[:-1], bounds_list))

        for i, (lo, hi) in enumerate(layer_ranges):
            if shaded:
                ax.axhspan(lo, hi, alpha=0.045, color="#000000", zorder=0)
            ax.axhline(lo, color="#d0d0d0", linewidth=0.5, zorder=1)
            shaded = not shaded
        ax.axhline(layer_ranges[-1][1], color="#d0d0d0", linewidth=0.5, zorder=1)

        if is_left:
            for i, (lo, hi) in enumerate(layer_ranges):
                if i < len(layer_names):
                    ax.text(-0.22, (lo + hi) / 2, layer_names[i],
                            transform=ax.get_yaxis_transform(),
                            ha="center", va="center", fontsize=11,
                            color="#777777", fontweight="bold", clip_on=False)

        # Xenium
        xen_profile = compute_mean_profile(xenium_df, ct, class_filter)
        if xen_profile is not None:
            ax.plot(xen_profile["mean"], xen_profile["depth_midpoint"],
                    color=XENIUM_COLOR, linewidth=2.0, zorder=5)
            ax.fill_betweenx(xen_profile["depth_midpoint"],
                             xen_profile["mean"] - xen_profile["sem"],
                             xen_profile["mean"] + xen_profile["sem"],
                             color=XENIUM_COLOR, alpha=0.15, zorder=4)

        # MERFISH
        mer_profile = compute_mean_profile(merfish_df, ct, class_filter)
        if mer_profile is not None:
            ax.plot(mer_profile["mean"], mer_profile["depth_midpoint"],
                    color=MERFISH_COLOR, linewidth=2.0, zorder=5)
            ax.fill_betweenx(mer_profile["depth_midpoint"],
                             mer_profile["mean"] - mer_profile["sem"],
                             mer_profile["mean"] + mer_profile["sem"],
                             color=MERFISH_COLOR, alpha=0.15, zorder=4)

        ax.set_ylim(1.0, 0.0)
        ax.set_xlim(left=0)
        ax.tick_params(labelsize=12, length=3, width=0.5)
        if is_left:
            ax.set_ylabel("Cortical Depth", fontsize=16, labelpad=45)
        else:
            ax.tick_params(axis="y", labelleft=False)
        if is_bottom:
            ax.set_xlabel("% of class", fontsize=14)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(0.6)
        ax.spines["bottom"].set_linewidth(0.6)
        ax.spines["left"].set_color("#888888")
        ax.spines["bottom"].set_color("#888888")
        ax.set_title(ct, fontsize=15, fontweight="bold", pad=6)

    legend_elements = [
        Line2D([0], [0], color=XENIUM_COLOR, lw=2.5,
               label=f"Xenium Control (n={xenium_df['donor'].nunique()})"),
        Line2D([0], [0], color=MERFISH_COLOR, lw=2.5,
               label=f"MERFISH (n={merfish_df['donor'].nunique()})"),
    ]
    fig.legend(handles=legend_elements, loc="upper right",
              fontsize=14, framealpha=0.95, edgecolor="#cccccc",
              bbox_to_anchor=(0.97, 0.99))
    fig.suptitle(f"{class_label} — Xenium vs MERFISH Depth Profiles",
                 fontsize=20, fontweight="bold", x=0.13, ha="left", y=1.0)

    outpath = os.path.join(DEPTH_DIR, filename)
    fig.savefig(outpath, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {outpath}")
    plt.close(fig)


def plot_boundary_derivation(merfish_raw, cells_raw, mer_bounds, xen_bounds):
    """Show how boundaries are derived from excitatory marker crossovers."""
    bin_edges = np.linspace(0, 1, 201)
    bin_mids = (bin_edges[:-1] + bin_edges[1:]) / 2

    layer_colors = {"L2/3": "#56B4E9", "L4": "#009E73",
                    "L5": "#CC79A7", "L6": "#D55E00"}

    fig, axes = plt.subplots(1, 3, figsize=(22, 8))

    # Panel 1: MERFISH pairwise proportions
    # Panel 2: Xenium pairwise proportions
    for panel_idx, (ax, data, depth_col, type_col, title, bounds) in enumerate([
        (axes[0], merfish_raw, "depth", "subclass",
         f"MERFISH (n={merfish_raw['donor'].nunique()})", mer_bounds),
        (axes[1], cells_raw, "predicted_norm_depth", "subclass_label",
         f"Xenium Controls", xen_bounds),
    ]):
        pairs = [
            ("L2/3", ["L2/3 IT"], "L4", ["L4 IT"]),
            ("L4", ["L4 IT"], "L5", ["L5 IT", "L5 ET", "L5/6 NP"]),
            ("L5", ["L5 IT", "L5 ET", "L5/6 NP"],
             "L6", ["L6 CT", "L6 IT", "L6 IT Car3", "L6b"]),
        ]

        for l1_name, l1_types, l2_name, l2_types in pairs:
            a = data[data[type_col].isin(l1_types)][depth_col].dropna().clip(0, 1).values
            b = data[data[type_col].isin(l2_types)][depth_col].dropna().clip(0, 1).values
            ha, _ = np.histogram(a, bins=bin_edges)
            hb, _ = np.histogram(b, bins=bin_edges)
            ha_s = gaussian_filter1d(ha.astype(float), sigma=3)
            hb_s = gaussian_filter1d(hb.astype(float), sigma=3)
            total = ha_s + hb_s
            frac_a = np.where(total > 5, ha_s / total, np.nan)

            color = layer_colors[l1_name]
            ax.plot(frac_a, bin_mids, color=color, linewidth=2.0,
                    label=f"{l1_name} / ({l1_name}+{l2_name})")

        ax.axvline(0.5, color="#999999", linewidth=1.0, linestyle="--", alpha=0.5)

        # Mark boundaries
        boundary_names = ["L2-3/L4", "L4/L5", "L5/L6"]
        boundary_colors = ["#56B4E9", "#009E73", "#CC79A7"]
        for bname, bcolor in zip(boundary_names, boundary_colors):
            val = bounds.get(bname, np.nan)
            if not np.isnan(val):
                ax.axhline(val, color=bcolor, linewidth=2.0,
                          linestyle="-", alpha=0.7)
                ax.text(0.03, val + 0.01, f"{val:.3f}",
                       fontsize=11, color=bcolor, fontweight="bold",
                       va="bottom")

        # Also mark current boundaries as dashed
        for _, (lo, hi) in LAYER_BOUNDS_CURRENT.items():
            ax.axhline(lo, color="#333333", linewidth=0.8,
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

    # Panel 3: boundary comparison table
    ax = axes[2]
    ax.set_xlim(0, 10)
    ax.set_ylim(1.0, 0.0)

    layer_fill = {"L1": "#FFF3E0", "L2/3": "#E3F2FD", "L4": "#E8F5E9",
                  "L5": "#FFF9C4", "L6": "#FCE4EC"}
    layer_names = ["L1", "L2/3", "L4", "L5", "L6"]

    # Current boundaries (left)
    curr_bounds = [0.0, 0.10, 0.40, 0.55, 0.70, 0.90]
    for i, layer in enumerate(layer_names):
        lo, hi = curr_bounds[i], curr_bounds[i + 1]
        ax.axhspan(lo, hi, xmin=0.02, xmax=0.30,
                   alpha=0.4, color=layer_fill[layer])
        ax.text(1.6, (lo + hi) / 2, layer, fontsize=13,
               ha="center", va="center", fontweight="bold")
        ax.axhline(lo, xmin=0.02, xmax=0.30, color="#333333",
                  linewidth=1.0, zorder=5)
    ax.axhline(0.90, xmin=0.02, xmax=0.30, color="#333333",
              linewidth=1.0, zorder=5)
    ax.text(1.6, 0.96, "Current", fontsize=12, ha="center",
           fontweight="bold", color="#555555")

    # Proposed boundaries (right) — using Xenium values
    xen_b = xen_bounds
    prop_bounds = [0.0,
                   xen_b.get("L1/L2-3", 0.12),
                   xen_b.get("L2-3/L4", 0.47),
                   xen_b.get("L4/L5", 0.54),
                   xen_b.get("L5/L6", 0.70),
                   xen_b.get("L6/WM", 0.93)]
    for i, layer in enumerate(layer_names):
        lo, hi = prop_bounds[i], prop_bounds[i + 1]
        ax.axhspan(lo, hi, xmin=0.70, xmax=0.98,
                   alpha=0.4, color=layer_fill[layer])
        ax.text(8.4, (lo + hi) / 2, layer, fontsize=13,
               ha="center", va="center", fontweight="bold")
        ax.axhline(lo, xmin=0.70, xmax=0.98, color="#333333",
                  linewidth=1.0, zorder=5)
    ax.axhline(prop_bounds[-1], xmin=0.70, xmax=0.98, color="#333333",
              linewidth=1.0, zorder=5)
    ax.text(8.4, 0.96, "Proposed", fontsize=12, ha="center",
           fontweight="bold", color="#555555")

    # Arrows and shift labels
    for i in range(1, len(curr_bounds)):
        curr_val = curr_bounds[i]
        prop_val = prop_bounds[i]
        if abs(prop_val - curr_val) > 0.005:
            mid_y = (curr_val + prop_val) / 2
            ax.annotate("", xy=(7.0, prop_val), xytext=(3.0, curr_val),
                        arrowprops=dict(arrowstyle="->", color="#CC3311",
                                       linewidth=1.5,
                                       connectionstyle="arc3,rad=0.05"))
            shift = prop_val - curr_val
            ax.text(5.0, mid_y, f"{shift:+.3f}",
                   fontsize=11, ha="center", va="center",
                   color="#CC3311", fontweight="bold",
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
    print("Loading data...")

    # Xenium
    xenium_df = pd.read_csv(os.path.join(CRUMBLR_DIR,
                                          "crumblr_depth_input_subclass.csv"))
    xenium_ctrl = xenium_df[xenium_df["diagnosis"] == "Control"].copy()
    print(f"  Xenium Controls: {xenium_ctrl['donor'].nunique()} donors")

    # MERFISH profiles (pre-computed quantile bins)
    merfish_profiles = pd.read_csv(os.path.join(CRUMBLR_DIR,
                                                 "merfish_depth_profiles_subclass.csv"))
    print(f"  MERFISH profiles: {merfish_profiles['donor'].nunique()} donors")

    # Raw data for boundary computation
    print("  Loading raw MERFISH...")
    merfish_raw = load_merfish_cortical()
    merfish_raw["depth"] = merfish_raw["depth"].clip(0, 1)
    print(f"  MERFISH raw: {len(merfish_raw):,} cells")

    cells_raw = pd.read_csv(os.path.join(DEPTH_DIR, "cell_level_data.csv"))
    cells_ctrl = cells_raw[cells_raw["diagnosis"] == "Control"].copy()
    print(f"  Xenium raw Controls: {len(cells_ctrl):,} cells")

    # ── Compute boundaries ──
    print("\nComputing layer boundaries...")
    mer_bounds = compute_boundaries(merfish_raw, "depth", "subclass")
    xen_bounds = compute_boundaries(cells_ctrl, "predicted_norm_depth",
                                     "subclass_label")

    print("\n  Boundary        MERFISH    Xenium    Current")
    print("  " + "-" * 50)
    current = {"L1/L2-3": 0.10, "L2-3/L4": 0.40, "L4/L5": 0.55,
               "L5/L6": 0.70, "L6/WM": 0.90}
    for name in ["L1/L2-3", "L2-3/L4", "L4/L5", "L5/L6", "L6/WM"]:
        m = mer_bounds.get(name, np.nan)
        x = xen_bounds.get(name, np.nan)
        c = current[name]
        print(f"  {name:12s}  {m:8.4f}  {x:8.4f}  {c:8.2f}")

    # Save proposed boundaries
    rows = []
    for name in ["L1/L2-3", "L2-3/L4", "L4/L5", "L5/L6", "L6/WM"]:
        rows.append({
            "boundary": name,
            "current": current[name],
            "merfish": mer_bounds.get(name, np.nan),
            "xenium": xen_bounds.get(name, np.nan),
        })
    bounds_df = pd.DataFrame(rows)
    outpath = os.path.join(DEPTH_DIR, "proposed_layer_boundaries.csv")
    bounds_df.to_csv(outpath, index=False)
    print(f"\n  Saved: {outpath}")

    # Proposed full layer boundaries (using Xenium values)
    proposed = {
        "L1": (0.0, xen_bounds["L1/L2-3"]),
        "L2/3": (xen_bounds["L1/L2-3"], xen_bounds["L2-3/L4"]),
        "L4": (xen_bounds["L2-3/L4"], xen_bounds["L4/L5"]),
        "L5": (xen_bounds["L4/L5"], xen_bounds["L5/L6"]),
        "L6": (xen_bounds["L5/L6"], xen_bounds["L6/WM"]),
    }
    print("\n  Proposed layer boundaries:")
    for layer, (lo, hi) in proposed.items():
        print(f"    {layer}: {lo:.4f} - {hi:.4f}")

    # ── Common subclasses ──
    common = sorted(set(xenium_ctrl["celltype"]) & set(merfish_profiles["celltype"]))
    glut = sorted([ct for ct in common if infer_class(ct) == "Glut"])
    gaba = sorted([ct for ct in common if infer_class(ct) == "GABA"])
    nn = sorted([ct for ct in common if infer_class(ct) in ("NN", "Other")])

    # ── Profile comparisons ──
    print("\nPlotting neuronal comparison...")
    plot_comparison_profiles(xenium_ctrl, merfish_profiles,
                              glut + gaba, "Neuronal Subclasses",
                              class_filter="neuronal",
                              filename="xenium_vs_merfish_neuronal.png",
                              proposed_bounds=xen_bounds)

    print("Plotting non-neuronal comparison...")
    plot_comparison_profiles(xenium_ctrl, merfish_profiles,
                              nn, "Non-Neuronal Subclasses",
                              class_filter="nonneuronal",
                              filename="xenium_vs_merfish_nonneuronal.png",
                              proposed_bounds=xen_bounds)

    # ── Boundary derivation figure ──
    print("Plotting boundary derivation...")
    plot_boundary_derivation(merfish_raw, cells_ctrl, mer_bounds, xen_bounds)

    print("\nDone!")


if __name__ == "__main__":
    main()
