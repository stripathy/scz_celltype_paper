#!/usr/bin/env python3
"""
4-panel depth comparison: MERFISH vs Xenium.

Row 1: SEA-AD MERFISH section -- predicted depth (left), predicted layer (right)
       with outline showing the region with manual depth annotations
Row 2: Xenium Br8667 -- predicted depth (left), predicted layer (right)

Illustrates how depth labels were transferred from MERFISH (sparse, manual
annotations in a small region) to Xenium (dense, predicted for entire section).

Layer assignment: depth bins for all cells, with Vascular cells identified
via OOD spatial domain classification (Xenium only).

Output: output/presentation/slide_depth_comparison.png
"""

import os
import sys
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from scipy.spatial import ConvexHull

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BG_COLOR, H5AD_DIR, MERFISH_PATH, PRESENTATION_DIR,
    LAYER_COLORS, LAYER_ORDER, SUBCLASS_CONF_THRESH,
)

OUT_DIR = PRESENTATION_DIR

MERFISH_SECTION = "H21.33.031.CX24.MTG.02.007.1.01.01"
XENIUM_SAMPLE = "Br8667"

BG = BG_COLOR

# Depth colormap: viridis (0=pia/surface, 1=WM/deep)
DEPTH_CMAP = cm.viridis


def load_merfish_section():
    """Load MERFISH section with spatial coords, manual + predicted depth/layer."""
    print(f"Loading MERFISH section {MERFISH_SECTION}...")
    merfish = ad.read_h5ad(MERFISH_PATH, backed="r")
    obs = merfish.obs
    mask = obs["Section"] == MERFISH_SECTION
    sec = obs[mask].copy()

    coords = merfish.obsm["spatial"][mask.values]
    sec["x"] = coords[:, 0]
    sec["y"] = coords[:, 1]

    sec["manual_depth"] = sec["Normalized depth from pia"].astype(float)
    sec["pred_depth"] = sec["predicted_norm_depth"].astype(float)
    sec["pred_layer"] = sec["predicted_layer"].astype(str).fillna("")
    sec["manual_layer"] = sec["Layer annotation"].astype(str)
    sec.loc[sec["manual_layer"] == "nan", "manual_layer"] = ""

    # Build composite: manual where available, predicted elsewhere
    has_manual = ~np.isnan(sec["manual_depth"].values)
    composite_depth = sec["pred_depth"].values.copy()
    composite_depth[has_manual] = sec["manual_depth"].values[has_manual]

    composite_layer = sec["pred_layer"].values.copy()
    manual_layer_vals = sec["manual_layer"].values
    has_manual_layer = manual_layer_vals != ""
    composite_layer[has_manual_layer] = manual_layer_vals[has_manual_layer]

    sec["composite_depth"] = composite_depth
    sec["composite_layer"] = composite_layer
    sec["has_manual"] = has_manual
    return sec


def load_xenium_sample():
    """Load Xenium sample with spatial coords, predicted depth, and layer."""
    print(f"Loading Xenium {XENIUM_SAMPLE}...")
    fpath = os.path.join(H5AD_DIR, f"{XENIUM_SAMPLE}_annotated.h5ad")
    adata = ad.read_h5ad(fpath, backed="r")

    has_corr = "corr_qc_pass" in adata.obs.columns
    cols = ["sample_id", "predicted_norm_depth", "layer",
            "spatial_domain", "qc_pass"]
    if has_corr:
        cols.append("corr_qc_pass")
    else:
        cols.append("subclass_label_confidence")
    obs = adata.obs[cols].copy()
    coords = adata.obsm["spatial"]
    obs["x"] = coords[:, 0]
    obs["y"] = coords[:, 1]
    obs = obs[obs["qc_pass"] == True].copy()

    # Correlation QC or legacy confidence filter
    if has_corr:
        obs = obs[obs["corr_qc_pass"] == True]
    else:
        obs = obs[obs["subclass_label_confidence"].astype(float) >= SUBCLASS_CONF_THRESH]
    obs["depth"] = obs["predicted_norm_depth"].astype(float)
    return obs


def draw_manual_region_outline(ax, x, y, has_manual, color="#ffffff",
                                linewidth=2.0, alpha=0.8):
    """Draw a convex hull outline around cells with manual annotations."""
    pts = np.column_stack([x[has_manual], y[has_manual]])
    if len(pts) < 3:
        return
    hull = ConvexHull(pts)
    hull_pts = pts[hull.vertices]
    # Close the polygon
    hull_pts = np.vstack([hull_pts, hull_pts[0]])
    ax.plot(hull_pts[:, 0], hull_pts[:, 1], color=color,
            linewidth=linewidth, alpha=alpha, linestyle="--", zorder=10)


def plot_depth_panel(ax, x, y, depth, title, has_depth_mask=None, s=0.8,
                     manual_mask=None):
    """Plot cells colored by continuous depth (viridis).

    If manual_mask is provided, draws an outline around the manually
    annotated region.
    """
    if has_depth_mask is None:
        has_depth_mask = ~np.isnan(depth)

    no_depth = ~has_depth_mask
    if no_depth.sum() > 0:
        ax.scatter(x[no_depth], y[no_depth], s=s * 0.3, c="#333333",
                   alpha=0.15, rasterized=True, linewidths=0, zorder=1)

    # Plot cells with depth colored by viridis
    if has_depth_mask.sum() > 0:
        norm = Normalize(vmin=0, vmax=1)
        colors = DEPTH_CMAP(norm(np.clip(depth[has_depth_mask], 0, 1)))
        ax.scatter(x[has_depth_mask], y[has_depth_mask], s=s,
                   c=colors, alpha=0.7, rasterized=True, linewidths=0,
                   zorder=3)

    # Draw outline around manually annotated region
    if manual_mask is not None:
        draw_manual_region_outline(ax, x, y, manual_mask)

    n_with = has_depth_mask.sum()
    n_total = len(x)
    ax.text(0.02, 0.98,
            f"{n_with:,} / {n_total:,} cells\nwith depth ({n_with/n_total*100:.0f}%)",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=13, color="#dddddd",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#333333",
                      edgecolor="#555555", alpha=0.85))

    ax.set_title(title, fontsize=18, fontweight="bold", color="white", pad=6)
    ax.set_facecolor(BG)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_layer_panel(ax, x, y, layers, title, s=0.8, manual_mask=None):
    """Plot cells colored by discrete layer assignment."""
    # Plot cells without layer as dim gray
    no_layer = (layers == "") | layers.isna() if hasattr(layers, "isna") else (layers == "")
    has_layer = ~no_layer

    if isinstance(no_layer, np.ndarray):
        pass
    else:
        no_layer = no_layer.values
        has_layer = has_layer.values

    if no_layer.sum() > 0:
        ax.scatter(x[no_layer], y[no_layer], s=s * 0.3, c="#333333",
                   alpha=0.15, rasterized=True, linewidths=0, zorder=1)

    # Plot each layer
    layer_vals = layers.values if hasattr(layers, "values") else layers
    for lname, color in LAYER_COLORS.items():
        mask = (layer_vals == lname) & has_layer
        if mask.sum() > 0:
            ax.scatter(x[mask], y[mask], s=s, c=[color],
                       alpha=0.7, rasterized=True, linewidths=0, zorder=3)

    # Draw outline around manually annotated region
    if manual_mask is not None:
        draw_manual_region_outline(ax, x, y, manual_mask)

    n_with = has_layer.sum()
    n_total = len(x)
    ax.text(0.02, 0.98,
            f"{n_with:,} / {n_total:,} cells\nwith layer ({n_with/n_total*100:.0f}%)",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=13, color="#dddddd",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#333333",
                      edgecolor="#555555", alpha=0.85))

    ax.set_title(title, fontsize=18, fontweight="bold", color="white", pad=6)
    ax.set_facecolor(BG)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def main():
    merfish = load_merfish_section()
    xenium = load_xenium_sample()

    # 2 rows x 2 cols: depth | layer
    fig = plt.figure(figsize=(20, 18), facecolor=BG)
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           wspace=0.05, hspace=0.18,
                           left=0.02, right=0.90, top=0.90, bottom=0.06)

    # --- Row 1: MERFISH (manual where available, predicted elsewhere) ---
    mx = merfish["x"].values
    my = merfish["y"].values
    m_depth = merfish["composite_depth"].values
    m_has_depth = ~np.isnan(m_depth)
    m_layer = merfish["composite_layer"]
    m_has_manual = merfish["has_manual"].values

    ax1 = fig.add_subplot(gs[0, 0])
    plot_depth_panel(ax1, mx, my, m_depth,
                     "Depth (manual inside outline, predicted outside)",
                     has_depth_mask=m_has_depth, s=1.0,
                     manual_mask=m_has_manual)

    ax2 = fig.add_subplot(gs[0, 1])
    plot_layer_panel(ax2, mx, my, m_layer,
                     "Layer (manual inside outline, predicted outside)",
                     s=1.0, manual_mask=m_has_manual)

    # --- Row 2: Xenium ---
    xx = xenium["x"].values
    xy = xenium["y"].values
    x_depth = xenium["depth"].values
    x_has_depth = ~np.isnan(x_depth)

    ax3 = fig.add_subplot(gs[1, 0])
    plot_depth_panel(ax3, xx, xy, x_depth, "Predicted depth (all cells)",
                     has_depth_mask=x_has_depth, s=1.0)

    x_layer = xenium["layer"]
    ax4 = fig.add_subplot(gs[1, 1])
    plot_layer_panel(ax4, xx, xy, x_layer,
                     "Predicted layer (depth bins + Vascular OOD)", s=1.0)

    # Row labels
    n_manual = m_has_manual.sum()
    n_total_m = len(mx)
    fig.text(0.01, 0.92,
             f"SEA-AD MERFISH (reference)  --  "
             f"manual annotations inside outline ({n_manual:,}/{n_total_m:,}, "
             f"{n_manual/n_total_m*100:.0f}%), predicted outside",
             fontsize=16, fontweight="bold", color="#dddddd", ha="left")
    fig.text(0.01, 0.48, f"Xenium ({XENIUM_SAMPLE}, Control)", fontsize=20,
             fontweight="bold", color="#dddddd", ha="left")

    # Shared colorbar for depth (right side) -- inverted so pia (0) is at top
    cbar_ax = fig.add_axes([0.91, 0.52, 0.015, 0.36])
    norm = Normalize(vmin=0, vmax=1)
    sm = cm.ScalarMappable(cmap=DEPTH_CMAP, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=cbar_ax)
    cbar.ax.invert_yaxis()
    cbar.set_label("Normalized depth\n(0 = pia, 1 = WM)", fontsize=13,
                    color="white", labelpad=10)
    cbar.ax.tick_params(colors="white", labelsize=11)

    # Shared layer legend (bottom right)
    legend_elements = [
        Line2D([0], [0], marker="o", color=BG,
               markerfacecolor=LAYER_COLORS[l], markersize=12,
               label=l, linewidth=0)
        for l in LAYER_ORDER
    ]
    legend_elements.append(
        Line2D([0], [0], color="white", linewidth=2, linestyle="--",
               label="Manual annotation region", alpha=0.8)
    )
    fig.legend(handles=legend_elements, loc="lower right",
               fontsize=13, frameon=False, labelcolor="white",
               bbox_to_anchor=(0.98, 0.02), ncol=1)

    outpath = os.path.join(OUT_DIR, "slide_depth_comparison.png")
    plt.savefig(outpath, dpi=150, facecolor=BG)
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == "__main__":
    main()
