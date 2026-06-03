"""
Spatial plotting utilities for Xenium cell type visualizations.

Provides rasterized color-blended images for large spatial datasets,
layer zone visualizations, depth prediction panels, and combined
summary figures.

Design principles:
- Large, legible text (fontsize 16+ labels, 20+ titles)
- Dark backgrounds for rasterized spatial plots
- Consistent panel sizes within multi-panel figures
- Always include sample sizes in titles/legends
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgb
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter
from collections import Counter


def build_raster(x, y, colors_per_cell, bin_size=20,
                 bg_color=(0.07, 0.07, 0.15)):
    """
    Build a rasterized RGB image from per-cell RGB color arrays.

    Each pixel is the average color of all cells in that bin.
    Much faster and more visible than scatter plots for large datasets.

    Parameters
    ----------
    x, y : np.ndarray
        Spatial coordinates.
    colors_per_cell : np.ndarray
        (n_cells, 3) RGB array, values in [0, 1].
    bin_size : float
        Pixel size in spatial units (um).
    bg_color : tuple
        RGB background for empty pixels.

    Returns
    -------
    img : np.ndarray (ny, nx, 3)
    extent : list [x_min, x_max, y_max, y_min] for imshow
    """
    x_edges = np.arange(x.min() - bin_size, x.max() + 2 * bin_size, bin_size)
    y_edges = np.arange(y.min() - bin_size, y.max() + 2 * bin_size, bin_size)
    nx = len(x_edges) - 1
    ny = len(y_edges) - 1

    img_r = np.zeros((ny, nx))
    img_g = np.zeros((ny, nx))
    img_b = np.zeros((ny, nx))
    counts = np.zeros((ny, nx))

    xi = np.clip(np.digitize(x, x_edges) - 1, 0, nx - 1)
    yi = np.clip(np.digitize(y, y_edges) - 1, 0, ny - 1)

    for i in range(len(x)):
        r, g, b = colors_per_cell[i]
        img_r[yi[i], xi[i]] += r
        img_g[yi[i], xi[i]] += g
        img_b[yi[i], xi[i]] += b
        counts[yi[i], xi[i]] += 1

    mask = counts > 0
    img_r[mask] /= counts[mask]
    img_g[mask] /= counts[mask]
    img_b[mask] /= counts[mask]
    img_r[~mask] = bg_color[0]
    img_g[~mask] = bg_color[1]
    img_b[~mask] = bg_color[2]

    img = np.stack([img_r, img_g, img_b], axis=-1)
    extent = [x_edges[0], x_edges[-1], y_edges[-1], y_edges[0]]
    return img, extent


def build_color_image(x, y, labels, color_map, bin_size=20,
                      bg_color=(0.07, 0.07, 0.15)):
    """
    Build a rasterized RGB image from spatial cell labels and a color map.

    Parameters
    ----------
    x, y : np.ndarray
        Spatial coordinates.
    labels : np.ndarray
        Cell type label per cell.
    color_map : dict
        {label: hex_color} or {label: (r,g,b)} mapping.
    bin_size : float
        Pixel size in spatial units (um).
    bg_color : tuple
        RGB background for empty pixels.

    Returns
    -------
    img : np.ndarray (ny, nx, 3)
    extent : list
    """
    colors = np.array([
        to_rgb(color_map.get(l, '#CCCCCC')) for l in labels
    ])
    return build_raster(x, y, colors, bin_size=bin_size, bg_color=bg_color)


def plot_spatial_celltype(ax, x, y, labels, color_map, bin_size=20,
                          title="", legend=True, max_legend=15):
    """
    Plot a rasterized spatial cell type map on a matplotlib axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    x, y : np.ndarray
    labels : np.ndarray
    color_map : dict
    bin_size : float
    title : str
    legend : bool
    max_legend : int
    """
    img, extent = build_color_image(x, y, labels, color_map, bin_size=bin_size)
    ax.imshow(img, extent=extent, aspect="equal",
              interpolation="nearest", origin="upper")

    if legend:
        counts = Counter(labels)
        show = sorted(counts.keys(), key=lambda k: -counts[k])[:max_legend]
        handles = [
            Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=to_rgb(color_map.get(s, '#CCC')),
                   markersize=14, label=f'{s} ({counts[s]:,})')
            for s in show
        ]
        fontsize = 14 if len(show) <= 10 else 12 if len(show) <= 15 else 9
        ax.legend(handles=handles, loc="upper right", fontsize=fontsize,
                  framealpha=0.9, markerscale=1.2)

    ax.set_title(title, fontsize=22, fontweight="bold")
    ax.set_xlabel("x (um)", fontsize=16)
    ax.set_ylabel("y (um)", fontsize=16)
    ax.tick_params(labelsize=13)


def plot_summary(adata, output_path, model_bundle=None, dpi=150):
    """
    Generate a 2x3 summary figure for a single annotated sample.

    Row 1: Subclass | Cluster | Class
    Row 2: Predicted depth (continuous) | Discrete layers | Depth histogram

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated sample with 'subclass_label', 'supertype_label',
        'class_label', 'predicted_norm_depth' in .obs,
        spatial coords in .obsm['spatial'].
    output_path : str
        Path to save figure.
    model_bundle : dict, optional
        Depth model bundle (for title annotation).
    dpi : int
        Output resolution.
    """
    from depth_model import LAYER_BINS, LAYER_COLORS

    coords = adata.obsm['spatial']
    x, y = coords[:, 0], coords[:, 1]
    n_cells = adata.shape[0]
    sample_id = adata.obs['sample_id'].iloc[0] if 'sample_id' in adata.obs else ''

    subclass = adata.obs['subclass_label'].values.astype(str)
    cluster = adata.obs['supertype_label'].values.astype(str)
    class_label = adata.obs['class_label'].values.astype(str)

    # Color maps
    palette_main = np.vstack([
        plt.cm.tab20(np.linspace(0, 1, 20)),
        plt.cm.Set3(np.linspace(0, 1, 12))
    ])
    palette_cluster = np.vstack([
        plt.cm.tab20(np.linspace(0, 1, 20)),
        plt.cm.tab20b(np.linspace(0, 1, 20)),
        plt.cm.tab20c(np.linspace(0, 1, 20)),
        plt.cm.Set1(np.linspace(0, 1, 9)),
        plt.cm.Set2(np.linspace(0, 1, 8)),
        plt.cm.Set3(np.linspace(0, 1, 12)),
    ])

    def _make_cmap(labels, pal):
        unique = sorted(set(labels))
        return {u: pal[i % len(pal)][:3] for i, u in enumerate(unique)}

    sub_cmap = _make_cmap(subclass, palette_main)
    clust_cmap = _make_cmap(cluster, palette_cluster)
    class_cmap = {
        'Neuronal: Glutamatergic': (0.2, 0.6, 0.9),
        'Neuronal: GABAergic': (0.9, 0.3, 0.2),
        'Non-neuronal and Non-neural': (0.3, 0.8, 0.3),
        'Glutamatergic': (0.2, 0.6, 0.9),
        'GABAergic': (0.9, 0.3, 0.2),
        'Non-neuronal': (0.3, 0.8, 0.3),
    }

    fig, axes = plt.subplots(2, 3, figsize=(54, 36), facecolor='white')

    # Row 1: Cell type annotations
    # Subclass
    colors_sub = np.array([sub_cmap.get(s, (0.5, 0.5, 0.5)) for s in subclass])
    img, ext = build_raster(x, y, colors_sub)
    axes[0, 0].imshow(img, extent=ext, aspect='equal', interpolation='nearest', origin='upper')
    counts_sub = Counter(subclass)
    top = sorted(counts_sub.keys(), key=lambda k: -counts_sub[k])[:15]
    handles = [Line2D([0], [0], marker='o', color='w',
                      markerfacecolor=sub_cmap.get(s, (0.5, 0.5, 0.5)),
                      markersize=14, label=f'{s} ({counts_sub[s]:,})') for s in top]
    axes[0, 0].legend(handles=handles, loc='upper right', fontsize=14, framealpha=0.9)
    axes[0, 0].set_title('Subclass Labels', fontsize=22, fontweight='bold')
    axes[0, 0].set_xlabel('x (um)', fontsize=16)
    axes[0, 0].set_ylabel('y (um)', fontsize=16)
    axes[0, 0].tick_params(labelsize=13)

    # Cluster
    colors_cl = np.array([clust_cmap.get(c, (0.5, 0.5, 0.5)) for c in cluster])
    img, ext = build_raster(x, y, colors_cl)
    axes[0, 1].imshow(img, extent=ext, aspect='equal', interpolation='nearest', origin='upper')
    counts_cl = Counter(cluster)
    top_cl = sorted(counts_cl.keys(), key=lambda k: -counts_cl[k])[:15]
    handles = [Line2D([0], [0], marker='o', color='w',
                      markerfacecolor=clust_cmap.get(c, (0.5, 0.5, 0.5)),
                      markersize=14, label=f'{c} ({counts_cl[c]:,})') for c in top_cl]
    axes[0, 1].legend(handles=handles, loc='upper right', fontsize=12, framealpha=0.9)
    axes[0, 1].set_title('Cluster Labels', fontsize=22, fontweight='bold')
    axes[0, 1].set_xlabel('x (um)', fontsize=16)
    axes[0, 1].set_ylabel('y (um)', fontsize=16)
    axes[0, 1].tick_params(labelsize=13)

    # Class
    colors_cls = np.array([class_cmap.get(c, (0.5, 0.5, 0.5)) for c in class_label])
    img, ext = build_raster(x, y, colors_cls)
    axes[0, 2].imshow(img, extent=ext, aspect='equal', interpolation='nearest', origin='upper')
    counts_cls = Counter(class_label)
    handles = [Line2D([0], [0], marker='o', color='w',
                      markerfacecolor=class_cmap.get(c, (0.5, 0.5, 0.5)),
                      markersize=16, label=f'{c} ({counts_cls[c]:,})')
               for c in class_cmap if c in counts_cls]
    axes[0, 2].legend(handles=handles, loc='upper right', fontsize=16, framealpha=0.9)
    axes[0, 2].set_title('Class Labels', fontsize=22, fontweight='bold')
    axes[0, 2].set_xlabel('x (um)', fontsize=16)
    axes[0, 2].set_ylabel('y (um)', fontsize=16)
    axes[0, 2].tick_params(labelsize=13)

    # Row 2: Depth predictions
    if 'predicted_norm_depth' in adata.obs.columns:
        pred_depth = adata.obs['predicted_norm_depth'].values.astype(float)

        # Continuous depth
        cmap_d = plt.cm.viridis
        vmin, vmax = pred_depth.min(), pred_depth.max()
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        colors_d = np.array([cmap_d(norm(d))[:3] for d in pred_depth])
        img, ext = build_raster(x, y, colors_d)
        axes[1, 0].imshow(img, extent=ext, aspect='equal', interpolation='nearest', origin='upper')
        sm = plt.cm.ScalarMappable(cmap=cmap_d, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=axes[1, 0], shrink=0.6, pad=0.02)
        cbar.set_label('Predicted normalized depth\n(0 = pia, 1 = WM)', fontsize=16)
        cbar.ax.tick_params(labelsize=14)
        axes[1, 0].set_title(f'Predicted Normalized Depth\n(range: {vmin:.3f} to {vmax:.3f})',
                             fontsize=22, fontweight='bold')
        axes[1, 0].set_xlabel('x (um)', fontsize=16)
        axes[1, 0].set_ylabel('y (um)', fontsize=16)
        axes[1, 0].tick_params(labelsize=13)

        # Discrete layers
        colors_lay = np.zeros((n_cells, 3))
        cell_layers = np.full(n_cells, '', dtype=object)
        for lname, (lo, hi) in LAYER_BINS.items():
            mask = (pred_depth >= lo) & (pred_depth < hi)
            colors_lay[mask] = LAYER_COLORS[lname]
            cell_layers[mask] = lname
        img, ext = build_raster(x, y, colors_lay)
        axes[1, 1].imshow(img, extent=ext, aspect='equal', interpolation='nearest', origin='upper')
        layer_counts = Counter(cell_layers)
        handles = [Line2D([0], [0], marker='o', color='w',
                          markerfacecolor=LAYER_COLORS[l], markersize=16,
                          label=f'{l} ({layer_counts.get(l, 0):,})')
                   for l in LAYER_BINS]
        axes[1, 1].legend(handles=handles, loc='upper right', fontsize=16, framealpha=0.9)
        axes[1, 1].set_title('Discrete Layers\n(from predicted depth)', fontsize=22, fontweight='bold')
        axes[1, 1].set_xlabel('x (um)', fontsize=16)
        axes[1, 1].set_ylabel('y (um)', fontsize=16)
        axes[1, 1].tick_params(labelsize=13)

        # Depth histogram
        ax_h = axes[1, 2]
        ax_h.set_facecolor('white')
        ax_h.hist(pred_depth, bins=80, color='steelblue', edgecolor='none', alpha=0.8)
        boundaries = [0.1225, 0.4696, 0.5443, 0.7079, 0.9275]
        for b in boundaries:
            ax_h.axvline(b, color='red', linestyle='--', linewidth=2, alpha=0.8)
        for lname, (lo, hi) in LAYER_BINS.items():
            lo_c = max(lo, pred_depth.min() - 0.05)
            hi_c = min(hi, pred_depth.max() + 0.05)
            ax_h.axvspan(lo_c, hi_c, alpha=0.15, color=LAYER_COLORS[lname])
            mid = (max(lo, pred_depth.min()) + min(hi, pred_depth.max())) / 2
            ax_h.text(mid, ax_h.get_ylim()[1] * 0.5, lname,
                     ha='center', va='center', fontsize=18,
                     fontweight='bold', color=LAYER_COLORS[lname], alpha=0.7)
        ax_h.set_xlabel('Predicted normalized depth', fontsize=18)
        ax_h.set_ylabel('Number of cells', fontsize=18)
        ax_h.set_title(f'Depth Distribution\n(n={n_cells:,} cells)', fontsize=22, fontweight='bold')
        ax_h.tick_params(labelsize=14)
    else:
        for j in range(3):
            axes[1, j].text(0.5, 0.5, 'No depth predictions',
                           transform=axes[1, j].transAxes,
                           ha='center', va='center', fontsize=20)

    fig.suptitle(f'{sample_id} — Xenium Spatial Transcriptomics Summary\n'
                 f'n = {n_cells:,} cells | SEA-AD label transfer | MERFISH-trained depth model',
                 fontsize=28, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_path}")


# Keep legacy functions for backwards compatibility
def plot_layer_zones(ax, layer_result, x_range, y_range, title="Cortical layers"):
    """Plot dominant cortical layer zones (density-based, legacy)."""
    from layers import LAYER_NAMES, LAYER_COLORS as LEGACY_COLORS
    densities = layer_result["densities"]
    contours = layer_result["contours"]
    dominant = densities["dominant"]
    x_centers = densities["x_centers"]
    y_centers = densities["y_centers"]

    cmap_colors = ["#000000"] + [LEGACY_COLORS[ln] for ln in LAYER_NAMES]
    cmap = ListedColormap(cmap_colors)
    ax.pcolormesh(x_centers, y_centers, dominant,
                  cmap=cmap, vmin=-1, vmax=len(LAYER_NAMES) - 1, alpha=0.7)
    for ln, contour_list in contours.items():
        for cx, cy in contour_list:
            ax.plot(cx, cy, color="black", linewidth=1.5, alpha=0.7)
    ax.set_xlim(x_range[0] - 100, x_range[1] + 100)
    ax.set_ylim(y_range[1] + 100, y_range[0] - 100)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=20, fontweight="bold")
    ax.set_xlabel("x (um)", fontsize=16)
    ax.set_ylabel("y (um)", fontsize=16)


def plot_representative(adata, colors, output_path, dpi=200):
    """Legacy 3-panel figure (subclass | cluster | layer zones)."""
    from layers import segment_layers
    coords = adata.obsm["spatial"]
    x, y = coords[:, 0], coords[:, 1]
    subclass = adata.obs["subclass_label"].values
    cluster = adata.obs["supertype_label"].values
    sample_id = adata.obs["sample_id"].iloc[0] if "sample_id" in adata.obs else ""
    subclass_colors = colors.get("subclass_label", {})
    cluster_colors = colors.get("supertype_label", {})
    layer_result = segment_layers(adata)

    fig, axes = plt.subplots(1, 3, figsize=(48, 16), facecolor="white")
    plot_spatial_celltype(axes[0], x, y, subclass, subclass_colors, title="Subclass")
    plot_spatial_celltype(axes[1], x, y, cluster, cluster_colors, title="Cluster")
    plot_layer_zones(axes[2], layer_result,
                     x_range=(x.min(), x.max()), y_range=(y.min(), y.max()))
    fig.suptitle(f"{sample_id} — Xenium spatial cell type map (SEA-AD label transfer)",
                 fontsize=22, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")
