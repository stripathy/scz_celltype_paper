#!/usr/bin/env python3
"""
Visualize spatial doublet cells with mixed marker expression.

Creates publication-quality figures showing:
1. Marker expression heatmap: individual doublet cells vs normal Glut/GABA cells
2. Spatial zoom: doublet cells in spatial context with transcript dots
3. Total UMI comparison: doublets have ~2x counts (evidence of two cells)

Output: output/presentation/doublet_*.png
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BG_COLOR, H5AD_DIR, PRESENTATION_DIR

# ── Marker genes for doublet detection ──
GABA_MARKERS = ['GAD1', 'GAD2', 'SLC32A1', 'SST', 'PVALB', 'VIP', 'LAMP5']
GLUT_MARKERS = ['CUX2', 'RORB', 'GRIN2A', 'THEMIS']
GABA_TRIPLE = ['SST', 'PVALB', 'LAMP5']
ALL_MARKERS = GABA_MARKERS + GLUT_MARKERS

BG = BG_COLOR
SAMPLE = 'Br8667'  # Default viewer sample, good doublet diversity

OUT_DIR = PRESENTATION_DIR


def load_sample_data(sample_id):
    """Load sample h5ad and extract doublet info + marker expression."""
    fpath = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    print(f"Loading {fpath}...")
    adata = ad.read_h5ad(fpath)

    # Get marker expression (raw counts)
    X = adata.X
    if sparse.issparse(X):
        X = X.toarray()

    gene_names = list(adata.var_names)
    marker_idx = {g: gene_names.index(g) for g in ALL_MARKERS if g in gene_names}

    # Build marker expression DataFrame
    marker_expr = pd.DataFrame(
        {g: X[:, idx] for g, idx in marker_idx.items()},
        index=adata.obs_names
    )

    # Attach metadata
    marker_expr['doublet_suspect'] = adata.obs['doublet_suspect'].values
    marker_expr['doublet_type'] = adata.obs['doublet_type'].astype(str).values
    marker_expr['corr_class'] = adata.obs['corr_class'].astype(str).values
    marker_expr['corr_subclass'] = adata.obs['corr_subclass'].astype(str).values
    marker_expr['total_counts'] = adata.obs['total_counts'].values
    marker_expr['x'] = adata.obs['x_centroid'].values if 'x_centroid' in adata.obs else adata.obsm['spatial'][:, 0]
    marker_expr['y'] = adata.obs['y_centroid'].values if 'y_centroid' in adata.obs else adata.obsm['spatial'][:, 1]
    marker_expr['qc_pass'] = adata.obs['qc_pass'].values

    return marker_expr, adata


def plot_marker_heatmap(marker_expr, outpath):
    """
    Figure 1: Heatmap showing marker expression in individual doublet cells
    vs randomly sampled normal Glutamatergic and GABAergic cells.
    """
    # Get cells
    glut_gaba = marker_expr[
        (marker_expr['doublet_type'] == 'Glut+GABA') &
        (marker_expr['qc_pass'] == True)
    ]
    gaba_gaba = marker_expr[
        (marker_expr['doublet_type'] == 'GABA+GABA') &
        (marker_expr['qc_pass'] == True)
    ]
    normal_glut = marker_expr[
        (marker_expr['corr_class'] == 'Glutamatergic') &
        (marker_expr['doublet_suspect'] == False) &
        (marker_expr['qc_pass'] == True)
    ]
    normal_gaba = marker_expr[
        (marker_expr['corr_class'] == 'GABAergic') &
        (marker_expr['doublet_suspect'] == False) &
        (marker_expr['qc_pass'] == True)
    ]

    print(f"  Glut+GABA doublets: {len(glut_gaba)}")
    print(f"  GABA+GABA doublets: {len(gaba_gaba)}")
    print(f"  Normal Glut: {len(normal_glut)}")
    print(f"  Normal GABA: {len(normal_gaba)}")

    # Sample cells for display
    np.random.seed(42)
    n_show = 15  # cells per category
    n_doublet_show = min(20, len(glut_gaba))
    n_gaba_doublet_show = min(10, len(gaba_gaba))

    sampled_glut = normal_glut.sample(n=n_show, random_state=42)
    sampled_gaba = normal_gaba.sample(n=n_show, random_state=42)
    sampled_gg = glut_gaba.sample(n=n_doublet_show, random_state=42)
    sampled_gabagaba = gaba_gaba.sample(n=n_gaba_doublet_show, random_state=42)

    # Combine in order
    display_cells = pd.concat([sampled_glut, sampled_gg, sampled_gaba, sampled_gabagaba])
    labels = (
        ['Normal Glut'] * n_show +
        ['Glut+GABA\nDoublet'] * n_doublet_show +
        ['Normal GABA'] * n_show +
        ['GABA+GABA\nDoublet'] * n_gaba_doublet_show
    )

    # Extract marker columns only
    marker_genes = [g for g in ALL_MARKERS if g in marker_expr.columns]
    heatmap_data = display_cells[marker_genes].values.astype(float)

    # Log-transform for visualization: log2(counts + 1)
    heatmap_log = np.log2(heatmap_data + 1)

    # ── Plot ──
    fig = plt.figure(figsize=(14, 10), facecolor=BG)
    gs = gridspec.GridSpec(1, 2, width_ratios=[30, 1], wspace=0.03)

    ax = fig.add_subplot(gs[0])
    cax = fig.add_subplot(gs[1])

    # Custom colormap: dark background → blue → yellow
    cmap = mcolors.LinearSegmentedColormap.from_list(
        'expr', ['#1a1a2e', '#0f3460', '#16a085', '#f1c40f', '#e74c3c'], N=256
    )

    im = ax.imshow(heatmap_log, aspect='auto', cmap=cmap, interpolation='nearest')

    # Gene labels (x-axis)
    ax.set_xticks(range(len(marker_genes)))
    ax.set_xticklabels(marker_genes, fontsize=13, color='white', rotation=45, ha='right',
                       fontweight='bold')

    # Category separators and labels
    boundaries = [0, n_show, n_show + n_doublet_show,
                  n_show + n_doublet_show + n_show,
                  n_show + n_doublet_show + n_show + n_gaba_doublet_show]

    cat_labels = ['Normal\nGlutamatergic', 'Glut+GABA\nDoublets',
                  'Normal\nGABAergic', 'GABA+GABA\nDoublets']
    cat_colors = ['#00ADF8', '#FF6B6B', '#F05A28', '#FF6BFF']

    for i in range(len(boundaries) - 1):
        mid = (boundaries[i] + boundaries[i + 1]) / 2
        ax.text(-0.7, mid, cat_labels[i], fontsize=12, color=cat_colors[i],
                ha='right', va='center', fontweight='bold',
                transform=ax.get_yaxis_transform())

    # Draw separator lines
    for b in boundaries[1:-1]:
        ax.axhline(y=b - 0.5, color='white', linewidth=2, alpha=0.6)

    # Draw vertical separator between GABA and Glut markers
    gaba_end = len(GABA_MARKERS) - 0.5
    ax.axvline(x=gaba_end, color='#888888', linewidth=1.5, linestyle='--', alpha=0.6)
    ax.text(gaba_end / 2, -1.5, 'GABAergic markers', fontsize=11, color='#F05A28',
            ha='center', va='bottom', fontweight='bold')
    ax.text(gaba_end + (len(GLUT_MARKERS)) / 2, -1.5, 'Glut markers', fontsize=11,
            color='#00ADF8', ha='center', va='bottom', fontweight='bold')

    ax.set_yticks([])
    ax.set_title(f'Marker Expression in Spatial Doublets ({SAMPLE})',
                 fontsize=20, fontweight='bold', color='white', pad=25)
    ax.set_facecolor(BG)
    ax.tick_params(colors='white')
    for spine in ax.spines.values():
        spine.set_color('#555555')

    # Colorbar
    cb = plt.colorbar(im, cax=cax)
    cb.set_label('log₂(counts + 1)', fontsize=12, color='white')
    cb.ax.yaxis.set_tick_params(color='white')
    cb.ax.yaxis.set_ticklabels(cb.ax.yaxis.get_ticklabels(), color='white', fontsize=10)

    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


def plot_marker_barplot(marker_expr, outpath):
    """
    Figure 2: Mean marker expression (% cells expressing) across cell categories.
    Clean bar chart showing the mixed expression signature.
    """
    marker_genes = [g for g in ALL_MARKERS if g in marker_expr.columns]

    categories = {
        'Normal Glut': marker_expr[
            (marker_expr['corr_class'] == 'Glutamatergic') &
            (marker_expr['doublet_suspect'] == False) &
            (marker_expr['qc_pass'] == True)
        ],
        'Glut+GABA Doublets': marker_expr[
            (marker_expr['doublet_type'] == 'Glut+GABA') &
            (marker_expr['qc_pass'] == True)
        ],
        'Normal GABA': marker_expr[
            (marker_expr['corr_class'] == 'GABAergic') &
            (marker_expr['doublet_suspect'] == False) &
            (marker_expr['qc_pass'] == True)
        ],
        'GABA+GABA Doublets': marker_expr[
            (marker_expr['doublet_type'] == 'GABA+GABA') &
            (marker_expr['qc_pass'] == True)
        ],
    }

    # Compute % cells expressing each marker (count > 0)
    pct_data = {}
    for cat_name, cat_df in categories.items():
        pcts = []
        for g in marker_genes:
            pct = (cat_df[g] > 0).mean() * 100
            pcts.append(pct)
        pct_data[cat_name] = pcts

    # ── Plot ──
    fig, ax = plt.subplots(figsize=(16, 7), facecolor=BG)

    x = np.arange(len(marker_genes))
    width = 0.2
    colors = ['#00ADF8', '#FF6B6B', '#F05A28', '#FF6BFF']

    for i, (cat_name, pcts) in enumerate(pct_data.items()):
        offset = (i - 1.5) * width
        bars = ax.bar(x + offset, pcts, width, label=cat_name,
                      color=colors[i], alpha=0.85, edgecolor='white', linewidth=0.5)

    # Vertical separator between GABA and Glut markers
    gaba_end = len(GABA_MARKERS) - 0.5
    ax.axvline(x=gaba_end, color='#888888', linewidth=2, linestyle='--', alpha=0.5)
    ax.text(gaba_end / 2, ax.get_ylim()[1] * 0.95, 'GABAergic markers',
            fontsize=13, color='#F05A28', ha='center', va='top', fontweight='bold')
    ax.text(gaba_end + (len(GLUT_MARKERS)) / 2, ax.get_ylim()[1] * 0.95,
            'Glutamatergic markers', fontsize=13, color='#00ADF8',
            ha='center', va='top', fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(marker_genes, fontsize=14, color='white', rotation=45,
                       ha='right', fontweight='bold')
    ax.set_ylabel('% cells expressing marker', fontsize=16, color='white')
    ax.set_title(f'Marker Detection Rates: Doublets vs Normal Cells ({SAMPLE})',
                 fontsize=20, fontweight='bold', color='white', pad=15)

    ax.legend(fontsize=13, loc='upper right', facecolor='#333333',
              edgecolor='#555555', labelcolor='white')

    ax.set_facecolor(BG)
    ax.tick_params(colors='white', labelsize=12)
    for spine in ax.spines.values():
        spine.set_color('#555555')
    ax.grid(axis='y', alpha=0.2, color='#555555')
    ax.set_ylim(0, 105)

    plt.tight_layout()
    plt.savefig(outpath, dpi=200, facecolor=BG)
    plt.close()
    print(f"  Saved: {outpath}")


def plot_umi_comparison(marker_expr, outpath):
    """
    Figure 3: Total UMI count distributions for doublets vs normal cells.
    Shows doublets have ~2x counts (evidence of two captured cells).
    """
    categories = {
        'Normal Glut': marker_expr[
            (marker_expr['corr_class'] == 'Glutamatergic') &
            (marker_expr['doublet_suspect'] == False) &
            (marker_expr['qc_pass'] == True)
        ]['total_counts'],
        'Glut+GABA\nDoublets': marker_expr[
            (marker_expr['doublet_type'] == 'Glut+GABA') &
            (marker_expr['qc_pass'] == True)
        ]['total_counts'],
        'Normal GABA': marker_expr[
            (marker_expr['corr_class'] == 'GABAergic') &
            (marker_expr['doublet_suspect'] == False) &
            (marker_expr['qc_pass'] == True)
        ]['total_counts'],
        'GABA+GABA\nDoublets': marker_expr[
            (marker_expr['doublet_type'] == 'GABA+GABA') &
            (marker_expr['qc_pass'] == True)
        ]['total_counts'],
    }

    fig, ax = plt.subplots(figsize=(10, 7), facecolor=BG)

    colors = ['#00ADF8', '#FF6B6B', '#F05A28', '#FF6BFF']
    positions = range(len(categories))
    labels = list(categories.keys())

    bp = ax.boxplot(
        [v.values for v in categories.values()],
        positions=positions,
        widths=0.6,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color='white', linewidth=2),
        whiskerprops=dict(color='#888888'),
        capprops=dict(color='#888888'),
    )

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
        patch.set_edgecolor('white')

    # Add median labels
    for i, (cat_name, vals) in enumerate(categories.items()):
        median = vals.median()
        mean = vals.mean()
        n = len(vals)
        ax.text(i, median + ax.get_ylim()[1] * 0.02,
                f'median={median:.0f}\nn={n:,}',
                ha='center', va='bottom', fontsize=11, color='white',
                fontweight='bold')

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=14, color='white', fontweight='bold')
    ax.set_ylabel('Total UMI counts per cell', fontsize=16, color='white')
    ax.set_title(f'Doublets Have Elevated UMI Counts ({SAMPLE})',
                 fontsize=20, fontweight='bold', color='white', pad=15)

    ax.set_facecolor(BG)
    ax.tick_params(colors='white', labelsize=12)
    for spine in ax.spines.values():
        spine.set_color('#555555')
    ax.grid(axis='y', alpha=0.2, color='#555555')

    plt.tight_layout()
    plt.savefig(outpath, dpi=200, facecolor=BG)
    plt.close()
    print(f"  Saved: {outpath}")


def plot_spatial_doublet_zoom(marker_expr, outpath):
    """
    Figure 4: Spatial zoom showing doublet cells in their tissue context.
    Highlights doublet cells among normal cells in a small region.
    """
    # Find a good region with multiple doublets
    doublets = marker_expr[marker_expr['doublet_suspect'] == True]
    if len(doublets) == 0:
        print("  No doublets found, skipping spatial plot")
        return

    # Find a cluster of doublets for the zoom
    # Pick a doublet and show its neighborhood
    np.random.seed(42)

    # Find a region with multiple doublets close together
    from scipy.spatial import KDTree

    doublet_coords = doublets[['x', 'y']].values
    tree = KDTree(doublet_coords)

    # For each doublet, count how many other doublets are within 200 microns
    best_idx = 0
    best_count = 0
    for i in range(min(500, len(doublet_coords))):
        neighbors = tree.query_ball_point(doublet_coords[i], r=200)
        if len(neighbors) > best_count:
            best_count = len(neighbors)
            best_idx = i

    center_x = doublet_coords[best_idx, 0]
    center_y = doublet_coords[best_idx, 1]
    radius = 300  # microns

    # Get all cells in this region
    region = marker_expr[
        (marker_expr['x'] > center_x - radius) &
        (marker_expr['x'] < center_x + radius) &
        (marker_expr['y'] > center_y - radius) &
        (marker_expr['y'] < center_y + radius) &
        (marker_expr['qc_pass'] == True)
    ]

    print(f"  Spatial zoom: {len(region)} cells in {2*radius}µm region, "
          f"{region['doublet_suspect'].sum()} doublets")

    # ── Plot ──
    fig, axes = plt.subplots(1, 3, figsize=(24, 8), facecolor=BG)

    # Panel A: Colored by class
    class_colors = {
        'Glutamatergic': '#00ADF8',
        'GABAergic': '#F05A28',
        'Non-neuronal': '#808080',
    }

    ax = axes[0]
    for cls, color in class_colors.items():
        cells = region[(region['corr_class'] == cls) & (region['doublet_suspect'] == False)]
        ax.scatter(cells['x'], cells['y'], c=color, s=25, alpha=0.5,
                   edgecolors='none', label=cls)

    # Overlay doublets with distinctive markers
    glut_gaba = region[region['doublet_type'] == 'Glut+GABA']
    gaba_gaba = region[region['doublet_type'] == 'GABA+GABA']
    ax.scatter(glut_gaba['x'], glut_gaba['y'], c='#FF6B6B', s=120,
               marker='*', edgecolors='white', linewidths=0.8, zorder=10,
               label=f'Glut+GABA doublet (n={len(glut_gaba)})')
    ax.scatter(gaba_gaba['x'], gaba_gaba['y'], c='#FF6BFF', s=120,
               marker='*', edgecolors='white', linewidths=0.8, zorder=10,
               label=f'GABA+GABA doublet (n={len(gaba_gaba)})')

    ax.set_title('Cell Class + Doublets', fontsize=16, fontweight='bold', color='white')
    ax.legend(fontsize=10, loc='lower left', facecolor='#333333',
              edgecolor='#555555', labelcolor='white', markerscale=0.8)

    # Panel B: GABA marker score
    marker_genes = [g for g in GABA_MARKERS if g in region.columns]
    gaba_score = (region[marker_genes] > 0).sum(axis=1)

    ax = axes[1]
    sc = ax.scatter(region['x'], region['y'], c=gaba_score.values,
                    cmap='YlOrRd', s=25, alpha=0.7, edgecolors='none',
                    vmin=0, vmax=7)
    # Highlight doublets
    ax.scatter(glut_gaba['x'], glut_gaba['y'], facecolors='none',
               edgecolors='white', s=120, linewidths=1.5, zorder=10)
    ax.scatter(gaba_gaba['x'], gaba_gaba['y'], facecolors='none',
               edgecolors='cyan', s=120, linewidths=1.5, zorder=10)

    cb = plt.colorbar(sc, ax=ax, shrink=0.7)
    cb.set_label('# GABA markers detected', fontsize=12, color='white')
    cb.ax.yaxis.set_tick_params(color='white')
    plt.setp(cb.ax.yaxis.get_ticklabels(), color='white')

    ax.set_title('GABAergic Marker Score', fontsize=16, fontweight='bold', color='white')

    # Panel C: Glut marker score
    marker_genes_g = [g for g in GLUT_MARKERS if g in region.columns]
    glut_score = (region[marker_genes_g] > 0).sum(axis=1)

    ax = axes[2]
    sc2 = ax.scatter(region['x'], region['y'], c=glut_score.values,
                     cmap='YlGnBu', s=25, alpha=0.7, edgecolors='none',
                     vmin=0, vmax=4)
    ax.scatter(glut_gaba['x'], glut_gaba['y'], facecolors='none',
               edgecolors='white', s=120, linewidths=1.5, zorder=10)
    ax.scatter(gaba_gaba['x'], gaba_gaba['y'], facecolors='none',
               edgecolors='cyan', s=120, linewidths=1.5, zorder=10)

    cb2 = plt.colorbar(sc2, ax=ax, shrink=0.7)
    cb2.set_label('# Glut markers detected', fontsize=12, color='white')
    cb2.ax.yaxis.set_tick_params(color='white')
    plt.setp(cb2.ax.yaxis.get_ticklabels(), color='white')

    ax.set_title('Glutamatergic Marker Score', fontsize=16, fontweight='bold', color='white')

    # Style all axes
    for ax in axes:
        ax.set_facecolor(BG)
        ax.set_aspect('equal')
        ax.tick_params(colors='white', labelsize=10)
        for spine in ax.spines.values():
            spine.set_color('#555555')
        ax.set_xlabel('x (µm)', fontsize=12, color='white')
        ax.set_ylabel('y (µm)', fontsize=12, color='white')
        # Invert y for spatial convention
        ax.invert_yaxis()

    fig.suptitle(f'Spatial Context of Doublet Cells ({SAMPLE})',
                 fontsize=22, fontweight='bold', color='white', y=1.02)

    plt.tight_layout()
    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


def main():
    print(f"=== Doublet Visualization ({SAMPLE}) ===\n")

    marker_expr, adata = load_sample_data(SAMPLE)

    n_doublets = marker_expr['doublet_suspect'].sum()
    n_gg = (marker_expr['doublet_type'] == 'Glut+GABA').sum()
    n_gabagaba = (marker_expr['doublet_type'] == 'GABA+GABA').sum()
    print(f"\n  Total doublets: {n_doublets} ({n_gg} Glut+GABA, {n_gabagaba} GABA+GABA)")
    print()

    # Figure 1: Marker expression heatmap
    print("Figure 1: Marker expression heatmap...")
    plot_marker_heatmap(
        marker_expr,
        os.path.join(OUT_DIR, 'doublet_marker_heatmap.png')
    )

    # Figure 2: Marker detection rate barplot
    print("\nFigure 2: Marker detection rate barplot...")
    plot_marker_barplot(
        marker_expr,
        os.path.join(OUT_DIR, 'doublet_marker_barplot.png')
    )

    # Figure 3: UMI count comparison
    print("\nFigure 3: UMI count comparison...")
    plot_umi_comparison(
        marker_expr,
        os.path.join(OUT_DIR, 'doublet_umi_comparison.png')
    )

    # Figure 4: Spatial zoom
    print("\nFigure 4: Spatial zoom...")
    plot_spatial_doublet_zoom(
        marker_expr,
        os.path.join(OUT_DIR, 'doublet_spatial_zoom.png')
    )

    print("\n=== Done ===")


if __name__ == '__main__':
    main()
