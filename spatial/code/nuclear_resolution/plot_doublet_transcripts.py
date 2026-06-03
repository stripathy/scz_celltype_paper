#!/usr/bin/env python3
"""
Visualize individual transcript molecules within doublet vs normal cells.

Shows cell + nucleus boundary polygons with per-molecule dots colored by
marker type (GABAergic vs Glutamatergic), demonstrating that doublet cells
contain transcripts from both neuronal classes within a single segmented cell.

Output: output/presentation/doublet_transcript_examples.png
"""

import os
import sys
import json
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BG_COLOR, H5AD_DIR, PRESENTATION_DIR

# ── Config ──
SAMPLE = 'Br8667'
BG = BG_COLOR
OUT_DIR = PRESENTATION_DIR

GABA_VIZ_GENES = ['GAD1', 'GAD2', 'SLC32A1', 'SST', 'PVALB', 'VIP', 'LAMP5']
GLUT_VIZ_GENES = ['CUX2', 'RORB', 'GRIN2A', 'THEMIS']

TRANSCRIPT_DIR = f'/Users/shreejoy/Github/SCZ_Xenium/output/viewer/transcripts/{SAMPLE}'
BOUNDARY_PATH = f'/Users/shreejoy/Github/SCZ_Xenium/data/raw/GSM9223485_{SAMPLE}-cell_boundaries.csv.gz'
NUCLEUS_PATH = f'/Users/shreejoy/Github/SCZ_Xenium/data/raw/GSM9223485_{SAMPLE}-nucleus_boundaries.csv.gz'
H5AD_PATH = os.path.join(H5AD_DIR, f'{SAMPLE}_annotated.h5ad')

GABA_COLOR = '#FF4444'
GLUT_COLOR = '#44AAFF'
OTHER_COLOR = '#333333'
PAD = 18  # microns around cell

# Per-gene colors for GABA+GABA legend
GABA_GENE_COLORS = {
    'GAD1': '#FF4444', 'GAD2': '#FF7744', 'SLC32A1': '#FF9944',
    'SST': '#FF44FF', 'PVALB': '#FF0000', 'VIP': '#FFAA00', 'LAMP5': '#FF6688',
}
GLUT_GENE_COLORS = {
    'CUX2': '#44AAFF', 'RORB': '#2288DD', 'GRIN2A': '#66CCFF', 'THEMIS': '#0066CC',
}


def load_transcript_index():
    with open(os.path.join(TRANSCRIPT_DIR, 'gene_index.json')) as f:
        return json.load(f)


def load_gene_transcripts(gene, index_info):
    gene_file = os.path.join(TRANSCRIPT_DIR, f'{gene}.json')
    if not os.path.exists(gene_file):
        return np.array([]), np.array([])
    with open(gene_file) as f:
        data = json.load(f)
    x = np.array(data['x'], dtype=np.float64) * index_info['x_scale'] + index_info['x_offset']
    y = np.array(data['y'], dtype=np.float64) * index_info['y_scale'] + index_info['y_offset']
    return x, y


def load_boundaries(path, label="boundaries"):
    print(f"Loading {label}...")
    df = pd.read_csv(path, compression='gzip')
    print(f"  {len(df):,} vertices, {df['cell_id'].nunique():,} cells")
    return df


def find_example_cells(adata):
    obs = adata.obs
    X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    gene_names = list(adata.var_names)

    def _score(genes):
        s = np.zeros(X.shape[0])
        for g in genes:
            if g in gene_names:
                s += (X[:, gene_names.index(g)] > 0).astype(float)
        return s

    gaba_score = _score(GABA_VIZ_GENES)
    glut_score = _score(GLUT_VIZ_GENES)

    # Glut+GABA doublets — highest combined score
    mask = (obs['doublet_type'].astype(str) == 'Glut+GABA') & (obs['qc_pass'] == True)
    idx = np.where(mask.values)[0]
    best_doublets = idx[np.argsort(-(gaba_score[idx] + glut_score[idx]))][:20]

    # Pure normal Glut (0 GABA markers, 3+ Glut)
    mask = ((obs['corr_class'].astype(str) == 'Glutamatergic') &
            (~obs['doublet_suspect']) & (obs['qc_pass'] == True))
    idx = np.where(mask.values)[0]
    pure_glut = idx[(gaba_score[idx] == 0) & (glut_score[idx] >= 3)]

    # Pure normal GABA (0 Glut markers, 4+ GABA)
    mask = ((obs['corr_class'].astype(str) == 'GABAergic') &
            (~obs['doublet_suspect']) & (obs['qc_pass'] == True))
    idx = np.where(mask.values)[0]
    pure_gaba = idx[(glut_score[idx] == 0) & (gaba_score[idx] >= 4)]

    # GABA+GABA doublets
    mask = (obs['doublet_type'].astype(str) == 'GABA+GABA') & (obs['qc_pass'] == True)
    gabagaba = np.where(mask.values)[0]

    return best_doublets, pure_glut, pure_gaba, gabagaba


def plot_cell_panel(ax, cell_idx, adata, cell_bnd, nuc_bnd,
                    gene_transcripts, title, title_color='white'):
    """Plot a single cell with cell + nucleus boundaries and transcripts."""
    obs = adata.obs.iloc[cell_idx]
    cell_id = obs.name

    # Get boundaries by direct cell_id match
    cell_verts = cell_bnd[cell_bnd['cell_id'] == cell_id][['vertex_x', 'vertex_y']].values
    nuc_verts = nuc_bnd[nuc_bnd['cell_id'] == cell_id][['vertex_x', 'vertex_y']].values

    # Fallback to centroid proximity if cell_id doesn't match
    cx = obs.get('x_centroid', None)
    cy = obs.get('y_centroid', None)
    if cx is None and 'spatial' in adata.obsm:
        cx, cy = adata.obsm['spatial'][cell_idx]

    if len(cell_verts) == 0 and cx is not None:
        nearby = cell_bnd[
            (cell_bnd['vertex_x'].between(cx - 50, cx + 50)) &
            (cell_bnd['vertex_y'].between(cy - 50, cy + 50))
        ]
        best_id, best_d = None, 1e9
        for cid in nearby['cell_id'].unique():
            v = nearby[nearby['cell_id'] == cid][['vertex_x', 'vertex_y']].values
            d = np.hypot(*(v.mean(0) - [cx, cy]))
            if d < best_d:
                best_d, best_id = d, cid
        if best_id is not None and best_d < 20:
            cell_verts = cell_bnd[cell_bnd['cell_id'] == best_id][['vertex_x', 'vertex_y']].values
            nuc_verts = nuc_bnd[nuc_bnd['cell_id'] == best_id][['vertex_x', 'vertex_y']].values

    # Determine view region
    if len(cell_verts) > 0:
        xmin, xmax = cell_verts[:, 0].min() - PAD, cell_verts[:, 0].max() + PAD
        ymin, ymax = cell_verts[:, 1].min() - PAD, cell_verts[:, 1].max() + PAD
    else:
        xmin, xmax = cx - 30, cx + 30
        ymin, ymax = cy - 30, cy + 30

    # ── Transcripts ──
    # Other genes (very dim background)
    for gene, (gx, gy) in gene_transcripts.items():
        if gene in GABA_VIZ_GENES + GLUT_VIZ_GENES:
            continue
        m = (gx >= xmin) & (gx <= xmax) & (gy >= ymin) & (gy <= ymax)
        if m.sum() > 0:
            ax.scatter(gx[m], gy[m], c=OTHER_COLOR, s=2, alpha=0.12, zorder=1)

    # GABA marker transcripts (per-gene colors)
    for gene in GABA_VIZ_GENES:
        if gene in gene_transcripts:
            gx, gy = gene_transcripts[gene]
            m = (gx >= xmin) & (gx <= xmax) & (gy >= ymin) & (gy <= ymax)
            if m.sum() > 0:
                c = GABA_GENE_COLORS.get(gene, GABA_COLOR)
                ax.scatter(gx[m], gy[m], c=c, s=18, alpha=0.85,
                           edgecolors='white', linewidths=0.2, zorder=5)

    # Glut marker transcripts (per-gene colors)
    for gene in GLUT_VIZ_GENES:
        if gene in gene_transcripts:
            gx, gy = gene_transcripts[gene]
            m = (gx >= xmin) & (gx <= xmax) & (gy >= ymin) & (gy <= ymax)
            if m.sum() > 0:
                c = GLUT_GENE_COLORS.get(gene, GLUT_COLOR)
                ax.scatter(gx[m], gy[m], c=c, s=18, alpha=0.85,
                           edgecolors='white', linewidths=0.2, zorder=5)

    # ── Neighboring cell boundaries ──
    nearby_cells = cell_bnd[
        (cell_bnd['vertex_x'].between(xmin, xmax)) &
        (cell_bnd['vertex_y'].between(ymin, ymax))
    ]
    for cid in nearby_cells['cell_id'].unique():
        if cid == cell_id:
            continue
        v = cell_bnd[cell_bnd['cell_id'] == cid][['vertex_x', 'vertex_y']].values
        if len(v) > 2:
            mx, my = v.mean(0)
            if xmin < mx < xmax and ymin < my < ymax:
                ax.add_patch(Polygon(v, closed=True, fill=False,
                                     edgecolor='#777777', linewidth=1.0, alpha=0.55, zorder=2))
    # Neighboring nucleus boundaries
    nearby_nuc = nuc_bnd[
        (nuc_bnd['vertex_x'].between(xmin, xmax)) &
        (nuc_bnd['vertex_y'].between(ymin, ymax))
    ]
    for cid in nearby_nuc['cell_id'].unique():
        if cid == cell_id:
            continue
        v = nuc_bnd[nuc_bnd['cell_id'] == cid][['vertex_x', 'vertex_y']].values
        if len(v) > 2:
            mx, my = v.mean(0)
            if xmin < mx < xmax and ymin < my < ymax:
                ax.add_patch(Polygon(v, closed=True, fill=False,
                                     edgecolor='#666666', linewidth=0.7,
                                     linestyle='--', alpha=0.4, zorder=2))

    # ── Focus cell boundary (solid white) ──
    if len(cell_verts) > 0:
        ax.add_patch(Polygon(cell_verts, closed=True, fill=False,
                             edgecolor='white', linewidth=2.5, alpha=0.95, zorder=8))

    # ── Focus nucleus boundary (dashed yellow) ──
    if len(nuc_verts) > 0:
        ax.add_patch(Polygon(nuc_verts, closed=True, fill=False,
                             edgecolor='#FFD700', linewidth=1.8,
                             linestyle='--', alpha=0.9, zorder=9))

    # ── Per-gene legend (show which markers are present in this cell) ──
    X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    gn = list(adata.var_names)
    legend_items = []
    for g in GABA_VIZ_GENES:
        if g in gn and X[cell_idx, gn.index(g)] > 0:
            c = GABA_GENE_COLORS.get(g, GABA_COLOR)
            legend_items.append(
                Line2D([0], [0], marker='o', color='w', markerfacecolor=c,
                       markersize=7, label=g, linestyle='None'))
    for g in GLUT_VIZ_GENES:
        if g in gn and X[cell_idx, gn.index(g)] > 0:
            c = GLUT_GENE_COLORS.get(g, GLUT_COLOR)
            legend_items.append(
                Line2D([0], [0], marker='o', color='w', markerfacecolor=c,
                       markersize=7, label=g, linestyle='None'))
    if legend_items:
        ax.legend(handles=legend_items, loc='upper right', fontsize=7,
                  facecolor='#1a1a2e', edgecolor='#444444', labelcolor='white',
                  framealpha=0.85, handletextpad=0.3, borderpad=0.3)

    # ── Title with cell ID and coordinates ──
    subclass = str(obs.get('corr_subclass', obs.get('subclass_label', '?')))
    ax.set_title(title, fontsize=13, fontweight='bold', color=title_color, pad=6)

    # Info line: cell_id, section, coordinates, UMI
    x_coord = cx if cx is not None else (cell_verts[:, 0].mean() if len(cell_verts) > 0 else 0)
    y_coord = cy if cy is not None else (cell_verts[:, 1].mean() if len(cell_verts) > 0 else 0)
    info = (f"{SAMPLE} | {cell_id}\n"
            f"{subclass} | ({x_coord:.0f}, {y_coord:.0f}) µm | UMI: {int(obs['total_counts'])}")
    ax.text(0.5, -0.01, info, transform=ax.transAxes, ha='center', va='top',
            fontsize=7.5, color='#bbbbbb', linespacing=1.4)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymax, ymin)
    ax.set_aspect('equal')
    ax.set_facecolor('#0d0d1a')
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    for spine in ax.spines.values():
        spine.set_color('#333333')

    # Scale bar (bottom-left)
    bar_len = 10
    bx = xmin + 2
    by = ymax - 2
    ax.plot([bx, bx + bar_len], [by, by], color='white', linewidth=2, zorder=10)
    ax.text(bx + bar_len / 2, by - 1, '10 µm', color='white',
            fontsize=7, ha='center', va='top', zorder=10)


def main():
    print(f"=== Doublet Transcript Visualization ({SAMPLE}) ===\n")

    adata = ad.read_h5ad(H5AD_PATH)
    print(f"  {adata.n_obs:,} cells loaded")

    best_doublets, pure_glut, pure_gaba, gabagaba = find_example_cells(adata)
    print(f"  Doublet candidates: {len(best_doublets)} Glut+GABA, {len(gabagaba)} GABA+GABA")
    print(f"  Pure Glut: {len(pure_glut)}, Pure GABA: {len(pure_gaba)}")

    index_info = load_transcript_index()
    gene_transcripts = {}
    for gene in GABA_VIZ_GENES + GLUT_VIZ_GENES:
        x, y = load_gene_transcripts(gene, index_info)
        if len(x) > 0:
            gene_transcripts[gene] = (x, y)
            print(f"  {gene}: {len(x):,} transcripts")

    cell_bnd = load_boundaries(BOUNDARY_PATH, "cell boundaries")
    nuc_bnd = load_boundaries(NUCLEUS_PATH, "nucleus boundaries")

    # ── 2 x 3 figure ──
    np.random.seed(42)
    fig, axes = plt.subplots(2, 3, figsize=(18, 13), facecolor=BG)

    cells = [
        (pure_glut[np.random.randint(0, min(100, len(pure_glut)))],
         'Normal Glutamatergic', '#44AAFF'),
        (best_doublets[0], 'Glut+GABA Doublet', '#FF6B6B'),
        (best_doublets[1], 'Glut+GABA Doublet', '#FF6B6B'),
        (pure_gaba[np.random.randint(0, min(100, len(pure_gaba)))],
         'Normal GABAergic', '#F05A28'),
        (gabagaba[0] if len(gabagaba) > 0 else None,
         'GABA+GABA Doublet', '#FF6BFF'),
        (best_doublets[2], 'Glut+GABA Doublet', '#FF6B6B'),
    ]

    for i, (idx, title, color) in enumerate(cells):
        r, c = divmod(i, 3)
        if idx is None:
            axes[r, c].set_visible(False)
            continue
        print(f"  Panel {i+1}: {title} (cell {adata.obs_names[idx]})...")
        plot_cell_panel(axes[r, c], idx, adata, cell_bnd, nuc_bnd,
                        gene_transcripts, title, color)

    # ── Legend: per-gene colors + boundary types ──
    legend_elements = []
    # GABA genes
    for g in GABA_VIZ_GENES:
        legend_elements.append(
            Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=GABA_GENE_COLORS[g],
                   markersize=8, label=g, linestyle='None'))
    # Glut genes
    for g in GLUT_VIZ_GENES:
        legend_elements.append(
            Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=GLUT_GENE_COLORS[g],
                   markersize=8, label=g, linestyle='None'))
    # Boundary types
    legend_elements.extend([
        Line2D([0], [0], color='white', linewidth=2, label='Cell boundary'),
        Line2D([0], [0], color='#FFD700', linewidth=1.5, linestyle='--',
               label='Nucleus boundary'),
        Line2D([0], [0], color='#777777', linewidth=0.8, label='Neighbor cells'),
    ])
    fig.legend(handles=legend_elements, loc='lower center',
               ncol=7, fontsize=9.5, facecolor='#1a1a2e', edgecolor='#555555',
               labelcolor='white', bbox_to_anchor=(0.5, 0.005),
               handletextpad=0.3, columnspacing=1.0)

    fig.suptitle(
        'Transcript Molecules Within Individual Cells: Doublets vs Normal',
        fontsize=19, fontweight='bold', color='white', y=0.97
    )

    plt.tight_layout(rect=[0, 0.04, 1, 0.95])
    outpath = os.path.join(OUT_DIR, 'doublet_transcript_examples.png')
    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == '__main__':
    main()
