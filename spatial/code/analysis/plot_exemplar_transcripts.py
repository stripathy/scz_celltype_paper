#!/usr/bin/env python3
"""
Visualize transcript molecules in exemplar well-classified cells.

Shows cells from diverse subclasses with their characteristic marker
transcripts, plus cell and nucleus boundaries, demonstrating that the
300-gene Xenium panel provides clear marker-based cell type identity
at the single-molecule level.

Output: output/presentation/exemplar_transcript_classification.png
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

# Exemplar subclasses with key markers
EXEMPLARS = [
    {
        'subclass': 'L2/3 IT', 'class': 'Glutamatergic',
        'markers': ['CUX2', 'LAMP5'],
        'negative_markers': ['SST', 'PVALB', 'GAD1', 'GAD2'],
        'color': '#B1EC30', 'title': 'L2/3 IT',
    },
    {
        'subclass': 'L4 IT', 'class': 'Glutamatergic',
        'markers': ['RORB', 'GRIN2A'],
        'negative_markers': ['SST', 'PVALB', 'GAD1', 'GAD2'],
        'color': '#6DE843', 'title': 'L4 IT',
    },
    {
        'subclass': 'L6b', 'class': 'Glutamatergic',
        'markers': ['THEMIS', 'GRIN2A'],
        'negative_markers': ['CUX2', 'GAD1', 'SST', 'PVALB'],
        'color': '#007B10', 'title': 'L6b',
    },
    {
        'subclass': 'Pvalb', 'class': 'GABAergic',
        'markers': ['PVALB', 'GAD1', 'GAD2'],
        'negative_markers': ['SST', 'VIP', 'CUX2', 'RORB'],
        'color': '#D93137', 'title': 'Pvalb',
    },
    {
        'subclass': 'Sst', 'class': 'GABAergic',
        'markers': ['SST', 'GAD1', 'GAD2'],
        'negative_markers': ['PVALB', 'VIP', 'CUX2', 'RORB'],
        'color': '#FF8B94', 'title': 'Sst',
    },
    {
        'subclass': 'Vip', 'class': 'GABAergic',
        'markers': ['VIP', 'GAD1', 'GAD2'],
        'negative_markers': ['SST', 'PVALB', 'CUX2', 'RORB'],
        'color': '#FF6B6B', 'title': 'Vip',
    },
]

ALL_VIZ_GENES = list(set(
    g for ex in EXEMPLARS for g in ex['markers'] + ex['negative_markers']
))

# Per-gene colors for transcript dots
GENE_COLORS = {
    'CUX2': '#44AAFF', 'RORB': '#2288DD', 'GRIN2A': '#66CCFF',
    'THEMIS': '#0066CC', 'GAD1': '#FF4444', 'GAD2': '#FF7744',
    'SLC32A1': '#FF9944', 'SST': '#FF44FF', 'PVALB': '#FF0000',
    'VIP': '#FFAA00', 'LAMP5': '#FF6688',
}

TRANSCRIPT_DIR = f'/Users/shreejoy/Github/SCZ_Xenium/output/viewer/transcripts/{SAMPLE}'
BOUNDARY_PATH = f'/Users/shreejoy/Github/SCZ_Xenium/data/raw/GSM9223485_{SAMPLE}-cell_boundaries.csv.gz'
NUCLEUS_PATH = f'/Users/shreejoy/Github/SCZ_Xenium/data/raw/GSM9223485_{SAMPLE}-nucleus_boundaries.csv.gz'
H5AD_PATH = os.path.join(H5AD_DIR, f'{SAMPLE}_annotated.h5ad')

PAD = 18


def load_transcript_index():
    with open(os.path.join(TRANSCRIPT_DIR, 'gene_index.json')) as f:
        return json.load(f)


def load_gene_transcripts(gene, index_info):
    path = os.path.join(TRANSCRIPT_DIR, f'{gene}.json')
    if not os.path.exists(path):
        return np.array([]), np.array([])
    with open(path) as f:
        d = json.load(f)
    x = np.array(d['x'], dtype=np.float64) * index_info['x_scale'] + index_info['x_offset']
    y = np.array(d['y'], dtype=np.float64) * index_info['y_scale'] + index_info['y_offset']
    return x, y


def load_boundaries(path, label="boundaries"):
    print(f"Loading {label}...")
    df = pd.read_csv(path, compression='gzip')
    print(f"  {len(df):,} vertices, {df['cell_id'].nunique():,} cells")
    return df


def find_exemplar_cell(adata, X_dense, gene_names, ex):
    obs = adata.obs
    mask = (
        (obs['corr_subclass'].astype(str) == ex['subclass']) &
        (~obs['doublet_suspect']) & (obs['qc_pass'] == True)
    )
    candidates = np.where(mask.values)[0]
    if len(candidates) == 0:
        return None

    pos = np.zeros(len(candidates))
    for g in ex['markers']:
        if g in gene_names:
            gi = gene_names.index(g)
            pos += (X_dense[candidates, gi] > 0).astype(float)
            pos += np.log1p(X_dense[candidates, gi]) * 0.5

    neg = np.zeros(len(candidates))
    for g in ex['negative_markers']:
        if g in gene_names:
            gi = gene_names.index(g)
            neg += (X_dense[candidates, gi] > 0).astype(float)

    score = pos - neg * 5
    top = candidates[score >= np.percentile(score, 95)]
    if len(top) == 0:
        top = candidates[np.argsort(-score)[:10]]

    # Among top, pick decent total counts
    counts = adata.obs.iloc[top]['total_counts'].values
    return top[np.argmax(counts)]


def plot_exemplar_panel(ax, cell_idx, adata, cell_bnd, nuc_bnd,
                        gene_transcripts, ex):
    obs = adata.obs.iloc[cell_idx]
    cell_id = obs.name

    cell_verts = cell_bnd[cell_bnd['cell_id'] == cell_id][['vertex_x', 'vertex_y']].values
    nuc_verts = nuc_bnd[nuc_bnd['cell_id'] == cell_id][['vertex_x', 'vertex_y']].values

    cx = obs.get('x_centroid', None)
    cy = obs.get('y_centroid', None)

    # Fallback centroid proximity match
    if len(cell_verts) == 0 and cx is not None:
        nearby = cell_bnd[
            cell_bnd['vertex_x'].between(cx - 50, cx + 50) &
            cell_bnd['vertex_y'].between(cy - 50, cy + 50)
        ]
        best_id, best_d = None, 1e9
        for cid in nearby['cell_id'].unique():
            v = nearby[nearby['cell_id'] == cid][['vertex_x', 'vertex_y']].values
            d = np.hypot(*(v.mean(0) - [cx, cy]))
            if d < best_d:
                best_d, best_id = d, cid
        if best_id and best_d < 20:
            cell_verts = cell_bnd[cell_bnd['cell_id'] == best_id][['vertex_x', 'vertex_y']].values
            nuc_verts = nuc_bnd[nuc_bnd['cell_id'] == best_id][['vertex_x', 'vertex_y']].values

    if len(cell_verts) > 0:
        xmin, xmax = cell_verts[:, 0].min() - PAD, cell_verts[:, 0].max() + PAD
        ymin, ymax = cell_verts[:, 1].min() - PAD, cell_verts[:, 1].max() + PAD
    else:
        xmin, xmax = cx - 30, cx + 30
        ymin, ymax = cy - 30, cy + 30

    # ── Background transcripts (dim) ──
    for gene, (gx, gy) in gene_transcripts.items():
        if gene in ex['markers'] + ex['negative_markers']:
            continue
        m = (gx >= xmin) & (gx <= xmax) & (gy >= ymin) & (gy <= ymax)
        if m.sum() > 0:
            ax.scatter(gx[m], gy[m], c='#2a2a2a', s=2, alpha=0.1, zorder=1)

    # ── Negative markers (dim, small) ──
    for gene in ex['negative_markers']:
        if gene in gene_transcripts:
            gx, gy = gene_transcripts[gene]
            m = (gx >= xmin) & (gx <= xmax) & (gy >= ymin) & (gy <= ymax)
            if m.sum() > 0:
                c = GENE_COLORS.get(gene, '#666666')
                ax.scatter(gx[m], gy[m], c=c, s=8, alpha=0.25, zorder=3)

    # ── Positive markers (bright, labeled in legend) ──
    legend_items = []
    for gene in ex['markers']:
        if gene in gene_transcripts:
            gx, gy = gene_transcripts[gene]
            m = (gx >= xmin) & (gx <= xmax) & (gy >= ymin) & (gy <= ymax)
            if m.sum() > 0:
                c = GENE_COLORS.get(gene, '#FFFFFF')
                ax.scatter(gx[m], gy[m], c=c, s=22, alpha=0.9,
                           edgecolors='white', linewidths=0.2, zorder=6)
                legend_items.append(
                    Line2D([0], [0], marker='o', color='w', markerfacecolor=c,
                           markersize=7, label=gene, linestyle='None'))

    # ── Neighboring boundaries ──
    nearby_c = cell_bnd[
        cell_bnd['vertex_x'].between(xmin, xmax) &
        cell_bnd['vertex_y'].between(ymin, ymax)
    ]
    for cid in nearby_c['cell_id'].unique():
        if cid == cell_id:
            continue
        v = cell_bnd[cell_bnd['cell_id'] == cid][['vertex_x', 'vertex_y']].values
        if len(v) > 2:
            mx, my = v.mean(0)
            if xmin < mx < xmax and ymin < my < ymax:
                ax.add_patch(Polygon(v, closed=True, fill=False,
                                     edgecolor='#777777', lw=1.0, alpha=0.55, zorder=2))

    nearby_n = nuc_bnd[
        nuc_bnd['vertex_x'].between(xmin, xmax) &
        nuc_bnd['vertex_y'].between(ymin, ymax)
    ]
    for cid in nearby_n['cell_id'].unique():
        if cid == cell_id:
            continue
        v = nuc_bnd[nuc_bnd['cell_id'] == cid][['vertex_x', 'vertex_y']].values
        if len(v) > 2:
            mx, my = v.mean(0)
            if xmin < mx < xmax and ymin < my < ymax:
                ax.add_patch(Polygon(v, closed=True, fill=False,
                                     edgecolor='#666666', lw=0.7,
                                     linestyle='--', alpha=0.4, zorder=2))

    # ── Focus cell boundary (solid, subclass color) ──
    if len(cell_verts) > 0:
        ax.add_patch(Polygon(cell_verts, closed=True, fill=False,
                             edgecolor=ex['color'], linewidth=2.2, alpha=0.95, zorder=8))

    # ── Focus nucleus boundary (dashed yellow) ──
    if len(nuc_verts) > 0:
        ax.add_patch(Polygon(nuc_verts, closed=True, fill=False,
                             edgecolor='#FFD700', linewidth=1.8,
                             linestyle='--', alpha=0.9, zorder=9))

    # ── Title + info with cell ID and coordinates ──
    ax.set_title(ex['title'], fontsize=14, fontweight='bold',
                 color=ex['color'], pad=6)
    x_coord = cx if cx is not None else (cell_verts[:, 0].mean() if len(cell_verts) > 0 else 0)
    y_coord = cy if cy is not None else (cell_verts[:, 1].mean() if len(cell_verts) > 0 else 0)
    info = (f"{SAMPLE} | {cell_id}\n"
            f"({x_coord:.0f}, {y_coord:.0f}) µm | UMI: {int(obs['total_counts'])}")
    ax.text(0.5, -0.01, info, transform=ax.transAxes, ha='center', va='top',
            fontsize=7.5, color='#bbbbbb', linespacing=1.4)

    if legend_items:
        ax.legend(handles=legend_items, loc='upper right', fontsize=8,
                  facecolor='#1a1a2e', edgecolor='#444444', labelcolor='white',
                  framealpha=0.85, handletextpad=0.3, borderpad=0.3)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymax, ymin)
    ax.set_aspect('equal')
    ax.set_facecolor('#0d0d1a')
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    for spine in ax.spines.values():
        spine.set_color('#333333')

    # Scale bar
    bx, by = xmin + 2, ymax - 2
    ax.plot([bx, bx + 10], [by, by], color='white', linewidth=2, zorder=10)
    ax.text(bx + 5, by - 1, '10 µm', color='white', fontsize=7,
            ha='center', va='top', zorder=10)


def main():
    print(f"=== Exemplar Cell Transcript Visualization ({SAMPLE}) ===\n")

    adata = ad.read_h5ad(H5AD_PATH)
    X_dense = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    gene_names = list(adata.var_names)

    index_info = load_transcript_index()
    gene_transcripts = {}
    for gene in ALL_VIZ_GENES:
        x, y = load_gene_transcripts(gene, index_info)
        if len(x) > 0:
            gene_transcripts[gene] = (x, y)
            print(f"  {gene}: {len(x):,} transcripts")

    cell_bnd = load_boundaries(BOUNDARY_PATH, "cell boundaries")
    nuc_bnd = load_boundaries(NUCLEUS_PATH, "nucleus boundaries")

    print("\nFinding exemplar cells...")
    cell_indices = []
    for ex in EXEMPLARS:
        idx = find_exemplar_cell(adata, X_dense, gene_names, ex)
        cell_indices.append(idx)
        if idx is not None:
            print(f"  {ex['subclass']}: {adata.obs_names[idx]}, "
                  f"UMI={int(adata.obs.iloc[idx]['total_counts'])}")

    # ── 2 x 3 figure ──
    fig, axes = plt.subplots(2, 3, figsize=(18, 13), facecolor=BG)

    for i, (ex, idx) in enumerate(zip(EXEMPLARS, cell_indices)):
        r, c = divmod(i, 3)
        if idx is not None:
            print(f"  Plotting {ex['subclass']}...")
            plot_exemplar_panel(axes[r, c], idx, adata, cell_bnd, nuc_bnd,
                                gene_transcripts, ex)
        else:
            axes[r, c].set_visible(False)

    # ── Legend ──
    legend_elements = [
        Line2D([0], [0], color='gray', linewidth=2, label='Cell boundary'),
        Line2D([0], [0], color='#FFD700', linewidth=1.5, linestyle='--',
               label='Nucleus boundary'),
        Line2D([0], [0], color='#555555', linewidth=0.8, label='Neighbor cells'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3,
               fontsize=11, facecolor='#1a1a2e', edgecolor='#555555',
               labelcolor='white', bbox_to_anchor=(0.5, 0.005))

    fig.suptitle(
        'Transcript-Level Cell Type Identity: Exemplar Cells',
        fontsize=19, fontweight='bold', color='white', y=0.97
    )

    plt.tight_layout(rect=[0, 0.04, 1, 0.95])
    outpath = os.path.join(OUT_DIR, 'exemplar_transcript_classification.png')
    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == '__main__':
    main()
