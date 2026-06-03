#!/usr/bin/env python3
"""
Ghost Cells Prototype: Identify unsegmented cells from unassigned transcripts.

In Xenium brain tissue, ~40–55% of transcripts fall outside segmented cell
boundaries. Some arise from truly missing cells (DAPI too faint/out of focus),
others from cytoplasmic/neuritic RNA or diffusion artifacts. This prototype
quantifies the unassigned fraction, maps spatial hotspots, clusters unassigned
transcripts into candidate "ghost cells" via DBSCAN, and visualizes exemplars.

Parts:
  0. Data loading: transcripts (skip MT genes), cell boundaries → STRtree
  1. Assignment fraction: point-in-polygon → per-gene & overall stats → CSV
  2. Spatial density maps: 2×2 panel comparing transcript vs cell density
  3. DBSCAN clustering: find spatially coherent unassigned transcript clusters
  4. Exemplar visualization: 2×3 panel of candidate ghost cells with boundaries

Sample: Br5400 (medium-sized, clean cortical coverage)
Output: output/ghost_cells/

Usage:
    python3 -u ghost_cells_prototype.py
"""

import os
import sys
import json
import time
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
from shapely import Polygon, STRtree, points
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BG_COLOR, H5AD_DIR, PRESENTATION_DIR, style_dark_axis

# Import nucleus polygon loader from nuclear_resolution
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "nuclear_resolution"))
from nuclear_counts import load_nucleus_polygons

# ══════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
SAMPLE = 'Br5400'
CELL_BOUNDARY_PATH = os.path.join(
    BASE_DIR, 'data/raw/GSM9223474_Br5400-cell_boundaries.csv.gz')
NUCLEUS_BOUNDARY_PATH = os.path.join(
    BASE_DIR, 'data/raw/GSM9223474_Br5400-nucleus_boundaries.csv.gz')
TRANSCRIPT_DIR = os.path.join(
    BASE_DIR, f'output/viewer/transcripts/{SAMPLE}')
H5AD_PATH = os.path.join(H5AD_DIR, f'{SAMPLE}_annotated.h5ad')
OUT_DIR = os.path.join(BASE_DIR, 'output/ghost_cells')

BG = BG_COLOR

# Genes to skip (mitochondrial + mitochondrial rRNA pseudogenes)
SKIP_PREFIXES = ('MT-', 'MTRNR2L')

# Point-in-polygon query
CHUNK_SIZE = 500_000

# DBSCAN parameters
DBSCAN_EPS = 8.0           # microns
DBSCAN_MIN_SAMPLES = 30
DBSCAN_MAX_POINTS = 5_000_000  # subsample if more unassigned than this

# Density grid
BIN_SIZE = 50.0            # microns

# Ghost cluster filters
MAX_CLUSTER_RADIUS = 25.0  # microns — plausible cell-sized

# Marker genes for ghost cluster classification
MARKER_GENES = {
    'Glutamatergic': ['CUX2', 'RORB', 'THEMIS', 'SLC17A7', 'GRIN2A',
                      'SATB2', 'SLC17A6'],
    'GABAergic': ['GAD1', 'GAD2', 'SST', 'PVALB', 'VIP', 'LAMP5',
                  'SLC32A1', 'SNCG'],
    'Astrocyte': ['GFAP', 'AQP4', 'SLC1A2', 'SLC1A3', 'ALDH1L1'],
    'Oligodendrocyte': ['MBP', 'MOBP', 'MOG', 'PLP1', 'OPALIN'],
    'Microglia': ['CSF1R', 'CX3CR1', 'AIF1', 'P2RY12', 'TMEM119'],
    'OPC': ['PDGFRA', 'VCAN', 'CSPG4'],
    'Endothelial': ['CLDN5', 'FLT1', 'PECAM1'],
    'VLMC': ['COL1A1', 'COL1A2', 'PDGFRB'],
}

# Per-gene colors for transcript dots in exemplar panels
GENE_COLORS = {
    'CUX2': '#44AAFF', 'RORB': '#2288DD', 'GRIN2A': '#66CCFF',
    'THEMIS': '#0066CC', 'GAD1': '#FF4444', 'GAD2': '#FF7744',
    'SLC32A1': '#FF9944', 'SST': '#FF44FF', 'PVALB': '#FF0000',
    'VIP': '#FFAA00', 'LAMP5': '#FF6688', 'GFAP': '#33FF33',
    'AQP4': '#22DD22', 'MBP': '#AAAAFF', 'MOBP': '#8888EE',
    'PLP1': '#9999FF', 'MOG': '#7777DD', 'CSF1R': '#FFFF44',
    'PDGFRA': '#FF88FF', 'CLDN5': '#88FFFF', 'SLC17A7': '#55BBFF',
    'SATB2': '#3399FF', 'OPALIN': '#BBBBFF',
}
def get_gene_color(gene_name):
    """Deterministic gene color: manual overrides for markers, hash-based for rest."""
    if gene_name in GENE_COLORS:
        return GENE_COLORS[gene_name]
    # Deterministic hash → bright, saturated HSV color
    h = hash(gene_name) % 360
    s = 0.7 + (hash(gene_name + '_s') % 20) / 100.0  # 0.70–0.89
    v = 0.8 + (hash(gene_name + '_v') % 20) / 100.0  # 0.80–0.99
    rgb = matplotlib.colors.hsv_to_rgb([h / 360.0, s, v])
    return '#{:02x}{:02x}{:02x}'.format(
        int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))


# ══════════════════════════════════════════════════════════════════════
# PART 0: DATA LOADING
# ══════════════════════════════════════════════════════════════════════

def load_gene_index(transcript_dir):
    """Parse gene_index.json for coordinate transform + gene list."""
    with open(os.path.join(transcript_dir, 'gene_index.json')) as f:
        return json.load(f)


def load_cell_polygons(boundary_path):
    """Load cell boundary CSV → Shapely polygons + STRtree.

    Adapted from nuclear_counts.py:load_nucleus_polygons().

    Returns
    -------
    cell_ids : list of str/int
    polygons : list of Polygon
    tree : STRtree
    centroids : np.ndarray of shape (n_cells, 2)
    """
    print("Loading cell boundaries...", flush=True)
    t0 = time.time()
    df = pd.read_csv(boundary_path, compression='gzip')
    n_cells_raw = df['cell_id'].nunique()
    print(f"  {len(df):,} vertices, {n_cells_raw:,} cells", flush=True)

    grouped = df.groupby('cell_id', sort=False)
    cell_ids = []
    polygons = []
    centroids_list = []
    n_invalid = 0

    for cell_id, grp in grouped:
        verts = grp[['vertex_x', 'vertex_y']].values
        if len(verts) < 3:
            n_invalid += 1
            continue
        try:
            poly = Polygon(verts)
            if not poly.is_valid or poly.area <= 0:
                poly = poly.buffer(0)
            if poly.is_valid and poly.area > 0:
                cell_ids.append(cell_id)
                polygons.append(poly)
                centroids_list.append(verts.mean(axis=0))
            else:
                n_invalid += 1
        except Exception:
            n_invalid += 1

    centroids = np.array(centroids_list)

    print(f"  Building STRtree ({len(polygons):,} polygons)...", flush=True)
    tree = STRtree(polygons)

    elapsed = time.time() - t0
    print(f"  {len(polygons):,} valid polygons ({n_invalid} invalid) "
          f"in {elapsed:.1f}s", flush=True)
    return cell_ids, polygons, tree, centroids


def load_all_transcripts(transcript_dir, gene_index):
    """Load all non-MT gene transcript coords from per-gene JSONs.

    Returns
    -------
    x : np.ndarray (float32)
    y : np.ndarray (float32)
    gene_indices : np.ndarray (uint16) — index into gene_names
    gene_names : list of str
    """
    x_scale = gene_index['x_scale']
    y_scale = gene_index['y_scale']
    x_offset = gene_index['x_offset']
    y_offset = gene_index['y_offset']

    all_x, all_y, all_gi = [], [], []
    gene_names = []
    total_skipped = 0
    total_loaded = 0

    genes = gene_index['genes']
    n_genes = len(genes)
    t0 = time.time()

    for i, ginfo in enumerate(genes):
        gene = ginfo['gene']
        n_total = ginfo['n']

        # Skip mitochondrial genes
        if any(gene.startswith(prefix) for prefix in SKIP_PREFIXES):
            total_skipped += n_total
            continue

        gene_idx = len(gene_names)
        gene_names.append(gene)

        # Load transcript coordinates
        gene_path = os.path.join(transcript_dir, ginfo['file'])
        with open(gene_path) as f:
            gdata = json.load(f)

        gx = np.array(gdata['x'], dtype=np.float32) * x_scale + x_offset
        gy = np.array(gdata['y'], dtype=np.float32) * y_scale + y_offset
        all_x.append(gx)
        all_y.append(gy)
        all_gi.append(np.full(len(gx), gene_idx, dtype=np.uint16))
        total_loaded += len(gx)

        if (i + 1) % 50 == 0 or i == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            print(f"  [{i+1:3d}/{n_genes}] {gene:12s}: {len(gx):>10,} tx "
                  f"[{rate:.0f} genes/sec]", flush=True)

    x = np.concatenate(all_x)
    y = np.concatenate(all_y)
    gene_indices = np.concatenate(all_gi)

    elapsed = time.time() - t0
    print(f"\n  Loaded {total_loaded:,} transcripts from {len(gene_names)} genes "
          f"({elapsed:.1f}s)")
    print(f"  Skipped {total_skipped:,} mitochondrial transcripts "
          f"({100*total_skipped/(total_loaded+total_skipped):.1f}%)")
    return x, y, gene_indices, gene_names


def classify_transcripts(x, y, tree, chunk_size=CHUNK_SIZE):
    """Chunked STRtree query: which transcripts fall inside any cell polygon?

    Returns boolean array is_assigned (True = inside a cell boundary).
    """
    n = len(x)
    is_assigned = np.zeros(n, dtype=bool)
    t0 = time.time()
    n_chunks = (n + chunk_size - 1) // chunk_size

    print(f"  Classifying {n:,} transcripts in {n_chunks} chunks...",
          flush=True)

    for ci, start in enumerate(range(0, n, chunk_size)):
        end = min(start + chunk_size, n)
        coords = np.column_stack([
            x[start:end].astype(np.float64),
            y[start:end].astype(np.float64),
        ])
        pts = points(coords)
        hits = tree.query(pts, predicate='intersects')

        if len(hits[0]) > 0:
            is_assigned[start + hits[0]] = True

        if (ci + 1) % 20 == 0 or ci == 0:
            elapsed = time.time() - t0
            pct_done = 100 * end / n
            rate = end / elapsed if elapsed > 0 else 0
            print(f"    chunk {ci+1}/{n_chunks} "
                  f"({pct_done:.0f}%, {rate/1e6:.1f}M tx/sec)", flush=True)

    elapsed = time.time() - t0
    n_assigned = is_assigned.sum()
    n_unassigned = n - n_assigned
    print(f"\n  Done in {elapsed:.0f}s")
    print(f"  Assigned:   {n_assigned:>12,} ({100*n_assigned/n:.1f}%)")
    print(f"  Unassigned: {n_unassigned:>12,} ({100*n_unassigned/n:.1f}%)")
    return is_assigned


# ══════════════════════════════════════════════════════════════════════
# PART 1: ASSIGNMENT FRACTION
# ══════════════════════════════════════════════════════════════════════

def compute_assignment_stats(x, y, gene_indices, gene_names, is_assigned):
    """Compute per-gene and overall transcript assignment fractions."""
    print("\n" + "=" * 60)
    print("PART 1: Transcript Assignment Fraction")
    print("=" * 60)

    records = []
    for gi, gene in enumerate(gene_names):
        mask = gene_indices == gi
        n_total = mask.sum()
        n_assigned = (mask & is_assigned).sum()
        n_unassigned = n_total - n_assigned
        frac = n_unassigned / n_total if n_total > 0 else 0
        records.append({
            'gene': gene,
            'n_total': int(n_total),
            'n_assigned': int(n_assigned),
            'n_unassigned': int(n_unassigned),
            'frac_unassigned': round(float(frac), 4),
        })

    # Overall
    n_total = len(x)
    n_assigned = int(is_assigned.sum())
    n_unassigned = n_total - n_assigned
    records.append({
        'gene': 'ALL',
        'n_total': n_total,
        'n_assigned': n_assigned,
        'n_unassigned': n_unassigned,
        'frac_unassigned': round(n_unassigned / n_total, 4),
    })

    df = pd.DataFrame(records)
    csv_path = os.path.join(OUT_DIR, f'assignment_stats_{SAMPLE}.csv')
    df.to_csv(csv_path, index=False)
    print(f"\n  Saved: {csv_path}")

    # Summary
    overall = df[df['gene'] == 'ALL'].iloc[0]
    print(f"\n  Overall: {overall['n_assigned']:,} / {overall['n_total']:,} "
          f"assigned ({100*(1-overall['frac_unassigned']):.1f}%)")

    df_genes = df[df['gene'] != 'ALL'].sort_values('frac_unassigned',
                                                      ascending=False)
    print(f"\n  Top 10 genes by unassigned fraction:")
    for _, row in df_genes.head(10).iterrows():
        print(f"    {row['gene']:15s}: {100*row['frac_unassigned']:5.1f}% "
              f"unassigned ({row['n_unassigned']:>10,} / {row['n_total']:>10,})")

    print(f"\n  Bottom 10 genes by unassigned fraction:")
    for _, row in df_genes.tail(10).iterrows():
        print(f"    {row['gene']:15s}: {100*row['frac_unassigned']:5.1f}% "
              f"unassigned ({row['n_unassigned']:>10,} / {row['n_total']:>10,})")

    return df


# ══════════════════════════════════════════════════════════════════════
# PART 2: SPATIAL DENSITY MAPS
# ══════════════════════════════════════════════════════════════════════

def plot_density_comparison(x, y, is_assigned, centroids):
    """2×2 density comparison: all tx, cell, unassigned tx, hotspot ratio."""
    print("\n" + "=" * 60)
    print("PART 2: Spatial Density Maps")
    print("=" * 60)

    # Grid bins
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()
    x_bins = np.arange(x_min, x_max + BIN_SIZE, BIN_SIZE)
    y_bins = np.arange(y_min, y_max + BIN_SIZE, BIN_SIZE)

    # 1. All transcript density
    all_density, _, _ = np.histogram2d(x, y, bins=[x_bins, y_bins])
    # 2. Unassigned transcript density
    ux, uy = x[~is_assigned], y[~is_assigned]
    unassigned_density, _, _ = np.histogram2d(ux, uy, bins=[x_bins, y_bins])
    # 3. Cell density (from centroids)
    cell_density, _, _ = np.histogram2d(
        centroids[:, 0], centroids[:, 1], bins=[x_bins, y_bins])
    # 4. Ghost hotspot ratio
    ghost_ratio = unassigned_density / (cell_density + 1)

    extent = [x_bins[0], x_bins[-1], y_bins[-1], y_bins[0]]

    fig, axes = plt.subplots(2, 2, figsize=(16, 14), facecolor=BG)

    panels = [
        (all_density.T, 'hot', 'All transcript density', None),
        (cell_density.T, 'viridis', 'Cell density', None),
        (unassigned_density.T, 'inferno', 'Unassigned transcript density', None),
        (ghost_ratio.T, 'magma', 'Ghost hotspot ratio\n(unassigned tx / cell)',
         np.percentile(ghost_ratio[ghost_ratio > 0], 95) if (ghost_ratio > 0).any() else 1),
    ]

    for ax, (data, cmap, title, vmax) in zip(axes.flat, panels):
        # Use log scale for density panels, linear for ratio
        if vmax is None:
            im = ax.imshow(data, extent=extent, aspect='equal', cmap=cmap,
                           norm=LogNorm(vmin=max(1, data[data > 0].min()) if (data > 0).any() else 1,
                                        vmax=data.max() if data.max() > 0 else 1),
                           interpolation='nearest')
        else:
            im = ax.imshow(data, extent=extent, aspect='equal', cmap=cmap,
                           vmin=0, vmax=vmax, interpolation='nearest')
        ax.set_title(title, fontsize=18, fontweight='bold', color='white',
                     pad=10)
        ax.set_facecolor('#0a0a0a')
        ax.tick_params(colors='#888888', labelsize=10)
        for spine in ax.spines.values():
            spine.set_color('#444444')
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(colors='#aaaaaa', labelsize=9)

    # Scale bar on first panel
    ax0 = axes[0, 0]
    sb_x = x_min + 100
    sb_y = y_max - 100
    ax0.plot([sb_x, sb_x + 500], [sb_y, sb_y], color='white', lw=3)
    ax0.text(sb_x + 250, sb_y - 50, '500 µm', color='white', fontsize=12,
             ha='center', va='top')

    fig.suptitle(
        f'Ghost Cell Density Analysis — {SAMPLE}\n'
        f'{BIN_SIZE:.0f} µm bins',
        fontsize=22, fontweight='bold', color='white', y=0.98
    )

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    outpath = os.path.join(OUT_DIR, f'density_comparison_{SAMPLE}.png')
    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")

    # Print summary stats
    n_bins_with_ghosts = (ghost_ratio > 10).sum()
    print(f"  Grid: {len(x_bins)-1} × {len(y_bins)-1} bins at {BIN_SIZE:.0f} µm")
    print(f"  Bins with ghost ratio > 10: {n_bins_with_ghosts}")
    print(f"  Max ghost ratio: {ghost_ratio.max():.1f}")
    print(f"  Median ghost ratio (nonzero): "
          f"{np.median(ghost_ratio[ghost_ratio > 0]):.1f}")


# ══════════════════════════════════════════════════════════════════════
# PART 3: DBSCAN CLUSTERING
# ══════════════════════════════════════════════════════════════════════

def find_ghost_clusters(x, y, gene_indices, gene_names, is_assigned):
    """DBSCAN on unassigned transcripts → ghost cell candidate table."""
    print("\n" + "=" * 60)
    print("PART 3: DBSCAN Ghost Cell Clustering")
    print("=" * 60)

    # Extract unassigned coordinates
    unassigned_mask = ~is_assigned
    ux = x[unassigned_mask]
    uy = y[unassigned_mask]
    ugenes = gene_indices[unassigned_mask]
    n_unassigned = len(ux)
    print(f"  Unassigned transcripts: {n_unassigned:,}")

    # Subsample if needed
    if n_unassigned > DBSCAN_MAX_POINTS:
        print(f"  Subsampling to {DBSCAN_MAX_POINTS:,} for DBSCAN...")
        rng = np.random.RandomState(42)
        idx = rng.choice(n_unassigned, DBSCAN_MAX_POINTS, replace=False)
        idx.sort()
        ux_sub = ux[idx]
        uy_sub = uy[idx]
        ugenes_sub = ugenes[idx]
    else:
        ux_sub = ux
        uy_sub = uy
        ugenes_sub = ugenes

    # Run DBSCAN
    print(f"  Running DBSCAN (eps={DBSCAN_EPS}, min_samples={DBSCAN_MIN_SAMPLES}) "
          f"on {len(ux_sub):,} points...", flush=True)
    t0 = time.time()
    coords = np.column_stack([ux_sub, uy_sub])
    db = DBSCAN(eps=DBSCAN_EPS, min_samples=DBSCAN_MIN_SAMPLES,
                n_jobs=-1, algorithm='ball_tree')
    labels = db.fit_predict(coords)
    elapsed = time.time() - t0

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = (labels == -1).sum()
    print(f"  DBSCAN: {n_clusters:,} clusters, {n_noise:,} noise points "
          f"({elapsed:.0f}s)")

    if n_clusters == 0:
        print("  No clusters found — try reducing eps or min_samples")
        return pd.DataFrame(), ux_sub, uy_sub, ugenes_sub, labels

    # Characterize each cluster
    all_marker_genes = set()
    for genes in MARKER_GENES.values():
        all_marker_genes.update(genes)
    gene_name_set = set(gene_names)

    records = []
    for cl in range(n_clusters):
        mask = labels == cl
        cx = ux_sub[mask]
        cy = uy_sub[mask]
        cg = ugenes_sub[mask]
        n_tx = int(mask.sum())

        centroid_x = float(cx.mean())
        centroid_y = float(cy.mean())
        radius = float(np.sqrt(((cx - centroid_x)**2 +
                                 (cy - centroid_y)**2).max()))

        # Gene composition
        gene_counts = {}
        for gi in cg:
            gname = gene_names[gi]
            gene_counts[gname] = gene_counts.get(gname, 0) + 1

        top_genes_sorted = sorted(gene_counts.items(), key=lambda x: -x[1])
        top5 = top_genes_sorted[:5]
        top_genes_str = ', '.join(f'{g}({n})' for g, n in top5)

        # Classify by marker gene overlap
        best_type = 'Unknown'
        best_score = 0
        top10_names = set(g for g, _ in top_genes_sorted[:10])

        for cell_type, markers in MARKER_GENES.items():
            markers_in_panel = [m for m in markers if m in gene_name_set]
            if not markers_in_panel:
                continue
            overlap = len(top10_names & set(markers_in_panel))
            # Also count fraction of cluster transcripts from this type's markers
            marker_tx = sum(gene_counts.get(m, 0) for m in markers_in_panel)
            score = overlap + 0.5 * (marker_tx / n_tx if n_tx > 0 else 0)
            if score > best_score:
                best_score = score
                best_type = cell_type

        # Marker confidence: fraction of cluster transcripts from assigned type
        if best_type != 'Unknown':
            type_markers = [m for m in MARKER_GENES[best_type]
                           if m in gene_name_set]
            marker_tx = sum(gene_counts.get(m, 0) for m in type_markers)
            marker_conf = marker_tx / n_tx if n_tx > 0 else 0
        else:
            marker_conf = 0

        n_unique_genes = len(gene_counts)

        records.append({
            'cluster_id': cl,
            'centroid_x': round(centroid_x, 1),
            'centroid_y': round(centroid_y, 1),
            'n_transcripts': n_tx,
            'n_genes_unique': n_unique_genes,
            'radius_um': round(radius, 1),
            'top_genes': top_genes_str,
            'inferred_type': best_type,
            'marker_confidence': round(marker_conf, 3),
        })

    df = pd.DataFrame(records)

    # Filter to cell-sized clusters
    df_cell_sized = df[df['radius_um'] <= MAX_CLUSTER_RADIUS].copy()
    print(f"\n  Cell-sized clusters (radius ≤ {MAX_CLUSTER_RADIUS} µm): "
          f"{len(df_cell_sized):,} / {len(df):,}")

    # Save full table
    csv_path = os.path.join(OUT_DIR, f'ghost_clusters_{SAMPLE}.csv')
    df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    # Summary
    print(f"\n  Cluster size distribution (all):")
    print(f"    Median:  {df['n_transcripts'].median():.0f} transcripts, "
          f"{df['radius_um'].median():.1f} µm radius")
    print(f"    Max:     {df['n_transcripts'].max()} transcripts, "
          f"{df['radius_um'].max():.1f} µm radius")

    if len(df_cell_sized) > 0:
        type_counts = df_cell_sized['inferred_type'].value_counts()
        print(f"\n  Cell-sized cluster types:")
        for t, n in type_counts.items():
            print(f"    {t:20s}: {n:>5,}")

    return df, ux_sub, uy_sub, ugenes_sub, labels


# ══════════════════════════════════════════════════════════════════════
# PART 3b: NUCLEAR FRACTION VALIDATION
# ══════════════════════════════════════════════════════════════════════

def compute_nuclear_fractions(x, y, gene_indices, gene_names, is_assigned):
    """Classify assigned transcripts as nuclear/cytoplasmic using nucleus polygons.

    Returns per-gene nuclear fraction dict and saves CSV.
    """
    print("\n" + "=" * 60)
    print("PART 3b: Nuclear Fraction Validation")
    print("=" * 60)

    # Load nucleus polygons
    nuc_ids, nuc_polys, _ = load_nucleus_polygons(NUCLEUS_BOUNDARY_PATH)
    print(f"  Building nucleus STRtree ({len(nuc_polys):,} polygons)...",
          flush=True)
    t0 = time.time()
    nuc_tree = STRtree(nuc_polys)
    print(f"  STRtree built in {time.time() - t0:.1f}s", flush=True)

    # Classify assigned transcripts against nucleus polygons
    assigned_mask = is_assigned
    ax = x[assigned_mask]
    ay = y[assigned_mask]
    a_genes = gene_indices[assigned_mask]
    n_assigned_total = len(ax)
    print(f"  Classifying {n_assigned_total:,} assigned transcripts "
          f"against nucleus polygons...", flush=True)

    is_nuclear = np.zeros(n_assigned_total, dtype=bool)
    t0 = time.time()
    n_chunks = (n_assigned_total + CHUNK_SIZE - 1) // CHUNK_SIZE

    for ci, start in enumerate(range(0, n_assigned_total, CHUNK_SIZE)):
        end = min(start + CHUNK_SIZE, n_assigned_total)
        coords = np.column_stack([
            ax[start:end].astype(np.float64),
            ay[start:end].astype(np.float64),
        ])
        pts = points(coords)
        hits = nuc_tree.query(pts, predicate='intersects')
        if len(hits[0]) > 0:
            is_nuclear[start + hits[0]] = True

        if (ci + 1) % 20 == 0 or ci == 0:
            elapsed = time.time() - t0
            pct_done = 100 * end / n_assigned_total
            rate = end / elapsed if elapsed > 0 else 0
            print(f"    chunk {ci+1}/{n_chunks} "
                  f"({pct_done:.0f}%, {rate/1e6:.1f}M tx/sec)", flush=True)

    elapsed = time.time() - t0
    n_nuclear = is_nuclear.sum()
    print(f"\n  Done in {elapsed:.0f}s")
    print(f"  Nuclear:     {n_nuclear:>12,} "
          f"({100*n_nuclear/n_assigned_total:.1f}% of assigned)")
    print(f"  Cytoplasmic: {n_assigned_total - n_nuclear:>12,} "
          f"({100*(n_assigned_total - n_nuclear)/n_assigned_total:.1f}% of assigned)")

    # Compute per-gene nuclear fraction (among assigned transcripts)
    gene_nuc_fracs = {}
    records = []
    for gi, gene in enumerate(gene_names):
        gene_mask = a_genes == gi
        n_gene_assigned = gene_mask.sum()
        if n_gene_assigned == 0:
            gene_nuc_fracs[gene] = 0.0
            continue
        n_gene_nuclear = (gene_mask & is_nuclear).sum()
        frac = float(n_gene_nuclear) / n_gene_assigned
        gene_nuc_fracs[gene] = frac
        records.append({
            'gene': gene,
            'n_assigned': int(n_gene_assigned),
            'n_nuclear': int(n_gene_nuclear),
            'n_cytoplasmic': int(n_gene_assigned - n_gene_nuclear),
            'nuclear_fraction': round(frac, 4),
        })

    nf_df = pd.DataFrame(records).sort_values('nuclear_fraction',
                                                ascending=False)
    csv_path = os.path.join(OUT_DIR,
                             f'gene_nuclear_fractions_{SAMPLE}.csv')
    nf_df.to_csv(csv_path, index=False)
    print(f"\n  Saved: {csv_path}")

    # Print top/bottom genes
    print(f"\n  Top 10 genes by nuclear fraction:")
    for _, row in nf_df.head(10).iterrows():
        print(f"    {row['gene']:15s}: {100*row['nuclear_fraction']:5.1f}% "
              f"nuclear ({row['n_nuclear']:>8,} / {row['n_assigned']:>8,})")

    print(f"\n  Bottom 10 genes by nuclear fraction:")
    for _, row in nf_df.tail(10).iterrows():
        print(f"    {row['gene']:15s}: {100*row['nuclear_fraction']:5.1f}% "
              f"nuclear ({row['n_nuclear']:>8,} / {row['n_assigned']:>8,})")

    return gene_nuc_fracs


def score_ghost_clusters_nuclear(cluster_df, gene_nuc_fracs, gene_names,
                                  ux_sub, uy_sub, ugenes_sub, labels):
    """Add nuclear_score to each ghost cluster from ALL transcripts (not just top 5).

    Uses the raw DBSCAN labels + unassigned transcript data to compute
    accurate weighted-average nuclear fraction for each cluster.
    """
    print("\n  Scoring ghost clusters by nuclear fraction "
          f"(using all genes per cluster)...")

    # Build gene_name → nuclear_frac lookup by gene index
    gene_nf_by_idx = np.zeros(len(gene_names), dtype=np.float64)
    for gi, gn in enumerate(gene_names):
        gene_nf_by_idx[gi] = gene_nuc_fracs.get(gn, 0.0)

    nuclear_scores = []
    n_clusters = len(cluster_df)

    for cl in range(n_clusters):
        mask = labels == cl
        cg = ugenes_sub[mask]
        if len(cg) == 0:
            nuclear_scores.append(0.0)
            continue
        # Weighted average: sum(nuclear_frac_gene * count) / total_count
        # Since each transcript has weight 1, this is just mean(nuclear_frac)
        score = float(gene_nf_by_idx[cg].mean())
        nuclear_scores.append(round(score, 4))

    cluster_df['nuclear_score'] = nuclear_scores

    # Summary by type
    cell_sized = cluster_df[cluster_df['radius_um'] <= MAX_CLUSTER_RADIUS]
    if len(cell_sized) > 0:
        print(f"\n  Nuclear score by inferred type (cell-sized clusters):")
        for t, grp in cell_sized.groupby('inferred_type'):
            print(f"    {t:20s}: median={grp['nuclear_score'].median():.3f}, "
                  f"mean={grp['nuclear_score'].mean():.3f}, n={len(grp)}")

    return cluster_df


# ══════════════════════════════════════════════════════════════════════
# PART 4: EXEMPLAR VISUALIZATION
# ══════════════════════════════════════════════════════════════════════

def add_nearest_cell_distance(cluster_df, centroids):
    """Add nearest_cell_dist column using cKDTree on cell centroids."""
    print("  Computing nearest-cell distance for each cluster...", flush=True)
    cell_tree = cKDTree(centroids)

    dists = []
    for _, row in cluster_df.iterrows():
        d, _ = cell_tree.query([row['centroid_x'], row['centroid_y']], k=1)
        dists.append(round(float(d), 1))
    cluster_df['nearest_cell_dist'] = dists
    return cluster_df


def select_exemplar_clusters(cluster_df, adata, centroids):
    """Pick 6 diverse ghost clusters: 3 low-density + 3 type-diverse."""
    if len(cluster_df) == 0:
        return []

    # Filter to cell-sized with enough transcripts
    df = cluster_df[cluster_df['radius_um'] <= MAX_CLUSTER_RADIUS].copy()
    df = df[df['n_transcripts'] >= 50].copy()

    if len(df) == 0:
        df = cluster_df.nlargest(6, 'n_transcripts')
        return df.to_dict('records')

    # Approximate cortical depth for each cluster centroid
    if adata is not None and 'predicted_norm_depth' in adata.obs.columns:
        coords = adata.obsm['spatial']
        depths = adata.obs['predicted_norm_depth'].values
        for idx in df.index:
            cx, cy = df.loc[idx, 'centroid_x'], df.loc[idx, 'centroid_y']
            dists = (coords[:, 0] - cx)**2 + (coords[:, 1] - cy)**2
            nearest = np.argmin(dists)
            df.loc[idx, 'approx_depth'] = depths[nearest]

    selected = []
    used_ids = set()

    # ── Group A: 3 low-density exemplars ──
    # Prefer clusters far from existing cells, with higher nuclear scores
    if 'nearest_cell_dist' in df.columns:
        low_density = df.sort_values(
            ['nearest_cell_dist', 'nuclear_score' if 'nuclear_score' in df.columns
             else 'n_transcripts'],
            ascending=[False, False])
        for _, row in low_density.iterrows():
            if row['cluster_id'] not in used_ids and len(selected) < 3:
                rec = row.to_dict()
                rec['_density_tag'] = 'LOW DENSITY'
                selected.append(rec)
                used_ids.add(row['cluster_id'])

    # ── Group B: 3 type-diverse exemplars ──
    types_seen = set()
    for _, row in df.sort_values('n_transcripts', ascending=False).iterrows():
        t = row['inferred_type']
        if (row['cluster_id'] not in used_ids and t not in types_seen
                and len(selected) < 6):
            rec = row.to_dict()
            rec['_density_tag'] = 'HIGH DENSITY' if row.get(
                'nearest_cell_dist', 0) < 20 else ''
            selected.append(rec)
            types_seen.add(t)
            used_ids.add(row['cluster_id'])

    # Fill remaining slots
    if len(selected) < 6:
        for _, row in df.sort_values('n_transcripts', ascending=False).iterrows():
            if row['cluster_id'] not in used_ids and len(selected) < 6:
                rec = row.to_dict()
                rec['_density_tag'] = ''
                selected.append(rec)
                used_ids.add(row['cluster_id'])

    return selected


def plot_ghost_exemplars(exemplars, x, y, gene_indices, gene_names,
                         is_assigned, cell_bnd_df, nuc_bnd_df):
    """2×3 panel of ghost cell candidates with boundaries overlaid."""
    print("\n" + "=" * 60)
    print("PART 4: Exemplar Ghost Cell Visualization")
    print("=" * 60)

    if not exemplars:
        print("  No exemplars to plot.")
        return

    n_panels = min(len(exemplars), 6)
    nrows = 2
    ncols = 3
    fig, axes = plt.subplots(nrows, ncols, figsize=(20, 14), facecolor=BG)

    PAD = 40  # µm around cluster centroid

    for i, ex in enumerate(exemplars[:n_panels]):
        r, c = divmod(i, ncols)
        ax = axes[r, c]
        cx, cy = ex['centroid_x'], ex['centroid_y']

        xmin, xmax = cx - PAD, cx + PAD
        ymin, ymax = cy - PAD, cy + PAD

        # ── Assigned transcripts (dim background) ──
        mask_region = ((x >= xmin) & (x <= xmax) &
                       (y >= ymin) & (y <= ymax))
        assigned_in_region = mask_region & is_assigned
        if assigned_in_region.sum() > 0:
            ax.scatter(x[assigned_in_region], y[assigned_in_region],
                       c='#2a2a2a', s=2, alpha=0.15, zorder=1,
                       rasterized=True)

        # ── Unassigned transcripts (colored by gene) ──
        unassigned_in_region = mask_region & ~is_assigned
        if unassigned_in_region.sum() > 0:
            for gi_val in np.unique(gene_indices[unassigned_in_region]):
                gene = gene_names[gi_val]
                gene_mask = unassigned_in_region & (gene_indices == gi_val)
                if gene_mask.sum() == 0:
                    continue
                color = get_gene_color(gene)
                ax.scatter(x[gene_mask], y[gene_mask],
                           c=color, s=18, alpha=0.85, zorder=5,
                           edgecolors='white', linewidths=0.15,
                           rasterized=True)

        # ── Cell boundaries (nearby) ──
        nearby_c = cell_bnd_df[
            cell_bnd_df['vertex_x'].between(xmin, xmax) &
            cell_bnd_df['vertex_y'].between(ymin, ymax)
        ]
        for cid in nearby_c['cell_id'].unique():
            v = cell_bnd_df[cell_bnd_df['cell_id'] == cid][
                ['vertex_x', 'vertex_y']].values
            if len(v) > 2:
                mx, my = v.mean(0)
                if xmin < mx < xmax and ymin < my < ymax:
                    ax.add_patch(MplPolygon(
                        v, closed=True, fill=False,
                        edgecolor='#888888', lw=1.2, alpha=0.6, zorder=3))

        # ── Nucleus boundaries (nearby, dashed) ──
        nearby_n = nuc_bnd_df[
            nuc_bnd_df['vertex_x'].between(xmin, xmax) &
            nuc_bnd_df['vertex_y'].between(ymin, ymax)
        ]
        for cid in nearby_n['cell_id'].unique():
            v = nuc_bnd_df[nuc_bnd_df['cell_id'] == cid][
                ['vertex_x', 'vertex_y']].values
            if len(v) > 2:
                mx, my = v.mean(0)
                if xmin < mx < xmax and ymin < my < ymax:
                    ax.add_patch(MplPolygon(
                        v, closed=True, fill=False,
                        edgecolor='#666666', lw=0.8, linestyle='--',
                        alpha=0.45, zorder=2))

        # ── Circle around ghost cluster centroid ──
        radius = ex.get('radius_um', 10)
        circle = plt.Circle((cx, cy), radius, fill=False,
                             edgecolor='#FF4444', linewidth=2.0,
                             linestyle=':', alpha=0.8, zorder=7)
        ax.add_patch(circle)

        # ── Top gene legend ──
        top_str = ex.get('top_genes', '')
        top_pairs = []
        for part in top_str.split(', ')[:4]:
            if '(' in part:
                gname = part.split('(')[0]
                color = get_gene_color(gname)
                top_pairs.append((gname, color))
        if top_pairs:
            legend_items = [
                Line2D([0], [0], marker='o', color='w',
                       markerfacecolor=c, markersize=7,
                       label=g, linestyle='None')
                for g, c in top_pairs
            ]
            ax.legend(handles=legend_items, loc='upper right', fontsize=8,
                      facecolor='#1a1a2e', edgecolor='#444444',
                      labelcolor='white', framealpha=0.85,
                      handletextpad=0.3, borderpad=0.3)

        # ── Title ──
        inferred = ex.get('inferred_type', 'Unknown')
        n_tx = ex.get('n_transcripts', 0)
        density_tag = ex.get('_density_tag', '')
        title_suffix = f'  [{density_tag}]' if density_tag else ''
        ax.set_title(f'{inferred} ghost ({n_tx} tx){title_suffix}',
                     fontsize=14, fontweight='bold', color='#FF6666', pad=6)

        # ── Info text ──
        depth_str = ''
        if 'approx_depth' in ex:
            depth_str = f' | depth={ex["approx_depth"]:.2f}'
        nuc_str = ''
        if 'nuclear_score' in ex:
            nuc_str = f' | nuc={ex["nuclear_score"]:.3f}'
        dist_str = ''
        if 'nearest_cell_dist' in ex:
            dist_str = f' | cell_dist={ex["nearest_cell_dist"]:.0f}µm'
        info = (f'({cx:.0f}, {cy:.0f}) µm | r={radius:.0f} µm'
                f'{depth_str}{nuc_str}{dist_str}')
        ax.text(0.5, -0.01, info, transform=ax.transAxes, ha='center',
                va='top', fontsize=8, color='#bbbbbb')

        # ── Scale bar ──
        sb_x, sb_y = xmin + 2, ymax - 2
        ax.plot([sb_x, sb_x + 10], [sb_y, sb_y], color='white', lw=2,
                zorder=10)
        ax.text(sb_x + 5, sb_y - 1, '10 µm', color='white', fontsize=7,
                ha='center', va='top', zorder=10)

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymax, ymin)
        ax.set_aspect('equal')
        ax.set_facecolor('#0d0d1a')
        ax.tick_params(left=False, bottom=False, labelleft=False,
                       labelbottom=False)
        for spine in ax.spines.values():
            spine.set_color('#333333')

    # Hide unused panels
    for i in range(n_panels, nrows * ncols):
        r, c = divmod(i, ncols)
        axes[r, c].set_visible(False)

    # Figure legend
    legend_elements = [
        Line2D([0], [0], color='#888888', linewidth=1.5,
               label='Cell boundary'),
        Line2D([0], [0], color='#666666', linewidth=1, linestyle='--',
               label='Nucleus boundary'),
        plt.Circle((0, 0), 0.1, fill=False, edgecolor='#FF4444',
                    linewidth=2, linestyle=':', label='Ghost cluster'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3,
               fontsize=11, facecolor='#1a1a2e', edgecolor='#555555',
               labelcolor='white', bbox_to_anchor=(0.5, 0.005))

    fig.suptitle(
        f'Ghost Cell Candidates — {SAMPLE}\n'
        f'Unassigned transcript clusters with no overlapping cell boundary',
        fontsize=19, fontweight='bold', color='white', y=0.97
    )

    plt.tight_layout(rect=[0, 0.04, 1, 0.93])
    outpath = os.path.join(OUT_DIR, f'ghost_cell_exemplars_{SAMPLE}.png')
    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


# ══════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    print("=" * 60)
    print(f"GHOST CELLS PROTOTYPE — {SAMPLE}")
    print("=" * 60)

    # ── Part 0: Load data ──
    print("\n" + "=" * 60)
    print("PART 0: Data Loading")
    print("=" * 60)

    gene_index = load_gene_index(TRANSCRIPT_DIR)
    cell_ids, polygons, tree, centroids = load_cell_polygons(
        CELL_BOUNDARY_PATH)
    x, y, gene_indices, gene_names = load_all_transcripts(
        TRANSCRIPT_DIR, gene_index)

    # ── Part 1: Assignment fraction ──
    is_assigned = classify_transcripts(x, y, tree)
    stats_df = compute_assignment_stats(x, y, gene_indices, gene_names,
                                         is_assigned)

    # ── Part 2: Density maps ──
    plot_density_comparison(x, y, is_assigned, centroids)

    # ── Part 3: DBSCAN ──
    result = find_ghost_clusters(x, y, gene_indices, gene_names, is_assigned)
    if isinstance(result, tuple):
        cluster_df, ux_sub, uy_sub, ugenes_sub, dbscan_labels = result
    else:
        cluster_df = result
        ux_sub = uy_sub = ugenes_sub = dbscan_labels = None

    # ── Part 3b: Nuclear fraction validation ──
    gene_nuc_fracs = compute_nuclear_fractions(
        x, y, gene_indices, gene_names, is_assigned)

    if len(cluster_df) > 0 and dbscan_labels is not None:
        cluster_df = score_ghost_clusters_nuclear(
            cluster_df, gene_nuc_fracs, gene_names,
            ux_sub, uy_sub, ugenes_sub, dbscan_labels)

        # Add nearest-cell distance
        cluster_df = add_nearest_cell_distance(cluster_df, centroids)

        # Re-save cluster CSV with new columns
        csv_path = os.path.join(OUT_DIR, f'ghost_clusters_{SAMPLE}.csv')
        cluster_df.to_csv(csv_path, index=False)
        print(f"  Updated: {csv_path}")

    # ── Part 4: Exemplar visualization ──
    if len(cluster_df) > 0:
        print("\nLoading boundary DataFrames for visualization...")
        cell_bnd_df = pd.read_csv(CELL_BOUNDARY_PATH, compression='gzip')
        nuc_bnd_df = pd.read_csv(NUCLEUS_BOUNDARY_PATH, compression='gzip')

        print("Loading h5ad for depth info...")
        adata = ad.read_h5ad(H5AD_PATH, backed='r')

        exemplars = select_exemplar_clusters(cluster_df, adata, centroids)
        print(f"  Selected {len(exemplars)} exemplar clusters")

        plot_ghost_exemplars(exemplars, x, y, gene_indices, gene_names,
                             is_assigned, cell_bnd_df, nuc_bnd_df)

        # ── Print strategic summary ──
        print("\n" + "=" * 60)
        print("STRATEGIC SUMMARY")
        print("=" * 60)
        cell_sized = cluster_df[cluster_df['radius_um'] <= MAX_CLUSTER_RADIUS]
        if len(cell_sized) > 0 and 'nuclear_score' in cell_sized.columns:
            print(f"\n  Cell-sized clusters: {len(cell_sized):,}")
            print(f"  Nuclear score distribution:")
            print(f"    25th pctile: {cell_sized['nuclear_score'].quantile(0.25):.3f}")
            print(f"    Median:      {cell_sized['nuclear_score'].median():.3f}")
            print(f"    75th pctile: {cell_sized['nuclear_score'].quantile(0.75):.3f}")
            # High-confidence ghost cells: high nuclear score + far from cells
            if 'nearest_cell_dist' in cell_sized.columns:
                high_conf = cell_sized[
                    (cell_sized['nuclear_score'] > cell_sized['nuclear_score'].median())
                    & (cell_sized['nearest_cell_dist'] > 20)
                ]
                print(f"\n  High-confidence candidates "
                      f"(nuclear_score > median, cell_dist > 20 µm): "
                      f"{len(high_conf):,} ({100*len(high_conf)/len(cell_sized):.1f}%)")
                print(f"  Type breakdown:")
                for t, n in high_conf['inferred_type'].value_counts().items():
                    print(f"    {t:20s}: {n:>5,}")
    else:
        print("\nNo clusters found — skipping exemplar visualization.")

    total = time.time() - t_start
    print(f"\n{'=' * 60}")
    print(f"DONE — {total:.0f}s ({total/60:.1f} min)")
    print(f"Output: {OUT_DIR}")
    print(f"{'=' * 60}")


if __name__ == '__main__':
    main()
