#!/usr/bin/env python3
"""
Compute per-gene fraction of transcripts inside vs outside cell boundaries
for the Br8667 exemplar sample.

Strategy:
1. Load all 67,688 cell boundary polygons from the viewer export
2. Build a union geometry of all cells (or use STRtree for fast spatial query)
3. For each of the 300 genes, load transcript coordinates
4. Test which transcripts fall inside any cell polygon
5. Report per-gene % in-cell vs % ambient (extracellular)

Uses Shapely 2.x vectorized operations for speed.
"""

import os
import sys
import json
import time
import numpy as np
import pandas as pd
from shapely import Polygon, MultiPolygon, STRtree, points, prepare
from shapely.ops import unary_union
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
BOUNDARY_PATH = os.path.join(BASE_DIR, "output", "viewer", "boundaries", "Br8667.json")
TRANSCRIPT_DIR = os.path.join(BASE_DIR, "output", "viewer", "transcripts")
GENE_INDEX_PATH = os.path.join(TRANSCRIPT_DIR, "gene_index.json")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "plots")

t0 = time.time()

# ── 1. Load cell boundaries ──
print("Loading cell boundaries...")
with open(BOUNDARY_PATH, 'r') as f:
    bdata = json.load(f)

n_cells = bdata['n_cells']
verts_per_cell = bdata['verts_per_cell']
bx_offset = bdata['x_offset']
by_offset = bdata['y_offset']
bx_scale = bdata['x_scale']
by_scale = bdata['y_scale']
bx = np.array(bdata['bx'], dtype=np.float64)
by = np.array(bdata['by'], dtype=np.float64)

# Convert to real coordinates (micrometers)
bx_um = bx * bx_scale + bx_offset
by_um = by * by_scale + by_offset

print(f"  {n_cells:,} cells, {verts_per_cell} vertices each")
print(f"  Boundary coord range: X=[{bx_um.min():.1f}, {bx_um.max():.1f}], Y=[{by_um.min():.1f}, {by_um.max():.1f}]")

# Build Shapely polygons
print("Building cell polygons...")
cell_polys = []
n_valid = 0
for i in range(n_cells):
    start = i * verts_per_cell
    end = start + verts_per_cell
    coords = list(zip(bx_um[start:end], by_um[start:end]))
    try:
        poly = Polygon(coords)
        if poly.is_valid and poly.area > 0:
            cell_polys.append(poly)
            n_valid += 1
        else:
            # Try to fix
            poly = poly.buffer(0)
            if poly.is_valid and poly.area > 0:
                cell_polys.append(poly)
                n_valid += 1
    except Exception:
        pass

    if (i + 1) % 20000 == 0:
        print(f"  Built {i+1:,}/{n_cells:,} polygons ({n_valid:,} valid)")

print(f"  {n_valid:,} valid cell polygons out of {n_cells:,}")
print(f"  Polygon construction: {time.time()-t0:.1f}s")

# ── 2. Build spatial index ──
print("\nBuilding STRtree spatial index...")
t1 = time.time()
tree = STRtree(cell_polys)
print(f"  STRtree built in {time.time()-t1:.1f}s")

# ── 3. Load gene index ──
with open(GENE_INDEX_PATH, 'r') as f:
    gindex = json.load(f)

tx_x_offset = gindex['x_offset']
tx_y_offset = gindex['y_offset']
tx_x_scale = gindex['x_scale']
tx_y_scale = gindex['y_scale']

genes = gindex['genes']
n_genes = len(genes)
print(f"\n{n_genes} genes to process")
print(f"  Transcript coord range offsets: X={tx_x_offset}, Y={tx_y_offset}")

# ── 4. Process each gene ──
print("\nProcessing genes...")
results = []
total_in = 0
total_out = 0

for gi, ginfo in enumerate(genes):
    gene_name = ginfo['gene']
    n_transcripts = ginfo['n']

    # Load transcript coordinates
    gene_path = os.path.join(TRANSCRIPT_DIR, ginfo['file'])
    with open(gene_path, 'r') as f:
        gdata = json.load(f)

    tx = np.array(gdata['x'], dtype=np.float64)
    ty = np.array(gdata['y'], dtype=np.float64)

    # Convert to real coordinates (micrometers)
    tx_um = tx * tx_x_scale + tx_x_offset
    ty_um = ty * tx_y_scale + tx_y_offset

    # Create Shapely points
    pts = points(np.column_stack([tx_um, ty_um]))

    # Query STRtree: for each point, find candidate polygons whose bounding box contains it
    # Then do exact containment test
    # Use query_nearest or query method
    # STRtree.query(geometry, predicate) returns indices of geometries that satisfy predicate

    # Process in chunks to manage memory (shapely can handle large batches)
    chunk_size = 500000
    n_inside = 0

    for c_start in range(0, len(pts), chunk_size):
        c_end = min(c_start + chunk_size, len(pts))
        chunk_pts = pts[c_start:c_end]

        # query with 'intersects' predicate returns (input_idx, tree_idx) pairs
        hits = tree.query(chunk_pts, predicate='intersects')

        # Count unique input points that hit at least one polygon
        if len(hits[0]) > 0:
            unique_hits = np.unique(hits[0])
            n_inside += len(unique_hits)

    n_outside = n_transcripts - n_inside
    pct_inside = 100.0 * n_inside / n_transcripts if n_transcripts > 0 else 0

    total_in += n_inside
    total_out += n_outside

    results.append({
        'gene': gene_name,
        'n_transcripts': n_transcripts,
        'n_inside': n_inside,
        'n_outside': n_outside,
        'pct_inside': pct_inside,
        'pct_outside': 100.0 - pct_inside,
    })

    if (gi + 1) % 10 == 0 or gi == 0 or gi == n_genes - 1:
        elapsed = time.time() - t0
        print(f"  [{gi+1:3d}/{n_genes}] {gene_name:20s}: {n_transcripts:>10,} total, "
              f"{pct_inside:5.1f}% in cells | {elapsed:.0f}s")

# ── 5. Summary ──
df = pd.DataFrame(results)
df = df.sort_values('pct_inside', ascending=False).reset_index(drop=True)

total_transcripts = total_in + total_out
print(f"\n{'='*80}")
print(f"SUMMARY: Br8667 transcript capture")
print(f"{'='*80}")
print(f"Total transcripts: {total_transcripts:,}")
print(f"Inside cells:      {total_in:,} ({100*total_in/total_transcripts:.1f}%)")
print(f"Outside cells:     {total_out:,} ({100*total_out/total_transcripts:.1f}%)")

print(f"\n{'='*80}")
print(f"TOP 20 GENES: Highest % inside cells")
print(f"{'='*80}")
print(f"{'Rank':<5} {'Gene':<20} {'N_total':>12} {'%_inside':>10} {'%_outside':>10}")
print("-" * 60)
for i in range(min(20, len(df))):
    r = df.iloc[i]
    print(f"{i+1:<5} {r['gene']:<20} {r['n_transcripts']:>12,} {r['pct_inside']:>9.1f}% {r['pct_outside']:>9.1f}%")

print(f"\n{'='*80}")
print(f"BOTTOM 20 GENES: Lowest % inside cells (most 'ambient')")
print(f"{'='*80}")
print(f"{'Rank':<5} {'Gene':<20} {'N_total':>12} {'%_inside':>10} {'%_outside':>10}")
print("-" * 60)
for i in range(min(20, len(df))):
    r = df.iloc[-(i+1)]
    print(f"{len(df)-i:<5} {r['gene']:<20} {r['n_transcripts']:>12,} {r['pct_inside']:>9.1f}% {r['pct_outside']:>9.1f}%")

# Distribution stats
print(f"\n{'='*80}")
print(f"DISTRIBUTION OF % INSIDE CELLS (across {len(df)} genes)")
print(f"{'='*80}")
print(f"  Mean:   {df['pct_inside'].mean():.1f}%")
print(f"  Median: {df['pct_inside'].median():.1f}%")
print(f"  Std:    {df['pct_inside'].std():.1f}%")
print(f"  Min:    {df['pct_inside'].min():.1f}% ({df.iloc[-1]['gene']})")
print(f"  Max:    {df['pct_inside'].max():.1f}% ({df.iloc[0]['gene']})")
print(f"  Q25:    {df['pct_inside'].quantile(0.25):.1f}%")
print(f"  Q75:    {df['pct_inside'].quantile(0.75):.1f}%")

# Save results
out_csv = os.path.join(OUTPUT_DIR, "transcript_capture_by_gene.csv")
df.to_csv(out_csv, index=False)
print(f"\nResults saved to: {out_csv}")

elapsed = time.time() - t0
print(f"\nTotal elapsed: {elapsed:.0f}s ({elapsed/60:.1f} min)")
