"""
Nuclear-only transcript counting via point-in-polygon spatial queries.

Assigns each transcript molecule to the nucleus polygon it falls within,
producing a nuclear-only count matrix. Used by optional nuclear doublet
resolution (see 04_run_nuclear_doublet_resolution.py).

Key functions:
  load_nucleus_polygons()       — Parse boundary CSV → Shapely polygons
  build_nuclear_count_matrix()  — STRtree spatial query → sparse count matrix
  build_nuclear_adata()         — Wrap counts + metadata into AnnData

Requires shapely >= 2.0 for vectorized STRtree.query() with predicate support.
"""

import os
import json
import time
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from shapely import Polygon, STRtree, points


# ═══════════════════════════════════════════════════════════════════════
# 1. NUCLEUS POLYGON LOADING
# ═══════════════════════════════════════════════════════════════════════

def load_nucleus_polygons(nuc_path):
    """Load nucleus boundary CSV and build Shapely polygons per cell.

    Parameters
    ----------
    nuc_path : str
        Path to nucleus boundary CSV (gzipped). Expected columns:
        cell_id, vertex_x, vertex_y.

    Returns
    -------
    cell_ids : list of str
        Unique cell IDs in order
    polygons : list of shapely.Polygon
        One valid polygon per cell (same order as cell_ids)
    cell_id_to_poly_idx : dict
        {cell_id: index in polygons list}
    """
    print("  Loading nucleus boundaries...", flush=True)
    t0 = time.time()
    df = pd.read_csv(nuc_path, compression='gzip')
    print(f"    {len(df):,} vertices, {df['cell_id'].nunique():,} unique cells",
          flush=True)

    # Group by cell_id preserving order of first appearance
    grouped = df.groupby('cell_id', sort=False)

    cell_ids = []
    polygons = []
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
            else:
                n_invalid += 1
        except Exception:
            n_invalid += 1

    cell_id_to_poly_idx = {cid: i for i, cid in enumerate(cell_ids)}

    elapsed = time.time() - t0
    print(f"    {len(polygons):,} valid polygons ({n_invalid} invalid) "
          f"in {elapsed:.1f}s", flush=True)
    return cell_ids, polygons, cell_id_to_poly_idx


# ═══════════════════════════════════════════════════════════════════════
# 2. NUCLEAR COUNT MATRIX BUILDING
# ═══════════════════════════════════════════════════════════════════════

def build_nuclear_count_matrix(polygons, cell_ids, cell_id_to_poly_idx,
                                adata_ref, transcript_dir, gene_index,
                                chunk_size=500_000):
    """Build nuclear-only count matrix via point-in-polygon queries.

    For each gene, loads transcript coordinates from per-gene JSON files,
    creates Shapely points, queries the STRtree spatial index, and accumulates
    per-cell counts via np.bincount.

    Parameters
    ----------
    polygons : list of Polygon
        Nucleus polygons (from load_nucleus_polygons).
    cell_ids : list of str
        Cell IDs matching polygon order.
    cell_id_to_poly_idx : dict
        {cell_id: polygon index}.
    adata_ref : AnnData
        Reference h5ad for canonical cell/gene ordering.
    transcript_dir : str
        Directory containing per-gene JSON files + gene_index.json.
    gene_index : dict
        Parsed gene_index.json (has x_scale, y_scale, x_offset, y_offset, genes).
    chunk_size : int
        Transcripts per STRtree query batch (default 500K).

    Returns
    -------
    count_matrix : scipy.sparse.csr_matrix
        (n_cells, n_genes) nuclear transcript counts.
    gene_stats : list of dict
        Per-gene counting statistics.
    """
    print("  Building STRtree spatial index...", flush=True)
    t0 = time.time()
    tree = STRtree(polygons)
    print(f"    STRtree built in {time.time() - t0:.1f}s", flush=True)

    # Mapping: polygon index → h5ad row index
    obs_names = list(adata_ref.obs_names)
    obs_name_to_row = {name: i for i, name in enumerate(obs_names)}

    # Build poly_idx → h5ad_row lookup array
    n_polys = len(polygons)
    poly_to_h5ad = np.full(n_polys, -1, dtype=np.int64)
    n_mapped = 0
    for cid, pidx in cell_id_to_poly_idx.items():
        if cid in obs_name_to_row:
            poly_to_h5ad[pidx] = obs_name_to_row[cid]
            n_mapped += 1
    print(f"    Mapped {n_mapped:,}/{n_polys:,} polygons to h5ad rows", flush=True)

    n_cells = adata_ref.n_obs
    n_genes = adata_ref.n_vars
    gene_names = list(adata_ref.var_names)
    gene_name_to_col = {g: i for i, g in enumerate(gene_names)}

    # Dense count matrix (n_cells × n_genes × 4 bytes)
    count_matrix = np.zeros((n_cells, n_genes), dtype=np.int32)

    # Transcript coordinate conversion params
    x_scale = gene_index['x_scale']
    y_scale = gene_index['y_scale']
    x_offset = gene_index['x_offset']
    y_offset = gene_index['y_offset']

    gene_stats = []
    total_nuclear = 0
    total_transcripts = 0

    n_genes_total = len(gene_index['genes'])
    print(f"  Processing {n_genes_total} genes...", flush=True)
    t_start = time.time()

    for gi, ginfo in enumerate(gene_index['genes']):
        gene_name = ginfo['gene']
        n_total = ginfo['n']

        # Skip genes not in h5ad
        if gene_name not in gene_name_to_col:
            continue
        gene_col = gene_name_to_col[gene_name]

        # Load transcript coordinates
        gene_path = os.path.join(transcript_dir, ginfo['file'])
        with open(gene_path) as f:
            gdata = json.load(f)

        tx = np.array(gdata['x'], dtype=np.float64) * x_scale + x_offset
        ty = np.array(gdata['y'], dtype=np.float64) * y_scale + y_offset

        # Process in chunks
        all_h5ad_rows = []

        for c_start in range(0, len(tx), chunk_size):
            c_end = min(c_start + chunk_size, len(tx))
            coords = np.column_stack([tx[c_start:c_end], ty[c_start:c_end]])
            pts = points(coords)

            # Query STRtree: returns (point_idx, polygon_idx) pairs
            hits = tree.query(pts, predicate='intersects')

            if len(hits[0]) > 0:
                poly_indices = hits[1]
                h5ad_rows = poly_to_h5ad[poly_indices]
                # Filter out unmapped polygons (-1)
                valid = h5ad_rows >= 0
                if valid.any():
                    all_h5ad_rows.append(h5ad_rows[valid])

        # Accumulate counts for this gene
        n_nuclear_gene = 0
        if all_h5ad_rows:
            all_rows = np.concatenate(all_h5ad_rows)
            counts = np.bincount(all_rows, minlength=n_cells)
            count_matrix[:, gene_col] = counts.astype(np.int32)
            n_nuclear_gene = int(counts.sum())

        pct = 100.0 * n_nuclear_gene / n_total if n_total > 0 else 0
        gene_stats.append({
            'gene': gene_name,
            'n_total_transcripts': n_total,
            'n_nuclear': n_nuclear_gene,
            'pct_nuclear': round(pct, 2),
        })

        total_nuclear += n_nuclear_gene
        total_transcripts += n_total

        if (gi + 1) % 50 == 0 or gi == 0:
            elapsed = time.time() - t_start
            rate = (gi + 1) / elapsed * 60
            print(f"    [{gi+1:3d}/{n_genes_total}] {gene_name:12s}: "
                  f"{n_nuclear_gene:>10,} / {n_total:>10,} nuclear ({pct:5.1f}%) "
                  f"[{rate:.0f} genes/min]", flush=True)

    elapsed = time.time() - t_start
    pct_total = (100.0 * total_nuclear / total_transcripts
                 if total_transcripts > 0 else 0)
    print(f"    Done in {elapsed:.0f}s ({elapsed/60:.1f} min) | "
          f"Total: {total_nuclear:,}/{total_transcripts:,} nuclear "
          f"({pct_total:.1f}%)", flush=True)

    # Convert to sparse
    count_sparse = sparse.csr_matrix(count_matrix)
    return count_sparse, gene_stats


# ═══════════════════════════════════════════════════════════════════════
# 3. ANNDATA CONSTRUCTION
# ═══════════════════════════════════════════════════════════════════════

def build_nuclear_adata(nuc_counts, adata_ref):
    """Build AnnData from nuclear count matrix with metadata from reference.

    Parameters
    ----------
    nuc_counts : scipy.sparse.csr_matrix
        (n_cells, n_genes) nuclear transcript counts.
    adata_ref : AnnData
        Reference h5ad (whole-cell) for obs/var/obsm metadata.

    Returns
    -------
    adata_nuc : AnnData
        Nuclear count AnnData with added columns:
        nuclear_total_counts, nuclear_n_genes, nuclear_fraction.
    """
    adata_nuc = ad.AnnData(
        X=nuc_counts,
        obs=adata_ref.obs.copy(),
        var=adata_ref.var.copy(),
    )
    if 'spatial' in adata_ref.obsm:
        adata_nuc.obsm['spatial'] = adata_ref.obsm['spatial'].copy()

    # Add nuclear-specific metadata
    nuc_total = np.array(nuc_counts.sum(axis=1)).flatten()
    nuc_ngenes = np.array((nuc_counts > 0).sum(axis=1)).flatten()
    wc_total = adata_ref.obs['total_counts'].values.astype(float)

    adata_nuc.obs['nuclear_total_counts'] = nuc_total.astype(int)
    adata_nuc.obs['nuclear_n_genes'] = nuc_ngenes.astype(int)
    adata_nuc.obs['nuclear_fraction'] = np.where(
        wc_total > 0, nuc_total / wc_total, 0.0).astype(np.float32)

    return adata_nuc
