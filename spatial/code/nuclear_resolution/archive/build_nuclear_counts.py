#!/usr/bin/env python3
"""
Build a nuclear-only count matrix for Br8667 by assigning each transcript
molecule to the nucleus polygon it falls within.

Hypothesis: Spatial doublets arise from cytoplasmic mRNA spillover between
neighboring cells. Nuclear-only counts should be much cleaner.

Approach:
1. Parse nucleus boundary polygons from CSV (68,976 nuclei)
2. Build Shapely STRtree spatial index
3. For each gene (300), load transcript coordinates, query STRtree
4. Accumulate per-cell nuclear counts into a count matrix
5. Save as AnnData with same obs/var structure as whole-cell h5ad

Output:
  output/h5ad/Br8667_nuclear_counts.h5ad
  output/presentation/nuclear_counting_stats.csv
"""

import os
import sys
import json
import time
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from shapely import Polygon, STRtree, points
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import H5AD_DIR, PRESENTATION_DIR

# ── Config ──
SAMPLE = 'Br8667'
H5AD_PATH = os.path.join(H5AD_DIR, f'{SAMPLE}_annotated.h5ad')
NUCLEUS_PATH = os.path.expanduser(
    f'~/Github/SCZ_Xenium/data/raw/GSM9223485_{SAMPLE}-nucleus_boundaries.csv.gz')
TRANSCRIPT_DIR = os.path.expanduser(
    f'~/Github/SCZ_Xenium/output/viewer/transcripts/{SAMPLE}')
GENE_INDEX_PATH = os.path.join(TRANSCRIPT_DIR, 'gene_index.json')
OUTPUT_H5AD = os.path.join(H5AD_DIR, f'{SAMPLE}_nuclear_counts.h5ad')
OUTPUT_STATS = os.path.join(PRESENTATION_DIR, 'nuclear_counting_stats.csv')

CHUNK_SIZE = 500_000  # transcripts per STRtree query batch


def load_nucleus_polygons(nuc_path):
    """Load nucleus boundary CSV and build Shapely polygons per cell.

    Returns
    -------
    cell_ids : list of str
        Unique cell IDs in order
    polygons : list of shapely.Polygon
        One valid polygon per cell (same order as cell_ids)
    cell_id_to_poly_idx : dict
        {cell_id: index in polygons list}
    """
    print("Loading nucleus boundaries...")
    t0 = time.time()
    df = pd.read_csv(nuc_path, compression='gzip')
    print(f"  {len(df):,} vertices, {df['cell_id'].nunique():,} unique cells")

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

        if len(cell_ids) % 20000 == 0 and len(cell_ids) > 0:
            print(f"  Built {len(cell_ids):,} polygons...")

    cell_id_to_poly_idx = {cid: i for i, cid in enumerate(cell_ids)}

    elapsed = time.time() - t0
    print(f"  {len(polygons):,} valid polygons ({n_invalid} invalid) in {elapsed:.1f}s")
    return cell_ids, polygons, cell_id_to_poly_idx


def build_nuclear_count_matrix(polygons, cell_ids, cell_id_to_poly_idx,
                                adata_ref, transcript_dir, gene_index):
    """Build nuclear-only count matrix via point-in-polygon queries.

    Parameters
    ----------
    polygons : list of Polygon
    cell_ids : list of str
    cell_id_to_poly_idx : dict
    adata_ref : AnnData
        Reference h5ad for canonical cell/gene ordering
    transcript_dir : str
    gene_index : dict
        Parsed gene_index.json

    Returns
    -------
    count_matrix : scipy.sparse.csr_matrix
        (n_cells, n_genes) nuclear transcript counts
    gene_stats : list of dict
        Per-gene counting statistics
    """
    print("\nBuilding STRtree spatial index...")
    t0 = time.time()
    tree = STRtree(polygons)
    print(f"  STRtree built in {time.time() - t0:.1f}s")

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
    print(f"  Mapped {n_mapped:,}/{n_polys:,} polygons to h5ad rows")

    n_cells = adata_ref.n_obs
    n_genes = adata_ref.n_vars
    gene_names = list(adata_ref.var_names)
    gene_name_to_col = {g: i for i, g in enumerate(gene_names)}

    # Dense count matrix (68,976 × 300 × 4 bytes = ~79 MB)
    count_matrix = np.zeros((n_cells, n_genes), dtype=np.int32)

    # Transcript coordinate conversion params
    x_scale = gene_index['x_scale']
    y_scale = gene_index['y_scale']
    x_offset = gene_index['x_offset']
    y_offset = gene_index['y_offset']

    gene_stats = []
    total_nuclear = 0
    total_transcripts = 0

    print(f"\nProcessing {len(gene_index['genes'])} genes...")
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
        n_nuclear_gene = 0
        all_h5ad_rows = []

        for c_start in range(0, len(tx), CHUNK_SIZE):
            c_end = min(c_start + CHUNK_SIZE, len(tx))
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

        if (gi + 1) % 25 == 0 or gi == 0:
            elapsed = time.time() - t_start
            rate = (gi + 1) / elapsed * 60
            print(f"  [{gi+1:3d}/{len(gene_index['genes'])}] {gene_name:12s}: "
                  f"{n_nuclear_gene:>10,} / {n_total:>10,} nuclear ({pct:5.1f}%) "
                  f"[{rate:.0f} genes/min]")

    elapsed = time.time() - t_start
    pct_total = 100.0 * total_nuclear / total_transcripts if total_transcripts > 0 else 0
    print(f"\n  Done in {elapsed:.0f}s ({elapsed/60:.1f} min)")
    print(f"  Total: {total_nuclear:,} / {total_transcripts:,} nuclear ({pct_total:.1f}%)")

    # Convert to sparse
    count_sparse = sparse.csr_matrix(count_matrix)
    return count_sparse, gene_stats


def main():
    t_global = time.time()
    print("=" * 70)
    print(f"Nuclear Count Matrix Builder ({SAMPLE})")
    print("=" * 70)

    # 1. Load nucleus polygons
    cell_ids, polygons, cell_id_to_poly_idx = load_nucleus_polygons(NUCLEUS_PATH)

    # 2. Load reference h5ad
    print(f"\nLoading reference h5ad: {H5AD_PATH}")
    adata_ref = ad.read_h5ad(H5AD_PATH)
    print(f"  {adata_ref.n_obs:,} cells × {adata_ref.n_vars} genes")

    # Check overlap
    ref_ids = set(adata_ref.obs_names)
    poly_ids = set(cell_ids)
    overlap = ref_ids & poly_ids
    print(f"  Cell ID overlap: {len(overlap):,} / {len(ref_ids):,} h5ad cells "
          f"({100*len(overlap)/len(ref_ids):.1f}%)")

    # 3. Load gene index
    with open(GENE_INDEX_PATH) as f:
        gene_index = json.load(f)
    print(f"  Gene index: {gene_index['n_genes']} genes, "
          f"{gene_index['total_transcripts']:,} total transcripts")

    # 4. Build nuclear count matrix
    nuc_counts, gene_stats = build_nuclear_count_matrix(
        polygons, cell_ids, cell_id_to_poly_idx,
        adata_ref, TRANSCRIPT_DIR, gene_index)

    # 5. Create nuclear AnnData
    print("\nBuilding nuclear AnnData...")
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

    # Summary stats
    qc_mask = adata_nuc.obs['qc_pass'].values.astype(bool)
    nf = adata_nuc.obs.loc[qc_mask, 'nuclear_fraction']
    nt = adata_nuc.obs.loc[qc_mask, 'nuclear_total_counts']
    print(f"\n  QC-pass cells: {qc_mask.sum():,}")
    print(f"  Nuclear fraction: median={nf.median():.3f}, "
          f"mean={nf.mean():.3f}, IQR=[{nf.quantile(0.25):.3f}, {nf.quantile(0.75):.3f}]")
    print(f"  Nuclear UMI: median={nt.median():.0f}, "
          f"mean={nt.mean():.0f}, min={nt.min()}, max={nt.max()}")
    print(f"  Cells with <50 nuclear UMI: {(nt < 50).sum():,} "
          f"({100*(nt < 50).sum()/qc_mask.sum():.1f}%)")

    # 6. Save
    print(f"\nSaving: {OUTPUT_H5AD}")
    adata_nuc.write_h5ad(OUTPUT_H5AD)
    print(f"  Size: {os.path.getsize(OUTPUT_H5AD) / 1e6:.1f} MB")

    # 7. Save per-gene stats
    stats_df = pd.DataFrame(gene_stats)
    stats_df = stats_df.sort_values('pct_nuclear', ascending=False)
    stats_df.to_csv(OUTPUT_STATS, index=False)
    print(f"Saved: {OUTPUT_STATS}")

    # Print top/bottom genes by nuclear fraction
    print("\n  Top 10 most nuclear genes:")
    for _, row in stats_df.head(10).iterrows():
        print(f"    {row['gene']:12s}: {row['pct_nuclear']:5.1f}% "
              f"({row['n_nuclear']:>10,} / {row['n_total_transcripts']:>10,})")
    print("\n  Bottom 10 least nuclear genes:")
    for _, row in stats_df.tail(10).iterrows():
        print(f"    {row['gene']:12s}: {row['pct_nuclear']:5.1f}% "
              f"({row['n_nuclear']:>10,} / {row['n_total_transcripts']:>10,})")

    elapsed = time.time() - t_global
    print(f"\nTotal elapsed: {elapsed:.0f}s ({elapsed/60:.1f} min)")


if __name__ == '__main__':
    main()
