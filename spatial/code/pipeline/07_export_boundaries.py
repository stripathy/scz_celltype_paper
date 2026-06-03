#!/usr/bin/env python3
"""
Step 7: Export cell and nucleus boundary polygons for the viewer.

Reads the raw cell_boundaries.csv.gz and nucleus_boundaries.csv.gz files
for ALL cells (including QC-fail) and writes compact JSON files that the
browser viewer can load on demand when zoomed in.

Each cell has exactly 25 boundary vertices (Xenium default for both cell
and nucleus segmentation). Vertices are quantized to uint16 for compact storage.

Output:
  - output/viewer/boundaries/<sample_id>.json          (cell boundaries)
  - output/viewer/boundaries/<sample_id>_nucleus.json  (nucleus boundaries)
"""

import os
import sys
import csv
import gzip
import json
import time
import glob
import numpy as np
import anndata as ad

# Pipeline config
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_config import RAW_DIR, H5AD_DIR, VIEWER_DIR

OUTPUT_DIR = os.path.join(VIEWER_DIR, "boundaries")
VERTS_PER_CELL = 25
QUANT_RESOLUTION = 0.2  # micrometers per integer step


def _read_boundary_csv(boundary_path, cell_id_set):
    """Read boundary CSV and return dict: cell_id -> [(x, y), ...]."""
    cells = {}
    with gzip.open(boundary_path, 'rt') as f:
        reader = csv.DictReader(f)
        for row in reader:
            cid = row['cell_id']
            if cid not in cell_id_set:
                continue
            if cid not in cells:
                cells[cid] = []
            cells[cid].append((float(row['vertex_x']), float(row['vertex_y'])))
    return cells


def _build_quantized_arrays(cells, cell_id_order, x_min, y_min):
    """Build flat quantized coordinate arrays aligned to cell order."""
    bx = []
    by = []
    n_with_boundary = 0

    for cid in cell_id_order:
        if cid in cells:
            verts = cells[cid]
            n_with_boundary += 1
        else:
            # No boundary data — use centroid placeholder (will appear as a dot)
            verts = [(0, 0)]

        # Pad or truncate to exactly VERTS_PER_CELL
        while len(verts) < VERTS_PER_CELL:
            verts.append(verts[-1])
        verts = verts[:VERTS_PER_CELL]

        for vx, vy in verts:
            bx.append(int(round((vx - x_min) / QUANT_RESOLUTION)))
            by.append(int(round((vy - y_min) / QUANT_RESOLUTION)))

    return bx, by, n_with_boundary


def _write_boundary_json(sample_id, n_cells, n_with_boundary,
                         x_min, y_min, bx, by, fpath, boundary_type="cell"):
    """Write boundary JSON file."""
    data = {
        "sample_id": sample_id,
        "boundary_type": boundary_type,
        "n_cells": n_cells,
        "n_with_boundary": n_with_boundary,
        "verts_per_cell": VERTS_PER_CELL,
        "x_offset": float(x_min),
        "y_offset": float(y_min),
        "x_scale": QUANT_RESOLUTION,
        "y_scale": QUANT_RESOLUTION,
        "bx": bx,
        "by": by,
    }

    with open(fpath, 'w') as f:
        json.dump(data, f, separators=(',', ':'))

    return os.path.getsize(fpath) / 1024 / 1024


def process_sample(sample_id, raw_dir, h5ad_dir, output_dir):
    """Export cell and nucleus boundary polygons for one sample."""

    # Find cell boundary file
    cell_boundary_files = glob.glob(
        os.path.join(raw_dir, f"*{sample_id}-cell_boundaries.csv.gz"))
    if not cell_boundary_files:
        print(f"  {sample_id}: no cell boundary file found, skipping")
        return None

    # Find nucleus boundary file (optional)
    nucleus_boundary_files = glob.glob(
        os.path.join(raw_dir, f"*{sample_id}-nucleus_boundaries.csv.gz"))
    has_nucleus = len(nucleus_boundary_files) > 0

    # Read h5ad to get QC-pass cell IDs in order
    h5ad_path = os.path.join(h5ad_dir, f"{sample_id}_annotated.h5ad")
    if not os.path.exists(h5ad_path):
        print(f"  {sample_id}: no h5ad file found, skipping")
        return None

    adata = ad.read_h5ad(h5ad_path)
    # Include ALL cells (not just QC-pass) so the viewer can render
    # boundaries for QC-failed cells too (shown as white)
    all_cell_ids = adata.obs_names.tolist()
    all_set = set(all_cell_ids)
    n_qc = len(all_cell_ids)

    # ── Cell boundaries ────────────────────────────────────────────
    cell_polys = _read_boundary_csv(cell_boundary_files[0], all_set)

    # Compute global coordinate range (from cell boundaries)
    all_x = [v[0] for vs in cell_polys.values() for v in vs]
    all_y = [v[1] for vs in cell_polys.values() for v in vs]
    x_min = min(all_x)
    y_min = min(all_y)
    x_max = max(all_x)
    y_max = max(all_y)

    # Verify uint16 range
    x_range_q = (x_max - x_min) / QUANT_RESOLUTION
    y_range_q = (y_max - y_min) / QUANT_RESOLUTION
    assert x_range_q < 65535, f"X range too large for uint16: {x_range_q}"
    assert y_range_q < 65535, f"Y range too large for uint16: {y_range_q}"

    bx, by, n_cell_matched = _build_quantized_arrays(
        cell_polys, all_cell_ids, x_min, y_min)

    cell_fpath = os.path.join(output_dir, f"{sample_id}.json")
    cell_size = _write_boundary_json(
        sample_id, n_qc, n_cell_matched, x_min, y_min, bx, by,
        cell_fpath, boundary_type="cell")

    # ── Nucleus boundaries ─────────────────────────────────────────
    nucleus_size = 0
    n_nuc_matched = 0
    if has_nucleus:
        nuc_polys = _read_boundary_csv(nucleus_boundary_files[0], all_set)

        nbx, nby, n_nuc_matched = _build_quantized_arrays(
            nuc_polys, all_cell_ids, x_min, y_min)

        nuc_fpath = os.path.join(output_dir, f"{sample_id}_nucleus.json")
        nucleus_size = _write_boundary_json(
            sample_id, n_qc, n_nuc_matched, x_min, y_min, nbx, nby,
            nuc_fpath, boundary_type="nucleus")

    nuc_str = f", {n_nuc_matched:,} nuclei ({nucleus_size:.1f} MB)" if has_nucleus else ""
    print(f"  {sample_id}: {n_qc:,} cells, {n_cell_matched:,} cell boundaries "
          f"({cell_size:.1f} MB){nuc_str}")

    return cell_size + nucleus_size


def main():
    t0 = time.time()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Find all annotated h5ad files
    fnames = sorted(f for f in os.listdir(H5AD_DIR) if f.endswith("_annotated.h5ad"))
    sample_ids = [f.replace("_annotated.h5ad", "") for f in fnames]

    print(f"Exporting cell + nucleus boundaries for {len(sample_ids)} samples")
    print(f"Output: {OUTPUT_DIR}/\n")

    total_size = 0
    n_exported = 0

    for sid in sample_ids:
        fsize = process_sample(sid, RAW_DIR, H5AD_DIR, OUTPUT_DIR)
        if fsize is not None:
            total_size += fsize
            n_exported += 1

    elapsed = time.time() - t0
    print(f"\nDone! {n_exported} samples exported in {elapsed:.1f}s")
    print(f"Total size: {total_size:.1f} MB")
    print(f"Output: {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
