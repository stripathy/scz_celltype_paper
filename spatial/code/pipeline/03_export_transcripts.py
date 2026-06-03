#!/usr/bin/env python3
"""
Step 3: Export per-gene transcript molecule coordinates.

Reads Xenium transcripts.zarr (or .zarr.zip) files for one or more samples,
extracts (x, y) coordinates for each gene, and writes compact JSON files.
Used by optional nuclear doublet resolution (see nuclear_resolution/) and
by the browser-based viewer.

Output (per sample):
  - output/viewer/transcripts/<SAMPLE_ID>/gene_index.json
  - output/viewer/transcripts/<SAMPLE_ID>/<GENE>.json

Usage:
  python 03_export_transcripts.py                # process all available samples
  python 03_export_transcripts.py Br8667 Br6032  # process specific samples
"""

import os
import sys
import json
import time
import glob
import re
import numpy as np

# Pipeline config
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_config import RAW_DIR, VIEWER_DIR

# ── Configuration ──
TRANSCRIPTS_BASE_DIR = os.path.join(VIEWER_DIR, "transcripts")

# Quantization: encode coords as uint16 integers for compact JSON
# Resolution = 0.2 um (uint16 max = 65535 -> covers 0-13107 um range, plenty)
QUANT_RESOLUTION = 0.2  # micrometers per integer step

# Map sample IDs to their zarr paths (auto-discovered)
def discover_zarr_files():
    """Find all transcript zarr files (both .zarr dirs and .zarr.zip archives)."""
    samples = {}

    # Search for .zarr directories
    for p in glob.glob(os.path.join(RAW_DIR, "**/*-transcripts.zarr"), recursive=True):
        m = re.search(r'(Br\d+)-transcripts\.zarr$', p)
        if m:
            samples[m.group(1)] = p

    # Search for .zarr.zip files (override dirs if both exist — zips are newer)
    for p in glob.glob(os.path.join(RAW_DIR, "**/*-transcripts.zarr.zip"), recursive=True):
        m = re.search(r'(Br\d+)-transcripts\.zarr\.zip$', p)
        if m:
            samples[m.group(1)] = p

    return samples


def process_sample(sample_id, zarr_path):
    """Process a single sample's transcript zarr into per-gene JSON files."""
    t0 = time.time()

    import zarr

    output_dir = os.path.join(TRANSCRIPTS_BASE_DIR, sample_id)

    # Check if already processed
    index_path = os.path.join(output_dir, "gene_index.json")
    if os.path.exists(index_path):
        print(f"  [{sample_id}] Already exported — skipping (delete {index_path} to re-export)")
        return

    print(f"\n{'='*70}")
    print(f"  Processing {sample_id}: {zarr_path}")
    print(f"{'='*70}")

    # Open zarr — handle both directory and zip stores
    is_zip = zarr_path.endswith('.zarr.zip')
    if is_zip:
        store = zarr.storage.ZipStore(zarr_path, mode='r')
    else:
        store = zarr.storage.LocalStore(zarr_path)
    z = zarr.open_group(store, mode='r')

    gene_names = list(z.attrs['gene_names'])
    total_rnas = z.attrs['number_rnas']
    print(f"  Total RNAs: {total_rnas:,}")
    print(f"  Total gene entries: {len(gene_names)}")

    # Identify real genes vs control probes
    real_gene_mask = [
        'Unassigned' not in g and 'BLANK' not in g and 'NegControl' not in g
        for g in gene_names
    ]
    real_gene_indices = [i for i, is_real in enumerate(real_gene_mask) if is_real]
    real_gene_names = [gene_names[i] for i in real_gene_indices]
    print(f"  Real genes: {len(real_gene_names)} / {len(gene_names)}")

    # List level-0 grid cells (finest resolution)
    # For zip stores, use zarr API to list groups instead of filesystem
    if is_zip:
        grid0 = z['grids']['0']
        grid_cells = sorted([
            k for k in grid0.keys()
            if ',' in k
        ])
    else:
        grid0_path = os.path.join(zarr_path, 'grids', '0')
        grid_cells = sorted([
            g for g in os.listdir(grid0_path)
            if ',' in g and not g.startswith('.')
        ])
    print(f"  Level-0 grid cells: {len(grid_cells)}")

    # ── Pass 1: Read all transcripts, collect per-gene coordinates ──
    print(f"  Reading all transcripts from {len(grid_cells)} grid cells...")

    gene_idx_to_list_idx = {}
    for list_idx, gene_idx in enumerate(real_gene_indices):
        gene_idx_to_list_idx[gene_idx] = list_idx

    gene_x_coords = [[] for _ in real_gene_indices]
    gene_y_coords = [[] for _ in real_gene_indices]

    total_read = 0
    total_real = 0

    for ci, gc in enumerate(grid_cells):
        cell = z['grids']['0'][gc]
        gi = np.asarray(cell['gene_identity']).ravel()
        loc = np.asarray(cell['location'])

        total_read += len(gi)

        for gene_idx in np.unique(gi):
            gene_idx_int = int(gene_idx)
            if gene_idx_int not in gene_idx_to_list_idx:
                continue

            list_idx = gene_idx_to_list_idx[gene_idx_int]
            mask = gi == gene_idx
            n_hits = mask.sum()
            total_real += n_hits

            gene_x_coords[list_idx].append(loc[mask, 0])
            gene_y_coords[list_idx].append(loc[mask, 1])

        if (ci + 1) % 100 == 0 or ci == len(grid_cells) - 1:
            elapsed = time.time() - t0
            print(f"    [{ci+1}/{len(grid_cells)}] {total_read:,} transcripts read "
                  f"({total_real:,} real) | {elapsed:.1f}s")

    # ── Concatenate and compute global spatial extent ──
    print(f"  Concatenating per-gene coordinates...")
    t1 = time.time()

    gene_data = []
    global_x_min = np.inf
    global_x_max = -np.inf
    global_y_min = np.inf
    global_y_max = -np.inf

    for list_idx, gene_idx in enumerate(real_gene_indices):
        name = gene_names[gene_idx]

        if len(gene_x_coords[list_idx]) == 0:
            gene_data.append((name, np.array([], dtype=np.float32),
                              np.array([], dtype=np.float32)))
            continue

        x = np.concatenate(gene_x_coords[list_idx]).astype(np.float32)
        y = np.concatenate(gene_y_coords[list_idx]).astype(np.float32)

        if len(x) > 0:
            global_x_min = min(global_x_min, x.min())
            global_x_max = max(global_x_max, x.max())
            global_y_min = min(global_y_min, y.min())
            global_y_max = max(global_y_max, y.max())

        gene_data.append((name, x, y))

    del gene_x_coords, gene_y_coords

    print(f"  Spatial extent: X=[{global_x_min:.1f}, {global_x_max:.1f}] um, "
          f"Y=[{global_y_min:.1f}, {global_y_max:.1f}] um")
    print(f"  Concatenation took {time.time() - t1:.1f}s")

    # ── Quantize and write per-gene JSON files ──
    os.makedirs(output_dir, exist_ok=True)

    x_offset = float(global_x_min)
    y_offset = float(global_y_min)
    x_scale = QUANT_RESOLUTION
    y_scale = QUANT_RESOLUTION

    # Verify uint16 range is sufficient
    x_range_quant = (global_x_max - global_x_min) / QUANT_RESOLUTION
    y_range_quant = (global_y_max - global_y_min) / QUANT_RESOLUTION
    assert x_range_quant < 65535, f"X range too large for uint16: {x_range_quant}"
    assert y_range_quant < 65535, f"Y range too large for uint16: {y_range_quant}"
    print(f"  Quantized range: X=0-{int(x_range_quant)}, Y=0-{int(y_range_quant)} "
          f"(uint16 max=65535, OK)")

    print(f"  Writing per-gene JSON files to {output_dir}/")
    t2 = time.time()

    gene_index = []
    total_bytes = 0
    genes_written = 0

    for name, x, y in gene_data:
        n = len(x)
        if n == 0:
            continue

        x_quant = np.round((x - x_offset) / x_scale).astype(np.uint16)
        y_quant = np.round((y - y_offset) / y_scale).astype(np.uint16)

        gene_json = {
            "gene": name,
            "n": int(n),
            "x": x_quant.tolist(),
            "y": y_quant.tolist(),
        }

        fpath = os.path.join(output_dir, f"{name}.json")
        with open(fpath, 'w') as f:
            json.dump(gene_json, f, separators=(',', ':'))

        fsize = os.path.getsize(fpath)
        total_bytes += fsize
        genes_written += 1

        gene_index.append({
            "gene": name,
            "n": int(n),
            "file": f"{name}.json",
            "size_kb": round(fsize / 1024, 1),
        })

    gene_index.sort(key=lambda g: g['gene'])

    index_json = {
        "sample_id": sample_id,
        "n_genes": genes_written,
        "total_transcripts": int(total_real),
        "x_offset": x_offset,
        "y_offset": y_offset,
        "x_scale": x_scale,
        "y_scale": y_scale,
        "genes": gene_index,
    }

    with open(index_path, 'w') as f:
        json.dump(index_json, f, indent=2)

    # Close zip store if applicable
    if is_zip:
        store.close()

    elapsed = time.time() - t0
    print(f"\n  Done! {genes_written} gene files written in {elapsed:.1f}s")
    print(f"  Total file size: {total_bytes / 1024 / 1024:.1f} MB")
    print(f"  Gene index: {index_path}")

    gene_index.sort(key=lambda g: -g['n'])
    print(f"\n  Top 10 most abundant genes:")
    for g in gene_index[:10]:
        print(f"    {g['gene']:20s} {g['n']:>10,} molecules  {g['size_kb']:>8.1f} KB")


def main():
    # Discover all available zarr files
    all_samples = discover_zarr_files()
    print(f"Discovered {len(all_samples)} transcript zarr file(s):")
    for sid, path in sorted(all_samples.items()):
        print(f"  {sid}: {path}")

    # If command-line args given, filter to those samples
    if len(sys.argv) > 1:
        requested = sys.argv[1:]
        samples_to_process = {}
        for sid in requested:
            if sid in all_samples:
                samples_to_process[sid] = all_samples[sid]
            else:
                print(f"WARNING: No zarr file found for {sid}")
        if not samples_to_process:
            print("ERROR: No valid samples to process.")
            sys.exit(1)
    else:
        samples_to_process = all_samples

    print(f"\nWill process {len(samples_to_process)} sample(s): {', '.join(sorted(samples_to_process.keys()))}")

    os.makedirs(TRANSCRIPTS_BASE_DIR, exist_ok=True)

    for sample_id in sorted(samples_to_process.keys()):
        process_sample(sample_id, samples_to_process[sample_id])

    print(f"\n{'='*70}")
    print(f"All done! Transcript data exported to {TRANSCRIPTS_BASE_DIR}/")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
