#!/usr/bin/env python3
"""
Export MERSCOPE transcript data to viewer-compatible JSON format.

Reads MERSCOPE barcodes CSV files (decoded transcript positions) and converts
them to the compact per-gene JSON format used by the spatial viewer.

Output format matches the Xenium viewer transcript format:
  - transcripts/{sample_id}/gene_index.json
  - transcripts/{sample_id}/{GENE}.json (per gene)

Usage:
    python3 merscope_export_transcripts.py
"""

import os
import sys
import json
import time
import numpy as np
import pandas as pd

# Paths
BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
TRANSCRIPT_DIR = os.path.join(BASE_DIR, "data", "merscope_4k_probe_testing", "merscope_transcripts")
VIEWER_DIR = os.path.join(BASE_DIR, "output", "merscope_viewer", "transcripts")

# Map from CSV filename prefix to our sample ID format
SAMPLES = [
    {
        "csv": os.path.join(TRANSCRIPT_DIR, "H18.06.006.MTG.4000.expand.rep1.barcodes.csv.gz"),
        "sample_id": "H18.06.006_MTG_4000_rep1",
    },
    {
        "csv": os.path.join(TRANSCRIPT_DIR, "H18.06.006.MTG.250.expand.rep1.barcodes.csv.gz"),
        "sample_id": "H18.06.006_MTG_250_rep1",
    },
]


def process_sample(csv_path, sample_id):
    """Process a single MERSCOPE transcript CSV into viewer JSON files."""
    t0 = time.time()
    out_dir = os.path.join(VIEWER_DIR, sample_id)
    os.makedirs(out_dir, exist_ok=True)

    print(f"\nProcessing {sample_id}...")
    print(f"  Reading {csv_path}...")

    # Read CSV — only columns we need
    df = pd.read_csv(csv_path, usecols=["gene_name", "global_x", "global_y"])
    print(f"  {len(df):,} transcripts, {df['gene_name'].nunique()} genes")

    # Coordinate ranges
    x_min, x_max = df["global_x"].min(), df["global_x"].max()
    y_min, y_max = df["global_y"].min(), df["global_y"].max()
    print(f"  x: [{x_min:.1f}, {x_max:.1f}], y: [{y_min:.1f}, {y_max:.1f}]")

    # Choose quantization scale (0.2 µm per step, matching Xenium viewer)
    scale = 0.2
    x_offset = float(x_min)
    y_offset = float(y_min)

    # Quantize coordinates
    df["qx"] = np.round((df["global_x"].values - x_offset) / scale).astype(np.int32)
    df["qy"] = np.round((df["global_y"].values - y_offset) / scale).astype(np.int32)

    # Group by gene and write per-gene JSON files
    gene_groups = df.groupby("gene_name")
    gene_index = []
    n_written = 0

    for gene, group in gene_groups:
        n = len(group)
        gene_file = f"{gene}.json"
        gene_data = {
            "gene": gene,
            "n": n,
            "x": group["qx"].tolist(),
            "y": group["qy"].tolist(),
        }
        gene_path = os.path.join(out_dir, gene_file)
        with open(gene_path, "w") as f:
            json.dump(gene_data, f, separators=(",", ":"))

        file_size_kb = os.path.getsize(gene_path) / 1024
        gene_index.append({
            "gene": gene,
            "n": n,
            "file": gene_file,
            "size_kb": round(file_size_kb, 1),
        })
        n_written += 1
        if n_written % 500 == 0:
            print(f"  ... {n_written}/{len(gene_groups)} genes")

    # Sort gene index by name
    gene_index.sort(key=lambda g: g["gene"])

    # Write gene_index.json
    index_data = {
        "sample_id": sample_id,
        "n_genes": len(gene_index),
        "total_transcripts": len(df),
        "x_offset": x_offset,
        "y_offset": y_offset,
        "x_scale": scale,
        "y_scale": scale,
        "genes": gene_index,
    }
    index_path = os.path.join(out_dir, "gene_index.json")
    with open(index_path, "w") as f:
        json.dump(index_data, f, indent=2)

    total_size = sum(g["size_kb"] for g in gene_index) / 1024
    elapsed = time.time() - t0
    print(f"  Done: {n_written} genes, {len(df):,} transcripts")
    print(f"  Total size: {total_size:.1f} MB, Time: {elapsed:.1f}s")
    print(f"  Output: {out_dir}/")


def main():
    os.makedirs(VIEWER_DIR, exist_ok=True)

    for sample in SAMPLES:
        csv_path = sample["csv"]
        if not os.path.exists(csv_path):
            print(f"Skipping {sample['sample_id']}: {csv_path} not found")
            continue
        process_sample(csv_path, sample["sample_id"])


if __name__ == "__main__":
    main()
