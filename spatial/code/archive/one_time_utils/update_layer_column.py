#!/usr/bin/env python3
"""
Quick re-annotation: update the 'layer' column in all h5ad files
to use depth bins + Vascular OOD override (no Extra-cortical layer).

Does NOT re-run spatial domain clustering — just rewrites the layer column
using existing spatial_domain and predicted_norm_depth columns.
"""

import os
import sys
import time
import numpy as np
import anndata as ad
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "modules"))
from depth_model import assign_discrete_layers

H5AD_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "output", "h5ad")
H5AD_DIR = os.path.abspath(H5AD_DIR)

fnames = sorted(f for f in os.listdir(H5AD_DIR) if f.endswith("_annotated.h5ad"))
print(f"Found {len(fnames)} h5ad files")
print(f"Updating layer column: depth bins + Vascular OOD override")
print("=" * 60)

t0 = time.time()
for fname in fnames:
    sample_id = fname.replace("_annotated.h5ad", "")
    path = os.path.join(H5AD_DIR, fname)

    adata = ad.read_h5ad(path)

    # Get QC mask
    if "qc_pass" in adata.obs.columns:
        qc_mask = adata.obs["qc_pass"].values.astype(bool)
    else:
        qc_mask = np.ones(adata.n_obs, dtype=bool)

    # Assign depth-based layers for all cells
    depth = adata.obs["predicted_norm_depth"].values
    depth_layers = assign_discrete_layers(depth)

    # Build new hybrid layer: depth bins everywhere, only Vascular from OOD
    new_layer = np.full(adata.n_obs, "Unassigned", dtype=object)
    new_layer[qc_mask] = depth_layers[qc_mask]

    # Override Vascular from spatial_domain
    if "spatial_domain" in adata.obs.columns:
        vasc_mask = adata.obs["spatial_domain"].values == "Vascular"
        new_layer[vasc_mask] = "Vascular"

    # Non-QC cells stay Unassigned
    new_layer[~qc_mask] = "Unassigned"

    # Count changes
    old_layer = adata.obs["layer"].values.copy()
    changed = (old_layer != new_layer).sum()
    old_ec = (old_layer == "Extra-cortical").sum()

    adata.obs["layer"] = new_layer
    adata.write_h5ad(path)

    # Layer distribution
    c = Counter(new_layer[qc_mask])
    layer_str = "  ".join(
        f"{l}:{c.get(l, 0):,}"
        for l in ["L1", "L2/3", "L4", "L5", "L6", "WM", "Vascular"]
    )

    print(
        f"{sample_id}: {changed:,} changed (was {old_ec:,} Extra-cortical) | {layer_str}"
    )

elapsed = time.time() - t0
print(f"\nDone in {elapsed:.0f}s")
