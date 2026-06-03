#!/usr/bin/env python3
"""
Step 6: Export viewer JSON files from annotated h5ad files.

Reads per-sample h5ad files, filters to QC-pass cells only, and writes
compact JSON files for the interactive Xenium spatial viewer.

Output:
  - output/viewer/<sample_id>.json  (one per sample)
  - output/viewer/index.json        (global index with sample list + color palettes)
"""

import os
import sys
import json
import time
import numpy as np
import anndata as ad

# Pipeline config
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_config import H5AD_DIR, VIEWER_DIR, METADATA_PATH, MODULES_DIR

# Modules
sys.path.insert(0, MODULES_DIR)
from metadata import get_diagnosis_map

# ── Layer categories (depth bins + Vascular from BANKSY) ──
LAYER_CATS = ["L1", "L2/3", "L4", "L5", "L6", "WM", "Vascular"]


# ── Color palettes ──
# Official SEA-AD colors (from MERFISH .uns['Subclass_colors'])
SUBCLASS_COLORS = {
    "L2/3 IT": "#B1EC30",
    "L4 IT": "#00E5E5",
    "L5 IT": "#50B2AD",
    "L5 ET": "#0D5B78",
    "L5/6 NP": "#3E9E64",
    "L6 IT": "#A19922",
    "L6 IT Car3": "#5100FF",
    "L6 CT": "#2D8CB8",
    "L6b": "#7044AA",
    "Sst": "#FF9900",
    "Sst Chodl": "#B1B10C",
    "Pvalb": "#D93137",
    "Vip": "#A45FBF",
    "Lamp5": "#DA808C",
    "Lamp5 Lhx6": "#935F50",
    "Sncg": "#DF70FF",
    "Pax6": "#71238C",
    "Chandelier": "#F641A8",
    "Astrocyte": "#665C47",
    "Oligodendrocyte": "#53776C",
    "OPC": "#374A45",
    "Microglia-PVM": "#94AF97",
    "Endothelial": "#8D6C62",
    "VLMC": "#697255",
    "QC Failed": "#ffffff",
}

CLASS_COLORS = {
    "Glutamatergic": "#00ADF8",
    "GABAergic": "#F05A28",
    "Non-neuronal": "#808080",
    "QC Failed": "#ffffff",
}

LAYER_COLORS = {
    "L1": "#E54C4C",
    "L2/3": "#4CCC4C",
    "L4": "#4C4CE5",
    "L5": "#E5E533",
    "L6": "#E57F19",
    "WM": "#999999",
    "Vascular": "#F24C99",
}

# Supertype → Subclass mapping for non-neuronal types whose supertype names
# don't match their subclass names (e.g., "Astro_1" → "Astrocyte")
SUPERTYPE_SUBCLASS_MAP = {
    "Astro": "Astrocyte",
    "Oligo": "Oligodendrocyte",
    "Micro-PVM": "Microglia-PVM",
    "Endo": "Endothelial",
    "Pericyte": "VLMC",
    "SMC": "VLMC",
}


def _build_supertype_colors(supertype_to_subclass, subclass_colors):
    """Build per-supertype colors by varying brightness of parent subclass color."""
    from collections import defaultdict

    def hex_to_rgb(h):
        h = h.lstrip('#')
        return [int(h[i:i+2], 16) for i in (0, 2, 4)]

    def rgb_to_hex(r, g, b):
        return '#{:02X}{:02X}{:02X}'.format(
            max(0, min(255, r)), max(0, min(255, g)), max(0, min(255, b)))

    def adjust_brightness(hex_color, factor):
        r, g, b = hex_to_rgb(hex_color)
        return rgb_to_hex(int(r * factor), int(g * factor), int(b * factor))

    # Group supertypes by parent subclass
    subclass_to_supertypes = defaultdict(list)
    for sup, sub in sorted(supertype_to_subclass.items()):
        subclass_to_supertypes[sub].append(sup)

    supertype_colors = {}
    for sub, sups in sorted(subclass_to_supertypes.items()):
        base_color = subclass_colors.get(sub, '#666666')
        n = len(sups)
        if n == 1:
            supertype_colors[sups[0]] = base_color
        else:
            for i, sup in enumerate(sorted(sups)):
                factor = 0.7 + (i / max(1, n - 1)) * 0.6
                supertype_colors[sup] = adjust_brightness(base_color, factor)

    return supertype_colors


def main():
    t0 = time.time()

    os.makedirs(VIEWER_DIR, exist_ok=True)

    # Load diagnosis mapping
    dx_map = get_diagnosis_map(METADATA_PATH)
    print(f"Loaded diagnosis map: {len(dx_map)} subjects")

    # Process each sample
    sample_index = []
    total_cells_all = 0
    total_qc_pass = 0
    # Track supertype -> subclass mapping across all samples
    global_supertype_to_subclass = {}

    fnames = sorted(f for f in os.listdir(H5AD_DIR) if f.endswith("_annotated.h5ad"))
    print(f"Found {len(fnames)} h5ad files\n")

    for fname in fnames:
        sample_id = fname.replace("_annotated.h5ad", "")
        path = os.path.join(H5AD_DIR, fname)

        print(f"Processing {sample_id}...", end=" ", flush=True)
        adata = ad.read_h5ad(path)
        n_total = adata.n_obs

        # Keep ALL cells — QC status is tracked as a field, not used for filtering
        adata = adata.copy()
        print(f"({n_total} total)", end=" ")

        total_cells_all += n_total
        total_qc_pass += n_total  # all cells exported now

        # Extract data
        x = adata.obsm["spatial"][:, 0].astype(np.float32)
        y = adata.obsm["spatial"][:, 1].astype(np.float32)
        depth = adata.obs["predicted_norm_depth"].values.astype(np.float32)

        # Use correlation-derived labels as primary (fall back to HANN)
        has_corr = "corr_subclass" in adata.obs.columns
        if has_corr:
            subclass = adata.obs["corr_subclass"].values.astype(str)
            supertype = adata.obs["corr_supertype"].values.astype(str)
            class_label = adata.obs["corr_class"].values.astype(str)
            hann_subclass = adata.obs["subclass_label"].values.astype(str)
        else:
            subclass = adata.obs["subclass_label"].values.astype(str)
            supertype = adata.obs["supertype_label"].values.astype(str)
            class_label = adata.obs["class_label"].values.astype(str)
            hann_subclass = subclass  # same as primary

        # Mapping confidence scores (from HANN bootstrap voting)
        # Quantize to uint8 (0-200 range, divide by 200 to recover) to save space
        # while preserving 0.5% precision
        conf_class = np.nan_to_num(adata.obs["class_label_confidence"].values.astype(np.float32), nan=0.0)
        conf_subclass = np.nan_to_num(adata.obs["subclass_label_confidence"].values.astype(np.float32), nan=0.0)
        conf_supertype = np.nan_to_num(adata.obs["supertype_label_confidence"].values.astype(np.float32), nan=0.0)
        conf_class_q = np.clip(np.round(conf_class * 200), 0, 200).astype(np.uint8)
        conf_subclass_q = np.clip(np.round(conf_subclass * 200), 0, 200).astype(np.uint8)
        conf_supertype_q = np.clip(np.round(conf_supertype * 200), 0, 200).astype(np.uint8)

        # Track supertype -> subclass mapping (before QC override)
        for sup_val, sub_val in zip(supertype, subclass):
            if sup_val not in global_supertype_to_subclass:
                global_supertype_to_subclass[sup_val] = sub_val

        # Compute qc_status: 0=pass, 1=fail_spatial_qc, 2=fail_low_margin, 3=fail_doublet
        qc_status = np.zeros(len(adata), dtype=np.uint8)
        has_qc_pass = "qc_pass" in adata.obs.columns
        has_corr_qc = has_corr and "corr_qc_pass" in adata.obs.columns
        has_doublet = "doublet_suspect" in adata.obs.columns

        if has_qc_pass:
            fail_spatial = ~adata.obs["qc_pass"].values.astype(bool)
            qc_status[fail_spatial] = 1

        if has_doublet:
            is_doublet = adata.obs["doublet_suspect"].values.astype(bool)
            # Only mark as doublet if not already failed spatial QC
            qc_status[(qc_status == 0) & is_doublet] = 3

        if has_corr_qc and has_qc_pass:
            passes_spatial = adata.obs["qc_pass"].values.astype(bool)
            fails_corr = ~adata.obs["corr_qc_pass"].values.astype(bool)
            # Low margin = passes spatial QC but fails corr_qc_pass AND not already doublet
            low_margin = passes_spatial & fails_corr
            if has_doublet:
                low_margin = low_margin & ~is_doublet
            qc_status[(qc_status == 0) & low_margin] = 2

        n_qc_fail = (qc_status > 0).sum()
        if n_qc_fail > 0:
            subclass = subclass.copy()
            supertype = supertype.copy()
            class_label = class_label.copy()
            qc_fail_mask = qc_status > 0
            subclass[qc_fail_mask] = "QC Failed"
            supertype[qc_fail_mask] = "QC Failed"
            class_label[qc_fail_mask] = "QC Failed"
        global_supertype_to_subclass["QC Failed"] = "QC Failed"

        # Round coordinates to save space
        x = np.round(x, 1)
        y = np.round(y, 1)
        depth = np.round(depth, 3)
        # Replace NaN with 0 for JSON compatibility (NaN is not valid JSON)
        depth = np.nan_to_num(depth, nan=0.0)

        # Build compact format: encode categories as integers
        subclass_cats = sorted(set(subclass))
        supertype_cats = sorted(set(supertype))
        class_cats = sorted(set(class_label))

        subclass_idx = {c: i for i, c in enumerate(subclass_cats)}
        supertype_idx = {c: i for i, c in enumerate(supertype_cats)}
        class_idx = {c: i for i, c in enumerate(class_cats)}

        # Layer column (from step 05)
        if "layer" in adata.obs.columns:
            layers = adata.obs["layer"].values.astype(str)
        else:
            # Fallback: assign from depth only
            from depth_model import assign_discrete_layers
            layers = assign_discrete_layers(depth)

        # Ensure all layer values are in LAYER_CATS
        layer_cats = list(LAYER_CATS)  # local copy to avoid mutation
        layer_idx = {c: i for i, c in enumerate(layer_cats)}
        for lbl in set(layers):
            if lbl not in layer_idx:
                layer_cats.append(lbl)
                layer_idx[lbl] = len(layer_cats) - 1

        # HANN subclass for tooltip comparison (integer-encoded)
        hann_subclass_cats = sorted(set(hann_subclass))
        hann_subclass_idx = {c: i for i, c in enumerate(hann_subclass_cats)}

        data = {
            "sample_id": sample_id,
            "diagnosis": dx_map.get(sample_id, "Unknown"),
            "n_cells": len(x),
            "x_range": [float(x.min()), float(x.max())],
            "y_range": [float(y.min()), float(y.max())],
            "subclass_cats": subclass_cats,
            "supertype_cats": supertype_cats,
            "class_cats": class_cats,
            "x": x.tolist(),
            "y": y.tolist(),
            "subclass": [subclass_idx[s] for s in subclass],
            "supertype": [supertype_idx[c] for c in supertype],
            "class": [class_idx[c] for c in class_label],
            "depth": depth.tolist(),
            "layer_cats": layer_cats,
            "layer": [layer_idx[lbl] for lbl in layers],
            # HANN mapping confidence (uint8 quantized: divide by 200 to get 0-1)
            "conf_class": conf_class_q.tolist(),
            "conf_subclass": conf_subclass_q.tolist(),
            "conf_supertype": conf_supertype_q.tolist(),
            # HANN subclass for tooltip comparison
            "hann_subclass_cats": hann_subclass_cats,
            "hann_subclass": [hann_subclass_idx[s] for s in hann_subclass],
        }

        # QC status: 0=pass, 1=fail_spatial, 2=fail_low_margin, 3=fail_doublet
        data["qc_status"] = qc_status.tolist()

        # Correlation margin (always export if available)
        if has_corr and "corr_subclass_margin" in adata.obs.columns:
            margin = np.nan_to_num(adata.obs["corr_subclass_margin"].values.astype(np.float32), nan=0.0)
            # Quantize margin: multiply by 1000, clamp to 0-255
            margin_q = np.clip(np.round(margin * 1000), 0, 255).astype(np.uint8)
            data["corr_margin"] = margin_q.tolist()

        # Write compact JSON (no whitespace)
        json_path = os.path.join(VIEWER_DIR, f"{sample_id}.json")
        with open(json_path, "w") as f:
            json.dump(data, f, separators=(",", ":"))

        file_size = os.path.getsize(json_path) / 1024 / 1024
        print(f"-> {file_size:.1f}MB")

        sample_index.append({
            "sample_id": sample_id,
            "diagnosis": dx_map.get(sample_id, "Unknown"),
            "n_cells": len(x),
        })

    # Build supertype colors from the observed supertype -> subclass mapping
    supertype_colors = _build_supertype_colors(
        global_supertype_to_subclass, SUBCLASS_COLORS)
    print(f"\nBuilt supertype colors: {len(supertype_colors)} supertypes "
          f"from {len(set(global_supertype_to_subclass.values()))} subclasses")

    # Write index.json
    index_data = {
        "samples": sample_index,
        "subclass_colors": SUBCLASS_COLORS,
        "supertype_colors": supertype_colors,
        "class_colors": CLASS_COLORS,
        "layer_colors": LAYER_COLORS,
        "supertype_subclass_map": SUPERTYPE_SUBCLASS_MAP,
    }

    index_path = os.path.join(VIEWER_DIR, "index.json")
    with open(index_path, "w") as f:
        json.dump(index_data, f, indent=2)

    # Build standalone HTML viewer
    try:
        from bundle_viewer import bundle_standalone_html
        html_path = os.path.join(VIEWER_DIR, "xenium_viewer_standalone.html")
        bundle_standalone_html(VIEWER_DIR, html_path)
        print(f"\nStandalone HTML viewer: {html_path}")
    except ImportError:
        print("\nNote: bundle_viewer module not available, skipping HTML build")
    except Exception as e:
        print(f"\nHTML viewer build failed: {e}")

    elapsed = time.time() - t0
    print(f"\nDone! Processed {len(fnames)} samples in {elapsed:.1f}s")
    print(f"  Total cells across all samples: {total_cells_all:,}")
    print(f"  QC-pass cells exported: {total_qc_pass:,}")
    print(f"  QC-fail cells removed: {total_cells_all - total_qc_pass:,}")
    print(f"  Output: {VIEWER_DIR}/")


if __name__ == "__main__":
    main()
