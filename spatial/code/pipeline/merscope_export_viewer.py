#!/usr/bin/env python3
"""
Export MERSCOPE viewer JSON files from annotated h5ad files.

Creates a separate viewer instance (output/merscope_viewer/) for the
MERSCOPE Fang et al. data, independent of the main Xenium viewer.

Adapts 07_export_viewer.py for MERSCOPE data:
  - Reads from output/merscope_h5ad/ instead of output/h5ad/
  - Sample metadata from MERSCOPE filenames (donor, region, panel size)
  - No diagnosis column (these are non-disease reference samples)
  - Prepares transcript directory stubs for later population

Output:
  - output/merscope_viewer/<sample_id>.json  (one per sample)
  - output/merscope_viewer/index.json        (global index + color palettes)
  - output/merscope_viewer/index.html        (copied from Xenium viewer)
  - output/merscope_viewer/transcripts/      (empty, for future transcript data)

Usage:
    python3 -u merscope_export_viewer.py
"""

import os
import sys
import json
import time
import shutil
import numpy as np
import anndata as ad

# Pipeline config
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_config import BASE_DIR, OUTPUT_DIR, MODULES_DIR, VIEWER_DIR

# Modules
sys.path.insert(0, MODULES_DIR)
from constants import SUBCLASS_TO_CLASS

# ── Paths ──
MERSCOPE_H5AD_DIR = os.path.join(OUTPUT_DIR, "merscope_h5ad")
MERSCOPE_VIEWER_DIR = os.path.join(OUTPUT_DIR, "merscope_viewer")

# ── Layer categories ──
LAYER_CATS = ["L1", "L2/3", "L4", "L5", "L6", "WM", "Vascular"]

# ── Color palettes (same as Xenium viewer) ──
SUBCLASS_COLORS = {
    "L2/3 IT": "#B1EC30", "L4 IT": "#00E5E5", "L5 IT": "#50B2AD",
    "L5 ET": "#0D5B78", "L5/6 NP": "#3E9E64", "L6 IT": "#A19922",
    "L6 IT Car3": "#5100FF", "L6 CT": "#2D8CB8", "L6b": "#7044AA",
    "Sst": "#FF9900", "Sst Chodl": "#B1B10C", "Pvalb": "#D93137",
    "Vip": "#A45FBF", "Lamp5": "#DA808C", "Lamp5 Lhx6": "#935F50",
    "Sncg": "#DF70FF", "Pax6": "#71238C", "Chandelier": "#F641A8",
    "Astrocyte": "#665C47", "Oligodendrocyte": "#53776C",
    "OPC": "#374A45", "Microglia-PVM": "#94AF97",
    "Endothelial": "#8D6C62", "VLMC": "#697255", "SMC": "#5A6349",
    "Pericyte": "#7A8267",
    "QC Failed": "#ffffff",
}

CLASS_COLORS = {
    "Glutamatergic": "#00ADF8", "GABAergic": "#F05A28",
    "Non-neuronal": "#808080", "QC Failed": "#ffffff",
}

LAYER_COLORS = {
    "L1": "#E54C4C", "L2/3": "#4CCC4C", "L4": "#4C4CE5",
    "L5": "#E5E533", "L6": "#E57F19", "WM": "#999999",
    "Vascular": "#F24C99",
}

SUPERTYPE_SUBCLASS_MAP = {
    "Astro": "Astrocyte", "Oligo": "Oligodendrocyte",
    "Micro-PVM": "Microglia-PVM", "Endo": "Endothelial",
    "Pericyte": "VLMC", "SMC": "VLMC",
}

# HANN class labels → standard names
HANN_CLASS_TO_CLASS = {
    'Neuronal: Glutamatergic': 'Glutamatergic',
    'Neuronal: GABAergic': 'GABAergic',
    'Non-neuronal and Non-neural': 'Non-neuronal',
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


def _parse_sample_metadata(sample_id):
    """Extract metadata from MERSCOPE sample ID.

    Format: H18.06.006_MTG_4000_rep1
    """
    parts = sample_id.split("_")
    donor = parts[0] if len(parts) >= 1 else sample_id
    region = parts[1] if len(parts) >= 2 else "Unknown"
    panel = parts[2] if len(parts) >= 3 else "Unknown"
    rep = parts[3] if len(parts) >= 4 else ""
    return {
        "donor": donor,
        "region": region,
        "panel_size": panel,
        "replicate": rep,
        "label": f"{donor} {region} {panel}g {rep}",
    }


def main():
    t0 = time.time()

    os.makedirs(MERSCOPE_VIEWER_DIR, exist_ok=True)
    os.makedirs(os.path.join(MERSCOPE_VIEWER_DIR, "transcripts"), exist_ok=True)

    # Copy index.html from Xenium viewer
    src_html = os.path.join(VIEWER_DIR, "index.html")
    dst_html = os.path.join(MERSCOPE_VIEWER_DIR, "index.html")
    if os.path.exists(src_html):
        shutil.copy2(src_html, dst_html)
        print(f"Copied viewer HTML from {src_html}")
    else:
        print(f"WARNING: {src_html} not found — viewer HTML not copied")

    # Process each sample
    sample_index = []
    global_supertype_to_subclass = {}

    fnames = sorted(f for f in os.listdir(MERSCOPE_H5AD_DIR)
                    if f.endswith("_annotated.h5ad"))
    print(f"Found {len(fnames)} h5ad files in {MERSCOPE_H5AD_DIR}\n")

    for fname in fnames:
        sample_id = fname.replace("_annotated.h5ad", "")
        path = os.path.join(MERSCOPE_H5AD_DIR, fname)
        meta = _parse_sample_metadata(sample_id)

        print(f"Processing {sample_id}...", end=" ", flush=True)
        adata = ad.read_h5ad(path)
        n_total = adata.n_obs

        # Filter to QC-pass cells
        if "qc_pass" in adata.obs.columns:
            qc_mask = adata.obs["qc_pass"].values.astype(bool)
            n_fail = (~qc_mask).sum()
            adata = adata[qc_mask].copy()
            print(f"({n_total} total, {n_fail} QC-fail, {adata.n_obs} kept)",
                  end=" ")
        else:
            print(f"(all {n_total})", end=" ")

        # Extract spatial coordinates
        x = adata.obsm["spatial"][:, 0].astype(np.float32)
        y = adata.obsm["spatial"][:, 1].astype(np.float32)

        # Depth (may not exist for all samples)
        if "predicted_norm_depth" in adata.obs.columns:
            depth = adata.obs["predicted_norm_depth"].values.astype(np.float32)
        else:
            depth = np.zeros(len(x), dtype=np.float32)

        # Cell type labels: prefer correlation classifier, fall back to HANN
        has_corr = "corr_subclass" in adata.obs.columns
        if has_corr:
            subclass = adata.obs["corr_subclass"].values.astype(str)
            supertype = adata.obs["corr_supertype"].values.astype(str)
            class_label = adata.obs["corr_class"].values.astype(str)
            hann_subclass = adata.obs["subclass_label"].values.astype(str)
        elif "subclass_label" in adata.obs.columns:
            subclass = adata.obs["subclass_label"].values.astype(str)
            supertype = adata.obs["supertype_label"].values.astype(str)
            # Map HANN class names to standard names
            class_raw = adata.obs["class_label"].values.astype(str)
            class_label = np.array([HANN_CLASS_TO_CLASS.get(c, c) for c in class_raw])
            hann_subclass = subclass.copy()
        else:
            # No annotations yet
            subclass = np.full(len(x), "Unassigned", dtype=object)
            supertype = np.full(len(x), "Unassigned", dtype=object)
            class_label = np.full(len(x), "Unassigned", dtype=object)
            hann_subclass = subclass.copy()

        # HANN confidence scores
        conf_cols = {
            "class": "class_label_confidence",
            "subclass": "subclass_label_confidence",
            "supertype": "supertype_label_confidence",
        }
        conf_q = {}
        for level, col in conf_cols.items():
            if col in adata.obs.columns:
                vals = np.nan_to_num(
                    adata.obs[col].values.astype(np.float32), nan=0.0)
            else:
                vals = np.zeros(len(x), dtype=np.float32)
            conf_q[level] = np.clip(
                np.round(vals * 200), 0, 200).astype(np.uint8)

        # Track supertype → subclass mapping
        for sup_val, sub_val in zip(supertype, subclass):
            if sup_val not in global_supertype_to_subclass:
                global_supertype_to_subclass[sup_val] = sub_val

        # QC flagging for correlation classifier
        qc_col = None
        if "corr_qc_pass" in adata.obs.columns:
            qc_col = "corr_qc_pass"
        n_qc_fail = 0
        if qc_col is not None:
            qc_fail_mask = ~adata.obs[qc_col].values.astype(bool)
            n_qc_fail = qc_fail_mask.sum()
            if n_qc_fail > 0:
                subclass = subclass.copy()
                supertype = supertype.copy()
                class_label = class_label.copy()
                subclass[qc_fail_mask] = "QC Failed"
                supertype[qc_fail_mask] = "QC Failed"
                class_label[qc_fail_mask] = "QC Failed"
        global_supertype_to_subclass["QC Failed"] = "QC Failed"

        # Round coordinates
        x = np.round(x, 1)
        y = np.round(y, 1)
        depth = np.round(depth, 3)
        depth = np.nan_to_num(depth, nan=0.0)

        # Build compact integer-encoded categories
        subclass_cats = sorted(set(subclass))
        supertype_cats = sorted(set(supertype))
        class_cats = sorted(set(class_label))

        subclass_idx = {c: i for i, c in enumerate(subclass_cats)}
        supertype_idx = {c: i for i, c in enumerate(supertype_cats)}
        class_idx = {c: i for i, c in enumerate(class_cats)}

        # Layer column
        if "layer" in adata.obs.columns:
            layers = adata.obs["layer"].values.astype(str)
        else:
            from depth_model import assign_discrete_layers
            layers = assign_discrete_layers(depth)

        layer_cats = list(LAYER_CATS)
        layer_idx = {c: i for i, c in enumerate(layer_cats)}
        for lbl in set(layers):
            if lbl not in layer_idx:
                layer_cats.append(lbl)
                layer_idx[lbl] = len(layer_cats) - 1

        # HANN subclass for tooltip comparison
        hann_subclass_cats = sorted(set(hann_subclass))
        hann_subclass_idx = {c: i for i, c in enumerate(hann_subclass_cats)}

        data = {
            "sample_id": sample_id,
            "diagnosis": f"{meta['panel_size']}g {meta['region']}",
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
            "conf_class": conf_q["class"].tolist(),
            "conf_subclass": conf_q["subclass"].tolist(),
            "conf_supertype": conf_q["supertype"].tolist(),
            "hann_subclass_cats": hann_subclass_cats,
            "hann_subclass": [hann_subclass_idx[s] for s in hann_subclass],
        }

        # QC flag
        if qc_col is not None:
            qc_flag = adata.obs[qc_col].values.astype(bool)
            data["corr_qc"] = [1 if q else 0 for q in qc_flag]

        # Correlation margin
        if has_corr and "corr_subclass_margin" in adata.obs.columns:
            margin = np.nan_to_num(
                adata.obs["corr_subclass_margin"].values.astype(np.float32),
                nan=0.0)
            margin_q = np.clip(
                np.round(margin * 1000), 0, 255).astype(np.uint8)
            data["corr_margin"] = margin_q.tolist()

        # Write compact JSON
        json_path = os.path.join(MERSCOPE_VIEWER_DIR, f"{sample_id}.json")
        with open(json_path, "w") as f:
            json.dump(data, f, separators=(",", ":"))

        file_size = os.path.getsize(json_path) / 1024 / 1024
        print(f"-> {file_size:.1f}MB")

        sample_index.append({
            "sample_id": sample_id,
            "diagnosis": f"{meta['panel_size']}g {meta['region']}",
            "n_cells": len(x),
        })

    # Build supertype colors
    supertype_colors = _build_supertype_colors(
        global_supertype_to_subclass, SUBCLASS_COLORS)
    print(f"\nBuilt supertype colors: {len(supertype_colors)} supertypes")

    # Write index.json
    index_data = {
        "samples": sample_index,
        "subclass_colors": SUBCLASS_COLORS,
        "supertype_colors": supertype_colors,
        "class_colors": CLASS_COLORS,
        "layer_colors": LAYER_COLORS,
        "supertype_subclass_map": SUPERTYPE_SUBCLASS_MAP,
    }

    index_path = os.path.join(MERSCOPE_VIEWER_DIR, "index.json")
    with open(index_path, "w") as f:
        json.dump(index_data, f, indent=2)

    elapsed = time.time() - t0
    total_cells = sum(s["n_cells"] for s in sample_index)
    print(f"\nDone! {len(fnames)} samples, {total_cells:,} cells in {elapsed:.1f}s")
    print(f"Output: {MERSCOPE_VIEWER_DIR}/")
    print(f"\nTo view: cd {MERSCOPE_VIEWER_DIR} && python3 -m http.server 8001")


if __name__ == "__main__":
    main()
