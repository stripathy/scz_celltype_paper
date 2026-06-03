#!/usr/bin/env python3
"""
Generate a gallery of representative curved cortex strips for both
Xenium and MERFISH data.

Creates two multi-panel figures:
  1. Xenium gallery (6 representative samples)
  2. MERFISH gallery (6 representative sections)

Each panel shows cells in selected strips colored by cortical layer,
with pia curves overlaid in white.
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import H5AD_DIR, SAMPLE_TO_DX, LAYER_COLORS, LAYER_ORDER

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "modules"))

from banksy_05_curved_strips import (
    extract_l1_points,
    fit_pia_curve,
    compute_normals,
    split_pia_at_folds,
    assign_cells_to_strips,
    score_strips,
    MIN_SEGMENT_L1_CELLS,
    STRIP_WIDTH_ALONG_PIA,
    STRIP_MAX_DEPTH_UM,
    PIA_SAMPLE_SPACING,
    OUT_DIR,
)
from scipy.interpolate import UnivariateSpline
from banksy_05c_merfish_strips import (
    prepare_merfish_section,
    MERFISH_PATH,
    MERFISH_OUT_DIR,
)

# ── Layer colors (same as pipeline, works for both Xenium and MERFISH) ──
# Extend with WM for display
DISPLAY_LAYER_COLORS = dict(LAYER_COLORS)
DISPLAY_LAYER_COLORS.setdefault("WM", "#888888")

DISPLAY_LAYER_ORDER = ["L1", "L2/3", "L4", "L5", "L6", "WM"]


def run_strip_pipeline(adata):
    """Run the full strip pipeline on an adata and return results.

    Expects adata to already have: banksy_is_l1, layer, banksy_domain,
    predicted_norm_depth, and spatial in obsm.
    """
    domains = adata.obs["banksy_domain"].values.astype(str)
    pred_depth = pd.to_numeric(adata.obs["predicted_norm_depth"], errors="coerce").values

    cortical_mask = domains == "Cortical"
    cortical_coords = adata.obsm["spatial"][cortical_mask]
    cortical_depths = pred_depth[cortical_mask]

    l1_clean, l1_segment_labels, l1_clean_idx = extract_l1_points(adata)
    if len(l1_clean) == 0:
        return None

    segments = np.unique(l1_segment_labels)
    pia_segments = []
    bank_counter = 0

    for seg_id in segments:
        seg_mask = l1_segment_labels == seg_id
        n_seg = seg_mask.sum()
        if n_seg < MIN_SEGMENT_L1_CELLS:
            continue
        seg_coords = l1_clean[seg_mask]
        spline_x, spline_y, t_range, method = fit_pia_curve(seg_coords, seg_id)
        if spline_x is None:
            continue
        n_sample = max(50, int((t_range[1] - t_range[0]) / 50))
        sample_t, sample_xy, normals, tangents = compute_normals(
            spline_x, spline_y, t_range, cortical_coords, cortical_depths,
            n_sample=n_sample
        )
        fitted_seg = {
            "segment_id": seg_id, "spline_x": spline_x, "spline_y": spline_y,
            "t_range": t_range, "sample_t": sample_t, "sample_xy": sample_xy,
            "normals": normals, "tangents": tangents, "method": method,
        }
        banks = split_pia_at_folds(fitted_seg)
        for b in banks:
            b["bank_id"] = bank_counter
            bank_counter += 1
        pia_segments.extend(banks)

    if not pia_segments:
        return None

    strip_width = STRIP_WIDTH_ALONG_PIA
    strip_ids, strip_seg, strip_bank, perp_depth, pia_arc = assign_cells_to_strips(
        adata.obsm["spatial"], domains, pia_segments, strip_width=strip_width
    )
    strip_scores, complete_ids, partial_ids = score_strips(strip_ids, pred_depth, domains)

    n_cortical = cortical_mask.sum()
    n_selected = sum((strip_ids == s).sum() for s in (complete_ids | partial_ids))
    coverage = n_selected / max(1, n_cortical) * 100

    # Compute strip boundary geometry for drawing
    # For each bank, compute the pia points + normals at strip boundaries
    strip_boundaries = []  # list of (global_strip_id, pia_xy, normal, bank_id)
    global_strip_offset = 0

    for seg_i, seg in enumerate(pia_segments):
        if seg["sample_xy"] is None or len(seg["sample_xy"]) < 2:
            continue

        t_min, t_max = seg["t_range"]
        arc_length = t_max - t_min
        n_strips_seg = max(1, int(np.ceil(arc_length / strip_width)))
        bank_id = seg.get("bank_id", seg_i)

        spline_x = seg["spline_x"]
        spline_y = seg["spline_y"]
        sample_t = seg["sample_t"]
        normals = seg["normals"]

        # For each strip boundary (including start and end of bank)
        for si in range(n_strips_seg + 1):
            t_boundary = t_min + si * strip_width
            t_boundary = min(t_boundary, t_max)

            # Get pia point at this boundary
            pia_pt = np.array([float(spline_x(t_boundary)),
                               float(spline_y(t_boundary))])

            # Interpolate normal at this boundary
            ni = np.searchsorted(sample_t, t_boundary)
            ni = np.clip(ni, 0, len(normals) - 1)
            normal = normals[ni]

            # Which global strip IDs are on each side of this boundary
            strip_left = global_strip_offset + si - 1 if si > 0 else -1
            strip_right = global_strip_offset + si if si < n_strips_seg else -1

            strip_boundaries.append({
                "pia_pt": pia_pt,
                "normal": normal,
                "strip_left": strip_left,
                "strip_right": strip_right,
                "bank_id": bank_id,
            })

        global_strip_offset += n_strips_seg

    return {
        "pia_segments": pia_segments,
        "strip_ids": strip_ids,
        "strip_bank": strip_bank,
        "complete_ids": complete_ids,
        "partial_ids": partial_ids,
        "n_cortical": n_cortical,
        "n_selected": n_selected,
        "coverage": coverage,
        "n_banks": len(set(strip_bank[strip_bank >= 0])),
        "n_complete": len(complete_ids),
        "n_partial": len(partial_ids),
        "strip_boundaries": strip_boundaries,
    }


def plot_strip_panel(ax, adata, result, title, show_all_cells=True):
    """Plot one panel: selected strips colored by layer + pia curves +
    complete strip boundaries drawn as perpendicular lines."""
    coords = adata.obsm["spatial"]
    x, y = coords[:, 0], coords[:, 1]
    layers = adata.obs["layer"].values.astype(str)
    strip_ids = result["strip_ids"]
    complete_ids = result["complete_ids"]
    partial_ids = result["partial_ids"]
    selected = complete_ids | partial_ids
    in_selected = np.array([s in selected for s in strip_ids])

    # Background: non-selected cells
    not_selected = ~in_selected
    if show_all_cells and not_selected.sum() > 0:
        ax.scatter(x[not_selected], y[not_selected], c="#1a1a1a", s=0.05,
                   alpha=0.15, rasterized=True)

    # Selected cells colored by layer
    for layer_name in DISPLAY_LAYER_ORDER:
        if layer_name not in DISPLAY_LAYER_COLORS:
            continue
        mask = in_selected & (layers == layer_name)
        if mask.sum() > 0:
            ax.scatter(x[mask], y[mask], c=[DISPLAY_LAYER_COLORS[layer_name]],
                       s=0.4, alpha=0.7, rasterized=True)

    # Pia curves in white
    for seg in result["pia_segments"]:
        if seg["sample_xy"] is not None:
            sx = seg["sample_xy"][:, 0]
            sy = seg["sample_xy"][:, 1]
            ax.plot(sx, sy, color="white", linewidth=1.8, zorder=5, alpha=0.9)

    # Draw complete strip boundaries
    # For each boundary that borders a complete strip, draw a perpendicular line
    # from the pia along the normal direction
    strip_boundaries = result.get("strip_boundaries", [])
    depth_line_len = 2000  # μm — how far to draw the boundary line from pia

    # Estimate actual max depth from the data for line length
    perp_depths = result.get("perp_depth", None)
    if perp_depths is None:
        # Use cells in complete strips to estimate depth extent
        in_complete = np.array([s in complete_ids for s in strip_ids])
        if in_complete.any():
            # rough depth estimate: distance from pia centroid
            pass

    drawn_boundaries = set()  # avoid double-drawing shared boundaries
    for bd in strip_boundaries:
        left_id = bd["strip_left"]
        right_id = bd["strip_right"]

        # Draw boundary if it borders a complete strip
        left_complete = left_id in complete_ids
        right_complete = right_id in complete_ids

        if not (left_complete or right_complete):
            continue

        # Create a hashable key to avoid duplicates
        bd_key = (round(bd["pia_pt"][0], 1), round(bd["pia_pt"][1], 1))
        if bd_key in drawn_boundaries:
            continue
        drawn_boundaries.add(bd_key)

        pia_pt = bd["pia_pt"]
        normal = bd["normal"]

        # Draw line from pia point along normal
        end_pt = pia_pt + normal * depth_line_len

        ax.plot([pia_pt[0], end_pt[0]], [pia_pt[1], end_pt[1]],
                color="#FFD700", linewidth=2.5, alpha=0.85, zorder=6,
                linestyle="-")

    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_facecolor("#0a0a0a")
    ax.set_xticks([])
    ax.set_yticks([])

    # Title with stats
    cov = result["coverage"]
    nc = result["n_complete"]
    np_ = result["n_partial"]
    nb = result["n_banks"]
    ax.set_title(f"{title}\n{nb} bank{'s' if nb!=1 else ''}, "
                 f"{nc}C+{np_}P strips, {cov:.0f}% coverage",
                 fontsize=12, fontweight="bold", color="white",
                 pad=6, backgroundcolor="#333333",
                 fontfamily="monospace")


def make_gallery(panels, suptitle, out_path, ncols=3):
    """Create a multi-panel gallery figure.

    panels: list of (adata, result, title)
    """
    nrows = (len(panels) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(8 * ncols, 7 * nrows),
                              facecolor="#111111")
    if nrows == 1:
        axes = axes.reshape(1, -1)

    for i, (adata, result, title) in enumerate(panels):
        row, col = divmod(i, ncols)
        ax = axes[row, col]
        plot_strip_panel(ax, adata, result, title)

    # Turn off unused axes
    for i in range(len(panels), nrows * ncols):
        row, col = divmod(i, ncols)
        axes[row, col].set_visible(False)

    # Shared legend
    legend_elements = [
        Patch(facecolor=DISPLAY_LAYER_COLORS.get(l, "#888"), label=l)
        for l in DISPLAY_LAYER_ORDER if l in DISPLAY_LAYER_COLORS
    ]
    legend_elements.append(
        plt.Line2D([0], [0], color="white", linewidth=2, label="Pia curve")
    )
    legend_elements.append(
        plt.Line2D([0], [0], color="#FFD700", linewidth=2.5, alpha=0.85,
                   linestyle="-", label="Strip boundary")
    )
    fig.legend(handles=legend_elements, loc="lower center", ncol=len(legend_elements),
               fontsize=13, framealpha=0.9, facecolor="#222222", labelcolor="white",
               edgecolor="#555555", handlelength=1.5, handleheight=1.2)

    fig.suptitle(suptitle, fontsize=22, fontweight="bold", color="white", y=0.98)
    fig.subplots_adjust(bottom=0.06, top=0.93, hspace=0.22, wspace=0.08)

    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="#111111")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ── Main ──────────────────────────────────────────────────────────────

def main():
    t0 = time.time()

    # ── XENIUM GALLERY ──
    print("=" * 60)
    print("XENIUM GALLERY")
    print("=" * 60)

    xenium_samples = ["Br2039", "Br5622", "Br2421", "Br5588", "Br5400", "Br1113"]
    xenium_panels = []

    for sample_id in xenium_samples:
        print(f"\nProcessing {sample_id}...")
        path = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
        adata = ad.read_h5ad(path)

        # QC filter
        if "hybrid_qc_pass" in adata.obs.columns:
            adata = adata[adata.obs["hybrid_qc_pass"].values.astype(bool)].copy()
        elif "qc_pass" in adata.obs.columns:
            adata = adata[adata.obs["qc_pass"].values.astype(bool)].copy()

        result = run_strip_pipeline(adata)
        if result is None:
            print(f"  SKIP: no strips found")
            continue

        dx = SAMPLE_TO_DX.get(sample_id, "?")
        title = f"{sample_id} ({dx})"
        xenium_panels.append((adata, result, title))
        print(f"  {result['n_banks']} banks, {result['n_complete']}C+{result['n_partial']}P, "
              f"{result['coverage']:.1f}%")

    if xenium_panels:
        out_path = os.path.join(OUT_DIR, "gallery_xenium_strips.png")
        make_gallery(xenium_panels, "Xenium SCZ — Curved Cortex Strips (Representative Samples)",
                     out_path, ncols=3)

    # ── MERFISH GALLERY ──
    print("\n" + "=" * 60)
    print("MERFISH GALLERY")
    print("=" * 60)

    print("Loading MERFISH data...")
    adata_full = ad.read_h5ad(MERFISH_PATH)
    print(f"  {adata_full.n_obs:,} cells loaded")

    # Pick representative sections: varied donors, high coverage, different sizes
    merfish_sections = [
        "H20.33.001.CX28.MTG.02.007.1.02.03",   # 100%, 3 banks
        "H20.33.015.CX24.MTG.02.007.1.03.01",   # 100%, 3 banks, large
        "H20.33.040.Cx25.MTG.02.007.1.01.04",   # 100%, 7+5 strips
        "H21.33.040.Cx22.MTG.02.007.3.03.01",   # 100%, 4 complete
        "H21.33.031.CX24.MTG.02.007.1.01.01",   # 93.8%, biggest section (128K)
        "H21.33.025.CX26.MTG.02.007.4.01.06",   # 100%, 3 banks
    ]

    merfish_panels = []

    for section_id in merfish_sections:
        print(f"\nProcessing {section_id}...")
        sec_mask = adata_full.obs["Section"] == section_id
        adata_sec = adata_full[sec_mask]

        adata_prep = prepare_merfish_section(adata_sec)
        result = run_strip_pipeline(adata_prep)

        if result is None:
            print(f"  SKIP: no strips found")
            continue

        # Short title: donor + section count
        donor = section_id.split(".")[0] + "." + section_id.split(".")[1] + "." + section_id.split(".")[2]
        n_cells = adata_prep.n_obs
        title = f"{donor} ({n_cells//1000}K cells)"
        merfish_panels.append((adata_prep, result, title))
        print(f"  {result['n_banks']} banks, {result['n_complete']}C+{result['n_partial']}P, "
              f"{result['coverage']:.1f}%")

    if merfish_panels:
        out_path = os.path.join(MERFISH_OUT_DIR, "gallery_merfish_strips.png")
        make_gallery(merfish_panels, "SEA-AD MERFISH — Curved Cortex Strips (Representative Sections)",
                     out_path, ncols=3)

    print(f"\nDone in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
