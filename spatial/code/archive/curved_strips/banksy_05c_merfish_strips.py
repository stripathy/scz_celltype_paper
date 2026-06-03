#!/usr/bin/env python3
"""
Apply curved cortex strip identification to complete SEA-AD MERFISH data.

Adapts the Xenium curved strip pipeline (banksy_05_curved_strips.py) to work
with the MERFISH data format (different column names, no BANKSY annotations).

Key differences from Xenium:
  - L1 detection: uses predicted_layer == 'L1' (no banksy_is_l1)
  - Domain classification: derived from predicted_layer + Subclass
  - Sections processed independently (69 sections, 27 donors)
  - No QC column — all cells used

Usage:
    python3 -u code/analysis/banksy_05c_merfish_strips.py --section <section_name>
    python3 -u code/analysis/banksy_05c_merfish_strips.py --all
    python3 -u code/analysis/banksy_05c_merfish_strips.py --all --save
"""

import os
import sys
import time
import argparse
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Project imports — reuse the core functions from the Xenium pipeline
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from banksy_05_curved_strips import (
    extract_l1_points,
    fit_pia_curve,
    compute_normals,
    split_pia_at_folds,
    assign_cells_to_strips,
    score_strips,
    plot_curved_strips,
    MIN_SEGMENT_L1_CELLS,
    STRIP_WIDTH_ALONG_PIA,
    OUT_DIR,
)

# ── Constants ─────────────────────────────────────────────────────────

MERFISH_PATH = os.path.expanduser(
    "~/Github/SCZ_Xenium/data/reference/SEAAD_MTG_MERFISH.2024-12-11.h5ad"
)
MERFISH_OUT_DIR = os.path.join(OUT_DIR, "merfish_strips")
os.makedirs(MERFISH_OUT_DIR, exist_ok=True)

# Vascular cell types in MERFISH
MERFISH_VASCULAR_TYPES = {"Endothelial", "VLMC"}


# ── Helper: Prepare MERFISH section for Xenium pipeline ──────────────

def prepare_merfish_section(adata_section):
    """Map MERFISH columns to what the Xenium curved strip pipeline expects.

    Creates:
      - banksy_is_l1: bool (from predicted_layer == 'L1')
      - layer: str (from predicted_layer, mapping WM → L6 for compatibility)
      - banksy_domain: str (Cortical / Vascular / Extra-cortical)
      - predicted_norm_depth: float (already exists)

    Parameters
    ----------
    adata_section : AnnData — subset of MERFISH data for one section

    Returns
    -------
    adata : AnnData — annotated copy ready for pipeline functions
    """
    adata = adata_section.copy()

    # 1. L1 detection
    pred_layer = adata.obs["predicted_layer"].astype(str).values
    adata.obs["banksy_is_l1"] = pred_layer == "L1"

    # 2. Layer column — map predicted_layer values
    # predicted_layer has: L1, L2/3, L4, L5, L6, WM
    adata.obs["layer"] = pred_layer

    # 3. Domain classification
    subclass = adata.obs["Subclass"].astype(str).values
    domains = np.full(len(adata), "Cortical", dtype=object)

    # Vascular: Endothelial + VLMC
    is_vascular = np.isin(subclass, list(MERFISH_VASCULAR_TYPES))
    domains[is_vascular] = "Vascular"

    # Extra-cortical: WM cells that aren't vascular
    is_wm = pred_layer == "WM"
    domains[is_wm & ~is_vascular] = "Extra-cortical"

    adata.obs["banksy_domain"] = domains

    # 4. predicted_norm_depth — already exists
    # Ensure it's numeric
    adata.obs["predicted_norm_depth"] = pd.to_numeric(
        adata.obs["predicted_norm_depth"], errors="coerce"
    )

    return adata


# ── Process one MERFISH section ──────────────────────────────────────

def process_merfish_section(adata_section, section_id, strip_width=STRIP_WIDTH_ALONG_PIA):
    """Full curved cortex strip pipeline for one MERFISH section.

    Parameters
    ----------
    adata_section : AnnData — raw MERFISH section data
    section_id : str — section identifier
    strip_width : float — strip width along pia in μm
    """
    t0 = time.time()

    # Extract donor from section name
    donor_id = section_id.split(".")[0] + "." + section_id.split(".")[1] + "." + section_id.split(".")[2]

    print(f"\n{'='*70}")
    print(f"CURVED CORTEX STRIPS (MERFISH): {section_id}")
    print(f"  Donor: {donor_id}, {adata_section.n_obs:,} cells")
    print(f"{'='*70}")

    # Prepare columns for pipeline
    adata = prepare_merfish_section(adata_section)
    print(f"  {adata.n_obs:,} cells")

    # Count domains
    domains = adata.obs["banksy_domain"].values
    n_cortical = (domains == "Cortical").sum()
    n_vascular = (domains == "Vascular").sum()
    n_extra = (domains == "Extra-cortical").sum()
    print(f"  Domains: {n_cortical:,} Cortical, {n_vascular:,} Vascular, {n_extra:,} Extra-cortical")

    # Step 1: Extract and clean L1 points
    print("\n  Step 1: Extracting L1 boundary...")
    l1_clean, l1_segment_labels, l1_clean_idx = extract_l1_points(adata)

    if len(l1_clean) == 0:
        print("  ERROR: No valid L1 points found. Cannot fit pia curve.")
        return _empty_result(section_id, donor_id, n_cortical, time.time() - t0)

    # Step 2+3: Fit pia curve and compute normals per segment
    print("\n  Step 2+3: Fitting pia curves and computing normals...")
    pred_depth = pd.to_numeric(adata.obs["predicted_norm_depth"], errors="coerce").values

    cortical_mask = domains == "Cortical"
    cortical_coords = adata.obsm["spatial"][cortical_mask]
    cortical_depths = pred_depth[cortical_mask]

    segments = np.unique(l1_segment_labels)
    pia_segments = []
    global_bank_counter = 0

    for seg_id in segments:
        seg_mask = l1_segment_labels == seg_id
        n_seg = seg_mask.sum()

        if n_seg < MIN_SEGMENT_L1_CELLS:
            print(f"    Skipping segment {seg_id}: only {n_seg} L1 cells")
            continue

        seg_coords = l1_clean[seg_mask]
        print(f"\n    Segment {seg_id}: {n_seg} L1 cells")

        # Fit pia curve
        spline_x, spline_y, t_range, method = fit_pia_curve(seg_coords, seg_id)

        if spline_x is None:
            print(f"    WARNING: Failed to fit pia curve for segment {seg_id}")
            continue

        # Compute normals
        n_sample = max(50, int((t_range[1] - t_range[0]) / 50))
        sample_t, sample_xy, normals, tangents = compute_normals(
            spline_x, spline_y, t_range, cortical_coords, cortical_depths,
            n_sample=n_sample
        )

        fitted_seg = {
            "segment_id": seg_id,
            "spline_x": spline_x,
            "spline_y": spline_y,
            "t_range": t_range,
            "sample_t": sample_t,
            "sample_xy": sample_xy,
            "normals": normals,
            "tangents": tangents,
            "method": method,
        }

        # Split at folds (gyral crowns / sulcal bottoms)
        banks = split_pia_at_folds(fitted_seg)
        for b in banks:
            b["bank_id"] = global_bank_counter
            global_bank_counter += 1
        if len(banks) > 1:
            print(f"    Split into {len(banks)} banks at folds")
            for b in banks:
                bl = b["t_range"][1] - b["t_range"][0]
                print(f"      Bank {b['bank_id']}: arc-length {bl:.0f} μm, "
                      f"{len(b['sample_t'])} sample points")
        else:
            bl = banks[0]["t_range"][1] - banks[0]["t_range"][0]
            print(f"    Single bank (bank {banks[0]['bank_id']}): arc-length {bl:.0f} μm")
        pia_segments.extend(banks)

    if not pia_segments:
        print("  ERROR: No valid pia curves fitted. Cannot define strips.")
        return _empty_result(section_id, donor_id, n_cortical, time.time() - t0)

    # Step 4: Assign cells to strips
    print(f"\n  Step 4: Assigning cells to strips (width={strip_width}μm)...")
    strip_ids, strip_segment, strip_bank, perp_depth_um, pia_arc_pos = assign_cells_to_strips(
        adata.obsm["spatial"], domains, pia_segments, strip_width=strip_width
    )

    # Step 5: Score strips
    print("\n  Step 5: Scoring strips...")
    strip_scores, complete_ids, partial_ids = score_strips(strip_ids, pred_depth, domains)

    # Step 6: Diagnostic figure
    print("\n  Step 6: Generating diagnostic figure...")
    try:
        fig = plot_curved_strips(adata, section_id, pia_segments, strip_ids,
                                 strip_scores, complete_ids, partial_ids, perp_depth_um)
        # Shorten filename for MERFISH (section IDs are very long)
        short_name = section_id.replace(".", "_")[:40]
        fig_path = os.path.join(MERFISH_OUT_DIR, f"curved_strips_{short_name}.png")
        fig.savefig(fig_path, dpi=120, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"  Saved: {fig_path}")
    except Exception as e:
        print(f"  WARNING: Figure generation failed: {e}")
        plt.close("all")

    elapsed = time.time() - t0
    n_selected = sum((strip_ids == s).sum() for s in (complete_ids | partial_ids))
    coverage = n_selected / max(1, n_cortical) * 100

    n_banks = len(set(strip_bank[strip_bank >= 0]))
    print(f"\n  SUMMARY:")
    print(f"    Banks (distinct cortical columns): {n_banks}")
    print(f"    Pia segments: {len(pia_segments)}")
    print(f"    Complete strips: {len(complete_ids)}")
    print(f"    Partial strips: {len(partial_ids)}")
    print(f"    Cells in selected strips: {n_selected:,} / {n_cortical:,} cortical "
          f"({coverage:.1f}%)")
    print(f"    Time: {elapsed:.0f}s")

    return {
        "section_id": section_id,
        "donor_id": donor_id,
        "n_cells": adata.n_obs,
        "n_cortical": n_cortical,
        "n_banks": n_banks,
        "n_pia_segments": len(pia_segments),
        "n_complete_strips": len(complete_ids),
        "n_partial_strips": len(partial_ids),
        "n_cells_selected": n_selected,
        "coverage_pct": coverage,
        "time_sec": elapsed,
        "strip_ids": strip_ids,
        "strip_bank": strip_bank,
        "complete_ids": complete_ids,
        "partial_ids": partial_ids,
    }


def _empty_result(section_id, donor_id, n_cortical, elapsed):
    """Return an empty result dict for failed sections."""
    return {
        "section_id": section_id,
        "donor_id": donor_id,
        "n_cells": 0,
        "n_cortical": n_cortical,
        "n_banks": 0,
        "n_pia_segments": 0,
        "n_complete_strips": 0,
        "n_partial_strips": 0,
        "n_cells_selected": 0,
        "coverage_pct": 0,
        "time_sec": elapsed,
        "strip_ids": None,
        "strip_bank": None,
        "complete_ids": set(),
        "partial_ids": set(),
    }


# ── CLI ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Apply curved cortex strips to SEA-AD MERFISH data")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--section", type=str,
                       help="Section name (e.g., H20.33.001.CX28.MTG.02.007.1.02.02)")
    group.add_argument("--all", action="store_true",
                       help="Process all 69 sections")
    parser.add_argument("--strip-width", type=float, default=STRIP_WIDTH_ALONG_PIA,
                        help=f"Strip width along pia in μm (default: {STRIP_WIDTH_ALONG_PIA})")
    parser.add_argument("--save", action="store_true",
                        help="Save strip columns back to MERFISH h5ad")
    args = parser.parse_args()

    # Load MERFISH data
    print(f"Loading MERFISH data from {MERFISH_PATH}...")
    t_load = time.time()
    adata_full = ad.read_h5ad(MERFISH_PATH)
    print(f"  Loaded {adata_full.n_obs:,} cells in {time.time()-t_load:.1f}s")

    sections = sorted(adata_full.obs["Section"].unique())
    print(f"  {len(sections)} sections available")

    if args.section:
        if args.section not in sections:
            print(f"  ERROR: Section '{args.section}' not found.")
            print(f"  Available sections: {sections[:5]}... ({len(sections)} total)")
            sys.exit(1)
        sections_to_process = [args.section]
    else:
        sections_to_process = sections

    print(f"\nProcessing {len(sections_to_process)} sections...")
    t_total = time.time()

    results = []
    for section_id in sections_to_process:
        try:
            # Subset to this section
            sec_mask = adata_full.obs["Section"] == section_id
            adata_section = adata_full[sec_mask]

            result = process_merfish_section(
                adata_section, section_id, strip_width=args.strip_width
            )
            if result:
                results.append(result)
        except Exception as e:
            print(f"\n  ERROR processing {section_id}: {e}")
            import traceback
            traceback.print_exc()
            results.append({
                "section_id": section_id,
                "donor_id": "",
                "n_cells": 0,
                "n_cortical": 0,
                "n_banks": 0,
                "n_pia_segments": 0,
                "n_complete_strips": 0,
                "n_partial_strips": 0,
                "n_cells_selected": 0,
                "coverage_pct": 0,
                "time_sec": 0,
            })

    elapsed_total = time.time() - t_total
    print(f"\n{'='*70}")
    print(f"BATCH SUMMARY ({len(results)} sections, {elapsed_total:.0f}s total)")
    print(f"{'='*70}")

    # Summary table
    df = pd.DataFrame([{
        "section_id": r["section_id"],
        "donor_id": r.get("donor_id", ""),
        "n_cells": r.get("n_cells", 0),
        "n_cortical": r.get("n_cortical", 0),
        "n_banks": r["n_banks"],
        "n_pia_segments": r["n_pia_segments"],
        "n_complete_strips": r["n_complete_strips"],
        "n_partial_strips": r["n_partial_strips"],
        "n_cells_selected": r["n_cells_selected"],
        "coverage_pct": r["coverage_pct"],
        "time_sec": r["time_sec"],
    } for r in results])

    print(df.to_string(index=False))

    # Summary stats
    valid = df[df["coverage_pct"] > 0]
    if len(valid) > 0:
        print(f"\nMean coverage: {valid['coverage_pct'].mean():.1f}%")
        print(f"Median coverage: {valid['coverage_pct'].median():.1f}%")
        print(f"Min coverage: {valid['coverage_pct'].min():.1f}% "
              f"({valid.loc[valid['coverage_pct'].idxmin(), 'section_id']})")
        print(f"Max coverage: {valid['coverage_pct'].max():.1f}% "
              f"({valid.loc[valid['coverage_pct'].idxmax(), 'section_id']})")

    # Save summary CSV
    csv_path = os.path.join(MERFISH_OUT_DIR, "merfish_strips_batch_summary.csv")
    df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    # Save strip columns back to h5ad if requested
    if args.save and len(results) > 0:
        print(f"\nSaving strip columns to MERFISH h5ad...")
        _save_merfish_results(adata_full, results, sections_to_process)


def _save_merfish_results(adata_full, results, sections_processed):
    """Save strip assignments back to the MERFISH h5ad file."""
    # Initialize columns for all cells
    n = adata_full.n_obs
    all_strip_ids = np.full(n, -1, dtype=int)
    all_strip_bank = np.full(n, -1, dtype=int)
    all_in_strip = np.zeros(n, dtype=bool)
    all_strip_tier = np.full(n, "", dtype=object)

    for result in results:
        if result.get("strip_ids") is None:
            continue

        section_id = result["section_id"]
        sec_mask = adata_full.obs["Section"] == section_id
        sec_indices = np.where(sec_mask.values)[0]

        strip_ids = result["strip_ids"]
        strip_bank = result["strip_bank"]
        complete_ids = result["complete_ids"]
        partial_ids = result["partial_ids"]

        # Make strip IDs globally unique by adding section offset
        # (each section's strip IDs start from 0)
        sec_offset = sec_indices[0] * 1000  # crude offset to avoid collisions

        for qi, fi in enumerate(sec_indices):
            sid = strip_ids[qi]
            all_strip_ids[fi] = sid + sec_offset if sid >= 0 else -1
            all_strip_bank[fi] = strip_bank[qi]
            if sid in complete_ids:
                all_in_strip[fi] = True
                all_strip_tier[fi] = "complete"
            elif sid in partial_ids:
                all_in_strip[fi] = True
                all_strip_tier[fi] = "partial"

    adata_full.obs["curved_strip_id"] = all_strip_ids
    adata_full.obs["in_curved_strip"] = all_in_strip
    adata_full.obs["curved_strip_tier"] = all_strip_tier
    adata_full.obs["curved_strip_bank"] = all_strip_bank

    out_path = os.path.join(MERFISH_OUT_DIR, "SEAAD_MTG_MERFISH_with_strips.h5ad")
    print(f"  Writing to {out_path}...")
    adata_full.write_h5ad(out_path)
    print(f"  Saved {n:,} cells with strip annotations")


if __name__ == "__main__":
    main()
