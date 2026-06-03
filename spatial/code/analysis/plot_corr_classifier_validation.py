#!/usr/bin/env python3
"""
Validation figure for the two-stage hierarchical correlation classifier.

6-panel figure comparing HANN (MapMyCells) vs correlation classifier:
  A. L6b in upper layers — bar chart (HANN vs Correlation)
  B. Per-sample L6b improvement (paired scatter, below diagonal = better)
  C. HANN-L6b reclassification (what did misplaced L6b become?)
  D. Correlation margin distribution (histogram, 1% threshold marked)
  E. Spatial map — representative sample, HANN vs Correlation L6b
  F. L6b depth distribution — HANN vs Correlation (overlapping histograms)

Output:
  - output/presentation/slide_corr_classifier_validation.png
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    H5AD_DIR, PRESENTATION_DIR, BG_COLOR, DX_COLORS,
    SAMPLE_TO_DX, EXCLUDE_SAMPLES, CORTICAL_LAYERS,
    LAYER_COLORS, LAYER_ORDER, style_dark_axis, ALL_CELL_COLOR,
    MARKER_SIZE_BG,
)

# SEA-AD subclass colors (same as viewer)
SUBCLASS_COLORS = {
    "L2/3 IT": "#B1EC30", "L4 IT": "#00E5E5", "L5 IT": "#50B2AD",
    "L5 ET": "#0D5B78", "L5/6 NP": "#3E9E64", "L6 IT": "#A19922",
    "L6 IT Car3": "#5100FF", "L6 CT": "#2D8CB8", "L6b": "#7044AA",
    "Sst": "#FF9900", "Sst Chodl": "#B1B10C", "Pvalb": "#D93137",
    "Vip": "#A45FBF", "Lamp5": "#DA808C", "Lamp5 Lhx6": "#935F50",
    "Sncg": "#DF70FF", "Pax6": "#71238C", "Chandelier": "#F641A8",
    "Astrocyte": "#665C47", "Oligodendrocyte": "#53776C", "OPC": "#374A45",
    "Microglia-PVM": "#94AF97", "Endothelial": "#8D6C62", "VLMC": "#697255",
}

# Layer depth bins (matching depth_model.py)
UPPER_LAYERS = {"L1", "L2/3", "L4"}


def load_sample_obs(sample_id):
    """Load obs DataFrame with both HANN and corr labels + depth."""
    fpath = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    adata = ad.read_h5ad(fpath, backed="r")

    cols = ["sample_id", "subclass_label", "qc_pass",
            "predicted_norm_depth", "layer"]
    has_corr = "corr_subclass" in adata.obs.columns
    if has_corr:
        cols += ["corr_subclass", "corr_supertype", "corr_subclass_margin",
                 "corr_qc_pass"]

    obs = adata.obs[cols].copy()
    coords = adata.obsm["spatial"]
    obs["x"] = coords[:, 0]
    obs["y"] = coords[:, 1]

    # Filter to hardware QC-pass
    obs = obs[obs["qc_pass"] == True].copy()
    return obs


def main():
    print("Loading data from all samples...")
    os.makedirs(PRESENTATION_DIR, exist_ok=True)

    all_obs = []
    samples = sorted(s for s in SAMPLE_TO_DX if s not in EXCLUDE_SAMPLES)
    for sid in samples:
        obs = load_sample_obs(sid)
        all_obs.append(obs)
        print(f"  {sid}: {len(obs):,} cells")

    df = pd.concat(all_obs, ignore_index=True)
    print(f"Total: {len(df):,} cells from {len(samples)} samples\n")

    has_corr = "corr_subclass" in df.columns
    if not has_corr:
        print("ERROR: corr_subclass not found in data. Run 02b first.")
        sys.exit(1)

    # ── Compute key metrics ──────────────────────────────────────────
    df["layer"] = df["layer"].astype(str)
    df["depth"] = df["predicted_norm_depth"].astype(float)

    # HANN L6b cells
    hann_l6b = df[df["subclass_label"] == "L6b"].copy()
    hann_l6b_upper = hann_l6b[hann_l6b["layer"].isin(UPPER_LAYERS)]
    pct_hann_upper = 100 * len(hann_l6b_upper) / len(hann_l6b)

    # Corr L6b cells
    corr_l6b = df[df["corr_subclass"] == "L6b"].copy()
    corr_l6b_upper = corr_l6b[corr_l6b["layer"].isin(UPPER_LAYERS)]
    pct_corr_upper = 100 * len(corr_l6b_upper) / len(corr_l6b)

    print(f"HANN L6b: {len(hann_l6b):,} total, {len(hann_l6b_upper):,} in upper layers ({pct_hann_upper:.1f}%)")
    print(f"Corr L6b: {len(corr_l6b):,} total, {len(corr_l6b_upper):,} in upper layers ({pct_corr_upper:.1f}%)")

    # ── Build figure ─────────────────────────────────────────────────
    fig = plt.figure(figsize=(22, 14), facecolor=BG_COLOR, dpi=150)
    gs = gridspec.GridSpec(2, 3, figure=fig,
                           wspace=0.30, hspace=0.35,
                           left=0.06, right=0.97, top=0.92, bottom=0.06)

    panel_labels = ['A', 'B', 'C', 'D', 'E', 'F']

    # ═══════════════════════════════════════════════════════════════════
    # Panel A: L6b in upper layers — bar chart
    # ═══════════════════════════════════════════════════════════════════
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.set_facecolor(BG_COLOR)

    bars = ax_a.bar([0, 1], [pct_hann_upper, pct_corr_upper],
                    color=["#ef5350", "#66bb6a"], width=0.5, edgecolor="white",
                    linewidth=0.5, alpha=0.85)
    ax_a.set_xticks([0, 1])
    ax_a.set_xticklabels(["HANN\n(MapMyCells)", "Correlation\nClassifier"],
                          fontsize=14, color="white")
    ax_a.set_ylabel("L6b in upper layers (%)", fontsize=16, color="white")
    ax_a.set_title("L6b misplacement", fontsize=18, color="white", fontweight="bold")
    ax_a.tick_params(colors="white", labelsize=13)
    ax_a.spines["bottom"].set_color("#555555")
    ax_a.spines["left"].set_color("#555555")
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    ax_a.set_ylim(0, max(pct_hann_upper, pct_corr_upper) * 1.3)

    # Value labels on bars
    for bar, val in zip(bars, [pct_hann_upper, pct_corr_upper]):
        ax_a.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                  f"{val:.1f}%", ha="center", va="bottom", fontsize=16,
                  color="white", fontweight="bold")

    # Arrow showing improvement
    improvement = pct_hann_upper - pct_corr_upper
    ax_a.annotate(f"−{improvement:.1f}pp",
                  xy=(0.5, pct_corr_upper + 1),
                  fontsize=14, color="#66bb6a", ha="center",
                  fontweight="bold")

    # ═══════════════════════════════════════════════════════════════════
    # Panel B: Per-sample L6b upper-layer % (paired scatter)
    # ═══════════════════════════════════════════════════════════════════
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.set_facecolor(BG_COLOR)

    per_sample_hann = []
    per_sample_corr = []
    per_sample_dx = []
    for sid in samples:
        sdf = df[df["sample_id"] == sid]
        h_l6b = sdf[sdf["subclass_label"] == "L6b"]
        c_l6b = sdf[sdf["corr_subclass"] == "L6b"]
        if len(h_l6b) > 0 and len(c_l6b) > 0:
            h_pct = 100 * h_l6b["layer"].isin(UPPER_LAYERS).sum() / len(h_l6b)
            c_pct = 100 * c_l6b["layer"].isin(UPPER_LAYERS).sum() / len(c_l6b)
            per_sample_hann.append(h_pct)
            per_sample_corr.append(c_pct)
            per_sample_dx.append(SAMPLE_TO_DX[sid])

    # Diagonal line (no improvement)
    max_val = max(max(per_sample_hann), max(per_sample_corr)) * 1.15
    ax_b.plot([0, max_val], [0, max_val], "--", color="#666666", linewidth=1, zorder=1)

    # Points colored by diagnosis
    for h, c, dx in zip(per_sample_hann, per_sample_corr, per_sample_dx):
        ax_b.scatter(h, c, color=DX_COLORS[dx], s=80, alpha=0.85,
                     edgecolors="white", linewidths=0.5, zorder=3)

    ax_b.set_xlabel("HANN: L6b in upper layers (%)", fontsize=14, color="white")
    ax_b.set_ylabel("Corr: L6b in upper layers (%)", fontsize=14, color="white")
    ax_b.set_title("Per-sample improvement", fontsize=18, color="white",
                    fontweight="bold")
    ax_b.set_xlim(0, max_val)
    ax_b.set_ylim(0, max_val)
    ax_b.tick_params(colors="white", labelsize=12)
    ax_b.spines["bottom"].set_color("#555555")
    ax_b.spines["left"].set_color("#555555")
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)

    # "Better" label below diagonal
    ax_b.text(max_val * 0.75, max_val * 0.25, "← better",
              fontsize=13, color="#66bb6a", ha="center", fontstyle="italic")

    # Legend
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=DX_COLORS["Control"],
               markersize=9, label="Control", linewidth=0),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=DX_COLORS["SCZ"],
               markersize=9, label="SCZ", linewidth=0),
    ]
    leg = ax_b.legend(handles=legend_handles, fontsize=12, loc="upper left",
                      framealpha=0.85, facecolor="#222222", edgecolor="#555555",
                      labelcolor="white")

    # ═══════════════════════════════════════════════════════════════════
    # Panel C: Reclassification of HANN-L6b cells
    # ═══════════════════════════════════════════════════════════════════
    ax_c = fig.add_subplot(gs[0, 2])
    ax_c.set_facecolor(BG_COLOR)

    # What did HANN-L6b cells become under correlation classifier?
    reclass = hann_l6b["corr_subclass"].value_counts()
    # Show top subclasses
    top_n = 10
    top_reclass = reclass.head(top_n)

    colors_bars = [SUBCLASS_COLORS.get(s, "#888888") for s in top_reclass.index]
    y_pos = np.arange(len(top_reclass))
    bars_c = ax_c.barh(y_pos, top_reclass.values, color=colors_bars,
                        edgecolor="white", linewidth=0.3, alpha=0.85)
    ax_c.set_yticks(y_pos)
    ax_c.set_yticklabels(top_reclass.index, fontsize=12, color="white")
    ax_c.invert_yaxis()
    ax_c.set_xlabel("Number of cells", fontsize=14, color="white")
    ax_c.set_title("HANN-L6b reclassified as...", fontsize=18, color="white",
                    fontweight="bold")
    ax_c.tick_params(colors="white", labelsize=12)
    ax_c.spines["bottom"].set_color("#555555")
    ax_c.spines["left"].set_color("#555555")
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)

    # Percentage labels
    total_hann_l6b = len(hann_l6b)
    for bar, val in zip(bars_c, top_reclass.values):
        pct = 100 * val / total_hann_l6b
        ax_c.text(bar.get_width() + total_hann_l6b * 0.01,
                  bar.get_y() + bar.get_height() / 2,
                  f"{pct:.1f}%", ha="left", va="center",
                  fontsize=11, color="#dddddd")

    # ═══════════════════════════════════════════════════════════════════
    # Panel D: Correlation margin distribution
    # ═══════════════════════════════════════════════════════════════════
    ax_d = fig.add_subplot(gs[1, 0])
    ax_d.set_facecolor(BG_COLOR)

    margins = df["corr_subclass_margin"].dropna().values
    qc_fail = df[df["corr_qc_pass"] == False]
    n_fail = len(qc_fail)
    pct_fail = 100 * n_fail / len(df)

    # Histogram
    bins = np.linspace(0, np.percentile(margins, 99.5), 80)
    ax_d.hist(margins, bins=bins, color="#4fc3f7", alpha=0.7,
              edgecolor="none", density=True)

    # Per-sample 1st percentile threshold (show range)
    thresholds = []
    for sid in samples:
        sdf = df[df["sample_id"] == sid]
        m = sdf["corr_subclass_margin"].dropna().values
        if len(m) > 0:
            thresholds.append(np.percentile(m, 1.0))

    thresh_mean = np.mean(thresholds)
    thresh_min = np.min(thresholds)
    thresh_max = np.max(thresholds)

    # Shade the QC-fail region
    ax_d.axvspan(0, thresh_max, color="#ef5350", alpha=0.15, zorder=0)
    ax_d.axvline(thresh_mean, color="#ef5350", linewidth=2, linestyle="--",
                 label=f"1st pctl (mean={thresh_mean:.3f})")

    ax_d.set_xlabel("Subclass correlation margin", fontsize=14, color="white")
    ax_d.set_ylabel("Density", fontsize=14, color="white")
    ax_d.set_title("Margin distribution", fontsize=18, color="white",
                    fontweight="bold")
    ax_d.tick_params(colors="white", labelsize=12)
    ax_d.spines["bottom"].set_color("#555555")
    ax_d.spines["left"].set_color("#555555")
    ax_d.spines["top"].set_visible(False)
    ax_d.spines["right"].set_visible(False)

    # QC-fail annotation
    ax_d.text(thresh_mean + 0.005, ax_d.get_ylim()[1] * 0.85,
              f"QC-fail: {n_fail:,}\n({pct_fail:.1f}%)",
              fontsize=13, color="#ef5350", fontweight="bold")

    leg_d = ax_d.legend(fontsize=11, loc="upper right",
                        framealpha=0.85, facecolor="#222222",
                        edgecolor="#555555", labelcolor="white")

    # ═══════════════════════════════════════════════════════════════════
    # Panel E: Spatial map — HANN vs Corr L6b for representative sample
    # ═══════════════════════════════════════════════════════════════════
    # Use sample with biggest L6b improvement
    best_improvement = 0
    best_sid = samples[0]
    for sid, h, c in zip(samples, per_sample_hann, per_sample_corr):
        if h - c > best_improvement:
            best_improvement = h - c
            best_sid = sid

    sdf = df[df["sample_id"] == best_sid].copy()

    # Left sub-panel: HANN L6b
    gs_e = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, 1],
                                            wspace=0.05)
    ax_e1 = fig.add_subplot(gs_e[0])
    ax_e2 = fig.add_subplot(gs_e[1])

    for ax_e, label_col, subtitle in [(ax_e1, "subclass_label", "HANN"),
                                       (ax_e2, "corr_subclass", "Corr")]:
        style_dark_axis(ax_e)

        # Background cells
        non_l6b = sdf[sdf[label_col] != "L6b"]
        ax_e.scatter(non_l6b["x"], non_l6b["y"], s=MARKER_SIZE_BG,
                     c=ALL_CELL_COLOR, alpha=0.15, rasterized=True,
                     linewidths=0, zorder=1)

        # L6b cells in deep layers (correct)
        l6b_cells = sdf[sdf[label_col] == "L6b"]
        l6b_deep = l6b_cells[~l6b_cells["layer"].isin(UPPER_LAYERS)]
        ax_e.scatter(l6b_deep["x"], l6b_deep["y"], s=6,
                     c=SUBCLASS_COLORS["L6b"], alpha=0.8, rasterized=True,
                     linewidths=0, zorder=3)

        # L6b cells in upper layers (misplaced — red)
        l6b_upper = l6b_cells[l6b_cells["layer"].isin(UPPER_LAYERS)]
        ax_e.scatter(l6b_upper["x"], l6b_upper["y"], s=12,
                     c="#ff1744", alpha=0.9, rasterized=True,
                     linewidths=0, zorder=4)

        n_deep = len(l6b_deep)
        n_up = len(l6b_upper)
        pct_up = 100 * n_up / (n_deep + n_up) if (n_deep + n_up) > 0 else 0
        ax_e.set_title(f"{subtitle}\n{n_up} in upper ({pct_up:.1f}%)",
                       fontsize=13, color="white")

    # Suptitle for panel E
    fig.text(gs[1, 1].get_position(fig).x0 + 0.01,
             gs[1, 1].get_position(fig).y1 + 0.01,
             f"L6b spatial: {best_sid}", fontsize=16, color="white",
             fontweight="bold", transform=fig.transFigure)

    # ═══════════════════════════════════════════════════════════════════
    # Panel F: L6b depth distribution — HANN vs Correlation
    # ═══════════════════════════════════════════════════════════════════
    ax_f = fig.add_subplot(gs[1, 2])
    ax_f.set_facecolor(BG_COLOR)

    hann_depth = hann_l6b["depth"].dropna().values
    corr_depth = corr_l6b["depth"].dropna().values

    depth_bins = np.linspace(0, 1, 60)
    ax_f.hist(hann_depth, bins=depth_bins, color="#ef5350", alpha=0.5,
              label=f"HANN L6b (n={len(hann_depth):,})", density=True,
              edgecolor="none")
    ax_f.hist(corr_depth, bins=depth_bins, color="#66bb6a", alpha=0.5,
              label=f"Corr L6b (n={len(corr_depth):,})", density=True,
              edgecolor="none")

    # MERFISH reference L6b depth
    merfish_l6b_depth = 0.800  # from prior analysis
    ax_f.axvline(merfish_l6b_depth, color="#ffd54f", linewidth=2, linestyle="--",
                 label=f"MERFISH L6b mean ({merfish_l6b_depth})")

    # Mean depth lines
    hann_mean = np.nanmean(hann_depth)
    corr_mean = np.nanmean(corr_depth)
    ax_f.axvline(hann_mean, color="#ef5350", linewidth=1.5, linestyle=":")
    ax_f.axvline(corr_mean, color="#66bb6a", linewidth=1.5, linestyle=":")

    ax_f.text(hann_mean - 0.02, ax_f.get_ylim()[1] * 0.5 if ax_f.get_ylim()[1] > 0 else 3.0,
              f"μ={hann_mean:.3f}", fontsize=11, color="#ef5350",
              ha="right", rotation=90)
    ax_f.text(corr_mean + 0.02, ax_f.get_ylim()[1] * 0.5 if ax_f.get_ylim()[1] > 0 else 3.0,
              f"μ={corr_mean:.3f}", fontsize=11, color="#66bb6a",
              ha="left", rotation=90)

    ax_f.set_xlabel("Predicted normalized depth", fontsize=14, color="white")
    ax_f.set_ylabel("Density", fontsize=14, color="white")
    ax_f.set_title("L6b depth shift", fontsize=18, color="white",
                    fontweight="bold")
    ax_f.tick_params(colors="white", labelsize=12)
    ax_f.spines["bottom"].set_color("#555555")
    ax_f.spines["left"].set_color("#555555")
    ax_f.spines["top"].set_visible(False)
    ax_f.spines["right"].set_visible(False)

    leg_f = ax_f.legend(fontsize=11, loc="upper left",
                        framealpha=0.85, facecolor="#222222",
                        edgecolor="#555555", labelcolor="white")

    # Layer depth guides
    layer_boundaries = [0.0, 0.1225, 0.4696, 0.5443, 0.7079, 0.9275, 1.0]
    layer_names = ["L1", "L2/3", "L4", "L5", "L6", "WM"]
    for i, (lb, ln) in enumerate(zip(layer_boundaries[:-1], layer_names)):
        mid = (lb + layer_boundaries[i + 1]) / 2
        ax_f.axvspan(lb, layer_boundaries[i + 1],
                     color=LAYER_COLORS.get(ln, (0.5, 0.5, 0.5)),
                     alpha=0.08, zorder=0)

    # ═══════════════════════════════════════════════════════════════════
    # Panel labels
    # ═══════════════════════════════════════════════════════════════════
    axes_for_labels = [ax_a, ax_b, ax_c, ax_d, ax_e1, ax_f]
    for ax, label in zip(axes_for_labels, panel_labels):
        ax.text(-0.08, 1.05, label, transform=ax.transAxes,
                fontsize=22, fontweight="bold", color="white",
                va="bottom", ha="right")

    # ── Suptitle ─────────────────────────────────────────────────────
    fig.suptitle("Correlation Classifier vs HANN (MapMyCells)",
                 fontsize=24, color="white", fontweight="bold", y=0.97)

    # ── Save ─────────────────────────────────────────────────────────
    out_path = os.path.join(PRESENTATION_DIR, "slide_corr_classifier_validation.png")
    fig.savefig(out_path, dpi=200, facecolor=BG_COLOR, bbox_inches="tight")
    print(f"\nSaved: {out_path}")
    plt.close()

    # ── Print summary stats ──────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"  HANN L6b in upper layers: {pct_hann_upper:.1f}% ({len(hann_l6b_upper):,}/{len(hann_l6b):,})")
    print(f"  Corr L6b in upper layers: {pct_corr_upper:.1f}% ({len(corr_l6b_upper):,}/{len(corr_l6b):,})")
    print(f"  Improvement: {improvement:.1f} percentage points")
    print(f"  HANN L6b mean depth: {hann_mean:.3f}")
    print(f"  Corr L6b mean depth: {corr_mean:.3f}")
    print(f"  MERFISH L6b ref depth: {merfish_l6b_depth}")
    print(f"  QC-fail cells: {n_fail:,} / {len(df):,} ({pct_fail:.1f}%)")
    print(f"  Margin 1st pctl range: [{thresh_min:.4f}, {thresh_max:.4f}]")

    # Agreement between HANN and Corr subclass
    agree = (df["subclass_label"] == df["corr_subclass"]).sum()
    print(f"  Subclass agreement: {100*agree/len(df):.1f}% ({agree:,}/{len(df):,})")


if __name__ == "__main__":
    main()
