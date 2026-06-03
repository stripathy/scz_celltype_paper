#!/usr/bin/env python3
"""
Supertype-level depth distributions: MERFISH vs Xenium.

For each subclass with >3 supertypes, plots paired violin distributions
(MERFISH green, Xenium orange) for every supertype. Supertypes are ordered
by median MERFISH depth (shallowest → deepest, left → right).

Y-axis: cortical depth (pia at top, WM at bottom).
One figure per cell class (Glutamatergic, GABAergic, Non-neuronal).

Output: output/presentation/supertype_depth_violins_{class}.png
"""

import os
import sys
import pickle
import argparse
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    H5AD_DIR, PRESENTATION_DIR, EXCLUDE_SAMPLES,
    CORTICAL_LAYERS, SUBCLASS_TO_CLASS, load_merfish_cortical, load_cells,
)

DEPTH_MODEL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "..", "output", "depth_model_normalized.pkl")

OUT_DIR = PRESENTATION_DIR
os.makedirs(OUT_DIR, exist_ok=True)

MIN_SUPERTYPES = 4  # only show subclasses with >= this many supertypes

# Subclass ordering within each class
GLUT_ORDER = ["L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP",
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b"]
GABA_ORDER = ["Lamp5", "Lamp5 Lhx6", "Sncg", "Vip", "Pax6",
              "Chandelier", "Pvalb", "Sst", "Sst Chodl"]
NN_ORDER = ["Astrocyte", "Oligodendrocyte", "OPC",
            "Microglia-PVM", "Endothelial", "VLMC"]

CLASS_SUBCLASS_ORDER = {
    "Glutamatergic": GLUT_ORDER,
    "GABAergic": GABA_ORDER,
    "Non-neuronal": NN_ORDER,
}

# Colors
MERFISH_COLOR = "#66BB6A"
XENIUM_COLOR = "#FF9800"

CLASS_BG = {
    "Glutamatergic": "#FFF8F0",
    "GABAergic": "#F0F8F0",
    "Non-neuronal": "#F0F4FA",
}

# Layer boundary positions and labels
LAYER_BOUNDARIES = [0.10, 0.40, 0.55, 0.70, 0.90]
LAYER_TICKS = [0.05, 0.25, 0.475, 0.625, 0.80, 0.95]
LAYER_NAMES = ["L1", "L2/3", "L4", "L5", "L6", "WM"]


def load_xenium_cortical():
    """Load all Xenium cortical cells with corr labels + depth via standard load_cells()."""
    print("Loading Xenium cortical cells...")
    h5ad_files = sorted(
        f for f in os.listdir(H5AD_DIR) if f.endswith("_annotated.h5ad")
    )
    rows = []
    for fname in h5ad_files:
        sample_id = fname.replace("_annotated.h5ad", "")
        if sample_id in EXCLUDE_SAMPLES:
            continue
        obs = load_cells(sample_id, cortical_only=True,
                         extra_obs_columns=["predicted_norm_depth"])
        # Rename back to corr_* for downstream compatibility within this script
        obs = obs.rename(columns={"subclass_label": "corr_subclass",
                                  "supertype_label": "corr_supertype"})
        rows.append(obs)
        print(f"  {sample_id}: {len(obs):,} cortical cells")

    df = pd.concat(rows, ignore_index=True)
    df["predicted_norm_depth"] = pd.to_numeric(df["predicted_norm_depth"],
                                                errors="coerce")
    print(f"  Total: {len(df):,} cortical cells from {len(rows)} samples")
    return df


def _shorten_supertype(sup, subclass):
    """Strip subclass prefix from supertype name for compact labels."""
    short = sup
    # Try common prefix patterns
    for prefix in [f"{subclass}_", f"{subclass} "]:
        if short.startswith(prefix):
            short = short[len(prefix):]
            break
    # For names like "Astro_1", "Oligo_2", strip the abbreviated prefix too
    abbreviated_prefixes = {
        "Astrocyte": "Astro_", "Oligodendrocyte": "Oligo_",
        "Microglia-PVM": "Micro-PVM_", "OPC": "OPC_",
    }
    if subclass in abbreviated_prefixes:
        pfx = abbreviated_prefixes[subclass]
        if sup.startswith(pfx):
            short = sup[len(pfx):]
    return short


def plot_class_figure(merfish_df, xenium_df, cell_class, subclass_order):
    """
    One figure for a cell class using GridSpec for variable row heights.
    Rows = subclasses, paired violins per supertype.
    """
    # Find qualifying subclasses
    qualifying = []
    sub_supertypes = {}
    for sub in subclass_order:
        if SUBCLASS_TO_CLASS.get(sub) != cell_class:
            continue
        m_sups = merfish_df.loc[merfish_df["subclass"] == sub, "supertype"].unique()
        if len(m_sups) >= MIN_SUPERTYPES:
            qualifying.append(sub)
            sub_supertypes[sub] = m_sups

    if not qualifying:
        print(f"  No qualifying subclasses for {cell_class}")
        return None

    n_rows = len(qualifying)

    # Row heights proportional to number of supertypes (more types = needs more width,
    # but since we use a fixed width, give more vertical breathing room to dense rows)
    n_sups_per_row = [len(sub_supertypes[sub]) for sub in qualifying]
    max_sups = max(n_sups_per_row)

    # Base row height scales with number of supertypes (min 3.2, max 5.5)
    row_heights = [max(3.2, 3.2 + (n - 4) * 0.15) for n in n_sups_per_row]

    fig_width = max(18, max_sups * 1.8 + 4)
    fig_height = sum(row_heights) + 1.5  # extra for suptitle

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = gridspec.GridSpec(n_rows, 1, figure=fig, height_ratios=row_heights,
                           hspace=0.45)

    bg_color = CLASS_BG.get(cell_class, "white")

    for row_idx, sub in enumerate(qualifying):
        ax = fig.add_subplot(gs[row_idx, 0])
        ax.set_facecolor(bg_color)

        # Get supertypes present in MERFISH for this subclass
        m_sub = merfish_df[merfish_df["subclass"] == sub]
        m_supertypes = m_sub["supertype"].unique()

        # Order by median MERFISH depth (shallowest first)
        sup_median_depth = {}
        for sup in m_supertypes:
            depths = m_sub.loc[m_sub["supertype"] == sup, "depth"].values
            sup_median_depth[sup] = np.median(depths) if len(depths) > 0 else 0.5
        ordered_sups = sorted(m_supertypes, key=lambda s: sup_median_depth[s])

        # Get Xenium data for this subclass
        x_sub = xenium_df[xenium_df["corr_subclass"].astype(str) == sub]

        n_sups = len(ordered_sups)

        # Spacing: each supertype pair occupies 1 unit, with pairs centered at 0,1,2,...
        pair_offset = 0.22
        violin_width = 0.38

        for i, sup in enumerate(ordered_sups):
            pos_m = i - pair_offset
            pos_x = i + pair_offset

            # MERFISH depths
            m_depths = m_sub.loc[m_sub["supertype"] == sup, "depth"].values
            # Xenium depths
            x_depths = x_sub.loc[
                x_sub["corr_supertype"].astype(str) == sup,
                "predicted_norm_depth"
            ].dropna().values

            # MERFISH violin
            if len(m_depths) >= 5:
                vp_m = ax.violinplot([m_depths], positions=[pos_m],
                                      showmedians=True, showextrema=False,
                                      widths=violin_width)
                for body in vp_m["bodies"]:
                    body.set_facecolor(MERFISH_COLOR)
                    body.set_alpha(0.75)
                    body.set_edgecolor("white")
                    body.set_linewidth(0.5)
                vp_m["cmedians"].set_color("#333333")
                vp_m["cmedians"].set_linewidth(1.5)
            elif len(m_depths) > 0:
                # Too few for violin — show individual points
                ax.scatter([pos_m] * len(m_depths), m_depths,
                           color=MERFISH_COLOR, s=8, alpha=0.6, zorder=3)

            # Xenium violin
            if len(x_depths) >= 5:
                vp_x = ax.violinplot([x_depths], positions=[pos_x],
                                      showmedians=True, showextrema=False,
                                      widths=violin_width)
                for body in vp_x["bodies"]:
                    body.set_facecolor(XENIUM_COLOR)
                    body.set_alpha(0.75)
                    body.set_edgecolor("white")
                    body.set_linewidth(0.5)
                vp_x["cmedians"].set_color("#333333")
                vp_x["cmedians"].set_linewidth(1.5)
            elif len(x_depths) > 0:
                ax.scatter([pos_x] * len(x_depths), x_depths,
                           color=XENIUM_COLOR, s=8, alpha=0.6, zorder=3)

        # Light vertical separators between supertype pairs
        for i in range(1, n_sups):
            ax.axvline(i - 0.5, color="#cccccc", linewidth=0.5, alpha=0.5)

        # X-axis: supertype labels
        tick_positions = np.arange(n_sups)
        tick_labels = []
        for sup in ordered_sups:
            n_m = len(m_sub[m_sub["supertype"] == sup])
            n_x = len(x_sub[x_sub["corr_supertype"].astype(str) == sup])
            short = _shorten_supertype(sup, sub)
            tick_labels.append(f"{short}\n({n_m:,} / {n_x:,})")

        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=12, ha="center",
                           linespacing=0.9)
        ax.set_xlim(-0.6, n_sups - 0.4)

        # Y-axis: depth (inverted — pia at top)
        ax.set_ylim(1.02, -0.02)  # inverted
        ax.set_yticks(LAYER_TICKS)
        ax.set_yticklabels(LAYER_NAMES, fontsize=14, fontweight="medium")

        # Layer boundary lines
        for boundary in LAYER_BOUNDARIES:
            ax.axhline(boundary, color="#999999", linewidth=0.7, alpha=0.5,
                       linestyle="--")

        # Only put y-label on middle row
        if row_idx == n_rows // 2:
            ax.set_ylabel("Cortical Depth", fontsize=16, fontweight="bold")

        # Subclass title — bold, left-aligned
        ax.set_title(f"{sub}  ({n_sups} supertypes)", fontsize=20,
                     fontweight="bold", loc="left", pad=8)

        # Remove top/right spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Legend on first row
    legend_elements = [
        Patch(facecolor=MERFISH_COLOR, alpha=0.75, edgecolor="white",
              label="MERFISH (manual depth)"),
        Patch(facecolor=XENIUM_COLOR, alpha=0.75, edgecolor="white",
              label="Xenium (predicted depth)"),
    ]
    fig.axes[0].legend(handles=legend_elements, loc="upper right", fontsize=14,
                        framealpha=0.95, edgecolor="#cccccc")

    # Footnote
    fig.text(0.01, -0.005,
             "Counts shown as (MERFISH / Xenium). Supertypes ordered by median MERFISH depth.",
             fontsize=12, color="#666666", ha="left", style="italic")

    fig.suptitle(f"Supertype Depth Distributions — {cell_class}",
                 fontsize=24, fontweight="bold")
    return fig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--test-only", action="store_true",
                        help="Use only held-out MERFISH test donors (stricter validation)")
    args = parser.parse_args()

    # Load data
    merfish_df = load_merfish_cortical()
    xenium_df = load_xenium_cortical()

    suffix = ""
    if args.test_only:
        # Filter to held-out test donors only
        with open(DEPTH_MODEL_PATH, "rb") as f:
            bundle = pickle.load(f)
        test_donors = [str(d) for d in bundle["test_donors"]]
        n_before = len(merfish_df)
        merfish_df = merfish_df[merfish_df["donor"].isin(test_donors)].copy()
        print(f"\n  --test-only: Filtering to {len(test_donors)} held-out donors")
        print(f"  {n_before:,} → {len(merfish_df):,} MERFISH cells")
        suffix = "_testonly"

    print(f"\nMERFISH: {len(merfish_df):,} cells, "
          f"{merfish_df['supertype'].nunique()} supertypes")
    print(f"Xenium:  {len(xenium_df):,} cells, "
          f"{xenium_df['corr_supertype'].nunique()} supertypes")

    # One figure per class
    for cell_class, subclass_order in CLASS_SUBCLASS_ORDER.items():
        print(f"\n{'='*60}")
        print(f"Plotting {cell_class}...")
        fig = plot_class_figure(merfish_df, xenium_df, cell_class, subclass_order)
        if fig is not None:
            tag = cell_class.lower().replace("-", "")
            out_path = os.path.join(OUT_DIR,
                                     f"supertype_depth_violins_{tag}{suffix}.png")
            fig.savefig(out_path, dpi=200, bbox_inches="tight")
            print(f"  Saved: {out_path}")
            plt.close(fig)


if __name__ == "__main__":
    main()
