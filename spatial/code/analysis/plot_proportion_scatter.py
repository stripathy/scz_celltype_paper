#!/usr/bin/env python3
"""
Unified proportion scatter: MERFISH vs Xenium cell type proportions.

Consolidates three former scripts into one with composable CLI flags:
  plot_predicted_proportion_scatter.py  (neurons-only, error bars, weighted corr)
  plot_cropped_proportions.py           (all cells, manual MERFISH, uncrop vs crop)
  plot_merfish_vs_xenium_proportions.py (all cells, subclass+supertype, uncrop only)

Usage examples (reproducing original outputs):

  # Reproduce plot_predicted_proportion_scatter.py
  python plot_proportion_scatter.py \
      --filter neurons --level supertype --merfish-subset all \
      --crop both --error-bars --weighted-dots

  # Reproduce plot_cropped_proportions.py
  python plot_proportion_scatter.py \
      --filter all --level supertype --merfish-subset manual \
      --crop both --no-error-bars --no-weighted-dots

  # Reproduce plot_merfish_vs_xenium_proportions.py
  python plot_proportion_scatter.py \
      --filter all --level both --merfish-subset all \
      --crop uncropped --no-error-bars --no-weighted-dots

Output: output/presentation/proportion_scatter_<tag>.png
        output/presentation/proportion_scatter_<tag>.csv
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BG_COLOR, H5AD_DIR, MERFISH_PATH, PRESENTATION_DIR,
    CONTROL_SAMPLES, CORTICAL_LAYERS, EXCLUDE_SAMPLES,
    CLASS_COLORS, classify_celltype, load_cells, load_merfish_cortical,
)

BG = BG_COLOR


# ──────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────

def is_neuron(ct):
    """Return True if *ct* is a neuronal cell type."""
    _, cls = classify_celltype(ct)
    return cls in ("Glutamatergic", "GABAergic")


def weighted_corr(x, y, w):
    """Weighted Pearson correlation."""
    w = w / w.sum()
    mx = np.sum(w * x)
    my = np.sum(w * y)
    dx = x - mx
    dy = y - my
    cov = np.sum(w * dx * dy)
    sx = np.sqrt(np.sum(w * dx**2))
    sy = np.sqrt(np.sum(w * dy**2))
    if sx == 0 or sy == 0:
        return np.nan
    return cov / (sx * sy)


def pct_formatter(val, pos):
    """Human-readable percentage tick formatter for log-scale axes."""
    pct = val * 100
    if pct >= 1:
        return f"{pct:.0f}%"
    elif pct >= 0.1:
        return f"{pct:.1f}%"
    elif pct >= 0.01:
        return f"{pct:.2f}%"
    else:
        return f"{pct:.3f}%"


def _per_donor_proportions(obs, donor_col, celltype_col):
    """Compute per-donor proportions, return long DataFrame.

    Returns columns: donor, celltype, proportion.
    """
    records = []
    for donor in obs[donor_col].unique():
        donor_df = obs[obs[donor_col] == donor]
        total = len(donor_df)
        counts = donor_df[celltype_col].value_counts()
        for ct, n in counts.items():
            records.append({"donor": donor, "celltype": str(ct),
                            "proportion": n / total})
    return pd.DataFrame(records)


# ──────────────────────────────────────────────────────────────────────
# MERFISH data loading
# ──────────────────────────────────────────────────────────────────────

def load_merfish_proportions(merfish_subset, level, return_std=False):
    """Load MERFISH proportions at the requested taxonomy level.

    Parameters
    ----------
    merfish_subset : str
        'all'    – all cells with Layer annotation in CORTICAL_LAYERS
                   (direct h5ad read, used by predicted / merfish_vs_xenium).
        'manual' – cells with manual depth annotation only
                   (load_merfish_cortical(), used by cropped_proportions).
    level : str
        'supertype' or 'subclass'.
    return_std : bool
        If True, also return per-celltype std across donors.

    Returns
    -------
    DataFrame with columns: celltype, merfish_mean  (and merfish_std if requested).
    """
    level_cap = level.capitalize()  # 'Supertype' or 'Subclass'

    if merfish_subset == "manual":
        print(f"  MERFISH manually annotated ({level_cap})...")
        obs = load_merfish_cortical()
        celltype_col = level  # 'supertype' or 'subclass'
        donor_col = "donor"
    else:
        print(f"  MERFISH all cortical ({level_cap})...")
        merfish = ad.read_h5ad(MERFISH_PATH, backed="r")
        obs = merfish.obs[["Donor ID", level_cap, "Layer annotation"]].copy()
        obs["Layer annotation"] = obs["Layer annotation"].astype(str)
        obs = obs[obs["Layer annotation"].isin(CORTICAL_LAYERS)]
        obs = obs.rename(columns={"Donor ID": "donor", level_cap: "celltype"})
        obs["celltype"] = obs["celltype"].astype(str)
        celltype_col = "celltype"
        donor_col = "donor"

    n_cells = len(obs)
    n_donors = obs[donor_col].nunique()
    print(f"    {n_cells:,} cells from {n_donors} donors")

    # Rename celltype column to 'celltype' if needed
    if celltype_col != "celltype":
        obs = obs.rename(columns={celltype_col: "celltype"})

    props = _per_donor_proportions(obs, donor_col, "celltype")
    stats = props.groupby("celltype")["proportion"].agg(
        ["mean", "std"]).reset_index()
    stats.columns = ["celltype", "merfish_mean", "merfish_std"]

    if not return_std:
        stats = stats[["celltype", "merfish_mean"]]

    return stats


# ──────────────────────────────────────────────────────────────────────
# Xenium data loading
# ──────────────────────────────────────────────────────────────────────

def load_xenium_proportions(level, crop, return_std=False):
    """Load Xenium control-sample proportions.

    Parameters
    ----------
    level : str
        'supertype' or 'subclass'.
    crop : str
        'uncropped' – all QC-pass cells.
        'cropped'   – cortical layers only (L1-L6).
    return_std : bool
        If True, also return per-celltype std across samples.

    Returns
    -------
    DataFrame with columns: celltype, xenium_mean (and xenium_std if requested).
    """
    col_map = {"supertype": "supertype_label", "subclass": "subclass_label"}
    use_col = col_map[level]

    cortical_only = (crop == "cropped")
    tag = "L1-L6 cropped" if cortical_only else "uncropped"
    print(f"  Xenium controls ({level}, {tag})...")

    controls = [s for s in CONTROL_SAMPLES if s not in EXCLUDE_SAMPLES]
    records = []
    total_cells = 0

    for sample_id in controls:
        obs = load_cells(sample_id, cortical_only=cortical_only)
        if cortical_only:
            obs = obs[obs["layer"].isin(CORTICAL_LAYERS)]
        total = len(obs)
        total_cells += total
        counts = obs[use_col].value_counts()
        for ct, n in counts.items():
            records.append({"donor": sample_id, "celltype": str(ct),
                            "proportion": n / total})

    print(f"    {total_cells:,} cells from {len(controls)} samples")

    df = pd.DataFrame(records)
    stats = df.groupby("celltype")["proportion"].agg(
        ["mean", "std"]).reset_index()
    stats.columns = ["celltype", "xenium_mean", "xenium_std"]

    if not return_std:
        stats = stats[["celltype", "xenium_mean"]]

    return stats


# ──────────────────────────────────────────────────────────────────────
# Scatter plot
# ──────────────────────────────────────────────────────────────────────

def plot_scatter(ax, df, x_col, y_col, title, xlabel, ylabel,
                 x_sd_col=None, y_sd_col=None,
                 use_error_bars=False, use_weighted_dots=False,
                 max_labels=15):
    """Plot MERFISH vs Xenium proportion scatter on log-log axes.

    Parameters
    ----------
    ax : matplotlib Axes
    df : DataFrame
        Must have columns: celltype, *x_col*, *y_col*, and optionally
        *x_sd_col*, *y_sd_col*.
    x_col, y_col : str
        Column names for x and y mean proportions.
    title, xlabel, ylabel : str
    x_sd_col, y_sd_col : str or None
        Column names for standard deviations (needed when use_error_bars=True).
    use_error_bars : bool
        Draw +/-SD error bars per point.
    use_weighted_dots : bool
        Scale dot size by geometric mean frequency and compute
        frequency-weighted Pearson r.
    max_labels : int
        Number of most-deviant cell type labels to show.
    """
    # Filter to valid rows
    valid_mask = df[x_col].notna() & df[y_col].notna()
    valid_mask &= (df[x_col] > 0) & (df[y_col] > 0)
    d = df[valid_mask].copy()

    x = d[x_col].values
    y = d[y_col].values
    colors = [classify_celltype(ct)[0] for ct in d["celltype"]]

    # Dot sizes
    if use_weighted_dots:
        geo_mean = np.sqrt(x * y)
        sizes = 30 + 800 * (geo_mean / geo_mean.max())
    else:
        sizes = 70

    # Error bars
    if use_error_bars and x_sd_col and y_sd_col:
        x_sd = d[x_sd_col].fillna(0).values
        y_sd = d[y_sd_col].fillna(0).values
        for i in range(len(x)):
            ax.errorbar(x[i], y[i], xerr=x_sd[i], yerr=y_sd[i],
                        fmt="none", ecolor=colors[i], alpha=0.3, linewidth=1.0,
                        capsize=0, zorder=2)

    ax.scatter(x, y, c=colors, s=sizes, alpha=0.8, edgecolors="white",
               linewidths=0.5, zorder=5)

    # Label most deviant points
    d["log_ratio"] = np.log2((y + 1e-7) / (x + 1e-7))
    d["abs_log_ratio"] = d["log_ratio"].abs()
    top = d.nlargest(max_labels, "abs_log_ratio")

    for _, row in top.iterrows():
        ax.annotate(row["celltype"],
                    (row[x_col], row[y_col]),
                    fontsize=9, color="#dddddd", alpha=0.9,
                    xytext=(5, 5), textcoords="offset points")

    # Diagonal line
    lo = min(x.min(), y.min()) * 0.3
    hi = max(x.max(), y.max()) * 3
    ax.plot([lo, hi], [lo, hi], "--", color="#888888", linewidth=1.5,
            alpha=0.6, zorder=1)

    # Correlations
    log_x = np.log10(x)
    log_y = np.log10(y)
    r_log, _ = pearsonr(log_x, log_y)
    rho, _ = spearmanr(x, y)

    if use_weighted_dots and use_error_bars:
        # Full stats block: weighted r, CV
        geo_mean = np.sqrt(x * y)
        r_weighted = weighted_corr(log_x, log_y, geo_mean)
        x_sd = d[x_sd_col].fillna(0).values if x_sd_col else np.zeros_like(x)
        y_sd = d[y_sd_col].fillna(0).values if y_sd_col else np.zeros_like(y)
        cv_x = x_sd / (x + 1e-10)
        cv_y = y_sd / (y + 1e-10)
        mean_cv_x = np.mean(cv_x[x_sd > 0]) if np.any(x_sd > 0) else 0
        mean_cv_y = np.mean(cv_y[y_sd > 0]) if np.any(y_sd > 0) else 0

        ax.text(0.04, 0.96,
                f"r = {r_log:.2f} (log)\n"
                f"r_w = {r_weighted:.2f} (freq-weighted)\n"
                f"\u03c1 = {rho:.2f}\n"
                f"n = {len(x)} types\n"
                f"CV: MERFISH={mean_cv_x:.2f}, Xenium={mean_cv_y:.2f}",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=13, color="#dddddd",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#333333",
                          edgecolor="#555555", alpha=0.85))
    else:
        ax.text(0.04, 0.96,
                f"r = {r_log:.2f} (log-scale)\n\u03c1 = {rho:.2f}\n"
                f"n = {len(x)} types",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=14, color="#dddddd",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#333333",
                          edgecolor="#555555", alpha=0.85))

    # Axis formatting
    ax.set_xscale("log")
    ax.set_yscale("log")
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_formatter(mticker.FuncFormatter(pct_formatter))
        axis.set_minor_formatter(mticker.NullFormatter())

    ax.set_xlabel(xlabel, fontsize=16, color="white")
    ax.set_ylabel(ylabel, fontsize=16, color="white")
    ax.set_title(title, fontsize=20, fontweight="bold", color="white", pad=10)
    ax.set_facecolor(BG)
    ax.tick_params(colors="white", labelsize=13)
    for spine in ax.spines.values():
        spine.set_color("#555555")
    ax.grid(True, alpha=0.2, color="#555555")


# ──────────────────────────────────────────────────────────────────────
# Panel assembly helpers
# ──────────────────────────────────────────────────────────────────────

def _build_legend_classes(cell_filter):
    """Build legend handles for the given cell filter."""
    if cell_filter == "neurons":
        entries = [
            ("Glutamatergic", CLASS_COLORS["Glutamatergic"]),
            ("GABAergic", CLASS_COLORS["GABAergic"]),
        ]
    else:
        entries = [
            ("Glutamatergic", CLASS_COLORS["Glutamatergic"]),
            ("GABAergic", CLASS_COLORS["GABAergic"]),
            ("Non-neuronal", CLASS_COLORS["Non-neuronal"]),
        ]
    return [
        Line2D([0], [0], marker="o", color=BG, markerfacecolor=color,
               markersize=12, label=label, linewidth=0)
        for label, color in entries
    ]


def _make_output_tag(args):
    """Build a filename tag from CLI args for easy identification."""
    parts = [args.filter, args.merfish_subset]
    if args.level == "both":
        parts.append("sub_sup")
    else:
        parts.append(args.level)
    parts.append(args.crop)
    if args.error_bars:
        parts.append("eb")
    if args.weighted_dots:
        parts.append("wd")
    return "_".join(parts)


def _merfish_xlabel(merfish_subset):
    """Build x-axis label based on MERFISH subset."""
    if merfish_subset == "manual":
        return "MERFISH (manually annotated, cortical)"
    else:
        return "MERFISH (all cells, L1-L6)"


def _xenium_ylabel(crop):
    """Build y-axis label based on crop setting."""
    if crop == "cropped":
        return "Xenium (L1-L6 only)"
    else:
        return "Xenium (all QC-pass cells)"


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Unified MERFISH vs Xenium proportion scatter plots.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--filter", choices=["neurons", "all"], default="all",
                        help="Cell type filter (default: all)")
    parser.add_argument("--level", choices=["supertype", "subclass", "both"],
                        default="supertype",
                        help="Taxonomy level (default: supertype)")
    parser.add_argument("--merfish-subset", choices=["all", "manual"],
                        default="all",
                        help="MERFISH reference subset (default: all)")
    parser.add_argument("--crop", choices=["uncropped", "cropped", "both"],
                        default="both",
                        help="Xenium crop setting (default: both)")
    parser.add_argument("--error-bars", action="store_true", default=False,
                        help="Show +/-SD error bars per point")
    parser.add_argument("--no-error-bars", dest="error_bars",
                        action="store_false")
    parser.add_argument("--weighted-dots", action="store_true", default=False,
                        help="Scale dot size by frequency, show weighted corr")
    parser.add_argument("--no-weighted-dots", dest="weighted_dots",
                        action="store_false")
    parser.add_argument("--max-labels", type=int, default=12,
                        help="Max cell type labels per panel (default: 12)")
    parser.add_argument("--output", type=str, default=None,
                        help="Override output filename (without directory)")
    args = parser.parse_args()

    out_tag = _make_output_tag(args)
    out_png = args.output or f"proportion_scatter_{out_tag}.png"
    out_png_path = os.path.join(PRESENTATION_DIR, out_png)
    out_csv_base = os.path.splitext(out_png)[0]

    need_std = args.error_bars or args.weighted_dots

    # Decide which levels to iterate over
    if args.level == "both":
        levels = ["subclass", "supertype"]
    else:
        levels = [args.level]

    # Decide which crops to iterate over
    if args.crop == "both":
        crops = ["uncropped", "cropped"]
    else:
        crops = [args.crop]

    # Build panel list: each panel is (level, crop) combination
    # If level==both, panels are side-by-side levels (crop is fixed).
    # If crop==both, panels are side-by-side crops (level is fixed).
    # If both are 'both', we'd get 4 panels — allowed but unusual.
    panels = [(lvl, crp) for lvl in levels for crp in crops]
    n_panels = len(panels)

    print(f"=== Proportion scatter: {n_panels} panel(s) ===")
    print(f"  filter={args.filter}, level={args.level}, "
          f"merfish={args.merfish_subset}, crop={args.crop}")
    print(f"  error_bars={args.error_bars}, weighted_dots={args.weighted_dots}")
    print()

    # --- Load data for each panel ---
    panel_data = []
    for lvl, crp in panels:
        print(f"--- Panel: {lvl} / {crp} ---")

        # MERFISH
        merfish = load_merfish_proportions(
            args.merfish_subset, lvl, return_std=need_std)

        # Xenium
        xenium = load_xenium_proportions(lvl, crp, return_std=need_std)

        # Merge
        merged = merfish.merge(xenium, on="celltype", how="inner")
        print(f"  Matched: {len(merged)} cell types")

        # Filter to neurons if requested
        if args.filter == "neurons":
            mask = merged["celltype"].apply(is_neuron)
            merged = merged[mask].copy()
            n_glut = sum(1 for ct in merged["celltype"]
                         if classify_celltype(ct)[1] == "Glutamatergic")
            n_gaba = sum(1 for ct in merged["celltype"]
                         if classify_celltype(ct)[1] == "GABAergic")
            print(f"  Neurons only: {len(merged)} types "
                  f"(Glut={n_glut}, GABA={n_gaba})")

        panel_data.append((lvl, crp, merged))
        print()

    # --- Plot ---
    fig_w = 11 * n_panels
    fig_h = 10
    fig, axes = plt.subplots(1, n_panels, figsize=(fig_w, fig_h), facecolor=BG,
                             squeeze=False)
    axes = axes.ravel()

    csv_paths = []
    for i, (lvl, crp, merged) in enumerate(panel_data):
        ax = axes[i]

        # Build title
        filter_tag = " (neurons only)" if args.filter == "neurons" else ""
        if n_panels > 1 and args.level == "both":
            title = f"{lvl.capitalize()} level{filter_tag}"
        elif n_panels > 1 and args.crop == "both":
            crop_tag = "uncropped" if crp == "uncropped" else "cropped to L1-L6"
            title = f"Xenium {crop_tag}{filter_tag}"
        else:
            crop_tag = "uncropped" if crp == "uncropped" else "cropped to L1-L6"
            title = f"{lvl.capitalize()} — Xenium {crop_tag}{filter_tag}"

        x_col = "merfish_mean"
        y_col = "xenium_mean"
        x_sd = "merfish_std" if need_std else None
        y_sd = "xenium_std" if need_std else None

        plot_scatter(
            ax, merged, x_col, y_col,
            title=title,
            xlabel=_merfish_xlabel(args.merfish_subset),
            ylabel=_xenium_ylabel(crp),
            x_sd_col=x_sd, y_sd_col=y_sd,
            use_error_bars=args.error_bars,
            use_weighted_dots=args.weighted_dots,
            max_labels=args.max_labels,
        )

        # Save per-panel CSV
        csv_name = f"{out_csv_base}_{lvl}_{crp}.csv"
        csv_path = os.path.join(PRESENTATION_DIR, csv_name)
        merged.to_csv(csv_path, index=False)
        csv_paths.append(csv_path)

    # Legend
    legend_handles = _build_legend_classes(args.filter)
    fig.legend(handles=legend_handles, loc="lower center",
               ncol=len(legend_handles), fontsize=14, frameon=False,
               labelcolor="white", bbox_to_anchor=(0.5, -0.01))

    plt.tight_layout(pad=2.0)
    plt.savefig(out_png_path, dpi=200, facecolor=BG, bbox_inches="tight")
    plt.close()

    print(f"Saved figure: {out_png_path}")
    for cp in csv_paths:
        print(f"Saved CSV:    {cp}")


if __name__ == "__main__":
    main()
