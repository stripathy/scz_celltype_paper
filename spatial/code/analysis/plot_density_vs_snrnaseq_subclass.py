#!/usr/bin/env python3
"""
Scatter plot: snRNAseq meta-analysis betas vs Xenium density log2FC at SUBCLASS level.

Two-panel layout: Neuronal (left) | Non-neuronal (right).
Uses Nicole's direct subclass-level crumblr results (not aggregated from supertype).
For density, aggregates supertype-level per-sample densities to subclass.

Input:
  output/density_analysis/density_per_sample_supertype.csv
  data/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts_subclass.csv (neuronal)
  data/nicole_scz_snrnaseq_betas/final_results_crumblr_7_nonN_cohorts_subclass.csv (non-neuronal)

Output:
  output/density_analysis/snrnaseq_vs_density_subclass.png
  output/density_analysis/snrnaseq_vs_density_subclass.csv
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import pearsonr, ttest_ind
from adjustText import adjust_text

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BASE_DIR, EXCLUDE_SAMPLES, INFER_CLASS_COLORS as CLASS_COLORS,
    infer_class, BG_COLOR,
)

DENSITY_RAW = os.path.join(BASE_DIR, "output", "density_analysis",
                            "density_per_sample_supertype.csv")
NICOLE_NEURONAL_PATH = os.path.join(BASE_DIR, "data", "nicole_scz_snrnaseq_betas",
                                     "final_results_crumblr_7_cohorts_subclass.csv")
NICOLE_NONNEURONAL_PATH = os.path.join(BASE_DIR, "data", "nicole_scz_snrnaseq_betas",
                                        "final_results_crumblr_7_nonN_cohorts_subclass.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "density_analysis")

BG = BG_COLOR


def supertype_to_subclass(st):
    """Map supertype name to subclass by stripping trailing _N."""
    parts = st.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return st


# Xenium density uses abbreviated subclass names; Nicole uses full names.
SUBCLASS_NAME_MAP = {
    "Astro": "Astrocyte",
    "Endo": "Endothelial",
    "Micro-PVM": "Microglia-PVM",
    "Oligo": "Oligodendrocyte",
}


def plot_panel(ax, panel_df, title, show_ylabel=True):
    """Plot a single density scatter panel — labels all points (subclass = few)."""
    ax.set_facecolor(BG)

    for _, row in panel_df.iterrows():
        color = CLASS_COLORS.get(row["class"], '#999999')
        # Horizontal (snRNAseq SE)
        if pd.notna(row.get("se_snrnaseq")):
            ax.plot([row["beta_snrnaseq"] - row["se_snrnaseq"],
                     row["beta_snrnaseq"] + row["se_snrnaseq"]],
                    [row["logFC_density"], row["logFC_density"]],
                    color=color, alpha=0.25, linewidth=1.5, zorder=3)
        # Vertical (density SE)
        if pd.notna(row.get("se_density")) and pd.notna(row["logFC_density"]):
            ax.plot([row["beta_snrnaseq"], row["beta_snrnaseq"]],
                    [row["logFC_density"] - row["se_density"],
                     row["logFC_density"] + row["se_density"]],
                    color=color, alpha=0.25, linewidth=1.5, zorder=3)

    for cls in ["Glut", "GABA", "NN", "Other"]:
        mask = panel_df["class"] == cls
        if mask.sum() == 0:
            continue
        sub = panel_df[mask]
        ax.scatter(sub["beta_snrnaseq"], sub["logFC_density"],
                   c=CLASS_COLORS[cls], s=120, alpha=0.8,
                   edgecolors="white", linewidth=0.8, zorder=5,
                   label=f"{cls} (n={mask.sum()})")

    valid = panel_df["beta_snrnaseq"].notna() & panel_df["logFC_density"].notna()
    m = panel_df[valid]
    if len(m) > 2:
        z = np.polyfit(m["beta_snrnaseq"], m["logFC_density"], 1)
        lim_x = max(abs(m["beta_snrnaseq"]).max(), 0.1) * 1.5
        x_line = np.linspace(-lim_x, lim_x, 100)
        ax.plot(x_line, np.polyval(z, x_line), color="#888888", alpha=0.5,
                linewidth=1.5, linestyle="-", zorder=1)

    ax.axhline(0, color="#555555", alpha=0.4, linewidth=0.8, zorder=1)
    ax.axvline(0, color="#555555", alpha=0.4, linewidth=0.8, zorder=1)

    # Label ALL points (subclass level = few)
    texts = []
    for _, row in panel_df.iterrows():
        txt = ax.text(row["beta_snrnaseq"], row["logFC_density"],
                      f"  {row['subclass']}",
                      fontsize=15, fontweight="bold",
                      color="white", alpha=0.95, zorder=10)
        texts.append(txt)

    if texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#aaaaaa",
                                   alpha=0.5, linewidth=0.7),
                    force_text=(0.9, 0.9), force_points=(0.6, 0.6),
                    expand_text=(1.3, 1.3), expand_points=(1.5, 1.5))

    if len(m) > 2:
        r_val, p_val = pearsonr(m["beta_snrnaseq"], m["logFC_density"])
        ax.text(0.97, 0.04,
                f"r = {r_val:.2f} (p = {p_val:.1e})\nn = {len(m)}",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=16, color="#dddddd",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#222222",
                          edgecolor="#555555", alpha=0.9))

    ax.set_xlabel("snRNAseq meta-analysis beta", fontsize=20, color="white")
    if show_ylabel:
        ax.set_ylabel("Xenium density logFC (cells/mm²)", fontsize=20, color="white")
    ax.set_title(title, fontsize=22, fontweight="bold", color="white", pad=10)

    ax.tick_params(colors="white", labelsize=16)
    for spine in ax.spines.values():
        spine.set_color("#555555")
    ax.grid(True, alpha=0.15, color="#555555")


def main():
    # ── Load per-sample supertype densities ──
    raw = pd.read_csv(DENSITY_RAW)
    raw = raw[raw["region"] == "cortical"]
    raw = raw[~raw["sample_id"].isin(EXCLUDE_SAMPLES)]
    raw["subclass"] = raw["supertype"].apply(supertype_to_subclass)
    print(f"Raw data: {len(raw)} rows, {raw['sample_id'].nunique()} samples, "
          f"{raw['subclass'].nunique()} subclasses")

    # Aggregate density to subclass level per sample
    sub_density = (raw.groupby(["sample_id", "diagnosis", "subclass"])
                   .agg(density_per_mm2=("density_per_mm2", "sum"),
                        count=("count", "sum"),
                        area_mm2=("area_mm2", "first"))
                   .reset_index())

    # Compute logFC per subclass
    results = []
    for sc in sorted(sub_density["subclass"].unique()):
        sc_data = sub_density[sub_density["subclass"] == sc]
        ctrl = sc_data[sc_data["diagnosis"] == "Control"]["density_per_mm2"]
        scz = sc_data[sc_data["diagnosis"] == "SCZ"]["density_per_mm2"]

        if len(ctrl) < 3 or len(scz) < 3:
            continue

        ctrl_mean = ctrl.mean()
        scz_mean = scz.mean()
        logFC = np.log(scz_mean / ctrl_mean) if ctrl_mean > 0 else np.nan

        _, pval = ttest_ind(np.log1p(ctrl), np.log1p(scz), equal_var=False)
        log_vals = np.log1p(sc_data["density_per_mm2"])
        se = log_vals.std() / np.sqrt(len(log_vals))

        results.append({
            "subclass": sc,
            "ctrl_mean_density": ctrl_mean,
            "scz_mean_density": scz_mean,
            "logFC_density": logFC,
            "pval_density": pval,
            "se_density": se,
            "n_ctrl": len(ctrl),
            "n_scz": len(scz),
            "class": infer_class(sc),
        })

    density_df = pd.DataFrame(results)
    # Normalize abbreviated subclass names to match Nicole's naming
    density_df["subclass"] = density_df["subclass"].replace(SUBCLASS_NAME_MAP)
    # Filter out SEAAD variants
    density_df = density_df[~density_df["subclass"].str.contains("SEAAD", na=False)]
    # Re-aggregate if multiple rows map to the same normalized name
    # (shouldn't happen but be safe)
    print(f"Density results: {len(density_df)} subclasses")

    # ── Load Nicole's direct subclass-level results ──
    dfs = []
    for path, stratum in [(NICOLE_NEURONAL_PATH, "neuronal"),
                           (NICOLE_NONNEURONAL_PATH, "non-neuronal")]:
        if not os.path.exists(path):
            print(f"  WARNING: {path} not found")
            continue
        df = pd.read_csv(path)
        df = df[~df["CellType"].str.contains("SEAAD", na=False)]
        df["analysis_stratum"] = stratum
        dfs.append(df)
    nicole = pd.concat(dfs, ignore_index=True)
    nicole = nicole.rename(columns={
        "CellType": "subclass", "estimate": "beta_snrnaseq",
        "pval": "pval_snrnaseq", "padj": "padj_snrnaseq",
    })
    nicole = nicole.rename(columns={"se": "se_snrnaseq"})
    print(f"Nicole subclass data: {len(nicole)} subclasses "
          f"({sum(nicole['analysis_stratum']=='neuronal')} neuronal, "
          f"{sum(nicole['analysis_stratum']=='non-neuronal')} non-neuronal)")

    # ── Merge ──
    merged = pd.merge(
        nicole[["subclass", "beta_snrnaseq", "se_snrnaseq", "pval_snrnaseq", "padj_snrnaseq"]],
        density_df, on="subclass", how="inner"
    )
    print(f"Merged: {len(merged)} shared subclasses")

    # Split
    neuronal = merged[merged["class"].isin(["Glut", "GABA"])].copy()
    nonneuronal = merged[merged["class"] == "NN"].copy()
    print(f"Neuronal: {len(neuronal)}, Non-neuronal: {len(nonneuronal)}")

    # Correlations
    valid = merged["logFC_density"].notna() & merged["beta_snrnaseq"].notna()
    m = merged[valid]
    r_all, p_all = pearsonr(m["beta_snrnaseq"], m["logFC_density"])
    r_neur, p_neur = pearsonr(neuronal["beta_snrnaseq"], neuronal["logFC_density"]) if len(neuronal) > 3 else (np.nan, np.nan)
    r_nn, p_nn = pearsonr(nonneuronal["beta_snrnaseq"], nonneuronal["logFC_density"]) if len(nonneuronal) > 2 else (np.nan, np.nan)
    print(f"All: r={r_all:.3f} (p={p_all:.1e}, n={len(m)})")
    print(f"Neuronal: r={r_neur:.3f} (p={p_neur:.1e}, n={len(neuronal)})")
    print(f"Non-neuronal: r={r_nn:.3f} (p={p_nn:.1e}, n={len(nonneuronal)})")

    # ── Two-panel figure ──
    fig, (ax_n, ax_nn) = plt.subplots(1, 2, figsize=(22, 10), facecolor=BG)

    plot_panel(ax_n, neuronal, title=f'Neuronal subclasses (n={len(neuronal)})',
              show_ylabel=True)
    plot_panel(ax_nn, nonneuronal, title=f'Non-neuronal subclasses (n={len(nonneuronal)})',
              show_ylabel=False)

    # Legend on left panel
    legend_elements = [
        Line2D([0], [0], marker="o", color=BG, markerfacecolor=CLASS_COLORS["Glut"],
               markersize=12, label=f"Glutamatergic (n={(neuronal['class']=='Glut').sum()})", linewidth=0),
        Line2D([0], [0], marker="o", color=BG, markerfacecolor=CLASS_COLORS["GABA"],
               markersize=12, label=f"GABAergic (n={(neuronal['class']=='GABA').sum()})", linewidth=0),
    ]
    leg = ax_n.legend(handles=legend_elements, loc="upper left", fontsize=15,
                      frameon=True, fancybox=True, framealpha=0.85,
                      edgecolor="#555555", labelcolor="white")
    leg.get_frame().set_facecolor("#222222")

    fig.suptitle("SCZ Subclass Density Effects: snRNAseq vs Xenium",
                 fontsize=26, fontweight="bold", color="white", y=1.02)

    plt.tight_layout(pad=2.0)
    out_png = os.path.join(OUTPUT_DIR, "snrnaseq_vs_density_subclass.png")
    plt.savefig(out_png, dpi=200, facecolor=BG, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {out_png}")

    # Save CSV
    out_csv = os.path.join(OUTPUT_DIR, "snrnaseq_vs_density_subclass.csv")
    merged.to_csv(out_csv, index=False)
    print(f"Saved: {out_csv}")

    # Print all subclasses
    print(f"\nAll subclasses (sorted by snRNAseq beta):")
    for _, row in merged.sort_values("beta_snrnaseq").iterrows():
        concordant = "Y" if (row["beta_snrnaseq"] * row["logFC_density"]) > 0 else "N"
        print(f"  {row['subclass']:20s} snRNA={row['beta_snrnaseq']:+.3f}  "
              f"density={row['logFC_density']:+.3f}  {concordant}")


if __name__ == "__main__":
    main()
