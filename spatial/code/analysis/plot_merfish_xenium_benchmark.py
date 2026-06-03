#!/usr/bin/env python3
"""
Benchmark cell typing: MERFISH vs Xenium (Corr Classifier vs Harmony).

Compares cell type proportions and depth distributions between:
  - SEA-AD MERFISH reference (manual depth annotations, cortical cells only)
  - Xenium correlation classifier (cortical cells only, pooled across all donors)
  - Xenium Harmony label transfer (cortical cells only, pooled across all donors)

All comparisons use cortical cells only (L1-L6) and pool across diagnosis
(no SCZ/Control split — this is purely a methods benchmark).

Produces:
  1. Subclass proportion scatter (MERFISH vs Xenium-Corr, MERFISH vs Xenium-Harmony)
  2. Supertype proportion scatter (same)
  3. Depth distribution violin plots per subclass (MERFISH vs Xenium-Corr vs Xenium-Harmony)
  4. Summary statistics CSV

Output: output/presentation/cell_typing_benchmark_*.png
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    H5AD_DIR, MERFISH_PATH, PRESENTATION_DIR, EXCLUDE_SAMPLES,
    CORTICAL_LAYERS, SUBCLASS_TO_CLASS, classify_celltype,
    load_merfish_cortical, load_cells, load_sample_adata, SUBCLASS_CONF_THRESH,
)

OUT_DIR = PRESENTATION_DIR
os.makedirs(OUT_DIR, exist_ok=True)

# Subclass ordering by class (for consistent plots)
GLUT_ORDER = ["L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP",
              "L6 CT", "L6 IT", "L6 IT Car3", "L6b"]
GABA_ORDER = ["Lamp5", "Lamp5 Lhx6", "Sncg", "Vip", "Pax6",
              "Chandelier", "Pvalb", "Sst", "Sst Chodl"]
NN_ORDER = ["Astrocyte", "Oligodendrocyte", "OPC",
            "Microglia-PVM", "Endothelial", "VLMC"]
SUBCLASS_ORDER = GLUT_ORDER + GABA_ORDER + NN_ORDER


# ── Data loading ──

def load_merfish_cortical_proportions(level="subclass"):
    """Load MERFISH cortical cells, compute per-donor proportions, return mean."""
    print("Loading MERFISH cortical cells (manual Layer annotation)...")
    df = load_merfish_cortical()  # uses config.py helper
    print(f"  {len(df):,} cortical cells with manual depth")

    col = level  # 'subclass' or 'supertype'
    records = []
    for donor, grp in df.groupby("donor"):
        total = len(grp)
        for ct, n in grp[col].value_counts().items():
            records.append({"donor": donor, "celltype": ct,
                            "proportion": n / total})
    prop_df = pd.DataFrame(records)
    mean_props = prop_df.groupby("celltype")["proportion"].mean().reset_index()
    mean_props.columns = ["celltype", "merfish_prop"]
    n_donors = df["donor"].nunique()
    print(f"  {n_donors} donors, {len(mean_props)} {level} types")
    return mean_props, df


def load_xenium_cortical_data():
    """Load all Xenium cortical cells via load_sample_adata() (standard QC).

    Returns a single DataFrame with all obs columns (including corr_* and
    harmony_* columns needed for the benchmark comparison).
    """
    print("Loading Xenium cortical cells...")
    h5ad_files = sorted(
        f for f in os.listdir(H5AD_DIR) if f.endswith("_annotated.h5ad")
    )

    dfs = []
    for fname in h5ad_files:
        sample_id = fname.replace("_annotated.h5ad", "")
        if sample_id in EXCLUDE_SAMPLES:
            continue
        adata = load_sample_adata(sample_id, cortical_only=True)
        dfs.append(adata.obs)
        print(f"  {sample_id}: {len(adata):,} cortical cells")

    df = pd.concat(dfs, ignore_index=True)
    print(f"  {len(df):,} cortical QC-pass cells from {len(dfs)} samples")
    return df


def compute_xenium_proportions(df, label_col, level_name):
    """Compute mean per-sample proportions for a given label column."""
    valid = df[label_col].notna()
    sub = df.loc[valid]

    records = []
    for sample_id, grp in sub.groupby("sample_id"):
        total = len(grp)
        for ct, n in grp[label_col].value_counts().items():
            records.append({"donor": sample_id, "celltype": str(ct),
                            "proportion": n / total})

    prop_df = pd.DataFrame(records)
    if len(prop_df) == 0:
        return pd.DataFrame(columns=["celltype", f"{level_name}_prop"])
    mean_props = prop_df.groupby("celltype")["proportion"].mean().reset_index()
    mean_props.columns = ["celltype", f"{level_name}_prop"]
    return mean_props


# ── Plotting ──

def plot_proportion_scatter(ax, merged, xcol, ycol, xlabel, ylabel, title,
                            max_labels=12):
    """Proportion scatter on log-log axes with correlation stats."""
    x = merged[xcol].values
    y = merged[ycol].values

    colors = [classify_celltype(ct)[0] for ct in merged["celltype"]]

    ax.scatter(x, y, c=colors, s=80, alpha=0.85, edgecolors="white",
               linewidths=0.5, zorder=5)

    # Label deviant points
    merged = merged.copy()
    merged["log_ratio"] = np.log2((y + 1e-7) / (x + 1e-7))
    merged["abs_lr"] = merged["log_ratio"].abs()
    top = merged.nlargest(max_labels, "abs_lr")
    for _, row in top.iterrows():
        ax.annotate(row["celltype"], (row[xcol], row[ycol]),
                    fontsize=9, alpha=0.8,
                    xytext=(5, 5), textcoords="offset points")

    # Diagonal
    pos_x = x[x > 0]
    pos_y = y[y > 0]
    if len(pos_x) > 0 and len(pos_y) > 0:
        lo = min(pos_x.min(), pos_y.min()) * 0.3
        hi = max(x.max(), y.max()) * 3
        ax.plot([lo, hi], [lo, hi], "--", color="#888888", lw=1.5, alpha=0.6, zorder=1)

    # Correlation
    valid = (x > 0) & (y > 0)
    if valid.sum() > 3:
        r, p = pearsonr(np.log10(x[valid]), np.log10(y[valid]))
        rho, _ = spearmanr(x[valid], y[valid])
        ax.text(0.04, 0.96,
                f"r = {r:.2f} (log)\n\u03c1 = {rho:.2f}\nn = {valid.sum()}",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=13, bbox=dict(boxstyle="round,pad=0.3",
                                       facecolor="white", alpha=0.7))

    ax.set_xscale("log")
    ax.set_yscale("log")

    def pct_fmt(val, pos):
        pct = val * 100
        if pct >= 1:
            return f"{pct:.0f}%"
        elif pct >= 0.1:
            return f"{pct:.1f}%"
        else:
            return f"{pct:.2f}%"

    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_formatter(mticker.FuncFormatter(pct_fmt))
        axis.set_minor_formatter(mticker.NullFormatter())

    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=16, fontweight="bold")
    ax.tick_params(labelsize=11)
    ax.grid(True, alpha=0.2)


def plot_depth_violins(merfish_df, xenium_df, subclasses=None):
    """
    Side-by-side violin plots of depth distribution per subclass.
    MERFISH (manual depth) vs Xenium Corr Classifier (predicted depth).
    """
    print("\n=== Depth Distribution Comparison ===")

    if subclasses is None:
        subclasses = [s for s in SUBCLASS_ORDER
                      if s in merfish_df["subclass"].values
                      and s in xenium_df["corr_subclass"].astype(str).values]

    fig, axes = plt.subplots(len(subclasses), 1, figsize=(14, 2.2 * len(subclasses)),
                              sharex=True)
    if len(subclasses) == 1:
        axes = [axes]

    for ax, sc in zip(axes, subclasses):
        # MERFISH
        m_depths = merfish_df.loc[merfish_df["subclass"] == sc, "depth"].values
        # Xenium corr
        x_mask = xenium_df["corr_subclass"].astype(str) == sc
        x_depths = xenium_df.loc[x_mask, "predicted_norm_depth"].dropna().values

        data = []
        labels = []
        colors_list = []
        if len(m_depths) > 0:
            data.append(m_depths)
            labels.append(f"MERFISH\n(n={len(m_depths):,})")
            colors_list.append("#66BB6A")
        if len(x_depths) > 0:
            data.append(x_depths)
            labels.append(f"Xenium Corr\n(n={len(x_depths):,})")
            colors_list.append("#FF9800")

        # Harmony if available
        if "harmony_subclass" in xenium_df.columns:
            h_mask = xenium_df["harmony_subclass"].astype(str) == sc
            h_depths = xenium_df.loc[h_mask, "predicted_norm_depth"].dropna().values
            if len(h_depths) > 0:
                data.append(h_depths)
                labels.append(f"Xenium Harmony\n(n={len(h_depths):,})")
                colors_list.append("#42A5F5")

        if len(data) == 0:
            ax.set_visible(False)
            continue

        vp = ax.violinplot(data, positions=range(len(data)), showmedians=True,
                           showextrema=False)
        for i, body in enumerate(vp["bodies"]):
            body.set_facecolor(colors_list[i])
            body.set_alpha(0.7)
        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(2)

        ax.set_yticks(range(len(data)))
        ax.set_yticklabels(labels, fontsize=10)
        ax.set_xlim(-0.1, 1.1)

        # Class color for subclass label
        cls = SUBCLASS_TO_CLASS.get(sc, "Unknown")
        bg = {"Glutamatergic": "#FFF3E0", "GABAergic": "#E8F5E9",
              "Non-neuronal": "#E3F2FD"}.get(cls, "white")
        ax.set_facecolor(bg)

        ax.set_title(sc, fontsize=13, fontweight="bold", loc="left", pad=4)
        ax.grid(axis="x", alpha=0.3)

    axes[-1].set_xlabel("Normalized Cortical Depth (0=pia, 1=WM)", fontsize=14)
    fig.suptitle("Depth Distribution by Subclass: MERFISH vs Xenium",
                 fontsize=18, fontweight="bold", y=1.01)
    plt.tight_layout()
    return fig


def plot_depth_summary(merfish_df, xenium_df):
    """
    Compact depth comparison: median depth per subclass, MERFISH vs Xenium.
    One scatter point per subclass, colored by class.
    """
    print("\n=== Median Depth Scatter ===")

    records = []
    for sc in SUBCLASS_ORDER:
        m_mask = merfish_df["subclass"] == sc
        x_mask = xenium_df["corr_subclass"].astype(str) == sc

        m_depth = merfish_df.loc[m_mask, "depth"].median() if m_mask.sum() > 50 else np.nan
        x_depth = xenium_df.loc[x_mask, "predicted_norm_depth"].median() if x_mask.sum() > 50 else np.nan

        h_depth = np.nan
        if "harmony_subclass" in xenium_df.columns:
            h_mask = xenium_df["harmony_subclass"].astype(str) == sc
            if h_mask.sum() > 50:
                h_depth = xenium_df.loc[h_mask, "predicted_norm_depth"].median()

        records.append({
            "subclass": sc,
            "merfish_depth": m_depth,
            "xenium_corr_depth": x_depth,
            "xenium_harmony_depth": h_depth,
            "class": SUBCLASS_TO_CLASS.get(sc, "Unknown"),
        })

    df = pd.DataFrame(records)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    class_colors = {"Glutamatergic": "#E65100", "GABAergic": "#2E7D32",
                    "Non-neuronal": "#1565C0"}

    # Panel 1: MERFISH vs Xenium Corr
    ax = axes[0]
    valid = df["merfish_depth"].notna() & df["xenium_corr_depth"].notna()
    sub = df[valid]
    colors = [class_colors.get(c, "gray") for c in sub["class"]]
    ax.scatter(sub["merfish_depth"], sub["xenium_corr_depth"],
               c=colors, s=100, alpha=0.85, edgecolors="white", lw=0.8, zorder=5)
    for _, row in sub.iterrows():
        ax.annotate(row["subclass"], (row["merfish_depth"], row["xenium_corr_depth"]),
                    fontsize=8, alpha=0.7, xytext=(4, 4), textcoords="offset points")
    ax.plot([0, 1], [0, 1], "--", color="gray", lw=1.5, alpha=0.5)
    if valid.sum() > 3:
        r, _ = pearsonr(sub["merfish_depth"], sub["xenium_corr_depth"])
        ax.text(0.04, 0.96, f"r = {r:.3f}", transform=ax.transAxes,
                fontsize=14, va="top", bbox=dict(facecolor="white", alpha=0.7))
    ax.set_xlabel("MERFISH Median Depth (manual annotation)", fontsize=13)
    ax.set_ylabel("Xenium Corr Classifier Median Depth (predicted)", fontsize=13)
    ax.set_title("Corr Classifier vs MERFISH", fontsize=15, fontweight="bold")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect("equal")
    ax.grid(alpha=0.2)

    # Panel 2: MERFISH vs Xenium Harmony
    ax = axes[1]
    valid2 = df["merfish_depth"].notna() & df["xenium_harmony_depth"].notna()
    sub2 = df[valid2]
    if len(sub2) > 0:
        colors2 = [class_colors.get(c, "gray") for c in sub2["class"]]
        ax.scatter(sub2["merfish_depth"], sub2["xenium_harmony_depth"],
                   c=colors2, s=100, alpha=0.85, edgecolors="white", lw=0.8, zorder=5)
        for _, row in sub2.iterrows():
            ax.annotate(row["subclass"],
                        (row["merfish_depth"], row["xenium_harmony_depth"]),
                        fontsize=8, alpha=0.7, xytext=(4, 4), textcoords="offset points")
        ax.plot([0, 1], [0, 1], "--", color="gray", lw=1.5, alpha=0.5)
        if valid2.sum() > 3:
            r2, _ = pearsonr(sub2["merfish_depth"], sub2["xenium_harmony_depth"])
            ax.text(0.04, 0.96, f"r = {r2:.3f}", transform=ax.transAxes,
                    fontsize=14, va="top", bbox=dict(facecolor="white", alpha=0.7))
    ax.set_xlabel("MERFISH Median Depth (manual annotation)", fontsize=13)
    ax.set_ylabel("Xenium Harmony Median Depth (predicted)", fontsize=13)
    ax.set_title("Harmony vs MERFISH", fontsize=15, fontweight="bold")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect("equal")
    ax.grid(alpha=0.2)

    # Legend
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c,
                       markersize=10, label=cls)
               for cls, c in class_colors.items()]
    axes[0].legend(handles=handles, fontsize=11, loc="lower right")

    fig.suptitle("Median Cortical Depth per Subclass: MERFISH Reference vs Xenium Methods",
                 fontsize=16, fontweight="bold", y=1.02)
    plt.tight_layout()
    return fig, df


def main():
    # ── Load data ──
    merfish_sub_props, merfish_df = load_merfish_cortical_proportions("subclass")
    merfish_sup_props, _ = load_merfish_cortical_proportions("supertype")
    xenium_df = load_xenium_cortical_data()

    # Xenium proportions: corr classifier
    xenium_corr_sub = compute_xenium_proportions(xenium_df, "corr_subclass", "corr")
    xenium_corr_sup = compute_xenium_proportions(xenium_df, "corr_supertype", "corr")

    # Xenium proportions: harmony
    has_harmony = "harmony_subclass" in xenium_df.columns
    if has_harmony:
        xenium_harm_sub = compute_xenium_proportions(xenium_df, "harmony_subclass", "harmony")
        xenium_harm_sup = compute_xenium_proportions(xenium_df, "harmony_supertype", "harmony")

    # ── Figure 1: Subclass proportion scatter ──
    print("\n=== Subclass Proportion Scatter ===")
    fig1, axes1 = plt.subplots(1, 2 if has_harmony else 1,
                                figsize=(18 if has_harmony else 9, 8))
    if not has_harmony:
        axes1 = [axes1]

    merged_corr_sub = merfish_sub_props.merge(xenium_corr_sub, on="celltype", how="inner")
    print(f"  Corr: {len(merged_corr_sub)} subclasses matched")
    plot_proportion_scatter(axes1[0], merged_corr_sub, "merfish_prop", "corr_prop",
                            "MERFISH proportion", "Xenium Corr Classifier proportion",
                            "Subclass: Corr Classifier vs MERFISH")

    if has_harmony:
        merged_harm_sub = merfish_sub_props.merge(xenium_harm_sub, on="celltype", how="inner")
        print(f"  Harmony: {len(merged_harm_sub)} subclasses matched")
        plot_proportion_scatter(axes1[1], merged_harm_sub, "merfish_prop", "harmony_prop",
                                "MERFISH proportion", "Xenium Harmony proportion",
                                "Subclass: Harmony vs MERFISH")

    fig1.suptitle("Cortical Cell Type Proportions: MERFISH vs Xenium",
                  fontsize=18, fontweight="bold", y=1.02)
    plt.tight_layout()
    path1 = os.path.join(OUT_DIR, "cell_typing_benchmark_subclass_proportions.png")
    fig1.savefig(path1, dpi=200, bbox_inches="tight")
    plt.close(fig1)
    print(f"  Saved: {path1}")

    # ── Figure 2: Supertype proportion scatter ──
    print("\n=== Supertype Proportion Scatter ===")
    fig2, axes2 = plt.subplots(1, 2 if has_harmony else 1,
                                figsize=(18 if has_harmony else 9, 8))
    if not has_harmony:
        axes2 = [axes2]

    merged_corr_sup = merfish_sup_props.merge(xenium_corr_sup, on="celltype", how="inner")
    print(f"  Corr: {len(merged_corr_sup)} supertypes matched")
    plot_proportion_scatter(axes2[0], merged_corr_sup, "merfish_prop", "corr_prop",
                            "MERFISH proportion", "Xenium Corr Classifier proportion",
                            "Supertype: Corr Classifier vs MERFISH")

    if has_harmony:
        merged_harm_sup = merfish_sup_props.merge(xenium_harm_sup, on="celltype", how="inner")
        print(f"  Harmony: {len(merged_harm_sup)} supertypes matched")
        plot_proportion_scatter(axes2[1], merged_harm_sup, "merfish_prop", "harmony_prop",
                                "MERFISH proportion", "Xenium Harmony proportion",
                                "Supertype: Harmony vs MERFISH")

    fig2.suptitle("Cortical Supertype Proportions: MERFISH vs Xenium",
                  fontsize=18, fontweight="bold", y=1.02)
    plt.tight_layout()
    path2 = os.path.join(OUT_DIR, "cell_typing_benchmark_supertype_proportions.png")
    fig2.savefig(path2, dpi=200, bbox_inches="tight")
    plt.close(fig2)
    print(f"  Saved: {path2}")

    # ── Figure 3: Median depth scatter ──
    fig3, depth_df = plot_depth_summary(merfish_df, xenium_df)
    path3 = os.path.join(OUT_DIR, "cell_typing_benchmark_depth_scatter.png")
    fig3.savefig(path3, dpi=200, bbox_inches="tight")
    plt.close(fig3)
    print(f"  Saved: {path3}")

    # ── Figure 4: Depth violin plots (excitatory types only — clearest depth signal) ──
    exc_types = [s for s in GLUT_ORDER
                 if s in merfish_df["subclass"].values
                 and s in xenium_df["corr_subclass"].astype(str).values]
    fig4 = plot_depth_violins(merfish_df, xenium_df, exc_types)
    path4 = os.path.join(OUT_DIR, "cell_typing_benchmark_depth_violins_excitatory.png")
    fig4.savefig(path4, dpi=200, bbox_inches="tight")
    plt.close(fig4)
    print(f"  Saved: {path4}")

    # ── Save summary CSV ──
    summary_records = []
    for _, row in merged_corr_sub.iterrows():
        rec = {"subclass": row["celltype"],
               "merfish_prop": row["merfish_prop"],
               "xenium_corr_prop": row["corr_prop"]}
        if has_harmony:
            harm_row = merged_harm_sub[merged_harm_sub["celltype"] == row["celltype"]]
            rec["xenium_harmony_prop"] = harm_row["harmony_prop"].values[0] if len(harm_row) > 0 else np.nan
        # Add depth info
        depth_row = depth_df[depth_df["subclass"] == row["celltype"]]
        if len(depth_row) > 0:
            rec["merfish_median_depth"] = depth_row["merfish_depth"].values[0]
            rec["xenium_corr_median_depth"] = depth_row["xenium_corr_depth"].values[0]
            if has_harmony:
                rec["xenium_harmony_median_depth"] = depth_row["xenium_harmony_depth"].values[0]
        summary_records.append(rec)

    summary = pd.DataFrame(summary_records)
    csv_path = os.path.join(OUT_DIR, "cell_typing_benchmark_summary.csv")
    summary.to_csv(csv_path, index=False)
    print(f"\n  Saved: {csv_path}")

    # ── Print summary stats ──
    print("\n" + "=" * 70)
    print("BENCHMARK SUMMARY (cortical cells only)")
    print("=" * 70)

    valid_cs = merged_corr_sub.dropna()
    r_corr_sub, _ = pearsonr(np.log10(valid_cs["merfish_prop"]),
                              np.log10(valid_cs["corr_prop"]))
    rho_corr_sub, _ = spearmanr(valid_cs["merfish_prop"], valid_cs["corr_prop"])
    print(f"  Subclass proportions (Corr vs MERFISH): r={r_corr_sub:.3f}, rho={rho_corr_sub:.3f}")

    if has_harmony:
        valid_hs = merged_harm_sub.dropna()
        if len(valid_hs) > 3:
            r_harm_sub, _ = pearsonr(np.log10(valid_hs["merfish_prop"] + 1e-8),
                                      np.log10(valid_hs["harmony_prop"] + 1e-8))
            rho_harm_sub, _ = spearmanr(valid_hs["merfish_prop"], valid_hs["harmony_prop"])
            print(f"  Subclass proportions (Harmony vs MERFISH): r={r_harm_sub:.3f}, rho={rho_harm_sub:.3f}")

    # Depth correlation
    dv = depth_df.dropna(subset=["merfish_depth", "xenium_corr_depth"])
    if len(dv) > 3:
        r_depth, _ = pearsonr(dv["merfish_depth"], dv["xenium_corr_depth"])
        print(f"  Median depth per subclass (Corr vs MERFISH): r={r_depth:.3f}")

    dv2 = depth_df.dropna(subset=["merfish_depth", "xenium_harmony_depth"])
    if len(dv2) > 3:
        r_depth2, _ = pearsonr(dv2["merfish_depth"], dv2["xenium_harmony_depth"])
        print(f"  Median depth per subclass (Harmony vs MERFISH): r={r_depth2:.3f}")


if __name__ == "__main__":
    main()
