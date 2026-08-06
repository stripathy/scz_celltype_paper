#!/usr/bin/env python
"""
generate_markdown_figures.py — Generate publication figures for the
rnaseq_specificity_enrichment_results.md document.

Creates:
  1. SEA-AD Manhattan scatter (webapp-inspired style)
  2. SST dual-hit figure (GWAS enrichment + composition change)
  3. SST composition 2-panel (focused on the key correlation)
  4. SST deep dive panel using RBH combined taxonomy enrichment + depth

All figures use cowplot-style aesthetics: white background, no grid,
half-open axes, large legible text.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from scipy import stats

from scz_celltype_enrichment.config import (
    TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR, FIGURE_DPI,
)
from scz_celltype_enrichment.utils import (
    classify_supertypes_by_class, extract_subclass,
)


# ── Cowplot-style global settings ──────────────────────────────────────────
def cowplot_style():
    mpl.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 1.5,
        "font.family": "sans-serif",
        "font.size": 14,
    })


CLASS_COLORS = {
    "GABAergic": "#D55E00",
    "Glutamatergic": "#0072B2",
    "Non-neuronal": "#009E73",
}

GABA_SUBCLASSES = [
    "Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb",
    "Sncg", "Sst", "Sst Chodl", "Vip",
]
GLUT_SUBCLASSES = [
    "L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP",
    "L6 CT", "L6 IT", "L6 IT Car3", "L6b",
]


def _seaad_sort_key(name):
    sc = extract_subclass(name)
    for i, p in enumerate(GABA_SUBCLASSES):
        if sc == p:
            return (0, i, name)
    for i, p in enumerate(GLUT_SUBCLASSES):
        if sc == p:
            return (1, i, name)
    return (2, 0, name)


# ══════════════════════════════════════════════════════════════════════════
# Figure 1: SEA-AD Manhattan Scatter (webapp-inspired)
# ══════════════════════════════════════════════════════════════════════════
def figure1_seaad_manhattan(enrich, output_path):
    """Manhattan-style scatter of all 137 SEA-AD supertypes."""
    cowplot_style()

    LS = 20; TS = 22; TKS = 15; LGS = 13; ANNOT = 10

    df = enrich.copy()
    df["logp"] = -np.log10(np.clip(df["p_value"], 1e-300, 1))
    cls_map = classify_supertypes_by_class(df["supertype"].tolist())
    df["cell_class"] = [cls_map.get(ct, "Other") for ct in df["supertype"]]
    df["subclass"] = df["supertype"].apply(extract_subclass)

    # Sort by subclass order
    df = df.sort_values("supertype", key=lambda s: s.map(_seaad_sort_key))
    df = df.reset_index(drop=True)
    df["x"] = range(len(df))

    fig, ax = plt.subplots(figsize=(22, 9))

    # Alternating subclass background bands
    subclasses = df["subclass"].values
    band_starts = [0]
    for i in range(1, len(subclasses)):
        if subclasses[i] != subclasses[i - 1]:
            band_starts.append(i)
    band_starts.append(len(df))
    for j in range(len(band_starts) - 1):
        if j % 2 == 0:
            ax.axvspan(band_starts[j] - 0.5, band_starts[j + 1] - 0.5,
                        color="#f5f5f5", zorder=0)

    # Plot points
    for cls in ["Non-neuronal", "Glutamatergic", "GABAergic"]:
        mask = df["cell_class"] == cls
        fdr_sig = df.loc[mask, "p_fdr"] < 0.05
        bonf_sig = df.loc[mask, "p_bonferroni"] < 0.05

        # Non-significant
        ns = mask & ~(df["p_fdr"] < 0.05)
        ax.scatter(df.loc[ns, "x"], df.loc[ns, "logp"],
                   c=CLASS_COLORS.get(cls, "gray"), s=30, alpha=0.35,
                   edgecolors="none", zorder=2)

        # FDR significant (open diamonds)
        fdr_only = mask & (df["p_fdr"] < 0.05) & ~(df["p_bonferroni"] < 0.05)
        ax.scatter(df.loc[fdr_only, "x"], df.loc[fdr_only, "logp"],
                   c=CLASS_COLORS.get(cls, "gray"), s=70, alpha=0.8,
                   marker="D", edgecolors="black", linewidth=0.5, zorder=3)

        # Bonferroni significant (filled diamonds, larger)
        bonf = mask & (df["p_bonferroni"] < 0.05)
        ax.scatter(df.loc[bonf, "x"], df.loc[bonf, "logp"],
                   c=CLASS_COLORS.get(cls, "gray"), s=120, alpha=0.9,
                   marker="D", edgecolors="black", linewidth=0.8, zorder=4)

    # Threshold lines
    bonf_thresh = -np.log10(0.05 / len(df))
    fdr_sig_types = df[df["p_fdr"] < 0.05]
    fdr_thresh = -np.log10(fdr_sig_types["p_value"].max()) if len(fdr_sig_types) > 0 else 2
    ax.axhline(bonf_thresh, color="#e74c3c", linestyle="--", linewidth=1.2,
               alpha=0.6, label=f"Bonferroni (p < {0.05/len(df):.1e})")
    ax.axhline(fdr_thresh, color="#f39c12", linestyle=":", linewidth=1.2,
               alpha=0.6, label="FDR < 0.05")

    # Label top 20 types
    top = df.nsmallest(20, "p_value")
    for _, row in top.iterrows():
        y_off = 5 if row["logp"] < 9 else -12
        ax.annotate(
            row["supertype"], (row["x"], row["logp"]),
            fontsize=ANNOT, fontweight="bold",
            color=CLASS_COLORS.get(row["cell_class"], "gray"),
            rotation=45, ha="left", va="bottom",
            xytext=(3, y_off), textcoords="offset points",
        )

    # Subclass labels along bottom
    for j in range(len(band_starts) - 1):
        mid = (band_starts[j] + band_starts[j + 1]) / 2 - 0.5
        sc_name = subclasses[band_starts[j]]
        if band_starts[j + 1] - band_starts[j] >= 3:
            ax.text(mid, -0.7, sc_name, ha="center", va="top",
                    fontsize=9, rotation=45, color="#666")

    # Legend
    legend_els = [
        Line2D([0], [0], marker="D", color="w", markerfacecolor="#D55E00",
               markersize=10, markeredgecolor="black", linewidth=0, label="GABAergic"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor="#0072B2",
               markersize=10, markeredgecolor="black", linewidth=0, label="Glutamatergic"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor="#009E73",
               markersize=10, markeredgecolor="black", linewidth=0, label="Non-neuronal"),
        Line2D([0], [0], color="#e74c3c", linestyle="--", linewidth=1.2, label="Bonferroni"),
        Line2D([0], [0], color="#f39c12", linestyle=":", linewidth=1.2, label="FDR 0.05"),
    ]
    ax.legend(handles=legend_els, fontsize=LGS, frameon=False, loc="upper right",
              ncol=2)

    n_fdr = (df["p_fdr"] < 0.05).sum()
    n_bonf = (df["p_bonferroni"] < 0.05).sum()
    ax.set_title(
        f"SCZ GWAS Enrichment Across 137 SEA-AD Cell Types\n"
        f"{n_bonf} Bonferroni-significant, {n_fdr} FDR-significant (41 GABAergic, 3 Glutamatergic)",
        fontsize=TS, fontweight="bold",
    )
    ax.set_xlabel("Cell type (grouped by subclass)", fontsize=LS)
    ax.set_ylabel(r"$-\log_{10}(p)$", fontsize=LS)
    ax.tick_params(labelsize=TKS)
    ax.set_xticks([])
    ax.set_xlim(-2, len(df) + 2)
    ax.set_ylim(-0.5, max(df["logp"]) * 1.15)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ══════════════════════════════════════════════════════════════════════════
# Figure 2: SST Dual-Hit (GWAS + Composition Side by Side)
# ══════════════════════════════════════════════════════════════════════════
def figure2_sst_dual_hit(comp_df, output_path):
    """Paired bar chart: GWAS enrichment + composition change per SST subtype."""
    cowplot_style()

    LS = 18; TS = 20; TKS = 14; ANNOT = 11

    sst = comp_df[comp_df["is_sst"]].copy()
    sst = sst.sort_values("gwas_logp", ascending=True)  # ascending for horizontal bars

    n = len(sst)
    y_pos = np.arange(n)

    fig, axes = plt.subplots(1, 3, figsize=(22, max(8, n * 0.5)),
                              gridspec_kw={"width_ratios": [2, 2, 1.5], "wspace": 0.35})

    # --- Panel A: GWAS enrichment ---
    ax = axes[0]
    colors = []
    for _, row in sst.iterrows():
        if row["p_fdr"] < 0.001:
            colors.append("#c0392b")  # dark red
        elif row["p_fdr"] < 0.05:
            colors.append("#D55E00")  # orange-red
        else:
            colors.append("#bdc3c7")  # gray
    ax.barh(y_pos, sst["gwas_logp"], color=colors, edgecolor="black",
            linewidth=0.5, height=0.7)

    # FDR and Bonferroni lines
    bonf_thresh = -np.log10(0.05 / 137)
    ax.axvline(bonf_thresh, color="#e74c3c", linestyle="--", linewidth=1.2,
               alpha=0.6, label="Bonferroni")
    fdr_line = -np.log10(0.05)
    ax.axvline(fdr_line, color="#f39c12", linestyle=":", linewidth=1.2,
               alpha=0.6, label="Nominal")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(sst["supertype"], fontsize=13, fontweight="bold")
    ax.set_xlabel(r"GWAS enrichment $-\log_{10}(p)$", fontsize=LS)
    ax.set_title("A. GWAS cell-type enrichment", fontsize=TS, fontweight="bold")
    ax.tick_params(labelsize=TKS)
    ax.legend(fontsize=11, frameon=False, loc="lower right")

    # --- Panel B: Composition change ---
    ax = axes[1]
    comp_colors = []
    for _, row in sst.iterrows():
        if row["comp_padj"] < 0.05:
            if row["comp_beta"] < 0:
                comp_colors.append("#2980b9")  # blue = depleted
            else:
                comp_colors.append("#e67e22")  # orange = increased
        else:
            comp_colors.append("#bdc3c7")
    ax.barh(y_pos, sst["comp_beta"], color=comp_colors, edgecolor="black",
            linewidth=0.5, height=0.7)
    ax.axvline(0, color="black", linewidth=0.8, alpha=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([""] * n)  # no y-labels (shared with panel A)
    ax.set_xlabel("Composition change (case − control)", fontsize=LS)
    ax.set_title("B. Case-control proportion change", fontsize=TS, fontweight="bold")
    ax.tick_params(labelsize=TKS)

    # Add significance stars
    for i, (_, row) in enumerate(sst.iterrows()):
        if row["comp_padj"] < 0.001:
            stars = "***"
        elif row["comp_padj"] < 0.01:
            stars = "**"
        elif row["comp_padj"] < 0.05:
            stars = "*"
        else:
            stars = ""
        if stars:
            x_pos = row["comp_beta"]
            offset = -0.02 if x_pos < 0 else 0.01
            ax.text(x_pos + offset, i, stars, fontsize=13, fontweight="bold",
                    va="center", ha="left" if x_pos >= 0 else "right",
                    color=comp_colors[i])

    # --- Panel C: Correlation scatter ---
    ax = axes[2]
    ax.scatter(sst["abs_comp_beta"], sst["gwas_logp"],
               c="#D55E00", s=80, alpha=0.8,
               edgecolors="black", linewidth=0.5, zorder=3)

    for _, row in sst.iterrows():
        ax.annotate(row["supertype"].replace("Sst_", "").replace("Sst Chodl_", "SstChodl_"),
                     (row["abs_comp_beta"], row["gwas_logp"]),
                     fontsize=9, xytext=(4, 3), textcoords="offset points")

    r, p = stats.spearmanr(sst["abs_comp_beta"], sst["gwas_logp"])
    ax.text(0.03, 0.97, f"Spearman r = {r:.2f}\np = {p:.2e}",
            transform=ax.transAxes, fontsize=12, va="top",
            fontstyle="italic", color="#555",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    ax.set_xlabel("|Composition change|", fontsize=LS - 2)
    ax.set_ylabel(r"GWAS $-\log_{10}(p)$", fontsize=LS - 2)
    ax.set_title("C. Convergence", fontsize=TS, fontweight="bold")
    ax.tick_params(labelsize=TKS)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ══════════════════════════════════════════════════════════════════════════
# Figure 3: Focused SST Composition Correlation (2-panel)
# ══════════════════════════════════════════════════════════════════════════
def figure3_sst_composition_correlation(comp_df, output_path):
    """2-panel: SST |composition| vs GWAS, and significance vs significance."""
    cowplot_style()

    LS = 20; TS = 22; TKS = 15; ANNOT = 11

    sst = comp_df[comp_df["is_sst"]]

    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    # --- Panel A: |comp_beta| vs GWAS -log10(p) ---
    ax = axes[0]

    # Color by dual significance
    for _, row in sst.iterrows():
        gwas_sig = row["p_fdr"] < 0.05
        comp_sig = row["comp_padj"] < 0.05
        if gwas_sig and comp_sig:
            color = "#c0392b"  # dual-hit: dark red
            size = 120
        elif gwas_sig:
            color = "#D55E00"  # GWAS only
            size = 80
        else:
            color = "#bdc3c7"  # neither
            size = 60
        ax.scatter(row["abs_comp_beta"], row["gwas_logp"],
                   c=color, s=size, alpha=0.85,
                   edgecolors="black", linewidth=0.5, zorder=3)

    for _, row in sst.iterrows():
        ax.annotate(row["supertype"], (row["abs_comp_beta"], row["gwas_logp"]),
                    fontsize=ANNOT, fontweight="bold",
                    xytext=(5, 3), textcoords="offset points")

    r, p = stats.spearmanr(sst["abs_comp_beta"], sst["gwas_logp"])
    ax.text(0.03, 0.97,
            f"Spearman r = {r:.2f}\np = {p:.2e}\nn = {len(sst)}",
            transform=ax.transAxes, fontsize=14, va="top",
            fontstyle="italic", color="#333",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#f8f8f8", edgecolor="#ccc"))

    # Legend
    legend_els = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#c0392b",
               markersize=12, markeredgecolor="black", label="GWAS + composition sig"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#D55E00",
               markersize=10, markeredgecolor="black", label="GWAS sig only"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#bdc3c7",
               markersize=8, markeredgecolor="black", label="Not GWAS sig"),
    ]
    ax.legend(handles=legend_els, fontsize=12, frameon=False, loc="center right")

    ax.set_xlabel("|Composition change| (case − control)", fontsize=LS)
    ax.set_ylabel(r"GWAS enrichment $-\log_{10}(p)$", fontsize=LS)
    ax.set_title("A. SST: effect size convergence", fontsize=TS, fontweight="bold")
    ax.tick_params(labelsize=TKS)

    # --- Panel B: significance vs significance ---
    ax = axes[1]

    for _, row in sst.iterrows():
        gwas_sig = row["p_fdr"] < 0.05
        comp_sig = row["comp_padj"] < 0.05
        if gwas_sig and comp_sig:
            color = "#c0392b"
            size = 120
        elif gwas_sig:
            color = "#D55E00"
            size = 80
        else:
            color = "#bdc3c7"
            size = 60
        ax.scatter(row["comp_logp"], row["gwas_logp"],
                   c=color, s=size, alpha=0.85,
                   edgecolors="black", linewidth=0.5, zorder=3)

    for _, row in sst.iterrows():
        ax.annotate(row["supertype"], (row["comp_logp"], row["gwas_logp"]),
                    fontsize=ANNOT, fontweight="bold",
                    xytext=(5, 3), textcoords="offset points")

    r, p = stats.spearmanr(sst["comp_logp"], sst["gwas_logp"])
    ax.text(0.03, 0.97,
            f"Spearman r = {r:.2f}\np = {p:.2e}\nn = {len(sst)}",
            transform=ax.transAxes, fontsize=14, va="top",
            fontstyle="italic", color="#333",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#f8f8f8", edgecolor="#ccc"))

    ax.set_xlabel(r"Composition change $-\log_{10}(p)$", fontsize=LS)
    ax.set_ylabel(r"GWAS enrichment $-\log_{10}(p)$", fontsize=LS)
    ax.set_title("B. SST: significance convergence", fontsize=TS, fontweight="bold")
    ax.tick_params(labelsize=TKS)

    plt.tight_layout(w_pad=3)
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ══════════════════════════════════════════════════════════════════════════
# Figure 4: SST Deep Dive — RBH Combined Enrichment + Cortical Depth
# ══════════════════════════════════════════════════════════════════════════
def figure4_sst_rbh_depth(rbh_enrich, layer_df, output_path):
    """
    3-panel SST deep dive using RBH combined taxonomy enrichment results
    intersected with SEA-AD MERFISH cortical depth data.

    A) Horizontal bar chart of SST subtypes by -log10(p) from RBH combined,
       colored by cortical layer (Upper/Middle/Deep)
    B) Scatter: cortical depth vs GWAS enrichment -log10(p)
    C) Summary: enrichment statistics and layer breakdown

    Parameters
    ----------
    rbh_enrich : pd.DataFrame
        RBH combined enrichment results (503 types).
    layer_df : pd.DataFrame
        SST layer classification with supertype, avg_depth, layer_group.
    output_path : str or Path
    """
    cowplot_style()

    LS = 18; TS = 20; TKS = 14; ANNOT = 11

    # Extract SST types from RBH combined
    sst = rbh_enrich[rbh_enrich["supertype"].str.startswith("Sst")].copy()
    sst["-log10p"] = -np.log10(sst["p_value"].clip(lower=1e-300))

    # Merge with layer info
    sst = sst.merge(
        layer_df[["supertype", "avg_depth", "layer_group"]],
        on="supertype", how="left"
    )
    sst = sst.sort_values("-log10p", ascending=True)

    layer_colors = {
        "Upper": "#E74C3C",
        "Middle": "#F39C12",
        "Deep": "#3498DB",
        "Unknown": "#95A5A6",
    }

    n = len(sst)
    y_pos = np.arange(n)

    fig, axes = plt.subplots(1, 3, figsize=(24, max(8, n * 0.5)),
                              gridspec_kw={"width_ratios": [2.5, 2, 1.5], "wspace": 0.35})

    # --- Panel A: SST subtypes by -log10(p), colored by layer ---
    ax = axes[0]
    colors = [layer_colors.get(lg, "#95A5A6") for lg in sst["layer_group"]]
    ax.barh(y_pos, sst["-log10p"], color=colors, edgecolor="black",
            linewidth=0.5, height=0.7)

    # FDR and Bonferroni lines (computed within 503-type context)
    bonf_thresh = -np.log10(0.05 / 503)
    ax.axvline(bonf_thresh, color="#e74c3c", linestyle="--", linewidth=1.2,
               alpha=0.6, label="Bonferroni (503 types)")
    fdr_types = sst[sst["p_fdr"] < 0.05]
    if len(fdr_types) > 0:
        fdr_line = -np.log10(fdr_types["p_value"].max())
        ax.axvline(fdr_line, color="#f39c12", linestyle=":", linewidth=1.2,
                   alpha=0.6, label="FDR < 0.05")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(sst["supertype"], fontsize=13, fontweight="bold")
    ax.set_xlabel(r"GWAS enrichment $-\log_{10}(p)$", fontsize=LS)
    ax.set_title("A. SST subtype enrichment\n(RBH combined taxonomy, 503 types)",
                 fontsize=TS, fontweight="bold")
    ax.tick_params(labelsize=TKS)
    ax.legend(fontsize=11, frameon=False, loc="lower right")

    # --- Panel B: Depth vs enrichment scatter ---
    ax = axes[1]
    valid = sst.dropna(subset=["avg_depth"])
    for lg in ["Upper", "Middle", "Deep"]:
        mask = valid["layer_group"] == lg
        if mask.any():
            ax.scatter(
                valid.loc[mask, "avg_depth"], valid.loc[mask, "-log10p"],
                c=layer_colors[lg], s=120, label=lg,
                edgecolors="black", linewidth=0.7, zorder=3,
            )

    # Label all points
    for _, row in valid.iterrows():
        ha = "left"
        x_off = 0.012
        # Adjust specific overlapping labels
        if row["supertype"] in ("Sst_22", "Sst_23"):
            ha = "right"
            x_off = -0.012
        ax.annotate(
            row["supertype"], (row["avg_depth"] + x_off, row["-log10p"]),
            fontsize=ANNOT, fontweight="bold", ha=ha, va="center",
        )

    # Spearman correlation: depth vs enrichment
    r, p = stats.spearmanr(valid["avg_depth"], valid["-log10p"])
    ax.text(0.03, 0.97,
            f"Spearman r = {r:.2f}\np = {p:.2e}\nn = {len(valid)}",
            transform=ax.transAxes, fontsize=12, va="top",
            fontstyle="italic", color="#333",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#f8f8f8", edgecolor="#ccc"))

    ax.set_xlabel("Cortical depth (0 = pial, 1 = WM)", fontsize=LS)
    ax.set_ylabel(r"GWAS enrichment $-\log_{10}(p)$", fontsize=LS)
    ax.set_title("B. Cortical depth vs enrichment", fontsize=TS, fontweight="bold")
    ax.tick_params(labelsize=TKS)

    # Layer background shading
    ax.axvspan(0, 0.35, color="#E74C3C", alpha=0.04, zorder=0)
    ax.axvspan(0.35, 0.55, color="#F39C12", alpha=0.04, zorder=0)
    ax.axvspan(0.55, 1.0, color="#3498DB", alpha=0.04, zorder=0)
    ax.text(0.17, ax.get_ylim()[0] + 0.3, "Upper", ha="center", fontsize=10,
            color="#E74C3C", alpha=0.6, fontweight="bold")
    ax.text(0.45, ax.get_ylim()[0] + 0.3, "Mid", ha="center", fontsize=10,
            color="#F39C12", alpha=0.6, fontweight="bold")
    ax.text(0.75, ax.get_ylim()[0] + 0.3, "Deep", ha="center", fontsize=10,
            color="#3498DB", alpha=0.6, fontweight="bold")

    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=c,
               markersize=12, markeredgecolor="black", label=l)
        for l, c in layer_colors.items() if l != "Unknown"
    ]
    ax.legend(handles=handles, fontsize=11, frameon=False, loc="upper right")

    # --- Panel C: Summary statistics ---
    ax = axes[2]
    ax.axis("off")

    n_bonf = (sst["p_bonferroni"] < 0.05).sum()
    n_fdr = (sst["p_fdr"] < 0.05).sum()
    n_total = len(sst)

    upper = valid[valid["layer_group"] == "Upper"]
    middle = valid[valid["layer_group"] == "Middle"]
    deep = valid[valid["layer_group"] == "Deep"]

    summary_lines = [
        "SST Enrichment Summary",
        "(RBH Combined Taxonomy)\n",
        f"Total SST supertypes:  {n_total}",
        f"FDR-significant:       {n_fdr}",
        f"Bonferroni-significant: {n_bonf}\n",
        f"Upper-layer types:  {len(upper)}",
        f"  Mean -log10(p): {upper['-log10p'].mean():.1f}",
        f"  FDR-sig: {(upper['p_fdr'] < 0.05).sum()}/{len(upper)}\n",
        f"Middle-layer types: {len(middle)}",
        f"  Mean -log10(p): {middle['-log10p'].mean():.1f}",
        f"  FDR-sig: {(middle['p_fdr'] < 0.05).sum()}/{len(middle)}\n",
        f"Deep-layer types:   {len(deep)}",
        f"  Mean -log10(p): {deep['-log10p'].mean():.1f}",
        f"  FDR-sig: {(deep['p_fdr'] < 0.05).sum()}/{len(deep)}\n",
        f"Depth-enrichment r: {r:.2f}",
        f"  p = {p:.2e}",
    ]
    summary_text = "\n".join(summary_lines)
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
            fontsize=13, verticalalignment="top", fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.6", facecolor="#f8f8f8", edgecolor="#ccc"))
    ax.set_title("C. Summary", fontsize=TS, fontweight="bold")

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════
def main():
    out_dir = FIGURES_DIR / "manuscript"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Generating manuscript figures for rnaseq_specificity_enrichment_results.md\n")

    # Load data
    seaad_enrich = pd.read_csv(str(TABLES_DIR / "seaad_magma_scz_enrichment_all_supertypes.csv"))
    rbh_enrich = pd.read_csv(str(TABLES_DIR / "rbh_combined_enrichment.csv"))
    comp_df = pd.read_csv(str(TABLES_DIR / "gwas_vs_casecontrol_composition.csv"))
    layer_df = pd.read_csv(str(TABLES_DIR / "sst_supertype_layer_info.csv"))

    # Figure 1
    print("Figure 1: SEA-AD Manhattan scatter...")
    figure1_seaad_manhattan(seaad_enrich, out_dir / "fig1_seaad_manhattan.png")

    # Figure 2
    print("Figure 2: SST dual-hit...")
    figure2_sst_dual_hit(comp_df, out_dir / "fig2_sst_dual_hit.png")

    # Figure 3
    print("Figure 3: SST composition correlation...")
    figure3_sst_composition_correlation(comp_df, out_dir / "fig3_sst_composition_correlation.png")

    # Figure 4
    print("Figure 4: SST RBH combined enrichment + depth...")
    figure4_sst_rbh_depth(rbh_enrich, layer_df, out_dir / "fig4_sst_rbh_depth.png")

    print(f"\nDone. All figures in: {out_dir}")
    print("\nFull figure list for the markdown:")
    print("  1. fig1_seaad_manhattan.png")
    print("  3. gene_driver_scatter_Sst_2.png (existing)")
    print("  4. rbh_combined_taxonomy_top_enrichments.png (existing)")
    print("  5. fig4_sst_rbh_depth.png (NEW — RBH combined SST + depth)")
    print("  6. fig2_sst_dual_hit.png")
    print("  7. fig3_sst_composition_correlation.png")


if __name__ == "__main__":
    main()
