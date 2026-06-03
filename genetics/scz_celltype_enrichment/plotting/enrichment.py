"""
enrichment.py — Core GWAS enrichment figures.

Figure 1: All supertypes enrichment bar plot
Figure 5: Shared vs unique gene signal (Jaccard + shared genes)
Figure 6: Conditional analysis results (marginal vs conditional)
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np

from ..config import (
    TITLE_FONTSIZE, LABEL_FONTSIZE, LEGEND_FONTSIZE, FIGURE_DPI, FDR_THRESHOLD,
)
from ..utils import classify_supertypes_by_class
from . import set_figure_style


# ============================================================================
# Figure 1: All supertypes enrichment bar plot
# ============================================================================

def plot_enrichment_all_supertypes(results_df, output_path):
    """
    Bar plot of -log10(p) for all ~137 supertypes, colored by cell class.

    Horizontal lines at Bonferroni and FDR significance thresholds.

    Parameters
    ----------
    results_df : pd.DataFrame
        Enrichment results with columns: supertype, p_value, p_fdr, p_bonferroni
    output_path : str or Path
        Path to save figure.
    """
    set_figure_style()

    df = results_df.sort_values("p_value").copy()
    df["-log10p"] = -np.log10(df["p_value"].clip(lower=1e-20))

    # Classify by cell class
    class_map = classify_supertypes_by_class(df["supertype"].tolist())
    df["class"] = df["supertype"].map(class_map)

    class_colors = {
        "GABAergic": "#E74C3C",
        "Glutamatergic": "#3498DB",
        "Non-neuronal": "#95A5A6",
    }
    df["color"] = df["class"].map(class_colors)

    fig, ax = plt.subplots(figsize=(24, 6))

    bars = ax.bar(range(len(df)), df["-log10p"], color=df["color"], width=0.8)

    # Significance thresholds
    n_types = len(df)
    bonf_line = -np.log10(0.05 / n_types)
    fdr_types = df[df["p_fdr"] < FDR_THRESHOLD]
    if len(fdr_types) > 0:
        fdr_line = -np.log10(fdr_types["p_value"].max())
    else:
        fdr_line = bonf_line

    ax.axhline(bonf_line, color="red", linestyle="--", linewidth=1.5,
               label=f"Bonferroni (p<0.05/{n_types})")
    ax.axhline(fdr_line, color="orange", linestyle="--", linewidth=1.5,
               label="FDR < 0.05")

    # Labels for top hits
    for i, row in df.head(20).iterrows():
        idx = df.index.get_loc(i)
        ax.text(idx, row["-log10p"] + 0.2, row["supertype"],
                rotation=90, ha="center", va="bottom", fontsize=8)

    ax.set_xlabel("Cell types (sorted by p-value)", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("-log10(p-value)", fontsize=LABEL_FONTSIZE)
    ax.set_title("SCZ GWAS Enrichment Across All Cell Types (MAGMA Gene Property Analysis)",
                 fontsize=TITLE_FONTSIZE)
    ax.set_xticks([])

    # Legend
    handles = [
        mpatches.Patch(color=c, label=l)
        for l, c in class_colors.items()
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=LEGEND_FONTSIZE)

    sns.despine()
    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ============================================================================
# Figure 5: Shared vs unique gene signal
# ============================================================================

def plot_shared_vs_unique(jaccard_df, shared_genes_df, output_path):
    """
    Multi-panel: Jaccard similarity heatmap + shared gene summary.

    Parameters
    ----------
    jaccard_df : pd.DataFrame
        Square Jaccard similarity matrix.
    shared_genes_df : pd.DataFrame
        Shared genes table with gene, n_types columns.
    output_path : str or Path
    """
    set_figure_style()

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Panel A: Jaccard heatmap
    ax = axes[0]
    sns.heatmap(jaccard_df, ax=ax, cmap="YlOrRd", vmin=0, vmax=0.5,
                annot=True, fmt=".2f", annot_kws={"size": 7},
                xticklabels=True, yticklabels=True,
                cbar_kws={"shrink": 0.6, "label": "Jaccard similarity"})
    ax.set_title("A) Gene Set Overlap (Top 200)", fontsize=TITLE_FONTSIZE - 2)
    ax.tick_params(axis="x", rotation=90, labelsize=8)
    ax.tick_params(axis="y", rotation=0, labelsize=8)

    # Panel B: Shared genes bar
    ax = axes[1]
    if len(shared_genes_df) > 0:
        top_shared = shared_genes_df.head(25)
        ax.barh(range(len(top_shared)), top_shared["n_types"],
                color="#E74C3C", alpha=0.7)
        ax.set_yticks(range(len(top_shared)))
        ax.set_yticklabels(top_shared["gene"], fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel("# Cell Types Shared", fontsize=LABEL_FONTSIZE)
        ax.set_title("B) Most Shared Genes", fontsize=TITLE_FONTSIZE - 2)
    else:
        ax.text(0.5, 0.5, "No shared genes found",
                ha="center", va="center", transform=ax.transAxes)
    sns.despine(ax=ax)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ============================================================================
# Figure 6: Conditional analysis results
# ============================================================================

def plot_conditional_results(conditional_df, output_path):
    """
    Scatter: marginal vs conditional -log10(p) for independent cell types.

    Parameters
    ----------
    conditional_df : pd.DataFrame
        From conditional.forward_selection, with marginal_p and conditional_p.
    output_path : str or Path
    """
    set_figure_style()

    df = conditional_df.copy()
    df["marginal_logp"] = -np.log10(df["marginal_p"].clip(lower=1e-20))
    df["conditional_logp"] = -np.log10(df["conditional_p"].clip(lower=1e-20))

    # Color by subclass
    subclasses = df["subclass"].unique()
    palette = sns.color_palette("husl", len(subclasses))
    color_map = dict(zip(subclasses, palette))

    fig, ax = plt.subplots(figsize=(10, 8))

    for sc in subclasses:
        mask = df["subclass"] == sc
        sub = df[mask]
        ax.scatter(sub["marginal_logp"], sub["conditional_logp"],
                  c=[color_map[sc]], s=100, label=sc,
                  edgecolors="k", linewidth=0.5, zorder=3)

    # Label all points
    for _, row in df.iterrows():
        ax.annotate(row["supertype"],
                   (row["marginal_logp"], row["conditional_logp"]),
                   fontsize=8, ha="left", va="bottom",
                   xytext=(3, 3), textcoords="offset points")

    # Reference line (y=x)
    max_val = max(df["marginal_logp"].max(), df["conditional_logp"].max()) * 1.1
    ax.plot([0, max_val], [0, max_val], "k--", alpha=0.3, linewidth=1)

    # Significance line
    ax.axhline(-np.log10(0.05), color="red", linestyle=":", alpha=0.5,
               label="p = 0.05")

    ax.set_xlabel("Marginal -log10(p)", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Conditional -log10(p)", fontsize=LABEL_FONTSIZE)
    ax.set_title("Independent Cell Types: Marginal vs Conditional Significance",
                 fontsize=TITLE_FONTSIZE)
    ax.legend(loc="upper left", fontsize=LEGEND_FONTSIZE - 1, ncol=2)

    sns.despine(ax=ax)
    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")
