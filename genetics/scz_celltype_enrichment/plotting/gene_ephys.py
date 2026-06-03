"""
gene_ephys.py — Gene-electrophysiology correlation figures.

Figure 7: Gene-ephys scatter plots (top correlated genes vs ephys feature)
Figure 8: HCN1 vs SAG relationship across supertypes
"""
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats

from ..config import (
    TITLE_FONTSIZE, LABEL_FONTSIZE, LEGEND_FONTSIZE, FIGURE_DPI,
)
from . import set_figure_style


# ============================================================================
# Figure 7: Gene-ephys scatter plots
# ============================================================================

def plot_gene_ephys_scatters(corr_df, feature_name, output_path, n_genes=9):
    """
    Grid of scatter plots: top correlated genes vs an ephys feature.

    Parameters
    ----------
    corr_df : pd.DataFrame
        Correlation results with gene, spearman_rho, pval_fdr.
    feature_name : str
        Name of ephys feature (e.g., 'SAG', 'TAU').
    output_path : str or Path
    n_genes : int
        Number of top genes to plot.
    """
    set_figure_style()

    top = corr_df.head(n_genes)
    ncols = 3
    nrows = (n_genes + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = axes.flatten() if hasattr(axes, "flatten") else [axes]

    for i, (_, row) in enumerate(top.iterrows()):
        if i >= len(axes):
            break
        ax = axes[i]
        rho = row["spearman_rho"]
        fdr = row.get("pval_fdr", row.get("pval", np.nan))
        direction = "+" if rho > 0 else "-"
        ax.set_title(
            f"{row['gene']}\n(rho={rho:.3f}, FDR={fdr:.1e})",
            fontsize=12
        )
        ax.text(0.5, 0.5, f"rho = {direction}{abs(rho):.3f}",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=16, fontweight="bold",
                color="red" if rho > 0 else "blue")
        ax.set_xlabel(f"{feature_name}", fontsize=12)
        ax.set_ylabel("Expression", fontsize=12)
        sns.despine(ax=ax)

    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(f"Top Genes Correlated with {feature_name}",
                 fontsize=TITLE_FONTSIZE, y=1.02)
    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ============================================================================
# Figure 8: HCN1 vs SAG relationship
# ============================================================================

def plot_hcn1_sag_relationship(ephys_by_supertype, specificity, output_path):
    """
    Scatter: HCN1 expression (from SEA-AD) vs mean SAG per supertype.

    Parameters
    ----------
    ephys_by_supertype : pd.DataFrame
        Mean ephys by supertype. Must have 'sag' column, index = supertype.
    specificity : pd.DataFrame
        Specificity matrix. Must have 'HCN1' in index.
    output_path : str or Path
    """
    set_figure_style()

    if "HCN1" not in specificity.index:
        print("  WARNING: HCN1 not found in specificity matrix, skipping plot")
        return

    # Get shared supertypes
    shared = sorted(set(ephys_by_supertype.index) & set(specificity.columns))
    if len(shared) < 3:
        print(f"  WARNING: Only {len(shared)} shared supertypes, skipping HCN1 plot")
        return

    hcn1 = specificity.loc["HCN1", shared].values
    sag = ephys_by_supertype.loc[shared, "mean_sag"].values

    rho, pval = stats.spearmanr(hcn1, sag)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(hcn1, sag, s=80, c="#E74C3C", edgecolors="k", linewidth=0.5, zorder=3)

    for i, st in enumerate(shared):
        ax.annotate(st, (hcn1[i], sag[i]), fontsize=8, ha="left", va="bottom",
                   xytext=(3, 3), textcoords="offset points")

    ax.set_xlabel("HCN1 Specificity (SEA-AD)", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Mean SAG", fontsize=LABEL_FONTSIZE)
    ax.set_title(f"HCN1 Expression vs SAG Across Supertypes\n"
                 f"(Spearman rho={rho:.3f}, p={pval:.1e})",
                 fontsize=TITLE_FONTSIZE)

    sns.despine(ax=ax)
    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")
