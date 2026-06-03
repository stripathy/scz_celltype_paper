"""
plotting — Publication-quality figures for SCZ cell-type enrichment analysis.

All figures follow consistent style:
  - Title fontsize >= 20, label fontsize >= 16, tick fontsize >= 13
  - DPI = 150 for print quality
  - sns.despine() on all axes
  - Large, legible text optimized for presentation/publication

Submodules:
  enrichment.py   — Core enrichment figures (bar plot, shared/unique, conditional)
  sst_spatial.py  — SST-focused figures (panels, layers, upper-layer, volcano)
  gene_ephys.py   — Gene-electrophysiology figures (scatters, HCN1-SAG)
"""
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for script usage
import matplotlib.pyplot as plt
import seaborn as sns

from ..config import (
    TITLE_FONTSIZE, LABEL_FONTSIZE, TICK_FONTSIZE,
    LEGEND_FONTSIZE, FIGURE_DPI,
)


def set_figure_style():
    """Set global matplotlib/seaborn style for all figures."""
    sns.set_style("white")
    plt.rcParams.update({
        "font.size": TICK_FONTSIZE,
        "axes.titlesize": TITLE_FONTSIZE,
        "axes.labelsize": LABEL_FONTSIZE,
        "xtick.labelsize": TICK_FONTSIZE,
        "ytick.labelsize": TICK_FONTSIZE,
        "legend.fontsize": LEGEND_FONTSIZE,
        "figure.dpi": FIGURE_DPI,
    })


# Re-export all plot functions for backward compatibility.
# Existing scripts can continue to use:
#   from scz_celltype_enrichment.plotting import plot_enrichment_all_supertypes
from .enrichment import (
    plot_enrichment_all_supertypes,
    plot_shared_vs_unique,
    plot_conditional_results,
)
from .sst_spatial import (
    plot_sst_enrichment_panel,
    plot_sst_pvalues_by_layer,
    plot_upper_layer_enrichment,
    plot_sst_layer_ephys,
    plot_depth_volcano,
)
from .gene_ephys import (
    plot_gene_ephys_scatters,
    plot_hcn1_sag_relationship,
)

__all__ = [
    "set_figure_style",
    "plot_enrichment_all_supertypes",
    "plot_shared_vs_unique",
    "plot_conditional_results",
    "plot_sst_enrichment_panel",
    "plot_sst_pvalues_by_layer",
    "plot_upper_layer_enrichment",
    "plot_sst_layer_ephys",
    "plot_depth_volcano",
    "plot_gene_ephys_scatters",
    "plot_hcn1_sag_relationship",
]
