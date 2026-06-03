"""
Shared constants and utility functions for SCZ Xenium presentation analysis.

This module centralizes all configuration that was previously duplicated across
~20 analysis scripts. Import from here instead of hardcoding constants.

Shared project-wide constants (SAMPLE_TO_DX, EXCLUDE_SAMPLES, SUBCLASS_TO_CLASS,
CLASS_COLORS, CORTICAL_LAYERS) live in code/modules/constants.py and are
re-exported here for backward compatibility.

NOTE on LAYER_COLORS: These are the PRESENTATION layer colors, intentionally
different from code/modules/depth_model.py's LAYER_COLORS (which uses different
hues for L5 and L6). Do not "fix" this — the divergence is deliberate.
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Import shared constants from modules/constants.py and re-export
# so all existing `from config import SAMPLE_TO_DX` etc. continue to work.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "modules"))
from constants import (  # noqa: F401 — re-exported
    SAMPLE_TO_DX, CONTROL_SAMPLES, SCZ_SAMPLES, EXCLUDE_SAMPLES,
    CORTICAL_LAYERS, CLASS_COLORS, SUBCLASS_TO_CLASS,
)

# ──────────────────────────────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────────────────────────────

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
H5AD_DIR = os.path.join(BASE_DIR, "output", "h5ad")
# SEA-AD snRNAseq reference — validation only (137K cells, 36K genes)
# NOT required for the core pipeline; used for proportion validation plots
SNRNASEQ_REF_PATH = os.path.join(BASE_DIR, "data", "reference",
                                  "seaad_mtg_snrnaseq_reference.h5ad")
# MERFISH spatial reference — required for depth model training (step 04)
# Also used for validation plots (proportion comparison, depth comparison)
MERFISH_PATH = os.path.join(BASE_DIR, "data", "reference",
                             "SEAAD_MTG_MERFISH.2024-12-11.h5ad")
METADATA_PATH = os.path.join(BASE_DIR, "data", "sample_metadata.xlsx")

# Reference availability flags — checked at import time so scripts can
# skip validation plots gracefully when reference data isn't available
HAS_MERFISH_REF = os.path.exists(MERFISH_PATH)
HAS_SNRNASEQ_REF = os.path.exists(SNRNASEQ_REF_PATH)
CRUMBLR_DIR = os.path.join(BASE_DIR, "output", "crumblr")
PRESENTATION_DIR = os.path.join(BASE_DIR, "output", "presentation")
MARKER_ANALYSIS_DIR = os.path.join(BASE_DIR, "output", "marker_analysis")

# Xenium panel metadata CSVs (downloaded from 10x)
PANEL_5K_PATH = os.path.expanduser(
    "~/Downloads/XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv")
PANEL_V1_PATH = os.path.expanduser(
    "~/Downloads/Xenium_hBrain_v1_metadata.csv")

# MERSCOPE data directory
MERSCOPE_H5AD_DIR = os.path.join(BASE_DIR, "output", "merscope_h5ad")

# ──────────────────────────────────────────────────────────────────────
# MapMyCells confidence filter
# ──────────────────────────────────────────────────────────────────────
# HANN subclass confidence filter — DISABLED.
# Previously set to 0.500, which removed ~25% of cells (concentrated in
# hard-to-distinguish deep-layer excitatory and rare interneuron types).
# Now set to 0.0 to keep all QC-pass cells with their raw HANN labels.
# Mapping quality is exposed in the viewer instead of hard-filtering.
SUBCLASS_CONF_THRESH = 0.0

# ──────────────────────────────────────────────────────────────────────
# Presentation colors
# ──────────────────────────────────────────────────────────────────────

BG_COLOR = "#0a0a0a"
DX_COLORS = {"Control": "#4fc3f7", "SCZ": "#ef5350"}

# Presentation layer colors — distinct from depth_model.py's palette
LAYER_COLORS = {
    "L1":       (0.9, 0.3, 0.3),
    "L2/3":     (0.3, 0.8, 0.3),
    "L4":       (0.3, 0.3, 0.9),
    "L5":       (0.9, 0.9, 0.2),
    "L6":       (0.9, 0.5, 0.1),
    "WM":       (0.6, 0.6, 0.6),
    "Vascular": (0.95, 0.3, 0.6),
}
LAYER_ORDER = ["L1", "L2/3", "L4", "L5", "L6", "WM", "Vascular"]
# CORTICAL_LAYERS imported from modules/constants

# ──────────────────────────────────────────────────────────────────────
# Vulnerable cell type groups
# ──────────────────────────────────────────────────────────────────────

SST_TYPES = ["Sst_2", "Sst_22", "Sst_25", "Sst_20", "Sst_3"]
SST_COLORS = {
    "Sst_2":  "#e41a1c",
    "Sst_22": "#377eb8",
    "Sst_25": "#984ea3",
    "Sst_20": "#ff7f00",
    "Sst_3":  "#ffff33",
}

L6B_TYPES = ["L6b_1", "L6b_2", "L6b_4"]
L6B_COLORS = {
    "L6b_1": "#1f78b4",
    "L6b_2": "#33a02c",
    "L6b_4": "#e31a1c",
}

# ──────────────────────────────────────────────────────────────────────
# Representative samples — closest to group medians for vulnerable Sst
# proportion, with good cortical coverage (>27K cells).
# Control median=1.59%: Br1113=1.66% (d=0.06), Br2719=1.53% (d=0.06)
# SCZ median=1.33%: Br5373=1.32% (d=0.008), Br6032=1.34% (d=0.008)
# ──────────────────────────────────────────────────────────────────────

REPRESENTATIVE_SAMPLES = [
    ("Br1113", "Control"),
    ("Br2719", "Control"),
    ("Br5373", "SCZ"),
    ("Br6032", "SCZ"),
]

# ──────────────────────────────────────────────────────────────────────
# Cell class classification
# ──────────────────────────────────────────────────────────────────────

# CLASS_COLORS and SUBCLASS_TO_CLASS imported from modules/constants

# Prefix-based inference (for crumblr/snRNAseq scripts)
GABA_PREFIXES = ['Sst', 'Pvalb', 'Vip', 'Lamp5', 'Sncg', 'Pax6', 'Chandelier']
GLUT_PREFIXES = ['L2/3', 'L4', 'L5', 'L6']
NN_PREFIXES = ['Astro', 'Oligo', 'OPC', 'Micro', 'Endo', 'VLMC', 'SMC', 'Pericyte']

# ──────────────────────────────────────────────────────────────────────
# Spatial plot defaults
# ──────────────────────────────────────────────────────────────────────

ALL_CELL_COLOR = "#444444"
MARKER_SIZE_BG = 0.3

# ──────────────────────────────────────────────────────────────────────
# Data loading utilities
# ──────────────────────────────────────────────────────────────────────


def load_cells(sample_id, cortical_only=False, extra_obs_columns=None,
               qc_mode='corr'):
    """Load QC-pass cells from one Xenium sample.

    Central data loader for all analysis scripts. Uses correlation classifier
    labels (corr_subclass, corr_supertype) if available, falling back to HANN
    labels. Applies the appropriate QC filter based on qc_mode.

    Parameters
    ----------
    sample_id : str
        Sample ID (e.g., "Br6389").
    cortical_only : bool
        If True, restrict to cortical layers (L1-L6) only.
    extra_obs_columns : list of str, optional
        Additional obs columns to include (e.g., ["predicted_norm_depth"]).
    qc_mode : str
        QC filtering strategy:
        - 'corr'   : use corr_qc_pass (default; spatial QC + margin filter + doublet exclusion)
        - 'hybrid' : use hybrid_qc_pass (nuclear doublet-resolved). Falls back to 'corr' if
                     hybrid_qc_pass column doesn't exist yet.

    Returns
    -------
    DataFrame with columns: sample_id, subclass_label, supertype_label,
    spatial_domain, layer, qc_pass, x, y, plus any extra_obs_columns.
    """
    fpath = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    adata = ad.read_h5ad(fpath, backed="r")

    # Determine which label columns to use
    has_corr = "corr_subclass" in adata.obs.columns
    subclass_col = "corr_subclass" if has_corr else "subclass_label"
    supertype_col = "corr_supertype" if has_corr else "supertype_label"

    cols = ["sample_id", subclass_col, supertype_col,
            "spatial_domain", "layer", "qc_pass"]
    if has_corr:
        cols.append("corr_qc_pass")
    else:
        cols.append("subclass_label_confidence")

    # Include hybrid QC columns if available and requested
    has_hybrid = "hybrid_qc_pass" in adata.obs.columns
    if has_hybrid:
        for hcol in ["hybrid_qc_pass", "nuclear_doublet_status"]:
            if hcol in adata.obs.columns and hcol not in cols:
                cols.append(hcol)

    if extra_obs_columns:
        for c in extra_obs_columns:
            if c not in cols and c in adata.obs.columns:
                cols.append(c)

    obs = adata.obs[cols].copy()

    # Rename to standard column names for downstream compatibility
    if has_corr:
        obs = obs.rename(columns={
            subclass_col: "subclass_label",
            supertype_col: "supertype_label",
        })

    coords = adata.obsm["spatial"]
    obs["x"] = coords[:, 0]
    obs["y"] = coords[:, 1]
    obs["layer"] = obs["layer"].astype(str)
    obs["subclass_label"] = obs["subclass_label"].astype(str)
    obs["supertype_label"] = obs["supertype_label"].astype(str)

    # Determine effective QC mode
    effective_mode = qc_mode
    if effective_mode == 'hybrid' and not has_hybrid:
        effective_mode = 'corr'  # fallback if optional nuclear resolution hasn't run

    # QC pass filter
    # For hybrid mode, hybrid_qc_pass IS the complete QC column (includes
    # high-UMI rescue, so we don't intersect with qc_pass which excludes them)
    if effective_mode == 'hybrid':
        mask = obs["hybrid_qc_pass"] == True
    else:
        mask = obs["qc_pass"] == True

    if cortical_only:
        mask = mask & (obs["spatial_domain"] == "Cortical")
    obs = obs[mask]

    # Apply classifier-specific QC filter (only needed for non-hybrid modes)
    if effective_mode != 'hybrid':
        if has_corr:
            obs = obs[obs["corr_qc_pass"] == True]
        else:
            obs = obs[obs["subclass_label_confidence"].astype(float)
                      >= SUBCLASS_CONF_THRESH]

    return obs.copy()


# Backward-compatible wrappers
def load_cortical(sample_id):
    """Load QC-pass cortical cells (L1-L6) with spatial coordinates."""
    return load_cells(sample_id, cortical_only=True)


def load_all_cells(sample_id):
    """Load all QC-pass cells with spatial coordinates and layer info."""
    return load_cells(sample_id, cortical_only=False)


def load_merfish_cortical():
    """Load MERFISH cortical cells with manual depth annotation.

    NOTE: This uses the MERFISH spatial reference (not snRNAseq) because it
    has spatial coordinates and manual depth annotations. Used ONLY for
    depth-related analyses and validation plots.

    Returns DataFrame with columns: donor, subclass, supertype, depth,
    layer_annotation. Only includes cells with manual depth annotations
    in cortical layers (L1-L6).

    Raises FileNotFoundError if the MERFISH reference is not available.
    Check HAS_MERFISH_REF before calling if you want to skip gracefully.
    """
    if not HAS_MERFISH_REF:
        raise FileNotFoundError(
            f"MERFISH reference not found at {MERFISH_PATH}. "
            "See data/README.md for download instructions."
        )
    adata = ad.read_h5ad(MERFISH_PATH, backed="r")
    obs = adata.obs[["Donor ID", "Subclass", "Supertype",
                      "Normalized depth from pia",
                      "Layer annotation"]].copy()

    # Filter to cells with manual depth annotation and cortical layers
    obs = obs.dropna(subset=["Normalized depth from pia"])
    obs["Layer annotation"] = obs["Layer annotation"].astype(str)
    obs = obs[obs["Layer annotation"].isin(CORTICAL_LAYERS)]

    obs = obs.rename(columns={
        "Donor ID": "donor",
        "Subclass": "subclass",
        "Supertype": "supertype",
        "Normalized depth from pia": "depth",
        "Layer annotation": "layer_annotation",
    })
    obs["subclass"] = obs["subclass"].astype(str)
    obs["supertype"] = obs["supertype"].astype(str)
    obs["depth"] = obs["depth"].astype(float)

    return obs


def load_snrnaseq_reference(level="Subclass", neurons_only=False):
    """Load Nicole's Sea-AD snRNAseq reference for proportion comparisons.

    This is the primary cell type reference (137K cells, 36K genes, 5 donors).
    Provides ground-truth cell type proportions from dissociated tissue.
    NOT required for the core pipeline — used only for validation/benchmarking.

    Parameters
    ----------
    level : str
        Taxonomy level: "Subclass" or "Supertype".
    neurons_only : bool
        If True, restrict to neuronal classes only.

    Returns
    -------
    DataFrame with columns: donor, celltype, class_label.

    Raises FileNotFoundError if the snRNAseq reference is not available.
    Check HAS_SNRNASEQ_REF before calling if you want to skip gracefully.
    """
    if not HAS_SNRNASEQ_REF:
        raise FileNotFoundError(
            f"snRNAseq reference not found at {SNRNASEQ_REF_PATH}. "
            "See data/README.md for download instructions. "
            "This file is only needed for validation plots, not the core pipeline."
        )
    adata = ad.read_h5ad(SNRNASEQ_REF_PATH, backed="r")
    obs = adata.obs[["donor_id", "Class", level]].copy()
    obs = obs.rename(columns={
        "donor_id": "donor",
        level: "celltype",
        "Class": "class_label",
    })
    obs["donor"] = obs["donor"].astype(str)
    obs["celltype"] = obs["celltype"].astype(str)
    obs["class_label"] = obs["class_label"].astype(str)

    # Map Nicole's class labels to our standard short names
    class_map = {
        "Neuronal: Glutamatergic": "Glutamatergic",
        "Neuronal: GABAergic": "GABAergic",
        "Non-neuronal and Non-neural": "Non-neuronal",
    }
    obs["class_label"] = obs["class_label"].map(class_map).fillna("Unknown")

    if neurons_only:
        obs = obs[obs["class_label"].isin(["Glutamatergic", "GABAergic"])]

    return obs


def compute_reference_proportions(level="Subclass", neurons_only=False):
    """Compute per-donor mean proportions from the snRNAseq reference.

    Parameters
    ----------
    level : str
        Taxonomy level: "Subclass" or "Supertype".
    neurons_only : bool
        If True, restrict to neuronal classes only.

    Returns
    -------
    DataFrame with columns: celltype, ref_mean, ref_std.
    """
    obs = load_snrnaseq_reference(level=level, neurons_only=neurons_only)

    records = []
    for donor in obs["donor"].unique():
        donor_df = obs[obs["donor"] == donor]
        total = len(donor_df)
        counts = donor_df["celltype"].value_counts()
        for ct, n in counts.items():
            records.append({"donor": donor, "celltype": ct,
                            "proportion": n / total})

    df = pd.DataFrame(records)
    stats = df.groupby("celltype")["proportion"].agg(
        ["mean", "std"]).reset_index()
    stats.columns = ["celltype", "ref_mean", "ref_std"]
    return stats


def load_sample_adata(sample_id, cortical_only=True, qc_mode='corr'):
    """Load a full AnnData object with standardized labels and QC filtering.

    Unlike load_cells(), this returns the full AnnData (including X matrix),
    suitable for pseudobulk DE or any analysis needing gene counts.

    Parameters
    ----------
    sample_id : str
        Sample ID (e.g., "Br6389").
    cortical_only : bool
        If True, restrict to spatial_domain == 'Cortical'.
    qc_mode : str
        QC filtering strategy (same as load_cells):
        - 'corr'   : use corr_qc_pass (default; spatial QC + margin filter + doublet exclusion)
        - 'hybrid' : use hybrid_qc_pass (nuclear doublet-resolved). Falls back to 'corr'
                     if hybrid_qc_pass doesn't exist.

    Returns
    -------
    AnnData with .obs columns standardized to subclass_label, supertype_label.
    """
    fpath = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    adata = ad.read_h5ad(fpath)

    has_corr = "corr_subclass" in adata.obs.columns
    has_hybrid = "hybrid_qc_pass" in adata.obs.columns

    effective_mode = qc_mode
    if effective_mode == 'hybrid' and not has_hybrid:
        effective_mode = 'corr'  # fallback if optional nuclear resolution hasn't run

    # Build QC mask
    # For hybrid mode, hybrid_qc_pass IS the complete QC column (includes
    # high-UMI rescue, so we don't intersect with qc_pass which excludes them)
    if effective_mode == 'hybrid':
        mask = adata.obs["hybrid_qc_pass"] == True
    else:
        mask = adata.obs["qc_pass"] == True

    if cortical_only:
        mask = mask & (adata.obs["spatial_domain"] == "Cortical")

    # Classifier-specific QC (only needed for non-hybrid modes)
    if effective_mode != 'hybrid':
        if has_corr:
            mask = mask & (adata.obs["corr_qc_pass"] == True)
        else:
            mask = mask & (adata.obs["subclass_label_confidence"].astype(float)
                           >= SUBCLASS_CONF_THRESH)

    adata = adata[mask].copy()

    # Standardize label columns
    if has_corr:
        adata.obs["subclass_label"] = adata.obs["corr_subclass"].astype(str)
        adata.obs["supertype_label"] = adata.obs["corr_supertype"].astype(str)
    else:
        adata.obs["subclass_label"] = adata.obs["subclass_label"].astype(str)
        adata.obs["supertype_label"] = adata.obs["supertype_label"].astype(str)

    adata.obs["sample_id"] = adata.obs["sample_id"].astype(str)

    return adata


# ──────────────────────────────────────────────────────────────────────
# Spatial visualization helpers
# ──────────────────────────────────────────────────────────────────────


def style_dark_axis(ax):
    """Apply standard dark presentation styling to a spatial axis."""
    ax.set_facecolor(BG_COLOR)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_inset(ax_main, df, loc="lower right"):
    """Draw a small full-tissue inset with layer-colored cells.

    Parameters
    ----------
    ax_main : matplotlib.axes.Axes
        The main axis to attach the inset to.
    df : DataFrame
        All QC-pass cells (from load_all_cells), must have x, y, layer columns.
    loc : str
        Location for the inset (e.g., "lower right", "lower left").
    """
    ax_inset = inset_axes(ax_main, width="30%", height="30%", loc=loc,
                           borderpad=0.5)
    ax_inset.set_facecolor("#111111")

    # Non-cortical cells dim
    non_cortical = df[~df["layer"].isin(CORTICAL_LAYERS)]
    if len(non_cortical) > 0:
        ax_inset.scatter(non_cortical["x"], non_cortical["y"],
                         s=0.05, c="#333333", alpha=0.3,
                         rasterized=True, linewidths=0)

    # Cortical cells colored by layer
    for layer in ["L1", "L2/3", "L4", "L5", "L6"]:
        layer_cells = df[df["layer"] == layer]
        if len(layer_cells) > 0:
            ax_inset.scatter(layer_cells["x"], layer_cells["y"],
                             s=0.08, c=[LAYER_COLORS[layer]], alpha=0.5,
                             rasterized=True, linewidths=0)

    ax_inset.set_aspect("equal")
    ax_inset.invert_yaxis()
    ax_inset.set_xticks([])
    ax_inset.set_yticks([])
    for spine in ax_inset.spines.values():
        spine.set_color("#555555")
        spine.set_linewidth(0.8)

    ax_inset.text(0.5, 0.02, "full section", transform=ax_inset.transAxes,
                  ha="center", va="bottom", fontsize=7, color="#aaaaaa",
                  fontstyle="italic")


def draw_layer_shading(ax, df, alpha=0.12):
    """Draw very faint layer-colored scatter for cortical context.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to draw on.
    df : DataFrame
        Cells to shade (should include WM for context). Must have x, y, layer.
    alpha : float
        Transparency for the shading dots (default 0.12).
    """
    for layer in LAYER_ORDER:
        layer_cells = df[df["layer"] == layer]
        if len(layer_cells) < 10:
            continue
        color = LAYER_COLORS[layer]
        ax.scatter(layer_cells["x"], layer_cells["y"],
                   s=0.8, c=[color], alpha=alpha,
                   rasterized=True, linewidths=0, zorder=1)


# ──────────────────────────────────────────────────────────────────────
# Statistics helpers
# ──────────────────────────────────────────────────────────────────────


def format_pval(p):
    """Format p-value for display on figures."""
    if p < 0.001:
        return f"p = {p:.1e}"
    elif p < 0.01:
        return f"p = {p:.3f}"
    elif p < 0.1:
        return f"p = {p:.2f}"
    else:
        return f"p = {p:.2f}"


# ──────────────────────────────────────────────────────────────────────
# Inferred-class colors (short names: Glut/GABA/NN/Other)
# ──────────────────────────────────────────────────────────────────────
# Used by crumblr volcano/effect plots and snRNAseq concordance scatters.
# These map the short class labels returned by infer_class() to hex colors.
# Distinct from CLASS_COLORS (which uses long names: Glutamatergic/GABAergic/Non-neuronal).
INFER_CLASS_COLORS = {
    'Glut': '#00ADF8',
    'GABA': '#F05A28',
    'NN': '#808080',
    'Other': '#999999',
}

# ──────────────────────────────────────────────────────────────────────
# Shared plotting utilities
# ──────────────────────────────────────────────────────────────────────


def pct_formatter(val, pos):
    """Format a proportion (0-1 scale) as a human-readable percentage string.

    Intended for use with matplotlib.ticker.FuncFormatter on log-scale axes
    displaying proportions. Used by proportion scatter plots.
    """
    pct = val * 100
    if pct >= 1:
        return f"{pct:.0f}%"
    elif pct >= 0.1:
        return f"{pct:.1f}%"
    elif pct >= 0.01:
        return f"{pct:.2f}%"
    else:
        return f"{pct:.3f}%"


def style_dark_boxplot(ax, ctrl_vals, scz_vals, ylabel=None, title=None,
                       pval=None, pval_label="", subtitle=None,
                       point_size=40, jitter_range=0.12, show_means=False):
    """Draw a dark-background boxplot + strip plot comparing Control vs SCZ.

    This is the shared pattern used by plot_aggregated_boxplots.py and
    plot_xenium_composition_boxplots.py. Draws paired boxplots with
    translucent colored boxes, jittered individual data points, and
    optional p-value annotation.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to draw on.
    ctrl_vals : array-like
        Values for the Control group.
    scz_vals : array-like
        Values for the SCZ group.
    ylabel : str, optional
        Y-axis label.
    title : str, optional
        Panel title (displayed above the axis).
    pval : float, optional
        P-value to annotate. If None, no annotation is drawn.
    pval_label : str
        Prefix for the p-value annotation (e.g., "Welch's t-test",
        "crumblr"). Empty string for no prefix.
    subtitle : str, optional
        Subtitle text placed just above the panel (e.g., cell type list).
    point_size : float
        Size of individual data points (default 40).
    jitter_range : float
        Half-width of horizontal jitter for strip points (default 0.12).
    show_means : bool
        If True, display group means alongside the p-value.

    Returns
    -------
    float or None
        The p-value passed in (for convenience in chaining).
    """
    ctrl_vals = np.asarray(ctrl_vals, dtype=float)
    scz_vals = np.asarray(scz_vals, dtype=float)

    positions = [0, 1]
    bp = ax.boxplot([ctrl_vals, scz_vals], positions=positions,
                    widths=0.5, patch_artist=True, showfliers=False,
                    medianprops=dict(color="white", linewidth=1.5),
                    whiskerprops=dict(color="#888888", linewidth=1.2),
                    capprops=dict(color="#888888", linewidth=1.2))

    for patch, color in zip(bp['boxes'],
                            [DX_COLORS["Control"], DX_COLORS["SCZ"]]):
        patch.set_facecolor(color)
        patch.set_alpha(0.4)
        patch.set_edgecolor("white")
        patch.set_linewidth(0.8)

    # Individual points with jitter
    n_ctrl = len(ctrl_vals)
    n_scz = len(scz_vals)
    for i, (vals, dx) in enumerate([(ctrl_vals, "Control"),
                                     (scz_vals, "SCZ")]):
        jitter = np.random.default_rng(42).uniform(
            -jitter_range, jitter_range, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals,
                   c=DX_COLORS[dx], s=point_size, alpha=0.85,
                   edgecolors="white", linewidths=0.5, zorder=5)

    ax.set_xticks(positions)
    ax.set_xticklabels([f"Ctrl\n(n={n_ctrl})", f"SCZ\n(n={n_scz})"],
                       fontsize=12, color="white")

    if title:
        ax.set_title(title, fontsize=15, fontweight="bold",
                     color="white", pad=6)

    # P-value annotation
    if pval is not None:
        pval_str = format_pval(pval)
        pval_color = ("white" if pval < 0.05
                      else ("#cccccc" if pval < 0.1 else "#888888"))
        pval_weight = "bold" if pval < 0.1 else "normal"

        label_parts = []
        if pval_label:
            label_parts.append(f"{pval_label} {pval_str}")
        else:
            label_parts.append(pval_str)
        if show_means:
            label_parts.append(
                f"Ctrl: {np.mean(ctrl_vals):.2f}, "
                f"SCZ: {np.mean(scz_vals):.2f}")

        ax.text(0.5, 0.02 if not show_means else 0.96,
                "\n".join(label_parts),
                transform=ax.transAxes,
                ha="center",
                va="bottom" if not show_means else "top",
                fontsize=10, color=pval_color, fontweight=pval_weight,
                fontstyle="italic",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="#1a1a1a",
                          edgecolor="#444444", alpha=0.85)
                if show_means else {})

    if subtitle:
        ax.text(0.5, 1.01, subtitle, transform=ax.transAxes,
                ha="center", va="bottom",
                fontsize=11, color="#aaaaaa", fontstyle="italic")

    if ylabel:
        ax.set_ylabel(ylabel, fontsize=15, color="white")

    ax.set_facecolor(BG_COLOR)
    ax.tick_params(colors="white", labelsize=11)
    for spine in ax.spines.values():
        spine.set_color("#555555")
    ax.grid(axis="y", alpha=0.15, color="#555555")

    return pval


def style_proportion_scatter(ax, x, y, celltypes, title="", xlabel="",
                             ylabel="", max_labels=15, show_corr=True,
                             classify_fn=None):
    """Draw a log-log proportion scatter with percentage tick formatting.

    This is the shared pattern used by plot_merfish_vs_xenium_proportions.py,
    plot_cropped_proportions.py, and plot_predicted_proportion_scatter.py.
    Plots cell type proportions on log-log axes with a diagonal identity line,
    class-colored dots, correlation statistics, and deviant-point labels.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to draw on.
    x, y : array-like
        Proportion values for the two platforms being compared.
    celltypes : array-like of str
        Cell type names (same length as x and y).
    title, xlabel, ylabel : str
        Axis labels and title.
    max_labels : int
        Maximum number of deviant points to label (default 15).
    show_corr : bool
        If True (default), show Pearson r (log-scale) and Spearman rho.
    classify_fn : callable, optional
        Function that takes a celltype name and returns (hex_color, class_name).
        Defaults to classify_celltype.
    """
    from scipy.stats import pearsonr, spearmanr

    if classify_fn is None:
        classify_fn = classify_celltype

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    celltypes = np.asarray(celltypes)

    colors = [classify_fn(ct)[0] for ct in celltypes]

    ax.scatter(x, y, c=colors, s=70, alpha=0.8, edgecolors="white",
               linewidths=0.5, zorder=5)

    # Label most deviant points (by log-ratio distance from diagonal)
    valid = (x > 0) & (y > 0)
    if valid.sum() > 0 and max_labels > 0:
        log_ratio = np.full_like(x, 0.0)
        log_ratio[valid] = np.abs(np.log2(
            (y[valid] + 1e-7) / (x[valid] + 1e-7)))
        top_idx = np.argsort(log_ratio)[-max_labels:]
        for idx in top_idx:
            if log_ratio[idx] > 0:
                ax.annotate(celltypes[idx],
                            (x[idx], y[idx]),
                            fontsize=9, color="#dddddd", alpha=0.9,
                            xytext=(5, 5), textcoords="offset points")

    # Diagonal identity line
    pos_x = x[x > 0]
    pos_y = y[y > 0]
    if len(pos_x) > 0 and len(pos_y) > 0:
        lo = min(pos_x.min(), pos_y.min()) * 0.3
        hi = max(x.max(), y.max()) * 3
        ax.plot([lo, hi], [lo, hi], "--", color="#888888", linewidth=1.5,
                alpha=0.6, zorder=1)

    # Correlation statistics
    if show_corr and valid.sum() > 3:
        r, p = pearsonr(np.log10(x[valid]), np.log10(y[valid]))
        rho, _ = spearmanr(x[valid], y[valid])
        ax.text(0.04, 0.96,
                f"r = {r:.2f} (log-scale)\n\u03c1 = {rho:.2f}"
                f"\nn = {valid.sum()} types",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=14, color="#dddddd",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#333333",
                          edgecolor="#555555", alpha=0.85))

    ax.set_xscale("log")
    ax.set_yscale("log")

    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_formatter(mticker.FuncFormatter(pct_formatter))
        axis.set_minor_formatter(mticker.NullFormatter())

    ax.set_xlabel(xlabel, fontsize=16, color="white")
    ax.set_ylabel(ylabel, fontsize=16, color="white")
    ax.set_title(title, fontsize=20, fontweight="bold", color="white", pad=10)
    ax.set_facecolor(BG_COLOR)
    ax.tick_params(colors="white", labelsize=13)
    for spine in ax.spines.values():
        spine.set_color("#555555")
    ax.grid(True, alpha=0.2, color="#555555")


def plot_volcano(ax, df, celltype_col="celltype", logfc_col="logFC",
                 pval_col="P.Value", fdr_col="FDR",
                 class_colors=None, classify_fn=None,
                 fdr_thresh_bold=0.05, fdr_thresh_light=0.10,
                 title="", xlabel="logFC (SCZ vs Control)",
                 ylabel="-log10(p-value)", dot_size=60):
    """Draw a volcano plot: logFC vs -log10(p-value), colored by cell class.

    This is the shared pattern used by plot_crumblr_results.py and
    plot_seaad_crumblr_results.py. Points are colored by inferred cell class,
    with FDR-significant types labeled.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to draw on.
    df : DataFrame
        Must contain columns for celltype, logFC, p-value, and FDR.
    celltype_col, logfc_col, pval_col, fdr_col : str
        Column names in df.
    class_colors : dict, optional
        Mapping of class name to hex color. Defaults to INFER_CLASS_COLORS.
    classify_fn : callable, optional
        Function mapping celltype -> class name. Defaults to infer_class.
    fdr_thresh_bold : float
        FDR threshold for bold labels (default 0.05).
    fdr_thresh_light : float
        FDR threshold for lighter labels (default 0.10).
    title, xlabel, ylabel : str
        Axis labels and title.
    dot_size : float
        Scatter dot size (default 60).
    """
    if class_colors is None:
        class_colors = INFER_CLASS_COLORS
    if classify_fn is None:
        classify_fn = infer_class

    df = df.copy()
    df['_class'] = df[celltype_col].apply(classify_fn)
    df['_nlog10p'] = -np.log10(df[pval_col])

    for cls in ['Glut', 'GABA', 'NN', 'Other']:
        mask = df['_class'] == cls
        if mask.sum() == 0:
            continue
        sub = df[mask]
        ax.scatter(sub[logfc_col], sub['_nlog10p'],
                   c=class_colors.get(cls, '#999999'), s=dot_size,
                   alpha=0.7, label=cls,
                   edgecolors='white', linewidth=0.5, zorder=5)

    # Label FDR-significant types
    for _, row in df.iterrows():
        if row[fdr_col] < fdr_thresh_bold:
            ax.annotate(row[celltype_col],
                        (row[logfc_col], row['_nlog10p']),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=9, fontweight='bold', alpha=0.9)
        elif row[fdr_col] < fdr_thresh_light:
            ax.annotate(row[celltype_col],
                        (row[logfc_col], row['_nlog10p']),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=8, alpha=0.7)

    # Reference lines
    if df[pval_col].min() < 0.05:
        ax.axhline(-np.log10(0.05), color='grey', linestyle='--',
                   alpha=0.5, linewidth=1)
        ax.text(ax.get_xlim()[1] * 0.95, -np.log10(0.05) + 0.05, 'p=0.05',
                ha='right', fontsize=9, color='grey', fontstyle='italic')
    ax.axvline(0, color='grey', alpha=0.3, linewidth=0.5)

    # Summary annotation
    n_fdr = (df[fdr_col] < fdr_thresh_bold).sum()
    ax.text(0.98, 0.98,
            f'n={len(df)} types\n{n_fdr} FDR<{fdr_thresh_bold}',
            transform=ax.transAxes, ha='right', va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(fontsize=11, loc='upper left')


# ──────────────────────────────────────────────────────────────────────
# Cell type classification
# ──────────────────────────────────────────────────────────────────────


def infer_class(celltype):
    """Infer broad cell class from celltype name using prefix matching.

    Returns 'GABA', 'Glut', 'NN', or 'Other'.
    """
    for p in GABA_PREFIXES:
        if celltype.startswith(p):
            return 'GABA'
    for p in GLUT_PREFIXES:
        if celltype.startswith(p):
            return 'Glut'
    for p in NN_PREFIXES:
        if celltype.startswith(p):
            return 'NN'
    return 'Other'


def classify_celltype(ct):
    """Return (hex_color, class_name) for a cell type name.

    Uses case-insensitive substring matching for robustness with
    both subclass and supertype names.
    """
    ct_lower = ct.lower()
    glut_markers = ["l2/3", "l4 ", "l4_", "l5 ", "l5_", "l5/6", "l6 ", "l6_",
                    "l6b", "it_", " it"]
    gaba_markers = ["sst", "pvalb", "vip", "lamp5", "sncg", "chandelier", "pax6"]
    if any(g in ct_lower for g in glut_markers):
        return CLASS_COLORS["Glutamatergic"], "Glutamatergic"
    elif any(g in ct_lower for g in gaba_markers):
        return CLASS_COLORS["GABAergic"], "GABAergic"
    else:
        return CLASS_COLORS["Non-neuronal"], "Non-neuronal"
