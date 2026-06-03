#!/usr/bin/env python3
"""
Cell & Nucleus Boundary Area by Cell Type.

Computes polygon areas for cell and nucleus boundaries across all 24 Xenium
samples, joins with cell type labels, and produces violin plots at the
subclass and supertype levels. Uses linear mixed effects models (random
intercept per sample) to test for Control vs SCZ differences in cell and
nucleus size, with FDR correction.

Output:
  output/presentation/boundary_areas_subclass.png
  output/presentation/boundary_areas_supertype_glut.png
  output/presentation/boundary_areas_supertype_gaba.png
  output/presentation/boundary_areas_supertype_nn.png
  output/presentation/boundary_areas_lme_forest.png
  output/presentation/cell_boundary_areas_by_subclass.csv
  output/presentation/cell_boundary_areas_by_supertype.csv
  output/presentation/cell_boundary_areas_lme_results.csv

Usage:
    python3 -u code/analysis/plot_cell_boundary_areas.py
"""

import os
import sys
import glob
import time
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (BG_COLOR, H5AD_DIR, PRESENTATION_DIR,
                    SUBCLASS_TO_CLASS, CLASS_COLORS, SAMPLE_TO_DX,
                    EXCLUDE_SAMPLES)

# ── SEA-AD subclass colors (from 07_export_viewer.py) ──
SUBCLASS_COLORS = {
    "L2/3 IT": "#B1EC30", "L4 IT": "#00E5E5", "L5 IT": "#50B2AD",
    "L5 ET": "#0D5B78", "L5/6 NP": "#3E9E64", "L6 IT": "#A19922",
    "L6 IT Car3": "#5100FF", "L6 CT": "#2D8CB8", "L6b": "#7044AA",
    "Sst": "#FF9900", "Sst Chodl": "#B1B10C", "Pvalb": "#D93137",
    "Vip": "#A45FBF", "Lamp5": "#DA808C", "Lamp5 Lhx6": "#935F50",
    "Sncg": "#DF70FF", "Pax6": "#71238C", "Chandelier": "#F641A8",
    "Astrocyte": "#665C47", "Oligodendrocyte": "#53776C",
    "OPC": "#374A45", "Microglia-PVM": "#94AF97",
    "Endothelial": "#8D6C62", "VLMC": "#697255",
}

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
RAW_DIR = os.path.join(BASE_DIR, "data/raw")
BG = BG_COLOR


# ══════════════════════════════════════════════════════════════════════
# DATA LOADING
# ══════════════════════════════════════════════════════════════════════

def find_boundary_file(sample_id, kind='cell'):
    """Find cell or nucleus boundary CSV for a sample."""
    pattern = os.path.join(RAW_DIR, f"*{sample_id}-{kind}_boundaries.csv.gz")
    matches = glob.glob(pattern)
    return matches[0] if matches else None


def _shoelace_area(verts):
    """Polygon area via the shoelace formula (no shapely dependency)."""
    x, y = verts[:, 0], verts[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def compute_polygon_areas(boundary_path):
    """Load boundary CSV → compute polygon area per cell_id.

    Returns dict {cell_id: area_um2}.
    """
    df = pd.read_csv(boundary_path, compression='gzip')
    grouped = df.groupby('cell_id', sort=False)

    areas = {}
    n_invalid = 0
    for cell_id, grp in grouped:
        verts = grp[['vertex_x', 'vertex_y']].values
        if len(verts) < 3:
            n_invalid += 1
            continue
        area = _shoelace_area(verts)
        if area > 0:
            areas[cell_id] = area
        else:
            n_invalid += 1

    return areas


def load_sample_areas(sample_id):
    """Load cell + nucleus areas and cell type labels for one sample.

    Returns DataFrame with columns:
        cell_id, sample_id, subclass_label, supertype_label,
        cell_area, nucleus_area
    """
    import anndata as ad

    cell_path = find_boundary_file(sample_id, 'cell')
    nuc_path = find_boundary_file(sample_id, 'nucleus')
    h5ad_path = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")

    if not cell_path or not nuc_path:
        print(f"  WARNING: Missing boundary files for {sample_id}")
        return pd.DataFrame()

    # Compute polygon areas
    t0 = time.time()
    cell_areas = compute_polygon_areas(cell_path)
    nuc_areas = compute_polygon_areas(nuc_path)
    t_areas = time.time() - t0

    # Load h5ad for cell type labels and QC
    adata = ad.read_h5ad(h5ad_path, backed='r')

    # Determine label columns
    has_corr = "corr_subclass" in adata.obs.columns
    sub_col = "corr_subclass" if has_corr else "subclass_label"
    sup_col = "corr_supertype" if has_corr else "supertype_label"

    # Determine QC column
    if "corr_qc_pass" in adata.obs.columns:
        qc_mask = adata.obs["corr_qc_pass"].values == True
    else:
        qc_mask = adata.obs["qc_pass"].values == True

    obs = adata.obs[[sub_col, sup_col]].copy()
    obs.columns = ['subclass_label', 'supertype_label']
    obs['qc_pass'] = qc_mask
    obs['cell_id'] = obs.index

    # Join with areas
    obs['cell_area'] = obs['cell_id'].map(cell_areas)
    obs['nucleus_area'] = obs['cell_id'].map(nuc_areas)
    obs['sample_id'] = sample_id
    obs['diagnosis'] = SAMPLE_TO_DX.get(sample_id, 'Unknown')

    # Filter: QC pass + have areas
    obs = obs[obs['qc_pass'] & obs['cell_area'].notna()
              & obs['nucleus_area'].notna()].copy()

    obs['subclass_label'] = obs['subclass_label'].astype(str)
    obs['supertype_label'] = obs['supertype_label'].astype(str)

    print(f"  {sample_id}: {len(obs):,} cells with areas "
          f"({t_areas:.1f}s for polygons)", flush=True)

    return obs[['cell_id', 'sample_id', 'diagnosis', 'subclass_label',
                'supertype_label', 'cell_area', 'nucleus_area']]


# ══════════════════════════════════════════════════════════════════════
# PLOTTING
# ══════════════════════════════════════════════════════════════════════

def style_violin_axis(ax, ylabel, title):
    """Standard dark styling for violin axes."""
    ax.set_facecolor('#111111')
    ax.set_ylabel(ylabel, fontsize=16, color='white', fontweight='bold')
    ax.set_title(title, fontsize=20, fontweight='bold', color='white', pad=12)
    ax.tick_params(axis='x', colors='white', labelsize=12)
    ax.tick_params(axis='y', colors='white', labelsize=13)
    for spine in ax.spines.values():
        spine.set_color('#555555')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


DX_COLORS = {"Control": "#4fc3f7", "SCZ": "#ef5350"}
DX_ORDER = ["Control", "SCZ"]


def plot_split_violins(ax, data_ctrl, data_scz, positions, colors,
                       area_col_label=''):
    """Draw split violins: Control (left half) vs SCZ (right half)."""
    from matplotlib.patches import PathPatch
    import matplotlib.path as mpath

    for i, pos in enumerate(positions):
        for dx_idx, (dx_data, dx_color, side) in enumerate([
            (data_ctrl[i], DX_COLORS['Control'], 'left'),
            (data_scz[i], DX_COLORS['SCZ'], 'right'),
        ]):
            if len(dx_data) < 5:
                continue
            parts = ax.violinplot([dx_data], positions=[pos],
                                   vert=False, showmedians=False,
                                   showextrema=False, widths=0.8)
            for body in parts['bodies']:
                # Clip to half
                vertices = body.get_paths()[0].vertices
                if side == 'left':
                    # Keep only bottom half (y <= pos)
                    vertices[:, 1] = np.clip(vertices[:, 1], -np.inf, pos)
                else:
                    # Keep only top half (y >= pos)
                    vertices[:, 1] = np.clip(vertices[:, 1], pos, np.inf)
                body.set_facecolor(dx_color)
                body.set_alpha(0.6)
                body.set_edgecolor('white')
                body.set_linewidth(0.3)

            # Median marker
            med = np.median(dx_data)
            y_offset = -0.15 if side == 'left' else 0.15
            ax.plot(med, pos + y_offset, 'o', color=dx_color,
                    markersize=4, zorder=5, markeredgecolor='white',
                    markeredgewidth=0.5)


def plot_subclass_violins(df):
    """Split violin plots of cell and nucleus areas by subclass, Control vs SCZ."""
    print("\n  Plotting subclass-level violins...", flush=True)

    # Order subclasses by median cell area (descending)
    medians = df.groupby('subclass_label')['cell_area'].median()
    order = medians.sort_values(ascending=True).index.tolist()
    counts = df['subclass_label'].value_counts()
    order = [s for s in order if counts.get(s, 0) >= 20]

    fig, axes = plt.subplots(1, 2, figsize=(22, max(10, len(order) * 0.55)),
                              facecolor=BG, sharey=True)

    for col_idx, (area_col, title) in enumerate([
        ('cell_area', 'Cell Boundary Area'),
        ('nucleus_area', 'Nucleus Boundary Area'),
    ]):
        ax = axes[col_idx]

        data_ctrl = []
        data_scz = []
        positions = []
        for i, subclass in enumerate(order):
            mask = df['subclass_label'] == subclass
            ctrl_vals = df.loc[mask & (df['diagnosis'] == 'Control'),
                               area_col].values
            scz_vals = df.loc[mask & (df['diagnosis'] == 'SCZ'),
                              area_col].values
            data_ctrl.append(ctrl_vals)
            data_scz.append(scz_vals)
            positions.append(i)

        if not positions:
            continue

        plot_split_violins(ax, data_ctrl, data_scz, positions,
                          [SUBCLASS_COLORS.get(s, '#888888') for s in order],
                          area_col)

        # Annotate n=
        for i, subclass in enumerate(order):
            n_ctrl = len(data_ctrl[i])
            n_scz = len(data_scz[i])
            ax.text(0.98, i, f' {n_ctrl:,} / {n_scz:,}',
                    va='center', ha='left', fontsize=8, color='#aaaaaa',
                    transform=ax.get_yaxis_transform())

        ax.set_yticks(range(len(order)))
        ax.set_yticklabels(order)
        ax.set_xlabel('Area (µm²)', fontsize=14, color='white')
        style_violin_axis(ax, '', title)

        all_vals = np.concatenate(data_ctrl + data_scz)
        if len(all_vals) > 0:
            ax.set_xlim(0, np.percentile(all_vals, 99))

    # Diagnosis legend
    from matplotlib.lines import Line2D
    legend_items = [
        Line2D([0], [0], marker='s', color='w', markerfacecolor=DX_COLORS['Control'],
               markersize=12, label='Control', linestyle='None'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor=DX_COLORS['SCZ'],
               markersize=12, label='SCZ', linestyle='None'),
    ]
    fig.legend(handles=legend_items, loc='lower center', ncol=2,
               fontsize=14, facecolor='#1a1a2e', edgecolor='#555555',
               labelcolor='white', bbox_to_anchor=(0.5, 0.005))

    n_ctrl = (df['diagnosis'] == 'Control').sum()
    n_scz = (df['diagnosis'] == 'SCZ').sum()
    fig.suptitle('Cell & Nucleus Boundary Areas by Subclass\n'
                 f'Control ({n_ctrl:,}) vs SCZ ({n_scz:,}) — '
                 f'n=Control/SCZ shown at right',
                 fontsize=20, fontweight='bold', color='white', y=0.98)

    plt.tight_layout(rect=[0, 0.03, 1, 0.94])
    outpath = os.path.join(PRESENTATION_DIR, 'boundary_areas_subclass.png')
    plt.savefig(outpath, dpi=150, facecolor=BG, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


def get_supertype_color(supertype, subclass, subtype_idx, n_subtypes):
    """Brightness-varied color for supertype based on parent subclass."""
    import matplotlib.colors as mcolors
    base = SUBCLASS_COLORS.get(subclass, '#888888')
    r, g, b = mcolors.hex2color(base)
    # Vary brightness from 0.7 to 1.3
    if n_subtypes > 1:
        factor = 0.7 + 0.6 * subtype_idx / (n_subtypes - 1)
    else:
        factor = 1.0
    r = min(1.0, r * factor)
    g = min(1.0, g * factor)
    b = min(1.0, b * factor)
    return (r, g, b)


def plot_supertype_violins(df, cell_class, class_label):
    """Split violin plots of cell and nucleus areas by supertype, Control vs SCZ."""
    class_df = df[df['subclass_label'].map(
        lambda s: SUBCLASS_TO_CLASS.get(s, '') == cell_class)].copy()

    if len(class_df) == 0:
        return

    # Only supertypes with n >= 50
    counts = class_df['supertype_label'].value_counts()
    valid_types = counts[counts >= 50].index
    class_df = class_df[class_df['supertype_label'].isin(valid_types)]

    if len(class_df) == 0:
        return

    # Order by median cell area
    medians = class_df.groupby('supertype_label')['cell_area'].median()
    order = medians.sort_values(ascending=True).index.tolist()

    # Build color map
    sup_to_sub = class_df.groupby('supertype_label')['subclass_label'].first()
    subclass_groups = {}
    for sup, sub in sup_to_sub.items():
        subclass_groups.setdefault(sub, []).append(sup)
    for sub in subclass_groups:
        subclass_groups[sub].sort()

    sup_colors = {}
    for sub, sups in subclass_groups.items():
        for idx, sup in enumerate(sups):
            sup_colors[sup] = get_supertype_color(sup, sub, idx, len(sups))

    n_types = len(order)
    fig_height = max(8, n_types * 0.42)
    fig, axes = plt.subplots(1, 2, figsize=(22, fig_height),
                              facecolor=BG, sharey=True)

    for col_idx, (area_col, title) in enumerate([
        ('cell_area', f'Cell Boundary Area — {class_label}'),
        ('nucleus_area', f'Nucleus Boundary Area — {class_label}'),
    ]):
        ax = axes[col_idx]

        data_ctrl = []
        data_scz = []
        positions = []
        for i, supertype in enumerate(order):
            mask = class_df['supertype_label'] == supertype
            ctrl_vals = class_df.loc[mask & (class_df['diagnosis'] == 'Control'),
                                     area_col].values
            scz_vals = class_df.loc[mask & (class_df['diagnosis'] == 'SCZ'),
                                    area_col].values
            data_ctrl.append(ctrl_vals)
            data_scz.append(scz_vals)
            positions.append(i)

        if not positions:
            continue

        plot_split_violins(ax, data_ctrl, data_scz, positions,
                          [sup_colors.get(s, '#888888') for s in order],
                          area_col)

        for i, supertype in enumerate(order):
            n_ctrl = len(data_ctrl[i])
            n_scz = len(data_scz[i])
            ax.text(0.98, i, f' {n_ctrl:,}/{n_scz:,}',
                    va='center', ha='left', fontsize=7, color='#aaaaaa',
                    transform=ax.get_yaxis_transform())

        ax.set_yticks(range(n_types))
        ax.set_yticklabels(order, fontsize=10)
        ax.set_xlabel('Area (µm²)', fontsize=14, color='white')
        style_violin_axis(ax, '', title)

        all_vals = np.concatenate(data_ctrl + data_scz)
        if len(all_vals) > 0:
            ax.set_xlim(0, np.percentile(all_vals, 99))

    # Diagnosis legend
    from matplotlib.lines import Line2D
    legend_items = [
        Line2D([0], [0], marker='s', color='w',
               markerfacecolor=DX_COLORS['Control'],
               markersize=10, label='Control', linestyle='None'),
        Line2D([0], [0], marker='s', color='w',
               markerfacecolor=DX_COLORS['SCZ'],
               markersize=10, label='SCZ', linestyle='None'),
    ]
    fig.legend(handles=legend_items, loc='lower center', ncol=2,
               fontsize=12, facecolor='#1a1a2e', edgecolor='#555555',
               labelcolor='white', bbox_to_anchor=(0.5, 0.003))

    n_ctrl = (class_df['diagnosis'] == 'Control').sum()
    n_scz = (class_df['diagnosis'] == 'SCZ').sum()
    fig.suptitle(f'Boundary Areas by Supertype — {class_label}\n'
                 f'Control ({n_ctrl:,}) vs SCZ ({n_scz:,}), '
                 f'supertypes with n ≥ 50',
                 fontsize=18, fontweight='bold', color='white', y=0.98)

    plt.tight_layout(rect=[0, 0.03, 1, 0.94])
    tag = cell_class.lower()[:4]
    outpath = os.path.join(PRESENTATION_DIR,
                            f'boundary_areas_supertype_{tag}.png')
    plt.savefig(outpath, dpi=150, facecolor=BG, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


def save_summary_csvs(df):
    """Save per-subclass and per-supertype area summary CSVs, stratified by dx."""
    for level, col in [('subclass', 'subclass_label'),
                        ('supertype', 'supertype_label')]:
        records = []
        for (name, dx), grp in df.groupby([col, 'diagnosis']):
            records.append({
                level: name,
                'diagnosis': dx,
                'n_cells': len(grp),
                'cell_area_median': round(grp['cell_area'].median(), 1),
                'cell_area_mean': round(grp['cell_area'].mean(), 1),
                'cell_area_std': round(grp['cell_area'].std(), 1),
                'nucleus_area_median': round(grp['nucleus_area'].median(), 1),
                'nucleus_area_mean': round(grp['nucleus_area'].mean(), 1),
                'nucleus_area_std': round(grp['nucleus_area'].std(), 1),
                'nc_ratio_median': round(
                    (grp['nucleus_area'] / grp['cell_area']).median(), 3),
            })
        out_df = pd.DataFrame(records).sort_values(
            ['cell_area_median', 'diagnosis'], ascending=[False, True])
        csv_path = os.path.join(PRESENTATION_DIR,
                                 f'cell_boundary_areas_by_{level}.csv')
        out_df.to_csv(csv_path, index=False)
        print(f"  Saved: {csv_path} ({len(out_df)} rows)")


# ══════════════════════════════════════════════════════════════════════
# MIXED EFFECTS MODELS
# ══════════════════════════════════════════════════════════════════════

def fit_lme_per_celltype(df, level_col='subclass_label',
                         min_cells_per_group=100,
                         min_median_cells_per_sample=None):
    """Fit linear mixed effects models for cell/nucleus area ~ diagnosis.

    For each cell type, fits:
        log(area) ~ C(diagnosis, Treatment(reference='Control')) + (1 | sample_id)

    Log-transform handles right-skewed area distributions. The SCZ coefficient
    represents the log-fold-change in area (SCZ vs Control).

    Parameters
    ----------
    df : DataFrame
        Pooled cell data (already excluding EXCLUDE_SAMPLES).
    level_col : str
        Column name for cell type grouping ('subclass_label' or
        'supertype_label').
    min_cells_per_group : int
        Minimum cells per diagnosis group to fit model.
    min_median_cells_per_sample : int or None
        If set, only test cell types where the median number of cells per
        sample is at least this value.

    Returns
    -------
    results_df : DataFrame
        One row per celltype × area_type with columns:
        celltype, area_type, beta_scz, se, pval, ci_lo, ci_hi,
        pct_change, n_ctrl, n_scz, n_samples_ctrl, n_samples_scz,
        converged.
    """
    level_short = level_col.replace('_label', '')
    print(f"\n  Fitting LME at {level_short} level "
          f"(log(area) ~ diagnosis + (1|sample))...", flush=True)
    t0 = time.time()

    # Pre-filter cell types by median cells per sample if requested
    celltypes = sorted(df[level_col].unique())
    if min_median_cells_per_sample is not None:
        per_sample = df.groupby([level_col, 'sample_id']).size().reset_index(
            name='n')
        medians = per_sample.groupby(level_col)['n'].median()
        passing = medians[medians >= min_median_cells_per_sample].index
        n_before = len(celltypes)
        celltypes = [ct for ct in celltypes if ct in passing]
        print(f"    Median-per-sample filter (≥{min_median_cells_per_sample}): "
              f"{len(celltypes)}/{n_before} types pass", flush=True)

    results = []

    for ct in celltypes:
        if ct in ('Unassigned', 'nan', ''):
            continue

        sub_df = df[df[level_col] == ct].copy()

        n_ctrl = (sub_df['diagnosis'] == 'Control').sum()
        n_scz = (sub_df['diagnosis'] == 'SCZ').sum()
        n_samples_ctrl = sub_df.loc[sub_df['diagnosis'] == 'Control',
                                     'sample_id'].nunique()
        n_samples_scz = sub_df.loc[sub_df['diagnosis'] == 'SCZ',
                                    'sample_id'].nunique()

        if n_ctrl < min_cells_per_group or n_scz < min_cells_per_group:
            continue

        if n_samples_ctrl < 3 or n_samples_scz < 3:
            continue

        for area_col, area_label in [('cell_area', 'cell'),
                                      ('nucleus_area', 'nucleus')]:
            sub_df['log_area'] = np.log(sub_df[area_col].clip(lower=1.0))
            sub_df['dx'] = pd.Categorical(
                sub_df['diagnosis'], categories=['Control', 'SCZ'])

            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    model = smf.mixedlm(
                        "log_area ~ dx",
                        data=sub_df,
                        groups=sub_df["sample_id"],
                    )
                    fit = model.fit(reml=True, method='lbfgs')

                    # Extract SCZ coefficient (dx[T.SCZ])
                    beta = fit.params['dx[T.SCZ]']
                    se = fit.bse['dx[T.SCZ]']
                    pval = fit.pvalues['dx[T.SCZ]']
                    ci_lo = beta - 1.96 * se
                    ci_hi = beta + 1.96 * se

                    # Convert log-scale beta to approximate percent change
                    pct_change = (np.exp(beta) - 1) * 100

                    converged = True

            except Exception as e:
                print(f"    WARNING: {ct}/{area_label} failed: {e}",
                      flush=True)
                beta = se = pval = ci_lo = ci_hi = pct_change = np.nan
                converged = False

            results.append({
                'celltype': ct,
                'area_type': area_label,
                'beta_scz': round(beta, 5) if not np.isnan(beta) else np.nan,
                'se': round(se, 5) if not np.isnan(se) else np.nan,
                'pval': pval,
                'ci_lo': round(ci_lo, 5) if not np.isnan(ci_lo) else np.nan,
                'ci_hi': round(ci_hi, 5) if not np.isnan(ci_hi) else np.nan,
                'pct_change_scz': round(pct_change, 2) if not np.isnan(pct_change) else np.nan,
                'n_ctrl': n_ctrl,
                'n_scz': n_scz,
                'n_samples_ctrl': n_samples_ctrl,
                'n_samples_scz': n_samples_scz,
                'converged': converged,
            })

    results_df = pd.DataFrame(results)

    if results_df.empty:
        print("    No cell types tested!", flush=True)
        return results_df

    # FDR correction (Benjamini-Hochberg) — separately for cell and nucleus
    for area_type in ['cell', 'nucleus']:
        mask = (results_df['area_type'] == area_type) & results_df['pval'].notna()
        if mask.sum() > 0:
            pvals = results_df.loc[mask, 'pval'].values
            _, fdr_vals, _, _ = multipletests(pvals, method='fdr_bh')
            results_df.loc[mask, 'fdr'] = fdr_vals
        else:
            results_df.loc[mask, 'fdr'] = np.nan

    elapsed = time.time() - t0
    n_tested = results_df['celltype'].nunique()
    n_sig_cell = ((results_df['area_type'] == 'cell') &
                  (results_df['fdr'] < 0.05)).sum()
    n_sig_nuc = ((results_df['area_type'] == 'nucleus') &
                 (results_df['fdr'] < 0.05)).sum()
    print(f"    Done in {elapsed:.1f}s — {n_tested} {level_short}s tested")
    print(f"    Significant (FDR<0.05): {n_sig_cell} cell area, "
          f"{n_sig_nuc} nucleus area", flush=True)

    return results_df


def compute_sample_level_means(df):
    """Compute sample-level mean areas per subclass for sanity checking.

    Returns DataFrame with one row per (subclass, sample_id) with mean areas.
    """
    grouped = df.groupby(['subclass_label', 'sample_id', 'diagnosis']).agg(
        cell_area_mean=('cell_area', 'mean'),
        nucleus_area_mean=('nucleus_area', 'mean'),
        n_cells=('cell_area', 'size'),
    ).reset_index()
    return grouped


def _infer_broad_class(celltype):
    """Infer broad class from celltype name (works for both subclass and supertype)."""
    # Direct lookup for subclass names
    if celltype in SUBCLASS_TO_CLASS:
        return SUBCLASS_TO_CLASS[celltype]
    # Prefix matching for supertype names (e.g., "L2/3 IT_1" → Glutamatergic)
    for subclass, broad in SUBCLASS_TO_CLASS.items():
        if celltype.startswith(subclass):
            return broad
    return 'Non-neuronal'


def plot_lme_forest(results_df, level='subclass', suffix=''):
    """Forest plot of LME effect sizes (SCZ vs Control) for cell and nucleus area.

    Two panels side-by-side: cell area (left), nucleus area (right).
    Cell types ordered by cell area effect size. Colored by broad class.
    Stars for FDR significance.

    Parameters
    ----------
    results_df : DataFrame
        From fit_lme_per_celltype. Must have 'celltype' column.
    level : str
        'subclass' or 'supertype' — affects title and filename.
    suffix : str
        Optional filename suffix (e.g., '_glut' for class-specific plots).
    """
    print(f"\n  Plotting LME forest plot ({level}{suffix})...", flush=True)

    n_types = results_df['celltype'].nunique()
    row_height = 0.45 if level == 'subclass' else 0.38
    fig_height = max(6, n_types * row_height)
    ylabel_fontsize = 13 if level == 'subclass' else 10
    marker_size = 8 if level == 'subclass' else 6
    pct_fontsize = 10 if level == 'subclass' else 8

    fig, axes = plt.subplots(1, 2, figsize=(20, fig_height),
                              facecolor=BG, sharey=True)

    for col_idx, (area_type, title) in enumerate([
        ('cell', 'Cell Boundary Area'),
        ('nucleus', 'Nucleus Boundary Area'),
    ]):
        ax = axes[col_idx]
        ax.set_facecolor('#111111')

        sub_df = results_df[results_df['area_type'] == area_type].copy()
        if sub_df.empty:
            continue

        # Order by effect size for cell area panel; use same order for nucleus
        if col_idx == 0:
            sub_df = sub_df.sort_values('beta_scz', ascending=True)
            order = sub_df['celltype'].tolist()
        else:
            sub_df = sub_df.set_index('celltype').loc[order].reset_index()

        positions = list(range(len(sub_df)))

        for i, (_, row) in enumerate(sub_df.iterrows()):
            celltype = row['celltype']
            beta = row['beta_scz']
            ci_lo = row['ci_lo']
            ci_hi = row['ci_hi']
            fdr = row['fdr']

            if np.isnan(beta):
                continue

            # Color by broad class
            broad_class = _infer_broad_class(celltype)
            color = CLASS_COLORS.get(broad_class, '#888888')

            # Horizontal error bar
            ax.plot([ci_lo, ci_hi], [i, i], color=color, linewidth=2,
                    alpha=0.8, solid_capstyle='round')
            # Point estimate
            ax.plot(beta, i, 'o', color=color, markersize=marker_size,
                    markeredgecolor='white', markeredgewidth=0.5, zorder=5)

            # FDR significance annotation
            if fdr < 0.001:
                star = '***'
            elif fdr < 0.01:
                star = '**'
            elif fdr < 0.05:
                star = '*'
            else:
                star = ''

            if star:
                ax.text(ci_hi + 0.003, i, star,
                        va='center', ha='left', fontsize=14, color='#ffd700',
                        fontweight='bold')

            # Percent change annotation on right
            pct = row['pct_change_scz']
            if not np.isnan(pct):
                pct_str = f"{pct:+.1f}%"
                ax.text(0.97, i, pct_str, va='center', ha='right',
                        fontsize=pct_fontsize, color='#cccccc',
                        transform=ax.get_yaxis_transform())

        # Zero line
        ax.axvline(0, color='#666666', linestyle='--', linewidth=1, alpha=0.7)

        ax.set_yticks(positions)
        ax.set_yticklabels(sub_df['celltype'].values)
        ax.set_xlabel('β (log-scale, SCZ vs Control)', fontsize=14, color='white')
        style_violin_axis(ax, '', f'{title}\nSCZ vs Control (LME)')
        ax.tick_params(axis='y', labelsize=ylabel_fontsize, colors='white')
        ax.tick_params(axis='x', labelsize=12, colors='white')

    # Class legend
    from matplotlib.lines import Line2D
    legend_items = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=CLASS_COLORS['Glutamatergic'],
               markersize=10, label='Glutamatergic', linestyle='None'),
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=CLASS_COLORS['GABAergic'],
               markersize=10, label='GABAergic', linestyle='None'),
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=CLASS_COLORS['Non-neuronal'],
               markersize=10, label='Non-neuronal', linestyle='None'),
    ]
    fig.legend(handles=legend_items, loc='lower center', ncol=3,
               fontsize=13, facecolor='#1a1a2e', edgecolor='#555555',
               labelcolor='white', bbox_to_anchor=(0.5, 0.003))

    level_label = level.capitalize()
    fig.suptitle(f'Mixed Effects Model: Cell Size by {level_label} (SCZ vs Control)\n'
                 'β = log-fold-change | random intercept: sample_id | '
                 '* FDR<0.05  ** FDR<0.01  *** FDR<0.001',
                 fontsize=16, fontweight='bold', color='white', y=0.98)

    plt.tight_layout(rect=[0, 0.05, 1, 0.93])
    outpath = os.path.join(PRESENTATION_DIR,
                            f'boundary_areas_lme_forest_{level}{suffix}.png')
    plt.savefig(outpath, dpi=150, facecolor=BG, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


def save_lme_results(results_df, level='subclass'):
    """Save LME results to CSV."""
    out_df = results_df.copy()
    # Round for readability
    for col in ['pval', 'fdr']:
        if col in out_df.columns:
            out_df[col] = out_df[col].apply(
                lambda x: f"{x:.2e}" if pd.notna(x) and x < 0.001 else
                (f"{x:.4f}" if pd.notna(x) else ''))

    csv_path = os.path.join(PRESENTATION_DIR,
                             f'cell_boundary_areas_lme_{level}.csv')
    out_df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path} ({len(out_df)} rows)")


def print_lme_summary(results_df, level='subclass'):
    """Print key LME results to console, sorted by p-value."""
    if results_df.empty:
        print(f"\n  No {level}-level LME results to print.")
        return

    print(f"\n  ── LME Results Summary ({level}) ──")
    for area_type in ['cell', 'nucleus']:
        print(f"\n  {area_type.upper()} AREA:")
        sub = results_df[results_df['area_type'] == area_type].sort_values('pval')
        for _, row in sub.iterrows():
            sig = ''
            if pd.notna(row['fdr']):
                if row['fdr'] < 0.001:
                    sig = ' ***'
                elif row['fdr'] < 0.01:
                    sig = ' **'
                elif row['fdr'] < 0.05:
                    sig = ' *'
            pval_str = f"{row['pval']:.2e}" if pd.notna(row['pval']) else "N/A"
            fdr_str = f"{row['fdr']:.3f}" if pd.notna(row['fdr']) else "N/A"
            pct_str = (f"{row['pct_change_scz']:+.1f}%"
                       if pd.notna(row['pct_change_scz']) else "N/A")
            name_width = 25 if level == 'supertype' else 20
            print(f"    {row['celltype']:{name_width}s}  β={row['beta_scz']:+.4f}  "
                  f"pct={pct_str:>7s}  p={pval_str}  FDR={fdr_str}{sig}")


def plot_sample_level_means(sample_means_df):
    """Dot plot of sample-level mean areas per subclass, colored by diagnosis.

    Each dot is one sample's mean for that subclass. This validates
    that the LME results are driven by consistent sample-level differences,
    not just cell-count imbalances.
    """
    print("\n  Plotting sample-level mean areas...", flush=True)

    # Exclude rare subclasses
    counts = sample_means_df.groupby('subclass_label')['sample_id'].nunique()
    valid = counts[counts >= 10].index
    plot_df = sample_means_df[sample_means_df['subclass_label'].isin(valid)].copy()

    if plot_df.empty:
        return

    # Order by overall median cell area
    medians = plot_df.groupby('subclass_label')['cell_area_mean'].median()
    order = medians.sort_values(ascending=True).index.tolist()

    fig, axes = plt.subplots(1, 2, figsize=(20, max(9, len(order) * 0.5)),
                              facecolor=BG, sharey=True)

    for col_idx, (area_col, title) in enumerate([
        ('cell_area_mean', 'Mean Cell Area per Sample'),
        ('nucleus_area_mean', 'Mean Nucleus Area per Sample'),
    ]):
        ax = axes[col_idx]
        ax.set_facecolor('#111111')

        for i, subclass in enumerate(order):
            mask = plot_df['subclass_label'] == subclass
            for dx, color in DX_COLORS.items():
                dx_mask = mask & (plot_df['diagnosis'] == dx)
                vals = plot_df.loc[dx_mask, area_col].values
                jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(vals))
                offset = -0.15 if dx == 'Control' else 0.15
                ax.scatter(vals, np.full(len(vals), i) + offset + jitter * 0.3,
                          s=30, c=color, alpha=0.7, edgecolors='white',
                          linewidths=0.3, zorder=5)

                # Group mean bar
                if len(vals) > 0:
                    mean_val = vals.mean()
                    ax.plot([mean_val, mean_val],
                           [i + offset - 0.12, i + offset + 0.12],
                           color=color, linewidth=2.5, alpha=0.9, zorder=6)

        ax.set_yticks(range(len(order)))
        ax.set_yticklabels(order)
        ax.set_xlabel('Area (µm²)', fontsize=14, color='white')
        style_violin_axis(ax, '', title)

    # Legend
    from matplotlib.lines import Line2D
    legend_items = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=DX_COLORS['Control'],
               markersize=10, label='Control', linestyle='None'),
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=DX_COLORS['SCZ'],
               markersize=10, label='SCZ', linestyle='None'),
    ]
    fig.legend(handles=legend_items, loc='lower center', ncol=2,
               fontsize=13, facecolor='#1a1a2e', edgecolor='#555555',
               labelcolor='white', bbox_to_anchor=(0.5, 0.003))

    fig.suptitle('Sample-Level Mean Boundary Areas by Subclass\n'
                 'Each dot = one sample, bars = group mean',
                 fontsize=16, fontweight='bold', color='white', y=0.98)

    plt.tight_layout(rect=[0, 0.04, 1, 0.93])
    outpath = os.path.join(PRESENTATION_DIR, 'boundary_areas_sample_means.png')
    plt.savefig(outpath, dpi=150, facecolor=BG, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


# ══════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()
    os.makedirs(PRESENTATION_DIR, exist_ok=True)

    print("=" * 60)
    print("CELL & NUCLEUS BOUNDARY AREAS BY CELL TYPE")
    print("=" * 60)

    # Load all samples
    all_dfs = []
    samples = sorted(SAMPLE_TO_DX.keys())
    for i, sample_id in enumerate(samples):
        print(f"\n[{i+1}/{len(samples)}] {sample_id}:", flush=True)
        sample_df = load_sample_areas(sample_id)
        if len(sample_df) > 0:
            all_dfs.append(sample_df)

    df = pd.concat(all_dfs, ignore_index=True)
    print(f"\n  Total: {len(df):,} cells across {len(all_dfs)} samples")

    # Quick summary
    print(f"\n  Overall cell area: median={df['cell_area'].median():.0f} µm², "
          f"mean={df['cell_area'].mean():.0f} µm²")
    print(f"  Overall nucleus area: median={df['nucleus_area'].median():.0f} µm², "
          f"mean={df['nucleus_area'].mean():.0f} µm²")

    # ── Figures ──
    print("\n" + "=" * 60)
    print("PLOTTING")
    print("=" * 60)

    plot_subclass_violins(df)

    for cell_class, label in [
        ('Glutamatergic', 'Glutamatergic'),
        ('GABAergic', 'GABAergic'),
        ('Non-neuronal', 'Non-neuronal'),
    ]:
        plot_supertype_violins(df, cell_class, label)

    # ── CSVs ──
    save_summary_csvs(df)

    # ── Mixed effects models ──
    print("\n" + "=" * 60)
    print("MIXED EFFECTS MODELS")
    print("=" * 60)

    # Exclude WM-outlier sample (EXCLUDE_SAMPLES from config)
    df_model = df[~df['sample_id'].isin(EXCLUDE_SAMPLES)].copy()
    print(f"  Modeling on {len(df_model):,} cells "
          f"({df_model['sample_id'].nunique()} samples) "
          f"after excluding {EXCLUDE_SAMPLES}")

    # ── Subclass-level LME ──
    lme_subclass = fit_lme_per_celltype(
        df_model, level_col='subclass_label', min_cells_per_group=100)

    print_lme_summary(lme_subclass, 'subclass')
    plot_lme_forest(lme_subclass, level='subclass')
    save_lme_results(lme_subclass, level='subclass')

    # ── Supertype-level LME ──
    # Filter to supertypes with median ≥ 20 cells per sample
    lme_supertype = fit_lme_per_celltype(
        df_model, level_col='supertype_label',
        min_cells_per_group=100,
        min_median_cells_per_sample=20)

    print_lme_summary(lme_supertype, 'supertype')
    save_lme_results(lme_supertype, level='supertype')

    # Supertype forest plots — split by broad class for readability
    if not lme_supertype.empty:
        # All supertypes in one big plot
        plot_lme_forest(lme_supertype, level='supertype')

        # Per-class plots
        for class_name, class_tag in [('Glutamatergic', '_glut'),
                                       ('GABAergic', '_gaba'),
                                       ('Non-neuronal', '_nn')]:
            class_types = lme_supertype[
                lme_supertype['celltype'].apply(
                    lambda ct: _infer_broad_class(ct) == class_name
                )
            ]
            if not class_types.empty and class_types['celltype'].nunique() >= 2:
                plot_lme_forest(class_types, level='supertype',
                               suffix=class_tag)

    # ── Sample-level means (sanity check) ──
    print("\n" + "=" * 60)
    print("SAMPLE-LEVEL MEANS")
    print("=" * 60)
    sample_means = compute_sample_level_means(df_model)
    plot_sample_level_means(sample_means)

    total = time.time() - t_start
    print(f"\n{'=' * 60}")
    print(f"DONE — {total:.0f}s ({total/60:.1f} min)")
    print(f"{'=' * 60}")


if __name__ == '__main__':
    main()
