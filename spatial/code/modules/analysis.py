"""
Depth-stratified cell type proportion analysis for SCZ vs Control comparison.

Computes cell type proportions within cortical depth strata, runs
Mann-Whitney U tests with FDR correction, and validates against
MERFISH reference proportions.

IMPORTANT CAVEATS (from validation analysis):
- Sst is ~5-10x over-represented in Xenium relative to MERFISH reference,
  likely due to label transfer artifacts. Any Sst findings should be
  interpreted with extreme caution.
- L6 IT, L5 IT, Microglia-PVM are systematically under-represented.
- Oligodendrocyte is over-represented in deep layers.
- See xenium_vs_merfish_proportions.csv for the full comparison.
"""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from depth_model import DEPTH_STRATA, LAYER_BINS


def compute_proportions(adata, depth_col='predicted_norm_depth',
                        celltype_col='subclass_label',
                        sample_col='sample_id',
                        diagnosis_map=None,
                        strata=None, min_cells=100):
    """
    Compute cell type proportions per sample within each depth stratum.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated dataset with depth predictions and cell type labels.
    depth_col : str
        Column in .obs with predicted normalized depth.
    celltype_col : str
        Column in .obs with cell type labels.
    sample_col : str
        Column in .obs with sample IDs.
    diagnosis_map : dict, optional
        {sample_id: diagnosis} mapping.
    strata : dict, optional
        {stratum_name: (lower, upper)}. Defaults to DEPTH_STRATA.
    min_cells : int
        Minimum cells in a stratum to include.

    Returns
    -------
    pd.DataFrame
        Columns: sample, diagnosis, stratum, celltype, proportion,
        count, n_stratum.
    """
    if strata is None:
        strata = DEPTH_STRATA

    samples = sorted(adata.obs[sample_col].unique())
    celltypes = sorted(adata.obs[celltype_col].unique())

    records = []
    for sid in samples:
        mask = adata.obs[sample_col] == sid
        depths = adata.obs.loc[mask, depth_col].values.astype(float)
        labels = adata.obs.loc[mask, celltype_col].values.astype(str)
        dx = diagnosis_map.get(sid, 'Unknown') if diagnosis_map else 'Unknown'

        for sname, (lo, hi) in strata.items():
            stratum_mask = (depths >= lo) & (depths < hi)
            n_stratum = stratum_mask.sum()
            if n_stratum < min_cells:
                continue

            labels_in_stratum = labels[stratum_mask]
            counts = pd.Series(labels_in_stratum).value_counts()

            for ct in celltypes:
                records.append({
                    'sample': sid,
                    'diagnosis': dx,
                    'stratum': sname,
                    'celltype': ct,
                    'proportion': counts.get(ct, 0) / n_stratum,
                    'count': counts.get(ct, 0),
                    'n_stratum': n_stratum,
                })

    return pd.DataFrame(records)


def test_case_control(proportions_df, case_label='SCZ',
                      control_label='Control'):
    """
    Run Mann-Whitney U tests comparing case vs control proportions
    for each cell type × depth stratum combination.

    Parameters
    ----------
    proportions_df : pd.DataFrame
        Output from compute_proportions().
    case_label, control_label : str
        Diagnosis labels.

    Returns
    -------
    pd.DataFrame
        Test results with columns: stratum, celltype, scz_mean, ctrl_mean,
        log2fc, pval, padj, n_scz, n_ctrl.
    """
    results = []
    strata = proportions_df['stratum'].unique()
    celltypes = proportions_df['celltype'].unique()

    for stratum in strata:
        for ct in celltypes:
            mask = ((proportions_df['stratum'] == stratum) &
                    (proportions_df['celltype'] == ct))
            sub = proportions_df[mask]

            scz = sub[sub['diagnosis'] == case_label]['proportion'].values
            ctrl = sub[sub['diagnosis'] == control_label]['proportion'].values

            if len(scz) < 3 or len(ctrl) < 3:
                continue
            # Skip very rare types
            if scz.mean() < 0.001 and ctrl.mean() < 0.001:
                continue

            stat, pval = mannwhitneyu(scz, ctrl, alternative='two-sided')
            results.append({
                'stratum': stratum,
                'celltype': ct,
                'scz_mean': scz.mean(),
                'ctrl_mean': ctrl.mean(),
                'log2fc': np.log2((scz.mean() + 1e-6) / (ctrl.mean() + 1e-6)),
                'U_stat': stat,
                'pval': pval,
                'n_scz': len(scz),
                'n_ctrl': len(ctrl),
            })

    df = pd.DataFrame(results)
    if len(df) > 0:
        _, padj, _, _ = multipletests(df['pval'].values, method='fdr_bh')
        df['padj'] = padj
    return df


def flag_outlier_samples(adata, depth_col='predicted_norm_depth',
                         sample_col='sample_id',
                         merfish_depth=None,
                         min_layer_frac=0.05,
                         max_wm_frac=0.35):
    """
    Flag samples with missing cortical layers or extreme depth distributions.

    Criteria:
    1. Any cortical layer (L2/3, L4, L5, L6) with < min_layer_frac of cells
    2. WM fraction > max_wm_frac
    3. Wasserstein distance to MERFISH > 2x median (if merfish_depth provided)

    Parameters
    ----------
    adata : anndata.AnnData
        Dataset with predicted depth.
    depth_col : str
        Depth column name.
    sample_col : str
        Sample ID column name.
    merfish_depth : np.ndarray, optional
        MERFISH reference depth values for comparison.
    min_layer_frac : float
        Minimum fraction threshold for cortical layers.
    max_wm_frac : float
        Maximum WM fraction before flagging.

    Returns
    -------
    dict
        {sample_id: [list of flag reasons]} for flagged samples.
        Samples not flagged are not included.
    """
    from scipy.stats import wasserstein_distance

    cortical_layers = ['L2/3', 'L4', 'L5', 'L6']
    samples = sorted(adata.obs[sample_col].unique())

    # Compute per-sample metrics
    w_dists = {}
    layer_fracs = {}
    for sid in samples:
        mask = adata.obs[sample_col] == sid
        d = adata.obs.loc[mask, depth_col].values.astype(float)

        # Layer fractions
        fracs = {}
        for lname, (lo, hi) in LAYER_BINS.items():
            fracs[lname] = ((d >= lo) & (d < hi)).sum() / len(d)
        layer_fracs[sid] = fracs

        if merfish_depth is not None:
            w_dists[sid] = wasserstein_distance(d, merfish_depth)

    # Flag
    median_w = np.median(list(w_dists.values())) if w_dists else 0
    flags = {}
    for sid in samples:
        reasons = []
        for l in cortical_layers:
            if layer_fracs[sid][l] < min_layer_frac:
                reasons.append(
                    f"{l} underrepresented ({layer_fracs[sid][l]:.1%})")
        if layer_fracs[sid].get('WM', 0) > max_wm_frac:
            reasons.append(f"excess WM ({layer_fracs[sid]['WM']:.1%})")
        if w_dists and w_dists[sid] > 2 * median_w:
            reasons.append(f"high Wasserstein ({w_dists[sid]:.3f})")
        if reasons:
            flags[sid] = reasons

    return flags


def validate_against_merfish(xenium_props, merfish_adata,
                             strata=None, depth_col='Normalized depth from pia',
                             subclass_col='Subclass',
                             donor_col='Donor ID'):
    """
    Compare Xenium cell type proportions against MERFISH reference.

    Computes per-donor proportions in MERFISH for each depth stratum,
    then compares mean proportions to identify systematic biases
    in the Xenium label transfer.

    Parameters
    ----------
    xenium_props : pd.DataFrame
        Output from compute_proportions() for Xenium controls.
    merfish_adata : anndata.AnnData
        SEA-AD MERFISH dataset.
    strata : dict, optional
        Depth strata. Defaults to DEPTH_STRATA.
    depth_col : str
        Depth column in MERFISH .obs.
    subclass_col : str
        Subclass column in MERFISH .obs.
    donor_col : str
        Donor column in MERFISH .obs.

    Returns
    -------
    pd.DataFrame
        Comparison table with columns: stratum, celltype, merfish_mean,
        xenium_mean, diff, log2_ratio, abs_diff.
    """
    if strata is None:
        strata = DEPTH_STRATA

    depth = merfish_adata.obs[depth_col].values.astype(float)
    has_depth = ~np.isnan(depth)
    subclass = merfish_adata.obs[subclass_col].values.astype(str)
    donors = merfish_adata.obs[donor_col].values.astype(str)
    all_types = sorted(set(subclass))

    # MERFISH per-donor proportions
    merfish_records = []
    for donor in np.unique(donors[has_depth]):
        d_mask = (donors == donor) & has_depth
        d_depth = depth[d_mask]
        d_sub = subclass[d_mask]

        for sname, (lo, hi) in strata.items():
            s_mask = (d_depth >= lo) & (d_depth < hi)
            n = s_mask.sum()
            if n < 50:
                continue
            counts = pd.Series(d_sub[s_mask]).value_counts()
            for ct in all_types:
                merfish_records.append({
                    'donor': donor, 'stratum': sname,
                    'celltype': ct,
                    'proportion': counts.get(ct, 0) / n,
                })

    df_mer = pd.DataFrame(merfish_records)
    mer_means = (df_mer.groupby(['stratum', 'celltype'])['proportion']
                 .mean().reset_index()
                 .rename(columns={'proportion': 'merfish_mean'}))

    # Xenium means
    xen_means = (xenium_props.groupby(['stratum', 'celltype'])['proportion']
                 .mean().reset_index()
                 .rename(columns={'proportion': 'xenium_mean'}))

    # Merge
    comp = pd.merge(mer_means, xen_means, on=['stratum', 'celltype'], how='outer')
    comp = comp.dropna(subset=['merfish_mean', 'xenium_mean'])
    comp['diff'] = comp['xenium_mean'] - comp['merfish_mean']
    comp['abs_diff'] = comp['diff'].abs()
    comp['log2_ratio'] = np.log2(
        (comp['xenium_mean'] + 1e-5) / (comp['merfish_mean'] + 1e-5))

    return comp.sort_values('abs_diff', ascending=False)
