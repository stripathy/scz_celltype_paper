import os
import sys
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import SAMPLE_TO_DX, BASE_DIR, SUBCLASS_CONF_THRESH

sample_to_dx = SAMPLE_TO_DX
OUTPUT_DIR = os.path.join(BASE_DIR, "output")

# Depth strata
depth_strata = {
    'L2/3': (0.10, 0.30),
    'L4': (0.30, 0.45),
    'L5': (0.45, 0.65),
    'L6': (0.65, 0.85),
}

# Load data
print("Loading Xenium data...")
xenium = ad.read_h5ad(os.path.join(OUTPUT_DIR, 'all_samples_annotated.h5ad'))
print(f"Loaded: {xenium.shape[0]:,} cells")

# Remove Br2039 outlier
outlier = 'Br2039'
mask_keep = xenium.obs['sample_id'] != outlier
xenium = xenium[mask_keep].copy()
print(f"After removing {outlier}: {xenium.shape[0]:,} cells")

# Apply bottom-1% subclass confidence filter
if 'subclass_label_confidence' in xenium.obs.columns:
    conf = xenium.obs['subclass_label_confidence'].astype(float)
    conf_mask = conf >= SUBCLASS_CONF_THRESH
    n_drop = (~conf_mask).sum()
    xenium = xenium[conf_mask].copy()
    print(f"Confidence filter (subclass >= {SUBCLASS_CONF_THRESH}): "
          f"dropped {n_drop:,} cells, {xenium.shape[0]:,} remaining")

samples = sorted(xenium.obs['sample_id'].unique())
print(f"Samples: {len(samples)} (12 Control, 11 SCZ)")

# ===================================================================
# 1. SUBCLASS-LEVEL ANALYSIS
# ===================================================================
print("\n" + "="*80)
print("SUBCLASS-LEVEL DEPTH-STRATIFIED ANALYSIS")
print("="*80)

subclasses = sorted(xenium.obs['subclass_label'].unique())
print(f"Subclass types: {len(subclasses)}")

# Compute proportions: for each sample × stratum, proportion of each subclass
# among all cells in that stratum
all_subclass_props = []

for sid in samples:
    sample_mask = xenium.obs['sample_id'] == sid
    depths = xenium.obs.loc[sample_mask, 'predicted_norm_depth'].values.astype(float)
    sub_labels = xenium.obs.loc[sample_mask, 'subclass_label'].values.astype(str)
    dx = sample_to_dx[sid]
    
    for stratum_name, (lo, hi) in depth_strata.items():
        stratum_mask = (depths >= lo) & (depths < hi)
        n_stratum = stratum_mask.sum()
        
        if n_stratum < 100:
            continue
        
        sub_in_stratum = sub_labels[stratum_mask]
        counts = pd.Series(sub_in_stratum).value_counts()
        
        for sc in subclasses:
            prop = counts.get(sc, 0) / n_stratum
            all_subclass_props.append({
                'sample': sid, 'diagnosis': dx, 'stratum': stratum_name,
                'celltype': sc, 'proportion': prop, 'count': counts.get(sc, 0),
                'n_stratum': n_stratum, 'level': 'subclass'
            })

df_sub = pd.DataFrame(all_subclass_props)

# Statistical tests: for each subclass × stratum, Mann-Whitney SCZ vs Control
test_results_sub = []
for stratum_name in depth_strata:
    for sc in subclasses:
        mask = (df_sub['stratum'] == stratum_name) & (df_sub['celltype'] == sc)
        sub_df = df_sub[mask]
        
        scz_vals = sub_df[sub_df['diagnosis'] == 'SCZ']['proportion'].values
        ctrl_vals = sub_df[sub_df['diagnosis'] == 'Control']['proportion'].values
        
        if len(scz_vals) >= 3 and len(ctrl_vals) >= 3:
            stat, pval = mannwhitneyu(scz_vals, ctrl_vals, alternative='two-sided')
            test_results_sub.append({
                'stratum': stratum_name, 'celltype': sc,
                'scz_mean': scz_vals.mean(), 'ctrl_mean': ctrl_vals.mean(),
                'scz_median': np.median(scz_vals), 'ctrl_median': np.median(ctrl_vals),
                'log2fc': np.log2((scz_vals.mean() + 1e-6) / (ctrl_vals.mean() + 1e-6)),
                'U_stat': stat, 'pval': pval,
                'n_scz': len(scz_vals), 'n_ctrl': len(ctrl_vals),
                'level': 'subclass'
            })

df_tests_sub = pd.DataFrame(test_results_sub)

# FDR correction across ALL tests
_, padj, _, _ = multipletests(df_tests_sub['pval'].values, method='fdr_bh')
df_tests_sub['padj'] = padj

# Print significant results
print(f"\nTotal tests: {len(df_tests_sub)}")
sig = df_tests_sub[df_tests_sub['padj'] < 0.05].sort_values('padj')
print(f"FDR-significant (padj < 0.05): {len(sig)}")
if len(sig) > 0:
    print(f"\n{'Stratum':<8} {'Subclass':<25} {'SCZ_mean':>9} {'Ctrl_mean':>10} {'log2FC':>8} {'pval':>10} {'padj':>10}")
    for _, r in sig.iterrows():
        print(f"{r['stratum']:<8} {r['celltype']:<25} {r['scz_mean']:>9.4f} {r['ctrl_mean']:>10.4f} "
              f"{r['log2fc']:>8.2f} {r['pval']:>10.4f} {r['padj']:>10.4f}")

# Nominal significance
nom = df_tests_sub[df_tests_sub['pval'] < 0.05].sort_values('pval')
print(f"\nNominally significant (p < 0.05): {len(nom)}")
if len(nom) > 0:
    print(f"\n{'Stratum':<8} {'Subclass':<25} {'SCZ_mean':>9} {'Ctrl_mean':>10} {'log2FC':>8} {'pval':>10} {'padj':>10}")
    for _, r in nom.iterrows():
        direction = "↑SCZ" if r['log2fc'] > 0 else "↓SCZ"
        print(f"{r['stratum']:<8} {r['celltype']:<25} {r['scz_mean']:>9.4f} {r['ctrl_mean']:>10.4f} "
              f"{r['log2fc']:>8.2f} {r['pval']:>10.4f} {r['padj']:>10.4f}  {direction}")

# ===================================================================
# 2. CLUSTER-LEVEL ANALYSIS
# ===================================================================
print("\n" + "="*80)
print("SUPERTYPE-LEVEL DEPTH-STRATIFIED ANALYSIS")
print("="*80)

clusters = sorted(xenium.obs['supertype_label'].unique())
print(f"Cluster types: {len(clusters)}")

all_cluster_props = []
for sid in samples:
    sample_mask = xenium.obs['sample_id'] == sid
    depths = xenium.obs.loc[sample_mask, 'predicted_norm_depth'].values.astype(float)
    clust_labels = xenium.obs.loc[sample_mask, 'supertype_label'].values.astype(str)
    dx = sample_to_dx[sid]
    
    for stratum_name, (lo, hi) in depth_strata.items():
        stratum_mask = (depths >= lo) & (depths < hi)
        n_stratum = stratum_mask.sum()
        
        if n_stratum < 100:
            continue
        
        cl_in_stratum = clust_labels[stratum_mask]
        counts = pd.Series(cl_in_stratum).value_counts()
        
        for cl in clusters:
            prop = counts.get(cl, 0) / n_stratum
            all_cluster_props.append({
                'sample': sid, 'diagnosis': dx, 'stratum': stratum_name,
                'celltype': cl, 'proportion': prop, 'count': counts.get(cl, 0),
                'n_stratum': n_stratum, 'level': 'cluster'
            })

df_clust = pd.DataFrame(all_cluster_props)

test_results_clust = []
for stratum_name in depth_strata:
    for cl in clusters:
        mask = (df_clust['stratum'] == stratum_name) & (df_clust['celltype'] == cl)
        sub_df = df_clust[mask]
        
        scz_vals = sub_df[sub_df['diagnosis'] == 'SCZ']['proportion'].values
        ctrl_vals = sub_df[sub_df['diagnosis'] == 'Control']['proportion'].values
        
        if len(scz_vals) >= 3 and len(ctrl_vals) >= 3:
            # Skip cell types that are essentially absent in this stratum
            if scz_vals.mean() < 0.001 and ctrl_vals.mean() < 0.001:
                continue
            stat, pval = mannwhitneyu(scz_vals, ctrl_vals, alternative='two-sided')
            test_results_clust.append({
                'stratum': stratum_name, 'celltype': cl,
                'scz_mean': scz_vals.mean(), 'ctrl_mean': ctrl_vals.mean(),
                'scz_median': np.median(scz_vals), 'ctrl_median': np.median(ctrl_vals),
                'log2fc': np.log2((scz_vals.mean() + 1e-6) / (ctrl_vals.mean() + 1e-6)),
                'U_stat': stat, 'pval': pval,
                'n_scz': len(scz_vals), 'n_ctrl': len(ctrl_vals),
                'level': 'cluster'
            })

df_tests_clust = pd.DataFrame(test_results_clust)
_, padj_c, _, _ = multipletests(df_tests_clust['pval'].values, method='fdr_bh')
df_tests_clust['padj'] = padj_c

sig_c = df_tests_clust[df_tests_clust['padj'] < 0.05].sort_values('padj')
print(f"\nTotal tests: {len(df_tests_clust)}")
print(f"FDR-significant (padj < 0.05): {len(sig_c)}")
if len(sig_c) > 0:
    print(f"\n{'Stratum':<8} {'Cluster':<35} {'SCZ_mean':>9} {'Ctrl_mean':>10} {'log2FC':>8} {'pval':>10} {'padj':>10}")
    for _, r in sig_c.iterrows():
        print(f"{r['stratum']:<8} {r['celltype']:<35} {r['scz_mean']:>9.4f} {r['ctrl_mean']:>10.4f} "
              f"{r['log2fc']:>8.2f} {r['pval']:>10.4f} {r['padj']:>10.4f}")

nom_c = df_tests_clust[df_tests_clust['pval'] < 0.05].sort_values('pval')
print(f"\nNominally significant (p < 0.05): {len(nom_c)}")
if len(nom_c) > 0:
    print(f"\n{'Stratum':<8} {'Cluster':<35} {'SCZ_mean':>9} {'Ctrl_mean':>10} {'log2FC':>8} {'pval':>10} {'padj':>10}")
    for _, r in nom_c.head(30).iterrows():
        direction = "↑SCZ" if r['log2fc'] > 0 else "↓SCZ"
        print(f"{r['stratum']:<8} {r['celltype']:<35} {r['scz_mean']:>9.4f} {r['ctrl_mean']:>10.4f} "
              f"{r['log2fc']:>8.2f} {r['pval']:>10.4f} {r['padj']:>10.4f}  {direction}")

# ===================================================================
# 3. SST-FOCUSED ANALYSIS
# ===================================================================
print("\n" + "="*80)
print("SST-FOCUSED ANALYSIS (all depth strata)")
print("="*80)

sst_clusters = [c for c in clusters if 'Sst' in c or 'SST' in c]
print(f"SST clusters: {sst_clusters}")

# For SST, also compute proportion relative to all GABAergic neurons
print("\nSST subclass proportion (of all cells) per stratum:")
sst_sub_tests = df_tests_sub[df_tests_sub['celltype'] == 'Sst']
for _, r in sst_sub_tests.iterrows():
    direction = "↑SCZ" if r['log2fc'] > 0 else "↓SCZ"
    print(f"  {r['stratum']}: SCZ={r['scz_mean']:.4f}, Ctrl={r['ctrl_mean']:.4f}, "
          f"log2FC={r['log2fc']:.2f}, p={r['pval']:.4f}, padj={r['padj']:.4f} {direction}")

print("\nSST cluster-level results:")
sst_clust_tests = df_tests_clust[df_tests_clust['celltype'].str.contains('Sst', case=False, na=False)]
for _, r in sst_clust_tests.sort_values('pval').iterrows():
    direction = "↑SCZ" if r['log2fc'] > 0 else "↓SCZ"
    print(f"  {r['stratum']:<8} {r['celltype']:<30} SCZ={r['scz_mean']:.4f}, Ctrl={r['ctrl_mean']:.4f}, "
          f"log2FC={r['log2fc']:.2f}, p={r['pval']:.4f}, padj={r['padj']:.4f} {direction}")

# ===================================================================
# 4. SAVE ALL RESULTS
# ===================================================================
# Combine tests
all_tests = pd.concat([df_tests_sub, df_tests_clust], ignore_index=True)
all_tests.to_csv(os.path.join(OUTPUT_DIR, 'depth_stratified_tests.csv'), index=False)

# Save proportions
all_props = pd.concat([df_sub, df_clust], ignore_index=True)
all_props.to_csv(os.path.join(OUTPUT_DIR, 'depth_stratified_proportions.csv'), index=False)

print(f"\nSaved: depth_stratified_tests.csv ({len(all_tests)} tests)")
print(f"Saved: depth_stratified_proportions.csv ({len(all_props)} rows)")
print("\nDone!")
