import os
import sys
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR, MERFISH_PATH, CONTROL_SAMPLES, SUBCLASS_CONF_THRESH

controls = CONTROL_SAMPLES

depth_strata = {
    'L2/3': (0.10, 0.30),
    'L4': (0.30, 0.45),
    'L5': (0.45, 0.65),
    'L6': (0.65, 0.85),
}

# ===================================================================
# 1. Compute MERFISH reference proportions per depth stratum
# ===================================================================
print("Loading MERFISH data...")
merfish = ad.read_h5ad(MERFISH_PATH)

# Use normalized depth from pia for MERFISH (ground truth depth)
merfish_depth = merfish.obs['Normalized depth from pia'].values.astype(float)
has_depth = ~np.isnan(merfish_depth)
merfish_sub = merfish.obs['Subclass'].values.astype(str)

# Also get MERFISH donors (SEA-AD are all controls)
merfish_donors = merfish.obs['Donor ID'].values.astype(str) if 'Donor ID' in merfish.obs.columns else merfish.obs['Specimen Barcode'].values.astype(str)

# MERFISH per-donor proportions (only depth-annotated cells)
merfish_records = []
unique_donors = np.unique(merfish_donors[has_depth])
print(f"MERFISH donors with depth annotations: {len(unique_donors)}")

all_merfish_subclasses = sorted(set(merfish_sub))
print(f"MERFISH subclasses: {len(all_merfish_subclasses)}")

for donor in unique_donors:
    donor_mask = (merfish_donors == donor) & has_depth
    d = merfish_depth[donor_mask]
    s = merfish_sub[donor_mask]
    
    for stratum_name, (lo, hi) in depth_strata.items():
        stratum_mask = (d >= lo) & (d < hi)
        n_stratum = stratum_mask.sum()
        if n_stratum < 50:
            continue
        
        s_in_stratum = s[stratum_mask]
        counts = pd.Series(s_in_stratum).value_counts()
        
        for sc in all_merfish_subclasses:
            prop = counts.get(sc, 0) / n_stratum
            merfish_records.append({
                'donor': donor, 'stratum': stratum_name,
                'celltype': sc, 'proportion': prop,
                'count': counts.get(sc, 0), 'n_stratum': n_stratum,
                'source': 'MERFISH'
            })

df_merfish = pd.DataFrame(merfish_records)
print(f"MERFISH records: {len(df_merfish)}")

# MERFISH mean proportions per stratum (average across donors)
merfish_means = df_merfish.groupby(['stratum', 'celltype'])['proportion'].agg(['mean', 'std', 'count']).reset_index()
merfish_means.columns = ['stratum', 'celltype', 'merfish_mean', 'merfish_std', 'merfish_n_donors']

del merfish  # free memory

# ===================================================================
# 2. Compute Xenium control proportions per depth stratum
# ===================================================================
print("\nLoading Xenium data...")
xenium = ad.read_h5ad(os.path.join(BASE_DIR, 'output', 'all_samples_annotated.h5ad'))

# Apply bottom-1% subclass confidence filter
if 'subclass_label_confidence' in xenium.obs.columns:
    conf = xenium.obs['subclass_label_confidence'].astype(float)
    conf_mask = conf >= SUBCLASS_CONF_THRESH
    n_drop = (~conf_mask).sum()
    xenium = xenium[conf_mask].copy()
    print(f"Confidence filter: dropped {n_drop:,} cells, {xenium.shape[0]:,} remaining")

# Keep only controls
ctrl_mask = xenium.obs['sample_id'].isin(controls)
xenium_ctrl = xenium[ctrl_mask].copy()
print(f"Xenium control cells: {xenium_ctrl.shape[0]:,} from {len(controls)} samples")

all_xenium_subclasses = sorted(xenium_ctrl.obs['subclass_label'].unique())
print(f"Xenium subclasses: {len(all_xenium_subclasses)}")

xenium_records = []
for sid in controls:
    sample_mask = xenium_ctrl.obs['sample_id'] == sid
    depths = xenium_ctrl.obs.loc[sample_mask, 'predicted_norm_depth'].values.astype(float)
    sub_labels = xenium_ctrl.obs.loc[sample_mask, 'subclass_label'].values.astype(str)
    
    for stratum_name, (lo, hi) in depth_strata.items():
        stratum_mask = (depths >= lo) & (depths < hi)
        n_stratum = stratum_mask.sum()
        if n_stratum < 50:
            continue
        
        s_in_stratum = sub_labels[stratum_mask]
        counts = pd.Series(s_in_stratum).value_counts()
        
        for sc in all_xenium_subclasses:
            prop = counts.get(sc, 0) / n_stratum
            xenium_records.append({
                'donor': sid, 'stratum': stratum_name,
                'celltype': sc, 'proportion': prop,
                'count': counts.get(sc, 0), 'n_stratum': n_stratum,
                'source': 'Xenium'
            })

df_xenium = pd.DataFrame(xenium_records)
print(f"Xenium records: {len(df_xenium)}")

# Xenium mean proportions
xenium_means = df_xenium.groupby(['stratum', 'celltype'])['proportion'].agg(['mean', 'std', 'count']).reset_index()
xenium_means.columns = ['stratum', 'celltype', 'xenium_mean', 'xenium_std', 'xenium_n_donors']

# ===================================================================
# 3. Merge and compare
# ===================================================================
# Find common cell types
common_types = sorted(set(all_merfish_subclasses) & set(all_xenium_subclasses))
merfish_only = sorted(set(all_merfish_subclasses) - set(all_xenium_subclasses))
xenium_only = sorted(set(all_xenium_subclasses) - set(all_merfish_subclasses))

print(f"\nCommon subclasses: {len(common_types)}")
print(f"MERFISH-only: {merfish_only}")
print(f"Xenium-only: {xenium_only}")

# Merge on stratum + celltype
comparison = pd.merge(merfish_means, xenium_means, on=['stratum', 'celltype'], how='outer')
comparison = comparison[comparison['celltype'].isin(common_types)]
comparison = comparison.dropna(subset=['merfish_mean', 'xenium_mean'])

# Compute deviation metrics
comparison['diff'] = comparison['xenium_mean'] - comparison['merfish_mean']
comparison['abs_diff'] = comparison['diff'].abs()
comparison['log2_ratio'] = np.log2((comparison['xenium_mean'] + 1e-5) / (comparison['merfish_mean'] + 1e-5))
comparison['abs_log2_ratio'] = comparison['log2_ratio'].abs()

# Print summary
print("\n" + "="*80)
print("COMPARISON: Xenium Control vs MERFISH Reference")
print("="*80)

for st in depth_strata:
    sub = comparison[comparison['stratum'] == st].sort_values('abs_diff', ascending=False)
    r_val, _ = pearsonr(sub['merfish_mean'], sub['xenium_mean'])
    rho, _ = spearmanr(sub['merfish_mean'], sub['xenium_mean'])
    print(f"\n--- {st} --- (Pearson r={r_val:.3f}, Spearman ρ={rho:.3f})")
    print(f"{'Celltype':<25} {'MERFISH':>8} {'Xenium':>8} {'Diff':>8} {'log2R':>8}")
    for _, r in sub.head(10).iterrows():
        flag = " ***" if r['abs_diff'] > 0.03 else (" **" if r['abs_diff'] > 0.02 else (" *" if r['abs_diff'] > 0.01 else ""))
        print(f"{r['celltype']:<25} {r['merfish_mean']:>8.4f} {r['xenium_mean']:>8.4f} {r['diff']:>+8.4f} {r['log2_ratio']:>+8.2f}{flag}")

# Overall biggest deviations
print("\n" + "="*80)
print("TOP 20 LARGEST DEVIATIONS (absolute difference)")
print("="*80)
top_dev = comparison.sort_values('abs_diff', ascending=False).head(20)
print(f"{'Stratum':<8} {'Celltype':<25} {'MERFISH':>8} {'Xenium':>8} {'Diff':>8} {'log2R':>8}")
for _, r in top_dev.iterrows():
    print(f"{r['stratum']:<8} {r['celltype']:<25} {r['merfish_mean']:>8.4f} {r['xenium_mean']:>8.4f} {r['diff']:>+8.4f} {r['log2_ratio']:>+8.2f}")

# Save comparison
comparison.to_csv(os.path.join(BASE_DIR, 'output', 'xenium_vs_merfish_proportions.csv'), index=False)

# ===================================================================
# 4. FIGURES
# ===================================================================

# Figure 1: Scatter plots per stratum (MERFISH mean vs Xenium mean)
fig, axes = plt.subplots(1, 4, figsize=(36, 9), facecolor='white')

for j, st in enumerate(depth_strata):
    ax = axes[j]
    sub = comparison[comparison['stratum'] == st].copy()
    
    r_val, _ = pearsonr(sub['merfish_mean'], sub['xenium_mean'])
    rho, _ = spearmanr(sub['merfish_mean'], sub['xenium_mean'])
    
    # Plot identity line
    max_val = max(sub['merfish_mean'].max(), sub['xenium_mean'].max()) * 1.1
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.4, linewidth=2)
    
    # Scatter with labels for outliers
    for _, row in sub.iterrows():
        color = 'red' if row['abs_diff'] > 0.015 else 'steelblue'
        size = 120 if row['abs_diff'] > 0.015 else 60
        ax.scatter(row['merfish_mean'], row['xenium_mean'], s=size, color=color,
                  alpha=0.7, edgecolor='black', linewidth=0.5, zorder=5)
        if row['abs_diff'] > 0.015:
            ax.annotate(row['celltype'], (row['merfish_mean'], row['xenium_mean']),
                       fontsize=11, fontweight='bold', color='red',
                       xytext=(5, 5), textcoords='offset points')
    
    ax.set_xlabel('MERFISH proportion (SEA-AD reference)', fontsize=16)
    ax.set_ylabel('Xenium proportion (Control samples)', fontsize=16)
    ax.set_title(f'{st}\nPearson r = {r_val:.3f} | Spearman ρ = {rho:.3f}', fontsize=18, fontweight='bold')
    ax.tick_params(labelsize=13)
    ax.set_xlim(-0.005, max_val)
    ax.set_ylim(-0.005, max_val)

fig.suptitle('Cell Type Proportions: Xenium Controls vs MERFISH Reference\n'
             'Each point = one subclass (mean across donors per depth stratum)\n'
             'Red = deviation > 1.5% | Dashed = identity line',
             fontsize=22, fontweight='bold', y=1.04)
plt.tight_layout()
out1 = os.path.join(BASE_DIR, 'output', 'xenium_vs_merfish_scatter.png')
plt.savefig(out1, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"\nSaved: {out1}")

# Figure 2: Deviation heatmap (celltype × stratum)
fig, ax = plt.subplots(figsize=(14, 16), facecolor='white')

# Build matrix
celltypes_sorted = sorted(common_types)
diff_matrix = np.zeros((len(celltypes_sorted), len(list(depth_strata.keys()))))
strata_list = list(depth_strata.keys())

for i, ct in enumerate(celltypes_sorted):
    for j, st in enumerate(strata_list):
        row = comparison[(comparison['celltype'] == ct) & (comparison['stratum'] == st)]
        if len(row) == 1:
            diff_matrix[i, j] = row['diff'].values[0]

vmax = max(abs(diff_matrix.min()), abs(diff_matrix.max()), 0.05)
im = ax.imshow(diff_matrix, aspect='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)

for i in range(len(celltypes_sorted)):
    for j in range(len(strata_list)):
        val = diff_matrix[i, j]
        text_color = 'white' if abs(val) > vmax * 0.6 else 'black'
        star = '***' if abs(val) > 0.03 else ('**' if abs(val) > 0.02 else ('*' if abs(val) > 0.01 else ''))
        ax.text(j, i, f'{val:+.3f}\n{star}' if star else f'{val:+.3f}',
               ha='center', va='center', fontsize=11,
               fontweight='bold' if star else 'normal', color=text_color)

ax.set_xticks(range(len(strata_list)))
ax.set_xticklabels(strata_list, fontsize=16, fontweight='bold')
ax.set_yticks(range(len(celltypes_sorted)))
ax.set_yticklabels(celltypes_sorted, fontsize=13)
ax.set_xlabel('Cortical Depth Stratum', fontsize=18, fontweight='bold')

cbar = plt.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
cbar.set_label('Xenium − MERFISH (proportion difference)', fontsize=14)
cbar.ax.tick_params(labelsize=12)

ax.set_title('Xenium vs MERFISH: Proportion Deviation per Subclass & Depth\n'
             '(Xenium Controls mean − MERFISH mean)\n'
             '* > 1%, ** > 2%, *** > 3% absolute deviation',
             fontsize=18, fontweight='bold', pad=20)
plt.tight_layout()
out2 = os.path.join(BASE_DIR, 'output', 'xenium_vs_merfish_heatmap.png')
plt.savefig(out2, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: {out2}")

# Figure 3: Boxplot comparison for biggest outlier cell types
# Show per-donor distributions side by side
top_outliers = comparison.groupby('celltype')['abs_diff'].max().sort_values(ascending=False).head(8).index.tolist()
print(f"\nTop 8 outlier cell types: {top_outliers}")

fig, axes = plt.subplots(len(top_outliers), len(strata_list), figsize=(28, 5 * len(top_outliers)), facecolor='white')

for i, ct in enumerate(top_outliers):
    for j, st in enumerate(strata_list):
        ax = axes[i, j]
        
        # MERFISH per-donor
        merfish_vals = df_merfish[(df_merfish['celltype'] == ct) & (df_merfish['stratum'] == st)]['proportion'].values
        # Xenium per-donor (controls only)
        xenium_vals = df_xenium[(df_xenium['celltype'] == ct) & (df_xenium['stratum'] == st)]['proportion'].values
        
        if len(merfish_vals) > 0 and len(xenium_vals) > 0:
            bp = ax.boxplot([merfish_vals, xenium_vals], labels=['MERFISH\n(SEA-AD)', 'Xenium\n(Control)'],
                           patch_artist=True, widths=0.6,
                           medianprops=dict(color='black', linewidth=2))
            bp['boxes'][0].set_facecolor('#2CA02C'); bp['boxes'][0].set_alpha(0.5)
            bp['boxes'][1].set_facecolor('#FF7F0E'); bp['boxes'][1].set_alpha(0.5)
            
            ax.scatter(np.ones(len(merfish_vals)) + np.random.normal(0, 0.04, len(merfish_vals)),
                      merfish_vals, color='#1B7A1B', s=35, alpha=0.7, zorder=5)
            ax.scatter(np.ones(len(xenium_vals)) * 2 + np.random.normal(0, 0.04, len(xenium_vals)),
                      xenium_vals, color='#CC6600', s=35, alpha=0.7, zorder=5)
            
            # Annotation
            diff = xenium_vals.mean() - merfish_vals.mean()
            ax.set_title(f'{st}\nΔ = {diff:+.4f}', fontsize=14,
                        fontweight='bold' if abs(diff) > 0.015 else 'normal',
                        color='red' if abs(diff) > 0.015 else 'black')
        
        if j == 0:
            ax.set_ylabel(f'{ct}\nproportion', fontsize=14, fontweight='bold')
        ax.tick_params(labelsize=11)
        ax.set_xlim(0.3, 2.7)

fig.suptitle('Per-Donor Cell Type Proportions: MERFISH (SEA-AD) vs Xenium (Controls)\n'
             'Top 8 most deviant cell types across depth strata',
             fontsize=22, fontweight='bold', y=1.01)
plt.tight_layout()
out3 = os.path.join(BASE_DIR, 'output', 'xenium_vs_merfish_boxplots.png')
plt.savefig(out3, dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: {out3}")

print("\nDone!")
