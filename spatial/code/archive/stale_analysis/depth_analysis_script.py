#!/usr/bin/env python3
"""
Depth-from-pia analysis: compare Xenium predicted depth distributions to MERFISH
reference. Flags sections with missing layers and evaluates resampling feasibility.

Output:
  output/depth_resampling_analysis.csv
  output/depth_resampling_analysis.png
"""

import os
import sys
from collections import Counter
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp, wasserstein_distance

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR, MERFISH_PATH, SAMPLE_TO_DX, SUBCLASS_CONF_THRESH

OUTPUT_DIR = os.path.join(BASE_DIR, "output")
H5AD_PATH = os.path.join(OUTPUT_DIR, "all_samples_annotated.h5ad")

LAYER_BINS = {
    'L1': (-np.inf, 0.1), 'L2/3': (0.1, 0.3), 'L4': (0.3, 0.45),
    'L5': (0.45, 0.65), 'L6': (0.65, 0.85), 'WM': (0.85, np.inf),
}
CORTICAL_LAYERS = ['L1', 'L2/3', 'L4', 'L5', 'L6']

# ============================================================
# 1. Load MERFISH depth distribution (cropped/annotated cells)
# ============================================================
print("Loading MERFISH data...")
merfish = ad.read_h5ad(MERFISH_PATH)
merfish_depth = merfish.obs['Normalized depth from pia'].values.astype(float)
merfish_has_depth = ~np.isnan(merfish_depth)
merfish_depth_annotated = merfish_depth[merfish_has_depth]
print(f"MERFISH annotated cells: {len(merfish_depth_annotated):,}")
print(f"MERFISH depth range: {merfish_depth_annotated.min():.3f} to {merfish_depth_annotated.max():.3f}")
print(f"MERFISH depth mean: {merfish_depth_annotated.mean():.3f}, std: {merfish_depth_annotated.std():.3f}")

print("\nMERFISH layer distribution (annotated cells):")
merfish_layer_fracs = {}
for lname, (lo, hi) in LAYER_BINS.items():
    mask = (merfish_depth_annotated >= lo) & (merfish_depth_annotated < hi)
    frac = mask.sum() / len(merfish_depth_annotated)
    merfish_layer_fracs[lname] = frac
    print(f"  {lname}: {mask.sum():,} ({frac:.1%})")

del merfish  # free memory

# ============================================================
# 2. Load all Xenium data
# ============================================================
print("\nLoading combined Xenium data...")
xenium = ad.read_h5ad(H5AD_PATH)
print(f"Xenium: {xenium.shape[0]:,} cells")

# Bottom-1% subclass confidence filter
if 'subclass_label_confidence' in xenium.obs.columns:
    conf = xenium.obs['subclass_label_confidence'].astype(float)
    conf_mask = conf >= SUBCLASS_CONF_THRESH
    n_drop = (~conf_mask).sum()
    xenium = xenium[conf_mask].copy()
    print(f"Confidence filter: dropped {n_drop:,} cells, {xenium.shape[0]:,} remaining")

# ============================================================
# 3. Per-sample depth analysis
# ============================================================
print("\n" + "="*80)
print("PER-SAMPLE DEPTH ANALYSIS")
print("="*80)

samples = sorted(xenium.obs['sample_id'].unique())

results = []
for sid in samples:
    mask = xenium.obs['sample_id'] == sid
    depths = xenium.obs.loc[mask, 'predicted_norm_depth'].values.astype(float)
    n_cells = len(depths)
    dx = SAMPLE_TO_DX.get(sid, 'Unknown')

    # Layer fractions
    layer_fracs = {}
    layer_counts = {}
    for lname, (lo, hi) in LAYER_BINS.items():
        lmask = (depths >= lo) & (depths < hi)
        frac = lmask.sum() / n_cells
        layer_fracs[lname] = frac
        layer_counts[lname] = lmask.sum()

    # Compare to MERFISH distribution
    ks_stat, ks_p = ks_2samp(depths, merfish_depth_annotated)
    w_dist = wasserstein_distance(depths, merfish_depth_annotated)

    # Define "covers a layer" as having >= 2% of cells in that layer
    missing_layers = [l for l in LAYER_BINS if layer_fracs[l] < 0.02]
    present_layers = [l for l in LAYER_BINS if layer_fracs[l] >= 0.02]

    p5, p95 = np.percentile(depths, [5, 95])

    results.append({
        'sample': sid, 'diagnosis': dx, 'n_cells': n_cells,
        'depth_mean': depths.mean(), 'depth_std': depths.std(),
        'depth_min': depths.min(), 'depth_max': depths.max(),
        'p5': p5, 'p95': p95,
        'ks_stat': ks_stat, 'ks_p': ks_p, 'wasserstein': w_dist,
        'n_layers_present': len(present_layers),
        'missing_layers': ', '.join(missing_layers) if missing_layers else 'none',
        'present_layers': ', '.join(present_layers),
        **{f'frac_{l}': layer_fracs[l] for l in LAYER_BINS},
        **{f'count_{l}': layer_counts[l] for l in LAYER_BINS},
    })

df = pd.DataFrame(results)
df = df.sort_values('wasserstein')

print("\nPer-sample summary (sorted by Wasserstein distance to MERFISH):")
print(f"{'Sample':<10} {'Dx':<8} {'n_cells':>8} {'mean':>6} {'std':>6} {'p5':>6} {'p95':>6} {'W_dist':>7} {'Missing':>20}")
for _, r in df.iterrows():
    print(f"{r['sample']:<10} {r['diagnosis']:<8} {r['n_cells']:>8,} {r['depth_mean']:>6.3f} {r['depth_std']:>6.3f} "
          f"{r['p5']:>6.3f} {r['p95']:>6.3f} {r['wasserstein']:>7.3f} {r['missing_layers']:>20}")

# ============================================================
# 4. Flag sections with missing layers
# ============================================================
print("\n" + "="*80)
print("LAYER COVERAGE ANALYSIS")
print("="*80)

print(f"\nMERFISH layer fractions (target):")
for l, f in merfish_layer_fracs.items():
    print(f"  {l}: {f:.1%}")

print(f"\nXenium layer fractions per sample:")
print(f"{'Sample':<10} {'Dx':<8} " + " ".join(f"{l:>7}" for l in LAYER_BINS))
for _, r in df.sort_values('sample').iterrows():
    fracs = " ".join(f"{r[f'frac_{l}']:>7.1%}" for l in LAYER_BINS)
    flag = " *** OUTLIER" if r['missing_layers'] != 'none' else ""
    print(f"{r['sample']:<10} {r['diagnosis']:<8} {fracs}{flag}")

outliers = df[df['missing_layers'] != 'none']
print(f"\n*** FLAGGED OUTLIER SECTIONS (missing layers with <2% of cells): ***")
if len(outliers) > 0:
    for _, r in outliers.iterrows():
        print(f"  {r['sample']} ({r['diagnosis']}): missing {r['missing_layers']}")
else:
    print("  None — all sections have >=2% in all layers")

print(f"\nStricter flagging (<5% in any cortical layer L1-L6, WM excluded):")
for _, r in df.sort_values('sample').iterrows():
    missing_strict = [l for l in CORTICAL_LAYERS if r[f'frac_{l}'] < 0.05]
    if missing_strict:
        frac_strs = ', '.join(f"{l}={r[f'frac_{l}']:.1%}" for l in missing_strict)
        print(f"  {r['sample']} ({r['diagnosis']}): under-represented: {', '.join(missing_strict)} "
              f"(fracs: {frac_strs})")

# ============================================================
# 5. Resampling feasibility analysis
# ============================================================
print("\n" + "="*80)
print("RESAMPLING FEASIBILITY ANALYSIS")
print("="*80)

n_bins = 20
bin_edges = np.linspace(0, 1, n_bins + 1)

merfish_hist, _ = np.histogram(merfish_depth_annotated, bins=bin_edges)
merfish_props = merfish_hist / merfish_hist.sum()

print(f"\nMERFISH target distribution ({n_bins} bins from 0 to 1):")
for i in range(n_bins):
    print(f"  [{bin_edges[i]:.2f}, {bin_edges[i+1]:.2f}): {merfish_props[i]:.3f} ({merfish_hist[i]:,} cells)")

print(f"\nPer-sample resampling budget:")
print(f"{'Sample':<10} {'Dx':<8} {'n_orig':>8} {'n_resamp':>8} {'retain%':>8} {'bottleneck_bin':>15}")

resamp_results = []
for _, r in df.sort_values('sample').iterrows():
    sid = r['sample']
    mask = xenium.obs['sample_id'] == sid
    depths = xenium.obs.loc[mask, 'predicted_norm_depth'].values.astype(float)
    n_orig = len(depths)

    sample_hist, _ = np.histogram(depths, bins=bin_edges)

    max_n = np.inf
    bottleneck = -1
    for i in range(n_bins):
        if merfish_props[i] > 0.001:
            possible = sample_hist[i] / merfish_props[i]
            if possible < max_n:
                max_n = possible
                bottleneck = i

    max_n = int(max_n)
    retain_pct = max_n / n_orig * 100
    bn_label = f"[{bin_edges[bottleneck]:.2f}-{bin_edges[bottleneck+1]:.2f})"

    resamp_results.append({
        'sample': sid, 'diagnosis': r['diagnosis'],
        'n_orig': n_orig, 'n_resamp': max_n, 'retain_pct': retain_pct,
        'bottleneck_bin': bn_label, 'bottleneck_idx': bottleneck
    })
    print(f"{sid:<10} {r['diagnosis']:<8} {n_orig:>8,} {max_n:>8,} {retain_pct:>7.1f}% {bn_label:>15}")

resamp_df = pd.DataFrame(resamp_results)
print(f"\nOverall resampling summary:")
print(f"  Total original cells: {resamp_df['n_orig'].sum():,}")
print(f"  Total after resampling: {resamp_df['n_resamp'].sum():,}")
print(f"  Overall retention: {resamp_df['n_resamp'].sum() / resamp_df['n_orig'].sum() * 100:.1f}%")
print(f"  Min cells after resampling: {resamp_df['n_resamp'].min():,} ({resamp_df.loc[resamp_df['n_resamp'].idxmin(), 'sample']})")
print(f"  Max cells after resampling: {resamp_df['n_resamp'].max():,}")
print(f"  Median retention: {resamp_df['retain_pct'].median():.1f}%")

for dx in ['SCZ', 'Control']:
    sub = resamp_df[resamp_df['diagnosis'] == dx]
    print(f"\n  {dx}: {len(sub)} samples")
    print(f"    Original: {sub['n_orig'].sum():,} cells (mean {sub['n_orig'].mean():,.0f})")
    print(f"    After resampling: {sub['n_resamp'].sum():,} cells (mean {sub['n_resamp'].mean():,.0f})")
    print(f"    Retention: {sub['n_resamp'].sum() / sub['n_orig'].sum() * 100:.1f}%")

# ============================================================
# 6. Save analysis table
# ============================================================
out_csv = os.path.join(OUTPUT_DIR, 'depth_resampling_analysis.csv')
df.to_csv(out_csv, index=False)
print(f"\nSaved analysis table: {out_csv}")

# ============================================================
# 7. Generate comprehensive figure
# ============================================================
fig, axes = plt.subplots(3, 2, figsize=(28, 36))

# (0,0): MERFISH vs all Xenium depth histogram overlay
axes[0,0].hist(merfish_depth_annotated, bins=80, density=True, alpha=0.6, color='black', label='MERFISH (annotated)')
for _, r in df.iterrows():
    mask = xenium.obs['sample_id'] == r['sample']
    d = xenium.obs.loc[mask, 'predicted_norm_depth'].values.astype(float)
    color = 'red' if r['diagnosis'] == 'SCZ' else 'blue'
    axes[0,0].hist(d, bins=80, density=True, alpha=0.08, color=color, histtype='step', linewidth=1)
axes[0,0].hist([], color='red', alpha=0.5, label='Xenium SCZ')
axes[0,0].hist([], color='blue', alpha=0.5, label='Xenium Control')
axes[0,0].set_xlabel('Normalized depth from pia', fontsize=16)
axes[0,0].set_ylabel('Density', fontsize=16)
axes[0,0].set_title('Depth Distributions: MERFISH vs Xenium', fontsize=20, fontweight='bold')
axes[0,0].legend(fontsize=14)
axes[0,0].tick_params(labelsize=13)
for b in [0.1, 0.3, 0.45, 0.65, 0.85]:
    axes[0,0].axvline(b, color='gray', linestyle=':', alpha=0.5)

# (0,1): Per-sample layer fractions heatmap
layer_names = list(LAYER_BINS.keys())
samples_sorted = df.sort_values('sample')['sample'].values
frac_matrix = np.array([[df.loc[df['sample']==s, f'frac_{l}'].values[0] for l in layer_names] for s in samples_sorted])
display_labels = [f"{s} ({SAMPLE_TO_DX.get(s, '?')})" for s in samples_sorted]

im = axes[0,1].imshow(frac_matrix, aspect='auto', cmap='YlOrRd', vmin=0, vmax=0.4)
axes[0,1].set_xticks(range(len(layer_names)))
axes[0,1].set_xticklabels(layer_names, fontsize=14, fontweight='bold')
axes[0,1].set_yticks(range(len(samples_sorted)))
axes[0,1].set_yticklabels(display_labels, fontsize=11)
for i in range(len(samples_sorted)):
    for j in range(len(layer_names)):
        val = frac_matrix[i, j]
        axes[0,1].text(j, i, f'{val:.0%}', ha='center', va='center', fontsize=10,
                      fontweight='bold' if val < 0.02 else 'normal',
                      color='white' if val > 0.2 else 'black')
plt.colorbar(im, ax=axes[0,1], shrink=0.8, label='Fraction of cells')
axes[0,1].set_title('Layer Fractions per Sample\n(vs MERFISH target)', fontsize=20, fontweight='bold')

# (1,0): Wasserstein distance by diagnosis
scz_w = df[df['diagnosis']=='SCZ']['wasserstein'].values
ctrl_w = df[df['diagnosis']=='Control']['wasserstein'].values
axes[1,0].boxplot([ctrl_w, scz_w], labels=['Control', 'SCZ'])
axes[1,0].set_ylabel('Wasserstein distance to MERFISH', fontsize=16)
axes[1,0].set_title('Distributional Distance to MERFISH\nby Diagnosis', fontsize=20, fontweight='bold')
axes[1,0].tick_params(labelsize=13)
for i, (vals, color) in enumerate([(ctrl_w, 'blue'), (scz_w, 'red')]):
    axes[1,0].scatter(np.ones_like(vals) * (i+1) + np.random.normal(0, 0.05, len(vals)), vals,
                     color=color, alpha=0.6, s=60, zorder=5)

# (1,1): Resampling retention by sample
resamp_sorted = resamp_df.sort_values('retain_pct', ascending=True)
colors_bar = ['red' if d == 'SCZ' else 'blue' for d in resamp_sorted['diagnosis']]
axes[1,1].barh(range(len(resamp_sorted)), resamp_sorted['retain_pct'].values, color=colors_bar, alpha=0.7)
axes[1,1].set_yticks(range(len(resamp_sorted)))
axes[1,1].set_yticklabels([f"{r['sample']} ({r['diagnosis']})" for _, r in resamp_sorted.iterrows()], fontsize=11)
axes[1,1].set_xlabel('% cells retained after resampling', fontsize=16)
axes[1,1].set_title('Resampling Retention per Sample\n(matching MERFISH depth distribution)', fontsize=20, fontweight='bold')
axes[1,1].axvline(50, color='gray', linestyle='--', alpha=0.5)
axes[1,1].tick_params(labelsize=13)

# (2,0): Resampled cell counts by diagnosis
axes[2,0].clear()
dx_groups = {'Control': resamp_df[resamp_df['diagnosis']=='Control'].sort_values('n_resamp'),
             'SCZ': resamp_df[resamp_df['diagnosis']=='SCZ'].sort_values('n_resamp')}
all_sorted = pd.concat([dx_groups['Control'], dx_groups['SCZ']])
bar_colors = ['blue' if d == 'Control' else 'red' for d in all_sorted['diagnosis']]
axes[2,0].barh(range(len(all_sorted)), all_sorted['n_resamp'].values, color=bar_colors, alpha=0.7)
axes[2,0].set_yticks(range(len(all_sorted)))
axes[2,0].set_yticklabels([f"{r['sample']} ({r['diagnosis']})" for _, r in all_sorted.iterrows()], fontsize=11)
axes[2,0].set_xlabel('Cells after resampling', fontsize=16)
axes[2,0].set_title('Cell Count After Depth-Matched Resampling', fontsize=20, fontweight='bold')
axes[2,0].tick_params(labelsize=13)

# (2,1): Bottleneck bins
bn_counts = Counter(resamp_df['bottleneck_bin'].values)
bn_labels = sorted(bn_counts.keys())
axes[2,1].bar(range(len(bn_labels)), [bn_counts[b] for b in bn_labels], color='steelblue')
axes[2,1].set_xticks(range(len(bn_labels)))
axes[2,1].set_xticklabels(bn_labels, rotation=45, ha='right', fontsize=12)
axes[2,1].set_ylabel('Number of samples', fontsize=16)
axes[2,1].set_xlabel('Depth bin', fontsize=16)
axes[2,1].set_title('Resampling Bottleneck Bins\n(which depth bins limit sample size)', fontsize=20, fontweight='bold')
axes[2,1].tick_params(labelsize=13)

fig.suptitle('Depth-Matched Resampling Analysis\nMERFISH Target vs Xenium Sections',
             fontsize=26, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])

out_fig = os.path.join(OUTPUT_DIR, 'depth_resampling_analysis.png')
plt.savefig(out_fig, dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print(f"\nSaved figure: {out_fig}")
print("\nDone!")
