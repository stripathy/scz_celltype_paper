#!/usr/bin/env python3
"""
Calibrate correlation classifier margin threshold using MERFISH ground truth.

Uses 1.89M SEA-AD MERFISH cells with Allen Institute ground-truth labels
to empirically measure how corr_subclass_margin relates to misclassification
rate, then identifies principled absolute thresholds.

Steps:
  1. Load MERFISH reclassified data (has both ground-truth and pipeline margin)
  2. Build margin-vs-accuracy curves (overall and by class)
  3. Identify candidate thresholds at specific accuracy levels
  4. Cross-validate stability across sections
  5. Apply candidates to Xenium and compare to current 1% percentile

Output: output/presentation/margin_calibration_*.png
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
sys.path.insert(0, os.path.join(BASE_DIR, "code", "analysis"))
from config import SUBCLASS_TO_CLASS, PRESENTATION_DIR

MERFISH_RECLASSIFIED = os.path.join(BASE_DIR, "output", "merfish_benchmark",
                                     "merfish_reclassified.h5ad")
H5AD_DIR = os.path.join(BASE_DIR, "output", "h5ad")
OUT_DIR = PRESENTATION_DIR


def load_merfish_data():
    """Load MERFISH obs with ground truth and pipeline margin."""
    print("Loading MERFISH reclassified data...")
    adata = ad.read_h5ad(MERFISH_RECLASSIFIED, backed='r')
    cols = ['Subclass', 'Supertype', 'Class', 'Section', 'Donor ID',
            'corr_subclass', 'corr_subclass_corr', 'corr_subclass_margin',
            'corr_supertype', 'corr_qc_pass', 'doublet_suspect']
    df = adata.obs[cols].copy()
    adata.file.close()
    print(f"  {len(df):,} cells loaded")

    # Compute correctness
    df['correct_subclass'] = df['corr_subclass'] == df['Subclass']
    df['correct_supertype'] = df['corr_supertype'] == df['Supertype']

    # Infer class from corr_subclass
    df['corr_class_simple'] = df['corr_subclass'].map(
        lambda x: SUBCLASS_TO_CLASS.get(x, 'Other'))

    # Ground truth class
    class_map = {
        'Neuronal: Glutamatergic': 'Glut',
        'Neuronal: GABAergic': 'GABA',
        'Non-neuronal and Non-neural': 'NN',
    }
    df['gt_class'] = df['Class'].astype(str).map(class_map).fillna('Other')

    return df


def margin_accuracy_curve(df, n_bins=50):
    """Compute accuracy as a function of margin, using equal-count bins."""
    margins = df['corr_subclass_margin'].values
    correct = df['correct_subclass'].values

    # Equal-count bins for smooth curve
    sorted_idx = np.argsort(margins)
    bin_size = len(sorted_idx) // n_bins
    bins = []
    for i in range(n_bins):
        start = i * bin_size
        end = (i + 1) * bin_size if i < n_bins - 1 else len(sorted_idx)
        idx = sorted_idx[start:end]
        bins.append({
            'margin_low': margins[idx].min(),
            'margin_high': margins[idx].max(),
            'margin_median': np.median(margins[idx]),
            'accuracy': correct[idx].mean(),
            'n_cells': len(idx),
        })
    return pd.DataFrame(bins)


def cumulative_accuracy_below(df, thresholds):
    """For each threshold, compute accuracy of cells BELOW that margin."""
    margins = df['corr_subclass_margin'].values
    correct = df['correct_subclass'].values
    results = []
    for t in thresholds:
        mask = margins < t
        n = mask.sum()
        if n > 0:
            acc = correct[mask].mean()
        else:
            acc = np.nan
        results.append({
            'threshold': t,
            'n_below': n,
            'pct_below': 100 * n / len(margins),
            'accuracy_below': acc,
            'n_above': (~mask).sum(),
            'accuracy_above': correct[~mask].mean() if (~mask).sum() > 0 else np.nan,
        })
    return pd.DataFrame(results)


def find_threshold_for_accuracy(df, target_accuracy):
    """Find the margin threshold below which accuracy drops to target."""
    margins = df['corr_subclass_margin'].values
    correct = df['correct_subclass'].values

    # Sort by margin
    sorted_idx = np.argsort(margins)
    sorted_margins = margins[sorted_idx]
    sorted_correct = correct[sorted_idx]

    # Sliding window: cumulative accuracy from lowest margin upward
    cum_correct = np.cumsum(sorted_correct)
    cum_count = np.arange(1, len(sorted_correct) + 1)
    cum_accuracy = cum_correct / cum_count

    # Find where cumulative accuracy (of cells below threshold) equals target
    # We want: accuracy of cells with margin < threshold = target_accuracy
    # Search from high margin downward
    for i in range(len(cum_accuracy) - 1, 0, -1):
        if cum_accuracy[i] <= target_accuracy:
            return sorted_margins[i]
    return sorted_margins[0]


# ══════════════════════════════════════════════════════════════════════
# STEP 1 & 2: MERFISH calibration curves
# ══════════════════════════════════════════════════════════════════════

def step1_calibration_curves(df):
    """Build and plot margin-vs-accuracy curves."""
    print("\n=== Step 1: Margin calibration curves ===")

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # Panel A: Overall accuracy vs margin (binned)
    ax = axes[0, 0]
    curve = margin_accuracy_curve(df, n_bins=80)
    ax.scatter(curve['margin_median'], curve['accuracy'], s=20, c='#4477AA', alpha=0.8)
    ax.set_xlabel('Correlation margin (best - 2nd best)', fontsize=13)
    ax.set_ylabel('Subclass accuracy (vs Allen ground truth)', fontsize=13)
    ax.set_title('A) Accuracy vs Margin (MERFISH, equal-count bins)', fontsize=14, fontweight='bold')
    ax.axhline(0.8, color='red', linestyle='--', alpha=0.5, label='80% accuracy')
    ax.axhline(0.7, color='orange', linestyle='--', alpha=0.5, label='70% accuracy')
    ax.set_ylim(0.3, 1.02)
    ax.legend(fontsize=11)
    ax.grid(alpha=0.2)

    # Panel B: Accuracy vs margin by cell class
    ax = axes[0, 1]
    class_colors = {'Glut': '#00ADF8', 'GABA': '#F05A28', 'NN': '#808080'}
    for cls, color in class_colors.items():
        sub = df[df['gt_class'] == cls]
        curve_cls = margin_accuracy_curve(sub, n_bins=40)
        ax.scatter(curve_cls['margin_median'], curve_cls['accuracy'],
                   s=20, c=color, alpha=0.7, label=f'{cls} (n={len(sub):,})')
    ax.set_xlabel('Correlation margin', fontsize=13)
    ax.set_ylabel('Subclass accuracy', fontsize=13)
    ax.set_title('B) Accuracy vs Margin by Cell Class', fontsize=14, fontweight='bold')
    ax.axhline(0.8, color='red', linestyle='--', alpha=0.3)
    ax.set_ylim(0.3, 1.02)
    ax.legend(fontsize=11)
    ax.grid(alpha=0.2)

    # Panel C: Accuracy-vs-retention tradeoff
    ax = axes[1, 0]
    thresholds = np.linspace(0, 0.3, 200)
    tradeoff = cumulative_accuracy_below(df, thresholds)

    # Plot: for each threshold, what's the accuracy of cells ABOVE it (kept cells)
    ax.plot(100 - tradeoff['pct_below'], tradeoff['accuracy_above'],
            color='#4477AA', linewidth=2)
    ax.set_xlabel('% cells retained', fontsize=13)
    ax.set_ylabel('Accuracy of retained cells', fontsize=13)
    ax.set_title('C) Accuracy-Retention Tradeoff', fontsize=14, fontweight='bold')
    ax.axhline(0.85, color='grey', linestyle=':', alpha=0.5)
    ax.set_xlim(85, 100.5)
    ax.set_ylim(0.8, 1.0)
    ax.grid(alpha=0.2)

    # Mark current 1% percentile equivalent
    pctl_1 = np.percentile(df['corr_subclass_margin'], 1)
    pctl_retention = 99.0
    pctl_acc = df.loc[df['corr_subclass_margin'] >= pctl_1, 'correct_subclass'].mean()
    ax.axvline(pctl_retention, color='red', linestyle='--', alpha=0.5,
               label=f'Current 1% pctl (margin={pctl_1:.3f})')
    ax.plot(pctl_retention, pctl_acc, 'ro', markersize=8)
    ax.legend(fontsize=11)

    # Panel D: Margin distribution with accuracy overlay
    ax = axes[1, 1]
    margins = df['corr_subclass_margin'].values
    bins = np.linspace(0, np.percentile(margins, 99.5), 100)
    ax.hist(margins, bins=bins, color='#4477AA', alpha=0.5, density=True, label='All cells')

    # Overlay: misclassified cells
    wrong = df[~df['correct_subclass']]
    ax.hist(wrong['corr_subclass_margin'].values, bins=bins, color='#CC6677',
            alpha=0.5, density=True, label='Misclassified cells')

    ax.set_xlabel('Correlation margin', fontsize=13)
    ax.set_ylabel('Density', fontsize=13)
    ax.set_title('D) Margin Distribution: All vs Misclassified', fontsize=14, fontweight='bold')
    ax.axvline(pctl_1, color='red', linestyle='--', linewidth=1.5,
               label=f'1st pctl = {pctl_1:.3f}')
    ax.legend(fontsize=11)
    ax.grid(alpha=0.2)

    plt.tight_layout()
    outpath = os.path.join(OUT_DIR, 'margin_calibration_curve.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")

    # Print key statistics
    print(f"\n  Overall subclass accuracy: {df['correct_subclass'].mean():.4f}")
    print(f"  1st percentile margin: {pctl_1:.4f}")
    print(f"  Accuracy of cells above 1st pctl: {pctl_acc:.4f}")

    # Find thresholds at various accuracy targets
    print(f"\n  Candidate thresholds (accuracy of cells BELOW threshold):")
    for target in [0.50, 0.60, 0.70, 0.80]:
        t = find_threshold_for_accuracy(df, target)
        n_below = (margins < t).sum()
        pct_below = 100 * n_below / len(margins)
        acc_above = df.loc[df['corr_subclass_margin'] >= t, 'correct_subclass'].mean()
        print(f"    Below-threshold acc ≤ {target:.0%}: margin = {t:.4f}, "
              f"excludes {n_below:,} ({pct_below:.1f}%), "
              f"retained accuracy = {acc_above:.4f}")

    return pctl_1


# ══════════════════════════════════════════════════════════════════════
# STEP 3: Cross-validate across sections
# ══════════════════════════════════════════════════════════════════════

def step3_cross_validate(df):
    """Check threshold stability across MERFISH sections."""
    print("\n=== Step 3: Cross-section stability ===")

    sections = df['Section'].unique()
    section_stats = []
    for sec in sorted(sections):
        sub = df[df['Section'] == sec]
        margins = sub['corr_subclass_margin'].values
        correct = sub['correct_subclass'].values

        # Per-section 1st percentile
        pctl_1 = np.percentile(margins, 1)

        # Find margin where accuracy of cells below drops to 70%
        sorted_idx = np.argsort(margins)
        sorted_m = margins[sorted_idx]
        sorted_c = correct[sorted_idx]
        cum_acc = np.cumsum(sorted_c) / np.arange(1, len(sorted_c) + 1)

        # Threshold where cumulative accuracy = 0.70
        t70 = sorted_m[0]
        for i in range(len(cum_acc) - 1, 0, -1):
            if cum_acc[i] <= 0.70:
                t70 = sorted_m[i]
                break

        section_stats.append({
            'section': sec,
            'donor': sub['Donor ID'].iloc[0],
            'n_cells': len(sub),
            'accuracy': correct.mean(),
            'pctl_1_margin': pctl_1,
            'median_margin': np.median(margins),
            'threshold_70pct': t70,
        })

    sdf = pd.DataFrame(section_stats)

    print(f"  Sections: {len(sdf)}")
    print(f"  1st pctl margin: mean={sdf['pctl_1_margin'].mean():.4f}, "
          f"std={sdf['pctl_1_margin'].std():.4f}, "
          f"range=[{sdf['pctl_1_margin'].min():.4f}, {sdf['pctl_1_margin'].max():.4f}]")
    print(f"  70%-acc threshold: mean={sdf['threshold_70pct'].mean():.4f}, "
          f"std={sdf['threshold_70pct'].std():.4f}, "
          f"range=[{sdf['threshold_70pct'].min():.4f}, {sdf['threshold_70pct'].max():.4f}]")
    print(f"  Median margin: mean={sdf['median_margin'].mean():.4f}, "
          f"std={sdf['median_margin'].std():.4f}")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    ax = axes[0]
    ax.hist(sdf['pctl_1_margin'], bins=20, color='#CC6677', alpha=0.7,
            label=f'1st pctl (CV={sdf["pctl_1_margin"].std()/sdf["pctl_1_margin"].mean():.2f})')
    ax.set_xlabel('Margin threshold', fontsize=13)
    ax.set_ylabel('Number of sections', fontsize=13)
    ax.set_title('A) Per-Section 1st Percentile\n(current approach)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)

    ax = axes[1]
    ax.hist(sdf['threshold_70pct'], bins=20, color='#4477AA', alpha=0.7,
            label=f'70%-acc threshold (CV={sdf["threshold_70pct"].std()/sdf["threshold_70pct"].mean():.2f})')
    ax.set_xlabel('Margin threshold', fontsize=13)
    ax.set_ylabel('Number of sections', fontsize=13)
    ax.set_title('B) Per-Section 70%-Accuracy Threshold\n(principled approach)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)

    ax = axes[2]
    ax.scatter(sdf['pctl_1_margin'], sdf['threshold_70pct'], s=30, alpha=0.7, c='#4477AA')
    ax.set_xlabel('1st percentile margin', fontsize=13)
    ax.set_ylabel('70%-accuracy threshold', fontsize=13)
    ax.set_title('C) Percentile vs Accuracy-Based\nThreshold Per Section', fontsize=14, fontweight='bold')
    # Add diagonal
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, 'k--', alpha=0.3)
    ax.grid(alpha=0.2)

    plt.tight_layout()
    outpath = os.path.join(OUT_DIR, 'margin_threshold_cross_validation.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")

    return sdf


# ══════════════════════════════════════════════════════════════════════
# STEP 4: Apply to Xenium data
# ══════════════════════════════════════════════════════════════════════

def step4_apply_to_xenium(candidate_thresholds):
    """Apply candidate absolute thresholds to Xenium and compare to 1% percentile."""
    print("\n=== Step 4: Apply to Xenium ===")

    import glob
    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))

    rows = []
    for fpath in h5ad_files:
        sid = os.path.basename(fpath).replace("_annotated.h5ad", "")
        adata = ad.read_h5ad(fpath, backed='r')
        obs = adata.obs

        # Only look at qc_pass cells (Step 01 passers)
        qc_mask = obs['qc_pass'].astype(bool).values
        margins = obs['corr_subclass_margin'].values
        current_pass = obs['corr_qc_pass'].astype(bool).values

        # Current 1st percentile (among qc_pass cells)
        qc_margins = margins[qc_mask]
        current_pctl = np.percentile(qc_margins, 1)

        row = {
            'sample': sid,
            'n_qc_pass': qc_mask.sum(),
            'current_pctl_threshold': current_pctl,
            'current_n_excluded': qc_mask.sum() - (qc_mask & current_pass & ~obs['doublet_suspect'].astype(bool).values).sum(),
        }

        for name, thresh in candidate_thresholds.items():
            pass_mask = qc_mask & (margins >= thresh) & (~obs['doublet_suspect'].astype(bool).values)
            row[f'n_pass_{name}'] = pass_mask.sum()
            row[f'n_excluded_{name}'] = qc_mask.sum() - pass_mask.sum()
            row[f'pct_excluded_{name}'] = 100 * (qc_mask.sum() - pass_mask.sum()) / qc_mask.sum()

        rows.append(row)
        adata.file.close()
        print(f"  {sid}: pctl={current_pctl:.4f}, n_qc_pass={row['n_qc_pass']:,}")

    xdf = pd.DataFrame(rows)

    # Summary
    print(f"\n  Xenium summary across {len(xdf)} samples:")
    print(f"  Current 1st pctl range: [{xdf['current_pctl_threshold'].min():.4f}, "
          f"{xdf['current_pctl_threshold'].max():.4f}]")
    for name in candidate_thresholds:
        total_excluded = xdf[f'n_excluded_{name}'].sum()
        total_qc = xdf['n_qc_pass'].sum()
        print(f"  Threshold '{name}' ({candidate_thresholds[name]:.4f}): "
              f"excludes {total_excluded:,} / {total_qc:,} "
              f"({100*total_excluded/total_qc:.2f}%)")

    return xdf


# ══════════════════════════════════════════════════════════════════════
# STEP 5: Crumblr sensitivity
# ══════════════════════════════════════════════════════════════════════

def step5_build_crumblr_inputs(candidate_thresholds):
    """Build crumblr inputs for each candidate threshold."""
    print("\n=== Step 5: Building crumblr inputs ===")

    import glob
    sys.path.insert(0, os.path.join(BASE_DIR, "code"))
    from modules.metadata import get_subject_info

    METADATA_PATH = os.path.join(BASE_DIR, "data", "sample_metadata.xlsx")
    CRUMBLR_DIR = os.path.join(BASE_DIR, "output", "crumblr")
    EXCLUDE_SAMPLES = set()  # No samples excluded

    meta = get_subject_info(METADATA_PATH).set_index("sample_id")

    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))

    for thresh_name, thresh_val in candidate_thresholds.items():
        print(f"\n  Building input for threshold '{thresh_name}' ({thresh_val:.4f})...")
        all_obs = []

        for fpath in h5ad_files:
            sid = os.path.basename(fpath).replace("_annotated.h5ad", "")
            if sid in EXCLUDE_SAMPLES:
                continue

            adata = ad.read_h5ad(fpath, backed='r')
            obs = adata.obs

            # Filter: qc_pass AND margin >= threshold AND not doublet AND cortical
            mask = (
                obs['qc_pass'].astype(bool) &
                (obs['corr_subclass_margin'] >= thresh_val) &
                (~obs['doublet_suspect'].astype(bool)) &
                (obs['spatial_domain'] == 'Cortical')
            )
            sub = obs.loc[mask, ['sample_id', 'subclass_label', 'supertype_label']].copy()
            all_obs.append(sub)
            adata.file.close()

        df = pd.concat(all_obs, ignore_index=True)
        print(f"    {len(df):,} cortical cells from {df['sample_id'].nunique()} samples")

        # Build count tables
        for level_col, level_name in [("subclass_label", "subclass"),
                                       ("supertype_label", "supertype")]:
            counts = df.groupby(["sample_id", level_col]).size().reset_index(name="count")
            totals = df.groupby("sample_id").size().reset_index(name="total")
            counts = counts.merge(totals, on="sample_id")
            counts = counts.rename(columns={"sample_id": "donor", level_col: "celltype"})
            counts = counts.merge(
                meta[["diagnosis", "sex", "age"]].reset_index(),
                left_on="donor", right_on="sample_id", how="left"
            ).drop(columns=["sample_id"])
            counts = counts.sort_values(["donor", "celltype"]).reset_index(drop=True)

            outpath = os.path.join(CRUMBLR_DIR,
                                   f"crumblr_input_{level_name}_margin_{thresh_name}.csv")
            counts.to_csv(outpath, index=False)
            print(f"    Saved: {os.path.basename(outpath)}")


# ══════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("MARGIN THRESHOLD CALIBRATION")
    print("=" * 60)

    # Step 1: Load MERFISH and build calibration curves
    df = load_merfish_data()
    pctl_1 = step1_calibration_curves(df)

    # Step 3: Cross-validate
    sdf = step3_cross_validate(df)

    # Define candidate thresholds based on MERFISH findings
    # We'll compute these from the data
    margins = df['corr_subclass_margin'].values
    correct = df['correct_subclass'].values

    # Find thresholds at various accuracy levels
    # For cells BELOW threshold, what's their accuracy?
    t50 = find_threshold_for_accuracy(df, 0.50)
    t60 = find_threshold_for_accuracy(df, 0.60)
    t70 = find_threshold_for_accuracy(df, 0.70)
    t80 = find_threshold_for_accuracy(df, 0.80)

    candidate_thresholds = {
        'permissive': t50,
        'moderate': t70,
        'strict': t80,
    }

    print(f"\n  Candidate thresholds:")
    for name, val in candidate_thresholds.items():
        n_below = (margins < val).sum()
        print(f"    {name}: margin = {val:.4f} "
              f"(would exclude {n_below:,} = {100*n_below/len(margins):.1f}% of MERFISH)")

    # Step 4: Apply to Xenium
    xdf = step4_apply_to_xenium(candidate_thresholds)

    # Step 5: Build crumblr inputs
    step5_build_crumblr_inputs(candidate_thresholds)

    # Save candidate thresholds for reference
    thresh_df = pd.DataFrame([
        {'name': name, 'threshold': val,
         'merfish_pct_excluded': 100 * (margins < val).sum() / len(margins),
         'merfish_acc_below': correct[margins < val].mean() if (margins < val).sum() > 0 else np.nan,
         'merfish_acc_above': correct[margins >= val].mean(),
         }
        for name, val in candidate_thresholds.items()
    ])
    thresh_df.to_csv(os.path.join(OUT_DIR, 'margin_candidate_thresholds.csv'), index=False)
    print(f"\n  Saved: margin_candidate_thresholds.csv")

    print("\nDone! Now run crumblr on each threshold variant:")
    for name in candidate_thresholds:
        print(f"  Rscript code/analysis/run_crumblr.R --suffix _margin_{name}")

    print("\nDone!")


if __name__ == '__main__':
    main()
