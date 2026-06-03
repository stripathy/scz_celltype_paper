"""
Depth-stratified scCODA compositional analysis: SCZ vs Control.

Runs scCODA within each cortical depth stratum (L2/3, L4, L5, L6)
at the subclass level, with Sex and Age as covariates.

Each cell type is tested as reference (following SEA-AD approach),
and results are aggregated.

Usage:
    python run_sccoda_stratified.py
"""

import os
import sys
import time
import pickle
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from datetime import datetime

# Monkey-patch scipy.signal.gaussian for scCODA compatibility
import scipy.signal
if not hasattr(scipy.signal, 'gaussian'):
    scipy.signal.gaussian = lambda M, std: np.exp(
        -0.5 * ((np.arange(0, M) - (M - 1) / 2) / std) ** 2
    )

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat

# ── Config imports ──────────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR, SAMPLE_TO_DX, METADATA_PATH, SUBCLASS_CONF_THRESH

H5AD_PATH = os.path.join(BASE_DIR, "output", "all_samples_annotated.h5ad")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "sccoda")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Depth strata ───────────────────────────────────────────────────
STRATA = {
    "L2/3": (0.10, 0.30),
    "L4":   (0.30, 0.45),
    "L5":   (0.45, 0.65),
    "L6":   (0.65, 0.85),
}


def load_metadata():
    """Load subject metadata from Excel file."""
    sys.path.insert(0, os.path.join(BASE_DIR, "code"))
    from modules.metadata import get_subject_info
    info = get_subject_info(METADATA_PATH)
    meta = {}
    for _, row in info.iterrows():
        sid = row['Brain Number']
        meta[sid] = {
            'Sex': row['Sex'],
            'Age': float(row['Age']) if row['Age'] else np.nan,
        }
    return meta


def build_count_table(obs_df, stratum_name, depth_range, metadata):
    """
    Build a cell count table for one depth stratum.

    Returns an AnnData with:
      .X = cell counts (samples x cell_types)
      .obs = covariates (diagnosis, Sex, Age)
      .var = cell type names
    """
    lo, hi = depth_range
    mask = (obs_df['predicted_norm_depth'] >= lo) & (obs_df['predicted_norm_depth'] < hi)
    stratum_df = obs_df.loc[mask]

    print(f"  {stratum_name}: {mask.sum():,} cells in depth [{lo}, {hi})")

    # Cross-tabulate: samples x subclasses
    ct = pd.crosstab(stratum_df['sample_id'], stratum_df['subclass_label'])

    # Remove cell types that are absent across all samples
    ct = ct.loc[:, ct.sum() > 0]

    # Remove cell types present in fewer than 50% of samples (scCODA struggles with very sparse types)
    presence = (ct > 0).mean()
    keep = presence[presence >= 0.5].index
    dropped = set(ct.columns) - set(keep)
    if dropped:
        print(f"    Dropped rare types (<50% presence): {dropped}")
    ct = ct[keep]

    print(f"    {ct.shape[0]} samples x {ct.shape[1]} cell types")
    print(f"    Cells per sample: {ct.sum(axis=1).min()}-{ct.sum(axis=1).max()}")

    # Build covariate DataFrame
    cov = pd.DataFrame(index=ct.index)
    cov['diagnosis'] = [SAMPLE_TO_DX[s] for s in ct.index]
    for sid in ct.index:
        if sid in metadata:
            cov.loc[sid, 'Sex'] = metadata[sid]['Sex']
            cov.loc[sid, 'Age'] = metadata[sid]['Age']
        else:
            cov.loc[sid, 'Sex'] = 'Unknown'
            cov.loc[sid, 'Age'] = np.nan

    # Convert Age to float and center it
    cov['Age'] = cov['Age'].astype(float)
    cov['Age_centered'] = cov['Age'] - cov['Age'].mean()

    # Make sure diagnosis is properly ordered (Control first as reference)
    # This is handled by patsy formula with Treatment encoding

    # Create AnnData for scCODA
    count_adata = ad.AnnData(
        X=ct.values.astype(np.float64),
        obs=cov,
        var=pd.DataFrame(index=ct.columns),
    )
    count_adata.obs_names = ct.index.tolist()

    return count_adata


def run_sccoda_all_references(count_adata, stratum_name, formula,
                               num_results=20000, num_burnin=5000,
                               min_acceptance=0.4, max_retries=3):
    """
    Run scCODA with each cell type as reference (following SEA-AD).

    Returns a DataFrame with results from all reference cell types.
    """
    cell_types = count_adata.var_names.tolist()
    all_results = []

    for ct in cell_types:
        print(f"    Reference: {ct}")

        for attempt in range(max_retries):
            try:
                model = mod.CompositionalAnalysis(
                    count_adata,
                    formula=formula,
                    reference_cell_type=ct,
                )
                result = model.sample_hmc(
                    num_results=num_results,
                    num_burnin=num_burnin,
                    verbose=False,
                )

                # Check acceptance rate
                accepted = result.sample_stats["is_accepted"].to_numpy()
                acc_rate = accepted.sum() / accepted.shape[1]

                if acc_rate < min_acceptance:
                    print(f"      Attempt {attempt+1}: acceptance {acc_rate:.2f} < {min_acceptance}, retrying...")
                    continue

                print(f"      Converged! acceptance={acc_rate:.2f}")

                # Extract results
                summary_df = result.summary_prepare()[1]
                summary_df = summary_df.reset_index()
                summary_df['reference_cell_type'] = ct
                summary_df['stratum'] = stratum_name
                summary_df['acceptance_rate'] = acc_rate
                all_results.append(summary_df)

                # Save individual result
                pkl_path = os.path.join(
                    OUTPUT_DIR, f"{stratum_name}_{ct.replace('/', '_')}_result.pkl"
                )
                with open(pkl_path, 'wb') as f:
                    pickle.dump(result, f)

                break

            except Exception as e:
                print(f"      Attempt {attempt+1} FAILED: {e}")
                if attempt == max_retries - 1:
                    print(f"      Skipping {ct} after {max_retries} failures")

    if all_results:
        return pd.concat(all_results, ignore_index=True)
    return pd.DataFrame()


def aggregate_results(results_df):
    """
    Aggregate scCODA results across all reference cell types.

    For each cell type and covariate, compute:
    - Fraction of reference choices where the effect was credible (inclusion prob > 0.9)
    - Mean log-fold change across all references
    - Mean inclusion probability
    """
    if results_df.empty:
        return pd.DataFrame()

    # Filter to the diagnosis effect (the one we care about most)
    # The covariate name should contain 'diagnosis' or 'SCZ'
    diag_mask = results_df['Covariate'].str.contains('diagnosis|SCZ', case=False, na=False)
    diag_results = results_df[diag_mask].copy()

    if diag_results.empty:
        print("  Warning: No diagnosis effects found in results")
        return results_df

    # For each cell type, aggregate across reference choices
    agg = diag_results.groupby(['stratum', 'Cell Type']).agg(
        mean_log2_fc=('log2-fold change', 'mean'),
        median_log2_fc=('log2-fold change', 'median'),
        mean_inclusion_prob=('Inclusion probability', 'mean'),
        frac_credible=('Final Parameter', lambda x: (x != 0).mean()),
        n_references=('reference_cell_type', 'nunique'),
        mean_expected_fc=('Expected fold change', 'mean'),
    ).reset_index()

    agg = agg.sort_values(['stratum', 'frac_credible'], ascending=[True, False])

    return agg


def plot_results(agg_df, all_results_df):
    """Generate summary figures."""

    if agg_df.empty:
        print("  No aggregated results to plot")
        return

    strata_names = sorted(agg_df['stratum'].unique())

    # ── Figure 1: Heatmap of inclusion probabilities ──
    fig, axes = plt.subplots(1, len(strata_names), figsize=(6*len(strata_names), 10),
                              squeeze=False)

    for idx, stratum in enumerate(strata_names):
        ax = axes[0, idx]
        sub = agg_df[agg_df['stratum'] == stratum].sort_values('mean_log2_fc')

        colors = ['#d73027' if fc > 0 else '#4575b4' for fc in sub['mean_log2_fc']]
        alphas = [max(0.3, min(1.0, ip)) for ip in sub['mean_inclusion_prob']]

        bars = ax.barh(range(len(sub)), sub['mean_log2_fc'], color=colors)
        for bar, alpha in zip(bars, alphas):
            bar.set_alpha(alpha)

        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(sub['Cell Type'], fontsize=12)
        ax.set_xlabel('Mean log2 FC (SCZ vs Control)', fontsize=14)
        ax.set_title(f'{stratum}\n(bar opacity = inclusion prob)', fontsize=16)
        ax.axvline(0, color='black', linewidth=0.5)

        # Mark credible effects
        for i, (_, row) in enumerate(sub.iterrows()):
            if row['frac_credible'] > 0.5:
                ax.text(row['mean_log2_fc'], i, ' *', fontsize=16, fontweight='bold',
                       va='center', color='black')

    plt.suptitle('scCODA: SCZ vs Control by Cortical Layer\n(* = credible in >50% of reference choices)',
                 fontsize=18, y=1.02)
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "sccoda_stratified_barplot.png"),
                dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved sccoda_stratified_barplot.png")

    # ── Figure 2: Dot plot (inclusion prob vs fold change) ──
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    for idx, stratum in enumerate(strata_names):
        ax = axes[idx // 2, idx % 2]
        sub = agg_df[agg_df['stratum'] == stratum]

        scatter = ax.scatter(
            sub['mean_log2_fc'], sub['mean_inclusion_prob'],
            s=100, c=['#d73027' if fc > 0 else '#4575b4' for fc in sub['mean_log2_fc']],
            edgecolors='black', linewidth=0.5, alpha=0.8
        )

        # Label points with high inclusion probability
        for _, row in sub.iterrows():
            if row['mean_inclusion_prob'] > 0.5 or abs(row['mean_log2_fc']) > 0.5:
                ax.annotate(row['Cell Type'], (row['mean_log2_fc'], row['mean_inclusion_prob']),
                           fontsize=9, ha='center', va='bottom')

        ax.axhline(0.9, color='red', linestyle='--', alpha=0.5, label='Inclusion=0.9')
        ax.axvline(0, color='black', linewidth=0.5)
        ax.set_xlabel('Mean log2 FC (SCZ vs Control)', fontsize=14)
        ax.set_ylabel('Mean Inclusion Probability', fontsize=14)
        ax.set_title(f'{stratum}', fontsize=16)
        ax.set_ylim(-0.05, 1.05)
        ax.legend(fontsize=10)

    plt.suptitle('scCODA: Effect Size vs Inclusion Probability\n(Red=increased in SCZ, Blue=decreased)',
                 fontsize=18)
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "sccoda_stratified_dotplot.png"),
                dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved sccoda_stratified_dotplot.png")

    # ── Figure 3: Summary table as figure ──
    credible = agg_df[agg_df['frac_credible'] > 0.5].copy()
    if not credible.empty:
        credible = credible.sort_values(['stratum', 'frac_credible'], ascending=[True, False])

        fig, ax = plt.subplots(figsize=(14, max(4, len(credible)*0.5 + 2)))
        ax.axis('off')

        table_data = credible[['stratum', 'Cell Type', 'mean_log2_fc',
                                'mean_inclusion_prob', 'frac_credible']].copy()
        table_data.columns = ['Layer', 'Cell Type', 'Mean log2 FC',
                              'Mean Incl. Prob', 'Frac Credible']
        table_data['Mean log2 FC'] = table_data['Mean log2 FC'].map('{:.3f}'.format)
        table_data['Mean Incl. Prob'] = table_data['Mean Incl. Prob'].map('{:.3f}'.format)
        table_data['Frac Credible'] = table_data['Frac Credible'].map('{:.2f}'.format)

        table = ax.table(
            cellText=table_data.values,
            colLabels=table_data.columns,
            loc='center',
            cellLoc='center',
        )
        table.auto_set_font_size(False)
        table.set_fontsize(12)
        table.scale(1.2, 1.8)

        ax.set_title('Credible Effects (>50% of reference choices)\nSCZ vs Control, depth-stratified',
                     fontsize=16, pad=20)

        fig.savefig(os.path.join(OUTPUT_DIR, "sccoda_credible_effects_table.png"),
                    dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved sccoda_credible_effects_table.png")


def main():
    t_start = time.time()

    # ── 1. Load obs data (don't need expression matrix) ──
    print("Loading annotated data (obs only)...")
    adata = ad.read_h5ad(H5AD_PATH, backed='r')
    obs_df = adata.obs.copy()
    obs_df.index = range(len(obs_df))  # Reset index to avoid duplicates
    del adata
    print(f"  {len(obs_df):,} cells, {obs_df['sample_id'].nunique()} samples")

    # Bottom-1% subclass confidence filter
    if 'subclass_label_confidence' in obs_df.columns:
        conf = obs_df['subclass_label_confidence'].astype(float)
        conf_mask = conf >= SUBCLASS_CONF_THRESH
        n_drop = (~conf_mask).sum()
        obs_df = obs_df[conf_mask]
        obs_df.index = range(len(obs_df))
        print(f"  Confidence filter: dropped {n_drop:,} cells, {len(obs_df):,} remaining")

    # ── 2. Load metadata ──
    print("\nLoading metadata...")
    metadata = load_metadata()
    print(f"  Metadata for {len(metadata)} subjects")
    for sid in sorted(SAMPLE_TO_DX.keys()):
        if sid in metadata:
            m = metadata[sid]
            print(f"    {sid}: {SAMPLE_TO_DX[sid]}, {m['Sex']}, Age={m['Age']}")

    # ── 3. Build count tables and run scCODA per stratum ──
    all_raw_results = []
    all_agg_results = []

    formula = "C(diagnosis, Treatment('Control')) + C(Sex, Treatment('F')) + Age_centered"
    print(f"\nFormula: {formula}")

    for stratum_name, depth_range in STRATA.items():
        print(f"\n{'='*60}")
        print(f"STRATUM: {stratum_name}")
        print(f"{'='*60}")

        # Build count table
        count_adata = build_count_table(obs_df, stratum_name, depth_range, metadata)

        # Print diagnosis balance
        dx_counts = count_adata.obs['diagnosis'].value_counts()
        print(f"    Diagnosis: {dict(dx_counts)}")

        # Run scCODA
        print(f"\n  Running scCODA (iterating over {count_adata.shape[1]} reference cell types)...")
        t_stratum = time.time()

        raw_results = run_sccoda_all_references(
            count_adata, stratum_name, formula,
            num_results=20000, num_burnin=5000,
            min_acceptance=0.4,
        )

        elapsed = time.time() - t_stratum
        print(f"\n  Stratum {stratum_name} done in {elapsed:.0f}s")

        if not raw_results.empty:
            all_raw_results.append(raw_results)

            # Aggregate
            agg = aggregate_results(raw_results)
            if not agg.empty:
                all_agg_results.append(agg)

                # Print top effects
                print(f"\n  Top effects for {stratum_name}:")
                top = agg.head(5)
                for _, row in top.iterrows():
                    print(f"    {row['Cell Type']}: log2FC={row['mean_log2_fc']:.3f}, "
                          f"incl_prob={row['mean_inclusion_prob']:.3f}, "
                          f"frac_credible={row['frac_credible']:.2f}")

    # ── 4. Combine and save results ──
    print(f"\n{'='*60}")
    print("SAVING RESULTS")
    print(f"{'='*60}")

    if all_raw_results:
        raw_df = pd.concat(all_raw_results, ignore_index=True)
        raw_path = os.path.join(OUTPUT_DIR, "sccoda_raw_results.csv")
        raw_df.to_csv(raw_path, index=False)
        print(f"  Raw results: {raw_path} ({len(raw_df)} rows)")

    if all_agg_results:
        agg_df = pd.concat(all_agg_results, ignore_index=True)
        agg_path = os.path.join(OUTPUT_DIR, "sccoda_aggregated_results.csv")
        agg_df.to_csv(agg_path, index=False)
        print(f"  Aggregated results: {agg_path} ({len(agg_df)} rows)")

        # ── 5. Generate figures ──
        print("\nGenerating figures...")
        plot_results(agg_df, raw_df)

    total = time.time() - t_start
    print(f"\nTotal time: {total:.0f}s")
    print(f"Output: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
