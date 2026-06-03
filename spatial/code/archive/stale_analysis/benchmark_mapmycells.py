#!/usr/bin/env python3
"""
Benchmark MapMyCells on snRNA-seq reference restricted to Xenium panel genes.

Takes a stratified subsample from the snRNA-seq reference (which has ground truth
labels), restricts to the 293 Xenium genes (converted to Ensembl IDs), and runs
through MapMyCells. This tests: 'how well does MapMyCells assign cell types when
limited to only the Xenium gene panel?'

Output:
  - output/benchmark/benchmark_results.csv (per-cell ground truth vs prediction)
  - Accuracy and bootstrapping_probability analysis at each taxonomy level
  - Threshold analysis: what accuracy do you get at each confidence cutoff?
"""

import os
import sys
import json
import tempfile
import numpy as np
import pandas as pd
import anndata as ad

# Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
STATS_PATH = os.path.join(BASE_DIR, "data", "reference",
                          "precomputed_stats.20231120.sea_ad.MTG.h5")
REF_PATH = os.path.join(BASE_DIR, "data", "reference",
                        "Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad")
GENE_MAP_PATH = os.path.join(BASE_DIR, "data", "reference",
                             "gene_symbol_to_ensembl.json")
OUT_DIR = os.path.join(BASE_DIR, "output", "benchmark")
os.makedirs(OUT_DIR, exist_ok=True)


def create_benchmark_dataset():
    """Subsample reference, restrict to Xenium genes, convert to Ensembl IDs."""
    print("Loading snRNA-seq reference...")
    ref = ad.read_h5ad(REF_PATH)
    print(f"  Shape: {ref.shape}")

    with open(GENE_MAP_PATH) as f:
        gene_mapping = json.load(f)

    # Subset to Xenium panel genes
    xenium_genes_in_ref = [g for g in gene_mapping.keys() if g in ref.var_names]
    ref_xenium = ref[:, xenium_genes_in_ref].copy()

    # Convert to Ensembl IDs
    ref_xenium.var_names = [gene_mapping[g] for g in xenium_genes_in_ref]
    ref_xenium.var_names_make_unique()

    # Stratified subsample: 200 cells per subclass
    np.random.seed(42)
    indices = []
    for subclass in ref_xenium.obs['subclass_label'].unique():
        mask = ref_xenium.obs['subclass_label'] == subclass
        sc_indices = np.where(mask)[0]
        n_sample = min(200, len(sc_indices))
        chosen = np.random.choice(sc_indices, n_sample, replace=False)
        indices.extend(chosen)

    ref_sub = ref_xenium[sorted(indices)].copy()
    print(f"  Benchmark dataset: {ref_sub.shape}")
    print(f"  Subclasses: {ref_sub.obs['subclass_label'].nunique()}")
    print(f"  Supertypes: {ref_sub.obs['cluster_label'].nunique()}")

    # Save ground truth
    ground_truth = ref_sub.obs[['class_label', 'subclass_label', 'cluster_label']].copy()
    ground_truth.columns = ['true_class', 'true_subclass', 'true_supertype']
    gt_path = os.path.join(OUT_DIR, 'benchmark_ground_truth.csv')
    ground_truth.to_csv(gt_path)

    # Save query h5ad
    query_path = os.path.join(OUT_DIR, 'benchmark_query.h5ad')
    ref_sub.write_h5ad(query_path)
    print(f"  Saved: {query_path}")

    return query_path, gt_path


def run_benchmark(query_path):
    """Run MapMyCells on the benchmark dataset."""
    from cell_type_mapper.cli.map_to_on_the_fly_markers import OnTheFlyMapper

    csv_path = os.path.join(OUT_DIR, 'benchmark_mapmycells.csv')
    json_path = os.path.join(OUT_DIR, 'benchmark_mapmycells.json')

    config = {
        "query_path": query_path,
        "extended_result_path": json_path,
        "csv_result_path": csv_path,
        "precomputed_stats": {"path": STATS_PATH},
        "type_assignment": {
            "normalization": "raw",
            "bootstrap_iteration": 100,
            "bootstrap_factor": 0.5,
            "algorithm": "hierarchical",
        },
        "n_processors": 1,
        "query_markers": {"n_per_utility": 30},
        "reference_markers": {"precomputed_path_list": None},
    }

    print("Running MapMyCells on benchmark...")
    runner = OnTheFlyMapper(args=[], input_data=config)
    runner.run()
    print("  Done!")

    return csv_path


def analyze_results(csv_path, gt_path):
    """Compare MapMyCells predictions to ground truth."""
    df = pd.read_csv(csv_path, comment='#')
    gt = pd.read_csv(gt_path, index_col=0)

    # Extract MapMyCells labels
    for level in ['class', 'subclass', 'supertype']:
        name_col = f'{level}_name'
        prob_col = f'{level}_bootstrapping_probability'
        if name_col in df.columns:
            gt[f'pred_{level}'] = df[name_col].values
            gt[f'{level}_prob'] = df[prob_col].values.astype(float)

    # Compare at each level
    print("\n" + "=" * 60)
    print("BENCHMARK: snRNA-seq reference → Xenium genes → MapMyCells")
    print("=" * 60)

    for level, true_col, pred_col in [
        ('Class', 'true_class', 'pred_class'),
        ('Subclass', 'true_subclass', 'pred_subclass'),
        ('Supertype', 'true_supertype', 'pred_supertype'),
    ]:
        if pred_col not in gt.columns:
            continue
        correct = gt[true_col] == gt[pred_col]
        accuracy = correct.mean()
        prob_col = f'{level.lower()}_prob'

        print(f"\n--- {level} ---")
        print(f"  Overall accuracy: {accuracy:.3f} ({correct.sum()}/{len(correct)})")
        print(f"  Bootstrapping prob: mean={gt[prob_col].mean():.3f}, "
              f"median={gt[prob_col].median():.3f}")

        # Score distribution for correct vs incorrect
        correct_probs = gt.loc[correct, prob_col]
        incorrect_probs = gt.loc[~correct, prob_col]
        print(f"  Correct cells:   mean prob={correct_probs.mean():.3f}, "
              f"median={correct_probs.median():.3f}")
        if len(incorrect_probs) > 0:
            print(f"  Incorrect cells: mean prob={incorrect_probs.mean():.3f}, "
                  f"median={incorrect_probs.median():.3f}, n={len(incorrect_probs)}")

        # Threshold analysis
        print(f"\n  Threshold analysis:")
        print(f"  {'Threshold':>10s} {'% pass':>8s} {'Accuracy':>10s} "
              f"{'N correct':>10s} {'N wrong':>10s}")
        for thresh in [0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            passing = gt[prob_col] >= thresh
            n_pass = passing.sum()
            if n_pass == 0:
                continue
            pass_correct = (gt.loc[passing, true_col] == gt.loc[passing, pred_col])
            acc_pass = pass_correct.mean()
            pct_pass = passing.mean() * 100
            n_right = pass_correct.sum()
            n_wrong = n_pass - n_right
            print(f"  {thresh:>10.1f} {pct_pass:>7.1f}% {acc_pass:>9.3f} "
                  f"{n_right:>10d} {n_wrong:>10d}")

    # Per-subclass accuracy
    print(f"\n--- Per-subclass accuracy ---")
    print(f"{'Subclass':>25s} {'N':>6s} {'Accuracy':>10s} {'Mean prob':>10s}")
    for sc in sorted(gt['true_subclass'].unique()):
        mask = gt['true_subclass'] == sc
        n = mask.sum()
        acc = (gt.loc[mask, 'true_subclass'] == gt.loc[mask, 'pred_subclass']).mean()
        prob = gt.loc[mask, 'subclass_prob'].mean()
        print(f"{sc:>25s} {n:>6d} {acc:>9.3f} {prob:>9.3f}")

    # Save full results
    results_path = os.path.join(OUT_DIR, 'benchmark_results.csv')
    gt.to_csv(results_path)
    print(f"\nSaved: {results_path}")

    return gt


def main():
    # Create benchmark dataset (if not already exists)
    query_path = os.path.join(OUT_DIR, 'benchmark_query.h5ad')
    gt_path = os.path.join(OUT_DIR, 'benchmark_ground_truth.csv')

    if not os.path.exists(query_path):
        query_path, gt_path = create_benchmark_dataset()
    else:
        print(f"Using existing benchmark dataset: {query_path}")

    # Run MapMyCells
    csv_path = os.path.join(OUT_DIR, 'benchmark_mapmycells.csv')
    if not os.path.exists(csv_path):
        csv_path = run_benchmark(query_path)
    else:
        print(f"Using existing MapMyCells results: {csv_path}")

    # Analyze
    gt = analyze_results(csv_path, gt_path)


if __name__ == "__main__":
    main()
