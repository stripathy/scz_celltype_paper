#!/usr/bin/env python3
"""
Step 2: Cell type annotation using Allen Institute's MapMyCells (HANN algorithm).

Uses bootstrapped HANN (Hierarchical Approximate Nearest Neighbor) mapping
against the SEA-AD MTG taxonomy to assign class, subclass, and supertype
labels to each cell. HANN performs a full-tree correlation walk per bootstrap
iteration, tallying votes at the leaf (supertype) level.

For each sample:
  1. Load existing h5ad (from steps 00+01, with expression + QC columns)
  2. Subset to qc_pass == True cells
  3. Convert gene symbols to Ensembl IDs (reference uses Ensembl)
  4. Run MapMyCells HANN mapping (100 bootstrap iterations)
  5. Parse HDF5 output: extract labels + confidence at class/subclass/supertype
  6. Write labels back to h5ad (QC-fail cells get 'Unassigned' + 0.0 confidence)

Requires:
  - Python 3.10+
  - cell_type_mapper >= 1.7.0:
      pip install cell_type_mapper@git+https://github.com/AllenInstitute/cell_type_mapper
  - SEA-AD MTG precomputed stats (~1-2 GB, download once):
      wget -P data/reference/ https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/
        mapmycells/SEAAD/20240831/precomputed_stats.20231120.sea_ad.MTG.h5

Usage:
    python3 -u 02_run_mapmycells.py [--sample SAMPLE_ID] [--n-workers N]
"""

import os
import sys
import time
import json
import argparse
import tempfile
import numpy as np
import pandas as pd
import anndata as ad
import h5py

# Pipeline config
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_config import (
    H5AD_DIR, OUTPUT_DIR, RAW_DIR, MODULES_DIR,
    PRECOMPUTED_STATS_PATH, GENE_MAPPING_PATH,
    TAXONOMY_LEVELS,
    MAPMYCELLS_BOOTSTRAP_ITER, MAPMYCELLS_BOOTSTRAP_FACTOR,
    MAPMYCELLS_N_PER_UTILITY,
)

# Modules
sys.path.insert(0, MODULES_DIR)
from loading import discover_samples


# ──────────────────────────────────────────────────────────────────────
# Gene mapping utilities
# ──────────────────────────────────────────────────────────────────────

def load_gene_mapping():
    """Load gene symbol -> Ensembl ID mapping from JSON."""
    with open(GENE_MAPPING_PATH) as f:
        return json.load(f)


def convert_genes_to_ensembl(adata, gene_mapping):
    """
    Convert gene symbols in adata.var_names to Ensembl IDs.

    Drops genes that can't be mapped. Returns a new AnnData with
    Ensembl IDs as var_names.
    """
    mappable = [g for g in adata.var_names if g in gene_mapping]
    adata_sub = adata[:, mappable].copy()
    adata_sub.var_names = [gene_mapping[g] for g in mappable]
    adata_sub.var_names_make_unique()
    return adata_sub


# ──────────────────────────────────────────────────────────────────────
# HANN taxonomy helpers
# ──────────────────────────────────────────────────────────────────────

def _load_taxonomy_tree():
    """Load the taxonomy tree from precomputed stats HDF5."""
    with h5py.File(PRECOMPUTED_STATS_PATH, 'r') as f:
        tree = json.loads(f['taxonomy_tree'][()].decode())
    return tree


def _build_taxonomy_lookups(tree):
    """
    Build lookups from the taxonomy tree:
      - supertype_id -> supertype_name
      - supertype_id -> subclass_name
      - supertype_id -> class_name

    The tree structure has:
      CCN..._CLAS: {class_id: [list of subclass_ids]}
      CCN..._SUBC: {subclass_id: [list of supertype_ids]}
      CCN..._SUPT: {supertype_id: [list of cell barcodes]}
    """
    hierarchy = tree['hierarchy']
    nm = tree['name_mapper']

    clas_level = hierarchy[0]
    subc_level = hierarchy[1]
    supt_level = hierarchy[2]

    clas_data = tree[clas_level]
    subc_data = tree[subc_level]

    supt_to_name = {sid: info['name']
                     for sid, info in nm[supt_level].items()}
    supt_to_subclass = {}
    supt_to_class = {}

    for clas_id, subc_ids in clas_data.items():
        clas_name = nm[clas_level][clas_id]['name']
        for subc_id in subc_ids:
            subc_name = nm[subc_level][subc_id]['name']
            for supt_id in subc_data.get(subc_id, []):
                supt_to_subclass[supt_id] = subc_name
                supt_to_class[supt_id] = clas_name

    return supt_to_name, supt_to_subclass, supt_to_class


# ──────────────────────────────────────────────────────────────────────
# MapMyCells runner
# ──────────────────────────────────────────────────────────────────────

def run_mapmycells_on_sample(query_h5ad_path, output_dir, n_processors=1):
    """
    Run MapMyCells HANN mapping on a single h5ad file.

    Returns hdf5_path for the results.
    """
    from cell_type_mapper.cli.map_to_on_the_fly_markers import OnTheFlyMapper

    hdf5_path = os.path.join(output_dir, "hann_output.h5")

    config = {
        "query_path": query_h5ad_path,
        "hdf5_result_path": hdf5_path,
        "precomputed_stats": {
            "path": PRECOMPUTED_STATS_PATH,
        },
        "type_assignment": {
            "normalization": "raw",
            "bootstrap_iteration": MAPMYCELLS_BOOTSTRAP_ITER,
            "bootstrap_factor": MAPMYCELLS_BOOTSTRAP_FACTOR,
            "algorithm": "hann",
        },
        "n_processors": n_processors,
        "query_markers": {
            "n_per_utility": MAPMYCELLS_N_PER_UTILITY,
        },
        "reference_markers": {
            "precomputed_path_list": None,
        },
    }

    runner = OnTheFlyMapper(args=[], input_data=config)
    runner.run()

    return hdf5_path


def parse_mapmycells_output(hdf5_path):
    """
    Parse HANN HDF5 output and extract human-readable labels + confidence.

    The HDF5 contains:
      - votes: (n_cells x n_clusters) vote counts across bootstrap iterations
      - correlation: (n_cells x n_clusters) average correlation per cluster
      - cluster_identifiers: cluster IDs (supertype-level taxonomy IDs)

    For each cell, assigns the cluster with the most votes (correlation as
    tiebreaker). Computes confidence at supertype, subclass, and class levels
    by aggregating votes up the taxonomy hierarchy.

    Returns a DataFrame with standard column names:
      class_label, subclass_label, supertype_label,
      class_label_confidence, subclass_label_confidence, supertype_label_confidence
    """
    tree = _load_taxonomy_tree()
    supt_to_name, supt_to_subclass, supt_to_class = _build_taxonomy_lookups(tree)

    with h5py.File(hdf5_path, 'r') as f:
        votes = f['votes'][:]
        correlation = f['correlation'][:]
        cluster_ids = [c.decode() if isinstance(c, bytes) else c
                       for c in f['cluster_identifiers'][:]]

    n_cells = votes.shape[0]

    # For each cell, pick cluster with most votes (correlation as tiebreaker)
    score = votes.astype(float) + correlation * 1e-6
    best_idx = np.argmax(score, axis=1)

    # Supertype confidence = fraction of votes for winning cluster
    max_votes = votes[np.arange(n_cells), best_idx]
    total_votes = votes.sum(axis=1)
    supt_confidence = np.where(total_votes > 0,
                                max_votes / total_votes,
                                0.0).astype(np.float32)

    best_cluster_ids = [cluster_ids[i] for i in best_idx]

    result = pd.DataFrame(index=range(n_cells))

    # Supertype (leaf level)
    result["supertype_label"] = [supt_to_name.get(cid, cid)
                                  for cid in best_cluster_ids]
    result["supertype_label_confidence"] = supt_confidence

    # Subclass (inferred from supertype via taxonomy)
    subclass_labels = [supt_to_subclass.get(cid, "Unknown")
                       for cid in best_cluster_ids]
    result["subclass_label"] = subclass_labels

    # Subclass confidence: sum votes for all supertypes in same subclass
    cluster_subclasses = [supt_to_subclass.get(cid, "") for cid in cluster_ids]
    unique_subclasses = sorted(set(cluster_subclasses))
    sub_to_idx = {s: i for i, s in enumerate(unique_subclasses)}

    n_subclasses = len(unique_subclasses)
    sub_votes = np.zeros((n_cells, n_subclasses), dtype=np.float32)
    for j, sc in enumerate(cluster_subclasses):
        if sc in sub_to_idx:
            sub_votes[:, sub_to_idx[sc]] += votes[:, j]

    subclass_confidence = np.zeros(n_cells, dtype=np.float32)
    for i in range(n_cells):
        sc = subclass_labels[i]
        if sc in sub_to_idx and total_votes[i] > 0:
            subclass_confidence[i] = sub_votes[i, sub_to_idx[sc]] / total_votes[i]
    result["subclass_label_confidence"] = subclass_confidence

    # Class (inferred from supertype via taxonomy)
    class_labels = [supt_to_class.get(cid, "Unknown")
                    for cid in best_cluster_ids]
    result["class_label"] = class_labels

    # Class confidence: sum votes for all supertypes in same class
    cluster_classes = [supt_to_class.get(cid, "") for cid in cluster_ids]
    unique_classes = sorted(set(cluster_classes))
    cls_to_idx = {c: i for i, c in enumerate(unique_classes)}
    n_classes = len(unique_classes)
    cls_votes = np.zeros((n_cells, n_classes), dtype=np.float32)
    for j, cc in enumerate(cluster_classes):
        if cc in cls_to_idx:
            cls_votes[:, cls_to_idx[cc]] += votes[:, j]

    class_confidence = np.zeros(n_cells, dtype=np.float32)
    for i in range(n_cells):
        cc = class_labels[i]
        if cc in cls_to_idx and total_votes[i] > 0:
            class_confidence[i] = cls_votes[i, cls_to_idx[cc]] / total_votes[i]
    result["class_label_confidence"] = class_confidence

    print(f"    HANN: {n_cells:,} cells, {len(cluster_ids)} clusters")
    print(f"    Mean confidence — supertype: {supt_confidence.mean():.3f}, "
          f"subclass: {subclass_confidence.mean():.3f}, "
          f"class: {class_confidence.mean():.3f}")

    return result


# ──────────────────────────────────────────────────────────────────────
# Per-sample processing
# ──────────────────────────────────────────────────────────────────────

def process_one_sample(sample_info, gene_mapping):
    """
    Run MapMyCells HANN label transfer on one sample.

    Reads the existing h5ad (with expression matrix + QC columns from
    steps 00 and 01), runs HANN mapping on QC-pass cells, and writes
    the label columns back.
    """
    sid = sample_info["sample_id"]
    h5ad_path = os.path.join(H5AD_DIR, f"{sid}_annotated.h5ad")
    t0 = time.time()

    try:
        # 1. Load existing h5ad
        print(f"  {sid}: Loading h5ad...")
        adata = ad.read_h5ad(h5ad_path)
        n_total = adata.shape[0]

        # 2. Get QC mask
        if "qc_pass" not in adata.obs.columns:
            print(f"  {sid}: WARNING - no qc_pass column, using all cells")
            qc_mask = np.ones(n_total, dtype=bool)
        else:
            qc_mask = adata.obs["qc_pass"].values.astype(bool)

        n_pass = int(qc_mask.sum())
        n_fail = n_total - n_pass
        print(f"  {sid}: {n_total:,} total, {n_pass:,} QC pass, "
              f"{n_fail:,} QC fail ({100*n_fail/n_total:.1f}%)")

        # 3. Subset to QC-pass cells and convert genes to Ensembl
        adata_pass = adata[qc_mask].copy()

        with tempfile.TemporaryDirectory() as tmpdir:
            adata_ensembl = convert_genes_to_ensembl(adata_pass, gene_mapping)
            n_mapped = adata_ensembl.shape[1]
            print(f"  {sid}: Mapped {n_mapped}/{adata_pass.shape[1]} genes to Ensembl IDs")

            # 4. Save temp h5ad for MapMyCells (it needs a file path)
            tmp_h5ad = os.path.join(tmpdir, f"{sid}_qcpass.h5ad")
            adata_ensembl.write_h5ad(tmp_h5ad)

            # 5. Run MapMyCells HANN
            print(f"  {sid}: Running MapMyCells HANN mapping...")
            t_map = time.time()
            hdf5_path = run_mapmycells_on_sample(
                tmp_h5ad, tmpdir, n_processors=1
            )
            print(f"  {sid}: MapMyCells done in {time.time()-t_map:.0f}s")

            # 6. Parse results
            labels_df = parse_mapmycells_output(hdf5_path)

            # Log summary
            for level in TAXONOMY_LEVELS:
                col = f"{level}_label"
                conf_col = f"{level}_label_confidence"
                if col in labels_df.columns:
                    n_unique = labels_df[col].nunique()
                    mean_conf = labels_df[conf_col].mean() if conf_col in labels_df.columns else 0
                    print(f"  {sid}: {level}: {n_unique} types, "
                          f"mean confidence={mean_conf:.3f}")

        # 7. Write labels into the full h5ad
        # Initialize all cells as Unassigned
        for level in TAXONOMY_LEVELS:
            adata.obs[f"{level}_label"] = "Unassigned"
            adata.obs[f"{level}_label_confidence"] = np.float32(0.0)

        # Copy labels for QC-pass cells
        pass_indices = np.where(qc_mask)[0]
        for col in labels_df.columns:
            if col in adata.obs.columns:
                adata.obs.iloc[pass_indices,
                               adata.obs.columns.get_loc(col)] = \
                    labels_df[col].values

        # 8. Save updated h5ad
        adata.write_h5ad(h5ad_path)

        elapsed = time.time() - t0
        print(f"  {sid}: Done in {elapsed:.0f}s")

        return {
            "sample_id": sid, "status": "success",
            "n_total": n_total, "n_pass": n_pass,
            "n_fail": n_fail, "time": elapsed,
        }

    except Exception as e:
        import traceback
        elapsed = time.time() - t0
        print(f"  {sid}: FAILED after {elapsed:.0f}s - {e}")
        traceback.print_exc()
        return {
            "sample_id": sid, "status": "error",
            "error": str(e), "time": elapsed,
        }


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Step 2: Run MapMyCells HANN label transfer on Xenium samples"
    )
    parser.add_argument("--sample", type=str, default=None,
                        help="Process only this sample ID (for testing)")
    parser.add_argument("--n-workers", type=int, default=1,
                        help="Number of MapMyCells internal processors")
    args = parser.parse_args()

    t_start = time.time()

    # Check prerequisites
    if not os.path.exists(PRECOMPUTED_STATS_PATH):
        print(f"ERROR: SEA-AD MTG precomputed stats not found at:")
        print(f"  {PRECOMPUTED_STATS_PATH}")
        print(f"\nDownload with:")
        print(f"  wget -P {os.path.dirname(PRECOMPUTED_STATS_PATH)} \\")
        print(f"    https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/"
              f"mapmycells/SEAAD/20240831/precomputed_stats.20231120.sea_ad.MTG.h5")
        sys.exit(1)

    if not os.path.exists(GENE_MAPPING_PATH):
        print(f"ERROR: Gene mapping not found at {GENE_MAPPING_PATH}")
        print("Run the gene mapping script first.")
        sys.exit(1)

    gene_mapping = load_gene_mapping()
    print(f"Loaded gene mapping: {len(gene_mapping)} symbol->Ensembl mappings")

    # Discover samples from raw data directory
    all_samples = discover_samples(RAW_DIR)
    print(f"Found {len(all_samples)} samples")

    if args.sample:
        all_samples = [s for s in all_samples if s["sample_id"] == args.sample]
        if not all_samples:
            print(f"ERROR: Sample '{args.sample}' not found")
            sys.exit(1)
        print(f"Processing single sample: {args.sample}")

    print(f"\nMethod: MapMyCells HANN mapping (SEA-AD MTG taxonomy)")
    print(f"Reference: {PRECOMPUTED_STATS_PATH}")
    print(f"{'='*60}")

    # Process sequentially (MapMyCells has its own internal parallelism)
    results = []
    for i, sample_info in enumerate(all_samples):
        print(f"\n[{i+1}/{len(all_samples)}] Processing {sample_info['sample_id']}...")
        result = process_one_sample(sample_info, gene_mapping=gene_mapping)
        results.append(result)

    # Summary
    total_time = time.time() - t_start
    n_ok = sum(1 for r in results if r["status"] == "success")
    total_cells = sum(r.get("n_total", 0) for r in results)
    total_pass = sum(r.get("n_pass", 0) for r in results)

    print(f"\n{'='*60}")
    print(f"SUMMARY - {total_time:.0f}s total")
    print(f"{'='*60}")
    print(f"  Success: {n_ok}/{len(results)}")
    if total_cells > 0:
        print(f"  Total cells: {total_cells:,}")
        print(f"  QC pass: {total_pass:,} ({100*total_pass/total_cells:.1f}%)")

    for r in results:
        status = "OK" if r["status"] == "success" else "FAIL"
        t = r.get("time", 0)
        n = r.get("n_total", 0)
        print(f"  {status} {r['sample_id']:10s}: {n:>6,} cells, {t:>4.0f}s")

    for r in results:
        if r["status"] == "error":
            print(f"  FAILED: {r['sample_id']}: {r.get('error', 'unknown')}")

    print(f"\nTotal time: {total_time:.0f}s")


if __name__ == "__main__":
    main()
