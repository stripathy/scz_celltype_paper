#!/usr/bin/env python3
"""
Build crumblr input CSVs for SEA-AD pseudoprogression analysis.

Creates long-format count tables for crumblr from three data sources:
  1. SEA-AD snRNAseq (89 donors, 1.38M cells)
  2. SEA-AD MERFISH original labels (27 donors, 1.89M cells)
  3. SEA-AD MERFISH reclassified by our pipeline (27 donors, QC-pass)

Each CSV contains: donor, celltype, count, total, CPS, sex, age
where CPS = Continuous Pseudo-progression Score (0.15-0.93).

Covariates match the original SEA-AD paper:
  - CPS (continuous pseudo-progression score) — primary variable of interest
  - Sex (Male/Female)
  - Age at Death (continuous)

Builds both "all cell types" and "neurons only" versions.
Neuron-only restricts to GABAergic + Glutamatergic subclasses and
recomputes totals as sum of neuronal cells per donor.

Output:
  output/crumblr_seaad/crumblr_input_{source}_{level}.csv           (all types)
  output/crumblr_seaad/crumblr_input_{source}_{level}_neurons.csv   (neurons only)
  where source = snrnaseq | merfish_orig | merfish_recl_qc
  and level = subclass | supertype

Usage:
    python3 -u build_seaad_crumblr_input.py
"""

import os
import sys
import pandas as pd
import numpy as np

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
REF_DIR = os.path.join(BASE_DIR, "data", "reference")
OUT_DIR = os.path.join(BASE_DIR, "output", "crumblr_seaad")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Neuronal subclass/supertype mapping ──────────────────────────────
# Non-neuronal subclasses to EXCLUDE for neuron-only analysis
NON_NEURONAL_SUBCLASSES = {
    'Astrocyte', 'Endothelial', 'Microglia-PVM', 'OPC',
    'Oligodendrocyte', 'VLMC',
}

# Supertype prefixes that map to non-neuronal classes
NON_NEURONAL_SUPERTYPE_PREFIXES = (
    'Astro', 'Endo', 'Micro-PVM', 'OPC', 'Oligo', 'VLMC',
    'SMC', 'Pericyte', 'Monocyte', 'Lymphocyte',
)


def load_donor_metadata():
    """Load donor-level metadata from the lean snRNAseq CSV."""
    obs = pd.read_csv(os.path.join(REF_DIR, "SEAAD_MTG_RNAseq_obs_lean.csv"),
                       index_col=0,
                       usecols=['exp_component_name', 'Donor ID',
                                'Continuous Pseudo-progression Score',
                                'Sex', 'Age at Death'])
    donor_meta = obs.groupby('Donor ID').agg({
        'Continuous Pseudo-progression Score': 'first',
        'Sex': 'first',
        'Age at Death': 'first',
    }).reset_index()
    donor_meta = donor_meta.rename(columns={
        'Donor ID': 'donor',
        'Continuous Pseudo-progression Score': 'CPS',
        'Sex': 'sex',
        'Age at Death': 'age',
    })
    # Drop donors without CPS
    donor_meta = donor_meta.dropna(subset=['CPS'])
    print(f"  Donor metadata: {len(donor_meta)} donors with CPS")
    return donor_meta


def is_neuronal_subclass(name):
    """Check if a subclass name is neuronal."""
    return name not in NON_NEURONAL_SUBCLASSES


def is_neuronal_supertype(name):
    """Check if a supertype name is neuronal."""
    return not name.startswith(NON_NEURONAL_SUPERTYPE_PREFIXES)


def build_crumblr_input(count_csv, level_col, donor_meta, source_name,
                         neurons_only=False):
    """Build crumblr input CSV from a pre-computed count file.

    If neurons_only=True, filters to neuronal cell types and recomputes
    totals as the sum of neuronal cells per donor.
    """
    df = pd.read_csv(count_csv)

    # Rename columns to standard format
    if level_col in df.columns:
        df = df.rename(columns={level_col: 'celltype', 'Donor ID': 'donor',
                                'n_cells': 'count', 'n_total': 'total'})
    elif 'Subclass' in df.columns:
        df = df.rename(columns={'Subclass': 'celltype', 'Donor ID': 'donor',
                                'n_cells': 'count', 'n_total': 'total'})
    elif 'Supertype' in df.columns:
        df = df.rename(columns={'Supertype': 'celltype', 'Donor ID': 'donor',
                                'n_cells': 'count', 'n_total': 'total'})
    elif 'Class' in df.columns:
        df = df.rename(columns={'Class': 'celltype', 'Donor ID': 'donor',
                                'n_cells': 'count', 'n_total': 'total'})

    # Keep only the columns we need
    df = df[['donor', 'celltype', 'count', 'total']].copy()

    # Filter to neurons only if requested
    if neurons_only:
        n_before_types = df['celltype'].nunique()
        # Determine level from the column names in original file
        if level_col in ('Subclass', 'subclass'):
            df = df[df['celltype'].apply(is_neuronal_subclass)]
        else:  # Supertype
            df = df[df['celltype'].apply(is_neuronal_supertype)]

        n_after_types = df['celltype'].nunique()
        print(f"      Neuron filter: {n_after_types}/{n_before_types} types kept")

        # Recompute totals as sum of neuronal cells per donor
        neuron_totals = df.groupby('donor')['count'].sum().reset_index()
        neuron_totals = neuron_totals.rename(columns={'count': 'neuron_total'})
        df = df.merge(neuron_totals, on='donor', how='left')
        df['total'] = df['neuron_total']
        df = df.drop(columns=['neuron_total'])

    # Merge with donor metadata
    df = df.merge(donor_meta, on='donor', how='inner')

    # Drop donors without metadata
    n_before = df['donor'].nunique()
    df = df.dropna(subset=['CPS', 'sex', 'age'])
    n_after = df['donor'].nunique()

    n_types = df['celltype'].nunique()
    suffix = " [neurons]" if neurons_only else ""
    print(f"    {source_name}{suffix}: {n_after} donors, {n_types} types, {len(df)} rows")
    if n_before != n_after:
        print(f"      (dropped {n_before - n_after} donors missing metadata)")

    return df


def main():
    print("=" * 70)
    print("Build crumblr inputs for SEA-AD pseudoprogression analysis")
    print("=" * 70)

    # Load donor metadata
    print("\nLoading donor metadata...")
    donor_meta = load_donor_metadata()

    # Define all input sources and levels
    sources = {
        'snrnaseq': {
            'subclass': os.path.join(REF_DIR, "SEAAD_snrnaseq_counts_by_subclass.csv"),
            'supertype': os.path.join(REF_DIR, "SEAAD_snrnaseq_counts_by_supertype.csv"),
        },
        'merfish_orig': {
            'subclass': os.path.join(REF_DIR, "SEAAD_merfish_original_counts_by_subclass.csv"),
            'supertype': os.path.join(REF_DIR, "SEAAD_merfish_original_counts_by_supertype.csv"),
        },
        'merfish_recl_qc': {
            'subclass': os.path.join(REF_DIR, "SEAAD_merfish_reclassified_qcpass_counts_by_subclass.csv"),
            'supertype': os.path.join(REF_DIR, "SEAAD_merfish_reclassified_qcpass_counts_by_supertype.csv"),
        },
        'merfish_depth': {
            'subclass': os.path.join(REF_DIR, "SEAAD_merfish_depth_annotated_counts_by_subclass.csv"),
            'supertype': os.path.join(REF_DIR, "SEAAD_merfish_depth_annotated_counts_by_supertype.csv"),
        },
    }

    for source_name, levels in sources.items():
        print(f"\n  Source: {source_name}")
        for level_name, csv_path in levels.items():
            if not os.path.exists(csv_path):
                print(f"    WARNING: {csv_path} not found, skipping")
                continue

            # Build all-types version
            df = build_crumblr_input(csv_path, level_name.capitalize(),
                                      donor_meta, f"{source_name}/{level_name}")
            df = df.sort_values(['donor', 'celltype']).reset_index(drop=True)
            out_path = os.path.join(OUT_DIR, f"crumblr_input_{source_name}_{level_name}.csv")
            df.to_csv(out_path, index=False)
            print(f"      -> {out_path}")

            # Build neurons-only version
            df_n = build_crumblr_input(csv_path, level_name.capitalize(),
                                        donor_meta, f"{source_name}/{level_name}",
                                        neurons_only=True)
            df_n = df_n.sort_values(['donor', 'celltype']).reset_index(drop=True)
            out_n = os.path.join(OUT_DIR, f"crumblr_input_{source_name}_{level_name}_neurons.csv")
            df_n.to_csv(out_n, index=False)
            print(f"      -> {out_n}")

    # Summary
    print(f"\n{'='*70}")
    print("Input files created:")
    for f in sorted(os.listdir(OUT_DIR)):
        if f.startswith('crumblr_input_'):
            fpath = os.path.join(OUT_DIR, f)
            df = pd.read_csv(fpath)
            print(f"  {f}: {df['donor'].nunique()} donors, {df['celltype'].nunique()} types")
    print("Done!")


if __name__ == "__main__":
    main()
