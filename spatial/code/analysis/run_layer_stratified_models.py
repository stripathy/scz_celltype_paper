#!/usr/bin/env python3
"""
Layer-stratified cell type proportion and density analysis.

Tests whether cell type proportions and densities differ between SCZ and
Control within each cortical layer, using per-layer OLS models:

    outcome ~ diagnosis + sex + scale(age) + scale(pmi)

Outcome transforms:
  - Proportion: logit(proportion)
  - Density: log(density + 1)

Multiple testing correction: global FDR (Benjamini-Hochberg) across all
subclass × layer × outcome tests.

Requires:
  output/depth_proportions/layer_counts.csv
  (from build_depth_proportion_input.py)

Output:
  output/depth_proportions/layer_model_results.csv

Usage:
    python3 -u run_layer_stratified_models.py
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR

OUTPUT_DIR = os.path.join(BASE_DIR, "output", "depth_proportions")
CORTICAL_LAYERS = ["L1", "L2/3", "L4", "L5", "L6"]

warnings.filterwarnings("ignore", category=RuntimeWarning)


def fit_layer_models(layer_df, outcome="proportion"):
    """Fit OLS models per (subclass x layer).

    For each subclass in each layer:
      outcome ~ diagnosis + sex + scale(age) + scale(pmi)

    Proportion: logit-transform; Density: log-transform.

    Returns DataFrame with one row per (subclass, layer) test.
    """
    subclasses = sorted(layer_df["subclass_label"].unique())
    results = []

    for layer in CORTICAL_LAYERS:
        layer_data = layer_df[layer_df["layer"] == layer]

        for sc in subclasses:
            sc_data = layer_data[layer_data["subclass_label"] == sc].copy()

            n_samples = sc_data["sample_id"].nunique()
            if n_samples < 10:
                continue

            n_scz = sc_data[sc_data["diagnosis"] == "SCZ"]["sample_id"].nunique()
            n_ctrl = sc_data[sc_data["diagnosis"] == "Control"]["sample_id"].nunique()
            if n_scz < 3 or n_ctrl < 3:
                continue

            # Transform outcome
            if outcome == "density_per_mm2":
                sc_data["y"] = np.log1p(sc_data[outcome])
            else:
                p = sc_data[outcome].clip(1e-6, 1 - 1e-6)
                sc_data["y"] = np.log(p / (1 - p))

            sc_data["age_z"] = ((sc_data["age"] - sc_data["age"].mean())
                                / sc_data["age"].std())
            sc_data["pmi_z"] = ((sc_data["pmi"] - sc_data["pmi"].mean())
                                / sc_data["pmi"].std())

            try:
                formula = ("y ~ C(diagnosis, Treatment('Control'))"
                           " + C(sex) + age_z + pmi_z")
                fit = smf.ols(formula, data=sc_data).fit()

                scz_key = [k for k in fit.params.index if "SCZ" in k]
                if scz_key:
                    coef = fit.params[scz_key[0]]
                    se = fit.bse[scz_key[0]]
                    pval = fit.pvalues[scz_key[0]]
                    tstat = fit.tvalues[scz_key[0]]
                else:
                    coef = se = pval = tstat = np.nan

                ctrl_mean = sc_data[sc_data["diagnosis"] == "Control"][outcome].mean()
                scz_mean = sc_data[sc_data["diagnosis"] == "SCZ"][outcome].mean()
                log2fc = (np.log2(scz_mean / ctrl_mean)
                          if ctrl_mean > 0 else np.nan)

                converged = True

            except Exception as e:
                coef = se = pval = tstat = log2fc = np.nan
                ctrl_mean = scz_mean = np.nan
                converged = False
                print(f"  {layer}/{sc}: failed ({e})")

            results.append({
                "layer": layer,
                "subclass": sc,
                "outcome": outcome,
                "n_samples": n_samples,
                "n_scz": n_scz,
                "n_ctrl": n_ctrl,
                "ctrl_mean": ctrl_mean,
                "scz_mean": scz_mean,
                "log2fc": log2fc,
                "coef": coef,
                "se": se,
                "tstat": tstat,
                "pval": pval,
                "converged": converged,
            })

    return pd.DataFrame(results)


def main():
    t0 = time.time()

    print("Loading layer-level counts...")
    layer_df = pd.read_csv(os.path.join(OUTPUT_DIR, "layer_counts.csv"))
    print(f"  {len(layer_df):,} rows, "
          f"{layer_df['sample_id'].nunique()} samples, "
          f"{layer_df['subclass_label'].nunique()} subclasses")

    all_results = []
    for outcome in ["proportion", "density_per_mm2"]:
        print(f"\n── Layer models: {outcome} ──")
        results = fit_layer_models(layer_df, outcome=outcome)
        all_results.append(results)

        # Quick summary
        n_nom = (results["pval"] < 0.05).sum()
        print(f"  {n_nom} nominally significant (p < 0.05) "
              f"of {len(results)} tests")

    layer_results = pd.concat(all_results, ignore_index=True)

    # Global FDR across ALL layer × subclass × outcome tests
    valid = layer_results["pval"].notna()
    if valid.any():
        _, fdr_vals, _, _ = multipletests(
            layer_results.loc[valid, "pval"], method="fdr_bh")
        layer_results.loc[valid, "fdr"] = fdr_vals

    out_path = os.path.join(OUTPUT_DIR, "layer_model_results.csv")
    layer_results.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")

    # Report top hits
    n_fdr05 = (layer_results.get("fdr", pd.Series()) < 0.05).sum()
    n_fdr10 = (layer_results.get("fdr", pd.Series()) < 0.10).sum()
    n_nom = (layer_results["pval"] < 0.05).sum()
    n_total = valid.sum()
    print(f"\nSummary: {n_fdr05} FDR<0.05, {n_fdr10} FDR<0.10, "
          f"{n_nom} nom p<0.05 (of {n_total} tests)")

    print(f"\n── Top 15 hits ──")
    top = layer_results.dropna(subset=["pval"]).sort_values("pval").head(15)
    for _, row in top.iterrows():
        sig = ("***" if row.get("fdr", 1) < 0.05 else
               ("**" if row.get("fdr", 1) < 0.10 else
                ("*" if row["pval"] < 0.05 else "")))
        print(f"  {row['layer']:5s} / {row['subclass']:25s} ({row['outcome']:15s}): "
              f"log2FC={row['log2fc']:+.3f}  p={row['pval']:.4g}  "
              f"FDR={row.get('fdr', np.nan):.4g} {sig}")

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s")


if __name__ == "__main__":
    main()
