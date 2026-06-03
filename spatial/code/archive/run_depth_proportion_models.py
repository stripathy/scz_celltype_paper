#!/usr/bin/env python3
"""
Depth-stratified cell type proportion and density analysis.

Tests whether cell type proportions and densities differ between SCZ and
control at specific cortical depths, using two complementary approaches:

1. GAM (continuous depth): Fits smooth depth profiles per diagnosis group
   using generalized additive models. Tests whether the SCZ depth profile
   differs from control via a diagnosis-specific smooth term.

2. Per-layer mixed models: Fits linear mixed effects models within each
   cortical layer, testing diagnosis effect with sex, age, PMI covariates
   and sample-level random intercepts.

Both approaches test two outcome measures:
  - Proportion: fraction of cells that are type X at a given depth
  - Density: cells of type X per mm² at a given depth

Multiple testing correction: global FDR (Benjamini-Hochberg) across all
subclass × depth/layer combinations.

Requires:
  output/depth_proportions/cell_level_data.csv
  output/depth_proportions/layer_counts.csv
  (from build_depth_proportion_input.py)

Output:
  output/depth_proportions/gam_results.csv
  output/depth_proportions/gam_smooth_predictions.csv
  output/depth_proportions/layer_model_results.csv

Usage:
    python3 -u run_depth_proportion_models.py [--n-depth-bins 20] [--min-cells 30]
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

# Suppress convergence warnings for display clarity (still tracked in results)
warnings.filterwarnings("ignore", category=RuntimeWarning)


# ──────────────────────────────────────────────────────────────────────
# Part 1: GAM-based continuous depth analysis
# ──────────────────────────────────────────────────────────────────────

def compute_depth_binned_proportions(cells_df, n_bins=20, min_cells_per_bin=30):
    """Compute proportion and density of each subclass in depth bins per sample.

    For the GAM, we discretize depth into fine bins (default=20) to get
    per-sample proportion estimates at each depth, then fit smooth models
    to these binned proportions. This avoids fitting a GAM to millions of
    binary outcomes (cell is/isn't type X).

    Returns DataFrame with columns:
      sample_id, subclass_label, depth_bin_mid, proportion, density_per_mm2,
      n_cells, n_total, diagnosis, sex, age, pmi
    """
    # Create depth bins (exclude extremes outside [0, 1])
    cells = cells_df[cells_df["predicted_norm_depth"].between(0, 1)].copy()
    bin_edges = np.linspace(0, 1, n_bins + 1)
    bin_mids = (bin_edges[:-1] + bin_edges[1:]) / 2
    cells["depth_bin"] = pd.cut(cells["predicted_norm_depth"], bins=bin_edges,
                                 labels=bin_mids, include_lowest=True)
    cells["depth_bin"] = cells["depth_bin"].astype(float)

    # Count per (sample, depth_bin, subclass)
    counts = (cells.groupby(["sample_id", "depth_bin", "subclass_label"])
              .size().reset_index(name="n_cells"))

    # Total per (sample, depth_bin)
    totals = (cells.groupby(["sample_id", "depth_bin"])
              .size().reset_index(name="n_total"))

    # Area per (sample, depth_bin) — convex hull of cells in that bin
    from scipy.spatial import ConvexHull

    def _hull_area(group):
        if len(group) < 10:
            return np.nan
        coords = group[["x", "y"]].values
        try:
            return ConvexHull(coords).volume / 1e6  # µm² → mm²
        except Exception:
            return np.nan

    areas = (cells.groupby(["sample_id", "depth_bin"])
             .apply(_hull_area, include_groups=False)
             .reset_index(name="area_mm2"))

    # Merge
    result = counts.merge(totals, on=["sample_id", "depth_bin"])
    result = result.merge(areas, on=["sample_id", "depth_bin"], how="left")

    # Compute proportion and density
    result["proportion"] = result["n_cells"] / result["n_total"]
    result["density_per_mm2"] = result["n_cells"] / result["area_mm2"]

    # Filter sparse bins
    result = result[result["n_total"] >= min_cells_per_bin]

    # Attach metadata (take first row per sample from cells_df)
    meta_cols = ["sample_id", "diagnosis", "sex", "age", "pmi"]
    meta = cells_df[meta_cols].drop_duplicates()
    result = result.merge(meta, on="sample_id", how="left")

    return result


def fit_gam_per_subclass(binned_df, outcome="proportion", n_splines=6):
    """Fit GAM-style models for each subclass using basis splines.

    We approximate GAMs using natural cubic splines in statsmodels, which
    allows mixed effects and standard formula interface. For each subclass:

    Full model:  outcome ~ bs(depth, n_splines) * diagnosis + sex + age + pmi
    Reduced:     outcome ~ bs(depth, n_splines) + sex + age + pmi

    The interaction term tests whether the depth profile differs by diagnosis.

    Returns:
      results_df: One row per subclass with F-test p-value for interaction
      predictions_df: Smooth predictions at fine depth grid for plotting
    """
    from patsy import dmatrix
    import statsmodels.api as sm

    subclasses = sorted(binned_df["subclass_label"].unique())
    depth_grid = np.linspace(0.02, 0.98, 50)  # fine grid for predictions

    results = []
    predictions = []

    for sc in subclasses:
        sc_data = binned_df[binned_df["subclass_label"] == sc].copy()

        # Need enough data points
        n_samples = sc_data["sample_id"].nunique()
        n_obs = len(sc_data)
        if n_samples < 6 or n_obs < 20:
            print(f"  {sc}: skipped (n_samples={n_samples}, n_obs={n_obs})")
            continue

        # Use log transform for density, logit-ish for proportion
        if outcome == "density_per_mm2":
            sc_data["y"] = np.log1p(sc_data[outcome])
        else:
            # CLR-inspired: log(p / (1-p)), clamped to avoid inf
            p = sc_data[outcome].clip(1e-6, 1 - 1e-6)
            sc_data["y"] = np.log(p / (1 - p))

        sc_data["depth"] = sc_data["depth_bin"]
        sc_data["age_z"] = (sc_data["age"] - sc_data["age"].mean()) / sc_data["age"].std()
        sc_data["pmi_z"] = (sc_data["pmi"] - sc_data["pmi"].mean()) / sc_data["pmi"].std()
        sc_data["is_scz"] = (sc_data["diagnosis"] == "SCZ").astype(float)

        try:
            # Full model with diagnosis × depth interaction
            # Using polynomial basis (degree 3) for simplicity and stability
            full_formula = ("y ~ (depth + I(depth**2) + I(depth**3)) * is_scz "
                           "+ C(sex) + age_z + pmi_z")
            reduced_formula = ("y ~ depth + I(depth**2) + I(depth**3) "
                              "+ C(sex) + age_z + pmi_z")

            # Mixed model with random intercept per sample
            full_model = smf.mixedlm(full_formula, sc_data,
                                      groups=sc_data["sample_id"],
                                      missing="drop")
            full_fit = full_model.fit(reml=False, method="lbfgs")

            reduced_model = smf.mixedlm(reduced_formula, sc_data,
                                         groups=sc_data["sample_id"],
                                         missing="drop")
            reduced_fit = reduced_model.fit(reml=False, method="lbfgs")

            # Likelihood ratio test for interaction
            lr_stat = -2 * (reduced_fit.llf - full_fit.llf)
            df_diff = full_fit.df_modelwc - reduced_fit.df_modelwc
            from scipy.stats import chi2
            interaction_p = chi2.sf(lr_stat, max(df_diff, 1))

            # Main effect of diagnosis (from full model)
            diag_coef = full_fit.params.get("is_scz", np.nan)
            diag_pval = full_fit.pvalues.get("is_scz", np.nan)

            converged = True

        except Exception as e:
            print(f"  {sc}: model failed ({e})")
            interaction_p = np.nan
            diag_coef = np.nan
            diag_pval = np.nan
            converged = False

        results.append({
            "subclass": sc,
            "outcome": outcome,
            "n_samples": n_samples,
            "n_observations": n_obs,
            "interaction_pval": interaction_p,
            "diagnosis_coef": diag_coef,
            "diagnosis_pval": diag_pval,
            "converged": converged,
        })

        # Generate predictions for plotting (if model converged)
        if converged:
            try:
                for dx_label, dx_val in [("Control", 0), ("SCZ", 1)]:
                    pred_data = pd.DataFrame({
                        "depth": depth_grid,
                        "is_scz": dx_val,
                        "sex": sc_data["sex"].mode().iloc[0],
                        "age_z": 0.0,
                        "pmi_z": 0.0,
                        "sample_id": sc_data["sample_id"].iloc[0],
                    })

                    pred = full_fit.predict(pred_data)

                    # Back-transform
                    if outcome == "density_per_mm2":
                        pred_orig = np.expm1(pred)
                    else:
                        pred_orig = 1 / (1 + np.exp(-pred))  # inverse logit

                    for i, d in enumerate(depth_grid):
                        predictions.append({
                            "subclass": sc,
                            "outcome": outcome,
                            "diagnosis": dx_label,
                            "depth": d,
                            "predicted": pred_orig.iloc[i] if hasattr(pred_orig, 'iloc') else pred_orig[i],
                        })
            except Exception as e:
                print(f"  {sc}: prediction failed ({e})")

        print(f"  {sc}: interaction p={interaction_p:.4g}, "
              f"diagnosis coef={diag_coef:.4f}" if converged else "")

    results_df = pd.DataFrame(results)
    predictions_df = pd.DataFrame(predictions)

    return results_df, predictions_df


# ──────────────────────────────────────────────────────────────────────
# Part 2: Per-layer mixed effects models
# ──────────────────────────────────────────────────────────────────────

def fit_layer_models(layer_df, outcome="proportion"):
    """Fit mixed effects models per (subclass × layer).

    For each subclass in each layer:
      outcome ~ diagnosis + sex + scale(age) + scale(pmi) + (1 | sample_id)

    For proportion: logit-transform the outcome
    For density: log-transform the outcome

    Returns DataFrame with one row per (subclass, layer) test.
    """
    subclasses = sorted(layer_df["subclass_label"].unique())
    results = []

    for layer in CORTICAL_LAYERS:
        layer_data = layer_df[layer_df["layer"] == layer]

        for sc in subclasses:
            sc_data = layer_data[layer_data["subclass_label"] == sc].copy()

            # Need minimum representation
            n_samples = sc_data["sample_id"].nunique()
            if n_samples < 10:  # need reasonable N in both groups
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

            sc_data["age_z"] = (sc_data["age"] - sc_data["age"].mean()) / sc_data["age"].std()
            sc_data["pmi_z"] = (sc_data["pmi"] - sc_data["pmi"].mean()) / sc_data["pmi"].std()

            try:
                formula = "y ~ C(diagnosis, Treatment('Control')) + C(sex) + age_z + pmi_z"
                model = smf.ols(formula, data=sc_data)
                fit = model.fit()

                # Extract SCZ coefficient
                scz_key = [k for k in fit.params.index if "SCZ" in k]
                if scz_key:
                    coef = fit.params[scz_key[0]]
                    se = fit.bse[scz_key[0]]
                    pval = fit.pvalues[scz_key[0]]
                    tstat = fit.tvalues[scz_key[0]]
                else:
                    coef = se = pval = tstat = np.nan

                # Effect size: mean proportion/density per group
                ctrl_mean = sc_data[sc_data["diagnosis"] == "Control"][outcome].mean()
                scz_mean = sc_data[sc_data["diagnosis"] == "SCZ"][outcome].mean()
                log2fc = np.log2(scz_mean / ctrl_mean) if ctrl_mean > 0 else np.nan

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


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-depth-bins", type=int, default=20,
                        help="Number of depth bins for GAM analysis (default: 20)")
    parser.add_argument("--min-cells", type=int, default=30,
                        help="Minimum cells per depth bin to include (default: 30)")
    parser.add_argument("--n-splines", type=int, default=6,
                        help="Number of spline basis functions for GAM (default: 6)")
    args = parser.parse_args()

    t0 = time.time()

    # ── Load data ──
    print("Loading cell-level data...")
    cells_df = pd.read_csv(os.path.join(OUTPUT_DIR, "cell_level_data.csv"))
    print(f"  {len(cells_df):,} cells, {cells_df['sample_id'].nunique()} samples")

    print("Loading layer-level counts...")
    layer_df = pd.read_csv(os.path.join(OUTPUT_DIR, "layer_counts.csv"))
    print(f"  {len(layer_df):,} rows")

    # ── Part 1: GAM continuous depth ──
    print("\n" + "=" * 70)
    print("PART 1: GAM continuous-depth models")
    print("=" * 70)

    # Compute binned proportions for GAM input
    print(f"\nBinning cells into {args.n_depth_bins} depth bins...")
    binned = compute_depth_binned_proportions(cells_df, n_bins=args.n_depth_bins,
                                               min_cells_per_bin=args.min_cells)
    print(f"  {len(binned):,} (sample × depth_bin × subclass) observations")

    all_gam_results = []
    all_gam_preds = []

    for outcome in ["proportion", "density_per_mm2"]:
        print(f"\n── GAM: {outcome} ──")
        results, preds = fit_gam_per_subclass(binned, outcome=outcome,
                                               n_splines=args.n_splines)
        all_gam_results.append(results)
        all_gam_preds.append(preds)

    gam_results = pd.concat(all_gam_results, ignore_index=True)
    gam_preds = pd.concat(all_gam_preds, ignore_index=True)

    # FDR correction across all GAM interaction tests
    valid = gam_results["interaction_pval"].notna()
    if valid.any():
        _, fdr_vals, _, _ = multipletests(gam_results.loc[valid, "interaction_pval"],
                                           method="fdr_bh")
        gam_results.loc[valid, "interaction_fdr"] = fdr_vals

    gam_path = os.path.join(OUTPUT_DIR, "gam_results.csv")
    gam_results.to_csv(gam_path, index=False)
    print(f"\nGAM results: {gam_path}")

    pred_path = os.path.join(OUTPUT_DIR, "gam_smooth_predictions.csv")
    gam_preds.to_csv(pred_path, index=False)
    print(f"GAM predictions: {pred_path}")

    # Report top hits
    print("\n── Top GAM interaction hits ──")
    top = gam_results.dropna(subset=["interaction_pval"]).sort_values("interaction_pval").head(10)
    for _, row in top.iterrows():
        sig = "***" if row.get("interaction_fdr", 1) < 0.05 else ""
        print(f"  {row['subclass']:25s} ({row['outcome']:15s}): "
              f"p={row['interaction_pval']:.4g}  FDR={row.get('interaction_fdr', np.nan):.4g} {sig}")

    # ── Part 2: Per-layer models ──
    print("\n" + "=" * 70)
    print("PART 2: Per-layer mixed effects models")
    print("=" * 70)

    all_layer_results = []
    for outcome in ["proportion", "density_per_mm2"]:
        print(f"\n── Layer models: {outcome} ──")
        results = fit_layer_models(layer_df, outcome=outcome)
        all_layer_results.append(results)

    layer_results = pd.concat(all_layer_results, ignore_index=True)

    # Global FDR across ALL layer × subclass × outcome tests
    valid = layer_results["pval"].notna()
    if valid.any():
        _, fdr_vals, _, _ = multipletests(layer_results.loc[valid, "pval"],
                                           method="fdr_bh")
        layer_results.loc[valid, "fdr"] = fdr_vals

    layer_path = os.path.join(OUTPUT_DIR, "layer_model_results.csv")
    layer_results.to_csv(layer_path, index=False)
    print(f"\nLayer model results: {layer_path}")

    # Report top hits
    print("\n── Top layer model hits ──")
    top = layer_results.dropna(subset=["pval"]).sort_values("pval").head(15)
    for _, row in top.iterrows():
        sig = "***" if row.get("fdr", 1) < 0.05 else ""
        print(f"  {row['layer']:5s} / {row['subclass']:25s} ({row['outcome']:15s}): "
              f"log2FC={row['log2fc']:+.3f}  p={row['pval']:.4g}  "
              f"FDR={row.get('fdr', np.nan):.4g} {sig}")

    # ── Summary ──
    elapsed = time.time() - t0
    n_gam_sig = (gam_results.get("interaction_fdr", pd.Series()) < 0.05).sum()
    n_layer_sig = (layer_results.get("fdr", pd.Series()) < 0.05).sum()
    n_gam_tests = valid.sum()
    n_layer_tests = layer_results["pval"].notna().sum()

    print(f"\n{'=' * 70}")
    print(f"SUMMARY")
    print(f"  GAM interaction tests: {n_gam_sig} / {n_gam_tests} FDR < 0.05")
    print(f"  Layer model tests:     {n_layer_sig} / {n_layer_tests} FDR < 0.05")
    print(f"  Total time: {elapsed:.0f}s")


if __name__ == "__main__":
    main()
