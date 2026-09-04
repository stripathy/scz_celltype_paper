#!/usr/bin/env python3
"""
Depth-model validation restricted to LOW-pathology SEA-AD donors.

Rationale: the Xenium SCZ cohort is not an AD cohort, so comparing against
high-pathology SEA-AD donors conflates platform/method consistency with AD
neurodegeneration. We restrict to the bottom-50%-CPS donors (Continuous
Pseudo-progression Score < 0.5, the lower half of the disease axis) for both
model building and the cross-platform comparisons.

This script answers the gating question -- does the depth model TOLERATE dropping
the high-pathology donors? -- via held-out donor CV WITHIN the low-CPS cohort:
every low-CPS donor is predicted by a model trained only on OTHER low-CPS donors.
Compare the out-of-fold R2 to the full-cohort R2 (0.889).

Outputs (output/depth_validation/lowcps/):
  lowcps_donors.csv, lowcps_summary.csv, lowcps_per_donor.csv,
  lowcps_per_subclass.csv, lowcps_predictions_sample.csv
"""
import os, sys, time
import numpy as np
import pandas as pd
import anndata as ad
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score, mean_absolute_error
from scipy.stats import pearsonr

sys.path.insert(0, "code/analysis"); sys.path.insert(0, "code")
from config import MERFISH_PATH, SUBCLASS_TO_CLASS
from modules.depth_model import build_neighborhood_features

K = 50
N_FOLDS = 5
CPS_MAX = 0.5            # bottom-50%-CPS cohort
CPS_COL = "Continuous Pseudo-progression Score"
GBR_KW = dict(n_estimators=300, max_depth=5, learning_rate=0.1,
              subsample=0.8, random_state=42, min_samples_leaf=20)
OUT = "output/depth_validation/lowcps"
os.makedirs(OUT, exist_ok=True)


def main():
    print("Loading MERFISH (backed)...", flush=True); t0 = time.time()
    a = ad.read_h5ad(MERFISH_PATH, backed="r")
    obs = a.obs
    depth = obs["Normalized depth from pia"].values.astype(float)
    cps = obs[CPS_COL].values.astype(float)
    donor = obs["Donor ID"].values.astype(str)
    subclass = obs["Subclass"].values.astype(str)
    section = obs["Section"].values.astype(str)
    coords = np.asarray(a.obsm["X_spatial_raw"])

    # donor-level CPS -> low-CPS cohort
    dcps = pd.Series(cps, index=donor).groupby(level=0).first()
    low_donors = sorted(dcps.index[dcps < CPS_MAX])
    pd.DataFrame({"donor": dcps.index, "CPS": dcps.values,
                  "set": np.where(dcps.values < CPS_MAX, "low-CPS (kept)", "high-CPS (dropped)")}
                 ).sort_values("CPS").to_csv(os.path.join(OUT, "lowcps_donors.csv"), index=False)
    print(f"  low-CPS cohort: {len(low_donors)} donors (CPS<{CPS_MAX}); dropped {len(dcps)-len(low_donors)} high-CPS")

    keep = (~np.isnan(depth)) & np.isin(donor, low_donors)
    print(f"  low-CPS depth-annotated cells: {keep.sum():,}  ({time.time()-t0:.0f}s)", flush=True)

    t1 = time.time()
    X = build_neighborhood_features(coords[keep], subclass[keep], sorted(set(subclass)),
                                    K=K, sections=section[keep])
    y = depth[keep]; g = donor[keep]; sub = subclass[keep]
    print(f"  features {X.shape}  ({time.time()-t1:.0f}s)", flush=True)

    oof = np.full(len(y), np.nan)
    for i, (tr, te) in enumerate(GroupKFold(n_splits=N_FOLDS).split(X, y, groups=g)):
        t2 = time.time()
        m = GradientBoostingRegressor(**GBR_KW).fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
        print(f"  fold {i+1}/{N_FOLDS}: train={len(tr):,} ({len(set(g[tr]))}d) test={len(te):,} "
              f"({len(set(g[te]))}d) R2={r2_score(y[te],oof[te]):.3f} ({time.time()-t2:.0f}s)", flush=True)

    R2 = r2_score(y, oof); MAE = mean_absolute_error(y, oof); R = pearsonr(y, oof)[0]
    sub_mean = pd.Series(y).groupby(pd.Series(sub)).transform("mean").values
    r2_sub = r2_score(y, sub_mean)
    within_r = float(pearsonr(y - sub_mean,
                              oof - pd.Series(oof).groupby(pd.Series(sub)).transform("mean").values)[0])
    print(f"\nLOW-CPS out-of-fold: R2={R2:.4f} MAE={MAE:.4f} r={R:.4f} (n={len(y):,}, {len(low_donors)} donors)")
    print(f"  subclass-only R2={r2_sub:.3f} | neighborhood R2={R2:.3f} | within-subclass r={within_r:.3f}")
    print(f"  (full-cohort reference: R2=0.889)")

    pd.DataFrame([{"metric": "oof_r2", "value": R2}, {"metric": "oof_mae", "value": MAE},
                  {"metric": "oof_pearson_r", "value": R}, {"metric": "n_cells", "value": int(len(y))},
                  {"metric": "n_donors", "value": len(low_donors)}, {"metric": "cps_max", "value": CPS_MAX},
                  {"metric": "subclass_only_r2", "value": round(r2_sub, 3)},
                  {"metric": "within_subclass_r", "value": round(within_r, 3)},
                  {"metric": "full_cohort_r2", "value": 0.889}]
                 ).to_csv(os.path.join(OUT, "lowcps_summary.csv"), index=False)

    rows = []
    for d in low_donors:
        mk = g == d
        rows.append({"donor": d, "CPS": round(float(dcps[d]), 3), "n_cells": int(mk.sum()),
                     "r2": r2_score(y[mk], oof[mk]), "mae": mean_absolute_error(y[mk], oof[mk])})
    pd.DataFrame(rows).sort_values("r2").to_csv(os.path.join(OUT, "lowcps_per_donor.csv"), index=False)

    rows = []
    for s in sorted(set(sub)):
        mk = sub == s
        if mk.sum() < 20:
            continue
        rows.append({"subclass": s, "klass": SUBCLASS_TO_CLASS.get(s, "Unknown"),
                     "manual_median": float(np.median(y[mk])), "pred_median": float(np.median(oof[mk])),
                     "n_cells": int(mk.sum()),
                     "within_r": float(pearsonr(y[mk], oof[mk])[0]) if mk.sum() > 5 else np.nan})
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "lowcps_per_subclass.csv"), index=False)

    df = pd.DataFrame({"manual": y, "pred": oof, "subclass": sub})
    df["klass"] = df["subclass"].map(lambda s: SUBCLASS_TO_CLASS.get(s, "Unknown"))
    df.sample(n=min(80000, len(df)), random_state=1).to_csv(
        os.path.join(OUT, "lowcps_predictions_sample.csv"), index=False)
    print(f"\nWrote low-CPS outputs to {OUT}", flush=True)


if __name__ == "__main__":
    main()
