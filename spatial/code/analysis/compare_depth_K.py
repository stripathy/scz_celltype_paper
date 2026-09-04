#!/usr/bin/env python3
"""
Retest K=50 vs K=100 neighborhood size for the depth model, in the current
low-CPS held-out-CV framework (GroupKFold by donor within the 13 low-CPS donors).
K changes only how many spatial neighbors are averaged into the 48-D subclass-
fraction features -- same feature dimension, same GBR cost -- so this is cheap.
"""
import os, sys, time
import numpy as np
import pandas as pd
import anndata as ad
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score, mean_absolute_error
from scipy.stats import pearsonr

_HERE = os.path.dirname(os.path.abspath(__file__))
_SPATIAL = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_SPATIAL, "code", "analysis"))
sys.path.insert(0, os.path.join(_SPATIAL, "code"))
from config import MERFISH_PATH
from modules.depth_model import build_neighborhood_features

N_FOLDS = 5
CPS_MAX = 0.5
GBR_KW = dict(n_estimators=300, max_depth=5, learning_rate=0.1,
              subsample=0.8, random_state=42, min_samples_leaf=20)


def main():
    a = ad.read_h5ad(MERFISH_PATH, backed="r")
    obs = a.obs
    depth = obs["Normalized depth from pia"].values.astype(float)
    cps = obs["Continuous Pseudo-progression Score"].values.astype(float)
    donor = obs["Donor ID"].values.astype(str)
    subc = obs["Subclass"].values.astype(str)
    section = obs["Section"].values.astype(str)
    coords = np.asarray(a.obsm["X_spatial_raw"])
    dcps = pd.Series(cps, index=donor).groupby(level=0).first()
    low = sorted(dcps.index[dcps < CPS_MAX])
    keep = (~np.isnan(depth)) & np.isin(donor, low)
    y = depth[keep]; g = donor[keep]; sub = subc[keep]
    co = coords[keep]; se = section[keep]; names = sorted(set(subc))
    print(f"low-CPS: {len(low)} donors, {keep.sum():,} cells\n")

    rows = []
    for K in [50, 100]:
        t0 = time.time()
        X = build_neighborhood_features(co, sub, names, K=K, sections=se)
        oof = np.full(len(y), np.nan)
        for tr, te in GroupKFold(n_splits=N_FOLDS).split(X, y, groups=g):
            oof[te] = GradientBoostingRegressor(**GBR_KW).fit(X[tr], y[tr]).predict(X[te])
        R2 = r2_score(y, oof); MAE = mean_absolute_error(y, oof); R = pearsonr(y, oof)[0]
        sm = pd.Series(y).groupby(pd.Series(sub)).transform("mean").values
        wr = float(pearsonr(y - sm, oof - pd.Series(oof).groupby(pd.Series(sub)).transform("mean").values)[0])
        rows.append({"K": K, "oof_r2": round(R2, 4), "oof_mae": round(MAE, 4),
                     "oof_r": round(R, 4), "within_subclass_r": round(wr, 3)})
        print(f"K={K:3d}: R2={R2:.4f}  MAE={MAE:.4f}  r={R:.4f}  within-r={wr:.3f}  ({time.time()-t0:.0f}s)", flush=True)
    out = os.path.join(_SPATIAL, "output/depth_validation/lowcps/K_comparison.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
