#!/usr/bin/env python3
"""
Donor-grouped cross-validation of the MERFISH cortical-depth model.

The deployed model (modules/depth_model.train_depth_model) holds out only the
3 SMALLEST donors (~3% of cells) as its test set. This script instead runs
GroupKFold by donor, so EVERY donor is predicted out-of-fold by a model trained
on a disjoint set of donors -- a proper "train on one subset of SEA-AD MERFISH,
test against held-out SEA-AD depths" generalization test.

Faithful to the deployed recipe: identical K=50 neighborhood-composition
features (build_neighborhood_features) and identical GradientBoostingRegressor
hyperparameters. Only the train/test partition changes (GroupKFold vs 3-smallest).

Outputs (output/depth_validation/):
  cv_summary.csv          overall out-of-fold R2 / MAE / Pearson r
  cv_per_donor.csv        per-donor held-out R2 / MAE / r (n=27) -> consistency
  cv_per_subclass.csv     manual vs cross-validated-predicted MEDIAN depth (per subclass)
  cv_per_layer.csv        out-of-fold MAE within each manual cortical layer
  cv_predictions_sample.csv  stratified ~80k-cell sample of (manual, pred) for density panel
"""
import os, sys, time
import numpy as np
import pandas as pd
import anndata as ad
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score, mean_absolute_error
from scipy.stats import pearsonr

HERE = os.path.dirname(os.path.abspath(__file__))
SPATIAL = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))   # .../spatial
sys.path.insert(0, os.path.join(SPATIAL, "code", "analysis"))
sys.path.insert(0, os.path.join(SPATIAL, "code"))
from config import MERFISH_PATH, SUBCLASS_TO_CLASS, CORTICAL_LAYERS
from modules.depth_model import build_neighborhood_features

K = 50
N_FOLDS = 5
GBR_KW = dict(n_estimators=300, max_depth=5, learning_rate=0.1,
              subsample=0.8, random_state=42, min_samples_leaf=20)
OUT = os.path.join(SPATIAL, "output", "depth_validation")
os.makedirs(OUT, exist_ok=True)


def main():
    print("Loading MERFISH (backed)...", flush=True); t0 = time.time()
    a = ad.read_h5ad(MERFISH_PATH, backed="r")
    obs = a.obs
    depth = obs["Normalized depth from pia"].values.astype(float)
    has = ~np.isnan(depth)                      # deployed model's exact cell set
    subclass = obs["Subclass"].values.astype(str)
    donors = obs["Donor ID"].values.astype(str)
    sections = obs["Section"].values.astype(str)
    layer = obs["Layer annotation"].astype(str).values
    coords = np.asarray(a.obsm["X_spatial_raw"])
    print(f"  total={a.n_obs:,}  has_depth={has.sum():,}  donors={len(set(donors[has]))}"
          f"  ({time.time()-t0:.0f}s)", flush=True)

    subclass_names = sorted(set(subclass))
    print(f"Building K={K} neighborhood features for {has.sum():,} cells (once)...", flush=True)
    t1 = time.time()
    X = build_neighborhood_features(coords[has], subclass[has], subclass_names,
                                    K=K, sections=sections[has])
    y = depth[has]; g = donors[has]; sub = subclass[has]; lay = layer[has]
    print(f"  features {X.shape}  ({time.time()-t1:.0f}s)", flush=True)

    # ---- GroupKFold by donor: out-of-fold predictions ----
    oof = np.full(len(y), np.nan)
    gkf = GroupKFold(n_splits=N_FOLDS)
    for i, (tr, te) in enumerate(gkf.split(X, y, groups=g)):
        t2 = time.time()
        m = GradientBoostingRegressor(**GBR_KW)
        m.fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
        print(f"  fold {i+1}/{N_FOLDS}: train={len(tr):,} ({len(set(g[tr]))} donors) "
              f"test={len(te):,} ({len(set(g[te]))} donors) "
              f"R2={r2_score(y[te], oof[te]):.3f} MAE={mean_absolute_error(y[te], oof[te]):.4f} "
              f"({time.time()-t2:.0f}s)", flush=True)

    R2 = r2_score(y, oof); MAE = mean_absolute_error(y, oof); R = pearsonr(y, oof)[0]
    print(f"\nOUT-OF-FOLD overall: R2={R2:.4f}  MAE={MAE:.4f}  r={R:.4f}  n={len(y):,}", flush=True)

    pd.DataFrame([{"metric": "oof_r2", "value": R2}, {"metric": "oof_mae", "value": MAE},
                  {"metric": "oof_pearson_r", "value": R}, {"metric": "n_cells", "value": len(y)},
                  {"metric": "n_folds", "value": N_FOLDS}, {"metric": "K", "value": K},
                  {"metric": "n_donors", "value": len(set(g))}]).to_csv(
                      os.path.join(OUT, "cv_summary.csv"), index=False)

    # ---- per-donor held-out metrics ----
    rows = []
    for d in sorted(set(g)):
        mk = g == d
        rows.append({"donor": d, "n_cells": int(mk.sum()),
                     "r2": r2_score(y[mk], oof[mk]),
                     "mae": mean_absolute_error(y[mk], oof[mk]),
                     "pearson_r": pearsonr(y[mk], oof[mk])[0]})
    pd.DataFrame(rows).sort_values("r2").to_csv(os.path.join(OUT, "cv_per_donor.csv"), index=False)

    # ---- per-subclass manual vs CV-predicted median depth ----
    rows = []
    for s in subclass_names:
        mk = sub == s
        if mk.sum() < 20:
            continue
        rows.append({"subclass": s, "klass": SUBCLASS_TO_CLASS.get(s, "Unknown"),
                     "manual_median": float(np.median(y[mk])),
                     "pred_median": float(np.median(oof[mk])),
                     "n_cells": int(mk.sum()),
                     "mae": mean_absolute_error(y[mk], oof[mk])})
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "cv_per_subclass.csv"), index=False)

    # ---- per-layer out-of-fold error ----
    rows = []
    for L in CORTICAL_LAYERS:
        mk = lay == L
        if mk.sum() < 20:
            continue
        rows.append({"layer": L, "n_cells": int(mk.sum()),
                     "manual_median": float(np.median(y[mk])),
                     "pred_median": float(np.median(oof[mk])),
                     "mae": mean_absolute_error(y[mk], oof[mk])})
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "cv_per_layer.csv"), index=False)

    # ---- stratified per-cell sample for the density scatter ----
    df = pd.DataFrame({"manual": y, "pred": oof, "subclass": sub, "donor": g})
    df["klass"] = df["subclass"].map(lambda s: SUBCLASS_TO_CLASS.get(s, "Unknown"))
    samp = df.sample(n=min(80000, len(df)), random_state=1)
    samp.to_csv(os.path.join(OUT, "cv_predictions_sample.csv"), index=False)

    print(f"\nWrote CV outputs to {OUT}", flush=True)


if __name__ == "__main__":
    main()
