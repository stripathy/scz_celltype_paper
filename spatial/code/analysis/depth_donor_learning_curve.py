#!/usr/bin/env python3
"""
Donor learning curve for the depth model (low-CPS cohort).

Answers: how many SEA-AD donors are needed to train a good depth model, and is
a small held-out set (3-4 donors) a stable test? For each training-set size k,
randomly pick k of the 13 low-CPS donors to train on and test on the remaining
(13-k); repeat over R seeds. Plot test R2 (mean +/- sd) vs k. The k=9 (test=4)
and k=10 (test=3) points are the user's proposed train/reserve setup.

Output: output/depth_validation/lowcps/learning_curve.csv
"""
import os, sys, time
import numpy as np
import pandas as pd
import anndata as ad
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import r2_score, mean_absolute_error

_HERE = os.path.dirname(os.path.abspath(__file__))            # .../spatial/code/analysis
_SPATIAL = os.path.dirname(os.path.dirname(_HERE))            # .../spatial
sys.path.insert(0, os.path.join(_SPATIAL, "code", "analysis"))
sys.path.insert(0, os.path.join(_SPATIAL, "code"))
from config import MERFISH_PATH
from modules.depth_model import build_neighborhood_features

K = 50
CPS_MAX = 0.5
KS = list(range(1, 11))      # train-set sizes 1..10
R = 5                         # random splits per size
GBR_KW = dict(n_estimators=300, max_depth=5, learning_rate=0.1,
              subsample=0.8, random_state=42, min_samples_leaf=20)
OUT = os.path.join(_SPATIAL, "output/depth_validation/lowcps")


def main():
    a = ad.read_h5ad(MERFISH_PATH, backed="r")
    obs = a.obs
    depth = obs["Normalized depth from pia"].values.astype(float)
    cps = obs["Continuous Pseudo-progression Score"].values.astype(float)
    donor = obs["Donor ID"].values.astype(str)
    subclass = obs["Subclass"].values.astype(str)
    section = obs["Section"].values.astype(str)
    coords = np.asarray(a.obsm["X_spatial_raw"])
    dcps = pd.Series(cps, index=donor).groupby(level=0).first()
    low = np.array(sorted(dcps.index[dcps < CPS_MAX]))
    keep = (~np.isnan(depth)) & np.isin(donor, low)
    print(f"low-CPS donors={len(low)}  cells={keep.sum():,}", flush=True)

    X = build_neighborhood_features(coords[keep], subclass[keep], sorted(set(subclass)),
                                    K=K, sections=section[keep])
    y = depth[keep]; g = donor[keep]
    print(f"features {X.shape}", flush=True)

    rows = []
    for k in KS:
        for seed in range(R):
            if k >= len(low):
                continue
            rng = np.random.default_rng(1000 * k + seed)
            tr_d = rng.choice(low, size=k, replace=False)
            te_d = np.setdiff1d(low, tr_d)
            tr = np.isin(g, tr_d); te = np.isin(g, te_d)
            t0 = time.time()
            m = GradientBoostingRegressor(**GBR_KW).fit(X[tr], y[tr])
            pr = m.predict(X[te])
            r2 = r2_score(y[te], pr); mae = mean_absolute_error(y[te], pr)
            rows.append({"k_train": k, "seed": seed, "n_test_donors": len(te_d),
                         "n_train_cells": int(tr.sum()), "n_test_cells": int(te.sum()),
                         "r2": r2, "mae": mae})
            print(f"  k={k} seed={seed}: train {len(tr_d)}d/{tr.sum():,}c -> "
                  f"test {len(te_d)}d R2={r2:.3f} ({time.time()-t0:.0f}s)", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "learning_curve.csv"), index=False)
    print("\nMean test R2 by train-set size:")
    print(df.groupby("k_train").agg(r2_mean=("r2", "mean"), r2_sd=("r2", "std"),
                                    mae_mean=("mae", "mean")).round(3).to_string())
    print(f"\nWrote {OUT}/learning_curve.csv", flush=True)


if __name__ == "__main__":
    main()
