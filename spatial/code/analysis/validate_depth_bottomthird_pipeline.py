#!/usr/bin/env python3
"""
Bottom-third-CPS depth pipeline (final validation design).

1. Bottom third of SEA-AD donors by CPS (9 lowest-pathology) -> split 6 train / 3 held-out.
2. Train depth model on the 6 training donors (faithful K=50 / GBR recipe).
3. Apply it to predict depth for ALL Xenium cortical cells (does NOT overwrite the
   deployed predicted_norm_depth -- writes a separate column for the comparison).
4. Compare Xenium against the 3 HELD-OUT MERFISH donors:
   - proportions (subclass + supertype), depth medians (subclass + supertype),
   - and the within-MERFISH held-out per-cell depth R2 (the clean depth validation).

Outputs:
  output/depth_validation/bottomthird/   depth-model held-out validation (+ cv_*-named for the figure)
  output/celltyping_supplement/data_bottomthird/   Xenium-vs-heldout comparison tables
  output/depth_validation/bottomthird/depth_model_bottomthird.pkl   the 6-donor model
"""
import os, sys, time, pickle
import numpy as np
import pandas as pd
import anndata as ad
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import r2_score, mean_absolute_error
from scipy.stats import pearsonr, spearmanr

_HERE = os.path.dirname(os.path.abspath(__file__))
_SPATIAL = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_SPATIAL, "code", "analysis"))
sys.path.insert(0, os.path.join(_SPATIAL, "code"))
from config import (MERFISH_PATH, H5AD_DIR, EXCLUDE_SAMPLES, CORTICAL_LAYERS,
                    SUBCLASS_TO_CLASS, load_cells, load_merfish_cortical)
from modules.depth_model import build_neighborhood_features, smooth_depth_spatial

K = 100               # deployed neighborhood size (K=100 > K=50 on held-out donors)
N_BOTTOM = 9          # bottom third of 27 donors
N_TEST = 3            # held-out donors (reserve), 6 train
MIN_CELLS = 20
GBR_KW = dict(n_estimators=300, max_depth=5, learning_rate=0.1,
              subsample=0.8, random_state=42, min_samples_leaf=20)
DV = os.path.join(_SPATIAL, "output/depth_validation/bottomthird")
CMP = os.path.join(_SPATIAL, "output/celltyping_supplement/data_bottomthird")
os.makedirs(DV, exist_ok=True); os.makedirs(CMP, exist_ok=True)


def mean_prop(df, donor_col, label_col):
    recs = []
    for d, grp in df.groupby(donor_col, observed=True):
        t = len(grp)
        if t:
            for ct, n in grp[label_col].value_counts().items():
                recs.append({"celltype": str(ct), "prop": n / t})
    return pd.DataFrame(recs).groupby("celltype")["prop"].mean()


def corr_block(m, xc, yc, log=True):
    x = m[xc].values.astype(float); y = m[yc].values.astype(float)
    out = {"n": int(len(m)), "pearson_raw": float(pearsonr(x, y)[0]), "spearman": float(spearmanr(x, y)[0])}
    if log:
        pos = (x > 0) & (y > 0)
        out["pearson_log10"] = float(pearsonr(np.log10(x[pos]), np.log10(y[pos]))[0]) if pos.sum() > 3 else np.nan
    return out


def main():
    print("Loading MERFISH obs+coords...", flush=True); t0 = time.time()
    a = ad.read_h5ad(MERFISH_PATH, backed="r")
    obs = a.obs
    depth = obs["Normalized depth from pia"].values.astype(float)
    cps = obs["Continuous Pseudo-progression Score"].values.astype(float)
    donor = obs["Donor ID"].values.astype(str)
    subc = obs["Subclass"].values.astype(str)
    supt = obs["Supertype"].values.astype(str)
    layer = obs["Layer annotation"].astype(str).values
    section = obs["Section"].values.astype(str)
    coords = np.asarray(a.obsm["X_spatial_raw"])
    sub_names = sorted(set(subc))

    # ---- bottom third by CPS, split train/test ----
    dcps = pd.Series(cps, index=donor).groupby(level=0).first().sort_values()
    bottom = list(dcps.index[:N_BOTTOM])
    test_donors = [bottom[i] for i in (1, 4, 7)]           # span CPS within the bottom third
    train_donors = [d for d in bottom if d not in test_donors]
    pd.DataFrame({"donor": bottom, "CPS": [round(float(dcps[d]), 3) for d in bottom],
                  "set": ["test" if d in test_donors else "train" for d in bottom]}
                 ).to_csv(os.path.join(DV, "bottomthird_donors.csv"), index=False)
    print(f"  bottom third (CPS<= {dcps[bottom[-1]]:.3f}): {N_BOTTOM} donors")
    print(f"  TRAIN ({len(train_donors)}): {train_donors}")
    print(f"  HELD-OUT TEST ({len(test_donors)}): {test_donors}  CPS={[round(float(dcps[d]),3) for d in test_donors]}")

    # ---- features on bottom-third depth-annotated cells; train on 6, test on 3 ----
    keep = (~np.isnan(depth)) & np.isin(donor, bottom)
    print(f"  building features for {keep.sum():,} bottom-third cells...", flush=True); t1 = time.time()
    X = build_neighborhood_features(coords[keep], subc[keep], sub_names, K=K, sections=section[keep])
    yk = depth[keep]; gk = donor[keep]; sk = subc[keep]
    print(f"  features {X.shape} ({time.time()-t1:.0f}s)", flush=True)
    tr = np.isin(gk, train_donors); te = np.isin(gk, test_donors)
    model = GradientBoostingRegressor(**GBR_KW).fit(X[tr], yk[tr])
    pred = model.predict(X[te])
    # match deployed model: post-hoc k=30 spatial smoothing, per section
    pred = smooth_depth_spatial(coords[keep][te], pred, k=30, sections=section[keep][te])
    yt, gt, st = yk[te], gk[te], sk[te]
    R2 = r2_score(yt, pred); MAE = mean_absolute_error(yt, pred); R = pearsonr(yt, pred)[0]
    sub_mean = pd.Series(yt).groupby(pd.Series(st)).transform("mean").values
    r2_sub = r2_score(yt, sub_mean)
    within_r = float(pearsonr(yt - sub_mean, pred - pd.Series(pred).groupby(pd.Series(st)).transform("mean").values)[0])
    print(f"\nHELD-OUT depth (3 MERFISH donors): R2={R2:.4f} MAE={MAE:.4f} r={R:.4f} "
          f"| subclass-only R2={r2_sub:.3f} within-r={within_r:.3f}", flush=True)

    pickle.dump({"model": model, "subclass_names": sub_names, "K": K,
                 "train_donors": train_donors, "test_donors": test_donors},
                open(os.path.join(DV, "depth_model_bottomthird.pkl"), "wb"))

    # cv_*-named outputs for plot_depth_validation_supplement.R
    pd.DataFrame([{"metric": "oof_r2", "value": R2}, {"metric": "oof_mae", "value": MAE},
                  {"metric": "oof_pearson_r", "value": R}, {"metric": "n_cells", "value": int(te.sum())},
                  {"metric": "n_donors", "value": len(test_donors)},
                  {"metric": "subclass_only_r2", "value": round(r2_sub, 3)},
                  {"metric": "within_subclass_r", "value": round(within_r, 3)}]
                 ).to_csv(os.path.join(DV, "cv_summary.csv"), index=False)
    pd.DataFrame([{"donor": d, "n_cells": int((gt == d).sum()), "r2": r2_score(yt[gt == d], pred[gt == d]),
                   "mae": mean_absolute_error(yt[gt == d], pred[gt == d])} for d in test_donors]
                 ).sort_values("r2").to_csv(os.path.join(DV, "cv_per_donor.csv"), index=False)
    pd.DataFrame([{"model": "Subclass identity only", "r2": round(r2_sub, 3)},
                  {"model": "Neighborhood model (held-out)", "r2": round(R2, 3)}]
                 ).to_csv(os.path.join(DV, "cv_decomposition.csv"), index=False)
    wrows = []
    for s in sorted(set(st)):
        mk = st == s
        if mk.sum() >= MIN_CELLS:
            wrows.append({"subclass": s, "klass": SUBCLASS_TO_CLASS.get(s, "Unknown"),
                          "within_r": float(pearsonr(yt[mk], pred[mk])[0]), "n": int(mk.sum())})
    pd.DataFrame(wrows).to_csv(os.path.join(DV, "cv_within_subclass.csv"), index=False)
    samp = pd.DataFrame({"manual": yt, "pred": pred, "subclass": st})
    samp["klass"] = samp["subclass"].map(lambda s: SUBCLASS_TO_CLASS.get(s, "Unknown"))
    samp.sample(n=min(80000, len(samp)), random_state=1).to_csv(os.path.join(DV, "cv_predictions_sample.csv"), index=False)

    # ---- apply model to Xenium (per sample; does NOT overwrite deployed depth) ----
    print("\nPredicting Xenium depth with the 6-donor model...", flush=True)
    xparts = []
    for sid in sorted(f.replace("_annotated.h5ad", "") for f in os.listdir(H5AD_DIR) if f.endswith("_annotated.h5ad")):
        if sid in EXCLUDE_SAMPLES:
            continue
        d = load_cells(sid, cortical_only=True, qc_mode="corr")
        feats = build_neighborhood_features(d[["x", "y"]].values, d["subclass_label"].astype(str).values,
                                            sub_names, K=K, sections=None)
        d = d[["sample_id", "subclass_label", "supertype_label"]].copy()
        d["xdepth"] = smooth_depth_spatial(d[["x", "y"]].values, model.predict(feats), k=30)
        xparts.append(d)
        print(f"  {sid}: {len(d):,} cells", flush=True)
    xen = pd.concat(xparts, ignore_index=True)
    print(f"  Xenium cortical cells with new depth: {len(xen):,}")

    # ---- held-out MERFISH reference (the 3 test donors, cortical) ----
    mref = load_merfish_cortical()
    mref = mref[mref["donor"].isin(test_donors)].copy()
    sup2sub = mref[["supertype", "subclass"]].drop_duplicates().set_index("supertype")["subclass"].to_dict()
    klass_of = lambda ct, lvl: SUBCLASS_TO_CLASS.get(ct if lvl == "subclass" else sup2sub.get(ct), "Unknown")
    print(f"  held-out MERFISH reference: {len(mref):,} cells, {mref['donor'].nunique()} donors")

    # ---- comparison: Xenium vs held-out-3 MERFISH ----
    stats = []
    for lvl, (mlab, xlab) in {"subclass": ("subclass", "subclass_label"),
                              "supertype": ("supertype", "supertype_label")}.items():
        mp = mean_prop(mref, "donor", mlab).rename("merfish_prop")
        xp = mean_prop(xen, "sample_id", xlab).rename("xenium_prop")
        mn = mref[mlab].value_counts().rename("n_merfish_cells"); xn = xen[xlab].value_counts().rename("n_xenium_cells")
        prop = pd.concat([mp, xp, mn, xn], axis=1).dropna(subset=["merfish_prop", "xenium_prop"])
        prop = prop[(prop.n_merfish_cells >= MIN_CELLS) & (prop.n_xenium_cells >= MIN_CELLS)].reset_index().rename(columns={"index": "celltype"})
        prop["klass"] = [klass_of(c, lvl) for c in prop.celltype]
        prop.sort_values("merfish_prop", ascending=False).to_csv(os.path.join(CMP, f"prop_{lvl}.csv"), index=False)
        cb = corr_block(prop, "merfish_prop", "xenium_prop"); stats.append({"panel": f"proportion_{lvl}", **cb})
        print(f"[prop {lvl}] n={cb['n']} raw={cb['pearson_raw']:.3f} log10={cb['pearson_log10']:.3f} Spearman={cb['spearman']:.3f}")

        md = mref.groupby(mlab, observed=True)["depth"].median().rename("merfish_depth")
        xd = xen.groupby(xlab, observed=True)["xdepth"].median().rename("xenium_depth")
        dep = pd.concat([md, mn, xd, xn], axis=1).dropna(subset=["merfish_depth", "xenium_depth"])
        dep = dep[(dep.n_merfish_cells >= MIN_CELLS) & (dep.n_xenium_cells >= MIN_CELLS)].reset_index().rename(columns={"index": "celltype"})
        dep["klass"] = [klass_of(c, lvl) for c in dep.celltype]
        dep.sort_values("merfish_depth").to_csv(os.path.join(CMP, f"depth_{lvl}.csv"), index=False)
        cb = corr_block(dep, "merfish_depth", "xenium_depth", log=False); stats.append({"panel": f"depth_{lvl}", **cb})
        print(f"[depth {lvl}] n={cb['n']} raw={cb['pearson_raw']:.3f} Spearman={cb['spearman']:.3f}")

    pd.DataFrame(stats).to_csv(os.path.join(CMP, "concordance_stats.csv"), index=False)
    print(f"\nDone. depth-model: {DV} | comparison: {CMP}  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
