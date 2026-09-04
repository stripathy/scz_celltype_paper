#!/usr/bin/env python3
"""
Cross-platform composition + depth comparison, MERFISH reference restricted to
the LOW-pathology cohort (CPS < 0.5) -- the SAME donor set used for the depth
model validation. Removes the AD-neurodegeneration confound and keeps one
coherent reference for both proportion and depth comparisons.

Same recipe as build_celltyping_validation_data.py (cortical L1-L6, corr labels,
proportions = mean across donors, depth = pooled-cell median); only the MERFISH
donor set changes (low-CPS only).

Outputs: output/celltyping_supplement/data_lowcps/
  prop_subclass.csv, prop_supertype.csv, depth_subclass.csv, depth_supertype.csv,
  concordance_stats.csv, lowcps_reference_donors.csv
"""
import os, sys
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, "code/analysis")
from config import (H5AD_DIR, EXCLUDE_SAMPLES, MERFISH_PATH, CORTICAL_LAYERS,
                    SUBCLASS_TO_CLASS, load_cells, load_merfish_cortical)

CPS_MAX = 0.5
MIN_CELLS = 20
OUT = "output/celltyping_supplement/data_lowcps"
os.makedirs(OUT, exist_ok=True)


def mean_prop(df, donor_col, label_col):
    recs = []
    for donor, grp in df.groupby(donor_col, observed=True):
        tot = len(grp)
        if tot:
            for ct, n in grp[label_col].value_counts().items():
                recs.append({"celltype": str(ct), "prop": n / tot})
    return pd.DataFrame(recs).groupby("celltype")["prop"].mean()


def corr_block(m, xc, yc):
    x = m[xc].values.astype(float); y = m[yc].values.astype(float)
    out = {"n": int(len(m)), "pearson_raw": float(pearsonr(x, y)[0]),
           "spearman": float(spearmanr(x, y)[0])}
    pos = (x > 0) & (y > 0)
    out["pearson_log10"] = float(pearsonr(np.log10(x[pos]), np.log10(y[pos]))[0]) if pos.sum() > 3 else np.nan
    return out


def main():
    # ---- low-CPS MERFISH reference ----
    a = ad.read_h5ad(MERFISH_PATH, backed="r")
    dcps = (pd.Series(a.obs["Continuous Pseudo-progression Score"].values.astype(float),
                      index=a.obs["Donor ID"].values.astype(str)).groupby(level=0).first())
    low_donors = sorted(dcps.index[dcps < CPS_MAX])
    pd.DataFrame({"donor": low_donors, "CPS": [round(float(dcps[d]), 3) for d in low_donors]}
                ).to_csv(os.path.join(OUT, "lowcps_reference_donors.csv"), index=False)
    print(f"Low-CPS MERFISH reference: {len(low_donors)} donors (CPS<{CPS_MAX})")

    mer = load_merfish_cortical()
    mer = mer[mer["donor"].isin(low_donors)].copy()
    print(f"  {len(mer):,} low-CPS cortical cells, {mer['donor'].nunique()} donors")
    sup2sub = mer[["supertype", "subclass"]].drop_duplicates().set_index("supertype")["subclass"].to_dict()

    def klass_of(ct, level):
        sub = ct if level == "subclass" else sup2sub.get(ct)
        return SUBCLASS_TO_CLASS.get(sub, "Unknown")

    # ---- Xenium (unchanged: all 24 donors, cortical, corr labels) ----
    sample_ids = sorted(f.replace("_annotated.h5ad", "") for f in os.listdir(H5AD_DIR)
                        if f.endswith("_annotated.h5ad"))
    xen = pd.concat([load_cells(s, cortical_only=True, qc_mode="corr",
                                extra_obs_columns=["predicted_norm_depth"])
                     for s in sample_ids if s not in EXCLUDE_SAMPLES], ignore_index=True)
    print(f"  {len(xen):,} Xenium cortical cells, {xen['sample_id'].nunique()} donors")

    stats = []
    for level, (mlab, xlab) in {"subclass": ("subclass", "subclass_label"),
                                "supertype": ("supertype", "supertype_label")}.items():
        # proportions (mean across donors)
        mp = mean_prop(mer, "donor", mlab).rename("merfish_prop")
        xp = mean_prop(xen, "sample_id", xlab).rename("xenium_prop")
        mn = mer[mlab].value_counts().rename("n_merfish_cells")
        xn = xen[xlab].value_counts().rename("n_xenium_cells")
        prop = pd.concat([mp, xp, mn, xn], axis=1).dropna(subset=["merfish_prop", "xenium_prop"])
        prop = prop[(prop.n_merfish_cells >= MIN_CELLS) & (prop.n_xenium_cells >= MIN_CELLS)].reset_index()
        prop = prop.rename(columns={"index": "celltype"})
        prop["klass"] = [klass_of(c, level) for c in prop.celltype]
        prop.sort_values("merfish_prop", ascending=False).to_csv(os.path.join(OUT, f"prop_{level}.csv"), index=False)
        cb = corr_block(prop, "merfish_prop", "xenium_prop"); stats.append({"panel": f"proportion_{level}", **cb})
        print(f"[prop {level}] n={cb['n']} Pearson(raw)={cb['pearson_raw']:.3f} log10={cb['pearson_log10']:.3f} Spearman={cb['spearman']:.3f}")

        # depth (pooled-cell median)
        md = mer.dropna(subset=["depth"]).groupby(mlab, observed=True)["depth"].median().rename("merfish_depth")
        mdn = mer.groupby(mlab, observed=True).size().rename("n_merfish_cells")
        xd = xen.dropna(subset=["predicted_norm_depth"]).groupby(xlab, observed=True)["predicted_norm_depth"].median().rename("xenium_depth")
        xdn = xen.groupby(xlab, observed=True).size().rename("n_xenium_cells")
        dep = pd.concat([md, mdn, xd, xdn], axis=1).dropna(subset=["merfish_depth", "xenium_depth"])
        dep = dep[(dep.n_merfish_cells >= MIN_CELLS) & (dep.n_xenium_cells >= MIN_CELLS)].reset_index()
        dep = dep.rename(columns={"index": "celltype"})
        dep["klass"] = [klass_of(c, level) for c in dep.celltype]
        dep.sort_values("merfish_depth").to_csv(os.path.join(OUT, f"depth_{level}.csv"), index=False)
        cb = corr_block(dep, "merfish_depth", "xenium_depth"); stats.append({"panel": f"depth_{level}", **cb})
        print(f"[depth {level}] n={cb['n']} Pearson(raw)={cb['pearson_raw']:.3f} Spearman={cb['spearman']:.3f}")

    pd.DataFrame(stats).to_csv(os.path.join(OUT, "concordance_stats.csv"), index=False)
    print(f"\nWrote low-CPS comparison tables to {OUT}")


if __name__ == "__main__":
    main()
