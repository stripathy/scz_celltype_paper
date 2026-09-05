#!/usr/bin/env python3
"""
Independent anchor for Xenium predicted depth: canonical cortical layer markers.

The depth model uses neighborhood SUBCLASS composition -- it never sees these
marker genes' expression as a depth target. So if Xenium predicted depth orders
known laminar markers correctly (upper markers shallow, deep markers deep), that
validates the Xenium depths using a signal independent of the MERFISH reference.

Outputs (output/depth_validation/):
  anchor_correlations.csv  per-marker Spearman r of expression vs predicted depth
  anchor_profiles.csv      mean z-scored expression per depth bin (for gradient plot)
"""
import os, sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sys.path.insert(0, "code/analysis")
from config import H5AD_DIR, EXCLUDE_SAMPLES

OUT = "output/depth_validation"
# canonical layer markers present in the 300-gene panel, ordered pia->WM
MARKERS = [("CUX2", "upper L2/3"), ("LAMP5", "upper L2/3"), ("CARTPT", "upper L2/3"),
           ("CALB1", "upper L2/3"), ("RORB", "L4"), ("PCP4", "L5"),
           ("NR4A2", "L6/L6b"), ("NXPH4", "L6/L6b")]
MGENES = [m for m, _ in MARKERS]
NBINS = 20


def main():
    samples = sorted(f.replace("_annotated.h5ad", "") for f in os.listdir(H5AD_DIR)
                     if f.endswith("_annotated.h5ad"))
    depth_all, expr_all = [], []
    for sid in samples:
        if sid in EXCLUDE_SAMPLES:
            continue
        a = ad.read_h5ad(os.path.join(H5AD_DIR, f"{sid}_annotated.h5ad"))
        ok = (a.obs.get("corr_qc_pass", pd.Series(True, index=a.obs.index)).astype(bool)
              & (a.obs["spatial_domain"].astype(str) == "Cortical")
              & a.obs["predicted_norm_depth"].notna())
        a = a[ok.values].copy()
        if a.n_obs == 0:
            continue
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        present = [g for g in MGENES if g in a.var_names]
        E = a[:, present].X
        E = np.asarray(E.todense()) if hasattr(E, "todense") else np.asarray(E)
        df = pd.DataFrame(E, columns=present)
        df["depth"] = a.obs["predicted_norm_depth"].values
        depth_all.append(df["depth"].values)
        expr_all.append(df)
        print(f"  {sid}: {a.n_obs:,} cortical cells", flush=True)

    big = pd.concat(expr_all, ignore_index=True)
    print(f"\nPooled cortical cells: {len(big):,}")

    # per-marker Spearman r vs predicted depth
    rows = []
    for g, lay in MARKERS:
        if g not in big.columns:
            continue
        r = big[[g, "depth"]].corr(method="spearman").iloc[0, 1]
        rows.append({"marker": g, "layer": lay, "spearman_r": round(float(r), 3),
                     "expected_sign": "-" if lay.startswith("upper") else
                                      ("+" if lay.startswith("L6") else "mid")})
    cor = pd.DataFrame(rows)
    cor.to_csv(os.path.join(OUT, "anchor_correlations.csv"), index=False)
    print("\nMarker vs predicted-depth Spearman r:")
    print(cor.to_string(index=False))

    # mean z-scored expression per depth bin (gradient profiles)
    big["bin"] = pd.cut(big["depth"].clip(0, 1), bins=np.linspace(0, 1, NBINS + 1),
                        labels=False, include_lowest=True)
    prof = []
    for g, lay in MARKERS:
        if g not in big.columns:
            continue
        z = (big[g] - big[g].mean()) / big[g].std()
        m = z.groupby(big["bin"]).mean()
        for b, v in m.items():
            prof.append({"marker": g, "layer": lay, "depth_bin": (b + 0.5) / NBINS,
                         "mean_z": float(v)})
    pd.DataFrame(prof).to_csv(os.path.join(OUT, "anchor_profiles.csv"), index=False)
    print(f"\nWrote anchor_correlations.csv + anchor_profiles.csv to {OUT}")


if __name__ == "__main__":
    main()
