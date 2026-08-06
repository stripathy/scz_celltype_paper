#!/usr/bin/env python3
"""
Add the Lieber/Kwon (spatialDLPFC_SCZ_XENIUM) layer + cell-type annotations onto
our cells, joined by raw cell position (sample, 1-indexed position), and compare
their coarse cell types against our SEA-AD subclasses.

Adds 4 obs columns to every per-sample h5ad AND the merged all_samples h5ad:
  lieber_spd, lieber_layer            (spatial-domain layer; spd07..spd04 -> L1/M..WM)
  lieber_cluster, lieber_celltype     (cell-type cluster -> coarse type)

STALE h5ads (as of 2026-07-27): the `lieber_celltype` column currently sitting in the
h5ads was written with the pre-correction CLUST2CT, so its MGE/CGE labels are transposed
(see the CLUST2CT comment below). subclass_vs_celltype.csv and Supp Fig S6 have been
corrected out-of-band; the h5ad columns have NOT. Rerun this script to bring them in line.
Nothing currently reads `lieber_celltype`, so this is latent, not active, breakage.

Outputs (output/depth_validation/lieber_layers/):
  subclass_vs_celltype.csv   confusion our corr_subclass x their cell type
  celltype_comparison.png    heatmap + class-level summary
"""
import os, sys, argparse
import numpy as np, pandas as pd
import anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__)); _SP = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_SP, "code", "analysis"))
from config import H5AD_DIR, EXCLUDE_SAMPLES, SUBCLASS_TO_CLASS
LAB = os.path.join(_SP, "output/depth_validation/lieber_layers")
MERGED = os.path.join(os.path.dirname(H5AD_DIR), "all_samples_annotated.h5ad")

SPD2LAYER = {"spd07": "L1/M", "spd06": "L2/3", "spd02": "L3/4", "spd05": "L5",
             "spd03": "L6", "spd01": "WMtz", "spd04": "WM"}
# Transcribed from the authors' FINAL object builder in LieberInstitute/
# spatialDLPFC_SCZ_XENIUM: code/analysis/07_cell_type_de/
# 04_create_SPE_with_cell_type_info.R (devel), which writes
# processed-data/07_cell_type_de/cleaned_spe_N24_with_cell_type_and_spds.RDS.
# This dict reproduces that case_when block exactly.
#
# Do NOT source these labels from 06_cell_type_clustering/
# 05_make_updated_cell_type_markers_heatmap_with_annotations.R (L47-48): that
# marker-heatmap script has MGE/CGE transposed (12 -> "MGE", 9 -> "CGE") relative
# to both the final object builder and their 04_explore_cell_type_labels.R
# (L53-54: 12 -> "In: VIP+, LAMP5", 9 -> "In: SST+, PVALB+"). We used script 05
# until 2026-07-27, which put the swap into our Supp Fig S6.
#
# Markers in these same cells confirm the final object is right: cluster 9 is
# SST+/PVALB+/VIP- (59%/47%/4% of cells) and is 91% our Sst+Pvalb+Chandelier =>
# MGE; cluster 12 is VIP+/CXCL14+/CALB2+ (53%/89%/70%) and is 98% our
# Vip+Lamp5+Pax6+Sncg => CGE. Independently confirmed by S.H. Kwon (2026-07-27).
CLUST2CT = {1: "Oligo", 6: "Oligo", 8: "Oligo", 10: "Ambig/Oligo", 7: "Mic",
            15: "Ambig/In/Endo", 18: "L5 Ex", 14: "L6 Ex", 11: "L4/5 Ex", 2: "L2/3 Ex",
            4: "Ast", 16: "Ast", 17: "Ast", 3: "Endo", 5: "Endo", 13: "Endo",
            9: "MGE", 12: "CGE"}

# ---- their per-cell label series (indexed by "{sample}_{rawpos}") ----
lay = pd.read_csv(os.path.join(LAB, "label_transfer_N24_k50_smoothed_labels.csv"))
lay.columns = ["cellid", "spd"]; spd_ser = lay.set_index("cellid")["spd"]
clu = pd.read_csv(os.path.join(LAB, "banksy_clustering_lambda0.1_res0.7.csv"))
clu.columns = ["row", "cluster", "cellid"]; clu_ser = clu.set_index("cellid")["cluster"]


def add_cols(obs, cellids):
    cid = pd.Series(cellids, index=obs.index)
    spd = cid.map(spd_ser)
    obs["lieber_spd"] = spd.fillna("unmatched").values
    obs["lieber_layer"] = spd.map(SPD2LAYER).fillna("unmatched").values
    cl = cid.map(clu_ser)
    obs["lieber_cluster"] = cl.astype("Int64").astype(str).replace("<NA>", "unmatched").values
    obs["lieber_celltype"] = np.where(cl.isna().values, "unmatched",
                                      cl.map(CLUST2CT).fillna("NA").values)
    return obs


def _audit_object(o, label):
    """Report whether depth / domain columns in an obs frame are current."""
    rep = {"object": label, "n_cells": len(o)}
    for c in ("predicted_norm_depth", "spatial_domain", "layer",
              "qc_pass", "corr_qc_pass"):
        rep[f"has_{c}"] = c in o.columns
    if {"predicted_norm_depth", "predicted_norm_depth_prev"} <= set(o.columns):
        d = o["predicted_norm_depth"].astype(float)
        dp = o["predicted_norm_depth_prev"].astype(float)
        rep["depth_is_current"] = not np.allclose(d, dp, equal_nan=True)
    if {"predicted_norm_depth", "predicted_norm_depth_raw"} <= set(o.columns):
        d = o["predicted_norm_depth"].astype(float)
        dr = o["predicted_norm_depth_raw"].astype(float)
        rep["depth_is_smoothed"] = not np.allclose(d, dr, equal_nan=True)
    if "spatial_domain" in o.columns:
        vc = o["spatial_domain"].astype(str).value_counts()
        rep["n_cortical"] = int(vc.get("Cortical", 0))
        if {"qc_pass", "corr_qc_pass"} <= set(o.columns):
            keep = o.qc_pass.astype(bool) & o.corr_qc_pass.astype(bool)
            rep["n_cortical_qc"] = int((keep & (o.spatial_domain.astype(str) == "Cortical")).sum())
    return rep


def stage():
    """Compute the corrected lieber_* columns and audit every object WITHOUT
    writing anything. Returns (staged, audits, diffs, comparison_df)."""
    staged, audits, diffs, coll = {}, [], [], []

    def _stage_one(path, cellids_fn, label):
        a = ad.read_h5ad(path, backed="r")
        o = a.obs.copy()
        a.file.close()
        audits.append(_audit_object(o, label))
        new = add_cols(o.copy(), cellids_fn(o))
        old_ct = o["lieber_celltype"].astype(str) if "lieber_celltype" in o.columns else None
        if old_ct is not None:
            chg = old_ct.values != new["lieber_celltype"].astype(str).values
            if chg.any():
                d = (pd.DataFrame({"old": old_ct.values[chg],
                                   "new": new["lieber_celltype"].astype(str).values[chg]})
                       .value_counts().rename("n").reset_index())
                d.insert(0, "object", label)
                diffs.append(d)
        staged[path] = new[["lieber_spd", "lieber_layer", "lieber_cluster",
                            "lieber_celltype"]]
        coll.append(new[["corr_subclass", "lieber_celltype", "lieber_cluster",
                         "lieber_layer", "spatial_domain"]].copy())
        return new

    for f in sorted(os.listdir(H5AD_DIR)):
        if not f.endswith("_annotated.h5ad"):
            continue
        sid = f.replace("_annotated.h5ad", "")
        if sid in EXCLUDE_SAMPLES:
            continue
        _stage_one(os.path.join(H5AD_DIR, f),
                   lambda o, sid=sid: [f"{sid}_{i+1}" for i in range(len(o))], sid)

    def _merged_ids(o):
        pos = o.groupby("sample_id", sort=False).cumcount()
        return (o["sample_id"].astype(str) + "_" + (pos + 1).astype(str)).values
    merged_obs = _stage_one(MERGED, _merged_ids, "MERGED")
    # the per-sample frames feed the confusion matrix (merged would double-count)
    df = pd.concat(coll[:-1], ignore_index=True)
    return staged, audits, diffs, df, merged_obs


def verify(audits, diffs, df):
    """Gate the write. Returns (ok, messages)."""
    msgs, ok = [], True

    # 1. every object must carry current depth + domain columns
    for a in audits:
        missing = [k[4:] for k in a if k.startswith("has_") and not a[k]]
        if missing:
            ok = False; msgs.append(f"FAIL  {a['object']}: missing {missing}")
        if a.get("depth_is_current") is False:
            ok = False; msgs.append(f"FAIL  {a['object']}: predicted_norm_depth == _prev "
                                    "(stale K=50 depths)")
        if a.get("depth_is_smoothed") is False:
            ok = False; msgs.append(f"FAIL  {a['object']}: depth == raw (smoothing not applied)")

    # 2. MGE/CGE must land on the right interneuron subclasses
    d = df[(df.lieber_celltype != "unmatched") & (df.corr_subclass.astype(str) != "nan")]
    ct = pd.crosstab(d.corr_subclass.astype(str), d.lieber_celltype.astype(str))
    for sub, want in [("Sst", "MGE"), ("Pvalb", "MGE"), ("Vip", "CGE"), ("Lamp5", "CGE")]:
        if sub in ct.index and {"MGE", "CGE"} <= set(ct.columns):
            got = ct.loc[sub, ["MGE", "CGE"]].idxmax()
            if got != want:
                ok = False; msgs.append(f"FAIL  {sub} -> {got}, expected {want} (MGE/CGE swap!)")
            else:
                msgs.append(f"ok    {sub} -> {want} "
                            f"({int(ct.loc[sub, want]):,} cells)")

    # 3. regenerated confusion matrix must match the corrected file on disk
    ref_p = os.path.join(LAB, "subclass_vs_celltype.csv")
    if os.path.exists(ref_p):
        ref = pd.read_csv(ref_p, index_col=0)
        common_i = [i for i in ref.index if i in ct.index]
        common_c = [c for c in ref.columns if c in ct.columns]
        same = np.allclose(ref.loc[common_i, common_c].values,
                           ct.loc[common_i, common_c].values)
        msgs.append(("ok    " if same else "WARN  ") +
                    "regenerated confusion matrix "
                    + ("matches" if same else "DIFFERS from")
                    + " subclass_vs_celltype.csv on disk")
    return ok, msgs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true",
                    help="write the staged columns into the h5ads "
                         "(default: stage and report only)")
    args = ap.parse_args()

    print("STAGING (no files written yet)...", flush=True)
    staged, audits, diffs, df, merged_obs = stage()

    print("\n--- object audit -------------------------------------------------")
    print(pd.DataFrame(audits).to_string(index=False))

    print("\n--- lieber_celltype changes that would be written ----------------")
    if diffs:
        dd = pd.concat(diffs, ignore_index=True)
        tot = dd.groupby(["old", "new"], as_index=False)["n"].sum().sort_values("n", ascending=False)
        print(tot.to_string(index=False))
        print(f"  total cells changed: {int(tot.n.sum()):,}")
    else:
        print("  none - h5ads already match the corrected mapping")

    print("\n--- verification -------------------------------------------------")
    ok, msgs = verify(audits, diffs, df)
    for m in msgs:
        print(" ", m)

    if not args.apply:
        print("\nDRY RUN. Nothing written. Re-run with --apply to write.")
    elif not ok:
        print("\nREFUSING TO WRITE: verification failed above.")
        sys.exit(1)
    else:
        print("\nAPPLYING to h5ads...", flush=True)
        for path, cols in staged.items():
            a = ad.read_h5ad(path)
            for c in cols.columns:
                a.obs[c] = cols[c].values
            a.write_h5ad(path)
            print(f"  wrote {os.path.basename(path)}", flush=True)

    _write_comparison(df)


def _write_comparison(df):
    """Confusion matrix + heatmap of our subclass vs their cell type."""
    df = df[(df.lieber_celltype != "unmatched") & (df.corr_subclass.astype(str) != "nan")].copy()
    df["our_class"] = df.corr_subclass.map(SUBCLASS_TO_CLASS).fillna("?")
    CT_ORDER = ["L2/3 Ex", "L4/5 Ex", "L5 Ex", "L6 Ex", "MGE", "CGE", "Ast", "Oligo",
                "Endo", "Mic", "Ambig/Oligo", "Ambig/In/Endo", "NA"]
    order = sorted(df.corr_subclass.unique(),
                   key=lambda s: ({"Glutamatergic": 0, "GABAergic": 1,
                                   "Non-neuronal": 2}.get(SUBCLASS_TO_CLASS.get(s, ""), 3), str(s)))
    ct = pd.crosstab(df.corr_subclass, df.lieber_celltype).reindex(index=order, columns=CT_ORDER).fillna(0)
    ct.to_csv(os.path.join(LAB, "subclass_vs_celltype.csv"))
    pct = ct.div(ct.sum(1), axis=0) * 100

    fig, ax = plt.subplots(figsize=(11, 11))
    im = ax.imshow(pct.values, cmap="magma", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(CT_ORDER))); ax.set_xticklabels(CT_ORDER, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=9)
    ax.set_xlabel("Lieber/Kwon cell type", fontsize=12); ax.set_ylabel("Our SEA-AD corr_subclass", fontsize=12)
    for i in range(len(order)):
        for j in range(len(CT_ORDER)):
            v = pct.values[i, j]
            if v >= 10:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                        color="white" if v < 55 else "black", fontsize=7)
    ax.set_title("Our subclass vs Lieber cell type (row %)\n", fontsize=13)
    fig.colorbar(im, fraction=0.035, label="% of our subclass")
    plt.tight_layout()
    fig.savefig(os.path.join(LAB, "celltype_comparison.png"), dpi=150, bbox_inches="tight")
    print("\nsaved celltype_comparison.png + subclass_vs_celltype.csv")
    inh = df[df.our_class == "GABAergic"]
    print("\nInterneuron subclass -> Lieber type (resolves MGE/CGE labeling):")
    x = pd.crosstab(inh.corr_subclass.astype(str), inh.lieber_celltype.astype(str))
    print(x.loc[(x[["MGE", "CGE"]].sum(1) > 0) if {"MGE", "CGE"} <= set(x.columns) else slice(None),
                [c for c in ["MGE", "CGE"] if c in x.columns]].to_string())


if __name__ == "__main__":
    main()
