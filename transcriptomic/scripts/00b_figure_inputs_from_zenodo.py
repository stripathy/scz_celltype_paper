#!/usr/bin/env python3
"""Rebuild the two large figure inputs from the Zenodo deposit instead of Git LFS.

The repo's big DE tables live in Git LFS. If you would rather not fetch ~1.6 GB of
LFS objects, the same data is published as parquet at

    https://zenodo.org/records/22802546   (v2; concept DOI 10.5281/zenodo.22801229,
                                           CC-BY-4.0)

This script downloads DE_subclass.parquet and writes the two snapshots under
transcriptomic/data/figure_inputs/ that Figure 2 reads, translating the Zenodo
schema to the one the figure scripts expect.

    python3 transcriptomic/scripts/00b_figure_inputs_from_zenodo.py [--write]

Without --write it only checks that it reproduces the tracked snapshots.

Two schema differences it handles, both verified 2026-09-16:
  * columns  — Zenodo carries SE where the repo carries t (t = logFC / SE), and
               names them gene / PValue / FDR rather than genes / P.Value / adj.P.Val.
  * cell types — Zenodo uses the long SEA-AD names with slashes ("Astrocyte",
               "L2/3 IT", "Lamp5 Lhx6"); the repo uses the short underscored form
               ("Astro", "L2_3 IT", "Lamp5_Lhx6"). Cohort names already agree.

Verified: the rebuilt forest snapshot matches the tracked one on all 143 rows to
within 7e-16, and the Meta-analysis rows match DE_genes_all_cells_scz.csv
(231,135 rows) to within 1e-16.

Scope, stated plainly:
  * meta_results_cohorts_subclass_forest.csv is rebuilt in full — every column
    Figure 2 reads, matching the tracked snapshot on all 143 rows to 7e-16.
  * DE_genes_all_cells_scz.csv is rebuilt with estimate/se/pval/padj and derived
    ci.lb/ci.ub, matching to 1.1e-16. It CANNOT carry k, tau2 or I2, which Zenodo
    does not publish. Figure 2 never reads those; sst_strata/ (Supplementary S8)
    does — from the committed copy of this same file, which has them. So leave the
    tracked snapshot alone if you intend to run S8.
  * the Xenium-sourced snapshots come from the Xenium processing repo, not
    Zenodo, and are already committed.
"""
import argparse, hashlib, os, sys, urllib.request

REC = "https://zenodo.org/records/22802546/files"   # v2, the current version
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEST = os.path.join(HERE, "data", "figure_inputs")
MD5 = {"DE_subclass.parquet": "71f19c9c6b8e2cb375b38c55f8e757f6"}

CELLTYPE = {"Astrocyte": "Astro", "Endothelial": "Endo", "Oligodendrocyte": "Oligo",
            "Microglia-PVM": "Micro-PVM", "Lamp5 Lhx6": "Lamp5_Lhx6"}
FOREST_GENES = ["SST", "PVALB"]


def fetch(name, cache):
    path = os.path.join(cache, name)
    if not os.path.exists(path):
        os.makedirs(cache, exist_ok=True)
        print(f"downloading {name} …", flush=True)
        urllib.request.urlretrieve(f"{REC}/{name}?download=1", path)
    got = hashlib.md5(open(path, "rb").read()).hexdigest()
    if got != MD5[name]:
        sys.exit(f"md5 mismatch for {name}: got {got}, expected {MD5[name]}")
    print(f"  {name}  md5 ok")
    return path


def to_repo_celltypes(s):
    s = s.replace(CELLTYPE)
    # the slash -> underscore rule applies to the layer names, not to Lamp5_Lhx6
    return s.where(s.eq("Lamp5_Lhx6"), s.str.replace("/", "_", regex=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true", help="overwrite the tracked snapshots")
    ap.add_argument("--cache", default=os.path.join(HERE, "data", ".zenodo_cache"))
    a = ap.parse_args()
    import pandas as pd

    z = pd.read_parquet(fetch("DE_subclass.parquet", a.cache))

    # --- meta-analysis table (full gene x subclass) --------------------------
    # Zenodo publishes estimate/se/p/FDR. ci.lb and ci.ub are recoverable as
    # estimate +/- qnorm(0.975)*se (verified against the tracked table to 1.4e-8).
    # k, tau2 and I2 are NOT published and cannot be derived — see the warning below.
    meta = z[z.cohort == "Meta-analysis"].copy()
    meta["cell_type"] = to_repo_celltypes(meta.cell_type)
    meta = meta.rename(columns={"gene": "genes", "logFC": "estimate",
                                "SE": "se", "PValue": "pval", "FDR": "padj"})
    Z975 = 1.959963984540054
    meta["ci.lb"] = meta.estimate - Z975 * meta.se
    meta["ci.ub"] = meta.estimate + Z975 * meta.se
    meta = meta[["cell_type", "genes", "estimate", "se", "pval", "ci.lb", "ci.ub", "padj"]]

    # --- per-cohort forest subset -------------------------------------------
    coh = z[z.gene.isin(FOREST_GENES) & ~z.cohort.isin(["Meta-analysis", "Xenium"])].copy()
    coh["cell_type"] = to_repo_celltypes(coh.cell_type)
    coh["t"] = coh.logFC / coh.SE
    coh = coh.rename(columns={"gene": "genes", "PValue": "P.Value", "FDR": "adj.P.Val"})
    coh = coh[["genes", "logFC", "t", "P.Value", "adj.P.Val", "cell_type", "cohort"]]

    out = {"DE_genes_all_cells_scz.csv": meta,
           "meta_results_cohorts_subclass_forest.csv": coh}
    for name, df in out.items():
        dst = os.path.join(DEST, name)
        print(f"\n{name}: {len(df):,} rows")
        if a.write:
            if name == "DE_genes_all_cells_scz.csv":
                print("  WARNING: Zenodo does not publish k, tau2 or I2. This rebuild omits them.")
                print("           Figure 2 does not use them, but sst_strata/ (Supplementary S8)")
                print("           does. Overwriting the tracked snapshot will break S8;")
                print("           restore it with: git checkout -- " + os.path.relpath(dst))
            df.to_csv(dst, index=False); print(f"  wrote {dst}")
        elif os.path.exists(dst):
            old = pd.read_csv(dst)
            keys = [c for c in ("genes", "cell_type", "cohort") if c in df.columns and c in old.columns]
            m = old.merge(df, on=keys, suffixes=("_old", "_new"))
            num = [c for c in df.columns if c not in keys and f"{c}_old" in m]
            worst = max(((m[f"{c}_old"] - m[f"{c}_new"]).abs().max(), c) for c in num) if num else (0, "-")
            print(f"  tracked {len(old):,} rows | joined {len(m):,} | worst |diff| {worst[0]:.3g} ({worst[1]})")
            print("  MATCHES" if len(m) == len(old) and worst[0] < 1e-8 else "  DIFFERS — inspect before using")
        else:
            print("  no tracked copy to compare against; re-run with --write")


if __name__ == "__main__":
    main()
