#!/usr/bin/env python3
"""
One-time migration: rename cluster_label → supertype_label in all Xenium h5ad files.

The MERFISH reference uses "Supertype" for the finest-level cell type annotation,
but the label transfer pipeline originally output this as "cluster_label" in the
Xenium h5ad files. This script aligns the naming convention.

Renames:
  cluster_label            → supertype_label
  cluster_label_confidence → supertype_label_confidence
"""

import os
import glob
import anndata as ad

H5AD_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "output", "h5ad")


def main():
    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))
    print(f"Found {len(h5ad_files)} h5ad files in {H5AD_DIR}")

    for fpath in h5ad_files:
        sample_id = os.path.basename(fpath).replace("_annotated.h5ad", "")
        adata = ad.read_h5ad(fpath)

        renamed = False

        if "cluster_label" in adata.obs.columns:
            adata.obs = adata.obs.rename(columns={
                "cluster_label": "supertype_label"
            })
            renamed = True

        if "cluster_label_confidence" in adata.obs.columns:
            adata.obs = adata.obs.rename(columns={
                "cluster_label_confidence": "supertype_label_confidence"
            })
            renamed = True

        if renamed:
            adata.write_h5ad(fpath)
            print(f"  {sample_id}: renamed cluster_label → supertype_label ✓")
        else:
            if "supertype_label" in adata.obs.columns:
                print(f"  {sample_id}: already has supertype_label, skipping")
            else:
                print(f"  {sample_id}: WARNING — no cluster_label or supertype_label found!")

    print("\nDone.")


if __name__ == "__main__":
    main()
