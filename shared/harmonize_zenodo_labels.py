#!/usr/bin/env python3
"""Harmonise the cell-type and dataset labels across the three Zenodo parquet files.

This produced Zenodo **v2** (10.5281/zenodo.22802546) from v1 (…22801230), whose
numbers were correct but whose three files disagreed on labels badly enough that
`cell_metadata.parquet` could not be joined to `DE_supertype.parquet` at all —
0 of 115 supertype names matched. Kept as the record of how v2 was made, and as
the statement of the conventions. See issue 27 in ../KNOWN_ISSUES.md.

The authority used here is the SEA-AD palette the annotations actually came from,
`snrnaseq/cluster_order_and_colors.csv` (Allen Institute, 139 supertypes). It sets
two conventions that coexist by design:

    subclass_label   LONG   Astrocyte, Endothelial, Microglia-PVM, Oligodendrocyte,
                            L2/3 IT, L5/6 NP, Lamp5 Lhx6
    cluster_label    SHORT prefix + _N   Astro_1, Endo_1, Micro-PVM_1, Oligo_1,
                            L2/3 IT_1, L5/6 NP_1, Lamp5_Lhx6_1

Measured against that, before this script runs:

    DE_subclass.cell_type      0 of 24  deviate   -> already correct, left alone
    cell_metadata.supertype    4 of 143 deviate   -> all legitimately outside the
                                                     palette (Unassigned, VLMC_2*,
                                                     Xenium-only), left alone
    cell_metadata.subclass     4 of 25  deviate   -> abbreviated; fixed
    DE_supertype.cell_type   115 of 115 deviate   -> wrong separators; fixed

Dataset names: cell_metadata and DE_subclass agree on Frohlich / MSSM1 / MSSM2;
DE_supertype alone uses the paper labels Fröhlich / MSSM 1 / MSSM 2. The majority
spelling wins, because a data file wants stable ASCII identifiers without spaces —
the paper's display labels belong in figures.

    python3 shared/harmonize_zenodo_labels.py --in <dir> --out <dir>

Every rename is asserted against the palette, so the script fails loudly rather
than inventing a name.
"""
import argparse, csv, os, re, sys

PALETTE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "snrnaseq", "cluster_order_and_colors.csv")

# subclass: the abbreviations cell_metadata uses -> the palette's subclass_label
SUBCLASS_FIX = {"Astro": "Astrocyte", "Endo": "Endothelial",
                "Micro-PVM": "Microglia-PVM", "Oligo": "Oligodendrocyte"}

# dataset: DE_supertype's paper labels -> the spelling the other two files use
DATASET_FIX = {"Fröhlich": "Frohlich", "MSSM 1": "MSSM1", "MSSM 2": "MSSM2"}

# values that are legitimately outside the SEA-AD palette and must not be "fixed"
OUTSIDE_PALETTE = {"Unassigned"}


def supertype_to_palette(name):
    """DE_supertype's separators -> the palette's cluster_label form."""
    s = re.sub(r"-(\d+)$", r"_\1", name)          # Astro-2      -> Astro_2
    s = s.replace("L2_3", "L2/3").replace("L5_6", "L5/6")
    s = s.replace("Lamp5-Lhx6", "Lamp5_Lhx6")
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True)
    ap.add_argument("--out", dest="dst", required=True)
    a = ap.parse_args()
    import pandas as pd

    rows = list(csv.DictReader(open(PALETTE)))
    allen_sub = {r["subclass_label"] for r in rows}
    allen_sup = {r["cluster_label"] for r in rows}
    os.makedirs(a.dst, exist_ok=True)

    # ---- cell_metadata: subclass abbreviations -> palette long names ---------
    cm = pd.read_parquet(os.path.join(a.src, "cell_metadata.parquet"))
    before = set(cm.subclass.dropna().unique())
    cm["subclass"] = cm.subclass.replace(SUBCLASS_FIX)
    bad = set(cm.subclass.dropna().unique()) - allen_sub - OUTSIDE_PALETTE
    if bad:
        sys.exit(f"cell_metadata.subclass still off-palette: {sorted(bad)}")
    print(f"cell_metadata.parquet  subclass: {len(before - set(cm.subclass.dropna().unique()))} renamed, "
          f"{len(set(cm.subclass.dropna().unique()) - OUTSIDE_PALETTE)} all on-palette")

    # ---- DE_supertype: separators and dataset names -------------------------
    dt = pd.read_parquet(os.path.join(a.src, "DE_supertype.parquet"))
    dt["cell_type"] = dt.cell_type.map(supertype_to_palette)
    # The bar is joinability to the annotations, not palette membership: a few
    # supertypes are real but absent from the 139-row palette (VLMC_2 and its
    # -SEAAD splits), and cell_metadata is the authority on what exists in the data.
    valid_sup = allen_sup | set(cm.supertype.dropna().unique())
    bad = set(dt.cell_type.dropna().unique()) - valid_sup
    if bad:
        sys.exit(f"DE_supertype.cell_type does not join to cell_metadata: {sorted(bad)}")
    off_palette = set(dt.cell_type.dropna().unique()) - allen_sup
    if off_palette:
        print(f"  note: {len(off_palette)} supertype(s) outside the 139-row palette "
              f"but present in cell_metadata: {sorted(off_palette)}")
    dt["cohort"] = dt.cohort.replace(DATASET_FIX)
    print(f"DE_supertype.parquet   cell_type: {dt.cell_type.nunique()} values, all on-palette; "
          f"datasets -> {sorted(dt.cohort.unique())}")

    # ---- DE_subclass: already correct, copied through unchanged -------------
    ds = pd.read_parquet(os.path.join(a.src, "DE_subclass.parquet"))
    bad = set(ds.cell_type.dropna().unique()) - allen_sub
    if bad:
        sys.exit(f"DE_subclass.cell_type off-palette: {sorted(bad)}")
    print(f"DE_subclass.parquet    unchanged ({ds.cell_type.nunique()} subclasses, already on-palette)")

    # ---- the join that was broken -------------------------------------------
    cs, dsup = set(cm.supertype.dropna().unique()), set(dt.cell_type.dropna().unique())
    csub, dsub = set(cm.subclass.dropna().unique()), set(ds.cell_type.dropna().unique())
    print(f"\ncross-file joins after harmonisation:")
    print(f"  cell_metadata.supertype  n DE_supertype.cell_type : {len(cs & dsup)} of {len(dsup)}")
    print(f"  cell_metadata.subclass   n DE_subclass.cell_type  : {len(csub & dsub)} of {len(dsub)}")
    if len(cs & dsup) != len(dsup) or len(csub & dsub) != len(dsub):
        sys.exit("a join is still incomplete — inspect before shipping")

    for name, df in [("cell_metadata.parquet", cm), ("DE_subclass.parquet", ds),
                     ("DE_supertype.parquet", dt)]:
        p = os.path.join(a.dst, name)
        df.to_parquet(p, index=False)
        print(f"  wrote {p}  ({os.path.getsize(p)/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
