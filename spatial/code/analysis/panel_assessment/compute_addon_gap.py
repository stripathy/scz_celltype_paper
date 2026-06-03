#!/usr/bin/env python3
"""
Compute how many add-on genes would be needed to bring Xenium v1
up to MERSCOPE 250-level supertype marker coverage.

Uses the within-subclass marker results from cross_platform_marker_adequacy.py
"""

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR, PANEL_V1_PATH

OUTPUT_DIR = os.path.join(BASE_DIR, "output", "marker_analysis")

def main():
    # Load coverage data
    cov = pd.read_csv(os.path.join(OUTPUT_DIR, "cross_platform_marker_coverage.csv"))

    # Load all within-subclass markers
    markers = pd.read_csv(os.path.join(OUTPUT_DIR, "within_subclass_markers_all.csv"))

    # Load Xenium v1 panel
    v1 = pd.read_csv(PANEL_V1_PATH)
    v1_genes = set(v1["Genes"].str.strip())

    # Load MERSCOPE 250 panel
    m250 = pd.read_csv(os.path.join(BASE_DIR, "data", "merscope_4k_probe_testing",
                                      "H18.06.006.MTG.250.expand.rep1.genes.csv"))
    m250_genes = set(m250.iloc[:, 0].str.strip())

    print(f"Xenium v1: {len(v1_genes)} genes")
    print(f"MERSCOPE 250: {len(m250_genes)} genes")
    print(f"Overlap: {len(v1_genes & m250_genes)} genes")

    # For each supertype, find markers that are in MERSCOPE 250 but NOT in Xenium v1
    # These are the genes that explain the coverage gap

    print("\n" + "=" * 80)
    print("GAP ANALYSIS: genes in MERSCOPE 250 that would improve Xenium v1 coverage")
    print("=" * 80)

    gap_genes = set()  # genes that would close the gap
    gap_records = []

    supertypes = sorted(markers["group"].unique())

    for st in supertypes:
        st_markers = markers[markers["group"] == st].copy()
        subclass = st_markers["subclass"].iloc[0]

        for topn in [10, 20, 50]:
            top = st_markers.head(topn)
            in_v1 = set(top[top["gene"].isin(v1_genes)]["gene"])
            in_m250 = set(top[top["gene"].isin(m250_genes)]["gene"])
            in_m250_not_v1 = in_m250 - v1_genes

            gap_records.append({
                "supertype": st,
                "subclass": subclass,
                "topN": topn,
                "n_in_v1": len(in_v1),
                "n_in_m250": len(in_m250),
                "n_gap": len(in_m250_not_v1),
                "gap_genes": "; ".join(sorted(in_m250_not_v1)),
            })

            if topn == 10:
                gap_genes.update(in_m250_not_v1)

    gap_df = pd.DataFrame(gap_records)

    # Summary for top-10
    top10 = gap_df[gap_df["topN"] == 10]
    print(f"\nTop-10 within-subclass markers:")
    print(f"  Xenium v1 mean coverage: {top10['n_in_v1'].mean():.1f}/10")
    print(f"  MERSCOPE 250 mean coverage: {top10['n_in_m250'].mean():.1f}/10")
    print(f"  Genes in M250 but not v1 (top10): {len(gap_genes)} unique genes")

    # Now compute: if we add these gap genes to v1, what's the new coverage?
    v1_plus_gap = v1_genes | gap_genes
    print(f"\n  If we add these {len(gap_genes)} genes to v1:")
    print(f"    New panel size: {len(v1_plus_gap)} genes")

    new_coverage = []
    for st in supertypes:
        st_markers = markers[markers["group"] == st]
        top10 = st_markers.head(10)
        n_new = top10["gene"].isin(v1_plus_gap).sum()
        n_orig = top10["gene"].isin(v1_genes).sum()
        new_coverage.append({"supertype": st, "v1_orig": n_orig,
                            "v1_plus_gap": n_new})

    new_cov_df = pd.DataFrame(new_coverage)
    print(f"    Mean top-10 coverage: {new_cov_df['v1_orig'].mean():.1f} → "
          f"{new_cov_df['v1_plus_gap'].mean():.1f}/10")

    # But that's only covering the M250 gap. Let's do a more systematic analysis:
    # For each supertype, greedily select markers until we reach M250's coverage level

    print("\n" + "=" * 80)
    print("GREEDY GENE SELECTION: minimum genes to match MERSCOPE 250 coverage")
    print("=" * 80)

    # Strategy: for each supertype, find how many top-N markers are in M250.
    # Then greedily select the top-ranked markers NOT in v1 until we match.

    needed_genes = set()
    per_type_needs = []

    for st in supertypes:
        st_markers = markers[markers["group"] == st]
        subclass = st_markers["subclass"].iloc[0]

        # M250 coverage at top-10
        top10 = st_markers.head(10)
        m250_cov = top10["gene"].isin(m250_genes).sum()
        v1_cov = top10["gene"].isin(v1_genes).sum()

        # Target: match m250_cov using top-50 markers
        # How many of top-50 markers are NOT in v1 but would help?
        top50 = st_markers.head(50)
        candidates = top50[~top50["gene"].isin(v1_genes)].copy()

        n_needed = max(0, m250_cov - v1_cov)
        selected = list(candidates.head(n_needed)["gene"])
        needed_genes.update(selected)

        per_type_needs.append({
            "supertype": st,
            "subclass": subclass,
            "v1_cov_top10": v1_cov,
            "m250_cov_top10": m250_cov,
            "deficit": n_needed,
            "genes_to_add": "; ".join(selected),
        })

    needs_df = pd.DataFrame(per_type_needs)
    needs_df.to_csv(os.path.join(OUTPUT_DIR, "addon_gap_analysis.csv"), index=False)

    print(f"\n  To match MERSCOPE 250 top-10 coverage per supertype:")
    print(f"    Total unique genes to add: {len(needed_genes)}")
    print(f"    New panel size: {len(v1_genes) + len(needed_genes)}")

    # Verify
    v1_plus = v1_genes | needed_genes
    verify = []
    for st in supertypes:
        st_markers = markers[markers["group"] == st]
        top10 = st_markers.head(10)
        n_new = top10["gene"].isin(v1_plus).sum()
        n_m250 = top10["gene"].isin(m250_genes).sum()
        verify.append({"supertype": st, "v1_plus": n_new, "m250": n_m250})
    verify_df = pd.DataFrame(verify)
    print(f"    Verified mean coverage: {verify_df['v1_plus'].mean():.1f}/10 "
          f"(MERSCOPE 250: {verify_df['m250'].mean():.1f}/10)")

    # Show supertypes with biggest deficits
    big_deficit = needs_df[needs_df["deficit"] > 0].sort_values("deficit", ascending=False)
    print(f"\n  Supertypes with deficit > 0: {len(big_deficit)}/{len(supertypes)}")
    print(f"\n  {'Supertype':20s} {'Sub':12s} {'v1':>4s} {'M250':>4s} {'Gap':>4s} Genes to add")
    print("  " + "-" * 90)
    for _, r in big_deficit.head(30).iterrows():
        print(f"  {r['supertype']:20s} {r['subclass']:12s} {r['v1_cov_top10']:4d} "
              f"{r['m250_cov_top10']:4d} {r['deficit']:4d} {r['genes_to_add']}")

    # ── Also compute: what if we want ALL supertypes at ≥3 markers? ──
    print("\n" + "=" * 80)
    print("MINIMUM GENES FOR ≥3 WITHIN-SUBCLASS MARKERS PER SUPERTYPE")
    print("=" * 80)

    threshold_genes = set()
    threshold_records = []

    for st in supertypes:
        st_markers = markers[markers["group"] == st]
        top50 = st_markers.head(50)

        current = top50["gene"].isin(v1_genes).sum()
        if current >= 3:
            threshold_records.append({
                "supertype": st,
                "current_v1": current,
                "deficit": 0,
                "genes_to_add": "",
            })
            continue

        needed = 3 - min(current, 3)
        candidates = top50[~top50["gene"].isin(v1_genes)]
        selected = list(candidates.head(needed)["gene"])
        threshold_genes.update(selected)
        threshold_records.append({
            "supertype": st,
            "current_v1": current,
            "deficit": needed,
            "genes_to_add": "; ".join(selected),
        })

    # But we need to check top-50 coverage, not just top-10!
    # Some supertypes might have 0 in top-10 but 3 in top-50
    # Recount using top-50
    for rec in threshold_records:
        st = rec["supertype"]
        st_markers = markers[markers["group"] == st]
        top50 = st_markers.head(50)
        rec["current_v1_top50"] = int(top50["gene"].isin(v1_genes).sum())

    thresh_df = pd.DataFrame(threshold_records)
    n_below = (thresh_df["current_v1_top50"] < 3).sum()
    print(f"\n  Supertypes with < 3 markers in v1 (top-50 window): {n_below}/{len(supertypes)}")
    print(f"  Genes needed to bring all to ≥3: {len(threshold_genes)}")
    print(f"  New panel size: {len(v1_genes) + len(threshold_genes)}")

    below3 = thresh_df[thresh_df["current_v1_top50"] < 3].sort_values("current_v1_top50")
    if len(below3) > 0:
        print(f"\n  {'Supertype':20s} {'top50_in_v1':>10s} {'deficit':>7s} Genes to add")
        print("  " + "-" * 70)
        for _, r in below3.iterrows():
            print(f"  {r['supertype']:20s} {r['current_v1_top50']:10d} "
                  f"{r['deficit']:7d} {r['genes_to_add']}")

    # ── Print the actual gene list ──
    print(f"\n" + "=" * 80)
    print(f"GENE LIST: {len(needed_genes)} genes to match MERSCOPE 250 coverage")
    print("=" * 80)

    # Count how many supertypes each gene helps
    gene_utility = {}
    for st in supertypes:
        st_markers = markers[markers["group"] == st]
        top10 = st_markers.head(10)
        for g in needed_genes:
            if g in top10["gene"].values:
                gene_utility[g] = gene_utility.get(g, 0) + 1

    sorted_genes = sorted(needed_genes, key=lambda g: -gene_utility.get(g, 0))
    print(f"\n  Gene (# supertypes it helps):")
    for g in sorted_genes:
        n = gene_utility.get(g, 0)
        print(f"    {g:20s} → {n} supertypes")


if __name__ == "__main__":
    main()
