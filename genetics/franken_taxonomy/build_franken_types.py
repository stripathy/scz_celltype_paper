#!/usr/bin/env python3
"""
Build franken_types.csv — the canonical 598-row provenance table for the
combined SEA-AD ⊕ Siletti ("Franken") taxonomy.

Each row is one candidate type (137 SEA-AD supertypes + 461 Siletti
clusters). Columns record the type's source, its RBH partner if any, and
whether it survived RBH dedup into the 503-type combined Franken taxonomy.

Inputs (relative to repo root):
  - results/tables/metaneighbor_full461_reciprocal_best_hits.csv
      The 95 reciprocal-best-hit pairs (output of step 09_metaneighbor_integration.py).
  - results/tables/metaneighbor_full461_all_matches.csv
      All 137 SEA-AD types with their best Siletti match (whether reciprocal or not).
  - results/tables/rbh_combined_enrichment.csv
      The 503-type combined-taxonomy enrichment table — used as the
      authoritative list of "what is in the combined Franken."

Output:
  - franken_taxonomy/franken_types.csv

Schema:
  type                  str   — type name (e.g. Pvalb_3, MGE_259, Misc_132)
  source                str   — "seaad" or "siletti"
  rbh_partner           str   — partner type name across datasets, if RBH (else "")
  rbh_mean_auroc        float — mean AUROC of the RBH pair (else NaN)
  best_match            str   — best cross-dataset match (filled even when not RBH)
  best_match_auroc      float — AUROC of the best cross-dataset match
  is_reciprocal         bool  — true iff this type is part of an RBH pair
  in_combined_franken   bool  — true iff this type is one of the 503 kept types
                                (= all SEA-AD; Siletti minus the 95 RBH-claimed)

Run from the repo root:
    python3 franken_taxonomy/build_franken_types.py
"""
from __future__ import annotations
import os
import sys
import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES = os.path.join(REPO_ROOT, "results", "tables")
RBH_PATH = os.path.join(TABLES, "metaneighbor_full461_reciprocal_best_hits.csv")
ALL_MATCHES_PATH = os.path.join(TABLES, "metaneighbor_full461_all_matches.csv")
COMBINED_PATH = os.path.join(TABLES, "rbh_combined_enrichment.csv")
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "franken_types.csv")


def main():
    rbh = pd.read_csv(RBH_PATH)
    all_matches = pd.read_csv(ALL_MATCHES_PATH)
    combined = pd.read_csv(COMBINED_PATH)

    rbh_seaad = set(rbh["seaad_type"])
    rbh_siletti_to_seaad = dict(zip(rbh["siletti_cluster"], rbh["seaad_type"]))
    rbh_seaad_to_siletti = dict(zip(rbh["seaad_type"], rbh["siletti_cluster"]))
    rbh_mean_auroc = dict(zip(rbh["seaad_type"], rbh["mean_auroc"]))
    rbh_claimed_siletti = set(rbh["siletti_cluster"])

    seaad_types = sorted(set(all_matches["seaad_type"]))
    siletti_clusters = sorted(set(all_matches["best_siletti"]) | rbh_claimed_siletti)
    siletti_extra = set(combined.loc[combined["source"] == "Siletti", "supertype"])
    siletti_clusters = sorted(set(siletti_clusters) | siletti_extra)

    sea_best = dict(zip(all_matches["seaad_type"], all_matches["best_siletti"]))
    sea_best_auroc = dict(zip(all_matches["seaad_type"], all_matches["best_auroc"]))

    in_combined = set(combined["supertype"])

    rows = []
    for t in seaad_types:
        is_recip = t in rbh_seaad
        rows.append({
            "type": t,
            "source": "seaad",
            "rbh_partner": rbh_seaad_to_siletti.get(t, ""),
            "rbh_mean_auroc": rbh_mean_auroc.get(t, float("nan")),
            "best_match": sea_best.get(t, ""),
            "best_match_auroc": sea_best_auroc.get(t, float("nan")),
            "is_reciprocal": is_recip,
            "in_combined_franken": t in in_combined,
        })
    for c in siletti_clusters:
        is_recip = c in rbh_claimed_siletti
        rows.append({
            "type": c,
            "source": "siletti",
            "rbh_partner": rbh_siletti_to_seaad.get(c, ""),
            "rbh_mean_auroc": (rbh_mean_auroc.get(rbh_siletti_to_seaad.get(c, ""),
                                                    float("nan"))
                                if is_recip else float("nan")),
            "best_match": "",
            "best_match_auroc": float("nan"),
            "is_reciprocal": is_recip,
            "in_combined_franken": c in in_combined,
        })

    df = pd.DataFrame(rows)
    df = df.sort_values(["source", "type"]).reset_index(drop=True)

    n_seaad = (df["source"] == "seaad").sum()
    n_siletti = (df["source"] == "siletti").sum()
    n_kept = df["in_combined_franken"].sum()
    n_dropped_siletti = ((df["source"] == "siletti") & ~df["in_combined_franken"]).sum()
    print(f"  SEA-AD types:           {n_seaad}")
    print(f"  Siletti clusters:       {n_siletti}")
    print(f"  In combined Franken:    {n_kept}")
    print(f"  Siletti dropped (RBH):  {n_dropped_siletti}")
    if n_kept != n_seaad + n_siletti - n_dropped_siletti:
        print("  WARN: combined-taxonomy size doesn't match SEA-AD + novel-Siletti")

    df.to_csv(OUT_PATH, index=False)
    print(f"\nWrote {OUT_PATH} ({len(df)} rows)")


if __name__ == "__main__":
    sys.exit(main())
