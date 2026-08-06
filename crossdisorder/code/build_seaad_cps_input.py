#!/usr/bin/env python3
"""
Build the crumblr input for the SEA-AD DLPFC (A9) pseudo-progression analysis
behind Figure 4j.

Companion to 14_seaad_dlpfc_crumblr.R, which fits the model.

What this does
--------------
Reads the SEA-AD A9 per-nucleus metadata table, restricts to neurons, and emits
a long-format per-donor x supertype count table with the donor covariates the
model needs.

Design choices, and why
-----------------------
* **A9 (DLPFC), not MTG.** Every schizophrenia dataset in this paper is frontal
  cortex, so putting the AD axis in DLPFC makes both axes of panel 4j regionally
  matched. The MTG version of this analysis reproduces it (rho = 0.79 vs 0.87
  here; the two regions agree with each other at rho = 0.95).

* **The 2024-02-13 release, not the 2026 re-annotation.** The 2026 DFC release
  annotates against an expanded taxonomy (152 supertypes, 15 of them novel
  `-SEAAD` disease states) that redistributes cells relative to the 137-supertype
  taxonomy used throughout this paper. The 2024-02-13 `Supertype` column is a
  strict subset of the MTG labels -- zero novel types -- and shares all 16 Sst
  supertypes with the rest of our analyses.

* **Neurons only.** SEA-AD nuclei were FACS-sorted to a fixed 70:30
  neuronal:non-neuronal ratio, so Gabitto et al. analyse the two compartments
  separately; the denominator here is all neuronal nuclei in a donor, which also
  matches the "% of all neurons" denominator on the schizophrenia side.

* **Severely affected donors are retained.** Gabitto et al. excluded these 9-11
  donors only from gene-expression tests (they show global transcriptional
  shutdown and lower library quality) and kept them for the compositional
  analyses, reporting consistent results either way. We follow that convention.
  Excluding them here moves the panel-4j correlation from rho = 0.87 to 0.77.

* **Neurotypical reference donors are dropped.** Three A9 donors (H18.30.002,
  H19.30.001, H19.30.002) are the Allen/BICCN reference brains -- no age, no
  PMI (the field reads "Reference"), and no CPS, because they are not on the AD
  pseudo-progression trajectory. 83 donors -> 80.

* **CPS comes from a committed lookup.** CPS is donor-constant and derived from
  quantitative neuropathology; it is absent from the A9 metadata. See
  data/seaad_dlpfc/README.md for provenance.

Usage
-----
    python3 build_seaad_cps_input.py [--metadata PATH]

Output
------
    results/intermediates/seaad_dlpfc/crumblr_input_supertype_neurons.csv
        donor, celltype, count, total, CPS, sex, age, pmi, severely_affected
"""

import argparse
import os

import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HERE = os.path.dirname(os.path.abspath(__file__))
COMP = os.path.dirname(HERE)               # crossdisorder/
DATA = os.path.join(COMP, "data")
OUT = os.path.join(COMP, "results")

# Default location; overridable with --metadata. See data/seaad_dlpfc/README.md
# for the download URL -- the file is ~1.3 GB and is not committed.
# Shared external reference, alongside the other cross-component inputs.
# Override with --metadata, or set SEAAD_METADATA.
DEFAULT_METADATA = os.environ.get(
    "SEAAD_METADATA",
    os.path.join(REPO, "shared", "seaad",
                 "SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv"))
CPS_LOOKUP = os.path.join(DATA, "seaad_donor_cps.csv")

NEURONAL_CLASSES = ("Neuronal: Glutamatergic", "Neuronal: GABAergic")

USECOLS = ["Donor ID", "Brain Region", "Sex", "Age at Death", "PMI",
           "Class", "Subclass", "Supertype", "Severely Affected Donor"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", default=DEFAULT_METADATA,
                    help="SEA-AD A9 per-nucleus metadata CSV (2024-02-13 release)")
    args = ap.parse_args()

    if not os.path.exists(args.metadata):
        raise SystemExit(
            f"Metadata not found: {args.metadata}\n"
            f"See {os.path.join(DATA, 'README.md')} for the download URL.")

    os.makedirs(OUT, exist_ok=True)

    print(f"Reading {args.metadata}")
    obs = pd.read_csv(args.metadata, usecols=USECOLS, low_memory=False)
    print(f"  {len(obs):,} nuclei, {obs['Donor ID'].nunique()} donors, "
          f"region(s): {sorted(obs['Brain Region'].unique())}")

    # ---- donor covariates -------------------------------------------------
    cps = pd.read_csv(CPS_LOOKUP)
    meta = (obs.groupby("Donor ID")
               .agg(sex=("Sex", "first"),
                    age=("Age at Death", "first"),
                    pmi=("PMI", "first"),
                    severely_affected=("Severely Affected Donor", "first"))
               .reset_index().rename(columns={"Donor ID": "donor"})
               .merge(cps, on="donor", how="left"))
    # SEA-AD codes the oldest bracket as "90+"; PMI reads "Reference" for the
    # neurotypical reference donors, which coerces to NaN and drops them below.
    meta["age"] = pd.to_numeric(
        meta.age.astype(str).str.replace("90+", "90", regex=False), errors="coerce")
    meta["pmi"] = pd.to_numeric(meta.pmi, errors="coerce")

    dropped = meta[meta.CPS.isna()].donor.tolist()
    if dropped:
        print(f"  dropping {len(dropped)} donor(s) without CPS "
              f"(neurotypical reference): {', '.join(sorted(dropped))}")

    # ---- neuronal counts --------------------------------------------------
    neurons = obs[obs["Class"].isin(NEURONAL_CLASSES)]
    print(f"  neuronal nuclei: {len(neurons):,} of {len(obs):,}")

    counts = (neurons.groupby(["Donor ID", "Supertype"], observed=True)
                     .size().reset_index(name="count")
                     .rename(columns={"Donor ID": "donor", "Supertype": "celltype"}))
    counts = counts.merge(
        counts.groupby("donor")["count"].sum().rename("total"), on="donor")

    df = (counts.merge(meta, on="donor", how="inner")
                .dropna(subset=["CPS", "sex", "age", "pmi"]))

    n_donors = df.donor.nunique()
    n_nuclei = int(df.groupby("donor")["total"].first().sum())
    print(f"  final: {n_donors} donors, {df.celltype.nunique()} neuronal supertypes, "
          f"{n_nuclei:,} neuronal nuclei")
    print(f"  severely affected retained: "
          f"{(df.groupby('donor')['severely_affected'].first() == 'Y').sum()} donors")

    out_path = os.path.join(OUT, "crumblr_input_supertype_neurons.csv")
    (df[["donor", "celltype", "count", "total", "CPS", "sex", "age", "pmi",
         "severely_affected"]]
       .sort_values(["donor", "celltype"])
       .to_csv(out_path, index=False))
    print(f"  -> {out_path}")


if __name__ == "__main__":
    main()
