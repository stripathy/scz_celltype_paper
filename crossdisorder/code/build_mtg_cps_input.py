"""
Per-donor neuronal supertype counts for the SEA-AD MTG cohort, as a robustness
counterpart to the DLPFC input behind Figure 4i.

Mirrors crossdisorder/code/build_seaad_cps_input.py exactly -- neurons only,
donor-level neuronal total as the denominator, neurotypical reference donors
dropped via their missing CPS -- but reads the MTG release instead of A9.

PMI is absent from the MTG metadata table, so it is joined by donor from the
DLPFC input. All 80 DLPFC donors are present in the MTG cohort, so the two runs
end up on the same donors with the same covariates and differ only in which
region was dissected.
"""
import os

import pandas as pd

XEN = os.environ.get("XENIUM_BASE", os.path.expanduser("~/Github/SCZ_Xenium"))
XD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MTG_OBS = (f"{XEN}/data/reference/"
           "SEAAD_MTG_RNAseq_obs_lean.csv")
DLPFC_IN = (f"{XD}/results/"
            "crumblr_input_supertype_neurons.csv")
OUT = (f"{XD}/results/"
       "crumblr_input_mtg_supertype_neurons.csv")

NEURONAL_CLASSES = ("Neuronal: Glutamatergic", "Neuronal: GABAergic")

obs = pd.read_csv(MTG_OBS, low_memory=False, usecols=[
    "Donor ID", "Class", "Supertype", "Continuous Pseudo-progression Score",
    "Sex", "Age at Death", "Severely Affected Donor"])
obs = obs.rename(columns={"Donor ID": "donor", "Supertype": "celltype",
                          "Continuous Pseudo-progression Score": "CPS",
                          "Sex": "sex", "Age at Death": "age",
                          "Severely Affected Donor": "severely_affected"})
print(f"  nuclei: {len(obs):,}   donors: {obs.donor.nunique()}")

meta = (obs.groupby("donor", observed=True)
        .agg(CPS=("CPS", "first"), sex=("sex", "first"), age=("age", "first"),
             severely_affected=("severely_affected", "first")).reset_index())
meta["CPS"] = pd.to_numeric(meta.CPS, errors="coerce")
# "90+" is recorded for the oldest donors; coerce as the DLPFC builder does
meta["age"] = pd.to_numeric(meta.age.astype(str).str.replace("+", "", regex=False),
                            errors="coerce")

dl = pd.read_csv(DLPFC_IN)
meta = meta.merge(dl[["donor", "pmi"]].drop_duplicates(), on="donor", how="left")

neurons = obs[obs["Class"].isin(NEURONAL_CLASSES)]
print(f"  neuronal nuclei: {len(neurons):,} of {len(obs):,}")
counts = (neurons.groupby(["donor", "celltype"], observed=True)
          .size().rename("count").reset_index())
counts["total"] = counts.groupby("donor")["count"].transform("sum")

df = (counts.merge(meta, on="donor", how="left")
      .dropna(subset=["CPS", "sex", "age", "pmi"]))
df.to_csv(OUT, index=False)
print(f"  final: {df.donor.nunique()} donors, {df.celltype.nunique()} neuronal "
      f"supertypes, {int(df['count'].sum()):,} nuclei")
print(f"  wrote {OUT}")
