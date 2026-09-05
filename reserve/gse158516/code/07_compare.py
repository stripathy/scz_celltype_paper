#!/usr/bin/env python3
"""
Is GSE158516 consistent with the paper's findings?

Tests, in the order the paper makes the claims:

  DE (Figure 2)
    1. SST down within Sst cells               (meta log2FC = -0.46, FDR = 0.049)
    2. PVALB not DE within Pvalb cells         (meta log2FC = -0.06, n.s.)
    3. Transcriptome-wide concordance of this cohort's per-subclass log2FCs with the
       7-dataset meta-analytic estimates -- the same comparison Fig. 2j makes for Xenium.
    4. The specific genes the paper names in Sst and Pvalb.

  Composition (Figure 3)
    5. The vulnerable Sst supertypes (Sst_2, Sst_3, Sst_20, Sst_22, Sst_25) depleted.
    6. L6b supertypes (L6b_1, L6b_2, L6b_4) increased.
    7. Per-supertype compositional effect correlated with the meta-analytic crumblr beta.
    8. Aggregate vulnerable-Sst proportion, control vs SCZ.

GSE158516 is an 8th, fully independent cohort: Reiner/Berrettini (UPenn) is not one of
the seven (Batiuk, HBCC, MSSM 1, MSSM 2, McLean, Multiome, Frohlich).
"""
import numpy as np
import pandas as pd
from scipy import stats

BASE = "/Users/shreejoy/Github/scz_celltype_paper"
OUTD = f"{BASE}/external_validation/gse158516/output"
SHARED = f"{BASE}/shared/snrnaseq_de"

VULNERABLE_SST = ["Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25"]
L6B = ["L6b_1", "L6b_2", "L6b_4"]

# SEA-AD subclass label -> the name used in the meta-analysis DE table
SUBCLASS_TO_META = {
    "Astrocyte": "Astro", "L2/3 IT": "L2_3 IT", "L5/6 NP": "L5_6 NP",
    "Microglia-PVM": "Micro-PVM", "Oligodendrocyte": "Oligo", "Endothelial": "Endo",
}


def hdr(t):
    print(f"\n{'=' * 78}\n{t}\n{'=' * 78}")


# ============================== DE =========================================
hdr("1-2. The two canonical interneuron markers")
de = pd.read_csv(f"{OUTD}/de_subclass.csv.gz")
de["ct_meta"] = de["cell_type"].map(lambda c: SUBCLASS_TO_META.get(c, c))
meta = pd.read_csv(f"{SHARED}/DE_genes_all_cells_scz.csv")

PAPER = {("SST", "Sst"): (-0.46, 0.049, "significant, down"),
         ("PVALB", "Pvalb"): (-0.06, 0.86, "not significant")}
rows = []
for (gene, ct), (mfc, mfdr, claim) in PAPER.items():
    r = de[(de["gene"] == gene) & (de["ct_meta"] == ct)]
    m = meta[(meta["genes"] == gene) & (meta["cell_type"] == ct)]
    rows.append(dict(gene=gene, cell_type=ct, paper_claim=claim,
                     meta_logFC=float(m["estimate"].iloc[0]) if len(m) else np.nan,
                     meta_FDR=float(m["padj"].iloc[0]) if len(m) else np.nan,
                     gse_logFC=float(r["logFC"].iloc[0]) if len(r) else np.nan,
                     gse_P=float(r["p"].iloc[0]) if len(r) else np.nan,
                     gse_FDR=float(r["FDR"].iloc[0]) if len(r) else np.nan))
marker_tab = pd.DataFrame(rows)
print(marker_tab.to_string(index=False))
marker_tab.to_csv(f"{OUTD}/compare_markers.csv", index=False)

hdr("3. Transcriptome-wide DE concordance with the 7-dataset meta-analysis")
j = de.merge(meta.rename(columns={"genes": "gene", "cell_type": "ct_meta",
                                  "estimate": "meta_logFC", "padj": "meta_FDR",
                                  "pval": "meta_P"}),
             on=["gene", "ct_meta"], how="inner")
print(f"gene x subclass pairs testable in both: {len(j):,} "
      f"({j['ct_meta'].nunique()} subclasses)")

for label, sub in [("all shared pairs", j),
                   ("meta FDR < 0.10", j[j["meta_FDR"] < 0.10]),
                   ("meta FDR < 0.05", j[j["meta_FDR"] < 0.05])]:
    if len(sub) < 10:
        print(f"{label:22s} n={len(sub)} -- too few to test")
        continue
    r, rp = stats.pearsonr(sub["meta_logFC"], sub["logFC"])
    rho, rhop = stats.spearmanr(sub["meta_logFC"], sub["logFC"])
    same = int(np.sum(np.sign(sub["meta_logFC"]) == np.sign(sub["logFC"])))
    sp = stats.binomtest(same, len(sub), 0.5).pvalue
    slope = np.polyfit(sub["meta_logFC"], sub["logFC"], 1)[0]
    print(f"{label:22s} n={len(sub):6d}  Pearson r={r:+.3f} (P={rp:.2g})  "
          f"Spearman rho={rho:+.3f}  slope={slope:+.2f}  "
          f"sign-concordant {same}/{len(sub)} = {same/len(sub):.1%} (P={sp:.2g})")
j.to_csv(f"{OUTD}/compare_de_joined.csv.gz", index=False, compression="gzip")

hdr("4. The specific genes named in the paper")
named = [("Sst", g) for g in ["AFG3L2", "NAT16", "SLC9A9", "STAC", "SMAD1", "ATP2B4", "DRD3", "KCTD4"]] + \
        [("Pvalb", g) for g in ["ANXA2", "NAT16", "VGF", "CIRBP", "SCN3A", "SMAD1", "TCAF2"]] + \
        [("Chandelier", "VGF")]
rows = []
for ct, g in named:
    r = de[(de["gene"] == g) & (de["ct_meta"] == ct)]
    m = meta[(meta["genes"] == g) & (meta["cell_type"] == ct)]
    if len(r) and len(m):
        rows.append(dict(cell_type=ct, gene=g,
                         meta_logFC=float(m["estimate"].iloc[0]), meta_FDR=float(m["padj"].iloc[0]),
                         gse_logFC=float(r["logFC"].iloc[0]), gse_P=float(r["p"].iloc[0]),
                         same_sign=np.sign(m["estimate"].iloc[0]) == np.sign(r["logFC"].iloc[0])))
nt = pd.DataFrame(rows)
print(nt.round(4).to_string(index=False))
if len(nt):
    print(f"\ndirectional agreement: {nt.same_sign.sum()}/{len(nt)}")
nt.to_csv(f"{OUTD}/compare_named_genes.csv", index=False)

# ========================== COMPOSITION ====================================
comp = pd.read_csv(f"{OUTD}/composition_crumblr.csv")
mc = pd.concat([pd.read_csv(f"{SHARED}/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv"),
                pd.read_csv(f"{SHARED}/nicole_scz_snrnaseq_betas/final_results_crumblr_7_nonN_cohorts.csv")])
mc_sub = pd.concat([pd.read_csv(f"{SHARED}/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts_subclass.csv"),
                    pd.read_csv(f"{SHARED}/nicole_scz_snrnaseq_betas/final_results_crumblr_7_nonN_cohorts_subclass.csv")])

hdr("5-6. The paper's headline supertypes")
st = comp[comp["level"] == "supertype"].set_index("cell_type")
mcs = mc.set_index("CellType")
rows = []
for t in VULNERABLE_SST + L6B:
    rows.append(dict(supertype=t,
                     meta_beta=mcs.loc[t, "estimate"] if t in mcs.index else np.nan,
                     meta_FDR=mcs.loc[t, "padj"] if t in mcs.index else np.nan,
                     gse_logFC=st.loc[t, "logFC"] if t in st.index else np.nan,
                     gse_P=st.loc[t, "P.Value"] if t in st.index else np.nan,
                     gse_FDR=st.loc[t, "FDR"] if t in st.index else np.nan,
                     prop_control=st.loc[t, "mean_prop_control"] if t in st.index else np.nan,
                     prop_scz=st.loc[t, "mean_prop_scz"] if t in st.index else np.nan))
head = pd.DataFrame(rows)
head["same_sign"] = np.sign(head.meta_beta) == np.sign(head.gse_logFC)
print(head.round(4).to_string(index=False))
print(f"\ndirectional agreement on the 8 headline supertypes: "
      f"{head.same_sign.sum()}/{head.same_sign.notna().sum()}")
head.to_csv(f"{OUTD}/compare_headline_supertypes.csv", index=False)

hdr("7. Compositional effect-size concordance across all matched cell types")
for level, mtab in [("supertype", mc), ("subclass", mc_sub)]:
    x = (comp[comp["level"] == level]
         .merge(mtab.rename(columns={"CellType": "cell_type", "estimate": "meta_beta",
                                     "padj": "meta_FDR"}), on="cell_type", how="inner"))
    for scope, sub in [("all", x), ("neuronal", x[x["compartment"] == "neuronal"])]:
        if len(sub) < 5:
            continue
        r, p = stats.pearsonr(sub["meta_beta"], sub["logFC"])
        rho, rp = stats.spearmanr(sub["meta_beta"], sub["logFC"])
        same = int(np.sum(np.sign(sub["meta_beta"]) == np.sign(sub["logFC"])))
        print(f"{level:10s} {scope:9s} n={len(sub):4d}  Pearson r={r:+.3f} (P={p:.3g})  "
              f"Spearman rho={rho:+.3f} (P={rp:.3g})  sign {same}/{len(sub)}")
    x.to_csv(f"{OUTD}/compare_composition_{level}.csv", index=False)

hdr("8. Aggregate vulnerable-Sst proportion, control vs SCZ")
counts = pd.read_csv(f"{OUTD}/counts_supertype.csv", index_col=0)
sm = pd.read_csv(f"{BASE}/external_validation/gse158516/data/sample_metadata.csv")
sm = sm.set_index("Deidentified ID").loc[counts.index]
comp_map = (pd.read_csv(f"{OUTD}/celltype_compartment.csv")
            .drop_duplicates("cell_type").set_index("cell_type")["compartment"])
neuronal = [c for c in counts.columns if comp_map.get(c) == "neuronal"]

for label, types in [("vulnerable Sst", VULNERABLE_SST), ("all Sst", [c for c in counts.columns if c.startswith("Sst_")]),
                     ("L6b", [c for c in counts.columns if c.startswith("L6b_")])]:
    present = [t for t in types if t in counts.columns]
    prop = counts[present].sum(1) / counts[neuronal].sum(1)
    a = prop[sm["Diagnosis"] == "Control"]; b = prop[sm["Diagnosis"] == "Schizophrenia"]
    t, p = stats.ttest_ind(a, b, equal_var=False)
    u, up = stats.mannwhitneyu(a, b)
    print(f"{label:16s} ({len(present)} types)  control {a.mean():.4%}  SCZ {b.mean():.4%}  "
          f"ratio {b.mean()/a.mean():.2f}  Welch P={p:.4f}  MWU P={up:.4f}")
