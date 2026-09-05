#!/usr/bin/env python3
"""
Check every numerical claim in README.md against the files in output/.

House rule: no statistic goes into a document unless it has been re-derived from the
source data. Any FAIL below means the README is wrong, not that the check is wrong.
"""
import re

import numpy as np
import pandas as pd
from scipy import stats

BASE = "/Users/shreejoy/Github/scz_celltype_paper/reserve/gse158516"
OUTD = f"{BASE}/output"
SHARED = "/Users/shreejoy/Github/scz_celltype_paper/shared/snrnaseq_de"

ok = fail = 0


def check(label, got, want, tol=5e-3, rel=False):
    global ok, fail
    good = (abs(got - want) <= tol * max(1e-12, abs(want))) if rel else (abs(got - want) <= tol)
    print(f"{'PASS' if good else 'FAIL'}  {label:58s} got {got:<14.6g} readme {want:<14.6g}")
    ok, fail = ok + good, fail + (not good)


meta = pd.read_csv(f"{BASE}/data/sample_metadata.csv")
qc = pd.read_csv(f"{OUTD}/qc_per_sample.csv")
d = qc.merge(meta, left_on="sample", right_on="Deidentified ID", how="left")
sub = d[d.in_paper]

print("--- cohort ---")
check("n donors in supplement", len(meta), 26, 0)
check("all male", (meta.Sex == "Male").sum(), 26, 0)
check("CT mean age", meta[meta.Diagnosis == "Control"].Age.mean(), 54.4, 0.05)
check("SCZ mean age", meta[meta.Diagnosis == "Schizophrenia"].Age.mean(), 41.7, 0.05)
ct, sz = (meta[meta.Diagnosis == g].Age for g in ["Control", "Schizophrenia"])
check("age Welch p (all 26)", stats.ttest_ind(ct, sz, equal_var=False).pvalue, 0.064, 5e-4)
check("CT PMI", meta[meta.Diagnosis == "Control"].PMI.mean(), 22.0, 0.05)
check("SCZ PMI", meta[meta.Diagnosis == "Schizophrenia"].PMI.mean(), 29.1, 0.05)
u = meta[meta.Age <= 70]
check("CT meeting <=70", (u.Diagnosis == "Control").sum(), 9, 0)
check("SCZ meeting <=70", (u.Diagnosis == "Schizophrenia").sum(), 12, 0)
check("fraction <=70 overall", (meta.Age <= 70).mean(), 0.808, 1e-3)
uct, usz = (u[u.Diagnosis == g].Age for g in ["Control", "Schizophrenia"])
check("age Welch p (<=70)", stats.ttest_ind(uct, usz, equal_var=False).pvalue, 1.000, 1e-3)

print("\n--- cell calling / QC ---")
cc = pd.read_csv(f"{OUTD}/cell_calling_comparison.csv")
ci = cc[cc.in_paper].merge(meta, left_on="sample", right_on="Deidentified ID")
ratio = ci.med_umi_topN / ci["Median UMI per Nuclei"]
check("median-UMI ratio ours/theirs", float(np.median(ratio)), 1.0001, 1e-3)
check("samples within 1%", int((100 * (ratio - 1).abs() < 1).sum()), 23, 0)
check("max % deviation", float((100 * (ratio - 1).abs()).max()), 3.2, 0.05)
check("nuclei called (26)", int(sub.n_called.sum()), 361681, 0)
check("== authors' reported total", int(meta["Number of Nuclei"].sum()), 361681, 0)
check("nuclei after QC (26)", int(sub.n_pass.sum()), 259187, 0)
check("nuclei after QC (32)", int(qc.n_pass.sum()), 294974, 0)
check("mean removal", sub.frac_removed.mean(), 0.285, 1e-3)
a, b = (sub[sub.Diagnosis == g] for g in ["Control", "Schizophrenia"])
check("removal CT", a.frac_removed.mean(), 0.282, 1e-3)
check("removal SCZ", b.frac_removed.mean(), 0.289, 1e-3)
check("removal Welch p", stats.ttest_ind(a.frac_removed, b.frac_removed, equal_var=False).pvalue, 0.80, 5e-3)
check("mito-fail Welch p", stats.ttest_ind(a.frac_fail_mito, b.frac_fail_mito, equal_var=False).pvalue, 0.71, 5e-3)

print("\n--- annotation ---")
o = pd.read_parquet("/Users/shreejoy/Github/shared_data/GSE158516/annotated/query_obs.parquet")
c = pd.read_parquet("/Users/shreejoy/Github/shared_data/GSE158516/annotated/query_obs_capped.parquet")
check("subclass conf >0.8", np.mean(o.subclass_confidence > 0.8), 0.939, 1e-3)
check("supertype conf median", float(o.supertype_confidence.median()), 0.75, 1e-3)
check("supertype conf >0.8", np.mean(o.supertype_confidence > 0.8), 0.431, 1e-3)
check("subclass agree prop-vs-capped",
      float((o.subclass.astype(str).values == c.subclass.astype(str).values).mean()), 0.956, 1e-3)
check("supertype agree prop-vs-capped",
      float((o.supertype.astype(str).values == c.supertype.astype(str).values).mean()), 0.677, 1e-3)
G = ["Sst", "Pvalb", "Vip", "Lamp5", "Sncg", "Chandelier", "Pax6", "Sst Chodl", "Lamp5 Lhx6"]
E = ["L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5/6 NP", "L6 CT", "L6b"]
gab = lambda s, col: 100 * s[col].isin(G).sum() / (s[col].isin(G).sum() + s[col].isin(E).sum())
check("GABA %% of neurons (proportional)", gab(o, "subclass"), 32.3, 0.05)
check("GABA %% of neurons (capped)", gab(c, "subclass"), 37.3, 0.05)

print("\n--- DE vs the meta-analysis ---")
mk = pd.read_csv(f"{OUTD}/compare_markers.csv").set_index("gene")
check("SST meta logFC", mk.loc["SST", "meta_logFC"], -0.458, 5e-4)
check("SST meta FDR", mk.loc["SST", "meta_FDR"], 0.049, 5e-4)
check("SST gse logFC", mk.loc["SST", "gse_logFC"], -0.175, 5e-4)
check("SST gse p", mk.loc["SST", "gse_P"], 0.61, 5e-3)
check("PVALB meta logFC", mk.loc["PVALB", "meta_logFC"], -0.056, 5e-4)
check("PVALB gse logFC", mk.loc["PVALB", "gse_logFC"], -0.021, 5e-4)
check("PVALB gse p", mk.loc["PVALB", "gse_P"], 0.94, 5e-3)

j = pd.read_csv(f"{OUTD}/compare_de_joined.csv.gz")
s = j[j.meta_FDR < 0.10]
check("n meta-sig pairs", len(s), 6864, 0)
r, p = stats.pearsonr(s.meta_logFC, s.logFC)
check("DE concordance r", r, 0.376, 5e-4)
same = int(np.sum(np.sign(s.meta_logFC) == np.sign(s.logFC)))
check("DE sign-concordance", same / len(s), 0.686, 1e-3)
check("DE sign-test P (log10)", np.log10(stats.binomtest(same, len(s), 0.5).pvalue), np.log10(3.2e-214), 0.05)
nt = pd.read_csv(f"{OUTD}/compare_named_genes.csv")
check("named genes concordant", int(nt.same_sign.sum()), 15, 0)
check("named genes total", len(nt), 16, 0)

print("\n--- composition ---")
h = pd.read_csv(f"{OUTD}/compare_headline_supertypes.csv")
check("headline supertypes concordant", int(h.same_sign.sum()), 8, 0)
for st, want in [("Sst_2", -0.364), ("Sst_20", -0.347), ("Sst_25", -0.216)]:
    check(f"{st} gse logFC (all 26)", float(h[h.supertype == st].gse_logFC.iloc[0]), want, 5e-4)
u70 = pd.read_csv(f"{OUTD}/under70_composition_crumblr.csv")
u70 = u70[u70.level == "supertype"].set_index("cell_type")
for st, want in [("Sst_2", -0.370), ("Sst_20", -0.338), ("Sst_25", -0.226)]:
    check(f"{st} gse logFC (<=70)", float(u70.loc[st, "logFC"]), want, 5e-4)

for level, want_r, want_p, want_n in [("subclass", 0.702, 8.0e-4, 19),
                                      ("supertype", 0.367, 8.0e-5, 110)]:
    x = pd.read_csv(f"{OUTD}/compare_composition_{level}.csv")
    n = x[x.compartment == "neuronal"]
    r, p = stats.pearsonr(n.meta_beta, n.logFC)
    check(f"{level} neuronal n", len(n), want_n, 0)
    check(f"{level} neuronal r", r, want_r, 5e-4)
    check(f"{level} neuronal P (log10)", np.log10(p), np.log10(want_p), 0.02)

comp = pd.read_csv(f"{OUTD}/composition_crumblr.csv")
sst = comp[(comp.level == "subclass") & (comp.cell_type == "Sst")].iloc[0]
check("subclass Sst logFC", sst.logFC, -0.010, 5e-4)
check("subclass Sst p", sst["P.Value"], 0.92, 5e-3)

print("\n--- Reiner reproduction ---")
rr = pd.read_csv(f"{OUTD}/reproduce_reiner.csv")
check("pooled concordant", int(rr.sign_concordant.sum()), 4647, 0)
check("pooled pairs", int(rr.n_genes.sum()), 4863, 0)
check("pooled fraction", rr.sign_concordant.sum() / rr.n_genes.sum(), 0.956, 1e-3)
sr = rr[rr.our_subclass == "Sst"].iloc[0]
check("Sst genes", int(sr.n_genes), 29, 0)
check("Sst concordant", int(sr.sign_concordant), 29, 0)
check("Sst spearman", sr.spearman, 0.81, 5e-3)
check("Sst median our p", sr.median_our_p, 0.0022, 5e-5)
check("Sst median bg p", sr.median_bg_p, 0.418, 5e-4)

print("\n--- under-70 DE ---")
u = pd.read_csv(f"{OUTD}/under70_de_subclass.csv.gz")
check("SST in Sst (<=70) logFC", float(u[(u.gene == "SST") & (u.cell_type == "Sst")].logFC.iloc[0]),
      -0.110, 5e-4)

print(f"\n{'=' * 72}\n{ok} passed, {fail} failed")
raise SystemExit(1 if fail else 0)
