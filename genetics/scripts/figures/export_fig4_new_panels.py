#!/usr/bin/env python3
"""
Extract panel data for the two new Figure 4 elements.

Panel 1 (volcano + violin): marker genes distinguishing SCZ-vulnerable Sst
supertypes from non-depleted Sst supertypes, in the SEA-AD *neurotypical*
reference (i.e. normative cell-identity differences, not disease DE).

  Cell-level Wilcoxon defines the volcano (standard marker-DE; the reference
  has only 5 donors, so a donor-paired test cannot reach genome-wide FDR).
  A donor-paired t-test (n = 5) is computed alongside so labelled genes can be
  reported with donor-level support -- the design used for the CALB1 violin.

Panel 2 (AD concordance): SCZ case-control crumblr beta vs the SEA-AD
Alzheimer's pseudo-progression (CPS) crumblr slope, across the 16 Sst supertypes.

Outputs -> genetics/results/figures/r_panels/
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats, sparse

REPO = "/Users/shreejoy/Github/scz_celltype_paper"
XEN = "/Users/shreejoy/Github/SCZ_Xenium"
# Canonical normative-expression source for ALL of Figure 4: the same SEA-AD
# neurotypical MTG reference (raw counts) that scripts/01_compute_specificity.py
# uses for MAGMA specificity and the per-supertype mean expression in panels d/e.
REF_H5AD = "/Users/shreejoy/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad"
SCZ_CRUMBLR = f"{REPO}/spatial/data/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv"
# AD axis = SEA-AD DLPFC (A9), region-matched to the frontal-cortex SCZ cohorts.
# Built in-repo by crossdisorder/code/build_seaad_cps_input.py +
# run_seaad_cps_crumblr.R; see crossdisorder/README.md. (Was the MTG run under
# SCZ_Xenium.)
AD_CRUMBLR = f"{REPO}/crossdisorder/results/crumblr_results_supertype_neurons.csv"
OUT = f"{REPO}/genetics/results/figures/r_panels"

# Group definitions -- match the compositional analysis (Fig 3) exactly.
VULNERABLE = ["Sst_25", "Sst_22", "Sst_2", "Sst_20", "Sst_3"]
NOT_DEPLETED = ["Sst_19", "Sst_9", "Sst_23", "Sst_11", "Sst_13",
                "Sst_1", "Sst_4", "Sst_5", "Sst_7", "Sst_10", "Sst_12"]

# Genes to guarantee a label on the volcano (story-critical), beyond top hits.
ALWAYS_LABEL = ["CALB1", "HCN1", "SST", "CALB2", "NPY", "RELN", "NOS1"]

os.makedirs(OUT, exist_ok=True)


def main():
    print("Loading SEA-AD neurotypical reference (backed)...")
    ad = sc.read_h5ad(REF_H5AD, backed="r")

    keep = ad.obs["Supertype"].astype(str).isin(VULNERABLE + NOT_DEPLETED).values
    print(f"  Sst cells: {keep.sum():,} of {ad.n_obs:,}")
    sub = ad[keep].to_memory()
    ad.file.close()

    sub.obs["group"] = np.where(
        sub.obs["Supertype"].astype(str).isin(VULNERABLE), "Vulnerable", "Not depleted")
    sub.obs["donor"] = sub.obs["donor_id"].astype(str)
    print(sub.obs.groupby(["group"], observed=True).size().to_string())
    print(f"  donors: {sub.obs['donor'].nunique()}")

    # Normalise exactly as the genetics pipeline does (and as Seurat
    # NormalizeData does): CP10K + log1p. Verified to reproduce the pipeline's
    # stored per-supertype means (HCN1 in Sst_25 = 2.4404).
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    X = sub.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    X = X.tocsc()
    genes = np.asarray(sub.var_names)
    is_v = (sub.obs["group"] == "Vulnerable").values
    donors = sub.obs["donor"].values
    uniq_donors = np.unique(donors)

    # --- expression filter: detected in >=10% of cells in either group ---
    det = np.asarray((X > 0).sum(axis=0)).ravel()
    det_v = np.asarray((X[is_v] > 0).sum(axis=0)).ravel() / max(is_v.sum(), 1)
    det_n = np.asarray((X[~is_v] > 0).sum(axis=0)).ravel() / max((~is_v).sum(), 1)
    expressed = (det_v >= 0.10) | (det_n >= 0.10)
    print(f"  genes passing 10% detection: {expressed.sum():,} of {len(genes):,}")

    idx = np.where(expressed)[0]
    Xe = X[:, idx]
    genes_e = genes[idx]

    # --- cell-level Wilcoxon (volcano) ---
    print("Cell-level Wilcoxon DE...")
    dense = np.asarray(Xe.todense())
    stat, pval = stats.mannwhitneyu(dense[is_v], dense[~is_v], axis=0, alternative="two-sided")

    # Seurat-style avg_log2FC on the CP10K (un-logged) scale
    mv = np.expm1(dense[is_v]).mean(axis=0)
    mn = np.expm1(dense[~is_v]).mean(axis=0)
    log2fc = np.log2((mv + 1.0) / (mn + 1.0))

    # BH FDR
    order = np.argsort(pval)
    ranked = pval[order]
    n = len(ranked)
    fdr_sorted = ranked * n / (np.arange(n) + 1)
    fdr_sorted = np.minimum.accumulate(fdr_sorted[::-1])[::-1]
    fdr = np.empty(n)
    fdr[order] = np.clip(fdr_sorted, 0, 1)

    # --- donor-level paired t-test (n = 5), same design as the violin ---
    print("Donor-level paired t-test (n=5 donors)...")
    dmean_v = np.vstack([dense[is_v & (donors == d)].mean(axis=0) for d in uniq_donors])
    dmean_n = np.vstack([dense[(~is_v) & (donors == d)].mean(axis=0) for d in uniq_donors])
    t_d, p_d = stats.ttest_rel(dmean_v, dmean_n, axis=0)

    # BH FDR on the donor-paired p-values (the threshold used to call genes in
    # the volcano and to annotate the violin).
    order_d = np.argsort(p_d)
    ranked_d = p_d[order_d]
    fdr_d_sorted = ranked_d * n / (np.arange(n) + 1)
    fdr_d_sorted = np.minimum.accumulate(fdr_d_sorted[::-1])[::-1]
    fdr_donor = np.empty(n)
    fdr_donor[order_d] = np.clip(fdr_d_sorted, 0, 1)

    res = pd.DataFrame({
        "gene": genes_e,
        "log2FC": log2fc,
        "pct_vulnerable": det_v[idx],
        "pct_not_depleted": det_n[idx],
        "mean_vulnerable": dense[is_v].mean(axis=0),
        "mean_not_depleted": dense[~is_v].mean(axis=0),
        "p_wilcoxon": pval,
        "fdr_wilcoxon": fdr,
        "donor_delta": (dmean_v - dmean_n).mean(axis=0),
        "p_donor_paired": p_d,
        "fdr_donor": fdr_donor,
    }).sort_values("p_wilcoxon")
    res["neg_log10_p"] = -np.log10(np.clip(res["p_wilcoxon"], 1e-300, None))
    res["always_label"] = res["gene"].isin(ALWAYS_LABEL)
    res.to_csv(f"{OUT}/panel_volcano_vulnerable_vs_notdepleted.csv", index=False)
    print(f"  wrote volcano table: {len(res):,} genes")

    for g in ALWAYS_LABEL:
        r = res[res["gene"] == g]
        if len(r):
            r = r.iloc[0]
            print(f"    {g:7s} log2FC={r['log2FC']:+.3f}  p_wilcox={r['p_wilcoxon']:.2e}  "
                  f"FDR={r['fdr_wilcoxon']:.2e}  donor-paired p={r['p_donor_paired']:.4f}")

    # --- per-cell expression for the CALB1 violin (+ SST for reference) ---
    viol = []
    for g in ["CALB1", "SST"]:
        if g not in set(genes):
            continue
        gi = int(np.where(genes == g)[0][0])
        vals = np.asarray(X[:, gi].todense()).ravel()
        viol.append(pd.DataFrame({
            "gene": g, "expression": vals,
            "group": sub.obs["group"].values, "donor": donors,
            "supertype": sub.obs["Supertype"].astype(str).values}))
    viol = pd.concat(viol, ignore_index=True)
    viol.to_csv(f"{OUT}/panel_violin_calb1_percell.csv", index=False)

    # donor-level means + paired p, for the violin annotation
    dsum = (viol.groupby(["gene", "donor", "group"], observed=True)["expression"]
                .mean().reset_index()
                .pivot_table(index=["gene", "donor"], columns="group", values="expression")
                .reset_index())
    fdr_lookup = dict(zip(res["gene"], res["fdr_donor"]))
    stats_rows = []
    for g, gg in dsum.groupby("gene"):
        gg = gg.dropna(subset=["Vulnerable", "Not depleted"])
        tt = stats.ttest_rel(gg["Vulnerable"], gg["Not depleted"])
        stats_rows.append({"gene": g, "n_donors": len(gg),
                           "mean_vulnerable": gg["Vulnerable"].mean(),
                           "mean_not_depleted": gg["Not depleted"].mean(),
                           "t": tt.statistic, "p_donor_paired": tt.pvalue,
                           "fdr_donor": fdr_lookup.get(g, float("nan"))})
        print(f"  violin {g}: n={len(gg)} donors, paired p={tt.pvalue:.4f}, "
              f"FDR={fdr_lookup.get(g, float('nan')):.4f}")
    dsum.to_csv(f"{OUT}/panel_violin_donor_means.csv", index=False)
    pd.DataFrame(stats_rows).to_csv(f"{OUT}/panel_violin_stats.csv", index=False)

    # --- AD concordance across the 16 Sst supertypes ---
    print("Building AD concordance panel data...")
    scz = (pd.read_csv(SCZ_CRUMBLR)
             .rename(columns={"CellType": "supertype", "estimate": "scz_beta",
                              "se": "scz_se", "padj": "scz_fdr"})
             [["supertype", "scz_beta", "scz_se", "scz_fdr"]])
    adc = (pd.read_csv(AD_CRUMBLR)
             .rename(columns={"celltype": "supertype", "logFC": "ad_slope",
                              "SE": "ad_se", "FDR": "ad_fdr"})
             [["supertype", "ad_slope", "ad_se", "ad_fdr"]])
    m = scz.merge(adc, on="supertype")
    m = m[m["supertype"].str.startswith("Sst_")].copy()
    m["group"] = np.where(m["supertype"].isin(VULNERABLE), "Vulnerable", "Not depleted")

    r, rp = stats.pearsonr(m["scz_beta"], m["ad_slope"])
    rho, sp = stats.spearmanr(m["scz_beta"], m["ad_slope"])
    conc = float((np.sign(m["scz_beta"]) == np.sign(m["ad_slope"])).mean() * 100)
    m.to_csv(f"{OUT}/panel_ad_concordance_sst.csv", index=False)
    pd.DataFrame([{"n": len(m), "pearson_r": r, "pearson_p": rp,
                   "spearman_rho": rho, "spearman_p": sp,
                   "pct_sign_concordant": conc}]).to_csv(
        f"{OUT}/panel_ad_concordance_stats.csv", index=False)
    print(f"  Sst n={len(m)}: r={r:.3f} (p={rp:.2g}), rho={rho:.3f}, {conc:.0f}% sign-concordant")
    record_provenance()
    print("\nDone.")


def record_provenance():
    """Record which upstream file each panel CSV came from, so the R renderer
    stops rather than drawing stale numbers. Source map (covering both
    exporters that write into r_panels/): r_panels_provenance.py."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from r_panels_provenance import record
    record(OUT)


if __name__ == "__main__":
    main()
