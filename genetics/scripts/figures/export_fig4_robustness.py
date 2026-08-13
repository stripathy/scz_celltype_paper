"""
Assemble the two comparison tables behind the Figure 4 robustness supplement.

Figure 4a and 4i each rest on a discretionary choice. Panel a depends on which
SCZ GWAS is used and which cortical region supplies the normotypic expression;
panel i depends on which region the AD compositional signal is measured in.
This collects the re-runs so the supplement can show all of them side by side.

Every re-run holds the taxonomy, the specificity recipe, the MAGMA settings and
the compositional model fixed. Only the named input changes -- in particular the
MTG enrichment reuses the published 501-type definition and RBH removals, so it
isolates the reference region rather than confounding it with a taxonomy change.

Inputs are produced by:
    build_spec_mtg_sametax.py            MTG expression on the published taxonomy
    magma --gene-covar                   the three .gsa.out files
    crossdisorder/code/build_mtg_cps_input.py + run_mtg_cps_crumblr.R
"""
import numpy as np
import pandas as pd
from scipy import stats

GEN = "/Users/shreejoy/Github/scz_celltype_paper/genetics"
XD = "/Users/shreejoy/Github/scz_celltype_paper/crossdisorder"
I = f"{GEN}/results/intermediates"
PANELS = f"{GEN}/results/figures/r_panels"
TABLES = f"{GEN}/results/tables"

RUNS = {
    "A9 + Bigdeli (published)": (f"{I}/T_a9rbh_bigdeli.gsa.out", f"{I}/namemap_a9rbh.csv"),
    "A9 + PGC3":                (f"{I}/T_a9_pgc3.gsa.out",       f"{I}/namemap_a9rbh.csv"),
    "MTG + Bigdeli":            (f"{I}/T_mtg_bigdeli.gsa.out",   f"{I}/namemap_mtg_sametax.csv"),
}
IS_SST = r"Sst_\d+$"


def enrichment(gsa, namemap):
    names = pd.read_csv(namemap).set_index("safe_name").cell_type
    g = pd.read_csv(gsa, sep=r"\s+", comment="#")
    g["supertype"] = g.VARIABLE.map(names)
    return g.dropna(subset=["supertype"]).set_index("supertype").P


# ---- panel a: enrichment vs SCZ depletion, under each configuration --------
comp = pd.read_csv(f"{PANELS}/panel_B_genetics_vs_depletion.csv")[
    ["supertype", "comp_beta", "color", "depleted_fdr20"]]
out = comp[comp.supertype.str.match(IS_SST)].copy()
for label, (gsa, nm) in RUNS.items():
    out[label] = out.supertype.map(-np.log10(enrichment(gsa, nm)))
    r = stats.spearmanr(out[label], -out.comp_beta)
    print(f"  {label:<26s} rho={r.statistic:+.3f}  p={r.pvalue:.4f}")
out.to_csv(f"{TABLES}/fig4_robustness_sst16.csv", index=False)

cols = list(RUNS)
print("\n  pairwise agreement of the enrichment vectors:")
for i in range(len(cols)):
    for j in range(i + 1, len(cols)):
        print(f"    {cols[i]} vs {cols[j]}: "
              f"{stats.spearmanr(out[cols[i]], out[cols[j]]).statistic:.3f}")

# ---- panel b: SCZ depletion vs AD depletion measured in MTG ----------------
pub = pd.read_csv(f"{PANELS}/panel_ad_concordance_sst.csv")
mtg = pd.read_csv(f"{XD}/results/crumblr_results_mtg_supertype_neurons.csv")
mtg = mtg[mtg.celltype.str.match(IS_SST)][["celltype", "logFC", "SE", "FDR"]]
mtg.columns = ["supertype", "ad_slope_mtg", "ad_se_mtg", "ad_fdr_mtg"]
ad = pub.merge(mtg, on="supertype")
ad.to_csv(f"{TABLES}/fig4_robustness_ad_mtg.csv", index=False)

for label, col in (("AD DLPFC (published)", "ad_slope"), ("AD MTG", "ad_slope_mtg")):
    r = stats.spearmanr(-ad.scz_beta, -ad[col])
    conc = 100 * np.mean(np.sign(-ad.scz_beta) == np.sign(-ad[col]))
    print(f"\n  SCZ vs {label:<22s} rho={r.statistic:+.3f}  p={r.pvalue:.2g}  "
          f"direction agreement {conc:.1f}%")
print(f"  AD DLPFC vs AD MTG slopes: "
      f"rho={stats.spearmanr(ad.ad_slope, ad.ad_slope_mtg).statistic:.3f}")
