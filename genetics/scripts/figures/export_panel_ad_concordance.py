#!/usr/bin/env python3
"""
SCZ vs Alzheimer's compositional concordance across the 16 Sst supertypes
(Figure 4i; also the SCZ axis of Supplementary Fig. S10b).

Extracted from export_fig4_new_panels.py. That script's other sections built
the volcano and CALB1 violin from the SEA-AD *MTG* reference; those panels were
rebuilt on the region-matched DLPFC (A9) reference by export_panels_dehi.py on
2026-08-12. Because the old script was one monolithic main(), re-running it to
refresh this panel would silently regress the volcano and violin to MTG. This
file does the one live piece and nothing else.

Inputs : shared/snrnaseq_de/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv
         crossdisorder/results/crumblr_results_supertype_neurons.csv
           (SEA-AD DLPFC pseudo-progression slopes; crossdisorder/README.md)
Outputs: r_panels/panel_ad_concordance_sst.csv
Read by: fig4_new_panels.R build_ad_concordance(); export_fig4_robustness.py
"""
import numpy as np
import pandas as pd
from scipy import stats

REPO = "/Users/shreejoy/Github/scz_celltype_paper"
SCZ_CRUMBLR = (f"{REPO}/shared/snrnaseq_de/nicole_scz_snrnaseq_betas"
               "/final_results_crumblr_7_cohorts.csv")
AD_CRUMBLR = f"{REPO}/crossdisorder/results/crumblr_results_supertype_neurons.csv"
OUT = f"{REPO}/genetics/results/figures/r_panels"

# The five depleted supertypes of Fig. 3 (FDR < 0.20), named exactly as in the
# compositional analysis so the two figures group cells the same way.
VULNERABLE = ["Sst_25", "Sst_22", "Sst_2", "Sst_20", "Sst_3"]

scz = (pd.read_csv(SCZ_CRUMBLR)
         .rename(columns={"CellType": "supertype", "estimate": "scz_beta",
                          "se": "scz_se", "padj": "scz_fdr"})
         [["supertype", "scz_beta", "scz_se", "scz_fdr"]])
adc = (pd.read_csv(AD_CRUMBLR)
         .rename(columns={"celltype": "supertype", "logFC": "ad_slope",
                          "SE": "ad_se", "FDR": "ad_fdr"})
         [["supertype", "ad_slope", "ad_se", "ad_fdr"]])

m = scz.merge(adc, on="supertype")
m = m[m.supertype.str.startswith("Sst_")].copy()
m["group"] = np.where(m.supertype.isin(VULNERABLE), "Vulnerable", "Not depleted")
m.to_csv(f"{OUT}/panel_ad_concordance_sst.csv", index=False)

r, rp = stats.pearsonr(m.scz_beta, m.ad_slope)
rho, sp = stats.spearmanr(m.scz_beta, m.ad_slope)
conc = float((np.sign(m.scz_beta) == np.sign(m.ad_slope)).mean() * 100)
print(f"panel i: n={len(m)} Sst supertypes | r={r:.3f} (p={rp:.2g}), "
      f"rho={rho:.3f} (p={sp:.2g}), {conc:.1f}% sign-concordant")
