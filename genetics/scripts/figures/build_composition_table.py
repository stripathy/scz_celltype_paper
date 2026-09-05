#!/usr/bin/env python3
"""
Per-supertype SCZ compositional effects, reshaped for the Figure 4 exporters.

Replaces scripts/13_gwas_vs_composition.py, which produced the same table by
inner-joining the crumblr meta-analysis against the retired 501-type combined
enrichment. That join never dropped a row (all 109 crumblr supertypes matched),
and every enrichment-derived column it added was overwritten downstream by
export_panels_abc.py from the SEA-AD-only MAGMA run. So the join contributed
nothing but a dependency on the retired taxonomy.

Input : shared/snrnaseq_de/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv
        (Endresz et al., 7-cohort crumblr meta-analysis; symlinked, not committed)
Output: genetics/results/tables/gwas_vs_casecontrol_composition.csv
Read by: export_panels_abc.py (panel a)
"""
import numpy as np
import pandas as pd

import os as _os

# Repo-relative, so the chain runs from any clone. External data that is not in
# the repo is still resolved by absolute path or an environment override below.
GEN = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
PAPER = _os.path.dirname(GEN)


CRUMBLR = (f"{GEN}/../shared/snrnaseq_de/nicole_scz_snrnaseq_betas"
           "/final_results_crumblr_7_cohorts.csv")
OUT = f"{GEN}/results/tables/gwas_vs_casecontrol_composition.csv"

comp = pd.read_csv(CRUMBLR)
out = pd.DataFrame({
    "supertype": comp.CellType,
    "comp_beta": comp.estimate,
    "se_y": comp.se,
    "comp_zval": comp.zval,
    "comp_pval": comp.pval,
    "comp_padj": comp.padj,
})
out["abs_comp_beta"] = out.comp_beta.abs()
out["comp_logp"] = -np.log10(np.clip(out.comp_pval, 1e-300, 1))
out["is_sst"] = out.supertype.str.startswith("Sst")
out.to_csv(OUT, index=False)

n_sst = int((out.supertype.str.startswith("Sst")
             & ~out.supertype.str.startswith("Sst Chodl")).sum())
print(f"wrote {OUT}\n  {len(out)} supertypes ({n_sst} Sst)")
