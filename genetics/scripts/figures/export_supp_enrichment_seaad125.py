"""
Table behind the supplementary "enrichment landscape" figure: SCZ common-variant
enrichment (MAGMA gene-property, Bigdeli 2026 GWAS, SEA-AD DLPFC expression
reference) for every one of the 125 SEA-AD DLPFC supertypes, annotated with
subclass grouping and the SEA-AD supertype palette used in Figs 1b/3a/4.

Supersedes export_supp_enrichment_501.py: on L. Duncan's advice the combined
SEA-AD + Siletti taxonomy was retired, so the landscape -- like Fig. 4a -- is
computed on the SEA-AD taxonomy alone (build_spec_seaad_only.py) and corrected
across the 125 supertypes actually tested.

Input:  results/intermediates/T_a9only_bigdeli.gsa.out (+ namemap_a9only.csv)
        data/seaad_supertype_colors.json
Output: results/tables/scz_enrichment_seaad125_bigdeli_dlpfc.csv
"""
import json
import re

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

import os as _os

# Repo-relative, so the chain runs from any clone. External data that is not in
# the repo is still resolved by absolute path or an environment override below.
GEN = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
PAPER = _os.path.dirname(GEN)


GSA = f"{GEN}/results/intermediates/T_a9only_bigdeli.gsa.out"
NAMEMAP = f"{GEN}/results/intermediates/namemap_a9only.csv"
COLORS = f"{GEN}/data/seaad_supertype_colors.json"
OUT = f"{GEN}/results/tables/scz_enrichment_seaad125_bigdeli_dlpfc.csv"

# SEA-AD subclass order (Fig. 3a order for neurons, then non-neuronal)
SEAAD_ORDER = ["Lamp5_Lhx6", "Lamp5", "Pax6", "Sncg", "Vip", "Sst Chodl", "Sst",
               "Pvalb", "Chandelier", "L2/3 IT", "L4 IT", "L5 IT", "L6 IT",
               "L6 IT Car3", "L5 ET", "L6 CT", "L6b", "L5/6 NP",
               "Astro", "Oligo", "OPC", "Micro-PVM", "Endo", "VLMC"]
CLASS = {**{s: "Inhibitory" for s in SEAAD_ORDER[:9]},
         **{s: "Excitatory" for s in SEAAD_ORDER[9:18]},
         **{s: "Non-neuronal" for s in SEAAD_ORDER[18:]}}
DEPLETED = ["Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25"]   # FDR < 0.20 in Fig. 3a


def seaad_subclass(name: str) -> str:
    base = re.sub(r"(-SEAAD)$", "", name)
    return re.sub(r"(_\d+)+$", "", base)


def main() -> None:
    names = pd.read_csv(NAMEMAP).set_index("safe_name").cell_type
    g = pd.read_csv(GSA, sep=r"\s+", comment="#")
    g["cell_type"] = g.VARIABLE.map(names)
    g = g.dropna(subset=["cell_type"]).copy()
    assert len(g) == 125, len(g)
    colors = json.load(open(COLORS))

    g["group"] = [seaad_subclass(n) for n in g.cell_type]
    bad = g[~g.group.isin(SEAAD_ORDER)]
    assert bad.empty, bad.cell_type.tolist()
    g["cell_class"] = g.group.map(CLASS)
    g["color"] = [colors.get(n, "#9e9e9e") for n in g.cell_type]
    g["neglog10p"] = -np.log10(g.P)
    g["fdr"] = multipletests(g.P, method="fdr_bh")[1]
    g["bonferroni"] = np.minimum(g.P * len(g), 1)
    g["depleted_scz"] = g.cell_type.isin(DEPLETED)
    out = g[["cell_type", "group", "cell_class", "color", "NGENES", "BETA", "SE",
             "P", "neglog10p", "fdr", "bonferroni", "depleted_scz"]].rename(
        columns={"NGENES": "n_genes", "BETA": "beta", "SE": "se", "P": "p"})
    out.to_csv(OUT, index=False)
    thr_bonf = -np.log10(0.05 / len(g))
    fdr_ok = out.p[out.fdr < 0.05]
    sig = out[out.fdr < 0.05]
    print(f"{len(out)} supertypes; Bonferroni threshold -log10P = {thr_bonf:.2f}; "
          f"n Bonferroni-sig = {(out.bonferroni < 0.05).sum()}; "
          f"n FDR<0.05 = {len(sig)} (max P at FDR<0.05 = {fdr_ok.max():.2e})")
    print(f"inhibitory among FDR-sig: {(sig.cell_class == 'Inhibitory').sum()}/{len(sig)} "
          f"({100 * (sig.cell_class == 'Inhibitory').mean():.0f}%)")
    print(out.sort_values('p').head(15)[['cell_type', 'group', 'neglog10p']].to_string(index=False))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
