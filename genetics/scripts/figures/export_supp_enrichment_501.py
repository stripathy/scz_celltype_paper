"""
Table behind the supplementary "enrichment landscape" figure: SCZ common-variant
enrichment (MAGMA gene-property, Bigdeli 2026 GWAS, SEA-AD DLPFC expression
reference) for every one of the 501 types in the combined SEA-AD + Siletti
taxonomy used in Figure 4, annotated with source (SEA-AD supertype vs novel
Siletti cluster), grouping (SEA-AD subclass or Siletti supercluster) and the
SEA-AD supertype palette used in Figs 1b/3a/4.

Input:  results/intermediates/T_a9rbh_bigdeli.gsa.out (+ namemap_a9rbh.csv)
        data/seaad_supertype_colors.json
Output: results/tables/scz_enrichment_501_bigdeli_dlpfc.csv
"""
import json
import re

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

GEN = "/Users/shreejoy/Github/scz_celltype_paper/genetics"
GSA = f"{GEN}/results/intermediates/T_a9rbh_bigdeli.gsa.out"
NAMEMAP = f"{GEN}/results/intermediates/namemap_a9rbh.csv"
COLORS = f"{GEN}/data/seaad_supertype_colors.json"
OUT = f"{GEN}/results/tables/scz_enrichment_501_bigdeli_dlpfc.csv"

# SEA-AD subclass order (Fig. 3a order for neurons, then non-neuronal)
SEAAD_ORDER = ["Lamp5_Lhx6", "Lamp5", "Pax6", "Sncg", "Vip", "Sst Chodl", "Sst",
               "Pvalb", "Chandelier", "L2/3 IT", "L4 IT", "L5 IT", "L6 IT",
               "L6 IT Car3", "L5 ET", "L6 CT", "L6b", "L5/6 NP",
               "Astro", "Oligo", "OPC", "Micro-PVM", "Endo", "VLMC"]
CLASS = {**{s: "Inhibitory" for s in SEAAD_ORDER[:9]},
         **{s: "Excitatory" for s in SEAAD_ORDER[9:18]},
         **{s: "Non-neuronal" for s in SEAAD_ORDER[18:]}}
# Siletti cluster prefix -> supercluster (Siletti et al. 2023 naming)
SILETTI = {
    "Amex": "Amygdala excitatory", "Astro": "Astrocyte", "Bgl": "Bergmann glia",
    "CA13": "Hippocampal CA1-3", "CA4": "Hippocampal CA4", "DG": "Hippocampal dentate gyrus",
    "CBI": "Cerebellar inhibitory", "CGE": "CGE interneuron", "MGE": "MGE interneuron",
    "LLC": "LAMP5-LHX6 and Chandelier", "COP": "Committed oligodendrocyte precursor",
    "Chrp": "Choroid plexus", "DLCT6b": "Deep-layer corticothalamic and 6b",
    "DLIT": "Deep-layer intratelencephalic", "DLNP": "Deep-layer near-projecting",
    "ULIT": "Upper-layer intratelencephalic", "EMSN": "Eccentric medium spiny neuron",
    "MSN": "Medium spiny neuron", "Epen": "Ependymal", "Fbl": "Fibroblast",
    "LRL": "Lower rhombic lip", "URL": "Upper rhombic lip", "Mgl": "Microglia",
    "Midi": "Midbrain-derived inhibitory", "Mmb": "Mammillary body", "Misc": "Miscellaneous",
    "L5ET": "Miscellaneous", "Bcell": "Miscellaneous", "Tcell": "Miscellaneous",
    "Nkcell": "Miscellaneous", "Mono": "Miscellaneous", "OPC": "Oligodendrocyte precursor",
    "Oligo": "Oligodendrocyte", "Per": "Vascular", "Vsmc": "Vascular", "VendA": "Vascular",
    "VendAC": "Vascular", "VendV": "Vascular", "VendVC": "Vascular", "VendPLVAP": "Vascular",
    "Splat": "Splatter", "Thex": "Thalamic excitatory",
}
DEPLETED = ["Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25"]   # FDR < 0.20 in Fig. 3a


def seaad_subclass(name: str) -> str:
    base = re.sub(r"(-SEAAD)$", "", name)
    base = re.sub(r"(_\d+)+$", "", base)
    return base


def main() -> None:
    names = pd.read_csv(NAMEMAP).set_index("safe_name").cell_type
    g = pd.read_csv(GSA, sep=r"\s+", comment="#")
    g["cell_type"] = g.VARIABLE.map(names)
    g = g.dropna(subset=["cell_type"]).copy()
    assert len(g) == 501, len(g)
    colors = json.load(open(COLORS))

    g["source"] = np.where(g.cell_type.isin(colors), "SEA-AD", "Siletti")
    g["group"] = [seaad_subclass(n) if s == "SEA-AD" else SILETTI[re.sub(r"_\d+$", "", n)]
                  for n, s in zip(g.cell_type, g.source)]
    bad = g[(g.source == "SEA-AD") & ~g.group.isin(SEAAD_ORDER)]
    assert bad.empty, bad.cell_type.tolist()
    g["cell_class"] = [CLASS.get(gr, "Siletti") if s == "SEA-AD" else "Siletti"
                       for gr, s in zip(g.group, g.source)]
    g["color"] = [colors.get(n, "#9e9e9e") for n in g.cell_type]
    g["neglog10p"] = -np.log10(g.P)
    g["fdr"] = multipletests(g.P, method="fdr_bh")[1]
    g["bonferroni"] = np.minimum(g.P * len(g), 1)
    g["depleted_scz"] = g.cell_type.isin(DEPLETED)
    out = g[["cell_type", "source", "group", "cell_class", "color", "NGENES", "BETA", "SE",
             "P", "neglog10p", "fdr", "bonferroni", "depleted_scz"]].rename(
        columns={"NGENES": "n_genes", "BETA": "beta", "SE": "se", "P": "p"})
    out.to_csv(OUT, index=False)
    thr_bonf = -np.log10(0.05 / len(g))
    fdr_ok = out.p[out.fdr < 0.05]
    print(f"{len(out)} types: {(out.source == 'SEA-AD').sum()} SEA-AD, {(out.source == 'Siletti').sum()} Siletti")
    print(f"Bonferroni threshold -log10P = {thr_bonf:.2f}; n Bonferroni-sig = {(out.bonferroni < 0.05).sum()}; "
          f"n FDR<0.05 = {(out.fdr < 0.05).sum()} (max P at FDR<0.05 = {fdr_ok.max():.2e})")
    print(out.sort_values('p').head(15)[['cell_type', 'source', 'group', 'neglog10p']].to_string(index=False))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
