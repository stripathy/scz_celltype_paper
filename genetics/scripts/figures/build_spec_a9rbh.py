"""
Specificity matrix using the A9-derived RBH taxonomy.

Identical to build_dlpfc_specificity.py from the merge step onward; the only
change is that redundancy between SEA-AD and Siletti is decided from the A9
expression matrix rather than MTG, so every step of the figure uses one
reference. The A9 pooling step is not repeated -- a9_supertype_log1p_mean.csv
is exactly what that code path wrote.
"""
import numpy as np
import pandas as pd
import h5py

W = "/Users/shreejoy/Github/scz_celltype_paper/genetics/results/intermediates"
SRC = "/Users/shreejoy/Github/scz_cell_type_enrichment"
GLOC = f"{SRC}/linking_cell_types_to_brain_phenotypes/Data/NCBI37.3.gene.loc.extendedMHCexcluded"

sea = pd.read_csv(f"{W}/a9_supertype_log1p_mean.csv", index_col=0)
print(f"A9 reference: {sea.shape[0]:,} genes x {sea.shape[1]} supertypes")

# ---------------- merge with Siletti, Duncan recipe ------------------------
ens_map = pd.read_csv(f"{W}/map_entrez_ensembl.csv", dtype=str)
sym_map = pd.read_csv(f"{W}/map_entrez_symbol.csv", dtype=str)
s2e = sym_map.merge(ens_map, on="ENTREZ")[["SYMBOL", "ENSEMBL"]]
s2e = s2e[~s2e.SYMBOL.duplicated(keep=False) & ~s2e.ENSEMBL.duplicated(keep=False)]
sea = sea[~sea.index.duplicated(keep=False)].join(
    s2e.set_index("SYMBOL"), how="inner").set_index("ENSEMBL")
print(f"SEA-AD DLPFC in ENSEMBL space: {sea.shape}")

with h5py.File(f"{W}/Siletti_L2-cluster-log1p_matrix.h5", "r") as f:
    mat, acc = f["matrix"][:], np.array([a.decode().split(".")[0] for a in f["Accession"][:]])
names = pd.read_csv(f"{W}/siletti_cluster_names.csv")
n2n = dict(zip(names.cluster_num, names.cluster_name))
sil = pd.DataFrame(mat, index=acc, columns=[n2n.get(i, f"Cluster{i}") for i in range(mat.shape[1])])
sil = sil[~sil.index.duplicated(keep=False)]

rbh = pd.read_csv(f"{W}/franken_rbh_A9.csv")
rbh = rbh[rbh.seaad_type.isin(sea.columns)]      # only merges whose SEA-AD side survives
comb = sea.join(sil, how="inner", lsuffix="_sea", rsuffix="_sil")
comb = comb[[c for c in comb.columns if c not in set(rbh.siletti_cluster)]]
print(f"combined taxonomy: {comb.shape[1]} types ({len(sea.columns)} DLPFC + "
      f"{comb.shape[1] - len(sea.columns)} Siletti), {comb.shape[0]:,} genes")

comb = comb.loc[comb.sum(axis=1) > 0]
comb = comb * (1000.0 / comb.sum(axis=0))
spec = comb.div(comb.sum(axis=1), axis=0)
gl = pd.read_csv(GLOC, sep=r"\s+", header=None,
                 names=["ENTREZ", "CHR", "START", "STOP", "STRAND", "SYMBOL"], dtype=str)
spec = spec.join(ens_map[ens_map.ENTREZ.isin(set(gl.ENTREZ))].set_index("ENSEMBL"), how="inner")
spec = spec[~spec.ENTREZ.duplicated(keep=False)]
order = list(comb.columns)
safe = {c: f"Type{i}" for i, c in enumerate(order)}
spec[["ENTREZ"] + order].rename(columns={**safe, "ENTREZ": "GENE"}).to_csv(
    f"{W}/spec_dlpfc_a9rbh.txt", sep="\t", index=False)
pd.DataFrame({"safe_name": [safe[c] for c in order], "cell_type": order}).to_csv(
    f"{W}/namemap_a9rbh.csv", index=False)
print(f"wrote spec_dlpfc.txt ({spec.shape[0]:,} genes x {len(order)} types)")
