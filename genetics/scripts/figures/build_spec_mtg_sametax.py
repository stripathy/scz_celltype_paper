"""
Specificity matrix with SEA-AD expression taken from MTG instead of A9/DLPFC,
holding the taxonomy fixed.

This is a robustness check on the *reference region*, not on the taxonomy. The
501-type definition, the RBH removals and every downstream step are identical to
the published build; the only change is that the SEA-AD half of the matrix is
populated from the MTG reference rather than A9. MTG is restricted to the same
125 supertypes A9 contributes, so the type list is unchanged and the two runs
are directly comparable cell type for cell type.
"""
import numpy as np
import pandas as pd
import h5py

GEN = "/Users/shreejoy/Github/scz_celltype_paper/genetics"
I = f"{GEN}/results/intermediates"
W = I
SRC = "/Users/shreejoy/Github/scz_cell_type_enrichment"
GLOC = f"{SRC}/linking_cell_types_to_brain_phenotypes/Data/NCBI37.3.gene.loc.extendedMHCexcluded"

# the published taxonomy: same 501 columns, in the same order
namemap = pd.read_csv(f"{I}/namemap_a9rbh.csv")
seaad_types = [c for c in namemap.cell_type if not c.startswith(("Mgl_", "Splat_"))]
a9_cols = set(pd.read_csv(f"{I}/a9_supertype_log1p_mean.csv", index_col=0, nrows=0).columns) \
    if __import__("os").path.exists(f"{I}/a9_supertype_log1p_mean.csv") else None

sea = pd.read_csv(f"{I}/seaad_supertype_log1p_mean.csv", index_col=0)
keep = [c for c in sea.columns if c in set(namemap.cell_type)]
sea = sea[keep]
print(f"MTG reference restricted to the published type list: {sea.shape[0]:,} genes x {sea.shape[1]} supertypes")

ens_map = pd.read_csv(f"{I}/map_entrez_ensembl.csv", dtype=str)
sym_map = pd.read_csv(f"{I}/map_entrez_symbol.csv", dtype=str)
s2e = sym_map.merge(ens_map, on="ENTREZ")[["SYMBOL", "ENSEMBL"]]
s2e = s2e[~s2e.SYMBOL.duplicated(keep=False) & ~s2e.ENSEMBL.duplicated(keep=False)]
sea = sea[~sea.index.duplicated(keep=False)].join(
    s2e.set_index("SYMBOL"), how="inner").set_index("ENSEMBL")

with h5py.File(f"{I}/Siletti_L2-cluster-log1p_matrix.h5", "r") as f:
    mat = f["matrix"][:]
    acc = np.array([a.decode().split(".")[0] for a in f["Accession"][:]])
names = pd.read_csv(f"{I}/siletti_cluster_names.csv")
n2n = dict(zip(names.cluster_num, names.cluster_name))
sil = pd.DataFrame(mat, index=acc,
                   columns=[n2n.get(i, f"Cluster{i}") for i in range(mat.shape[1])])
sil = sil[~sil.index.duplicated(keep=False)]

# identical RBH removals to the published build
rbh = pd.read_csv(f"{I}/franken_rbh_A9.csv")
rbh = rbh[rbh.seaad_type.isin(sea.columns)]
comb = sea.join(sil, how="inner", lsuffix="_sea", rsuffix="_sil")
comb = comb[[c for c in comb.columns if c not in set(rbh.siletti_cluster)]]
print(f"combined taxonomy: {comb.shape[1]} types, {comb.shape[0]:,} genes")

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
    f"{W}/spec_mtg_sametax.txt", sep="\t", index=False)
pd.DataFrame({"safe_name": [safe[c] for c in order], "cell_type": order}).to_csv(
    f"{W}/namemap_mtg_sametax.csv", index=False)
print(f"wrote spec_mtg_sametax.txt ({spec.shape[0]:,} genes x {len(order)} types)")
