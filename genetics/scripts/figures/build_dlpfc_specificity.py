"""
DLPFC (SEA-AD A9, 2024-02-13) specificity matrix, built exactly as the MTG one:
per-cell ln(1+raw UMIs) averaged within supertype, pooled over the three
neurotypical reference donors, then merged with Siletti through Duncan's recipe.

The three A9 reference donors (H18.30.002, H19.30.001, H19.30.002) are the same
individuals as three of the five MTG reference donors, so this is close to a
within-donor region swap rather than a change of cohort.
"""
import h5py
import numpy as np
import scipy.sparse as sp
import pandas as pd
import glob

W = "/Users/shreejoy/Github/scz_celltype_paper/genetics/results/intermediates"
SRC = "/Users/shreejoy/Github/scz_cell_type_enrichment"
GLOC = f"{SRC}/linking_cell_types_to_brain_phenotypes/Data/NCBI37.3.gene.loc.extendedMHCexcluded"
BLOCK = 8000

files = sorted(glob.glob("/Users/shreejoy/Downloads/*A9_RNAseq_final-nuclei.2024-02-13.h5ad"))
assert len(files) == 3, files

sums, counts, genes = {}, {}, None
for path in files:
    with h5py.File(path, "r") as f:
        d = f["obs/Supertype"]
        cats = np.array([x.decode() if isinstance(x, bytes) else x
                         for x in f[d.attrs["categories"]][:]])
        codes = d[:]
        if genes is None:
            genes = np.array([g.decode() for g in f["var/_index"][:]])
        n_cells, n_genes = f["X"].attrs["shape"]
        indptr = f["layers/UMIs/indptr"][:]
        donor = np.array([x.decode() if isinstance(x, bytes) else x
                          for x in f[f["obs/Donor ID"].attrs["categories"]][:]])[0]
        for i0 in range(0, n_cells, BLOCK):
            i1 = min(i0 + BLOCK, n_cells)
            p0, p1 = int(indptr[i0]), int(indptr[i1])
            blk = sp.csr_matrix(
                (np.log1p(f["layers/UMIs/data"][p0:p1].astype(np.float64)),
                 f["layers/UMIs/indices"][p0:p1], indptr[i0:i1 + 1] - p0),
                shape=(i1 - i0, n_genes))
            bc = codes[i0:i1]
            for c in np.unique(bc):
                if c < 0:
                    continue
                name = cats[c]
                s = np.asarray(blk[bc == c].sum(axis=0)).ravel()
                sums[name] = sums.get(name, 0) + s
                counts[name] = counts.get(name, 0) + int((bc == c).sum())
        print(f"  {donor}: {n_cells:,} nuclei, {len(np.unique(codes))} supertypes", flush=True)

sea = pd.DataFrame({k: sums[k] / counts[k] for k in sums}, index=genes)
print(f"\npooled DLPFC reference: {sea.shape[0]:,} genes x {sea.shape[1]} supertypes, "
      f"{sum(counts.values()):,} nuclei")

# keep only labels shared with the 137-supertype MTG taxonomy (drops -SEAAD states)
mtg_cols = pd.read_csv(f"{W}/seaad_supertype_log1p_mean.csv", index_col=0, nrows=0).columns
keep = [c for c in sea.columns if c in set(mtg_cols)]
print(f"supertypes shared with the MTG 137: {len(keep)}  "
      f"(dropped {sea.shape[1] - len(keep)}: {sorted(set(sea.columns) - set(keep))})")
sst = [c for c in keep if c.startswith("Sst_")]
print(f"  Sst_ supertypes: {len(sst)}; nuclei "
      f"min {min(counts[c] for c in sst)}, median {int(np.median([counts[c] for c in sst]))}, "
      f"max {max(counts[c] for c in sst)}")
sea = sea[keep]
sea.to_csv(f"{W}/a9_supertype_log1p_mean.csv")

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

rbh = pd.read_csv(f"{W}/franken_rbh_REBUILT.csv")
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
    f"{W}/spec_dlpfc.txt", sep="\t", index=False)
pd.DataFrame({"safe_name": [safe[c] for c in order], "cell_type": order}).to_csv(
    f"{W}/namemap_dlpfc.csv", index=False)
print(f"wrote spec_dlpfc.txt ({spec.shape[0]:,} genes x {len(order)} types)")
