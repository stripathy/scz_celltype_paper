"""
SEA-AD-only specificity matrices and their MAGMA gene-property runs.

Replaces the combined SEA-AD + Siletti (RBH) taxonomy of build_spec_a9rbh.py:
following reviewer feedback from L. Duncan, enrichment is presented on the
SEA-AD taxonomy itself, so specificity is computed across each region's own
supertypes with no whole-brain merge. The recipe is otherwise identical to
build_spec_a9rbh.py -- same log1p supertype means, same symbol -> ENSEMBL ->
ENTREZ mapping, same scale-to-1000 -> row-normalise -> gene.loc (MHC-excluded)
join order -- minus the Siletti inner join and RBH column removal. Dropping
the merge also drops its side effect on the gene universe (the Siletti join
kept only shared genes); a gene-universe-matched check showed that difference
is negligible (beta r = 0.999).

Writes, into results/intermediates/:
  spec_dlpfc_a9only.txt + namemap_a9only.csv    125 DLPFC supertypes (main)
  spec_mtg_only.txt     + namemap_mtgonly.csv   137 MTG supertypes (robustness)
  T_a9only_{bigdeli,pgc3}.gsa.out               MAGMA gene-property, one-sided
  T_mtgonly_{bigdeli,pgc3}.gsa.out              (direction=greater), per GWAS
"""
import subprocess

import pandas as pd

GEN = "/Users/shreejoy/Github/scz_celltype_paper/genetics"
W = f"{GEN}/results/intermediates"
SRC = "/Users/shreejoy/Github/scz_cell_type_enrichment"
GLOC = f"{SRC}/linking_cell_types_to_brain_phenotypes/Data/NCBI37.3.gene.loc.extendedMHCexcluded"
MAGMA = f"{GEN}/data/magma/magma_mac/magma"
GENES_RAW = {
    "bigdeli": f"{GEN}/data/gwas/magma_bigdeli/bigdeli.step2.genes.raw",
    "pgc3": f"{GEN}/data/gwas/magma_pgc3/pgc3.step2.genes.raw",
}
REGIONS = {  # tag -> (per-supertype ln(1+UMI) means, spec basename, namemap basename)
    "a9only": (f"{W}/a9_supertype_log1p_mean.csv", "spec_dlpfc_a9only", "namemap_a9only"),
    "mtgonly": (f"{W}/seaad_supertype_log1p_mean.csv", "spec_mtg_only", "namemap_mtgonly"),
}

ens_map = pd.read_csv(f"{W}/map_entrez_ensembl.csv", dtype=str)
sym_map = pd.read_csv(f"{W}/map_entrez_symbol.csv", dtype=str)
s2e = sym_map.merge(ens_map, on="ENTREZ")[["SYMBOL", "ENSEMBL"]]
s2e = s2e[~s2e.SYMBOL.duplicated(keep=False) & ~s2e.ENSEMBL.duplicated(keep=False)]
gl = pd.read_csv(GLOC, sep=r"\s+", header=None,
                 names=["ENTREZ", "CHR", "START", "STOP", "STRAND", "SYMBOL"], dtype=str)
loc_entrez = set(gl.ENTREZ)

for tag, (means_csv, spec_name, nm_name) in REGIONS.items():
    sea = pd.read_csv(means_csv, index_col=0)
    print(f"[{tag}] reference: {sea.shape[0]:,} genes x {sea.shape[1]} supertypes")
    sea = sea[~sea.index.duplicated(keep=False)].join(
        s2e.set_index("SYMBOL"), how="inner").set_index("ENSEMBL")

    mat = sea.loc[sea.sum(axis=1) > 0]
    mat = mat * (1000.0 / mat.sum(axis=0))
    spec = mat.div(mat.sum(axis=1), axis=0)
    spec = spec.join(ens_map[ens_map.ENTREZ.isin(loc_entrez)].set_index("ENSEMBL"),
                     how="inner")
    spec = spec[~spec.ENTREZ.duplicated(keep=False)]
    order = list(mat.columns)
    safe = {c: f"Type{i}" for i, c in enumerate(order)}
    spec[["ENTREZ"] + order].rename(columns={**safe, "ENTREZ": "GENE"}).to_csv(
        f"{W}/{spec_name}.txt", sep="\t", index=False)
    pd.DataFrame({"safe_name": [safe[c] for c in order], "cell_type": order}).to_csv(
        f"{W}/{nm_name}.csv", index=False)
    print(f"[{tag}] wrote {spec_name}.txt ({spec.shape[0]:,} genes x {len(order)} types)")

    for gwas, raw in GENES_RAW.items():
        out = f"{W}/T_{tag}_{gwas}"
        subprocess.run([MAGMA, "--gene-results", raw,
                        "--gene-covar", f"{W}/{spec_name}.txt",
                        "--model", "direction=greater", "--out", out],
                       check=True, capture_output=True, text=True)
        n = sum(1 for line in open(f"{out}.gsa.out") if not line.startswith("#")) - 1
        print(f"[{tag}] MAGMA x {gwas}: {out}.gsa.out ({n} types)")
