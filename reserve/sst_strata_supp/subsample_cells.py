#!/usr/bin/env python3
"""Cell-matched pseudobulks: the definitive test of the cell-count confound.

THE PROBLEM. The depletion groups are defined by cases having fewer cells, so
within the depleted group SCZ donors contribute ~30% fewer nuclei and ~36% fewer
counts to their pseudobulk than controls. Nuclei count predicts module scores
even among control donors alone, so some of the graded signal could be technical.
Read-depth matching (binomial thinning of control pseudobulks) left the result
intact, but thinning reproduces lost sequencing depth, not lost cells.

WHAT THIS DOES. Subsamples control donors' CELLS down to the case distribution,
within each cohort and stratum, then rebuilds the pseudobulks from scratch. That
removes the cell-count difference by construction rather than by modelling, and
it is the one correction the pseudobulks alone could not support.

Matching: donors are first restricted to those with >= MIN_CELLS nuclei in the
stratum (as the main pipeline does), then each control donor is assigned a target
cell count drawn from the case distribution at the control's own percentile rank,
so the whole distribution matches rather than just the mean. Targets therefore
never fall below MIN_CELLS and no donor is dropped by the matching itself. Cases
keep all their cells.

Output mirrors the layout of pseudobulk/donor_stratum/ so the existing R DE code
runs against it unchanged, one directory per random draw.

Usage:  python3 reserve/sst_strata_supp/subsample_cells.py [n_reps]
"""
import sys, os, warnings
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
warnings.filterwarnings("ignore")

N_REPS = int(sys.argv[1]) if len(sys.argv) > 1 else 5
MIN_CELLS = 10
CELLS = "transcriptomic/data/stratum_sst_cells_export"
OUT = "transcriptomic/results/sst_strata_gsea/supp/cellsub"
STRATA_CSV = "transcriptomic/results/sst_strata_gsea/strata_definition.csv"

strata = pd.read_csv(STRATA_CSV)
st_of = dict(zip(strata.CellType, strata.stratum))
cohorts = sorted(f.replace("_sst_cells.h5ad", "")
                 for f in os.listdir(CELLS) if f.endswith("_sst_cells.h5ad"))
print(f"cohorts: {cohorts}\nreplicates: {N_REPS}\n", flush=True)


def symbol_vector(var):
    """Resolve to gene symbols exactly as 02_stratum_de.R / 04_donor_pseudobulks.R do."""
    native = var.index.to_numpy().astype(object)
    sym = var["gene_symbol"].to_numpy().astype(object) if "gene_symbol" in var else native.copy()
    looks_sym = ~pd.Series(native).astype(str).str.startswith("ENSG").to_numpy()
    miss = pd.isna(sym)
    sym[miss & looks_sym] = native[miss & looks_sym]
    return sym


def pseudobulk(X, row_group, groups):
    """Sum rows of sparse X within each group. groups is the ordered group list."""
    idx = {g: i for i, g in enumerate(groups)}
    rows = np.fromiter((idx[g] for g in row_group), dtype=np.int64, count=len(row_group))
    M = sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, np.arange(len(rows)))),
                      shape=(len(groups), X.shape[0]))
    return np.asarray((M @ X).todense())


for co in cohorts:
    a = ad.read_h5ad(f"{CELLS}/{co}_sst_cells.h5ad")
    obs = a.obs.copy()
    obs["donor"] = obs["donor"].astype(str)
    obs["stratum"] = obs["supertype"].map(st_of)
    obs = obs.reset_index(drop=True)
    keep_cell = obs["stratum"].notna().to_numpy()

    sym = symbol_vector(a.var)
    ok_gene = pd.notna(sym)
    sym_ok = pd.Index(sym[ok_gene].astype(str))
    # collapse duplicate symbols by summing their columns, once per cohort
    uniq = sym_ok.unique()
    cidx = pd.Series(np.arange(len(uniq)), index=uniq)
    colmap = sp.csr_matrix(
        (np.ones(len(sym_ok), dtype=np.float32),
         (np.arange(len(sym_ok)), cidx.loc[sym_ok].to_numpy())),
        shape=(len(sym_ok), len(uniq)))
    X = a.X[:, ok_gene].astype(np.float32).tocsr()
    print(f"{co}: {X.shape[0]:,} nuclei, {len(uniq):,} symbols "
          f"({int((~ok_gene).sum()):,} unmapped dropped)", flush=True)

    # per donor x stratum cell inventory, restricted as the main pipeline does
    inv = (obs[keep_cell].groupby(["donor", "stratum", "diagnosis"], observed=True)
           .size().rename("n").reset_index())
    inv = inv[inv.n >= MIN_CELLS]

    for rep in range(1, N_REPS + 1):
        rng = np.random.default_rng(1000 + rep)
        sel_rows, grp_labels = [], []
        for st, blk in inv.groupby("stratum", observed=True):
            scz = blk[blk.diagnosis == "SCZ"].n.to_numpy()
            ctl = blk[blk.diagnosis == "Control"]
            if len(scz) < 3 or len(ctl) < 3:
                targets = dict(zip(blk.donor, blk.n))          # too few to match
            else:
                pr = (pd.Series(ctl.n.to_numpy()).rank(method="average").to_numpy()
                      / (len(ctl) + 1))
                tg = np.quantile(scz, pr)
                targets = dict(zip(blk[blk.diagnosis == "SCZ"].donor,
                                   blk[blk.diagnosis == "SCZ"].n))
                targets.update({d: int(min(n, np.floor(t)))
                                for d, n, t in zip(ctl.donor, ctl.n, tg)})
            for d, n in zip(blk.donor, blk.n):
                pool = np.where(keep_cell & (obs.donor == d).to_numpy()
                                & (obs.stratum == st).to_numpy())[0]
                k = int(targets[d])
                take = pool if k >= len(pool) else rng.choice(pool, size=k, replace=False)
                sel_rows.append(take)
                grp_labels.extend([f"{d}|{st}"] * len(take))
        sel_rows = np.concatenate(sel_rows)
        groups = sorted(set(grp_labels))
        pb = pseudobulk(X[sel_rows], grp_labels, groups) @ colmap.toarray()

        d = f"{OUT}/rep{rep}"
        os.makedirs(d, exist_ok=True)
        pd.DataFrame(pb.T, index=pd.Index(uniq, name="gene"), columns=groups) \
            .reset_index().to_parquet(f"{d}/{co}_donor_stratum_counts.parquet", index=False)
        meta = pd.DataFrame({"key": groups})
        meta[["donor", "stratum"]] = meta.key.str.split("|", expand=True)
        cnt = pd.Series(grp_labels).value_counts()
        meta["n_cells"] = meta.key.map(cnt).astype(int)
        meta = meta.merge(inv[["donor", "stratum", "diagnosis"]], on=["donor", "stratum"],
                          how="left")
        meta["cohort"] = co
        meta[["donor", "stratum", "n_cells", "diagnosis", "cohort"]] \
            .to_csv(f"{d}/{co}_donor_stratum_meta.csv", index=False)
        if rep == 1:
            chk = meta.groupby("diagnosis").n_cells.median()
            print(f"   rep1 median nuclei/donor after matching: "
                  f"{dict(chk.round(1))}", flush=True)
    del a, X
print("\ndone", flush=True)
