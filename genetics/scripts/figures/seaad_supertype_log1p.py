"""
SEA-AD supertype mean expression, computed Duncan-style:
per-cell ln(1+x) on RAW counts, then averaged within supertype.

This is the SEA-AD analogue of their Siletti_create_L2-log_dataset.py. The
existing in-house matrix instead averages CP10K-normalized values, which is a
different quantity; rebuilding it this way makes the SEA-AD and Siletti halves
of the Franken taxonomy commensurable and matches the published recipe.

log1p is applied to the sparse .data array only, which is exact because
log1p(0) = 0, so sparsity is preserved.
"""
import h5py
import numpy as np
import scipy.sparse as sp
import pandas as pd
import time

import os as _os

# Repo-relative, so the chain runs from any clone. External data that is not in
# the repo is still resolved by absolute path or an environment override below.
GEN = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
PAPER = _os.path.dirname(GEN)

H5 = _os.environ.get("SEAAD_MTG_H5AD",
     _os.path.expanduser("~/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad"))
OUT = f"{GEN}/results/intermediates/seaad_supertype_log1p_mean.csv"
BLOCK = 8000

f = h5py.File(H5, "r")
n_cells, n_genes = f["X"].attrs["shape"]
genes = np.array([g.decode() for g in f["var/_index"][:]])
cats = np.array([c.decode() for c in f["obs/Supertype/categories"][:]])
codes = f["obs/Supertype/codes"][:]
indptr = f["X/indptr"][:]
n_types = len(cats)
print(f"{n_cells:,} cells x {n_genes:,} genes; {n_types} supertypes", flush=True)

sums = np.zeros((n_types, n_genes), dtype=np.float64)
counts = np.bincount(codes, minlength=n_types).astype(np.float64)

t0 = time.time()
for i0 in range(0, n_cells, BLOCK):
    i1 = min(i0 + BLOCK, n_cells)
    p0, p1 = int(indptr[i0]), int(indptr[i1])
    block = sp.csr_matrix(
        (np.log1p(f["X/data"][p0:p1].astype(np.float64)),
         f["X/indices"][p0:p1],
         indptr[i0:i1 + 1] - p0),
        shape=(i1 - i0, n_genes))
    # indicator (n_types x block_cells) @ block -> per-type gene sums
    blk_codes = codes[i0:i1]
    ind = sp.csr_matrix(
        (np.ones(i1 - i0), (blk_codes, np.arange(i1 - i0))),
        shape=(n_types, i1 - i0))
    sums += np.asarray((ind @ block).todense())
    if (i0 // BLOCK) % 4 == 0:
        el = time.time() - t0
        print(f"  {i1:,}/{n_cells:,} cells | {el:.0f}s", flush=True)

f.close()
mean = sums / counts[:, None]
df = pd.DataFrame(mean.T, index=genes, columns=cats)
df.index.name = "Gene"
df.to_csv(OUT)
print(f"wrote {OUT}  ({df.shape[0]:,} genes x {df.shape[1]} supertypes) "
      f"in {(time.time()-t0)/60:.1f} min", flush=True)
