"""
Faithful re-implementation of Duncan et al.
Preprocessing_Siletti/create_matrices/Siletti_create_L2-log_dataset.py

Same arithmetic as the original — per cluster, mean over cells of ln(1+x):
    avg[:, k] = sum(log1p(matrix[:, cells_in_k]), axis=1) / n_cells_in_k

Two deviations, both performance-only and numerically neutral:
  * cells are perfectly grouped by cluster in this loom (461 contiguous runs),
    so each cluster is read as a contiguous slice instead of by fancy indexing;
  * each cluster is read in column blocks and accumulated, so peak memory stays
    bounded regardless of cluster size (the largest cluster would otherwise
    need ~100 GB), and clusters are processed by a worker pool.
"""
import h5py
import numpy as np
from multiprocessing import Pool
import time
import sys

LOOM = "/Users/shreejoy/Github/scz_cell_type_enrichment/data/adult_human_20221007.loom"
OUT = "/Users/shreejoy/Github/scz_celltype_paper/genetics/results/intermediates/Siletti_L2-cluster-log1p_matrix.h5"
BLOCK = 2000          # cells per read
N_WORKERS = 6         # disk-bound; more workers stop helping

_f = None


def _init():
    global _f
    _f = h5py.File(LOOM, "r")


def _cluster_mean(task):
    """Mean of log1p over one cluster's contiguous column range."""
    k, start, stop = task
    dset = _f["matrix"]
    n_genes = dset.shape[0]
    acc = np.zeros(n_genes, dtype=np.float64)
    for b0 in range(start, stop, BLOCK):
        b1 = min(b0 + BLOCK, stop)
        acc += np.log1p(dset[:, b0:b1]).sum(axis=1)
    return k, acc / (stop - start), stop - start


def main():
    with h5py.File(LOOM, "r") as f:
        clusters = f["col_attrs/Clusters"][:]
        accession = f["row_attrs/Accession"][:]
        gene = f["row_attrs/Gene"][:]
        n_genes = f["matrix"].shape[0]

    # contiguous [start, stop) per cluster
    edges = np.flatnonzero(np.r_[True, clusters[1:] != clusters[:-1], True])
    tasks = [(int(clusters[edges[i]]), int(edges[i]), int(edges[i + 1]))
             for i in range(len(edges) - 1)]
    n_clusters = int(clusters.max()) + 1
    assert len(tasks) == n_clusters, f"{len(tasks)} runs != {n_clusters} clusters"
    print(f"{n_genes:,} genes x {n_clusters} clusters; "
          f"cells {min(t[2]-t[1] for t in tasks):,}-{max(t[2]-t[1] for t in tasks):,} per cluster",
          flush=True)

    avg = np.zeros((n_genes, n_clusters), dtype=np.float64)
    t0 = time.time()
    # biggest clusters first so the pool drains evenly
    tasks.sort(key=lambda t: t[1] - t[2])
    with Pool(N_WORKERS, initializer=_init) as pool:
        for done, (k, col, n) in enumerate(pool.imap_unordered(_cluster_mean, tasks), 1):
            avg[:, k] = col
            if done % 20 == 0 or done == n_clusters:
                el = time.time() - t0
                print(f"  {done}/{n_clusters} clusters | {el/60:.1f} min | "
                      f"eta {el/done*(n_clusters-done)/60:.1f} min", flush=True)

    with h5py.File(OUT, "w") as w:
        w.create_dataset("matrix", data=avg)
        w.create_dataset("Cluster", data=np.array(
            [f"Cluster{i}".encode("UTF-8") for i in range(n_clusters)]))
        w.create_dataset("Accession", data=accession)
        w.create_dataset("Gene", data=gene)
    print(f"wrote {OUT} in {(time.time()-t0)/60:.1f} min", flush=True)


if __name__ == "__main__":
    sys.exit(main())
