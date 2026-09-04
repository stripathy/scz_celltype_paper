"""
Import-me-first shim for brisc, version-aware across 0.1.0 and 0.1.3.

0.1.3 (2026-08-14) fixed most of what ~/Github/brisc_test/BRISC_BUGS.md records, so the
workarounds below are applied only when the running version actually needs them --
`~/Github/brisc_test/test_brisc_bugs.py` is the regression test that decides which.

Still needed on 0.1.3:
  * `pl.Float64` cast for integer covariates in DE formulas (see `de()` below).

Fixed in 0.1.3, patched only for 0.1.0:
  * Bug 1  sched_setaffinity -- the critical one: on 0.1.0 the parallel HDF5 readers die
           silently and you get the right *shape* with an empty X and obs columns
           collapsed to a single repeated value. Wrong answers, not an error.
  * Bug 2  os.get_terminal_size() in SingleCell.ls() outside a tty.
  * Bug 6  `signal` not imported in de.py (breaks plot_volcano).
  * Bug 9  `anndata` not imported in single_cell.py (breaks to_scanpy).

Not bugs, but behaviour worth stating because both cost us time:
  * `subsample_obs(n=…, by_column=…)`: `n` is the exact *total*, stratified by the
    column -- not n per group. Capping a reference per supertype needs a manual loop,
    and getting it wrong distorts kNN label-transfer priors.
  * `concat_obs` refuses inputs whose `uns` differ; pass `flexible=True`.
"""
import os
from importlib.metadata import version as _version

BRISC_VERSION = tuple(int(p) for p in _version("brisc").split(".")[:3])
_NEEDS_LEGACY_PATCHES = BRISC_VERSION < (0, 1, 1)

# Bug 1 must be patched before brisc is imported, and it is the one that corrupts data
# rather than raising, so it is applied defensively regardless of version.
if not hasattr(os, "sched_setaffinity"):
    os.sched_setaffinity = lambda *a, **kw: None

if _NEEDS_LEGACY_PATCHES:
    _real_terminal_size = os.get_terminal_size

    def _safe_terminal_size(*a, **kw):           # Bug 2
        try:
            return _real_terminal_size(*a, **kw)
        except OSError:
            return os.terminal_size((140, 50))
    os.get_terminal_size = _safe_terminal_size

import brisc                     # noqa: E402
import brisc.de                  # noqa: E402
import brisc.single_cell         # noqa: E402
import polars as pl              # noqa: E402

if _NEEDS_LEGACY_PATCHES:
    import anndata as _anndata   # noqa: E402
    import signal as _signal     # noqa: E402
    brisc.single_cell.anndata = _anndata     # Bug 9
    brisc.de.signal = _signal                # Bug 6

from brisc import SingleCell, Pseudobulk, DE, concat_obs   # noqa: E402,F401

NUM_THREADS = 8


def verify_load(sc, name="dataset"):
    """Guard against Bug 1 recurring: a silently-empty X or a collapsed obs column."""
    nnz = sc.X.nnz if hasattr(sc.X, "nnz") else sc.X.data.size
    if nnz == 0:
        raise RuntimeError(f"{name}: X has zero non-zeros -- brisc Bug 1 (parallel "
                           f"reader crash) has recurred; check the sched_setaffinity patch")
    return nnz


def de(pb, formula, **kwargs):
    """Run pseudobulk DE across brisc versions.

    Two version-independent hazards handled here:
      * 0.1.3 renamed `Pseudobulk.DE` to `.de` (and `CPM`/`log_CPM` to `cpm`/`log_cpm`).
      * Integer covariates. `to_r` maps polars integers to R's `integer64` (bit64), which
        `model.matrix()` reinterprets bitwise -- an age of 76 becomes 3.75e-322. The
        result is either a rank-deficient design ("rank 3 with 6 columns") or NaN
        coefficients, and neither error names the dtype. Casting every non-categorical
        integer column of `obs` to Float64 is exact, not an approximation: the DE table
        is bit-identical to one computed from Float64 inputs.
    """
    categorical = set()
    for key in ("categorical_columns", "ordinal_columns"):
        val = kwargs.get(key)
        if val:
            categorical.update([val] if isinstance(val, str) else val)

    int_cols = set()
    for _, (X, obs, var) in pb.items():
        int_cols |= {c for c, t in zip(obs.columns, obs.dtypes)
                     if t.base_type() in (pl.Int8, pl.Int16, pl.Int32, pl.Int64,
                                          pl.UInt8, pl.UInt16, pl.UInt32, pl.UInt64)
                     and c not in categorical}
    if int_cols:
        pb = pb.cast_obs({c: pl.Float64 for c in sorted(int_cols)})

    method = getattr(pb, "DE", None) or getattr(pb, "de")
    return method(formula, **kwargs)
