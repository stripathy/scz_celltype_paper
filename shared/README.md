# shared/ — cross-component interfaces

Small, paper-internal assets that more than one component consumes. This is
**not** the large reference hub.

- [`snrnaseq_de/`](snrnaseq_de/) — the interface to the upstream snRNA-seq
  meta-analysis (DE + composition betas). See its README.
- `figure_inputs.R` / `figure_inputs.py` — the staleness guard. A figure's
  inputs are committed snapshots, so they can go stale silently; these record
  the upstream file behind each snapshot in a `MANIFEST.tsv` and stop a
  renderer that would otherwise draw old numbers.
- `verify_provenance.py` — re-checks every manifest in the repo at once, and
  exits non-zero if any set has drifted. The single command to run before
  trusting a figure.

```bash
python3 shared/verify_provenance.py            # report, exit 1 if stale
python3 shared/verify_provenance.py --quiet    # only print problems
```

**The large reference data (SEA-AD h5ads, the MERFISH atlas, GWAS summary
statistics) stays external**, in `~/Github/shared_data/` and each component's
git-ignored `data/`. Per `~/Github/DATA_LAYOUT.md`, consumers symlink into the
hub rather than copying, because it serves non-SCZ projects too. This repo
follows that convention and does not absorb it.
