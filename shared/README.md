# shared/ — cross-stream interfaces (NOT the big reference hub)

Small, paper-internal assets that ≥2 subdirs consume.

- [`snrnaseq_de/`](snrnaseq_de/) — the interface to the student's upstream
  pipeline (DE + composition betas). See its README.
- *(planned)* `franken_taxonomy/`, `gwas/` — dedup of the canonical cell-type
  taxonomy and PGC3 GWAS gene set, once `genetics/` is added.

**The large shared REFERENCE data (SEA-AD h5ads, etc.) stays in the EXTERNAL
`~/Github/shared_data/` hub** — it serves non-SCZ projects too. Per
`~/Github/DATA_LAYOUT.md`, consumers symlink into it with absolute paths rather
than copying. This monorepo follows that convention; it does **not** absorb
`shared_data/`.
