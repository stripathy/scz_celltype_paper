# Retiring `SCZ_Xenium` — parked decisions and plan

**Status: parked 2026-08-06.** No action taken beyond what is already committed.
This repo is being promoted to authoritative; `~/Github/SCZ_Xenium` will eventually
be retired. Nothing here is urgent — it is written so the work can be picked up
cold, months later, without re-deriving any of it.

---

## The decision that was made

`scz_celltype_paper` becomes the authoritative repo for everything, including the
Xenium processing pipeline (cell typing, depth inference, layer inference, spatial
domains). `SCZ_Xenium` is retired once its interactive cell browser is moved or
ported. **This repo will be made public** (timing TBD, likely at submission).

Rationale: the two code trees are already 100% duplicated with zero symlinks — 85
files identical, 5 diverged, 37 paper-only, 12 Xenium-only. The 15-file /
4,545-line cell-typing + depth + layer + spatial-domain core is already present and
git-tracked here. This was never a migration question; it is a question of which of
two existing copies is authoritative, and they have already drifted twice.

---

## Open decisions — these are yours, and nothing should proceed past them

### D1. Does the per-molecule transcript overlay survive publication?

This single decision determines whether 50 GB is a migration problem or a deletion.

- **If yes** — move `output/viewer/transcripts` (50 GB) + `boundaries` (739 MB) +
  `output/deploy` (5.4 GB) out of the repo tree to a plain data volume, and resolve
  `VIEWER_DIR` / `DEPLOY_DIR` from one env var read in `pipeline_config.py`.
- **If no** — publish `index.json` + 24 sample JSONs (98 MB) + boundaries (739 MB)
  and drop transcripts. **This needs zero code changes**: `checkTranscriptData()`
  already removes the Transcripts sidebar panel on a 404, and `loadBoundaryData()`
  uses `Promise.allSettled` and degrades to point rendering.

**The asymmetry that should drive the answer:** boundaries are 24/24 regenerable
from `data/raw` CSVs. Transcripts are only **19/24** — the zarrs for Br6032,
Br6389, Br6437, Br8433 and **Br8667** are absent, and Br8667 is the viewer's
hard-coded default sample. Regenerating those five means a GEO GSE307404
re-download plus a multi-hour zarr scan each (Br8667 alone indexes 281,894,732
transcripts). Only 98 MB is genuinely mandatory at runtime.

### D2. Public / citation strategy

`SCZ_Xenium` is public at `github.com/stripathy/SCZ_Xenium` and is the citable
artifact today. **Deleting the local 160 GB and deleting the GitHub repo are two
separate decisions — keep them decoupled regardless of which option is chosen.**

| Option | Effect | Cost |
|---|---|---|
| Make this repo public at submission + redirect notice on `SCZ_Xenium` | Cleanest going-forward story | Needs a content audit — exposes manuscript drafts and Nicole-owned material |
| Tag `SCZ_Xenium` at final state, set GitHub repo archived/read-only, delete only the local tree | Preserves every existing link and citation | Essentially zero — **lowest-regret if unsure** |
| Mint a Zenodo DOI from a tagged snapshot | Stable citation independent of GitHub; satisfies most data-availability policies | Small one-off |

Inputs not available to whoever wrote this: the target journal's data-availability
policy, and whether `sczxenium.netlify.app` already appears in a preprint or grant.

---

## Verified facts (do not re-derive)

- **Canonical object**: `SCZ_Xenium/output/all_samples_annotated.h5ad`, md5
  `763e6655fc55839f177d57aa98dc5453`, the reinstated 2026-04-01 state. Analysis set
  742,103 cortical → 356,313 neuronal / 385,790 non-neuronal. Corrected Lieber map
  cluster 9 = MGE, cluster 12 = CGE.
- **The viewer port is ~90% done.** All 17 front-end files (156 KB) are already at
  `spatial/output/viewer/`, verified byte-identical to `SCZ_Xenium`. `Makefile`,
  `package.json`, `playwright.config.js`, `tests/` and exporters 03/06/07 are
  already tracked. The only blocker is root `.gitignore:12` `spatial/output/`.
  The port is a `git add -f` on 17 paths — **1–3 hours, not 50 GB.**
- **Netlify**: manual `make deploy`, no CI/CD. siteId
  `78988b9e-d4c1-4dc6-8572-fe326d73bd2c`, recorded in exactly one untracked place
  (`SCZ_Xenium/output/deploy/.netlify/state.json`). Live URL believed to be
  `https://sczxenium.netlify.app` (asserted in `spatial-viewer-core/README.md`;
  never confirmed by a network call).
- **`core/` is vendored** from a third repo, `~/Github/spatial-viewer-core`, pinned
  at `VERSION = 525d6ee` (v0.3.0) while upstream is v0.5.0. Do **not** run
  `make sync-core` during the port — that conflates a move with a version bump.

### The QC divergence — settled, but for a subtler reason than first thought

`modules/cell_qc.py` and `pipeline/02b_run_correlation_classifier.py` differ between
repos; `SCZ_Xenium` is ahead. **Keep this repo's strict `cell_qc.py`** — it is the
pipeline that produced every published number.

An earlier analysis in this session claimed the relaxation was "verified inert"
because all 23,620 high-count cells carry `corr_qc_pass = False`. **That reasoning
was wrong.** `02b_run_correlation_classifier.py:238` assigns `corr_qc_pass = False`
as a *default* for cells that never got classified ("Initialize columns with
defaults for non-QC-pass cells"). Those cells failed strict `qc_pass` upstream, so
02b never evaluated them. Re-run the pipeline under relaxed QC and they get
classified and can enter the analysis set. The conclusion holds for the committed
object; it does **not** license adopting the change.

The two files are a **matched pair**. Adopting `cell_qc.py` without its 02b partner
creates a 25th "Unassigned" centroid absorbing **84,653 cells (6.5% of the analysis
set)** and drops Stage-1 agreement to 93.5%. Paired, centroids are bit-identical
(max diff 0.000e+00 over 24×300). If the newer version is ever adopted, adopt both,
and add an assert in 02b that raises if any centroid-building cell has
`subclass_label == 'Unassigned'` — that converts an unpaired adoption from silent
to loud.

---

## What actually gates deletion

Four things, none of which are the viewer:

1. **89 uncommitted paths in `SCZ_Xenium`**, including the entire 2026-07-31 April
   reinstatement. `HEAD == origin/main` at 2026-05-02, so **the public repo does not
   contain the current work**. `rm -rf` destroys it with no remote copy anywhere.
2. **~15 files exist in exactly one place** — untracked in `SCZ_Xenium` git *and*
   absent here. Highest loss-risk-per-hour item in the whole plan:
   - `code/pipeline/08_export_clean_h5ad.py` — sole producer of the
     collaborator-facing clean h5ad
   - `code/analysis/run_crumblr_prs.R` — sole producer of `output/crumblr/*_prs.csv`
     and `output/prs_analysis/`, the compositional-vs-polygenic-risk analysis
   - `code/analysis/qc_bias_{diagnose_high_counts_filter,failure_reason_breakdown,recover_high_marker_cells}.py`
     — the reviewer-anticipation control on whether QC preferentially discards
     high-SST/PVALB cells in SCZ
   - `code/analysis/plot_supertype_depth_violins_by_dx.py` — the
     diagnosis-dependent-misclassification check (this repo has only the non-by-dx variant)
   - `prototype_depth_smoothness{,_compare,_domain,_radius,_zoom}.py`,
     `prototype_sst_neighborhood_classifier.py` — methods provenance, belongs in
     `spatial/code/archive/`
3. **Single-copy docs**, including one citation that is *already* broken:
   `docs/pipeline_qc_audit.md` (20 KB) is cited by `spatial/methods_writeup.md`
   (lines 70, 95, 123, 330), `spatial/README.md:363` and
   `spatial/code/nuclear_resolution/README.md:11` — and is gitignored in **both**
   repos, so it is dead for every collaborator who clones today. Also
   `docs/{sst_vulnerability_interpretive_framework, L6b_findings_synthesis,
   dienel_robustness_analysis_plan, qc_bias_high_marker_cells_plan,
   SEAAD_MERFISH_dataset}.md` (68 KB of original interpretive writing, one co-owned
   with Nicole), `docs/papers/literature_review_calbindin_dbc.md` (46 KB, underlies
   the CALB1 narrative in Fig 4), and `archive/README.md` (the only written record
   of the June-11 depth-model regression and the unresolved L1 discrepancy,
   April 72,078 vs June 30,577).
4. **The live site is stale.** `output/deploy/Br*.json` are dated 2026-03-15 while
   `output/viewer/Br*.json` are 2026-05-01 and differ byte-for-byte — the published
   viewer serves **March cell data under May code**, i.e. it does not reflect the
   canonical April annotation. The deploy's `transcripts/` is a trimmed 6-sample
   subset with an incomplete MT-gene purge (9 `MT-*`/`MTRNR*` files vs 216 in
   `viewer/`). Redeploy from current data **before** deleting anything, or the
   published artifact can never be reconstructed.

---

## Plan

### Stage 3 — promote to authoritative (fully reversible, `SCZ_Xenium` untouched)

| | Step | Effort |
|---|---|---|
| 3.1 | Lock the QC pair; port only 02b's Unassigned-guard as a defensive no-op; add the loud assert; write the diff + empirical result into `spatial/code/archive/` as a note | 1–2 h |
| 3.2 | Reconcile the three conflicting `corr_qc_pass` totals — `methods_writeup.md:49` says 1,225,037, committed `pipeline_cell_flow.csv` sums to 1,297,413 with `hybrid_qc_pass` **exceeding** `qc_pass` (arithmetically impossible under strict QC), verified value from the canonical h5ads is **1,221,519**. Regenerate the CSV, fix the doc | 1–2 h |
| 3.3 | **Track the viewer front-end** — `git add -f` the 17 paths; move `spatial/.github/workflows/test.yml` to repo-root `.github/workflows/` (GitHub only reads root) with paths rewritten and `working-directory: spatial`; record the Netlify siteId in a tracked file | 1–3 h |
| 3.4 | **Rescue the single-copy code** (item 2 above) | 1 h |
| 3.5 | **Rescue the single-copy docs** (item 3 above); un-ignore `docs/pipeline_qc_audit.md` by narrowing `spatial/.gitignore`'s `docs/` rule to `docs/papers/` so copyrighted PDFs stay out | 1–2 h |
| 3.6 | Introduce the seam but **do not flip the default**: change to `BASE_DIR = os.environ.get('XENIUM_BASE', <script-relative spatial/>)` in `pipeline_config.py:17` and `config.py:37` (+ ~30 repeats), and set `XENIUM_BASE=~/Github/SCZ_Xenium` in your shell so nothing breaks yet | 2–3 h |

72 literal `Github/SCZ_Xenium` references exist (69 code, 3 docs); 50 are the bare
root constant, and because `spatial/` already mirrors the layout, one flip resolves
all 50 with no further edits. Full list was written to a scratch file during the
audit — regenerate with
`grep -rn 'Github/SCZ_Xenium' --include='*.py' --include='*.R' --include='*.sh' --include='*.md' .`

### Stage 4 — make this repo stand alone

| | Step | Effort |
|---|---|---|
| 4.1 | Move the ~4.6 GB must-keep payload into `spatial/` (same volume, instant): the canonical object, 24 per-sample h5ads + `correlation_centroids.pkl`, `PRISTINE.h5ad` (the only unaugmented provenance copy), the clean h5ad + guide, `depth_model*.pkl` (not bit-reproducible without the 3.1 GB SEA-AD reference), `output/de/` incl. `MANIFEST.tsv`, crumblr manifests + `*_prs.csv` + `_pre_qc_relax/`, and CSV-only from the `qc_bias`/`prs_analysis`/etc. dirs. Send `archive/2026-06-11/` (2.6 GB, superseded) to external storage, not into the repo | 2–3 h |
| 4.2 | Cut the three live symlinks (they dangle silently, no error until something reads a missing file): `shared/snrnaseq_de/nicole_scz_snrnaseq_betas` → `SCZ_Xenium/data/…` (22 MB, includes `sst_meta_de/` 19 CSVs that exist nowhere else and are untracked even there — input behind the Fig 4 SCZ–AD panel). ~~`spatial/output/de/de_results_{subclass,supertype}.csv`~~ — **done 2026-08-06, commit `fe72b5d`** | 1–2 h |
| 4.3 | Flip the `BASE_DIR` default, then **rehearse**: rename `SCZ_Xenium` → `SCZ_Xenium.hidden` and confirm a representative script from each seam runs green | 2 h |

### Stage 5 — the viewer gate

| | Step |
|---|---|
| 5.1 | Inventory what is actually live **while `SCZ_Xenium` is intact** — confirm the URL, the siteId, and exactly which samples and data vintage are published. After deletion this cannot be reconstructed |
| 5.2 | **Decision D1** |
| 5.3 | Redeploy from current data (better: re-run `06_export_viewer.py` from the canonical April h5ads so the published viewer matches the paper). If D1 = no, the deploy is ~840 MB and fully reproducible from h5ads + boundary CSVs, both 24/24 present |
| 5.4 | Execute D1 — relocate or retire the payload. `output/merscope_viewer/` (0.58 GB, older monolithic viewer for the 14 SEA-AD MERSCOPE sections, never deployed) can be regenerated from `output/merscope_h5ad/` (223 MB, present): migrate the exporters, drop the export |
| 5.5 | Fix or consciously accept three viewer defects — see below |

**Three viewer defects** (none block the port, all mislead if left):
1. `code/modules/bundle_viewer.py` can no longer build the offline standalone HTML.
   It string-matches an inline `<script>` block and a fetch-based `init()`/`loadSample()`
   deleted in the 2026-05-02 refactor; `html.index('<script>')` now raises
   `ValueError`, which `06_export_viewer.py` silently swallows as "HTML viewer build
   failed". The 25 MB file on disk is a stale pre-refactor artefact — **yet the
   README tells people to open it.** At minimum make the failure loud.
2. The Playwright workflow is almost certainly red — it serves `output/viewer/` over
   `http.server`, but `.gitignore` excludes the `*.json`, so a fresh CI checkout has
   code and zero data and `01_load.spec.js` times out waiting on `window.indexData`.
   Add a fixture step or mark the workflow manual.
3. Note in the README that `core/` is pinned at v0.3.0 while `spatial-viewer-core`
   is v0.5.0, so the divergence reads as deliberate rather than neglect.

### Stage 6 — deletion

| | Step |
|---|---|
| 6.1 | **Decision D2** |
| 6.2 | Commit or export the 89 dirty paths in `SCZ_Xenium` — **before any deletion.** Explicitly capture the SCZ_Xenium-ahead versions of `cell_qc.py` and `02b` as a diff note in `spatial/code/archive/`; do not overwrite this repo's files with them, but do not let the change vanish |
| 6.3 | Final checklist, then delete **the local tree only** |

**Pre-deletion checklist — all must be green:**
1. `SCZ_Xenium` renamed `.hidden` for 48 h with no script, figure rebuild, or notebook failing
2. Zero live `Github/SCZ_Xenium` grep hits; zero dangling symlinks
   (`find . -type l -exec test ! -e {} \; -print`)
3. md5 of the migrated canonical object still `763e6655fc55839f177d57aa98dc5453`
4. crumblr and DE outputs regenerate byte-identical from the migrated h5ads
   (the `MANIFEST.tsv` files make this mechanical)
5. The Netlify site serves the vintage decided in 5.3
6. D2 executed — tag pushed and/or DOI minted
7. 6.2 clean

Then delete the local directory. **Leave the GitHub repo alone unless D2 says
otherwise.** Reclaims ~160 GB.

---

## Provenance of this document

Written 2026-08-06 from a four-agent audit (viewer stack, migration gap, QC
divergence + path seams, synthesis). Claims marked "verified" were checked directly
against files on disk. The live Netlify URL was **not** confirmed by a network call.
Related: `DATA_FLOW.md` (object state and seams), `spatial/output/crumblr/README.md`
and `spatial/output/de/README.md` (canonical result sets).
