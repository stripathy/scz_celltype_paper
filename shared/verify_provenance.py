#!/usr/bin/env python3
"""Check every provenance manifest in the repo at once.

Each component records, beside its outputs, a MANIFEST.tsv naming the upstream
files those outputs were built from and their fingerprints at build time. This
script re-checks all of them and reports anything that has drifted — the single
command to run before trusting a figure or committing a result.

    python3 shared/verify_provenance.py            # report, exit 1 if stale
    python3 shared/verify_provenance.py --quiet    # only print problems

Exit status is 0 when everything is consistent, 1 when any set is stale, so it
works as a pre-commit or CI check.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
sys.path.insert(0, str(HERE))
from figure_inputs import problems, read_manifest  # noqa: E402

# The one place this repo reaches into the Xenium processing repo. Override with
# XENIUM_BASE when it is not checked out beside this one.
XEN = Path(os.environ.get("XENIUM_BASE", Path.home() / "Github/SCZ_Xenium")) / "output"

# snapshot/output dir -> (label, command that rebuilds it)
SETS = [
    (PAPER / "transcriptomic/data/figure_inputs",
     "Figure 2 inputs",
     "cd transcriptomic && Rscript scripts/00_refresh_figure_inputs.R"),
    (PAPER / "genetics/results/figures/r_panels",
     "Figure 4 panel CSVs",
     "re-run the genetics/scripts/figures/export_panel*.py exporters "
     "for the stale panels"),
    (XEN / "crumblr",
     "Xenium crumblr inputs",
     "cd spatial && python3 code/analysis/build_crumblr_input.py"),
    (XEN / "de",
     "Xenium DE pseudobulk",
     "cd spatial && python3 code/analysis/build_de_input.py"),
]


def main() -> int:
    quiet = "--quiet" in sys.argv
    stale = 0
    skipped = []
    for d, label, cmd in SETS:
        if not (d / "MANIFEST.tsv").exists():
            skipped.append(label)
            if not quiet:
                print(f"—  {label:<24} NOT CHECKED — no manifest at {d}")
            continue
        probs = problems(d)
        n = len({r["file"] for r in read_manifest(d)})
        if probs:
            stale += 1
            print(f"STALE  {label} — {len(probs)} problem(s)")
            for p in probs:
                print(p)
            print(f"       Fix: {cmd}\n")
        elif not quiet:
            print(f"ok     {label:<24} {n} files, all sources match")
    if stale:
        print(f"\n{stale} of {len(SETS)} sets are stale.")
        return 1
    if skipped:
        # Two of the four sets live in the external Xenium repo, so on a clone
        # they have no manifest to check. Saying "all consistent" there would be
        # a green light for checks that never ran.
        print(f"\n{len(SETS) - len(skipped)} of {len(SETS)} sets checked and consistent; "
              f"{len(skipped)} NOT CHECKED ({', '.join(skipped)}).")
        print("Set XENIUM_BASE to the Xenium processing repo to check those too.")
        return 0
    if not quiet:
        print("\nAll provenance manifests are consistent.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
