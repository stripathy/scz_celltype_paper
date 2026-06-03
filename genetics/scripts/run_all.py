#!/usr/bin/env python
"""
run_all.py — Run the complete SCZ cell-type enrichment pipeline.

Executes all analysis steps in order:
  01: Compute cell-type specificity from SEA-AD snRNA-seq (~5-10 min, cached)
  02: MAGMA-style gene property analysis (~30s)
  03: Conditional analysis / forward selection (~30s)
  04: Gene driver identification (~30s)
  05: Spatial layer-stratified SST analysis (~30s)
  06: Gene-electrophysiology correlations (~1-2 min)
  07: Regenerate all figures (~30s)
  08: Siletti enrichment (Conti vs reprocessed, ~5s)
  09: MetaNeighbor integration (SEA-AD ↔ Siletti, ~10s)
  10: RBH combined taxonomy (~3s)
  11: Gene driver scatter plots (~30s)
  12: ATAC-seq peak overlap (~5min, needs h5ad)
  13: GWAS vs case-control composition (~1s)

Usage:
    python scripts/run_all.py                # Steps 01-07 (SEA-AD core pipeline)
    python scripts/run_all.py --all          # Steps 01-13 (includes Siletti, ATAC, composition)
    python scripts/run_all.py --force        # Recompute everything including specificity
    python scripts/run_all.py --from 3       # Start from step 3
    python scripts/run_all.py --only 5       # Run only step 5
"""
import sys
import time
import subprocess
from pathlib import Path


SCRIPTS_DIR = Path(__file__).parent
CORE_STEPS = [
    ("01_compute_specificity.py", "Compute cell-type specificity"),
    ("02_magma_enrichment.py", "MAGMA-style gene property analysis"),
    ("03_conditional_analysis.py", "Conditional analysis / forward selection"),
    ("04_gene_drivers.py", "Gene driver identification"),
    ("05_spatial_layer_analysis.py", "Spatial layer-stratified SST analysis"),
    ("06_gene_ephys_correlations.py", "Gene-electrophysiology correlations"),
    ("07_make_all_figures.py", "Regenerate all figures"),
]

SILETTI_STEPS = [
    ("08_siletti_enrichment.py", "Siletti enrichment (Conti vs reprocessed)"),
    ("09_metaneighbor_integration.py", "MetaNeighbor integration (SEA-AD ↔ Siletti)"),
    ("10_combined_taxonomy.py", "RBH combined taxonomy"),
    ("11_gene_driver_scatter.py", "Gene driver scatter plots (top 10 types)"),
    ("12_atac_peak_overlap.py", "ATAC-seq peak overlap with FINEMAP credible sets"),
    ("13_gwas_vs_composition.py", "GWAS vs case-control composition"),
]


def run_step(script_name, description, extra_args=None):
    """Run a single pipeline step as a subprocess."""
    script_path = SCRIPTS_DIR / script_name
    cmd = [sys.executable, str(script_path)]
    if extra_args:
        cmd.extend(extra_args)

    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Script:  {script_name}")
    print(f"{'='*60}")

    t0 = time.time()
    result = subprocess.run(cmd, cwd=str(SCRIPTS_DIR.parent))

    elapsed = time.time() - t0
    if result.returncode != 0:
        print(f"\nERROR: {script_name} failed with return code {result.returncode}")
        print(f"Stopping pipeline.")
        sys.exit(result.returncode)

    print(f"\nCompleted {script_name} in {elapsed:.1f}s")
    return elapsed


def main():
    args = sys.argv[1:]

    # Parse arguments
    force = "--force" in args
    run_all_steps = "--all" in args
    start_from = 1
    only_step = None

    for i, arg in enumerate(args):
        if arg == "--from" and i + 1 < len(args):
            start_from = int(args[i + 1])
        if arg == "--only" and i + 1 < len(args):
            only_step = int(args[i + 1])

    STEPS = CORE_STEPS + (SILETTI_STEPS if run_all_steps or start_from > 7 or (only_step and only_step > 7) else [])

    print("=" * 60)
    print("SCZ Cell-Type Enrichment Pipeline")
    print("=" * 60)
    print(f"  Force recompute: {force}")
    print(f"  Start from step: {start_from}")
    print(f"  Include extended steps (08-13): {run_all_steps or start_from > 7}")
    if only_step:
        print(f"  Run only step:   {only_step}")
    print()

    total_t0 = time.time()
    timings = []

    for step_num, (script, desc) in enumerate(STEPS, start=1):
        if only_step and step_num != only_step:
            continue
        if step_num < start_from:
            print(f"  Skipping step {step_num}: {desc}")
            continue

        extra_args = ["--force"] if force and step_num == 1 else None
        elapsed = run_step(script, desc, extra_args=extra_args)
        timings.append((step_num, desc, elapsed))

    total_elapsed = time.time() - total_t0

    print(f"\n{'='*60}")
    print(f"PIPELINE COMPLETE")
    print(f"{'='*60}")
    print(f"\nTiming summary:")
    for step_num, desc, elapsed in timings:
        print(f"  Step {step_num:02d}: {elapsed:7.1f}s  {desc}")
    print(f"  {'─'*40}")
    print(f"  Total:  {total_elapsed:7.1f}s")

    print(f"\nResults saved to:")
    print(f"  Tables:  results/tables/")
    print(f"  Figures: results/figures/")


if __name__ == "__main__":
    main()
