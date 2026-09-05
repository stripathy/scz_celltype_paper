#!/usr/bin/env python3
"""Provenance for the Figure 4 panel CSVs in results/figures/r_panels/.

The R renderers (scz_sst_hcn1_story.R, fig4_new_panels.R) draw from those CSVs,
not from the analyses that produced them. That is what makes the figure
reproducible, and also what lets it go quietly stale: rerun an upstream table
and the figure keeps drawing the old numbers. This module records, for every
panel CSV, the upstream file it was derived from and that file's checksum, into
r_panels/MANIFEST.tsv. shared/figure_inputs.R reads the manifest and stops the
renderer when a source has moved on.

Five exporters write into r_panels/ (export_panels_abc, export_panel_d_genetrack,
export_panels_dehi, export_panels_ef, export_panel_ad_concordance), so the map
for all of them lives here — one map, one manifest, no chance of one exporter
dropping another's rows.

Run standalone to refresh provenance without re-exporting anything:
    python3 scripts/figures/r_panels_provenance.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

# Repo-relative, so the manifest can be refreshed from any clone.
PAPER = Path(__file__).resolve().parents[3]
GEN = PAPER / "genetics"
OUT = GEN / "results" / "figures" / "r_panels"

TABLES = GEN / "results" / "tables"
REFGENE = GEN / "data" / "gwas" / "ncbiRefSeq_hg38.txt.gz"
FINEMAP = GEN / "data" / "fine_mapping" / "pgc3_finemap_credible_sets.csv"

PATCHSEQ = GEN / "data" / "patchseq"
# Five exemplars since the 2026-08 rebuild: three depleted supertypes and two
# not-depleted, ordered by soma depth. Cells not in the upstream cache were
# pulled from DANDI 000636 / the patch-seq repo into data/patchseq so the chain
# no longer depends on anything outside this repo. 2026-08-31: Sst_20
# (1079568285) replaced the Sst_22 exemplar (907585117).
EXEMPLARS = ("1079568285", "819770858", "1037461069", "758996755", "797048104")
SWC = [PATCHSEQ / "swc" / f"{i}_upright.swc" for i in EXEMPLARS]
NWB = [PATCHSEQ / "nwb" / f"{i}.nwb" for i in EXEMPLARS]

XEN = Path(os.environ.get("XENIUM_BASE", os.path.expanduser("~/Github/SCZ_Xenium")))
REF_H5AD = Path(os.environ.get("SEAAD_MTG_H5AD",
    os.path.expanduser("~/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad")))
SCZ_CRUMBLR = PAPER / "shared/snrnaseq_de/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv"
AD_CRUMBLR = PAPER / "crossdisorder/results/crumblr_results_supertype_neurons.csv"

# Enrichment now comes from MAGMA run natively on the SEA-AD DLPFC taxonomy
# alone (125 supertypes; build_spec_seaad_only.py -- the combined SEA-AD +
# Siletti taxonomy was retired on L. Duncan's advice); the specificity matrix
# is regenerable, so the gsa output and the gene-level results stand in for it.
GSA = GEN / "results" / "intermediates" / "T_a9only_bigdeli.gsa.out"
BIGDELI_GENES = GEN / "data" / "gwas" / "magma_bigdeli" / "bigdeli.step2.genes.out"
BIGDELI_SS = GEN / "data" / "gwas" / "bigdeli_eur_scz_sum_stats.gz"
A9 = sorted(Path(os.environ.get("SEAAD_A9_DIR",
    os.path.expanduser("~/Downloads"))).glob(
    "*A9_RNAseq_final-nuclei.2024-02-13.h5ad"))
COMPOS = TABLES / "gwas_vs_casecontrol_composition.csv"
EPHYS = TABLES / "sst_supertype_ephys_summary.csv"

# panel CSV -> upstream file(s) it was derived from
SOURCES: dict[str, object] = {
    # --- export_panels_abc.py / export_panel_d_genetrack.py / export_panels_ef.py ---
    "panel_B_genetics_vs_depletion.csv": [GSA, COMPOS],
    "panel_C_gene_drivers.csv":          [GSA, BIGDELI_GENES],
    "panel_C_top_labels.csv":            [GSA, BIGDELI_GENES],
    "panel_C_meta.csv":                  [GSA, BIGDELI_GENES],
    "panel_D_credible_set.csv":          BIGDELI_SS,
    "panel_D_snps.csv":                  BIGDELI_SS,
    "panel_D_meta.csv":                  BIGDELI_SS,
    "panel_D_genes.csv":                 REFGENE,
    "panel_D_exons.csv":                 REFGENE,
    "panel_E_hcn1_vs_sag.csv":           EPHYS,
    "panel_F_morphology.csv":            SWC,
    "panel_F_meta.csv":                  SWC,
    "panel_G_traces.csv":                NWB,
    "panel_G_meta.csv":                  NWB,
    # --- export_panels_dehi.py / export_panel_ad_concordance.py ---
    "panel_volcano_vulnerable_vs_notdepleted.csv": A9,
    "panel_violin_calb1_percell.csv":             A9,
    "panel_violin_donor_means.csv":               A9,
    "panel_violin_stats.csv":                     A9,
    "panel_ad_concordance_sst.csv":               [SCZ_CRUMBLR, AD_CRUMBLR],
}

DESCRIPTIONS = {
    "panel_B_genetics_vs_depletion.csv": "SCZ GWAS enrichment vs compositional depletion",
    "panel_C_gene_drivers.csv":          "per-gene drivers of Sst_2 enrichment",
    "panel_D_credible_set.csv":          "Bigdeli SuSiE-R EUR credible set at the HCN1 locus (Suppl. Table 13)",
    "panel_D_genes.csv":                 "RefSeq gene track for the HCN1 locus",
    "panel_E_hcn1_vs_sag.csv":           "HCN1 expression vs patch-seq sag per supertype",
    "panel_F_morphology.csv":            "SWC reconstruction nodes, five depth-ordered Sst exemplars",
    "panel_G_traces.csv":                "hyperpolarising-step voltage traces, same five cells",
    "panel_volcano_vulnerable_vs_notdepleted.csv": "vulnerable vs not-depleted Sst marker DE",
    "panel_violin_calb1_percell.csv":    "per-cell CALB1 / SST expression",
    "panel_violin_stats.csv":            "donor-paired t-test per gene",
    "panel_ad_concordance_sst.csv":      "SCZ compositional beta vs SEA-AD CPS slope, Sst",
}

# Any panel CSV absent from SOURCES is reported by record() rather than silently
# skipped — an unverifiable manifest row would be worse than an admitted gap.


def record(out_dir=OUT) -> Path:
    """Write r_panels/MANIFEST.tsv for whichever panel CSVs currently exist."""
    sys.path.insert(0, str(PAPER / "shared"))
    from figure_inputs import write_manifest
    out_dir = Path(out_dir)
    present = {k: v for k, v in SOURCES.items() if (out_dir / k).exists()}
    missing = sorted(set(SOURCES) - set(present))
    path = write_manifest(out_dir, present, DESCRIPTIONS)
    print(f"  recorded {len(present)} panel CSVs -> {path}")
    if missing:
        print(f"  not present, so not recorded: {', '.join(missing)}")
    unmapped = sorted(p.name for p in out_dir.glob("*.csv") if p.name not in SOURCES)
    if unmapped:
        print(f"  NOTE: {len(unmapped)} panel CSVs have no recorded source "
              f"(add them to SOURCES): {', '.join(unmapped)}")
    return path


if __name__ == "__main__":
    record()
