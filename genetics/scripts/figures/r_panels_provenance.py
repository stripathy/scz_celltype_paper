#!/usr/bin/env python3
"""Provenance for the Figure 4 panel CSVs in results/figures/r_panels/.

The R renderers (scz_sst_hcn1_story.R, fig4_new_panels.R) draw from those CSVs,
not from the analyses that produced them. That is what makes the figure
reproducible, and also what lets it go quietly stale: rerun an upstream table
and the figure keeps drawing the old numbers. This module records, for every
panel CSV, the upstream file it was derived from and that file's checksum, into
r_panels/MANIFEST.tsv. shared/figure_inputs.R reads the manifest and stops the
renderer when a source has moved on.

Two exporters write into r_panels/ (export_for_R.py and
export_fig4_new_panels.py), so the map for both lives here — one map, one
manifest, no chance of one exporter dropping the other's rows.

Run standalone to refresh provenance without re-exporting anything:
    python3 scripts/figures/r_panels_provenance.py
"""
from __future__ import annotations

import sys
from pathlib import Path

PAPER = Path("/Users/shreejoy/Github/scz_celltype_paper")
GEN = PAPER / "genetics"
OUT = GEN / "results" / "figures" / "r_panels"

TABLES = GEN / "results" / "tables"
INTERM = Path("/Users/shreejoy/Github/scz_cell_type_enrichment/results/intermediates")
REFGENE = GEN / "data" / "gwas" / "ncbiRefSeq_hg38.txt.gz"
FINEMAP = GEN / "data" / "fine_mapping" / "pgc3_finemap_credible_sets.csv"

PATCHSEQ = GEN / "data" / "patchseq"
SWC = [PATCHSEQ / "swc" / f"{i}_upright.swc" for i in ("819770858", "758996755")]
NWB = [PATCHSEQ / "nwb" / f"{i}.nwb" for i in ("819770858", "758996755")]

XEN = Path("/Users/shreejoy/Github/SCZ_Xenium")
REF_H5AD = Path("/Users/shreejoy/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad")
SCZ_CRUMBLR = PAPER / "spatial/data/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv"
AD_CRUMBLR = PAPER / "crossdisorder/results/crumblr_results_supertype_neurons.csv"

ENRICH = TABLES / "rbh_combined_enrichment.csv"
COMPOS = TABLES / "gwas_vs_casecontrol_composition.csv"
EPHYS = TABLES / "sst_supertype_ephys_summary.csv"
SPEC = INTERM / "rbh_combined_specificity.csv"

# panel CSV -> upstream file(s) it was derived from
SOURCES: dict[str, object] = {
    # --- export_for_R.py ---
    "panel_A_enrichment.csv":            ENRICH,
    "panel_A_family_groups.csv":         ENRICH,
    "panel_A_thresholds.csv":            ENRICH,
    "panel_B_genetics_vs_depletion.csv": COMPOS,
    "panel_B_bd.csv":                    COMPOS,
    "panel_B_bigdeli.csv":               COMPOS,
    "panel_C_gene_drivers.csv":          SPEC,
    "panel_C_top_labels.csv":            SPEC,
    "panel_C_meta.csv":                  SPEC,
    "panel_D_credible_set.csv":          FINEMAP,
    "panel_D_snps.csv":                  FINEMAP,
    "panel_D_meta.csv":                  FINEMAP,
    "panel_D_genes.csv":                 REFGENE,
    "panel_D_exons.csv":                 REFGENE,
    "panel_E_hcn1_vs_sag.csv":           EPHYS,
    "panel_F_morphology.csv":            SWC,
    "panel_F_meta.csv":                  SWC,
    "panel_G_traces.csv":                NWB,
    "panel_G_meta.csv":                  NWB,
    # --- export_fig4_new_panels.py ---
    "panel_volcano_vulnerable_vs_notdepleted.csv": REF_H5AD,
    "panel_violin_calb1_percell.csv":             REF_H5AD,
    "panel_violin_donor_means.csv":               REF_H5AD,
    "panel_violin_stats.csv":                     REF_H5AD,
    "panel_ad_concordance_sst.csv":               [SCZ_CRUMBLR, AD_CRUMBLR],
    "panel_ad_concordance_stats.csv":             [SCZ_CRUMBLR, AD_CRUMBLR],
}

DESCRIPTIONS = {
    "panel_A_enrichment.csv":            "MAGMA gene-property enrichment per supertype",
    "panel_B_genetics_vs_depletion.csv": "SCZ GWAS enrichment vs compositional depletion",
    "panel_C_gene_drivers.csv":          "per-gene drivers of Sst_25 enrichment",
    "panel_D_credible_set.csv":          "PGC3 FINEMAP credible set at the HCN1 locus",
    "panel_D_genes.csv":                 "RefSeq gene track for the HCN1 locus",
    "panel_E_hcn1_vs_sag.csv":           "HCN1 expression vs patch-seq sag per supertype",
    "panel_F_morphology.csv":            "SWC reconstruction nodes, Sst_25 and Sst_5 exemplars",
    "panel_G_traces.csv":                "hyperpolarising-step voltage traces, same two cells",
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
