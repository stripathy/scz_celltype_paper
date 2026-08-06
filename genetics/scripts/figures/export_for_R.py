#!/usr/bin/env python3
"""
Export every data table the R figure script needs.

The R/ggplot script (`scz_sst_hcn1_story.R`) is the renderer; this Python
script does the heavy lifting (HDF5 reads, SWC parsing, NWB reads, liftover)
and emits flat CSVs the R script can ingest with read_csv().

Outputs land in `results/figures/r_panels/`:

  panel_A_enrichment.csv       137 SEA-AD supertypes, β + FDR + Bonferroni
  panel_A_family_groups.csv    family colorbar coordinates
  panel_B_genetics_vs_depletion.csv  16 Sst, β vs comp_beta, depleted flag
  panel_C_gene_drivers.csv     ~17K genes, specificity + -log10p in Sst_25
  panel_C_top_labels.csv       Top 15 driver labels with positions
  panel_D_snps.csv             GWAS SNPs in HCN1 hg38 window
  panel_D_genes.csv            RefSeq gene structures, hg38, long form (one row per exon)
  panel_D_meta.csv             Lead SNP info, HCN1 coords, window
  panel_E_hcn1_vs_sag.csv      16 Sst, HCN1 expression + mean sag
  panel_F_morphology.csv       SWC nodes for both example cells (long form)
  panel_F_meta.csv             Per-cell pia_dist, color, label
  panel_G_traces.csv           Voltage waveforms for both cells (long form)
"""
import os, sys, json, gzip
from pathlib import Path
import numpy as np
import pandas as pd

# Monorepo (paper) root — committed light inputs + panel-CSV outputs live here.
REPO = Path("/Users/shreejoy/Github/scz_celltype_paper/genetics")
DATA = REPO / "data"
TABLES = REPO / "results" / "tables"
OUT = REPO / "results" / "figures" / "r_panels"
OUT.mkdir(parents=True, exist_ok=True)

# ── External raw sources (LARGE / cross-repo) — referenced IN PLACE, deliberately
# NOT vendored into the paper repo (~325 MB + files in a third repo). Panels D/E
# and F/G read these; regenerating them requires these paths to exist locally.
# Panels A/B (and the two shipped figures) do NOT depend on this block.
# See scripts/figures/PORT_NOTES.md for the full provenance table.
ENRICH_REPO   = Path("/Users/shreejoy/Github/scz_cell_type_enrichment")
INTERM        = ENRICH_REPO / "results" / "intermediates"   # 43M mean-expr (panel E); 287M specificity (panel C, unused)
RAW_GWAS_DIR  = ENRICH_REPO / "data" / "gwas"               # 229M PGC3 sumstats (panel D locus)
MAGMA_REPO    = ENRICH_REPO / "linking_cell_types_to_brain_phenotypes"  # MAGMA out + gene.loc (panel C, unused)
PATCHSEQ_REPO = Path("/Users/shreejoy/Github/human_int_patch_seq")      # patchseq_builder code (parse_swc, orientation)
# The four raw files panels F/G need (2 SWC + 2 NWB, 54 MB) are held in the
# repo so those panels reproduce without the patch-seq repo checked out.
# Falls back to the upstream repo if the in-repo copies are absent.
PATCHSEQ_DATA = DATA / "patchseq"
if not PATCHSEQ_DATA.exists():
    PATCHSEQ_DATA = PATCHSEQ_REPO / "data"                              # swc/ and nwb/ differ upstream; see data/patchseq/README.md

# ────────────────────────────────────────────────────────────────────
# Common: SEA-AD supertype colors
# ────────────────────────────────────────────────────────────────────
with open(REPO / "data" / "seaad_supertype_colors.json") as f:
    SEAAD_COLORS = json.load(f)

# Family colors (for fallback + colorbar)
FAMILY_COLORS = {
    "Sst": "#FF9900", "Sst Chodl": "#B1B10C", "Pvalb": "#D93137",
    "Chandelier": "#F641A8", "Vip": "#A45FBF", "Sncg": "#DF70FF",
    "Pax6": "#71238C", "Lamp5": "#DA808C", "Lamp5 Lhx6": "#935F50",
    "L2/3 IT": "#B1EC30", "L4 IT": "#00E5E5", "L5 IT": "#50B2AD",
    "L6 IT": "#A19922", "L6 IT Car3": "#5100FF",
    "L5 ET": "#0D5B78", "L5/6 NP": "#3E9E64",
    "L6 CT": "#2D8CB8", "L6b": "#7044AA",
    "Astrocyte": "#665C47", "Oligodendrocyte": "#53776C",
    "OPC": "#374A45", "Microglia": "#94AF97",
    "Endothelial": "#8D6C62", "VLMC": "#697255",
    "Pericyte": "#A6678A", "Vsmc": "#7C4F8C",
    "Other": "#888888",
}
def family_of(s):
    aliases = {
        "Astro": "Astrocyte", "Oligo": "Oligodendrocyte",
        "Endo": "Endothelial", "Micro-PVM": "Microglia",
        "SMC": "Vsmc", "Lymphocyte": "Other", "Monocyte": "Other",
    }
    for fam in ["Sst Chodl", "L2/3 IT", "L4 IT", "L5 IT", "L6 IT Car3", "L6 IT",
                 "L5 ET", "L5/6 NP", "L6 CT", "L6b",
                 "Lamp5 Lhx6", "Lamp5", "Sst", "Pvalb", "Vip", "Sncg",
                 "Pax6", "Chandelier", "Astro", "Oligo", "OPC", "Micro-PVM",
                 "Endo", "VLMC", "SMC", "Pericyte", "Lymphocyte", "Monocyte"]:
        if s.startswith(fam):
            return aliases.get(fam, fam if fam in FAMILY_COLORS else "Other")
    return "Other"


# ────────────────────────────────────────────────────────────────────
# Panel A: SCZ enrichment landscape
# ────────────────────────────────────────────────────────────────────
print("Panel A …")
enrich = pd.read_csv(TABLES / "rbh_combined_enrichment.csv")
sea = enrich[enrich.source == "SEA-AD"].copy().reset_index(drop=True)
sea["family"] = sea.supertype.apply(family_of)
fam_order = [
    "Sst", "Sst Chodl", "Pvalb", "Chandelier", "Vip", "Sncg", "Lamp5",
    "Pax6", "Lamp5 Lhx6", "L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3",
    "L5 ET", "L5/6 NP", "L6 CT", "L6b",
    "Astrocyte", "Oligodendrocyte", "OPC", "Microglia",
    "Endothelial", "VLMC", "Vsmc", "Pericyte", "Other",
]
sea["family"] = pd.Categorical(sea.family, categories=fam_order, ordered=True)
# Sort within-family by p-value (most-significant first), matching the
# Duncan-style -log10(P) y-axis convention so bars decay left-to-right within
# each family group on the figure.
sea = sea.sort_values(["family", "p_value"], ascending=[True, True]).reset_index(drop=True)
sea["x"] = np.arange(len(sea))

BONF_P_THRESH = 0.05 / len(enrich)
sea["passes_bonf_503"] = sea.p_value < BONF_P_THRESH
sea["passes_fdr_05"] = sea.p_fdr < 0.05
sea["color"] = sea.supertype.map(SEAAD_COLORS).fillna(
    sea.family.astype(str).map(FAMILY_COLORS)).fillna("#888888")

sea[["supertype", "family", "x", "beta", "p_value", "p_fdr",
      "passes_bonf_503", "passes_fdr_05", "color"]].to_csv(
        OUT / "panel_A_enrichment.csv", index=False)

# Family groups for the colorbar
fam_groups = (sea.groupby("family", observed=True)
                 .x.agg(["min", "max", "size"]).reset_index())
fam_groups["color"] = fam_groups.family.astype(str).map(FAMILY_COLORS).fillna("#888888")
fam_groups["family"] = fam_groups.family.astype(str)
fam_groups.to_csv(OUT / "panel_A_family_groups.csv", index=False)

# Thresholds for h-lines
thresholds = pd.DataFrame([
    dict(label="FDR 0.05 (Franken 503)", beta=sea.loc[sea.passes_fdr_05, "beta"].min()),
    dict(label="Bonferroni 0.05/503", beta=sea.loc[sea.passes_bonf_503, "beta"].min()),
])
thresholds.to_csv(OUT / "panel_A_thresholds.csv", index=False)
print(f"  → {len(sea)} bars; Bonferroni-sig {sea.passes_bonf_503.sum()}")


# ────────────────────────────────────────────────────────────────────
# Panel B (was F): Genetic risk vs SCZ cell-abundance depletion
# ────────────────────────────────────────────────────────────────────
print("Panel B …")
comp = pd.read_csv(TABLES / "gwas_vs_casecontrol_composition.csv")
sst = (comp[comp.supertype.str.startswith("Sst")
            & ~comp.supertype.str.startswith("Sst Chodl")].copy())
sst["depleted_fdr20"] = (sst.comp_padj < 0.20) & (sst.comp_beta < 0)
sst["exemplar"] = sst.supertype.isin(["Sst_25", "Sst_5"])
sst["pt_size"] = 1.0 + 3.0 * np.clip(sst.comp_logp.values / 5, 0, 1)
sst["color"] = sst.supertype.map(SEAAD_COLORS).fillna("#888888")
# Add SCZ enrichment p-value / -log10(p) so Panel B can switch x-axis from
# β to -log10(P) to match Panel A's Duncan-style convention.
scz_p_lookup = dict(zip(enrich.supertype, enrich.p_value))
sst["scz_p"] = sst.supertype.map(scz_p_lookup)
sst["scz_neg_log10_p"] = -np.log10(sst.scz_p.clip(lower=1e-300))
sst[["supertype", "beta", "scz_p", "scz_neg_log10_p",
      "comp_beta", "comp_padj", "comp_logp",
      "depleted_fdr20", "exemplar", "pt_size", "color"]].to_csv(
        OUT / "panel_B_genetics_vs_depletion.csv", index=False)
print(f"  → {len(sst)} Sst types; {sst.depleted_fdr20.sum()} depleted at FDR<0.20")


# ────────────────────────────────────────────────────────────────────
# Panel C (was B): Gene drivers in Sst_25
# ────────────────────────────────────────────────────────────────────
print("Panel C …")
spec = pd.read_csv(INTERM / "rbh_combined_specificity.csv", index_col=0)
magma = pd.read_csv(MAGMA_REPO / "Example_results/"
                          "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.no_heading.step2.genes.out",
                    sep=r"\s+")
# Use the SAME MHC-excluded gene loc file the upstream MAGMA enrichment uses
# (Duncan/Trubetskoy convention). The MHC's long-range LD makes gene-level
# signals across ~190 genes simultaneously, including olfactory receptors with
# noise-driven 100% specificity. The non-MHC file has 19,175 genes.
gl = pd.read_csv(MAGMA_REPO / "Data/NCBI37.3.gene.loc.extendedMHCexcluded",
                 sep="\t", header=None,
                 names=["ENTREZ", "CHR", "START", "STOP", "STRAND", "SYMBOL"])
magma["symbol"] = magma.GENE.astype(int).map(dict(zip(gl.ENTREZ.astype(int), gl.SYMBOL)))

target = "Sst_25"
s = spec[target].rename("specificity")
m = (pd.merge(s.reset_index().rename(columns={"index": "symbol"}),
              magma[["symbol", "P"]].dropna(), on="symbol", how="inner"))

# Note: MHC exclusion is already handled by using the MHC-excluded gene-loc
# file above (matches the upstream Duncan/Trubetskoy MAGMA convention).
# So OR12D3, OR5V1, HLA-*, etc. never enter the symbol-mapped table.
m["neg_log10_p"] = -np.log10(m.P.clip(lower=1e-300))
from statsmodels.stats.multitest import multipletests
m["p_fdr"] = multipletests(m.P, method="fdr_bh")[1]
spec_q90 = m.specificity.quantile(0.90)
m["is_driver"] = (m.specificity >= spec_q90) & (m.p_fdr < 0.05)
m["is_hcn1"] = m.symbol == "HCN1"
m.to_csv(OUT / "panel_C_gene_drivers.csv", index=False)

# Top driver labels by RANK-SUM of (specificity, −log10(p)), with locus
# deduplication. Steps:
#   1. For each driver gene, compute rank by specificity (desc) and by
#      −log10(p) (desc); combined_rank = rank_spec + rank_p.
#   2. Annotate each gene with chr + midpoint from the gene-loc table.
#   3. Greedy locus-dedup loop: sort by combined_rank ascending; iterate
#      and accept a gene only if no already-accepted gene is within 1 Mb
#      on the same chromosome. Stop at 10 genes.
#   4. Always include HCN1 (highlighted in blue) even if it doesn't make
#      the deduplicated top 10.
LOCUS_DEDUP_KB = 1000   # 1 Mb window per labeled locus
drv = m[m.is_driver].copy()
drv["rank_spec"] = drv.specificity.rank(ascending=False, method="min")
drv["rank_p"]    = drv.neg_log10_p.rank(ascending=False, method="min")
drv["combined_rank"] = drv.rank_spec + drv.rank_p
# Bring in chr/start/stop from the (MHC-excluded) gene-loc table
gene_coords = gl.set_index("SYMBOL")[["CHR", "START", "STOP"]]
drv = drv.join(gene_coords, on="symbol", how="left")
drv["mid_bp"] = (drv.START + drv.STOP) / 2

# Greedy locus-deduped selection
drv_sorted = drv.sort_values("combined_rank").reset_index(drop=True)
accepted_loci = []   # list of (chr, mid_bp) for already-picked genes
labels_keep = []
for _, row in drv_sorted.iterrows():
    if pd.isna(row.CHR) or pd.isna(row.mid_bp):
        continue
    in_used_window = any(
        ch == row.CHR and abs(mb - row.mid_bp) < LOCUS_DEDUP_KB * 1000
        for ch, mb in accepted_loci
    )
    if in_used_window:
        continue
    accepted_loci.append((row.CHR, row.mid_bp))
    labels_keep.append(row.symbol)
    if len(labels_keep) >= 10:
        break

top_labels = drv[drv.symbol.isin(labels_keep)].copy()
top_labels["is_hcn1"] = top_labels.symbol == "HCN1"
# Always include HCN1 even if it didn't survive locus-dedup
if "HCN1" not in top_labels.symbol.values and "HCN1" in drv.symbol.values:
    hcn1_row = drv[drv.symbol == "HCN1"].copy()
    hcn1_row["is_hcn1"] = True
    top_labels = pd.concat([top_labels, hcn1_row], ignore_index=True)
    print(f"    HCN1 added explicitly outside top-10 "
          f"(rank-sum = {hcn1_row.combined_rank.iloc[0]:.0f})")
top_labels.to_csv(OUT / "panel_C_top_labels.csv", index=False)
print(f"    Top-10 driver labels by rank-sum (1 Mb locus-deduped): "
      f"{', '.join(labels_keep)}")

# Constants for axis truncation
# FDR 0.05 threshold in -log10(p) space: the y-value of the largest p that
# still passes BH FDR < 0.05. Used as the horizontal threshold line in
# Panel C (this is what defines drivers, NOT the conventional 5×10⁻⁸).
fdr_passing = m[m.p_fdr < 0.05]
if len(fdr_passing):
    fdr_p_cutoff = fdr_passing.P.max()
    fdr_neg_log10p_cutoff = -np.log10(fdr_p_cutoff)
else:
    fdr_neg_log10p_cutoff = float("nan")

pd.DataFrame([dict(spec_q90=spec_q90,
                    xlim_lo=spec_q90 * 0.4,
                    xlim_hi=0.05,
                    ylim_lo=1.3,
                    ylim_hi=20.0,   # fixed max so Panel C doesn't have dead space above HCN1
                    fdr_neg_log10p_cutoff=fdr_neg_log10p_cutoff,
                    target_type=target,
                    target_beta=enrich.loc[enrich.supertype == target, "beta"].iloc[0],
                    target_p_fdr=enrich.loc[enrich.supertype == target, "p_fdr"].iloc[0])
              ]).to_csv(OUT / "panel_C_meta.csv", index=False)
print(f"  → {len(m)} genes; {m.is_driver.sum()} drivers; HCN1 in driver set: "
      f"{m.loc[m.is_hcn1, 'is_driver'].iloc[0]}")


# ────────────────────────────────────────────────────────────────────
# Panel D: HCN1 locus zoom (Manhattan + gene track, hg38)
# ────────────────────────────────────────────────────────────────────
print("Panel D …")
from liftover import get_lifter
lifter = get_lifter("hg19", "hg38")

HCN1_START_HG38, HCN1_STOP_HG38 = 45_254_947, 45_696_380
# Asymmetric window so all 9 FINEMAP credible-set variants fit in frame.
# The most-downstream variant (rs10051592, PIP 0.002) sits ~165 kb 3' of
# HCN1, so the left pad is 200 kb. No credible-set variants extend past
# HCN1's 5' end (right side), so the right pad is just 50 kb.
WIN_PAD_LEFT  = 200_000
WIN_PAD_RIGHT = 50_000
win_lo = HCN1_START_HG38 - WIN_PAD_LEFT
win_hi = HCN1_STOP_HG38  + WIN_PAD_RIGHT

# GWAS SNPs (lift hg19 → hg38)
gwas_path = RAW_GWAS_DIR / "PGC3_SCZ_wave3.european.autosome.public.v3.vcf (1).tsv.gz"
rows = []
with gzip.open(gwas_path, "rt") as f:
    header = None
    for line in f:
        if line.startswith("##"): continue
        if header is None:
            header = line.rstrip("\n").split("\t")
            idx_chrom = header.index("CHROM")
            idx_pos = header.index("POS")
            idx_pval = header.index("PVAL")
            idx_id = header.index("ID")
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(idx_chrom, idx_pos, idx_pval): continue
        if p[idx_chrom] != "5": continue
        try:
            pos19 = int(p[idx_pos])
        except ValueError:
            continue
        if pos19 < win_lo - 50_000 or pos19 > win_hi + 50_000: continue
        try:
            pval = float(p[idx_pval])
        except ValueError:
            continue
        rows.append((pos19, pval, p[idx_id]))
df = pd.DataFrame(rows, columns=["pos_hg19", "p", "rsid"])
def _lift(pos):
    try:
        r = lifter.convert_coordinate("chr5", int(pos))
        return int(r[0][1]) if r else None
    except Exception:
        return None
df["pos_hg38"] = [_lift(p) for p in df.pos_hg19]
df = df.dropna(subset=["pos_hg38"])
df["pos_hg38"] = df.pos_hg38.astype(int)
df = df[(df.pos_hg38 >= win_lo) & (df.pos_hg38 <= win_hi)]
df["neg_log10_p"] = -np.log10(df.p.clip(lower=1e-300))
df["is_gwas_sig"] = df.p < 5e-8
df["pos_mb"] = df.pos_hg38 / 1e6
# ───── Fine-mapping (PGC3 FINEMAP credible sets) ─────
fm = pd.read_csv(REPO / "data" / "fine_mapping" / "pgc3_finemap_credible_sets.csv")
cs = fm[fm.locus == "rs10035564"].copy()  # the HCN1 credible set
# Lift credible-set positions hg19 → hg38
cs["pos_hg38"] = [_lift(p) for p in cs.pos]
cs = cs.dropna(subset=["pos_hg38"])
cs["pos_hg38"] = cs.pos_hg38.astype(int)
cs["pos_mb"] = cs.pos_hg38 / 1e6
cs["neg_log10_p"] = -np.log10(cs.gwas_p.clip(lower=1e-300))
cs["is_lead"] = cs.rsid == cs.loc[cs.pip.idxmax(), "rsid"]

# Classify each credible-set variant by HCN1 location (intron, UTR, downstream)
HCN1_EXON_STARTS = [45254947, 45267088, 45303598, 45353099, 45396491,
                     45461845, 45645184, 45695668]
HCN1_EXON_ENDS   = [45262810, 45267253, 45303839, 45353246, 45396710,
                     45462007, 45645608, 45696380]
HCN1_CDS_START, HCN1_CDS_END = 45261920, 45696093
def _classify_hcn1(pos38):
    if pos38 < HCN1_START_HG38:
        return f"3'_downstream"
    if pos38 > HCN1_STOP_HG38:
        return f"5'_upstream"
    for es, ee in zip(HCN1_EXON_STARTS, HCN1_EXON_ENDS):
        if es <= pos38 <= ee:
            if pos38 < HCN1_CDS_START: return "3'UTR"
            if pos38 > HCN1_CDS_END: return "5'UTR"
            return "CDS"
    return "intron"
cs["hcn1_location"] = cs.pos_hg38.apply(_classify_hcn1)
cs.to_csv(OUT / "panel_D_credible_set.csv", index=False)
print(f"    Fine-mapping credible set: {len(cs)} variants "
      f"(PIP sum = {cs.pip.sum():.3f})")
_inside_hcn1 = cs.hcn1_location.isin(["intron", "CDS", "5'UTR", "3'UTR"]).sum()
print(f"    Inside HCN1 body: {_inside_hcn1} of {len(cs)} variants")

# Tag each SNP with credible-set membership for the Manhattan render
cs_rsids = set(cs.rsid)
df["in_credible_set"] = df.rsid.isin(cs_rsids)
df["pip"] = df.rsid.map(dict(zip(cs.rsid, cs.pip))).fillna(0.0)
df[["rsid", "pos_hg38", "pos_mb", "p", "neg_log10_p", "is_gwas_sig",
     "in_credible_set", "pip"]].to_csv(OUT / "panel_D_snps.csv", index=False)

# Lead SNP + window meta
lead = df.loc[df.p.idxmin()]
# Second-best credible-set variant (the highest-PIP intronic one)
intronic = cs[cs.hcn1_location == "intron"].sort_values("pip", ascending=False)
second = intronic.iloc[0] if len(intronic) else cs.iloc[0]
meta_D = pd.DataFrame([dict(
    win_lo_hg38=win_lo, win_hi_hg38=win_hi,
    win_lo_mb=win_lo / 1e6, win_hi_mb=win_hi / 1e6,
    hcn1_start_mb=HCN1_START_HG38 / 1e6, hcn1_stop_mb=HCN1_STOP_HG38 / 1e6,
    lead_rsid=lead.rsid, lead_pos_mb=lead.pos_mb, lead_p=lead.p,
    lead_neglog10p=lead.neg_log10_p,
    lead_pip=cs.loc[cs.rsid == lead.rsid, "pip"].iloc[0] if lead.rsid in cs_rsids else np.nan,
    second_rsid=second.rsid, second_pos_mb=second.pos_mb,
    second_pip=second.pip, second_loc=second.hcn1_location,
    second_neglog10p=second.neg_log10_p,
    n_credible=len(cs), pip_sum=cs.pip.sum())])
meta_D.to_csv(OUT / "panel_D_meta.csv", index=False)

# Gene structures from hg38 RefSeq — emit one row per exon
import gzip as _gz
REFGENE = DATA / "gwas" / "ncbiRefSeq_hg38.txt.gz"
gene_rows = []
exon_rows = []
gene_seen = {}
with _gz.open(REFGENE, "rt") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 13: continue
        if parts[2] != "chr5": continue
        tx_start = int(parts[4])
        tx_end = int(parts[5])
        if tx_end < win_lo or tx_start > win_hi: continue
        symbol = parts[12]
        tx_len = tx_end - tx_start
        if symbol in gene_seen and gene_seen[symbol] >= tx_len: continue
        gene_seen[symbol] = tx_len
        # Re-emit clean rows
        cds_start = int(parts[6]); cds_end = int(parts[7])
        ex_s = [int(x) for x in parts[9].rstrip(",").split(",") if x]
        ex_e = [int(x) for x in parts[10].rstrip(",").split(",") if x]
        # Keep gene-level summary
        gene_rows.append(dict(
            symbol=symbol, strand=parts[3],
            start_mb=tx_start / 1e6, stop_mb=tx_end / 1e6,
            cds_start_mb=cds_start / 1e6, cds_end_mb=cds_end / 1e6,
            is_hcn1=(symbol == "HCN1"), tx_len=tx_len))
        # Per-exon rows
        for es, ee in zip(ex_s, ex_e):
            if ee < win_lo or es > win_hi: continue
            # Classify exon segment: 5'UTR / CDS / 3'UTR
            if cds_start == cds_end:
                # Non-coding
                exon_rows.append(dict(symbol=symbol, strand=parts[3],
                                       seg_type="UTR", is_hcn1=(symbol == "HCN1"),
                                       start_mb=max(es, win_lo) / 1e6,
                                       stop_mb=min(ee, win_hi) / 1e6))
            else:
                # 5' UTR
                if es < cds_start:
                    a = max(es, win_lo); b = min(ee, cds_start)
                    if b > a:
                        exon_rows.append(dict(symbol=symbol, strand=parts[3],
                                               seg_type="UTR",
                                               is_hcn1=(symbol == "HCN1"),
                                               start_mb=a / 1e6, stop_mb=b / 1e6))
                # CDS
                a = max(es, cds_start); b = min(ee, cds_end)
                if b > a:
                    exon_rows.append(dict(symbol=symbol, strand=parts[3],
                                           seg_type="CDS",
                                           is_hcn1=(symbol == "HCN1"),
                                           start_mb=max(a, win_lo) / 1e6,
                                           stop_mb=min(b, win_hi) / 1e6))
                # 3' UTR
                if ee > cds_end:
                    a = max(es, cds_end); b = min(ee, win_hi)
                    if b > a:
                        exon_rows.append(dict(symbol=symbol, strand=parts[3],
                                               seg_type="UTR",
                                               is_hcn1=(symbol == "HCN1"),
                                               start_mb=a / 1e6, stop_mb=b / 1e6))

# Lane assignment so genes don't overlap (sort by start, place in earliest free lane)
genes_df = pd.DataFrame(gene_rows).sort_values("start_mb").reset_index(drop=True)
lanes_end = []
genes_df["lane"] = 0
span_mb = (win_hi - win_lo) / 1e6
for i, g in genes_df.iterrows():
    label_extent = g.stop_mb + len(g.symbol) * span_mb * 0.008 + span_mb * 0.015
    placed = False
    for li, le in enumerate(lanes_end):
        if g.start_mb > le:
            genes_df.at[i, "lane"] = li
            lanes_end[li] = label_extent
            placed = True; break
    if not placed:
        genes_df.at[i, "lane"] = len(lanes_end)
        lanes_end.append(label_extent)
genes_df["n_lanes"] = max(genes_df.lane) + 1 if len(genes_df) else 1
genes_df.to_csv(OUT / "panel_D_genes.csv", index=False)

# Add lane info to each exon row
exon_df = pd.DataFrame(exon_rows)
lane_map = dict(zip(genes_df.symbol, genes_df.lane))
exon_df["lane"] = exon_df.symbol.map(lane_map)
exon_df.to_csv(OUT / "panel_D_exons.csv", index=False)
print(f"  → {len(df)} SNPs; {len(genes_df)} genes ({len(exon_df)} exon segments); lead = {lead.rsid}")


# ────────────────────────────────────────────────────────────────────
# Panel E (was C): HCN1 normalized expression vs sag, 16 Sst supertypes
# ────────────────────────────────────────────────────────────────────
print("Panel E …")
ephys = pd.read_csv(TABLES / "sst_supertype_ephys_summary.csv")
# Ported in-repo (Sst slice of the SEA-AD reference mean-expression matrix)
# so Figure 4 is reproducible without the external enrichment repo.
exp = pd.read_csv(REPO / "results" / "intermediates" /
                  "seaad_sst_supertype_mean_expression.csv", index_col=0)
hcn1_exp = exp.loc["HCN1"]
ephys["HCN1_expr"] = ephys.supertype.map(hcn1_exp.to_dict())
df_E = ephys.dropna(subset=["HCN1_expr", "mean_sag"]).copy()
# Add depletion status (FDR<0.20)
comp_lookup = comp.set_index("supertype")[["comp_padj", "comp_beta"]]
df_E["comp_padj"] = df_E.supertype.map(comp_lookup.comp_padj)
df_E["comp_beta"] = df_E.supertype.map(comp_lookup.comp_beta)
df_E["depleted_fdr20"] = (df_E.comp_padj < 0.20) & (df_E.comp_beta < 0)
df_E["exemplar"] = df_E.supertype.isin(["Sst_25", "Sst_5"])
df_E["color"] = df_E.supertype.map(SEAAD_COLORS).fillna("#888888")
df_E.to_csv(OUT / "panel_E_hcn1_vs_sag.csv", index=False)
print(f"  → {len(df_E)} Sst supertypes")


# ────────────────────────────────────────────────────────────────────
# Panel F (was M): Cell morphologies from SWC
# ────────────────────────────────────────────────────────────────────
print("Panel F …")
sys.path.insert(0, "/Users/shreejoy/Github/human_int_patch_seq")
from patchseq_builder.morphology.download import parse_swc
from patchseq_builder.morphology.orientation import INVERTED_SPECIMEN_IDS, flip_swc_y

MORPH_CELLS = [
    dict(specimen_id=819770858, supertype="Sst_25", layer="L2",
         swc_path=str(PATCHSEQ_DATA / "swc" / "819770858_upright.swc"),
         pia_dist_um=405.610339, color=SEAAD_COLORS["Sst_25"]),
    dict(specimen_id=758996755, supertype="Sst_5", layer="L4",
         swc_path=str(PATCHSEQ_DATA / "swc" / "758996755_upright.swc"),
         pia_dist_um=1500.0, color=SEAAD_COLORS["Sst_5"]),
]
all_segments = []
meta_F = []
for ci, cell in enumerate(MORPH_CELLS):
    nodes = parse_swc(cell["swc_path"])
    if cell["specimen_id"] in INVERTED_SPECIMEN_IDS:
        nodes = flip_swc_y(nodes)
    soma_xy = None
    for n in nodes.values():
        if n["type"] == 1:
            soma_xy = (n["x"], n["y"]); break
    if soma_xy is None:
        soma_xy = (np.mean([n["x"] for n in nodes.values()]),
                   np.mean([n["y"] for n in nodes.values()]))
    x_off = soma_xy[0]; y_soma = soma_xy[1]
    def to_um(n):
        return (n["x"] - x_off, cell["pia_dist_um"] + (y_soma - n["y"]))
    # Emit segments (one row per parent→child link)
    for nid, n in nodes.items():
        if n["parent"] < 0 or n["parent"] not in nodes: continue
        p = nodes[n["parent"]]
        px, py = to_um(p); cx, cy = to_um(n)
        all_segments.append(dict(cell=cell["supertype"],
                                  comp_type=n["type"],
                                  x0=px, y0=py, x1=cx, y1=cy))
    # Soma position
    sx, sy = to_um({"x": soma_xy[0], "y": soma_xy[1]})
    meta_F.append(dict(specimen_id=cell["specimen_id"],
                        supertype=cell["supertype"],
                        layer=cell["layer"],
                        pia_dist_um=cell["pia_dist_um"],
                        color=cell["color"],
                        soma_x_um=sx, soma_y_um=sy,
                        cell_order=ci))

pd.DataFrame(all_segments).to_csv(OUT / "panel_F_morphology.csv", index=False)
pd.DataFrame(meta_F).to_csv(OUT / "panel_F_meta.csv", index=False)
print(f"  → {len(all_segments)} morphology segments across {len(meta_F)} cells")


# ────────────────────────────────────────────────────────────────────
# Panel G (was E): Patch-seq voltage traces from NWB
# ────────────────────────────────────────────────────────────────────
print("Panel G …")
import pynwb
from pynwb import NWBHDF5IO

TRACES_CELLS = [
    dict(specimen_id=819770858, supertype="Sst_25", layer="L2",
         sag=0.564, color=SEAAD_COLORS["Sst_25"], target_pa=-100,
         nwb_path=str(PATCHSEQ_DATA / "nwb" / "819770858.nwb")),
    dict(specimen_id=758996755, supertype="Sst_5", layer="L4",
         sag=0.00196, color=SEAAD_COLORS["Sst_5"], target_pa=-30,
         nwb_path=str(PATCHSEQ_DATA / "nwb" / "758996755.nwb")),
]

def best_hyperpol_sweep(nwb_path, target_pa, amp_tol=25):
    candidates = []
    with NWBHDF5IO(nwb_path, "r", load_namespaces=True) as nio:
        nwb = nio.read()
        for k, v in nwb.acquisition.items():
            if not isinstance(v, pynwb.icephys.CurrentClampSeries): continue
            sd = getattr(v, "stimulus_description", "") or ""
            if not any(t in sd for t in ("X1PS_SubThresh", "X3LP_Rheo", "X4PS_SupraThresh")): continue
            s = nwb.stimulus.get(k.replace("_AD0", "_DA0"))
            if s is None: continue
            v_mv = np.asarray(v.data[:]) * v.conversion * 1000.0
            i_pa = np.asarray(s.data[:]) * s.conversion * 1e12
            rate = float(v.rate); t = np.arange(len(v_mv)) / rate
            mask = np.abs(i_pa) > 5
            if not mask.any(): continue
            edges = np.diff(mask.astype(int))
            starts = np.where(edges == 1)[0] + 1
            ends = np.where(edges == -1)[0] + 1
            if mask[0]: starts = np.r_[0, starts]
            if mask[-1]: ends = np.r_[ends, len(mask)]
            best = int(np.argmax(ends - starts))
            step_start = t[starts[best]]; step_end = t[ends[best] - 1]
            mm = np.zeros_like(mask); mm[starts[best]:ends[best]] = True
            step_amp = float(np.round(i_pa[mm].mean()))
            if step_amp >= 0: continue
            pre, post = 0.15, 0.30
            t0 = max(0, int((step_start - pre) * rate))
            t1 = min(len(t), int((step_end + post) * rate))
            v_seg = v_mv[t0:t1]; i_seg = i_pa[t0:t1]
            req_end = int((step_end + 0.10 - t[t0]) * rate)
            zero_mask = np.abs(v_seg) < 1e-9
            fz = int(np.argmax(zero_mask)) if zero_mask.any() else -1
            if zero_mask.any() and fz < req_end: continue
            if zero_mask.any(): v_seg = v_seg[:fz]; i_seg = i_seg[:fz]
            t_seg = t[t0:t0 + len(v_seg)] - step_start
            candidates.append(dict(step_amp_pa=step_amp,
                                    t=t_seg, v=v_seg, i=i_seg))
    if not candidates: return None
    return min(candidates, key=lambda c: abs(c["step_amp_pa"] - target_pa))

all_traces = []
trace_meta = []
for cell in TRACES_CELLS:
    sw = best_hyperpol_sweep(cell["nwb_path"], cell["target_pa"])
    if sw is None:
        print(f"  WARN: no sweep for {cell['specimen_id']}")
        continue
    # Decimate to keep CSV tractable (raw is 50 kHz; 1 kHz is enough for plotting)
    decim = max(1, int(len(sw["t"]) / 1500))
    t_ms = sw["t"][::decim] * 1000
    v_mv = sw["v"][::decim]
    i_pa = sw["i"][::decim]
    df_trace = pd.DataFrame(dict(
        cell=cell["supertype"], t_ms=t_ms, v_mV=v_mv, i_pA=i_pa))
    all_traces.append(df_trace)
    # Per-cell stats for sag annotation in R
    pre_mask = sw["t"] < -0.02
    step_mask = (sw["t"] > 0.05) & (sw["t"] < 0.7)
    late_mask = (sw["t"] > 0.7) & (sw["t"] < 1.0)
    V_base = float(sw["v"][pre_mask].mean()) if pre_mask.any() else np.nan
    V_peak_idx = int(np.argmin(sw["v"][step_mask]))
    V_peak_t_ms = float(sw["t"][step_mask][V_peak_idx] * 1000)
    V_peak = float(sw["v"][step_mask][V_peak_idx])
    V_steady = float(sw["v"][late_mask].mean()) if late_mask.any() else np.nan
    trace_meta.append(dict(
        specimen_id=cell["specimen_id"], supertype=cell["supertype"],
        layer=cell["layer"], sag=cell["sag"], color=cell["color"],
        step_amp_pa=sw["step_amp_pa"],
        V_base=V_base, V_peak=V_peak, V_peak_t_ms=V_peak_t_ms,
        V_steady=V_steady,
        label=f"{cell['specimen_id']} · {cell['supertype']} · {cell['layer']} "
              f"(sag = {cell['sag']:.2f}, step = {int(sw['step_amp_pa'])} pA)"))

pd.concat(all_traces, ignore_index=True).to_csv(OUT / "panel_G_traces.csv", index=False)
pd.DataFrame(trace_meta).to_csv(OUT / "panel_G_meta.csv", index=False)
print(f"  → {sum(len(t) for t in all_traces)} samples across {len(trace_meta)} cells")

print(f"\nAll CSVs written to: {OUT}")
print(f"Files: {sorted([f.name for f in OUT.glob('*.csv')])}")

# Record which upstream table each panel CSV came from, so the R renderer stops
# rather than drawing stale numbers. Source map: r_panels_provenance.py.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from r_panels_provenance import record as _record_provenance
_record_provenance(OUT)
