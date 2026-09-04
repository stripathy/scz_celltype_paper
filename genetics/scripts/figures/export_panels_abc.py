"""
Rebuild Figure 4 panels a, b and c on Bigdeli 2026 (European-ancestry).

  a  genetics-vs-depletion  <- Bigdeli EUR gene-property enrichment across the
                               125 SEA-AD DLPFC supertypes alone (the combined
                               SEA-AD + Siletti taxonomy was retired on
                               L. Duncan's advice; build_spec_seaad_only.py)
  b  Sst_2 gene drivers     <- same enrichment + Bigdeli EUR gene-level P
                               (Sst_2 is both the most SCZ-enriched Sst
                               supertype and the most significantly depleted,
                               so panel b names the type panel a puts at the
                               corner of the convergence; HCN1 is a top-decile
                               driver for Sst_2, Sst_3 and Sst_20 alike, so the
                               HCN1 argument does not rest on the choice)
  c  HCN1 locus zoom        <- Bigdeli EUR sumstats (already hg38, so the
                               hg19->hg38 liftover the PGC3 build needed is
                               dropped) + the SuSiE-R EUR credible set from
                               Bigdeli Supplementary Table 13, LOCUS319:
                               4 variants, cumulative PIP 0.985,
                               lead rs10035564 PIP 0.556.

The credible set is taken from the ancestry-matched SuSiE-R EUR analysis so it
matches the EUR sumstats used for enrichment. Their tables also report this SNP
under five other method/ancestry combinations spanning PIP 0.001-0.841; that
spread belongs in the supplement, not in a panel.
"""
import gzip
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

W = "/Users/shreejoy/Github/scz_celltype_paper/genetics/results/intermediates"
GEN = "/Users/shreejoy/Github/scz_celltype_paper/genetics"
SRC = "/Users/shreejoy/Github/scz_cell_type_enrichment"
GLOC = f"{SRC}/linking_cell_types_to_brain_phenotypes/Data/NCBI37.3.gene.loc.extendedMHCexcluded"
SUMSTATS = f"{GEN}/data/gwas/bigdeli_eur_scz_sum_stats.gz"
OUT = "/Users/shreejoy/Github/scz_celltype_paper/genetics/results/figures/r_panels"
TARGET = "Sst_2"

names = pd.read_csv(f"{W}/namemap_a9only.csv").set_index("safe_name").cell_type
gsa = pd.read_csv(f"{W}/T_a9only_bigdeli.gsa.out", sep=r"\s+", comment="#")
gsa["supertype"] = gsa.VARIABLE.map(names)
enrich = gsa.dropna(subset=["supertype"])[["supertype", "BETA", "P"]].rename(
    columns={"BETA": "beta", "P": "p_value"})
enrich["p_fdr"] = stats.false_discovery_control(enrich.p_value.values)

# ---------------------------------------------------------------- panel a
pa = pd.read_csv(f"{GEN}/results/figures/r_panels/panel_A_enrichment.csv")
SEAAD_COLORS = dict(zip(pa.supertype, pa.color))
comp = pd.read_csv(f"{GEN}/results/tables/gwas_vs_casecontrol_composition.csv")
sst = comp[comp.supertype.str.startswith("Sst")
           & ~comp.supertype.str.startswith("Sst Chodl")].copy()
sst["beta"] = sst.supertype.map(dict(zip(enrich.supertype, enrich.beta)))
sst["scz_p"] = sst.supertype.map(dict(zip(enrich.supertype, enrich.p_value)))
sst["scz_neg_log10_p"] = -np.log10(sst.scz_p.clip(lower=1e-300))
sst["depleted_fdr20"] = (sst.comp_padj < 0.20) & (sst.comp_beta < 0)
sst["exemplar"] = sst.supertype.isin(["Sst_25", "Sst_5"])
sst["pt_size"] = 1.0 + 3.0 * np.clip(sst.comp_logp.values / 5, 0, 1)
sst["color"] = sst.supertype.map(SEAAD_COLORS).fillna("#888888")
cols = ["supertype", "beta", "scz_p", "scz_neg_log10_p", "comp_beta", "comp_padj",
        "comp_logp", "depleted_fdr20", "exemplar", "pt_size", "color"]
sst[cols].to_csv(f"{OUT}/panel_B_genetics_vs_depletion.csv", index=False)
rho, p = stats.spearmanr(sst.scz_neg_log10_p, -sst.comp_beta)
print(f"panel a: n={len(sst)}  Spearman rho={rho:.3f}, p={p:.4f}")

# ---------------------------------------------------------------- panel b
gl = pd.read_csv(GLOC, sep=r"\s+", header=None,
                 names=["ENTREZ", "CHR", "START", "STOP", "STRAND", "SYMBOL"])
ent2sym = dict(zip(gl.ENTREZ.astype(int), gl.SYMBOL))
col = pd.read_csv(f"{W}/namemap_a9only.csv").query("cell_type == @TARGET").safe_name.iloc[0]
spec = pd.read_csv(f"{W}/spec_dlpfc_a9only.txt", sep="\t", usecols=["GENE", col])
spec["symbol"] = spec.GENE.astype(int).map(ent2sym)
magma = pd.read_csv(f"{GEN}/data/gwas/magma_bigdeli/bigdeli.step2.genes.out", sep=r"\s+")
magma["symbol"] = magma.GENE.astype(int).map(ent2sym)
m = (spec.dropna(subset=["symbol"])[["symbol", col]].rename(columns={col: "specificity"})
     .merge(magma[["symbol", "P"]].dropna(), on="symbol", how="inner"))
m["neg_log10_p"] = -np.log10(m.P.clip(lower=1e-300))
m["p_fdr"] = multipletests(m.P, method="fdr_bh")[1]
q90 = m.specificity.quantile(0.90)
m["is_driver"] = (m.specificity >= q90) & (m.p_fdr < 0.05)
m["is_hcn1"] = m.symbol == "HCN1"
m.to_csv(f"{OUT}/panel_C_gene_drivers.csv", index=False)

drv = m[m.is_driver].copy()
drv["rank_spec"] = drv.specificity.rank(ascending=False, method="min")
drv["rank_p"] = drv.neg_log10_p.rank(ascending=False, method="min")
drv["combined_rank"] = drv.rank_spec + drv.rank_p
drv = drv.join(gl.set_index("SYMBOL")[["CHR", "START", "STOP"]], on="symbol", how="left")
drv["mid_bp"] = (drv.START + drv.STOP) / 2
acc, keep = [], []
for _, r in drv.sort_values("combined_rank").reset_index(drop=True).iterrows():
    if pd.isna(r.CHR) or pd.isna(r.mid_bp):
        continue
    if any(c == r.CHR and abs(mb - r.mid_bp) < 1_000_000 for c, mb in acc):
        continue
    acc.append((r.CHR, r.mid_bp)); keep.append(r.symbol)
    if len(keep) >= 10:
        break
top = drv[drv.symbol.isin(keep)].copy()
top["is_hcn1"] = top.symbol == "HCN1"
if "HCN1" not in top.symbol.values and "HCN1" in drv.symbol.values:
    h = drv[drv.symbol == "HCN1"].copy(); h["is_hcn1"] = True
    top = pd.concat([top, h], ignore_index=True)
top.to_csv(f"{OUT}/panel_C_top_labels.csv", index=False)
fp = m[m.p_fdr < 0.05]
t = enrich[enrich.supertype == TARGET].iloc[0]
# log-x headroom above the most specific driver, as in the combined-taxonomy
# figure (0.05 over a 0.022 max); specificity shares are ~4x larger now that
# the denominator is 125 types instead of 501, so the cap scales with the data
xlim_hi = float(m.loc[m.is_driver, "specificity"].max()) * 2.25
pd.DataFrame([dict(spec_q90=q90, xlim_lo=q90 * 0.4, xlim_hi=xlim_hi, ylim_lo=1.3,
                   ylim_hi=20.0,
                   fdr_neg_log10p_cutoff=-np.log10(fp.P.max()) if len(fp) else np.nan,
                   target_type=TARGET, target_beta=t.beta, target_p_fdr=t.p_fdr)]
             ).to_csv(f"{OUT}/panel_C_meta.csv", index=False)
hc = m[m.symbol == "HCN1"].iloc[0]
print(f"panel b: {int(m.is_driver.sum())} drivers | HCN1 spec {hc.specificity:.5f} "
      f"(pct {100*(m.specificity < hc.specificity).mean():.1f}), P {hc.P:.3g}, driver={bool(hc.is_driver)}")
print("  labels:", ", ".join(keep))

# ---------------------------------------------------------------- panel c
CS = pd.DataFrame({
    "rsid": ["rs10035564", "rs4492120", "rs6859397", "rs16902086"],
    "pos_hg38": [45252398, 45187702, 45135948, 45285650],
    "pip": [0.556416, 0.314875, 0.069742, 0.043675]})
old_meta = pd.read_csv(f"{GEN}/results/figures/r_panels/panel_D_meta.csv")
WIN_LO, WIN_HI = int(old_meta.win_lo_hg38[0]), int(old_meta.win_hi_hg38[0])

rows = []
with gzip.open(SUMSTATS, "rt") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[0] != "5":
            continue
        pos = int(f[1])
        if pos < WIN_LO or pos > WIN_HI:
            continue
        rows.append((f[2], pos, float(f[5])))       # rsid, pos, -log10 p
snps = pd.DataFrame(rows, columns=["rsid", "pos_hg38", "neg_log10_p"])
snps["p"] = 10 ** (-snps.neg_log10_p)
snps["pos_mb"] = snps.pos_hg38 / 1e6
snps["is_gwas_sig"] = snps.p < 5e-8
snps["in_credible_set"] = snps.rsid.isin(CS.rsid)
snps["pip"] = snps.rsid.map(dict(zip(CS.rsid, CS.pip))).fillna(0.0)
snps[["rsid", "pos_hg38", "pos_mb", "p", "neg_log10_p",
      "is_gwas_sig", "in_credible_set", "pip"]].to_csv(f"{OUT}/panel_D_snps.csv", index=False)

cs = CS.merge(snps[["rsid", "p", "neg_log10_p"]], on="rsid", how="left")
cs["locus"] = "rs10035564"; cs["chr"] = 5; cs["pos"] = cs.pos_hg38
cs["pos_mb"] = cs.pos_hg38 / 1e6
cs["gwas_p"] = cs.p; cs["gwas_or"] = np.nan; cs["gene"] = np.nan
cs["is_lead"] = cs.rsid == "rs10035564"
hs, he = float(old_meta.hcn1_start_mb[0]) * 1e6, float(old_meta.hcn1_stop_mb[0]) * 1e6
cs["hcn1_location"] = np.where(cs.pos_hg38.between(hs, he), "intron",
                               np.where(cs.pos_hg38 < hs, "5'_upstream", "3'_downstream"))
cs[["locus", "rsid", "chr", "pos", "pip", "gene", "gwas_p", "gwas_or", "pos_hg38",
    "pos_mb", "neg_log10_p", "is_lead", "hcn1_location"]].to_csv(
    f"{OUT}/panel_D_credible_set.csv", index=False)

lead = cs[cs.is_lead].iloc[0]
second = cs[~cs.is_lead].nlargest(1, "pip").iloc[0]
meta = old_meta.copy()
for k, v in dict(lead_rsid=lead.rsid, lead_pos_mb=lead.pos_mb, lead_p=lead.gwas_p,
                 lead_neglog10p=lead.neg_log10_p, lead_pip=lead.pip,
                 second_rsid=second.rsid, second_pos_mb=second.pos_mb,
                 second_pip=second.pip, second_loc=second.hcn1_location,
                 second_neglog10p=second.neg_log10_p,
                 n_credible=len(cs), pip_sum=cs.pip.sum()).items():
    meta[k] = v
meta.to_csv(f"{OUT}/panel_D_meta.csv", index=False)
print(f"panel c: {len(snps):,} SNPs in window | credible set n={len(cs)}, "
      f"cumulative PIP={cs.pip.sum():.3f}, lead {lead.rsid} PIP={lead.pip:.3f}, "
      f"-log10P={lead.neg_log10_p:.2f}")
