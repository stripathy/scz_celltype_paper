#!/usr/bin/env Rscript
# Step 9 | Assert every number quoted in the Figure 5 text against its source.
#
# House rule: no statistic reaches a manuscript from memory or from reading a
# plot. This script re-derives each cited value and fails loudly on any mismatch,
# so a renumber cannot silently invalidate the draft. Run it after any pipeline
# change and before pasting text into the Doc.
#
# Reference: manuscript/drafting_history/figure5_section_DRAFT.md (VERSION 0 is the working base).
source("transcriptomic/scripts/fig5/_common.R")

fails <- 0L
chk <- function(label, cond, got = NULL) {
  ok <- isTRUE(cond)
  if (!ok) {
    fails <<- fails + 1L
    cat(sprintf("  FAIL  %s%s\n", label,
                if (!is.null(got)) sprintf("  [got: %s]", paste(got, collapse = ", ")) else ""))
  }
  invisible(ok)
}

strata <- load_strata()
g   <- read_csv(file.path(P$pb, "gsea_all_signatures.csv"),      show_col_types = FALSE)
sm  <- read_csv(file.path(P$pb, "stratum_meta_de.csv"),          show_col_types = FALSE)
im  <- read_csv(file.path(P$pb, "interaction_meta.csv"),         show_col_types = FALSE) |>
  filter(coef == "dxSCZ:stratumdepleted")
ig  <- read_csv(file.path(P$pb, "interaction_gsea.csv"),         show_col_types = FALSE)
ms  <- read_csv(file.path(P$pb, "module_score_tests.csv"),       show_col_types = FALSE)
xc  <- read_csv(file.path(P$out, "xenium_stratum_concordance.csv"), show_col_types = FALSE)

# ---- strata definition -------------------------------------------------------
ct <- cor.test(strata$estimate, strata$depth_xenium, method = "spearman")
chk("rho(estimate, depth) = 0.74", round(ct$estimate, 2) == 0.74, round(ct$estimate, 2))
chk("P = 0.0016", round(ct$p.value, 4) == 0.0016, round(ct$p.value, 4))
dm <- strata |> group_by(stratum) |>
  summarise(n = n(), d = round(mean(depth_xenium), 2)) |> arrange(match(stratum, STRATA_LEVELS))
chk("stratum n = 5, 6, 5", all(dm$n == c(5, 6, 5)), dm$n)
chk("mean depths 0.31, 0.49, 0.72", all(dm$d == c(0.31, 0.49, 0.72)), dm$d)
chk("depleted members", setequal(strata$CellType[strata$stratum == "depleted"],
                                c("Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25")))

# ---- gene-set burden (headline numbers are the FDR<0.10 tier) ----------------
burden <- g |> group_by(signature) |>
  summarise(n10 = sum(padj < 0.10), pct_dn = round(100 * mean(NES[padj < 0.10] < 0)),
            n05 = sum(padj < 0.05), .groups = "drop")
bn <- function(s, col) burden[[col]][burden$signature == s]
chk("depleted 130 sets at FDR<0.10", bn("depleted", "n10") == 130, bn("depleted", "n10"))
chk("depleted 88% down", bn("depleted", "pct_dn") == 88, bn("depleted", "pct_dn"))
chk("intermediate 64", bn("intermediate", "n10") == 64, bn("intermediate", "n10"))
chk("non_depleted 0", bn("non_depleted", "n10") == 0, bn("non_depleted", "n10"))
chk("Sst_subclass 55", bn("Sst_subclass", "n10") == 55, bn("Sst_subclass", "n10"))
chk("depleted up-sets 16 of 130",
    sum(g$padj < 0.10 & g$NES > 0 & g$signature == "depleted") == 16)

# ---- gene-level counts -------------------------------------------------------
gc <- sm |> group_by(stratum) |> summarise(n10 = sum(padj < 0.10)) |>
  arrange(match(stratum, c(STRATA_LEVELS, "all_sst")))
chk("genes FDR<0.10 = 19, 99, 35, 342", all(gc$n10 == c(19, 99, 35, 342)), gc$n10)

# ---- named gene sets ---------------------------------------------------------
nv <- function(pw, sg, field) { r <- g[g$pathway == pw & g$signature == sg, ]; r[[field]] }
chk("ribosomal subunit dep NES -2.6", round(nv("GOCC_RIBOSOMAL_SUBUNIT", "depleted", "NES"), 1) == -2.6)
chk("ribosomal subunit dep padj 1.3e-12",
    signif(nv("GOCC_RIBOSOMAL_SUBUNIT", "depleted", "padj"), 2) == 1.3e-12)
chk("cytoplasmic translation dep NES -2.3",
    round(nv("GOBP_CYTOPLASMIC_TRANSLATION", "depleted", "NES"), 1) == -2.3)
chk("cytoplasmic translation dep padj 6.8e-8",
    signif(nv("GOBP_CYTOPLASMIC_TRANSLATION", "depleted", "padj"), 2) == 6.8e-8)
chk("ribosomal subunit int NES -2.2", round(nv("GOCC_RIBOSOMAL_SUBUNIT", "intermediate", "NES"), 1) == -2.2)
chk("cytoplasmic translation int NES -1.9",
    round(nv("GOBP_CYTOPLASMIC_TRANSLATION", "intermediate", "NES"), 1) == -1.9)
chk("translation sets flat in non_depleted (NES >= -0.7, padj = 1)",
    all(round(c(nv("GOCC_RIBOSOMAL_SUBUNIT", "non_depleted", "NES"),
                nv("GOBP_CYTOPLASMIC_TRANSLATION", "non_depleted", "NES")), 1) >= -0.7) &&
    all(round(c(nv("GOCC_RIBOSOMAL_SUBUNIT", "non_depleted", "padj"),
                nv("GOBP_CYTOPLASMIC_TRANSLATION", "non_depleted", "padj")), 2) == 1))
ox <- g |> filter(pathway %in% BLOCKS$pathway[BLOCKS$block == "oxphos"])
chk("oxphos dep NES -1.9..-2.0 (block = the two sets panel d plots)",
    all(round(range(ox$NES[ox$signature == "depleted"]), 1) == c(-2.0, -1.9)))
chk("oxphos dep padj <= 0.0033", max(ox$padj[ox$signature == "depleted"]) <= 0.0033,
    signif(max(ox$padj[ox$signature == "depleted"]), 3))
chk("oxphos int not significant (min padj 0.64)",
    round(min(ox$padj[ox$signature == "intermediate"]), 2) == 0.64,
    round(min(ox$padj[ox$signature == "intermediate"]), 2))
ts <- g |> filter(pathway == "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING",
                  signature %in% STRATA_LEVELS)
chk("trans-synaptic NES -1.4..-1.5 in all strata", all(round(ts$NES, 1) %in% c(-1.4, -1.5)))
chk("K63 deubiquitination dep NES +2.15, padj 0.007",
    round(nv("GOBP_PROTEIN_K63_LINKED_DEUBIQUITINATION", "depleted", "NES"), 2) == 2.15 &&
    round(nv("GOBP_PROTEIN_K63_LINKED_DEUBIQUITINATION", "depleted", "padj"), 3) == 0.007)
chk("K63 positive in every stratum",
    all(g$NES[g$pathway == "GOBP_PROTEIN_K63_LINKED_DEUBIQUITINATION"] > 0))

# ---- named genes -------------------------------------------------------------
zt <- function(gn) round(sm$zval[match(paste(gn, STRATA_LEVELS), paste(sm$gene, sm$stratum))], 1)
for (spec in list(c("SST", "-3.4,-3.6,-2.1"), c("VGF", "-3.1,-2.9,-2.9"),
                  c("NMU", "-2.8,-2.7,-3.1"), c("RPL36", "-2.5,-1.7,0.5"),
                  c("EIF3G", "-1.8,-0.4,1.1"), c("NDUFS8", "-2.9,-0.1,0.5"),
                  c("COX5B", "-2.4,-0.7,0.9")))
  chk(sprintf("%s z = %s", spec[1], spec[2]),
      all(zt(spec[1]) == as.numeric(strsplit(spec[2], ",")[[1]])), zt(spec[1]))
chk("NMU subclass FDR 0.018",
    round(sm$padj[sm$gene == "NMU" & sm$stratum == "all_sst"], 3) == 0.018)

# ---- interaction -------------------------------------------------------------
chk("7 interaction genes at FDR<0.05", sum(im$padj < 0.05) == 7, sum(im$padj < 0.05))
chk("CELF2/GRIK4/QKI among them",
    all(c("CELF2", "GRIK4", "QKI") %in% im$gene[im$padj < 0.05]))
chk("SST interaction z -0.8", round(im$zval[im$gene == "SST"], 1) == -0.8)
chk("VGF interaction z -0.6", round(im$zval[im$gene == "VGF"], 1) == -0.6)
chk("COX5B interaction z -3.8, padj 0.13",
    round(im$zval[im$gene == "COX5B"], 1) == -3.8 &&
    round(im$padj[im$gene == "COX5B"], 2) == 0.13)
iv <- function(pw, cf, field) { r <- ig[ig$pathway == pw & ig$coef == cf, ]; r[[field]] }
chk("SRP interaction NES -2.7, padj 7.9e-10",
    round(iv("REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
             "dxSCZ:stratumdepleted", "NES"), 1) == -2.7 &&
    signif(iv("REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
              "dxSCZ:stratumdepleted", "padj"), 2) == 7.9e-10)
chk("elongation interaction padj 2.2e-9",
    signif(iv("REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
              "dxSCZ:stratumdepleted", "padj"), 2) == 2.2e-9)
# OxPhos interaction is TREND-LEVEL once Hallmark is excluded. Its former
# significance (FDR = 1.6e-4) rested entirely on HALLMARK_OXIDATIVE_PHOSPHORYLATION;
# the GO and Reactome oxphos sets were never significant in the interaction.
ox_i <- ig |> filter(coef == "dxSCZ:stratumdepleted",
                     pathway %in% BLOCKS$pathway[BLOCKS$block == "oxphos"])
chk("no oxphos block set reaches interaction FDR < 0.10",
    min(ox_i$padj) > 0.10, signif(min(ox_i$padj), 3))
chk("electron transport chain interaction trend, padj 0.053",
    signif(iv("GOBP_ELECTRON_TRANSPORT_CHAIN", "dxSCZ:stratumdepleted", "padj"), 2) == 0.053)
syn_i <- ig |> filter(coef == "dxSCZ:stratumdepleted",
                      pathway %in% BLOCKS$pathway[BLOCKS$block == "synaptic"])
chk("no synaptic gene set shows an interaction (all padj > 0.10)",
    min(syn_i$padj) > 0.10, signif(min(syn_i$padj), 3))
chk("synaptic interaction min padj 0.21", round(min(syn_i$padj), 2) == 0.21, round(min(syn_i$padj), 2))

# ---- module scores -----------------------------------------------------------
mss <- ms |> filter(test == "per_stratum")
gm <- function(mod, st, f) mss[[f]][mss$module == mod & mss$stratum == st]
chk("synaptic per-stratum -0.23, -0.18, -0.10",
    all(round(c(gm("synaptic", "depleted", "estimate"), gm("synaptic", "intermediate", "estimate"),
                gm("synaptic", "non_depleted", "estimate")), 2) == c(-0.23, -0.18, -0.10)))
chk("synaptic p 7.0e-8 / 1.6e-6 / 3.9e-3",
    all(signif(c(gm("synaptic", "depleted", "pval"), gm("synaptic", "intermediate", "pval"),
                 gm("synaptic", "non_depleted", "pval")), 2) == c(7.0e-8, 1.6e-6, 3.9e-3)))
chk("oxphos depleted -0.248, 7/7 cohorts negative",
    round(gm("oxphos", "depleted", "estimate"), 3) == -0.248 &&
    gm("oxphos", "depleted", "n_neg") == 7)
chk("translation depleted -0.216, 7/7 negative",
    round(gm("translation", "depleted", "estimate"), 3) == -0.216 &&
    gm("translation", "depleted", "n_neg") == 7)
msi <- ms |> filter(test == "within_donor_interaction")
chk("paired interaction p: oxphos .032, translation .028, synaptic 1.3e-4",
    round(msi$pval[msi$module == "oxphos"], 3) == 0.032 &&
    round(msi$pval[msi$module == "translation"], 3) == 0.028 &&
    signif(msi$pval[msi$module == "synaptic"], 2) == 1.3e-4)

# ---- coverage and Xenium -----------------------------------------------------
ds <- all_donor_meta() |> group_by(stratum) |> summarise(donors = n(), nuclei = sum(n_cells), med = median(n_cells)) |>
  arrange(match(stratum, STRATA_LEVELS))
chk("donors 395, 411, 342", all(ds$donors == c(395, 411, 342)), ds$donors)
chk("nuclei 44744, 42814, 13157", all(ds$nuclei == c(44744, 42814, 13157)), ds$nuclei)
chk("median nuclei 85, 90, 30.5", all(ds$med == c(85, 90, 30.5)), ds$med)
xz <- function(gn, f) round(xc[[f]][match(paste(gn, STRATA_LEVELS), paste(xc$gene, xc$stratum))], 1)
chk("Xenium SST z -3.2, -2.2, -2.5", all(xz("SST", "xen_z") == c(-3.2, -2.2, -2.5)), xz("SST", "xen_z"))
chk("Xenium VGF z -4.6, -3.3, -2.3", all(xz("VGF", "xen_z") == c(-4.6, -3.3, -2.3)), xz("VGF", "xen_z"))
mods <- modules(g)
panel <- unique(xc$gene)
chk("translation module 130 genes, 2 on panel",
    length(mods$translation) == 130 && sum(mods$translation %in% panel) == 2,
    sum(mods$translation %in% panel))
chk("oxphos module 66 genes, 2 on panel",
    length(mods$oxphos) == 66 && sum(mods$oxphos %in% panel) == 2,
    length(mods$oxphos))
chk("synaptic module 186 genes", length(mods$synaptic) == 186, length(mods$synaptic))

cat(sprintf("\n%s  %d checks failed\n", if (fails == 0) "PASS —" else "FAILURES —", fails))
if (fails > 0) quit(status = 1)
