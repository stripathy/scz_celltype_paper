suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr);library(stringr);library(purrr)})
OUT <- "transcriptomic/results/sst_strata_gsea"; POUT <- file.path(OUT,"pseudobulk")
fmt <- function(x,d=2) formatC(x, digits=d, format="f")
sci <- function(p) formatC(p, format="e", digits=1)

cat("################ 1. STRATA DEFINITION ################\n")
sd_ <- read_csv(file.path(OUT,"strata_definition.csv"), show_col_types=FALSE)
print(sd_ |> group_by(stratum) |> summarise(n=n(), members=paste(sort(CellType), collapse=", "),
      mean_depth=round(mean(depth_xenium),2), .groups="drop"))
ct <- cor.test(sd_$estimate, sd_$depth_xenium, method="spearman")
cat(sprintf("Spearman rho(crumblr est, median depth) = %.2f, p = %.3g, n = %d\n", ct$estimate, ct$p.value, nrow(sd_)))

cat("\n################ 2. GSEA BURDEN (panel b) ################\n")
g <- read_csv(file.path(OUT,"gsea_all_signatures.csv"), show_col_types=FALSE)
print(g |> mutate(dir=ifelse(NES<0,"down","up")) |> group_by(signature,dir) |>
  summarise(n05=sum(padj<0.05), n10=sum(padj<0.10), .groups="drop") |>
  pivot_wider(names_from=dir, values_from=c(n05,n10)))
gt <- g |> distinct(signature, pathway) |> count(signature); cat("gene sets tested per signature:\n"); print(gt)
cat("depleted FDR<0.05: total & pct down:\n")
dep05 <- g |> filter(signature=="depleted", padj<0.05); cat(sprintf("  n=%d, down=%d (%.0f%%)\n", nrow(dep05), sum(dep05$NES<0), 100*mean(dep05$NES<0)))
int05 <- g |> filter(signature=="intermediate", padj<0.05); cat(sprintf("  intermediate n=%d\n", nrow(int05)))
all05 <- g |> filter(signature=="all_sst", padj<0.05); cat(sprintf("  all_sst n=%d\n", nrow(all05)))

cat("\n################ 3. NAMED GENE SETS (panel d) ################\n")
SETS <- c("GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING","GOCC_SYNAPTIC_MEMBRANE",
  "GOMF_UBIQUITIN_CONJUGATING_ENZYME_BINDING","GOBP_PROTEIN_K63_LINKED_DEUBIQUITINATION",
  "GOCC_RIBOSOMAL_SUBUNIT","GOBP_CYTOPLASMIC_TRANSLATION",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT","GOBP_OXIDATIVE_PHOSPHORYLATION","HALLMARK_OXIDATIVE_PHOSPHORYLATION")
print(g |> filter(pathway %in% SETS) |> mutate(cell=sprintf("%.2f (%.3g)", NES, padj)) |>
  select(pathway, signature, cell) |> pivot_wider(names_from=signature, values_from=cell) |>
  arrange(match(pathway, SETS)), width=200)

cat("\n################ 4. UPREGULATED SETS, depleted, padj<0.05 ################\n")
up <- g |> filter(signature=="depleted", padj<0.05, NES>0) |> arrange(padj)
cat(sprintf("n=%d\n", nrow(up))); print(up |> select(pathway, NES, padj) |> mutate(NES=round(NES,2), padj=signif(padj,2)), n=15)

cat("\n################ 5. GENE-LEVEL META (stratum_meta_de) ################\n")
sm <- read_csv(file.path(POUT,"stratum_meta_de.csv"), show_col_types=FALSE)
cat("genes tested per stratum:\n"); print(sm |> count(stratum))
cat("genes at padj<0.05 / <0.10 per stratum:\n")
print(sm |> group_by(stratum) |> summarise(n05=sum(padj<0.05), n10=sum(padj<0.10)))
cat("FDR<0.10 genes by stratum:\n")
print(sm |> filter(padj<0.10) |> group_by(stratum) |> summarise(genes=paste(gene[order(padj)], collapse=", ")))
GENES <- c("SST","VGF","NMU","RASGRF1","RPL36","EIF3G","NDUFS8","COX5B","STAMBPL1","USP8","CALB1")
cat("named genes: z (padj) per stratum:\n")
print(sm |> filter(gene %in% GENES) |> mutate(cell=sprintf("%.1f (%.2g)", zval, padj)) |>
  select(gene, stratum, cell) |> pivot_wider(names_from=stratum, values_from=cell) |>
  arrange(match(gene, GENES)), width=200)

cat("\n################ 6. INTERACTION MODEL (genes) ################\n")
im <- read_csv(file.path(POUT,"interaction_meta.csv"), show_col_types=FALSE)
print(im |> group_by(coef) |> summarise(n05=sum(padj<0.05), n10=sum(padj<0.10), n_tested=n()))
cat("dep-vs-non FDR<0.05 genes:\n")
print(im |> filter(coef=="dxSCZ:stratumdepleted", padj<0.05) |> arrange(padj) |> select(gene, zval, padj))
cat("named genes, dep-vs-non interaction:\n")
print(im |> filter(coef=="dxSCZ:stratumdepleted", gene %in% c("SST","VGF","COX5B","NDUFS8","RPL36","EIF3G")) |>
  transmute(gene, z=round(zval,2), padj=signif(padj,2)))

cat("\n################ 7. INTERACTION GSEA ################\n")
ig <- read_csv(file.path(POUT,"interaction_gsea.csv"), show_col_types=FALSE)
ISETS <- c("REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION","GOBP_CYTOPLASMIC_TRANSLATION","GOCC_RIBOSOMAL_SUBUNIT",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION","GOBP_OXIDATIVE_PHOSPHORYLATION","REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING","GOCC_SYNAPTIC_MEMBRANE","GOBP_PROTEIN_K63_LINKED_DEUBIQUITINATION")
print(ig |> filter(pathway %in% ISETS, coef %in% c("dxSCZ:stratumdepleted","dxSCZ:stratumintermediate")) |>
  mutate(cell=sprintf("%.2f (%.2g)", NES, padj)) |> select(pathway, coef, cell) |>
  pivot_wider(names_from=coef, values_from=cell) |> arrange(match(pathway, ISETS)), width=220)
cat("top interaction sets (dep-vs-non) by padj:\n")
print(ig |> filter(coef=="dxSCZ:stratumdepleted") |> arrange(padj) |> head(8) |> transmute(pathway, NES=round(NES,2), padj=signif(padj,2)))

cat("\n################ 8. SIGN BIAS within modules (depleted stratum) ################\n")
le_union <- function(ps) g |> filter(pathway %in% ps, signature=="depleted") |> pull(leadingEdge) |>
  strsplit("|",fixed=TRUE) |> unlist() |> unique()
MODS <- list(translation=le_union(c("REACTOME_TRANSLATION","GOCC_RIBOSOMAL_SUBUNIT","GOBP_CYTOPLASMIC_TRANSLATION","REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION")),
  oxphos=le_union(c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","GOBP_OXIDATIVE_PHOSPHORYLATION","REACTOME_RESPIRATORY_ELECTRON_TRANSPORT")),
  synaptic=le_union(c("GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING","GOCC_SYNAPTIC_MEMBRANE","GOBP_NEUROTRANSMITTER_SECRETION","REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES")))
cat(sprintf("module sizes: %s\n", paste(names(MODS), lengths(MODS), collapse=", ")))
for (m in names(MODS)) for (st in c("depleted","non_depleted")) {
  z <- sm |> filter(stratum==st, gene %in% MODS[[m]]) |> pull(zval)
  bt <- binom.test(sum(z<0), length(z))
  cat(sprintf("  %-12s %-13s n=%3d  %%neg=%2.0f%%  binom p=%s\n", m, st, length(z), 100*mean(z<0), sci(bt$p.value)))
}

cat("\n################ 9. DONOR / NUCLEUS COUNTS ################\n")
st_of <- setNames(sd_$stratum, sd_$CellType)
gr <- map_dfr(Sys.glob("transcriptomic/data/stratum_pseudobulks_export/*_groups.csv"),
  ~read_csv(.x, show_col_types=FALSE, col_types=cols(donor=col_character(), .default=col_guess())))
sst <- gr |> filter(supertype %in% names(st_of)) |> mutate(stratum=st_of[supertype])
ds <- sst |> group_by(cohort, donor, diagnosis, stratum) |> summarise(n=sum(n_cells), .groups="drop")
cat(sprintf("cohorts: %d; donors with any Sst nuclei: %d (%s)\n", n_distinct(ds$cohort),
  n_distinct(paste(ds$cohort, ds$donor)),
  paste(capture.output(print(ds |> distinct(cohort,donor,diagnosis) |> count(diagnosis))), collapse=" ")))
print(ds |> group_by(stratum) |> summarise(donors_ge10=sum(n>=10), total_nuclei=sum(n), med_per_donor=median(n)))
dsd <- ds |> filter(n>=10) |> distinct(cohort,donor,diagnosis) ; cat("donors >=10 in >=1 stratum, by dx:\n"); print(dsd |> count(diagnosis))

cat("\n################ 10. XENIUM ################\n")
xc <- read_csv(file.path(OUT,"xenium_stratum_concordance.csv"), show_col_types=FALSE)
panel_genes <- unique(xc$gene); cat(sprintf("panel genes with stratum estimates: %d\n", length(panel_genes)))
print(xc |> filter(gene %in% c("SST","VGF")) |> transmute(gene, stratum, xen_z=round(xen_z,1)))
for (m in names(MODS)) cat(sprintf("  %-12s LE genes on panel: %d of %d\n", m, sum(MODS[[m]] %in% panel_genes), length(MODS[[m]])))
