#!/usr/bin/env bash
# run_magma_bigdeli.sh — gene-level MAGMA on the Bigdeli et al. 2026 SCZ EUR
# meta-analysis (MVP + PGC3 + Release7 + FinnGen), reproducing the EXACT settings
# of the original PGC3 run so that only the GWAS data changes:
#   - MAGMA v1.10            (same binary version)
#   - window 35,10           (35 kb upstream / 10 kb downstream)
#   - g1000_eur LD panel      (identical: 22,665,064 SNPs, 503 EUR individuals)
#   - NCBI37.3.gene.loc       (hg19, full — MHC dropped later in the Python join)
#   - per-variant effective N (Bigdeli NE column; mirrors PGC3's per-SNP NEFF)
#
# Build handling: Bigdeli sumstats are GRCh38, but the panel + gene-loc are hg19.
# We therefore annotate using the g1000_eur panel's OWN hg19 SNP positions and
# join the Bigdeli p-values by rsID (use=SNP,P). No liftOver needed; gene
# definitions are byte-identical to the PGC3 run.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT=/Users/shreejoy/Github/scz_celltype_paper/genetics
MAGMA="$ROOT/data/magma/magma_mac/magma"
PANEL="$ROOT/data/magma/g1000_eur/g1000_eur"
GENELOC="$HOME/Github/scz_cell_type_enrichment/linking_cell_types_to_brain_phenotypes/Data/NCBI37.3.gene.loc"
PVAL="$HERE/bigdeli.dedup.pval"
OUT="$HERE/bigdeli"

# 1) snp-loc from the panel .bim (hg19):  SNP  CHR  BP
awk '{print $2"\t"$1"\t"$4}' "${PANEL}.bim" > "$HERE/g1000_eur.snploc"

# 2) annotate SNPs to genes (same window as PGC3 step1)
"$MAGMA" --annotate window=35,10 \
  --snp-loc "$HERE/g1000_eur.snploc" \
  --gene-loc "$GENELOC" \
  --out "${OUT}.step1"

# 3) gene analysis (same as PGC3 step2, but ncol=N = per-variant effective N)
"$MAGMA" --bfile "$PANEL" \
  --pval "$PVAL" use=SNP,P ncol=N \
  --gene-annot "${OUT}.step1.genes.annot" \
  --out "${OUT}.step2"

echo "Done -> ${OUT}.step2.genes.out"
