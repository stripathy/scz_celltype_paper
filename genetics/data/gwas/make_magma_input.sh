#!/usr/bin/env bash
# make_magma_input.sh — build the MAGMA gene-analysis input from the Bigdeli
# et al. 2026 GWAS-VCF release.
#
#   SCZ_EUR_autosomes.bcf  ->  magma_bigdeli/bigdeli.dedup.pval
#
# run_magma_bigdeli.sh takes that .pval as a given input; this is the step that
# produces it. The BCF is the only source of the per-variant effective sample
# size (FORMAT/NE) that MAGMA is passed as ncol=N, so the analysis cannot be
# rebuilt from the LocusZoom export in this directory, which carries no N.
#
# Source: Synapse syn60527562 (Bigdeli et al. 2026 ancestry-specific and
# cross-ancestry meta-analyses), GWAS-VCF format.
#
# Three things happen here, in order:
#   1. rsID filter   — records whose ID is not an rs number are dropped
#                      (941,621 of 12,991,921), because MAGMA joins to the
#                      g1000_eur panel by rsID.
#   2. P from LP     — the VCF carries -log10(p) in FORMAT/LP.
#   3. rsID dedup    — 11,722 rsIDs appear on more than one record (multi-
#                      allelic sites and build-38 remappings). The smallest
#                      p-value per rsID is kept. Note this is the permissive
#                      choice; it affects 0.1% of variants.
#
# Output is 12,037,701 rows, sorted by rsID, with header SNP/P/N.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BCF="${1:-$HERE/SCZ_EUR_autosomes.bcf}"
OUT="${2:-$HERE/magma_bigdeli/bigdeli.dedup.pval}"

command -v bcftools >/dev/null || { echo "bcftools not found" >&2; exit 1; }
[ -f "$BCF" ] || { echo "missing $BCF — see README.md" >&2; exit 1; }

mkdir -p "$(dirname "$OUT")"

{
  printf 'SNP\tP\tN\n'
  bcftools query -f '%ID\t[%LP]\t[%NE]\n' "$BCF" \
    | awk -F'\t' '$1 ~ /^rs/ { printf "%s\t%g\t%g\n", $1, 10 ^ (-$2), $3 }' \
    | LC_ALL=C sort -k1,1 -k2,2g \
    | awk -F'\t' '!seen[$1]++'
} > "$OUT"

echo "wrote $OUT ($(($(wc -l < "$OUT") - 1)) variants)"
