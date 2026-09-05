#!/usr/bin/env python3
"""
Gene and exon structures for the HCN1 locus zoom (Figure 4c, gene track).

Extracted from export_for_R.py, whose other sections all built the retired
version of Figure 4 (501-type taxonomy, PGC3 sumstats with an hg19->hg38
liftover, the MTG expression reference, two patch-seq cells). Only this gene
track survived the 2026-08-31 rebuild, and re-running the old script would have
overwritten the live panels with retired data.

Run AFTER export_panels_abc.py, which writes the window into panel_D_meta.csv.

Input : genetics/data/gwas/ncbiRefSeq_hg38.txt.gz  (UCSC RefSeq, hg38; not
        committed -- see genetics/data/README.md for the download)
        genetics/results/figures/r_panels/panel_D_meta.csv  (the plot window)
Output: r_panels/panel_D_genes.csv, r_panels/panel_D_exons.csv
Read by: scz_sst_hcn1_story.R, build_hcn1_locus_plot()
"""
import gzip

import pandas as pd

import os as _os

# Repo-relative, so the chain runs from any clone. External data that is not in
# the repo is still resolved by absolute path or an environment override below.
GEN = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
PAPER = _os.path.dirname(GEN)


REFGENE = f"{GEN}/data/gwas/ncbiRefSeq_hg38.txt.gz"
OUT = f"{GEN}/results/figures/r_panels"

meta = pd.read_csv(f"{OUT}/panel_D_meta.csv")
win_lo, win_hi = int(meta.win_lo_hg38[0]), int(meta.win_hi_hg38[0])

# One row per gene (longest transcript wins) and one row per exon segment,
# each segment classified 5'UTR / CDS / 3'UTR so the track can draw thin UTRs
# and thick coding blocks.
gene_rows, exon_rows, gene_seen = [], [], {}
with gzip.open(REFGENE, "rt") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 13 or parts[2] != "chr5":
            continue
        tx_start, tx_end = int(parts[4]), int(parts[5])
        if tx_end < win_lo or tx_start > win_hi:
            continue
        symbol = parts[12]
        tx_len = tx_end - tx_start
        if symbol in gene_seen and gene_seen[symbol] >= tx_len:
            continue
        gene_seen[symbol] = tx_len
        cds_start, cds_end = int(parts[6]), int(parts[7])
        ex_s = [int(x) for x in parts[9].rstrip(",").split(",") if x]
        ex_e = [int(x) for x in parts[10].rstrip(",").split(",") if x]
        gene_rows.append(dict(
            symbol=symbol, strand=parts[3],
            start_mb=tx_start / 1e6, stop_mb=tx_end / 1e6,
            cds_start_mb=cds_start / 1e6, cds_end_mb=cds_end / 1e6,
            is_hcn1=(symbol == "HCN1"), tx_len=tx_len))

        def seg(a, b, seg_type):
            if b > a:
                exon_rows.append(dict(symbol=symbol, strand=parts[3],
                                      seg_type=seg_type,
                                      is_hcn1=(symbol == "HCN1"),
                                      start_mb=a / 1e6, stop_mb=b / 1e6))

        for es, ee in zip(ex_s, ex_e):
            if ee < win_lo or es > win_hi:
                continue
            if cds_start == cds_end:                     # non-coding transcript
                seg(max(es, win_lo), min(ee, win_hi), "UTR")
                continue
            if es < cds_start:                           # 5' UTR
                seg(max(es, win_lo), min(ee, cds_start), "UTR")
            a, b = max(es, cds_start), min(ee, cds_end)  # CDS
            seg(max(a, win_lo), min(b, win_hi), "CDS")
            if ee > cds_end:                             # 3' UTR
                seg(max(es, cds_end), min(ee, win_hi), "UTR")

# Lane assignment so gene labels don't collide: sort by start, place each gene
# in the earliest lane whose previous label has ended.
genes_df = pd.DataFrame(gene_rows).sort_values("start_mb").reset_index(drop=True)
genes_df["lane"] = 0
span_mb = (win_hi - win_lo) / 1e6
lanes_end = []
for i, g in genes_df.iterrows():
    label_extent = g.stop_mb + len(g.symbol) * span_mb * 0.008 + span_mb * 0.015
    for li, le in enumerate(lanes_end):
        if g.start_mb > le:
            genes_df.at[i, "lane"] = li
            lanes_end[li] = label_extent
            break
    else:
        genes_df.at[i, "lane"] = len(lanes_end)
        lanes_end.append(label_extent)
genes_df["n_lanes"] = max(genes_df.lane) + 1 if len(genes_df) else 1
genes_df.to_csv(f"{OUT}/panel_D_genes.csv", index=False)

exon_df = pd.DataFrame(exon_rows)
exon_df["lane"] = exon_df.symbol.map(dict(zip(genes_df.symbol, genes_df.lane)))
exon_df.to_csv(f"{OUT}/panel_D_exons.csv", index=False)
print(f"panel c gene track: {len(genes_df)} genes, {len(exon_df)} exon segments "
      f"in chr5:{win_lo:,}-{win_hi:,}")
