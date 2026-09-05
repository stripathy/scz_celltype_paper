#!/usr/bin/env python3
"""
Choose the GSE158516 cell-calling threshold from the barcode-rank curves.

Compares three candidate rules against the authors' own per-sample nuclei counts
(Supplementary Table 2), and draws the knee plots so the choice is made by looking
at the data rather than by assumption:

  ordmag   CellRanger v2: 99th pct of the top-N barcodes / 10. Reference point.
  topN     top-N barcodes by UMI, N = the authors' reported nuclei count.
  inflect  steepest-drop (inflection) point of the log-log rank/UMI curve, i.e. a
           rule that uses no information from the authors at all.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TOT = "/Users/shreejoy/Github/shared_data/GSE158516/barcode_totals"
META = "/Users/shreejoy/Github/scz_celltype_paper/reserve/gse158516/data/sample_metadata.csv"
OUTD = "/Users/shreejoy/Github/scz_celltype_paper/reserve/gse158516/output"
MIN_UMI_FLOOR, MIN_GENES = 500, 250

meta = pd.read_csv(META).set_index("Deidentified ID")
sids = sorted(f[:-4] for f in os.listdir(TOT) if f.endswith(".npz"))


def inflection(umi_sorted, lo=500, hi=50_000, w=250):
    """Steepest point of log10(UMI) vs log10(rank), searched over plausible ranks."""
    n = umi_sorted.size
    lo = max(lo, w + 1)
    hi = min(hi, n - w - 1)
    if hi <= lo:
        return n, umi_sorted[-1]
    y = np.log10(np.maximum(umi_sorted, 1))
    x = np.log10(np.arange(1, n + 1))
    ks = np.arange(lo, hi)                       # centred +/-w finite difference
    slope = (y[ks + w] - y[ks - w]) / (x[ks + w] - x[ks - w])
    k = int(ks[np.argmin(slope)])
    return k, umi_sorted[k]


rows = []
curves = {}
for sid in sids:
    z = np.load(f"{TOT}/{sid}.npz", allow_pickle=True)
    umi, ngene = z["umi"], z["ngene"]
    o = np.argsort(umi)[::-1]
    umi_s, ngene_s = umi[o], ngene[o]
    curves[sid] = umi_s

    reported = int(meta.loc[sid, "Number of Nuclei"]) if sid in meta.index else np.nan
    exp = int(reported) if np.isfinite(reported) else 12_000

    top = umi_s[:exp]
    thr_ordmag = max(np.percentile(top[top > 0], 99) / 10.0, MIN_UMI_FLOOR)
    n_ordmag = int(((umi_s >= thr_ordmag) & (ngene_s >= MIN_GENES)).sum())

    k_inf, thr_inf = inflection(umi_s)
    thr_inf = max(thr_inf, MIN_UMI_FLOOR)
    n_inf = int(((umi_s >= thr_inf) & (ngene_s >= MIN_GENES)).sum())

    thr_topN = umi_s[exp - 1]
    rows.append(dict(sample=sid, in_paper=sid in meta.index, reported=reported,
                     thr_ordmag=thr_ordmag, n_ordmag=n_ordmag,
                     thr_inflect=thr_inf, n_inflect=n_inf,
                     thr_topN=thr_topN, n_topN=exp,
                     med_umi_topN=float(np.median(umi_s[:exp])),
                     med_umi_inflect=float(np.median(umi_s[:n_inf])) if n_inf else np.nan))

df = pd.DataFrame(rows)
df.to_csv(f"{OUTD}/cell_calling_comparison.csv", index=False)
pd.set_option("display.width", 200)
print(df.round(0).to_string(index=False))

sub = df[df.in_paper]
for rule in ["ordmag", "inflect"]:
    r = np.corrcoef(sub["reported"], sub[f"n_{rule}"])[0, 1]
    ratio = (sub[f"n_{rule}"] / sub["reported"])
    print(f"\n{rule:8s} vs reported: r={r:.3f}  n/reported median={ratio.median():.2f} "
          f"(range {ratio.min():.2f}-{ratio.max():.2f})")
# does our recovered median UMI match theirs? (their Supp Table 2 column)
m = sub.merge(meta.reset_index()[["Deidentified ID", "Median UMI per Nuclei"]],
              left_on="sample", right_on="Deidentified ID")
for rule in ["topN", "inflect"]:
    r = np.corrcoef(m["Median UMI per Nuclei"], m[f"med_umi_{rule}"])[0, 1]
    print(f"median-UMI check, {rule:8s}: r={r:.3f}  ours/theirs median="
          f"{(m['med_umi_'+rule]/m['Median UMI per Nuclei']).median():.2f}")

# ---- knee plots -------------------------------------------------------------
fig, axes = plt.subplots(4, 8, figsize=(34, 17), sharex=True, sharey=True)
for ax, sid in zip(axes.ravel(), sids):
    u = curves[sid]
    r = df[df["sample"] == sid].iloc[0]
    ax.loglog(np.arange(1, u.size + 1), np.maximum(u, 1), lw=1.6, color="0.25")
    ax.axvline(r["n_inflect"], color="tab:red", lw=1.6, label="inflection")
    if np.isfinite(r["reported"]):
        ax.axvline(r["reported"], color="tab:blue", lw=1.6, ls="--", label="authors' N")
    ax.axvline(r["n_ordmag"], color="tab:green", lw=1.2, ls=":", label="ordmag")
    ax.set_title(f"{sid}{'' if r['in_paper'] else '  (dropped)'}", fontsize=15)
    ax.tick_params(labelsize=12)
axes[0, 0].legend(fontsize=12, loc="lower left")
fig.supxlabel("barcode rank", fontsize=20)
fig.supylabel("UMI count", fontsize=20)
fig.suptitle("GSE158516 barcode-rank curves (CellRanger raw matrices) with candidate cell-call thresholds",
             fontsize=22)
fig.tight_layout()
fig.savefig(f"{OUTD}/knee_plots.png", dpi=110)
print(f"\nwrote {OUTD}/knee_plots.png and cell_calling_comparison.csv")
