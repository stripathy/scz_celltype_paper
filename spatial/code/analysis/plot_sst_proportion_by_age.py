#!/usr/bin/env python3
"""
Plot total SST cell proportion (% of cortical cells) vs subject age.

One point per subject, colored by diagnosis. Includes linear regression
line and Pearson correlation for each group and overall.
"""

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from analysis.config import (
    BASE_DIR, METADATA_PATH, PRESENTATION_DIR,
    DX_COLORS, EXCLUDE_SAMPLES, SST_TYPES,
)
from modules.metadata import get_subject_info

# ── Load data ───────────────────────────────────────────────────────
prop_path = os.path.join(PRESENTATION_DIR, "xenium_composition_by_sample.csv")
df = pd.read_csv(prop_path)

# Filter to SST subtypes
sst = df[df["celltype"].isin(SST_TYPES)].copy()

# Total SST proportion per sample
sst_total = (sst.groupby(["sample_id", "diagnosis"])["proportion_pct"]
             .sum().reset_index()
             .rename(columns={"proportion_pct": "sst_pct"}))

# Merge with age
subject_info = get_subject_info(METADATA_PATH)
sst_total = sst_total.merge(subject_info[["sample_id", "age", "sex"]], on="sample_id")

# Exclude outlier samples
sst_total = sst_total[~sst_total["sample_id"].isin(EXCLUDE_SAMPLES)]

print(f"Samples: {len(sst_total)} ({(sst_total.diagnosis == 'Control').sum()} Ctrl, "
      f"{(sst_total.diagnosis == 'SCZ').sum()} SCZ)")
print(f"Age range: {sst_total.age.min():.1f} – {sst_total.age.max():.1f}")
print(f"SST proportion range: {sst_total.sst_pct.min():.2f}% – {sst_total.sst_pct.max():.2f}%")

# ── Statistics ──────────────────────────────────────────────────────
# Overall correlation
r_all, p_all = stats.pearsonr(sst_total["age"], sst_total["sst_pct"])
print(f"\nOverall:  r = {r_all:.3f}, p = {p_all:.3f} (n={len(sst_total)})")

# Per-group correlations
stats_text = []
for dx in ["Control", "SCZ"]:
    grp = sst_total[sst_total["diagnosis"] == dx]
    r, p = stats.pearsonr(grp["age"], grp["sst_pct"])
    print(f"{dx:8s}: r = {r:.3f}, p = {p:.3f} (n={len(grp)})")
    stats_text.append((dx, r, p, len(grp)))

# ── Plot ────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))

for dx in ["Control", "SCZ"]:
    grp = sst_total[sst_total["diagnosis"] == dx]
    ax.scatter(grp["age"], grp["sst_pct"],
               c=DX_COLORS[dx], s=90, edgecolors="white",
               linewidths=0.5, alpha=0.85, label=dx, zorder=3)

# Overall regression line
slope, intercept = np.polyfit(sst_total["age"], sst_total["sst_pct"], 1)
age_range = np.linspace(sst_total["age"].min() - 2, sst_total["age"].max() + 2, 100)
ax.plot(age_range, slope * age_range + intercept,
        color="#888888", linewidth=1.5, linestyle="--", alpha=0.7, zorder=2)

# Label each point with sample_id
for _, row in sst_total.iterrows():
    ax.annotate(row["sample_id"].replace("Br", ""),
                (row["age"], row["sst_pct"]),
                fontsize=7, color="#666666",
                textcoords="offset points", xytext=(5, 3))

# Annotation box
anno_lines = [f"Overall: r = {r_all:.3f}, p = {p_all:.3f}, n = {len(sst_total)}"]
for dx, r, p, n in stats_text:
    anno_lines.append(f"{dx}: r = {r:.3f}, p = {p:.3f}, n = {n}")
ax.text(0.02, 0.98, "\n".join(anno_lines),
        transform=ax.transAxes, fontsize=12, verticalalignment="top",
        fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="#cccccc", alpha=0.9))

ax.set_xlabel("Subject Age (years)", fontsize=16)
ax.set_ylabel("Total SST Proportion (% cortical cells)", fontsize=16)
ax.set_title("SST Cell Proportion vs. Subject Age", fontsize=20, fontweight="bold")
ax.tick_params(labelsize=14)
ax.legend(fontsize=14, framealpha=0.9)
ax.grid(True, alpha=0.2)

fig.tight_layout()

out_path = os.path.join(PRESENTATION_DIR, "sst_proportion_by_age.png")
fig.savefig(out_path, dpi=200, bbox_inches="tight")
print(f"\nSaved → {out_path}")
plt.close(fig)
