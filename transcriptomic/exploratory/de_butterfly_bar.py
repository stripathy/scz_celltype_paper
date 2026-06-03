"""
Butterfly bar chart of DE gene counts per cell type, split up/down, with two
FDR tiers distinguished: a light bar = FDR < 0.10, a darker overlaid bar =
FDR < 0.05 (subset). Up extends right (orange), down extends left (blue).

DATA PROVENANCE:
  Input  data/DE_genes_all_cells_scz.csv  — meta-analytic snRNA-seq DE
         (cols used: genes, cell_type, estimate, padj). Update this path if
         the meta-analysis is re-run.
  Output results/00_de_butterfly.png/.pdf
         results/00_de_counts_by_celltype.csv  (up/down counts at both tiers)
         results/00_de_up_padj0.1.csv, results/00_de_down_padj0.1.csv
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

INPUT      = "data/DE_genes_all_cells_scz.csv"
OUT_FIG    = "results/00_de_butterfly.png"
OUT_PDF    = "results/00_de_butterfly.pdf"
OUT_COUNTS = "results/00_de_counts_by_celltype.csv"
OUT_UP     = "results/00_de_up_padj0.1.csv"
OUT_DOWN   = "results/00_de_down_padj0.1.csv"

# Okabe-Ito up/down, with a light tint for the FDR<0.1 tier
UP_DARK,   UP_LIGHT   = "#D55E00", "#F2B58C"   # vermillion / light
DOWN_DARK, DOWN_LIGHT = "#0072B2", "#9FCAE6"   # blue / light

df = pd.read_csv(INPUT)
df["direction"] = np.where(df["estimate"] > 0, "up", "down")


def counts_at(thr):
    s = df[df["padj"] < thr]
    c = s.groupby(["cell_type", "direction"]).size().unstack(fill_value=0)
    for col in ("up", "down"):
        if col not in c.columns:
            c[col] = 0
    return c[["up", "down"]]


c10 = counts_at(0.10)
c05 = counts_at(0.05).reindex(c10.index, fill_value=0)  # FDR<0.05 is a subset of FDR<0.1

m = pd.DataFrame({
    "up_fdr10":   c10["up"],   "up_fdr05":   c05["up"],
    "down_fdr10": c10["down"], "down_fdr05": c05["down"],
})
m["total_fdr10"] = m["up_fdr10"] + m["down_fdr10"]
m = m.sort_values("total_fdr10", ascending=True)

# ---- Plot ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11, 10))
y = np.arange(len(m))
bar_kw = dict(height=0.74, edgecolor="black", linewidth=0.5)

# Light bars = FDR < 0.10 (full extent)
ax.barh(y,  m["up_fdr10"],   color=UP_LIGHT,   **bar_kw)
ax.barh(y, -m["down_fdr10"], color=DOWN_LIGHT, **bar_kw)
# Dark bars = FDR < 0.05 (overlaid subset, from 0)
ax.barh(y,  m["up_fdr05"],   color=UP_DARK,    **bar_kw)
ax.barh(y, -m["down_fdr05"], color=DOWN_DARK,  **bar_kw)

xmax = max(m["up_fdr10"].max(), m["down_fdr10"].max()) * 1.20
inner_min = xmax * 0.06   # only label the dark (FDR<0.05) count if the bar can fit it

for i, row in enumerate(m.itertuples()):
    # Outer label = FDR<0.10 total, at the bar tip
    if row.up_fdr10 > 0:
        ax.text(row.up_fdr10 + xmax * 0.012, i, str(int(row.up_fdr10)),
                va="center", ha="left", fontsize=12)
    if row.down_fdr10 > 0:
        ax.text(-row.down_fdr10 - xmax * 0.012, i, str(int(row.down_fdr10)),
                va="center", ha="right", fontsize=12)
    # Inner label = FDR<0.05 count, white, inside the dark bar when wide enough
    if row.up_fdr05 >= inner_min:
        ax.text(row.up_fdr05 - xmax * 0.012, i, str(int(row.up_fdr05)),
                va="center", ha="right", fontsize=10.5, color="white", fontweight="bold")
    if row.down_fdr05 >= inner_min:
        ax.text(-row.down_fdr05 + xmax * 0.012, i, str(int(row.down_fdr05)),
                va="center", ha="left", fontsize=10.5, color="white", fontweight="bold")

ax.set_yticks(y)
ax.set_yticklabels(m.index, fontsize=14)
ax.axvline(0, color="black", linewidth=0.9)
ax.set_xlim(-xmax, xmax)
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(int(v))}"))  # show |count|
ax.tick_params(axis="x", labelsize=13)
ax.set_xlabel("Number of DE genes  (← down   |   up →)", fontsize=17)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="y", length=0)

# Legend: 4 entries (direction x FDR tier)
handles = [
    mpatches.Patch(facecolor=UP_DARK,    edgecolor="black", label="Up, FDR < 0.05"),
    mpatches.Patch(facecolor=UP_LIGHT,   edgecolor="black", label="Up, FDR < 0.10"),
    mpatches.Patch(facecolor=DOWN_DARK,  edgecolor="black", label="Down, FDR < 0.05"),
    mpatches.Patch(facecolor=DOWN_LIGHT, edgecolor="black", label="Down, FDR < 0.10"),
]
ax.legend(handles=handles, loc="lower right", fontsize=12, frameon=False,
          title="white number = FDR < 0.05 count", title_fontsize=11)

tot10 = (int(m["up_fdr10"].sum()), int(m["down_fdr10"].sum()))
tot05 = (int(m["up_fdr05"].sum()), int(m["down_fdr05"].sum()))
ax.set_title(
    "Differentially expressed genes by cell type\n"
    f"FDR < 0.10: {tot10[0]} up / {tot10[1]} down      "
    f"FDR < 0.05: {tot05[0]} up / {tot05[1]} down",
    fontsize=18, pad=14,
)

plt.tight_layout()
plt.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
plt.savefig(OUT_PDF, bbox_inches="tight")
print(f"Saved figure: {OUT_FIG} / .pdf")

# ---- Tables ----------------------------------------------------------------
m.to_csv(OUT_COUNTS)
print(f"Saved counts table: {OUT_COUNTS}")

sig = df[df["padj"] < 0.10].copy()
sig[sig["direction"] == "up"].sort_values(["cell_type", "padj"]).to_csv(OUT_UP, index=False)
sig[sig["direction"] == "down"].sort_values(["cell_type", "padj"]).to_csv(OUT_DOWN, index=False)
print(f"Saved gene lists (FDR<0.1): {OUT_UP}, {OUT_DOWN}")
print(f"\nTotals  FDR<0.10: {tot10[0]} up / {tot10[1]} down   |   "
      f"FDR<0.05: {tot05[0]} up / {tot05[1]} down")
