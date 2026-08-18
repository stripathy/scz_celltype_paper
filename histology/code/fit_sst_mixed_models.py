#!/usr/bin/env python3
"""Fit the SST-density mixed models behind the RNAscope supplementary figure and
write their outputs as CSVs, so the figure script never hard-codes a statistic.

Model (per-frame counts, VIP-QC-filtered subjects), fitted three times:
    SST ~ diagnosis + age + sex + PMI + (1 | subject)
on (i) all counting frames, (ii) L2/3 frames only, (iii) L5/6 frames only.
Diagnosis is treatment-coded against Control, so each coefficient is the
difference in SST-positive cells per counting frame relative to controls.
Coefficients are also expressed per mm^2 (frame area = 0.110889 mm^2), which
is the unit plotted.

Reproduces the values quoted in the manuscript (SCZ: beta = -0.597, P = 0.071
overall; -0.762 / 0.126 in L2/3; -0.434 / 0.180 in L5/6).

Usage (from histology/):
    python3 code/fit_sst_mixed_models.py

Outputs (histology/results/):
    sst_mixed_models.csv    one row per (model, diagnosis) with beta, SE, P, CI
    sst_subject_density.csv per-subject mean SST density (all frames, L2/3, L5/6)
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data" / "sst_analysis_data.csv"
OUT_MODELS = ROOT / "results" / "sst_mixed_models.csv"
OUT_SUBJ = ROOT / "results" / "sst_subject_density.csv"

FOV_AREA_MM2 = 0.110889          # 333 x 333 um counting frame
FORMULA = "SST ~ C(diagnosis, Treatment('Control')) + age + sex_binary + PMI"
DIAGNOSES = ["SCHIZ", "Bipolar", "MDD"]
LABEL = {"SCHIZ": "SCZ", "Bipolar": "Bipolar", "MDD": "MDD"}


def fit(df: pd.DataFrame, model_label: str) -> list[dict]:
    df = df.reset_index(drop=True)
    m = smf.mixedlm(FORMULA, df, groups=df["subject"]).fit(reml=True)
    rows = []
    for dx in DIAGNOSES:
        k = f"C(diagnosis, Treatment('Control'))[T.{dx}]"
        b, se, p = m.params[k], m.bse[k], m.pvalues[k]
        rows.append(dict(
            model=model_label, diagnosis=LABEL[dx],
            n_subjects=df["subject"].nunique(), n_frames=len(df),
            beta_count=b, se_count=se, p=p,
            beta_density=b / FOV_AREA_MM2, se_density=se / FOV_AREA_MM2,
            ci_lo_density=(b - 1.96 * se) / FOV_AREA_MM2,
            ci_hi_density=(b + 1.96 * se) / FOV_AREA_MM2,
        ))
    return rows


def main() -> None:
    d = pd.read_csv(DATA)
    d = d[d["vip_pass"]].dropna(subset=["SST", "age", "sex_binary", "PMI"]).copy()
    d["diagnosis"] = d["diagnosis"].astype(str)
    print(f"VIP-QC-filtered frames: {len(d)}; subjects per diagnosis: "
          f"{d.groupby('diagnosis')['subject'].nunique().to_dict()}")

    rows = fit(d, "All layers")
    for layer in ["L2/3", "L5/6"]:
        rows += fit(d[d["layer"] == layer], layer)
    models = pd.DataFrame(rows)
    models.to_csv(OUT_MODELS, index=False)
    print(models[["model", "diagnosis", "beta_count", "se_count", "p", "beta_density"]]
          .round(3).to_string(index=False))

    # Per-subject densities for the boxplots (Control and SCZ only are plotted,
    # but all diagnoses are written so the file is reusable).
    d["density"] = d["SST"] / FOV_AREA_MM2
    all_layers = (d.groupby(["subject", "diagnosis"])["density"].mean()
                  .reset_index().assign(layer="All layers"))
    by_layer = (d.groupby(["subject", "diagnosis", "layer"])["density"].mean()
                .reset_index())
    subj = pd.concat([all_layers, by_layer], ignore_index=True)
    subj["diagnosis"] = subj["diagnosis"].map(lambda x: LABEL.get(x, x))
    subj.to_csv(OUT_SUBJ, index=False)
    ctrl_scz = subj[subj.diagnosis.isin(["Control", "SCZ"])]
    print(ctrl_scz.groupby(["layer", "diagnosis"])["density"]
          .agg(n="size", mean="mean", median="median").round(1))
    print(f"wrote {OUT_MODELS.name}, {OUT_SUBJ.name}")


if __name__ == "__main__":
    main()
