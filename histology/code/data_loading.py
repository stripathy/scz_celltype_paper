"""Load and merge cell count, diagnosis, and demographic data."""

import re
import numpy as np
import pandas as pd

from . import config


def load_cell_counts(csv_path=None):
    """Parse Cell_counts_NU.csv into a tidy per-site DataFrame.

    Returns DataFrame with columns:
        subject (str), section (str: 'L'/'R'), site (int: 1-20),
        layer (str: 'L23'/'L56'), ch488 (float), ch568 (float),
        area_um2 (float, only on first row per subject-section)
    """
    csv_path = csv_path or config.CELL_COUNTS_CSV
    raw = pd.read_csv(csv_path)

    rows = []
    for _, r in raw.iterrows():
        subj_str = str(r["Subject"]).strip()
        m = re.match(r"(\d+)\s*-\s*([LR])\s*-\s*(\d+)", subj_str)
        if not m:
            continue
        subj_id, section, site = m.group(1), m.group(2), int(m.group(3))

        # Apply ID corrections
        subj_id = config.ID_FIXES.get(subj_id, subj_id)

        layer = "L23" if site in config.L23_SITES else "L56"
        ch488 = pd.to_numeric(r["488 Channel"], errors="coerce")
        ch568 = pd.to_numeric(r["568 Channel"], errors="coerce")
        area = pd.to_numeric(r.get("Area (um^2)", np.nan), errors="coerce")

        rows.append({
            "subject": subj_id,
            "section": section,
            "site": site,
            "layer": layer,
            "ch488": ch488,
            "ch568": ch568,
            "area_um2": area,
        })

    return pd.DataFrame(rows)


def load_diagnosis_mapping(xlsx_path=None):
    """Load subject-to-diagnosis mapping from the authoritative Excel file.

    Returns dict: {subject_str: diagnosis_str}
    """
    xlsx_path = xlsx_path or config.DIAGNOSIS_XLSX
    df = pd.read_excel(xlsx_path, sheet_name="full cell counts")
    return dict(zip(df["Subject"].astype(str), df["Subject.Group"]))


def load_dwight_summary(xlsx_path=None):
    """Load Dwight's per-subject cell count summary from the Excel file.

    Returns DataFrame with columns: subject, PV, PYR_23, PYR_56, SST, VIP
    """
    xlsx_path = xlsx_path or config.DIAGNOSIS_XLSX
    df = pd.read_excel(xlsx_path, sheet_name="full cell counts")
    df = df[["Subject", "PV", "PYR_23", "PYR_56", "SST", "VIP"]].copy()
    df["Subject"] = df["Subject"].astype(str)
    return df.rename(columns={"Subject": "subject"})


def load_demographics(csv_path=None):
    """Load age, sex, PMI per subject from demographics CSV.

    Returns DataFrame with columns: subject, age, sex, PMI (one row per subject)
    """
    csv_path = csv_path or config.DEMOGRAPHICS_CSV
    df = pd.read_csv(csv_path)
    # CSV has multiple rows per subject (one per Cell.Type); deduplicate
    demo = (
        df.groupby("HU.")
        .first()
        .reset_index()
        [["HU.", "Age", "Sex", "PMI"]]
        .rename(columns={"HU.": "subject", "Age": "age", "Sex": "sex", "PMI": "PMI"})
    )
    demo["subject"] = demo["subject"].astype(str)
    return demo


def build_analysis_df(exclude_R=None, exclude_L=None):
    """Build the merged analysis DataFrame with one row per site.

    Applies exclusions, merges diagnosis and demographics, assigns cell type
    columns based on channel mapping.

    Returns two DataFrames: (df_R, df_L)
        df_R: R-section data with SST, VIP columns
        df_L: L-section data with PV, PYR columns
    Both have columns: subject, section, site, layer, diagnosis, age, sex,
                        sex_binary, <cell_type_1>, <cell_type_2>
    """
    exclude_R = exclude_R or config.EXCLUDE_R
    exclude_L = exclude_L or config.EXCLUDE_L

    # Load raw data
    counts = load_cell_counts()
    diag_map = load_diagnosis_mapping()
    demo = load_demographics()

    # Add diagnosis
    counts["diagnosis"] = counts["subject"].map(diag_map)

    # Split by section
    df_R = counts[counts["section"] == "R"].copy()
    df_L = counts[counts["section"] == "L"].copy()

    # Apply exclusions
    df_R = df_R[~df_R["subject"].isin(exclude_R)].copy()
    df_L = df_L[~df_L["subject"].isin(exclude_L)].copy()

    # Assign cell type columns (raw counts per FOV; FOV is constant 333² µm²)
    df_R["SST"] = df_R["ch488"]
    df_R["VIP"] = df_R["ch568"]
    df_L["PYR"] = df_L["ch488"]
    df_L["PV"] = df_L["ch568"]

    # Drop subjects with missing diagnosis
    df_R = df_R.dropna(subset=["diagnosis"])
    df_L = df_L.dropna(subset=["diagnosis"])

    # Drop subjects with all-missing cell counts
    for df in [df_R, df_L]:
        cell_cols = ["SST", "VIP"] if "SST" in df.columns else ["PV", "PYR"]
        subj_valid = df.groupby("subject")[cell_cols].apply(
            lambda g: g.notna().any().any()
        )
        valid_subjects = subj_valid[subj_valid].index
        df.drop(df[~df["subject"].isin(valid_subjects)].index, inplace=True)

    # Merge demographics
    dfs = [df_R, df_L]
    merged_dfs = []
    for df in dfs:
        df = df.drop(columns=[c for c in ["age", "sex", "sex_binary"] if c in df.columns],
                     errors="ignore")
        df = df.merge(demo[["subject", "age", "sex"]], on="subject", how="left")
        df["sex_binary"] = (df["sex"] == "M").astype(int)
        merged_dfs.append(df)
    df_R, df_L = merged_dfs

    # Print summary
    n_R = df_R["subject"].nunique()
    n_L = df_L["subject"].nunique()
    print(f"R-section (SST/VIP): {n_R} subjects, {len(df_R)} site-rows")
    print(f"L-section (PV/PYR): {n_L} subjects, {len(df_L)} site-rows")

    return df_R, df_L
