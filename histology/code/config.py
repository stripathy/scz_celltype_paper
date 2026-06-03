"""Constants, paths, exclusion lists, and channel mappings."""

from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────
PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"
STEREO_DIR = PROJECT_DIR / "dwight data 2"
COORD_DIR = PROJECT_DIR / "coordinates"

CELL_COUNTS_CSV = DATA_DIR / "Cell_counts_NU.csv"
DIAGNOSIS_XLSX = DATA_DIR / "full cell counts(Excel).xlsx"
DEMOGRAPHICS_CSV = DATA_DIR / "pTable with correct med info.csv"

# ── Subject ID corrections ─────────────────────────────────────────────
# Mismatches between Cell_counts_NU.csv and the authoritative Excel mapping
ID_FIXES = {"683": "863", "1159": "1157", "1449": "1444", "1381": "1391"}

# ── Exclusions ─────────────────────────────────────────────────────────
# R-section (SST/VIP): Dwight's exclusions + 1143 (missing data) + 1367 (bogus)
EXCLUDE_R = {"1088", "1143", "1153", "1188", "1226", "1240", "1341", "1367"}
# L-section (PV/PYR): 789 (Dwight's) + 1367 (bogus L data)
EXCLUDE_L = {"789", "1367"}

# ── Channel mapping ────────────────────────────────────────────────────
# Left section:  488 = PYR (SLC17A7), 568 = PV (PVALB)
# Right section: 488 = SST,           568 = VIP
CHANNEL_MAP = {
    ("L", "488"): "PYR",
    ("L", "568"): "PV",
    ("R", "488"): "SST",
    ("R", "568"): "VIP",
}

# ── Stereological parameters ──────────────────────────────────────────
FOV_SIDE_UM = 333  # field of view side length in µm
FOV_AREA_UM2 = FOV_SIDE_UM ** 2  # 110889 µm²
FOV_AREA_MM2 = FOV_AREA_UM2 / 1e6  # 0.110889 mm²

# Conversion factor: multiply raw count by this to get cells/mm²
COUNT_TO_DENSITY = 1.0 / FOV_AREA_MM2  # ≈ 9.018 cells/mm² per count

# ── Layer assignments ─────────────────────────────────────────────────
# Sites 1-10 = L2/3, sites 11-20 = L5/6
L23_SITES = set(range(1, 11))
L56_SITES = set(range(11, 21))

# ── Visualization ─────────────────────────────────────────────────────
DIAG_ORDER = ["Control", "Bipolar", "MDD", "SCHIZ"]
DIAG_COLORS = {
    "Control": "#4C72B0",
    "Bipolar": "#DD8452",
    "MDD": "#55A868",
    "SCHIZ": "#C44E52",
}
