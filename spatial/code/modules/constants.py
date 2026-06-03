"""
Shared constants for the SCZ Xenium project.

These constants are used by both pipeline scripts (code/pipeline/) and
analysis scripts (code/analysis/). Centralizing them here eliminates the
cross-layer dependency where pipeline scripts imported from analysis/config.py.

analysis/config.py re-exports these via `from modules.constants import *`,
so existing analysis scripts continue to work unchanged.
"""

# ──────────────────────────────────────────────────────────────────────
# Sample metadata
# ──────────────────────────────────────────────────────────────────────

# Hardcoded diagnosis mapping derived from code/modules/metadata.py.
# Avoids needing the Excel metadata file at import time.
SAMPLE_TO_DX = {
    'Br5400': 'Control', 'Br2039': 'SCZ', 'Br2719': 'Control',
    'Br1113': 'Control', 'Br5373': 'SCZ', 'Br5590': 'SCZ',
    'Br6432': 'Control', 'Br5314': 'Control', 'Br5436': 'Control',
    'Br8772': 'SCZ', 'Br8433': 'Control', 'Br5746': 'SCZ',
    'Br5588': 'SCZ', 'Br5973': 'SCZ', 'Br6032': 'SCZ',
    'Br6437': 'SCZ', 'Br5639': 'Control', 'Br6389': 'Control',
    'Br5622': 'Control', 'Br1139': 'SCZ', 'Br2421': 'SCZ',
    'Br5931': 'Control', 'Br6496': 'SCZ', 'Br8667': 'Control',
}

CONTROL_SAMPLES = sorted(k for k, v in SAMPLE_TO_DX.items() if v == "Control")
SCZ_SAMPLES = sorted(k for k, v in SAMPLE_TO_DX.items() if v == "SCZ")
EXCLUDE_SAMPLES = set()  # No samples excluded; Br2039 (WM-heavy) included — improves snRNAseq concordance

# ──────────────────────────────────────────────────────────────────────
# Cortical layers
# ──────────────────────────────────────────────────────────────────────

CORTICAL_LAYERS = {"L1", "L2/3", "L4", "L5", "L6"}

# ──────────────────────────────────────────────────────────────────────
# Cell class classification
# ──────────────────────────────────────────────────────────────────────

CLASS_COLORS = {
    "Glutamatergic": "#00ADF8",   # Allen Institute / SEA-AD reference
    "GABAergic":     "#F05A28",   # Allen Institute / SEA-AD reference
    "Non-neuronal":  "#808080",   # Allen Institute / SEA-AD reference
}

SUBCLASS_TO_CLASS = {
    # Glutamatergic
    "L2/3 IT": "Glutamatergic", "L4 IT": "Glutamatergic",
    "L5 IT": "Glutamatergic", "L5 ET": "Glutamatergic",
    "L5/6 NP": "Glutamatergic", "L6 IT": "Glutamatergic",
    "L6 IT Car3": "Glutamatergic", "L6 CT": "Glutamatergic",
    "L6b": "Glutamatergic",
    # GABAergic
    "Lamp5": "GABAergic", "Lamp5 Lhx6": "GABAergic",
    "Sncg": "GABAergic", "Vip": "GABAergic",
    "Pax6": "GABAergic", "Chandelier": "GABAergic",
    "Pvalb": "GABAergic", "Sst": "GABAergic",
    "Sst Chodl": "GABAergic",
    # Non-neuronal
    "Astrocyte": "Non-neuronal", "Oligodendrocyte": "Non-neuronal",
    "OPC": "Non-neuronal", "Microglia-PVM": "Non-neuronal",
    "Endothelial": "Non-neuronal", "VLMC": "Non-neuronal",
    "SMC": "Non-neuronal", "Pericyte": "Non-neuronal",
}
