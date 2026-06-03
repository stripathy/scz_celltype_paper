"""
Subject-level metadata loading for SCZ Xenium study.

Handles the non-standard Excel file format (requires manual XML parsing
due to non-standard namespace) and provides diagnosis mappings.
"""

import zipfile
import xml.etree.ElementTree as ET
import pandas as pd


def load_subject_metadata(xlsx_path):
    """
    Load subject metadata from the study Excel file.

    Uses manual XML parsing because the file has a non-standard
    namespace that causes openpyxl and pandas to fail.

    Parameters
    ----------
    xlsx_path : str
        Path to the metadata Excel file (sample_metadata.xlsx).

    Returns
    -------
    pd.DataFrame
        Subject-level metadata with columns including 'Brain Number (ID)',
        'PrimaryDx', 'Sex', 'Age', 'RIN', 'PMI', 'Xenium'.
    """
    with zipfile.ZipFile(xlsx_path) as z:
        # Parse shared strings
        with z.open('xl/sharedStrings.xml') as f:
            ss_root = ET.fromstring(f.read())
        ns = ss_root.tag.split('}')[0] + '}' if '}' in ss_root.tag else ''
        strings = []
        for si in ss_root.iter(f'{ns}t'):
            if si.text:
                strings.append(si.text)

        # Parse worksheet
        with z.open('xl/worksheets/sheet1.xml') as f:
            ws_root = ET.fromstring(f.read())
        ns_ws = ws_root.tag.split('}')[0] + '}' if '}' in ws_root.tag else ''

        rows_data = []
        for row in ws_root.iter(f'{ns_ws}row'):
            row_vals = []
            for cell in row.iter(f'{ns_ws}c'):
                typ = cell.get('t', '')
                v_elem = cell.find(f'{ns_ws}v')
                if v_elem is not None and v_elem.text:
                    if typ == 's':
                        val = strings[int(v_elem.text)]
                    else:
                        val = v_elem.text
                else:
                    val = ''
                row_vals.append(val)
            rows_data.append(row_vals)

    header = rows_data[0]
    df = pd.DataFrame(rows_data[1:], columns=header)
    return df


def get_diagnosis_map(xlsx_path):
    """
    Build a sample_id -> diagnosis mapping from metadata.

    Parameters
    ----------
    xlsx_path : str
        Path to the metadata Excel file.

    Returns
    -------
    dict
        {sample_id: 'SCZ' or 'Control'}, e.g. {'Br8667': 'Control'}.
    """
    df = load_subject_metadata(xlsx_path)

    # Find the brain ID and diagnosis columns
    id_col = [c for c in df.columns if 'Brain' in c or 'ID' in c]
    dx_col = [c for c in df.columns if 'Dx' in c or 'Diagnosis' in c]

    if not id_col or not dx_col:
        raise ValueError(f"Could not find ID/Dx columns in: {df.columns.tolist()}")

    id_col = id_col[0]
    dx_col = dx_col[0]

    mapping = {}
    for _, row in df.iterrows():
        sid = str(row[id_col]).strip()
        dx = str(row[dx_col]).strip()
        if not sid or sid == 'nan':
            continue
        if 'Schiz' in dx:
            mapping[sid] = 'SCZ'
        elif 'Control' in dx:
            mapping[sid] = 'Control'
    return mapping


def get_subject_info(xlsx_path):
    """
    Get full subject info for Xenium samples only.

    Returns
    -------
    pd.DataFrame
        One row per Xenium sample with columns: sample_id, diagnosis,
        sex, age, pmi.
    """
    df = load_subject_metadata(xlsx_path)

    # Filter to Xenium samples
    xenium_col = [c for c in df.columns if 'Xenium' in c]
    if xenium_col:
        df = df[df[xenium_col[0]].str.upper() == 'YES'].copy()

    id_col = [c for c in df.columns if 'Brain' in c][0]
    dx_col = [c for c in df.columns if 'Dx' in c][0]

    records = []
    for _, row in df.iterrows():
        sid = str(row[id_col]).strip()
        dx = 'SCZ' if 'Schiz' in str(row[dx_col]) else 'Control'
        pmi_val = row.get('PMI', '')
        records.append({
            'sample_id': sid,
            'diagnosis': dx,
            'sex': str(row.get('Sex', '')),
            'age': float(row.get('Age', 0)) if row.get('Age', '') else None,
            'pmi': float(pmi_val) if pmi_val not in ('', None, 'nan') else None,
        })
    return pd.DataFrame(records)
