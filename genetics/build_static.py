#!/usr/bin/env python3
"""
build_static.py — Generate a fully static version of the SCZ enrichment webapp.

Produces a self-contained static site in site/ that can be deployed to Netlify
or any static hosting. All data is pre-computed and baked into HTML/JSON files.

Routes:
  /               -> site/index.html  (main enrichment explorer)
  /drivers/{ct}   -> site/drivers/{ct}/index.html  (gene driver scatter per type)

Usage:
  python build_static.py          # Build into site/
  netlify deploy --dir=site       # Deploy to Netlify
"""
import csv
import json
import math
import os
import sys
import time
import shutil
from pathlib import Path
from urllib.parse import quote

import numpy as np
from statsmodels.stats.multitest import multipletests

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SITE_DIR = os.path.join(BASE_DIR, "site")

# ---------------------------------------------------------------------------
# Reuse data loading from app.py
# ---------------------------------------------------------------------------
sys.path.insert(0, BASE_DIR)
from app import prepare_data, DriverDataProvider, HTML_TEMPLATE, DRIVER_HTML_TEMPLATE


def build():
    t0 = time.time()

    # Clean and create output directory
    if os.path.exists(SITE_DIR):
        shutil.rmtree(SITE_DIR)
    os.makedirs(SITE_DIR)
    os.makedirs(os.path.join(SITE_DIR, "drivers"), exist_ok=True)

    # ── 1. Build main index.html ──────────────────────────────────────────
    print("Loading enrichment data...")
    app_data = prepare_data()
    data_json = json.dumps(app_data)
    n_types = len(app_data["points"])
    n_cond = len(app_data["conditional"])
    print(f"  {n_types} cell types, {n_cond} conditional results")

    print("Building index.html...")
    html = HTML_TEMPLATE.replace("%%DATA_JSON%%", data_json)
    html = html.replace("%%N_TYPES%%", str(n_types))
    html = html.replace("%%N_COND%%", str(n_cond))
    with open(os.path.join(SITE_DIR, "index.html"), "w") as f:
        f.write(html)
    print(f"  site/index.html ({len(html) // 1024} KB)")

    # ── 2. Build gene driver pages ────────────────────────────────────────
    print("Loading gene driver data (this may take a minute)...")
    driver_provider = DriverDataProvider()

    if driver_provider.available:
        cell_types = sorted(driver_provider.cell_types)
        print(f"  Building {len(cell_types)} driver pages...")

        # Write driver JSON data files separately (not embedded in HTML)
        # This avoids duplicating the ~10KB HTML template 598 times
        data_dir = os.path.join(SITE_DIR, "drivers", "_data")
        os.makedirs(data_dir, exist_ok=True)

        # Subsample background points to keep JSON files small
        # 2000 background points is plenty for visual context in a scatter
        MAX_BG = 1000
        rng = np.random.RandomState(42)

        for i, ct in enumerate(cell_types):
            ct_json_str = driver_provider.get_json(ct)
            if ct_json_str is None:
                continue

            # Parse, subsample background, re-serialize
            ct_data = json.loads(ct_json_str)
            n_bg = len(ct_data.get("bg_spec", []))
            if n_bg > MAX_BG:
                idx = sorted(rng.choice(n_bg, MAX_BG, replace=False))
                ct_data["bg_spec"] = [ct_data["bg_spec"][j] for j in idx]
                ct_data["bg_logp"] = [ct_data["bg_logp"][j] for j in idx]

            ct_json_str = json.dumps(ct_data, separators=(",", ":"))

            # URL-safe filename
            safe_name = ct.replace("/", "_")

            # Write JSON data file
            with open(os.path.join(data_dir, f"{safe_name}.json"), "w") as f:
                f.write(ct_json_str)

            # Write minimal HTML page that fetches the data
            ct_dir = os.path.join(SITE_DIR, "drivers", safe_name)
            os.makedirs(ct_dir, exist_ok=True)

            # Create a thin HTML wrapper that loads data via fetch
            loader_html = DRIVER_HTML_TEMPLATE.replace(
                "%%DRIVER_JSON%%",
                "null"  # placeholder — will be loaded via fetch
            ).replace("%%CELL_TYPE%%", ct)

            # The template has: const D = null;
            # Wrap all code after that in a function, then fetch + call it
            fetch_code = (
                'let D = null;\n'
                'fetch("/drivers/_data/' + safe_name + '.json")\n'
                '  .then(r => r.json())\n'
                '  .then(d => { D = d; initDriver(); })\n'
                '  .catch(e => console.error("Failed to load driver data", e));\n'
                'function initDriver() {'
            )
            loader_html = loader_html.replace("const D = null;", fetch_code)
            # Close the initDriver function before </script>
            loader_html = loader_html.replace(
                "</script>\n</body>",
                "}\n</script>\n</body>"
            )

            with open(os.path.join(ct_dir, "index.html"), "w") as f:
                f.write(loader_html)

            if (i + 1) % 100 == 0:
                print(f"    {i + 1}/{len(cell_types)} done...")

        print(f"  {len(cell_types)} driver pages built")

        # Update driver links in index.html to use URL-safe names
        index_path = os.path.join(SITE_DIR, "index.html")
        with open(index_path, "r") as f:
            index_html = f.read()

        # Replace the encodeURIComponent call with a version that also
        # replaces "/" with "_" to match our static file structure.
        # Original JS: encodeURIComponent(p.cell_type) + '" style=
        # Target JS:   encodeURIComponent(p.cell_type.replace(/\//g, '_')) + '/" style=
        index_html = index_html.replace(
            "encodeURIComponent(p.cell_type)",
            "encodeURIComponent(p.cell_type.replace(/\\//g, '_'))"
        )
        # Add trailing slash for directory-style URLs:
        # Original: + '" style=  ->  + '/" style=
        index_html = index_html.replace(
            """+ '" style="color:#e74c3c""",
            """+ '/" style="color:#e74c3c"""
        )

        with open(index_path, "w") as f:
            f.write(index_html)
    else:
        print("  WARNING: Gene driver data not available — skipping driver pages")

    # ── 3. Write _redirects for Netlify ───────────────────────────────────
    # Handle URL-encoded names (e.g., "L2%2F3 IT_2" -> "L2_3 IT_2")
    redirects = []
    if driver_provider.available:
        for ct in sorted(driver_provider.cell_types):
            safe_name = ct.replace("/", "_")
            if safe_name != ct:
                # Redirect the original URL-encoded path to the safe path
                encoded = quote(ct, safe="")
                redirects.append(f"/drivers/{encoded}/* /drivers/{safe_name}/:splat 200")
                redirects.append(f"/drivers/{ct}/* /drivers/{safe_name}/:splat 200")

    with open(os.path.join(SITE_DIR, "_redirects"), "w") as f:
        for r in redirects:
            f.write(r + "\n")

    # ── 4. Write netlify.toml-compatible headers ──────────────────────────
    # (Already in project root as netlify.toml)

    # ── 5. Summary ────────────────────────────────────────────────────────
    total_size = 0
    n_files = 0
    for root, dirs, files in os.walk(SITE_DIR):
        for fname in files:
            fpath = os.path.join(root, fname)
            total_size += os.path.getsize(fpath)
            n_files += 1

    elapsed = time.time() - t0
    print(f"\n{'=' * 50}")
    print(f"Static site built in {elapsed:.1f}s")
    print(f"  Output:  {SITE_DIR}/")
    print(f"  Files:   {n_files}")
    print(f"  Size:    {total_size / 1024 / 1024:.1f} MB")
    print(f"\nTo preview locally:")
    print(f"  cd site && python3 -m http.server 8080")
    print(f"\nTo deploy to Netlify:")
    print(f"  netlify deploy --dir=site --prod")


if __name__ == "__main__":
    os.chdir(BASE_DIR)
    build()
