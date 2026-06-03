#!/usr/bin/env python3
"""
SCZ Cell-Type Enrichment Explorer
Interactive web application for exploring schizophrenia cell-type enrichment results.
Run: python app.py
Then open: http://localhost:8050
"""

import csv
import json
import math
import os
import numpy as np
from http.server import HTTPServer, SimpleHTTPRequestHandler
from urllib.parse import urlparse, parse_qs
from statsmodels.stats.multitest import multipletests

PORT = int(os.environ.get("PORT", 8050))
BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def load_csv(path):
    """Load a CSV file and return list of dicts."""
    full = os.path.join(BASE_DIR, path)
    with open(full, "r") as f:
        reader = csv.DictReader(f)
        return [row for row in reader]


def load_json(path):
    """Load a JSON file."""
    full = os.path.join(BASE_DIR, path)
    with open(full, "r") as f:
        return json.load(f)


def prepare_data():
    """Load and merge all data files into JSON-serializable structures."""
    enrichment = load_csv("results/tables/rbh_combined_enrichment_app.csv")

    # Conditional analysis: load if available, otherwise empty
    cond_path = os.path.join(BASE_DIR, "results/tables/combined_conditional_analysis.csv")
    conditional = load_csv("results/tables/combined_conditional_analysis.csv") if os.path.exists(cond_path) else []

    correlation_path = os.path.join(BASE_DIR, "results/intermediates/siletti_cluster_seaad_correlation.csv")
    correlation = load_csv("results/intermediates/siletti_cluster_seaad_correlation.csv") if os.path.exists(correlation_path) else []

    classification = load_csv("results/intermediates/siletti_cluster_level_classification.csv")
    type_order = load_csv("results/intermediates/rbh_combined_type_order_info.csv")

    # Load region summary for Siletti types
    region_summary_path = os.path.join(BASE_DIR, "results/intermediates/siletti_cluster_region_summary.csv")
    region_summary = []
    if os.path.exists(region_summary_path):
        region_summary = load_csv("results/intermediates/siletti_cluster_region_summary.csv")

    # Load MetaNeighbor results (full 461 version)
    mn_seaad_path = os.path.join(BASE_DIR, "results/intermediates/rbh_seaad_to_siletti_metaneighbor.csv")
    mn_siletti_path = os.path.join(BASE_DIR, "results/intermediates/rbh_siletti_to_seaad_metaneighbor.csv")
    # Fall back to older files if RBH versions don't exist
    if not os.path.exists(mn_seaad_path):
        mn_seaad_path = os.path.join(BASE_DIR, "results/intermediates/seaad_to_siletti_metaneighbor.csv")
    if not os.path.exists(mn_siletti_path):
        mn_siletti_path = os.path.join(BASE_DIR, "results/intermediates/siletti_to_seaad_metaneighbor.csv")
    seaad_to_siletti = []
    siletti_to_seaad_mn = []
    if os.path.exists(mn_seaad_path):
        seaad_to_siletti = load_csv(mn_seaad_path)
    if os.path.exists(mn_siletti_path):
        siletti_to_seaad_mn = load_csv(mn_siletti_path)

    # Load SEA-AD supertype colors (also defines canonical order)
    seaad_colors_path = os.path.join(BASE_DIR, "data/seaad_supertype_colors.json")
    seaad_colors = {}
    if os.path.exists(seaad_colors_path):
        seaad_colors = load_json("data/seaad_supertype_colors.json")

    # Build lookup dicts
    corr_lookup = {r["siletti_cluster"]: r for r in correlation}
    class_lookup = {r["cluster_name"]: r for r in classification}
    region_lookup = {r["cluster_name"]: r for r in region_summary}
    seaad_corr_lookup = {r["seaad_supertype"]: r for r in seaad_to_siletti}
    siletti_mn_lookup = {r["siletti_cluster"]: r for r in siletti_to_seaad_mn}

    # Build ordered cell type list: use seaad_colors order for SEA-AD, then Siletti from type_order
    # SEA-AD order from the colors JSON
    seaad_ordered = list(seaad_colors.keys())
    # Find any SEA-AD enrichment types not in the color JSON and append them
    seaad_in_enrichment = {r["cell_type"] for r in enrichment if r["source"] == "SEA-AD"}
    seaad_missing = [ct for ct in seaad_in_enrichment if ct not in seaad_colors]
    seaad_ordered = seaad_ordered + sorted(seaad_missing)
    # Siletti order from existing type_order file
    siletti_from_order = [r["cell_type"] for r in type_order if r["source"] == "Siletti"]
    all_ordered = seaad_ordered + siletti_from_order
    # Give SEA-AD and Siletti roughly equal x-axis space
    # SEA-AD types are spaced so they span the same width as Siletti types
    n_seaad = len(seaad_ordered)
    n_siletti = len(siletti_from_order)
    order_map = {}
    if n_seaad > 0 and n_siletti > 0:
        # Scale SEA-AD spacing so both halves are equal width
        seaad_spacing = n_siletti / n_seaad  # e.g., ~2.5x wider spacing
        for i, ct in enumerate(seaad_ordered):
            order_map[ct] = i * seaad_spacing
        seaad_end = (n_seaad - 1) * seaad_spacing + seaad_spacing  # gap before Siletti
        for i, ct in enumerate(siletti_from_order):
            order_map[ct] = seaad_end + i
    else:
        order_map = {ct: i for i, ct in enumerate(all_ordered)}

    # Also build class/subclass maps from type_order
    class_map = {r["cell_type"]: r["class"] for r in type_order}
    subclass_map = {r["cell_type"]: r["subclass"] for r in type_order}
    # Add any SEA-AD types from colors that might not be in type_order
    for ct in seaad_ordered:
        if ct not in class_map:
            class_map[ct] = "SEA-AD"
            subclass_map[ct] = ""

    # Merge enrichment with ordering and class info
    points = []
    for row in enrichment:
        ct = row["cell_type"]
        p = float(row["p_value"])
        neglog10p = -math.log10(p) if p > 0 else 300
        beta = float(row["beta"])
        se = float(row["se"])
        t_stat = float(row["t_stat"])
        p_bonf = float(row["p_bonferroni"])
        p_fdr = float(row["p_fdr"])
        source = row["source"]
        supercluster = row["supercluster"]
        cell_class = class_map.get(ct, "Siletti")
        subclass = subclass_map.get(ct, "")
        x_order = order_map.get(ct, 9999)

        sig_bonf = p_bonf < 0.05
        sig_fdr = p_fdr < 0.05

        # Assign color: SEA-AD from JSON, Siletti default purple
        color = seaad_colors.get(ct, "#999999") if source == "SEA-AD" else "#8b5cf6"

        point = {
            "cell_type": ct,
            "source": source,
            "supercluster": supercluster,
            "cell_class": cell_class,
            "subclass": subclass,
            "color": color,
            "beta": round(beta, 4),
            "se": round(se, 4),
            "t_stat": round(t_stat, 4),
            "p_value": p,
            "p_fdr": p_fdr,
            "p_bonferroni": p_bonf,
            "neglog10p": round(neglog10p, 4),
            "x_order": x_order,
            "sig_bonf": sig_bonf,
            "sig_fdr": sig_fdr,
        }

        # Add SEA-AD -> Siletti match info (MetaNeighbor AUROC)
        if source == "SEA-AD":
            seaad_corr = seaad_corr_lookup.get(ct, {})
            point["best_siletti_match"] = seaad_corr.get("best_siletti_match", "")
            point["best_siletti_auroc"] = seaad_corr.get("best_siletti_auroc", "")
            point["second_siletti_match"] = seaad_corr.get("second_siletti_match", "")
            point["second_siletti_auroc"] = seaad_corr.get("second_siletti_auroc", "")

        # Add Siletti-specific info
        if source == "Siletti":
            # Strip "Siletti_" prefix to match keys in correlation/classification files
            siletti_key = ct.replace("Siletti_", "", 1) if ct.startswith("Siletti_") else ct
            corr_info = corr_lookup.get(siletti_key, {})
            class_info = class_lookup.get(siletti_key, {})
            # MetaNeighbor match (preferred)
            mn_info = siletti_mn_lookup.get(siletti_key, {})
            point["best_seaad_match"] = mn_info.get("best_seaad_match", corr_info.get("best_seaad_match", ""))
            point["best_seaad_auroc"] = mn_info.get("best_seaad_auroc", "")
            point["best_correlation"] = corr_info.get("best_correlation", "")
            point["frac_neocortical"] = class_info.get("frac_neocortical", "")
            point["n_total"] = class_info.get("n_total", "")
            point["n_neocortical"] = class_info.get("n_neocortical", "")
            point["n_subcortical"] = class_info.get("n_subcortical", "")
            point["include_siletti"] = class_info.get("include_siletti", "")
            # Region info
            region_info = region_lookup.get(siletti_key, {})
            point["top_regions"] = region_info.get("top_regions", "")
            point["top_divisions"] = region_info.get("top_divisions", "")

        points.append(point)

    # Sort by x_order for the plot
    points.sort(key=lambda p: p["x_order"])

    # Build group boundaries for x-axis labels
    # For SEA-AD: use subclass; for Siletti: use supercluster
    # Walk the ordered points and find contiguous runs of the same group
    group_bands = []
    current_group = None
    current_source = None
    band_start = None
    for pt in points:
        if pt["source"] == "SEA-AD":
            grp = pt["subclass"] or "Other"
        else:
            grp = pt["supercluster"] or "Other"
        key = (pt["source"], grp)
        if key != (current_source, current_group):
            if current_group is not None:
                group_bands.append({
                    "label": current_group,
                    "source": current_source,
                    "x_start": band_start,
                    "x_end": prev_x,
                })
            current_group = grp
            current_source = pt["source"]
            band_start = pt["x_order"]
        prev_x = pt["x_order"]
    # Final group
    if current_group is not None:
        group_bands.append({
            "label": current_group,
            "source": current_source,
            "x_start": band_start,
            "x_end": prev_x,
        })

    # Conditional analysis data
    cond_data = []
    for row in conditional:
        cond_data.append({
            "supertype": row["supertype"],
            "subclass": row["subclass"],
            "marginal_p": float(row["marginal_p"]),
            "marginal_t": float(row["marginal_t"]),
            "conditional_p": float(row["conditional_p"]),
            "conditional_t": float(row["conditional_t"]),
            "conditional_beta": float(row["conditional_beta"]),
            "source": row["source"],
            "supercluster": row["supercluster"],
        })

    # Compute thresholds
    n = len(points)
    bonf_threshold = -math.log10(0.05 / n) if n > 0 else 5
    # Find the FDR threshold: max p among significant FDR types
    fdr_pvals = [p["p_value"] for p in points if p["sig_fdr"]]
    fdr_threshold = -math.log10(max(fdr_pvals)) if fdr_pvals else 2

    return {
        "points": points,
        "conditional": cond_data,
        "bonf_threshold": round(bonf_threshold, 4),
        "fdr_threshold": round(fdr_threshold, 4),
        "n_types": n,
        "group_bands": group_bands,
    }


class DriverDataProvider:
    """Lazy-loading provider for gene driver scatter data.

    Loads specificity + GWAS data once, then computes per-cell-type JSON
    on demand. Handles both SEA-AD types (gene-symbol specificity) and
    Siletti clusters (ENTREZ-indexed specificity).
    """

    def __init__(self):
        import pandas as pd
        spec_path = os.path.join(BASE_DIR, "results/intermediates/seaad_supertype_specificity.csv")
        gwas_path = os.path.join(BASE_DIR, "results/intermediates/gwas_matched.csv")
        siletti_spec_path = os.path.join(BASE_DIR, "results/intermediates/reprocessed_siletti_specificity_entrez.csv")
        magma_path = os.path.join(BASE_DIR,
            "linking_cell_types_to_brain_phenotypes/Example_results/"
            "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.no_heading.step2.genes.out")

        self.available = os.path.exists(spec_path) and os.path.exists(gwas_path)
        if not self.available:
            return

        # SEA-AD: gene-symbol indexed, pre-merged with GWAS
        self.seaad_spec = pd.read_csv(spec_path, index_col=0)
        self.gwas = pd.read_csv(gwas_path)
        _, self.gwas_fdr, _, _ = multipletests(self.gwas["P"].values, method="fdr_bh")
        self.seaad_spec_matched = self.seaad_spec.loc[self.gwas["SYMBOL"].values]
        self.seaad_types = set(self.seaad_spec_matched.columns)

        # Siletti: ENTREZ-indexed, merge with MAGMA directly
        self.siletti_types = set()
        self.siletti_spec_matched = None
        self.siletti_gwas = None
        self.siletti_gwas_fdr = None
        self.siletti_gene_names = None

        if os.path.exists(siletti_spec_path) and os.path.exists(magma_path):
            siletti_spec = pd.read_csv(siletti_spec_path, index_col=0)
            magma_df = pd.read_csv(magma_path, sep=r"\s+")

            # Build ENTREZ → Symbol mapping for display
            gene_loc_path = os.path.join(BASE_DIR,
                "linking_cell_types_to_brain_phenotypes/Data/"
                "NCBI37.3.gene.loc.extendedMHCexcluded")
            gene_loc = pd.read_csv(gene_loc_path, sep="\t", header=None,
                names=["ENTREZ", "CHR", "START", "STOP", "STRAND", "SYMBOL"])
            entrez_to_sym = dict(zip(gene_loc["ENTREZ"].astype(int), gene_loc["SYMBOL"]))

            # Intersect Siletti specificity with MAGMA
            overlap = sorted(set(siletti_spec.index) & set(magma_df["GENE"]))
            self.siletti_spec_matched = siletti_spec.loc[overlap]
            self.siletti_gwas = magma_df.set_index("GENE").loc[overlap].reset_index()
            _, self.siletti_gwas_fdr, _, _ = multipletests(
                self.siletti_gwas["P"].values, method="fdr_bh")
            # Map ENTREZ to symbols for display
            self.siletti_gene_names = [
                entrez_to_sym.get(int(eid), str(int(eid)))
                for eid in self.siletti_gwas["GENE"].values
            ]
            self.siletti_types = set(self.siletti_spec_matched.columns)

        self.cell_types = self.seaad_types | self.siletti_types
        self._cache = {}

    def __contains__(self, ct):
        return self.available and ct in self.cell_types

    def get_json(self, ct):
        """Return JSON string for a cell type's driver data."""
        if ct in self._cache:
            return self._cache[ct]

        if ct in self.seaad_types:
            result = self._build_json_seaad(ct)
        elif ct in self.siletti_types:
            result = self._build_json_siletti(ct)
        else:
            return None

        self._cache[ct] = result
        return result

    def _build_json_seaad(self, ct):
        s = self.seaad_spec_matched[ct].values
        gwas = self.gwas
        gwas_fdr = self.gwas_fdr
        gene_names = gwas["SYMBOL"].values
        return self._build_json(ct, s, gwas, gwas_fdr, gene_names)

    def _build_json_siletti(self, ct):
        s = self.siletti_spec_matched[ct].values
        gwas = self.siletti_gwas
        gwas_fdr = self.siletti_gwas_fdr
        gene_names = self.siletti_gene_names
        return self._build_json(ct, s, gwas, gwas_fdr, gene_names)

    def _build_json(self, ct, s, gwas, gwas_fdr, gene_names):
        nonzero = s > 0
        spec_threshold = float(np.percentile(s[nonzero], 90)) if nonzero.sum() > 0 else 0
        logp = -np.log10(np.clip(gwas["P"].values, 1e-300, 1))

        genes = []
        bg_spec = []
        bg_logp = []
        n_total = len(gwas)

        for i in range(n_total):
            is_fdr = bool(gwas_fdr[i] < 0.05)
            is_hs = bool(s[i] > spec_threshold)

            if is_fdr or is_hs:
                genes.append({
                    "gene": str(gene_names[i]) if hasattr(gene_names, '__getitem__') else str(gene_names),
                    "spec": round(float(s[i]), 6),
                    "logp": round(float(logp[i]), 4),
                    "z": round(float(gwas.iloc[i]["ZSTAT"]), 3),
                    "p": float(gwas.iloc[i]["P"]),
                    "fdr": round(float(gwas_fdr[i]), 4),
                    "is_fdr": is_fdr,
                    "is_hs": is_hs,
                })
            else:
                bg_spec.append(round(float(s[i]), 5))
                bg_logp.append(round(float(logp[i]), 3))

        data = {
            "cell_type": ct,
            "genes": genes,
            "bg_spec": bg_spec,
            "bg_logp": bg_logp,
            "spec_threshold": round(spec_threshold, 6),
            "n_genes": n_total,
            "n_gwas_fdr": int((gwas_fdr < 0.05).sum()),
        }

        return json.dumps(data)


DRIVER_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Gene Drivers: %%CELL_TYPE%%</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
* { margin: 0; padding: 0; box-sizing: border-box; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background: #f5f6fa; color: #2c3e50; }
header { background: #1a1a2e; color: white; padding: 16px 24px; display: flex; align-items: center; justify-content: space-between; }
header h1 { font-size: 22px; font-weight: 600; }
header a { color: #8892b0; text-decoration: none; font-size: 14px; }
header a:hover { color: white; }
.container { max-width: 1400px; margin: 0 auto; padding: 16px; }
.plot-section { background: white; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); padding: 16px; margin-bottom: 16px; }
#driver-plot { width: 100%; height: 700px; }
.summary { display: flex; gap: 24px; flex-wrap: wrap; margin-bottom: 16px; }
.summary-card { background: white; border-radius: 8px; padding: 16px 24px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); text-align: center; }
.summary-card .number { font-size: 28px; font-weight: 700; color: #e74c3c; }
.summary-card .label { font-size: 13px; color: #666; }
.table-section { background: white; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); padding: 16px; }
table { width: 100%; border-collapse: collapse; font-size: 13px; }
table th { background: #f8f9fa; padding: 8px 6px; text-align: left; cursor: pointer; border-bottom: 2px solid #dee2e6; font-size: 12px; position: sticky; top: 0; }
table th:hover { background: #e9ecef; }
table td { padding: 6px; border-bottom: 1px solid #eee; }
table tr:hover { background: #f0f7ff; }
table tr.driver { background: #fff3f3; }
.table-wrap { max-height: 500px; overflow-y: auto; }
.controls { display: flex; gap: 12px; align-items: center; margin-bottom: 12px; }
.controls input { padding: 6px 10px; border: 1px solid #ddd; border-radius: 4px; font-size: 13px; width: 200px; }
.controls label { font-size: 13px; font-weight: 600; color: #555; }
.controls select { padding: 6px 10px; border: 1px solid #ddd; border-radius: 4px; font-size: 13px; }
</style>
</head>
<body>
<header>
  <div>
    <h1>Gene Drivers: %%CELL_TYPE%%</h1>
    <span style="color:#8892b0; font-size:13px">Genes with high cell-type specificity + GWAS FDR &lt; 0.05</span>
  </div>
  <a href="/">&larr; Back to enrichment explorer</a>
</header>
<div class="container">
  <div class="summary">
    <div class="summary-card"><div class="number" id="n-drivers">-</div><div class="label">Driver genes</div></div>
    <div class="summary-card"><div class="number" id="n-gwas">-</div><div class="label">GWAS FDR-sig genes</div></div>
    <div class="summary-card"><div class="number" id="n-hs">-</div><div class="label">High specificity genes</div></div>
    <div class="summary-card"><div class="number" style="color:#1a1a2e" id="n-total">-</div><div class="label">Total genes tested</div></div>
  </div>
  <div class="plot-section">
    <div id="driver-plot"></div>
  </div>
  <div class="table-section">
    <h2 style="margin-bottom:12px">Driver Genes <span id="table-count" style="background:#1a1a2e;color:white;padding:2px 8px;border-radius:10px;font-size:12px;margin-left:8px"></span></h2>
    <div class="controls">
      <label>Show:</label>
      <select id="gene-filter">
        <option value="drivers">Driver genes only</option>
        <option value="gwas_fdr">All GWAS FDR-sig</option>
        <option value="high_spec">All high specificity</option>
        <option value="all">All genes</option>
      </select>
      <label>Search:</label>
      <input type="text" id="gene-search" placeholder="Search gene name...">
    </div>
    <div class="table-wrap">
      <table id="gene-table">
        <thead><tr>
          <th data-col="gene">Gene</th>
          <th data-col="spec">Specificity</th>
          <th data-col="logp">-log10(p)</th>
          <th data-col="z">GWAS Z</th>
          <th data-col="p">GWAS p</th>
          <th data-col="fdr">GWAS FDR</th>
          <th data-col="score">Driver Score</th>
        </tr></thead>
        <tbody id="gene-tbody"></tbody>
      </table>
    </div>
  </div>
</div>
<script>
const D = %%DRIVER_JSON%%;
const genes = D.genes;        // relevant genes (FDR-sig or high-spec) with full metadata
const bgSpec = D.bg_spec;     // background gene x-coords
const bgLogp = D.bg_logp;     // background gene y-coords
const specThreshold = D.spec_threshold;

// Precompute derived fields
genes.forEach(g => {
  g.is_driver = g.is_fdr && g.is_hs;
  g.score = g.spec * g.logp;
});

const nDrivers = genes.filter(g => g.is_driver).length;
const nGwas = genes.filter(g => g.is_fdr).length;
const nHs = genes.filter(g => g.is_hs).length;
document.getElementById('n-drivers').textContent = nDrivers;
document.getElementById('n-gwas').textContent = D.n_gwas_fdr;
document.getElementById('n-hs').textContent = nHs;
document.getElementById('n-total').textContent = D.n_genes.toLocaleString();

// Find FDR line position
const fdrGenes = genes.filter(g => g.is_fdr);
const fdrLine = fdrGenes.length > 0 ? Math.min(...fdrGenes.map(g => g.logp)) : 2;

// Build scatter
function buildDriverPlot() {
  // Background (gray) genes — minimal data, no hover
  const bgTrace = {
    x: bgSpec, y: bgLogp,
    mode: 'markers', type: 'scatter',
    name: 'Other (' + bgSpec.length + ')',
    hoverinfo: 'skip',
    marker: { color: '#dddddd', size: 4, symbol: 'circle', opacity: 0.25 },
  };

  // Relevant gene categories
  const categories = [
    { label: 'High specificity only', filter: g => g.is_hs && !g.is_fdr, color: '#3498db', size: 8, symbol: 'circle', opacity: 0.5 },
    { label: 'GWAS FDR-sig only', filter: g => g.is_fdr && !g.is_hs, color: '#e67e22', size: 6, symbol: 'circle', opacity: 0.4 },
    { label: 'Driver (both)', filter: g => g.is_driver, color: '#e74c3c', size: 12, symbol: 'diamond', opacity: 0.85 },
  ];

  const traces = [bgTrace];
  categories.forEach(cat => {
    const subset = genes.filter(cat.filter);
    traces.push({
      x: subset.map(g => g.spec),
      y: subset.map(g => g.logp),
      text: subset.map(g =>
        '<b>' + g.gene + '</b><br>' +
        'Specificity: ' + g.spec.toFixed(4) + '<br>' +
        'GWAS Z: ' + g.z + '<br>' +
        'GWAS p: ' + (g.p < 0.001 ? g.p.toExponential(2) : g.p.toPrecision(3)) + '<br>' +
        'GWAS FDR: ' + (g.fdr < 0.001 ? g.fdr.toExponential(2) : g.fdr.toPrecision(3)) + '<br>' +
        'Driver score: ' + g.score.toFixed(4)
      ),
      mode: 'markers',
      type: 'scatter',
      name: cat.label + ' (' + subset.length + ')',
      hoverinfo: 'text',
      marker: { color: cat.color, size: cat.size, symbol: cat.symbol, opacity: cat.opacity,
                line: { width: cat.label.includes('Driver') ? 0.5 : 0, color: '#333' } },
    });
  });

  // Label top driver genes
  const topDrivers = genes.filter(g => g.is_driver).sort((a,b) => b.score - a.score).slice(0, 30);
  const annotations = topDrivers.map(g => ({
    x: g.spec, y: g.logp, text: '<b>' + g.gene + '</b>', showarrow: false,
    xanchor: 'left', yanchor: 'bottom',
    font: { size: 12, color: '#c0392b' },
    bgcolor: 'rgba(255,255,255,0.7)',
    borderpad: 1,
    xshift: 6, yshift: 4,
  }));

  // Also label high-specificity nominally significant genes
  const nomGenes = genes.filter(g => g.is_hs && !g.is_fdr && g.p < 0.05)
    .sort((a,b) => b.score - a.score).slice(0, 10);
  nomGenes.forEach(g => {
    annotations.push({
      x: g.spec, y: g.logp, text: g.gene, showarrow: false,
      xanchor: 'left', yanchor: 'bottom',
      font: { size: 11, color: '#2980b9' },
      bgcolor: 'rgba(255,255,255,0.6)',
      borderpad: 1,
      xshift: 6, yshift: 3,
    });
  });

  const xMax = Math.max(...genes.map(g => g.spec)) * 1.05;
  const yMax = Math.max(...genes.map(g => g.logp)) * 1.05;

  const layout = {
    xaxis: { title: D.cell_type + ' specificity score', zeroline: false },
    yaxis: { title: '-log10(GWAS p-value)', zeroline: false },
    shapes: [
      { type: 'line', x0: 0, x1: xMax, y0: fdrLine, y1: fdrLine, line: { color: '#e74c3c', width: 1.5, dash: 'dash' } },
      { type: 'line', x0: 0, x1: xMax, y0: -Math.log10(0.05 / D.n_genes), y1: -Math.log10(0.05 / D.n_genes), line: { color: '#8e44ad', width: 1.5, dash: 'dashdot' } },
      { type: 'line', x0: 0, x1: xMax, y0: -Math.log10(0.05), y1: -Math.log10(0.05), line: { color: '#f39c12', width: 1, dash: 'dot' } },
      { type: 'line', x0: specThreshold, x1: specThreshold, y0: 0, y1: yMax, line: { color: 'gray', width: 1, dash: 'dash' } },
    ],
    annotations: [
      { x: xMax*0.99, y: fdrLine, text: 'GWAS FDR 0.05', showarrow: false, xanchor: 'right', font: {size:11, color:'#e74c3c'} },
      { x: xMax*0.99, y: -Math.log10(0.05 / D.n_genes), text: 'Bonferroni', showarrow: false, xanchor: 'right', font: {size:11, color:'#8e44ad'} },
      { x: specThreshold, y: yMax*0.98, text: 'Top 10% spec', showarrow: false, xanchor: 'left', textangle: -90, font: {size:10, color:'gray'} },
      ...annotations
    ],
    showlegend: true,
    legend: { x: 0.01, y: 0.99 },
    margin: { t: 20, b: 60, l: 60, r: 20 },
    hovermode: 'closest',
    plot_bgcolor: '#fafafa',
  };

  Plotly.newPlot('driver-plot', traces, layout, { responsive: true });
}

// Gene table
let tblSortCol = 'score';
let tblSortAsc = false;

function fmtP(p) {
  if (p < 0.001) return p.toExponential(2);
  return p.toPrecision(3);
}

function buildGeneTable() {
  const filter = document.getElementById('gene-filter').value;
  const search = document.getElementById('gene-search').value.toLowerCase();

  let filtered = genes;
  if (filter === 'drivers') filtered = filtered.filter(g => g.is_driver);
  else if (filter === 'gwas_fdr') filtered = filtered.filter(g => g.is_fdr);
  else if (filter === 'high_spec') filtered = filtered.filter(g => g.is_hs);
  if (search) filtered = filtered.filter(g => g.gene.toLowerCase().includes(search));

  filtered = [...filtered].sort((a, b) => {
    let va = a[tblSortCol], vb = b[tblSortCol];
    if (typeof va === 'string') { va = va.toLowerCase(); vb = vb.toLowerCase(); }
    return tblSortAsc ? (va < vb ? -1 : 1) : (va > vb ? -1 : 1);
  });

  document.getElementById('table-count').textContent = filtered.length;

  let html = '';
  filtered.forEach(g => {
    const cls = g.is_driver ? 'driver' : '';
    html += '<tr class="' + cls + '">' +
      '<td><b>' + g.gene + '</b></td>' +
      '<td>' + g.spec.toFixed(4) + '</td>' +
      '<td>' + g.logp.toFixed(2) + '</td>' +
      '<td>' + g.z + '</td>' +
      '<td>' + fmtP(g.p) + '</td>' +
      '<td>' + fmtP(g.fdr) + '</td>' +
      '<td>' + g.score.toFixed(4) + '</td></tr>';
  });
  document.getElementById('gene-tbody').innerHTML = html;
}

document.querySelectorAll('#gene-table th').forEach(th => {
  th.addEventListener('click', () => {
    const col = th.dataset.col;
    if (col === tblSortCol) tblSortAsc = !tblSortAsc;
    else { tblSortCol = col; tblSortAsc = col === 'gene'; }
    buildGeneTable();
  });
});
document.getElementById('gene-filter').addEventListener('change', buildGeneTable);
document.getElementById('gene-search').addEventListener('input', buildGeneTable);

buildDriverPlot();
buildGeneTable();
</script>
</body>
</html>"""


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>SCZ Cell-Type Enrichment Explorer</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
* { margin: 0; padding: 0; box-sizing: border-box; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background: #f5f6fa; color: #2c3e50; }
header { background: #1a1a2e; color: white; padding: 16px 24px; display: flex; align-items: center; justify-content: space-between; }
header h1 { font-size: 22px; font-weight: 600; }
header .subtitle { font-size: 13px; color: #8892b0; }
.container { max-width: 1600px; margin: 0 auto; padding: 16px; }
.plot-section { background: white; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); padding: 16px; margin-bottom: 16px; }
.controls { display: flex; gap: 12px; align-items: center; flex-wrap: wrap; margin-bottom: 12px; }
.controls label { font-size: 13px; font-weight: 600; color: #555; }
.controls select, .controls input { padding: 6px 10px; border: 1px solid #ddd; border-radius: 4px; font-size: 13px; }
.controls input[type=text] { width: 200px; }
#manhattan-plot { width: 100%; height: 620px; }
.panels { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; margin-bottom: 16px; }
@media (max-width: 1000px) { .panels { grid-template-columns: 1fr; } }
.panel { background: white; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); padding: 16px; }
.panel h2 { font-size: 16px; margin-bottom: 12px; color: #1a1a2e; border-bottom: 2px solid #e8e8e8; padding-bottom: 8px; }
.detail-grid { display: grid; grid-template-columns: 140px 1fr; gap: 4px 12px; font-size: 13px; }
.detail-grid .label { font-weight: 600; color: #666; }
.detail-grid .value { color: #2c3e50; }
.detail-grid .value.sig { color: #e74c3c; font-weight: 600; }
table { width: 100%; border-collapse: collapse; font-size: 12px; }
table th { background: #f8f9fa; padding: 8px 6px; text-align: left; cursor: pointer; user-select: none; border-bottom: 2px solid #dee2e6; font-size: 11px; position: sticky; top: 0; }
table th:hover { background: #e9ecef; }
table td { padding: 6px; border-bottom: 1px solid #eee; }
table tr:hover { background: #f0f7ff; }
table tr.sig-bonf { background: #fff3f3; }
table tr.sig-fdr { background: #fff8f0; }
table tr.selected { background: #d4edff !important; }
.table-wrap { max-height: 450px; overflow-y: auto; }
.cond-table { font-size: 12px; }
.cond-table td, .cond-table th { padding: 5px 6px; }
.cond-table tr.bonf { background: #fff3f3; }
.stat-badge { display: inline-block; padding: 2px 8px; border-radius: 10px; font-size: 11px; font-weight: 600; }
.stat-badge.bonf { background: #e74c3c; color: white; }
.stat-badge.fdr { background: #f39c12; color: white; }
.stat-badge.ns { background: #bdc3c7; color: #666; }
.sort-arrow { font-size: 10px; margin-left: 3px; }
.region-bar { height: 18px; border-radius: 3px; display: flex; overflow: hidden; margin: 4px 0; }
.region-bar .neo { background: #3498db; }
.region-bar .sub { background: #e67e22; }
.region-bar .other { background: #95a5a6; }
.legend-row { display: flex; gap: 16px; font-size: 12px; align-items: center; flex-wrap: wrap; }
.legend-item { display: flex; align-items: center; gap: 4px; }
.legend-dot { width: 10px; height: 10px; border-radius: 50%; }
#detail-panel-content { min-height: 100px; }
.empty-state { color: #aaa; font-style: italic; padding: 20px; text-align: center; }
.count-badge { background: #1a1a2e; color: white; padding: 2px 8px; border-radius: 10px; font-size: 12px; margin-left: 8px; }
</style>
</head>
<body>
<header>
  <div>
    <h1>SCZ Cell-Type Enrichment Explorer</h1>
    <div class="subtitle">GWAS enrichment across %%N_TYPES%% cell types from SEA-AD and Siletti atlases</div>
  </div>
  <div class="legend-row">
    <div class="legend-item"><div class="legend-dot" style="background:#3c7abf"></div> Excitatory (SEA-AD)</div>
    <div class="legend-item"><div class="legend-dot" style="background:#ffb84f"></div> Sst</div>
    <div class="legend-item"><div class="legend-dot" style="background:#be5330"></div> Pvalb</div>
    <div class="legend-item"><div class="legend-dot" style="background:#e1b0f7"></div> Vip</div>
    <div class="legend-item"><div class="legend-dot" style="background:#f9b7c3"></div> Lamp5</div>
    <div class="legend-item"><div class="legend-dot" style="background:#2baf3a"></div> Non-neuronal (SEA-AD)</div>
    <div class="legend-item"><div class="legend-dot" style="background:#8b5cf6"></div> Siletti</div>
  </div>
</header>
<div class="container">
  <div class="plot-section">
    <div class="controls">
      <label>Source:</label>
      <select id="filter-source"><option value="all">All</option><option value="SEA-AD">SEA-AD</option><option value="Siletti">Siletti</option></select>
      <label>Significance:</label>
      <select id="filter-sig"><option value="all">All</option><option value="fdr">FDR &lt; 0.05</option><option value="bonf">Bonferroni &lt; 0.05</option></select>
      <label>Search:</label>
      <input type="text" id="filter-search" placeholder="Type cell type name...">
      <span id="count-display"></span>
    </div>
    <div id="manhattan-plot"></div>
  </div>

  <div class="panels">
    <div class="panel">
      <h2>Selected Cell Type Details</h2>
      <div id="detail-panel-content"><div class="empty-state">Click a point on the plot or a table row to see details</div></div>
    </div>
    <div class="panel">
      <h2>Independently Enriched Types <span class="count-badge">%%N_COND%%</span></h2>
      <div class="table-wrap">
        <table class="cond-table" id="cond-table">
          <thead><tr>
            <th>Cell Type</th><th>Source</th><th>Supercluster</th>
            <th>Marginal t</th><th>Marginal p</th>
            <th>Conditional t</th><th>Conditional p</th><th>Cond. Beta</th>
          </tr></thead>
          <tbody id="cond-tbody"></tbody>
        </table>
      </div>
    </div>
  </div>

  <div class="plot-section">
    <h2 style="margin-bottom:12px;">Full Results Table <span id="table-count" class="count-badge"></span></h2>
    <div class="controls">
      <label>Source:</label>
      <select id="tbl-source"><option value="all">All</option><option value="SEA-AD">SEA-AD</option><option value="Siletti">Siletti</option></select>
      <label>Significance:</label>
      <select id="tbl-sig"><option value="all">All</option><option value="fdr">FDR &lt; 0.05</option><option value="bonf">Bonferroni &lt; 0.05</option></select>
      <label>Search:</label>
      <input type="text" id="tbl-search" placeholder="Filter by name...">
    </div>
    <div class="table-wrap">
      <table id="results-table">
        <thead><tr>
          <th data-col="cell_type">Cell Type <span class="sort-arrow"></span></th>
          <th data-col="source">Source <span class="sort-arrow"></span></th>
          <th data-col="supercluster">Supercluster <span class="sort-arrow"></span></th>
          <th data-col="subclass">Subclass <span class="sort-arrow"></span></th>
          <th data-col="beta">Beta <span class="sort-arrow"></span></th>
          <th data-col="t_stat">t-stat <span class="sort-arrow"></span></th>
          <th data-col="p_value">p-value <span class="sort-arrow"></span></th>
          <th data-col="p_fdr">FDR <span class="sort-arrow"></span></th>
          <th data-col="p_bonferroni">Bonferroni <span class="sort-arrow"></span></th>
          <th data-col="sig">Sig</th>
        </tr></thead>
        <tbody id="results-tbody"></tbody>
      </table>
    </div>
  </div>
</div>

<script>
// Data injected from Python
const DATA = %%DATA_JSON%%;
const allPoints = DATA.points;
const condData = DATA.conditional;
const bonfThreshold = DATA.bonf_threshold;
const fdrThreshold = DATA.fdr_threshold;
const groupBands = DATA.group_bands;

const CLASS_COLORS = {
  'Neuronal: Glutamatergic': '#2563eb',
  'Neuronal: GABAergic': '#dc2626',
  'Non-neuronal and Non-neural': '#16a34a',
  'Siletti': '#8b5cf6'
};

function getColor(cls) { return CLASS_COLORS[cls] || '#8b5cf6'; }

function fmtP(p) {
  if (p === 0) return '0';
  if (p < 1e-300) return '<1e-300';
  if (p < 0.001) return p.toExponential(2);
  return p.toPrecision(3);
}

function filterPoints(source, sig, search) {
  let pts = allPoints;
  if (source !== 'all') pts = pts.filter(p => p.source === source);
  if (sig === 'fdr') pts = pts.filter(p => p.sig_fdr);
  if (sig === 'bonf') pts = pts.filter(p => p.sig_bonf);
  if (search) {
    const q = search.toLowerCase();
    pts = pts.filter(p => p.cell_type.toLowerCase().includes(q) || p.supercluster.toLowerCase().includes(q) || p.subclass.toLowerCase().includes(q));
  }
  return pts;
}

// Manhattan plot
function buildPlot() {
  const source = document.getElementById('filter-source').value;
  const sig = document.getElementById('filter-sig').value;
  const search = document.getElementById('filter-search').value;
  const pts = filterPoints(source, sig, search);
  document.getElementById('count-display').textContent = pts.length + ' types shown';

  // Group by source for legend, use per-point colors
  const groups = {};
  pts.forEach(p => {
    const grp = p.source;
    if (!groups[grp]) groups[grp] = { x: [], y: [], text: [], marker: { color: [], size: [], symbol: [] }, name: grp, type: 'scatter', mode: 'markers', hoverinfo: 'text' };
    const g = groups[grp];
    g.x.push(p.x_order);
    g.y.push(p.neglog10p);
    g.marker.color.push(p.color);
    let tooltip = '<b>' + p.cell_type + '</b><br>' +
      'Source: ' + p.source + '<br>' +
      (p.supercluster ? 'Supercluster: ' + p.supercluster + '<br>' : '') +
      (p.subclass && p.source === 'SEA-AD' ? 'Subclass: ' + p.subclass + '<br>' : '') +
      'Beta: ' + p.beta + '  |  t-stat: ' + p.t_stat + '<br>' +
      'p: ' + fmtP(p.p_value) + '  |  FDR: ' + fmtP(p.p_fdr);
    if (p.source === 'Siletti') {
      if (p.top_divisions) tooltip += '<br><b>Brain regions:</b> ' + p.top_divisions;
      if (p.n_total) tooltip += '<br>Cells: ' + p.n_total + ' (' + (parseFloat(p.frac_neocortical)*100).toFixed(0) + '% neocortical)';
      if (p.best_seaad_match) {
        let matchStr = p.best_seaad_match;
        if (p.best_seaad_auroc) matchStr += ' (AUROC=' + parseFloat(p.best_seaad_auroc).toFixed(3) + ')';
        else if (p.best_correlation) matchStr += ' (r=' + parseFloat(p.best_correlation).toFixed(3) + ')';
        tooltip += '<br>Best SEA-AD match: ' + matchStr;
      }
    }
    if (p.source === 'SEA-AD') {
      tooltip += '<br><b>Region:</b> Medial temporal gyrus (MTG)';
      if (p.best_siletti_match) tooltip += '<br>Best Siletti match: ' + p.best_siletti_match + ' (AUROC=' + parseFloat(p.best_siletti_auroc).toFixed(3) + ')';
    }
    g.text.push(tooltip);
    g.marker.size.push(p.sig_bonf ? 10 : (p.sig_fdr ? 7 : 4));
    g.marker.symbol.push(p.sig_bonf ? 'diamond' : (p.sig_fdr ? 'diamond-open' : 'circle'));
  });

  const traces = Object.values(groups);

  // Threshold lines
  const xRange = [Math.min(...pts.map(p => p.x_order)) - 5, Math.max(...pts.map(p => p.x_order)) + 5];
  const shapes = [
    { type: 'line', x0: xRange[0], x1: xRange[1], y0: bonfThreshold, y1: bonfThreshold, line: { color: '#e74c3c', width: 1.5, dash: 'dash' } },
    { type: 'line', x0: xRange[0], x1: xRange[1], y0: fdrThreshold, y1: fdrThreshold, line: { color: '#f39c12', width: 1.5, dash: 'dot' } }
  ];

  // Annotations for threshold labels
  const annotations = [
    { x: xRange[1], y: bonfThreshold, text: 'Bonferroni', showarrow: false, xanchor: 'right', font: { size: 11, color: '#e74c3c' } },
    { x: xRange[1], y: fdrThreshold, text: 'FDR 0.05', showarrow: false, xanchor: 'right', font: { size: 11, color: '#f39c12' } }
  ];

  // Build x-axis group bands and labels
  // Filter group bands to only include groups with visible points
  const visibleXOrders = new Set(pts.map(p => p.x_order));
  const filteredBands = groupBands.filter(band => {
    // Check if any visible point falls within this band's x range
    for (const x of visibleXOrders) {
      if (x >= band.x_start && x <= band.x_end) return true;
    }
    return false;
  });

  filteredBands.forEach((band, i) => {
    // Alternating background shading
    if (i % 2 === 0) {
      shapes.push({
        type: 'rect',
        xref: 'x', yref: 'paper',
        x0: band.x_start - 0.5, x1: band.x_end + 0.5,
        y0: 0, y1: 1,
        fillcolor: 'rgba(0,0,0,0.03)',
        line: { width: 0 },
        layer: 'below',
      });
    }
    // Vertical separator line between groups
    if (i > 0) {
      const prevBand = filteredBands[i - 1];
      const sepX = (prevBand.x_end + band.x_start) / 2;
      shapes.push({
        type: 'line',
        xref: 'x', yref: 'paper',
        x0: sepX, x1: sepX,
        y0: 0, y1: 1,
        line: { color: 'rgba(0,0,0,0.08)', width: 0.5 },
        layer: 'below',
      });
    }
    // Group label annotation below plot
    const midX = (band.x_start + band.x_end) / 2;
    // Abbreviate long labels
    let label = band.label;
    if (label.length > 14) {
      // Shorten common long names
      const abbrevs = {
        'Committed oligodendrocyte precursor': 'COP',
        'Oligodendrocyte precursor': 'OPC',
        'Oligodendrocyte': 'Oligo',
        'Deep-layer corticothalamic and 6b': 'DL CT/6b',
        'Deep-layer intratelencephalic': 'DL IT',
        'Hippocampal CA1-3': 'Hipp CA1-3',
        'Hippocampal CA4': 'Hipp CA4',
        'Hippocampal dentate gyrus': 'Hipp DG',
        'Medium spiny neuron': 'MSN',
        'Eccentric medium spiny neuron': 'eMSN',
        'LAMP5-LHX6 and Chandelier': 'LAMP5/Chand',
        'MGE interneuron': 'MGE IN',
        'CGE interneuron': 'CGE IN',
        'Upper rhombic lip': 'Upper RL',
        'Cerebellar inhibitory': 'Cereb Inh',
        'Lower rhombic lip': 'Lower RL',
        'Mammillary body': 'Mammillary',
        'Thalamic excitatory': 'Thal Exc',
        'Midbrain-derived inhibitory': 'Midbrain Inh',
        'Amygdala excitatory': 'Amyg Exc',
        'Choroid plexus': 'Choroid',
        'Microglia-PVM': 'Micro-PVM',
        'Non-neuronal and Non-neural': 'Non-neuronal',
        'Miscellaneous': 'Misc',
        'Bergmann glia': 'Bergmann',
      };
      label = abbrevs[label] || label.substring(0, 12);
    }
    annotations.push({
      x: midX,
      y: -0.02,
      yref: 'paper',
      text: label,
      showarrow: false,
      xanchor: 'right',
      yanchor: 'top',
      font: { size: 9, color: band.source === 'SEA-AD' ? '#444' : '#6d28d9' },
      textangle: -50,
    });
  });

  const layout = {
    xaxis: { title: '', showticklabels: false, zeroline: false },
    yaxis: { title: '-log10(p-value)', zeroline: false },
    shapes: shapes,
    annotations: annotations,
    showlegend: true,
    legend: { orientation: 'h', y: 1.12, x: 0.5, xanchor: 'center' },
    margin: { t: 40, b: 120, l: 60, r: 20 },
    hovermode: 'closest',
    plot_bgcolor: '#fafafa'
  };

  Plotly.react('manhattan-plot', traces, layout, { responsive: true });

  // Click handler
  const plotDiv = document.getElementById('manhattan-plot');
  plotDiv.removeAllListeners && plotDiv.removeAllListeners('plotly_click');
  plotDiv.on('plotly_click', function(data) {
    if (data.points.length > 0) {
      const xVal = data.points[0].x;
      const yVal = data.points[0].y;
      const match = pts.find(p => p.x_order === xVal && Math.abs(p.neglog10p - yVal) < 0.001);
      if (match) showDetail(match);
    }
  });

  // Resize markers on zoom — scale inversely with x-axis range
  const fullXRange = xRange[1] - xRange[0];
  // Store original sizes so we always scale from baseline
  const origSizes = traces.map(t => t.marker.size.slice());
  plotDiv.on('plotly_relayout', function(ed) {
    let visibleRange = fullXRange;
    if (ed['xaxis.range[0]'] !== undefined) {
      visibleRange = ed['xaxis.range[1]'] - ed['xaxis.range[0]'];
    }
    // ed['xaxis.autorange'] or no range keys = reset to full view → scale=1
    const scale = Math.max(1, Math.min(5, fullXRange / visibleRange));
    for (let i = 0; i < traces.length; i++) {
      Plotly.restyle('manhattan-plot', {
        'marker.size': [origSizes[i].map(s => Math.round(s * scale))]
      }, [i]);
    }
  });
}

function showDetail(p) {
  let html = '<div class="detail-grid">';
  html += '<div class="label">Cell Type</div><div class="value"><b>' + p.cell_type + '</b></div>';
  html += '<div class="label">Source</div><div class="value">' + p.source + '</div>';
  html += '<div class="label">Class</div><div class="value">' + p.cell_class + '</div>';
  html += '<div class="label">Subclass</div><div class="value">' + (p.subclass || '-') + '</div>';
  html += '<div class="label">Supercluster</div><div class="value">' + (p.supercluster || '-') + '</div>';
  html += '<div class="label">Beta</div><div class="value">' + p.beta + ' (SE: ' + p.se + ')</div>';
  html += '<div class="label">t-statistic</div><div class="value">' + p.t_stat + '</div>';
  html += '<div class="label">p-value</div><div class="value' + (p.sig_fdr ? ' sig' : '') + '">' + fmtP(p.p_value) + '</div>';
  html += '<div class="label">FDR</div><div class="value' + (p.sig_fdr ? ' sig' : '') + '">' + fmtP(p.p_fdr) + '</div>';
  html += '<div class="label">Bonferroni</div><div class="value' + (p.sig_bonf ? ' sig' : '') + '">' + fmtP(p.p_bonferroni) + '</div>';
  html += '<div class="label">Significance</div><div class="value">';
  if (p.sig_bonf) html += '<span class="stat-badge bonf">Bonferroni</span> ';
  else if (p.sig_fdr) html += '<span class="stat-badge fdr">FDR</span> ';
  else html += '<span class="stat-badge ns">Not significant</span>';
  html += '</div>';

  // Source-specific region info
  if (p.source === 'SEA-AD') {
    html += '<div class="label" style="grid-column:1/-1;margin-top:8px;border-top:1px solid #eee;padding-top:8px"><b>Region Information</b></div>';
    html += '<div class="label">Brain Region</div><div class="value">Medial Temporal Gyrus (MTG)</div>';
    html += '<div class="label">Atlas</div><div class="value">SEA-AD Reference (5 neurotypical donors)</div>';
    if (p.best_siletti_match) {
      html += '<div class="label">Best Siletti Match</div><div class="value">' + p.best_siletti_match + ' (AUROC=' + parseFloat(p.best_siletti_auroc).toFixed(3) + ')</div>';
      if (p.second_siletti_match) {
        html += '<div class="label">2nd Siletti Match</div><div class="value">' + p.second_siletti_match + ' (AUROC=' + parseFloat(p.second_siletti_auroc).toFixed(3) + ')</div>';
      }
    }
  }

  if (p.source === 'Siletti') {
    html += '<div class="label" style="grid-column:1/-1;margin-top:8px;border-top:1px solid #eee;padding-top:8px"><b>Region Information</b></div>';
    if (p.n_total) {
      html += '<div class="label">Total Cells</div><div class="value">' + parseInt(p.n_total).toLocaleString() + '</div>';
      const frac = parseFloat(p.frac_neocortical) || 0;
      const fracSub = p.n_subcortical && p.n_total ? (parseInt(p.n_subcortical) / parseInt(p.n_total)) : 0;
      const fracOther = Math.max(0, 1 - frac - fracSub);
      html += '<div class="label">Region Dist.</div><div class="value">';
      html += '<div class="region-bar" style="width:260px">';
      html += '<div class="neo" style="width:' + (frac*100) + '%" title="Neocortical"></div>';
      html += '<div class="sub" style="width:' + (fracSub*100) + '%" title="Subcortical"></div>';
      html += '<div class="other" style="width:' + (fracOther*100) + '%" title="Non-neocortical cortex"></div>';
      html += '</div>';
      html += '<span style="font-size:11px"><span style="color:#3498db">■</span> Neo: ' + (frac*100).toFixed(1) + '% | <span style="color:#e67e22">■</span> Subcort: ' + (fracSub*100).toFixed(1) + '% | <span style="color:#95a5a6">■</span> Other ctx: ' + (fracOther*100).toFixed(1) + '%</span>';
      html += '</div>';
    }
    if (p.top_divisions) {
      html += '<div class="label">Brain Divisions</div><div class="value" style="font-size:12px">' + p.top_divisions + '</div>';
    }
    if (p.top_regions) {
      html += '<div class="label">Top ROIs</div><div class="value" style="font-size:12px">' + p.top_regions + '</div>';
    }
    if (p.best_seaad_match) {
      let matchDetail = p.best_seaad_match;
      if (p.best_seaad_auroc) matchDetail += ' (AUROC=' + parseFloat(p.best_seaad_auroc).toFixed(3) + ')';
      else if (p.best_correlation) matchDetail += ' (r=' + parseFloat(p.best_correlation).toFixed(3) + ')';
      html += '<div class="label">Best SEA-AD Match</div><div class="value">' + matchDetail + '</div>';
    }
  }

  // Gene driver link (both SEA-AD and Siletti types)
  html += '<div class="label" style="grid-column:1/-1;margin-top:8px;border-top:1px solid #eee;padding-top:8px"><b>Gene Drivers</b></div>';
  html += '<div class="label">Explore</div><div class="value"><a href="/drivers/' + encodeURIComponent(p.cell_type) + '" style="color:#e74c3c;font-weight:600;text-decoration:none;font-size:14px">View gene driver scatter \u2192</a></div>';

  // Check if in conditional results
  const condMatch = condData.find(c => c.supertype === p.cell_type);
  if (condMatch) {
    html += '<div class="label" style="grid-column:1/-1;margin-top:8px;border-top:1px solid #eee;padding-top:8px"><b>Conditional Analysis</b></div>';
    html += '<div class="label">Marginal t / p</div><div class="value">' + condMatch.marginal_t.toFixed(3) + ' / ' + fmtP(condMatch.marginal_p) + '</div>';
    html += '<div class="label">Conditional t / p</div><div class="value">' + condMatch.conditional_t.toFixed(3) + ' / ' + fmtP(condMatch.conditional_p) + '</div>';
    html += '<div class="label">Conditional Beta</div><div class="value">' + condMatch.conditional_beta.toFixed(4) + '</div>';
    html += '<div class="label">Status</div><div class="value"><span class="stat-badge bonf">Independently Enriched</span></div>';
  }

  html += '</div>';
  document.getElementById('detail-panel-content').innerHTML = html;

  // Highlight row in table
  document.querySelectorAll('#results-tbody tr').forEach(tr => {
    tr.classList.toggle('selected', tr.dataset.ct === p.cell_type);
  });
}

// Conditional analysis table
function buildCondTable() {
  const tbody = document.getElementById('cond-tbody');
  tbody.innerHTML = '';
  condData.sort((a, b) => a.conditional_p - b.conditional_p).forEach(c => {
    const isBonf = c.conditional_p < (0.05 / DATA.n_types);
    const tr = document.createElement('tr');
    if (isBonf) tr.className = 'bonf';
    tr.innerHTML =
      '<td><b>' + c.supertype + '</b></td>' +
      '<td>' + c.source + '</td>' +
      '<td>' + c.supercluster + '</td>' +
      '<td>' + c.marginal_t.toFixed(2) + '</td>' +
      '<td>' + fmtP(c.marginal_p) + '</td>' +
      '<td>' + c.conditional_t.toFixed(2) + '</td>' +
      '<td>' + fmtP(c.conditional_p) + '</td>' +
      '<td>' + c.conditional_beta.toFixed(2) + '</td>';
    tr.style.cursor = 'pointer';
    tr.addEventListener('click', () => {
      const match = allPoints.find(p => p.cell_type === c.supertype);
      if (match) showDetail(match);
    });
    tbody.appendChild(tr);
  });
}

// Results table
let sortCol = 'p_value';
let sortAsc = true;

function buildTable() {
  const source = document.getElementById('tbl-source').value;
  const sig = document.getElementById('tbl-sig').value;
  const search = document.getElementById('tbl-search').value;
  let pts = filterPoints(source, sig, search);

  // Sort
  pts = [...pts].sort((a, b) => {
    let va = a[sortCol], vb = b[sortCol];
    if (typeof va === 'string') { va = va.toLowerCase(); vb = vb.toLowerCase(); }
    if (va < vb) return sortAsc ? -1 : 1;
    if (va > vb) return sortAsc ? 1 : -1;
    return 0;
  });

  document.getElementById('table-count').textContent = pts.length;

  const tbody = document.getElementById('results-tbody');
  // Build HTML in bulk for performance
  let html = '';
  pts.forEach(p => {
    let cls = '';
    if (p.sig_bonf) cls = 'sig-bonf';
    else if (p.sig_fdr) cls = 'sig-fdr';
    let badge = '<span class="stat-badge ns">NS</span>';
    if (p.sig_bonf) badge = '<span class="stat-badge bonf">Bonf</span>';
    else if (p.sig_fdr) badge = '<span class="stat-badge fdr">FDR</span>';
    html += '<tr class="' + cls + '" data-ct="' + p.cell_type + '">' +
      '<td><b>' + p.cell_type + '</b></td>' +
      '<td>' + p.source + '</td>' +
      '<td>' + p.supercluster + '</td>' +
      '<td>' + (p.subclass || '-') + '</td>' +
      '<td>' + p.beta.toFixed(2) + '</td>' +
      '<td>' + p.t_stat.toFixed(2) + '</td>' +
      '<td>' + fmtP(p.p_value) + '</td>' +
      '<td>' + fmtP(p.p_fdr) + '</td>' +
      '<td>' + fmtP(p.p_bonferroni) + '</td>' +
      '<td>' + badge + '</td></tr>';
  });
  tbody.innerHTML = html;

  // Click handlers for rows
  tbody.querySelectorAll('tr').forEach(tr => {
    tr.style.cursor = 'pointer';
    tr.addEventListener('click', () => {
      const ct = tr.dataset.ct;
      const match = allPoints.find(p => p.cell_type === ct);
      if (match) showDetail(match);
    });
  });

  // Update sort arrows
  document.querySelectorAll('#results-table th').forEach(th => {
    const col = th.dataset.col;
    const arrow = th.querySelector('.sort-arrow');
    if (arrow) arrow.textContent = col === sortCol ? (sortAsc ? ' \\u25B2' : ' \\u25BC') : '';
  });
}

// Sort headers
document.querySelectorAll('#results-table th[data-col]').forEach(th => {
  th.addEventListener('click', () => {
    const col = th.dataset.col;
    if (col === sortCol) sortAsc = !sortAsc;
    else { sortCol = col; sortAsc = col === 'p_value' || col === 'p_fdr' || col === 'p_bonferroni'; }
    buildTable();
  });
});

// Filter event listeners
['filter-source', 'filter-sig'].forEach(id => {
  document.getElementById(id).addEventListener('change', buildPlot);
});
document.getElementById('filter-search').addEventListener('input', buildPlot);
['tbl-source', 'tbl-sig'].forEach(id => {
  document.getElementById(id).addEventListener('change', buildTable);
});
document.getElementById('tbl-search').addEventListener('input', buildTable);

// Initial render
buildPlot();
buildCondTable();
buildTable();
</script>
</body>
</html>"""


class Handler(SimpleHTTPRequestHandler):
    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path == "/" or parsed.path == "/index.html":
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.end_headers()
            html = HTML_TEMPLATE.replace("%%DATA_JSON%%", data_json)
            html = html.replace("%%N_TYPES%%", str(len(app_data["points"])))
            html = html.replace("%%N_COND%%", str(len(app_data["conditional"])))
            self.wfile.write(html.encode("utf-8"))
        elif parsed.path.startswith("/drivers/"):
            # Gene driver scatter page for a specific cell type
            from urllib.parse import unquote
            cell_type = unquote(parsed.path.replace("/drivers/", "", 1)).rstrip("/")
            if cell_type in driver_provider:
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.end_headers()
                ct_json = driver_provider.get_json(cell_type)
                html = DRIVER_HTML_TEMPLATE.replace("%%DRIVER_JSON%%", ct_json)
                html = html.replace("%%CELL_TYPE%%", cell_type)
                self.wfile.write(html.encode("utf-8"))
            else:
                self.send_error(404, f"No driver data for cell type: {cell_type}")
        else:
            self.send_error(404)

    def log_message(self, format, *args):
        # Quieter logging
        pass


if __name__ == "__main__":
    print("Loading data...")
    app_data = prepare_data()
    data_json = json.dumps(app_data)
    print(f"  Loaded {len(app_data['points'])} cell types, {len(app_data['conditional'])} conditional results")
    print(f"  Bonferroni threshold: -log10(p) = {app_data['bonf_threshold']:.2f}")
    print(f"  FDR threshold: -log10(p) = {app_data['fdr_threshold']:.2f}")

    print("Loading gene driver data...")
    driver_provider = DriverDataProvider()
    if driver_provider.available:
        print(f"  Ready: {len(driver_provider.cell_types)} SEA-AD cell types (computed on demand)")
    else:
        print("  WARNING: Gene driver data not available (missing specificity or GWAS files)")

    server = HTTPServer(("0.0.0.0", PORT), Handler)
    print(f"\nServer running at http://localhost:{PORT}")
    print(f"  Gene driver pages at http://localhost:{PORT}/drivers/{{cell_type}}")
    print("Press Ctrl+C to stop")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        server.server_close()
