"""
sst_spatial.py — SST-focused and spatial analysis figures.

Figure 2: SST enrichment panel (4-panel)
Figure 3: SST subtypes by p-value, colored by layer
Figure 4: Upper-layer enrichment analysis
Figure 9: SST layer vs ephys relationships
Figure 10: SST depth volcano plot
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import pandas as pd

from ..config import (
    TITLE_FONTSIZE, LABEL_FONTSIZE, LEGEND_FONTSIZE, TICK_FONTSIZE,
    FIGURE_DPI, FDR_THRESHOLD,
)
from . import set_figure_style


# ============================================================================
# Figure 2: SST-focused enrichment panel
# ============================================================================

def plot_sst_enrichment_panel(enrichment_df, layer_df, gene_contributions, output_path):
    """
    4-panel SST-focused figure:
      A) All significant types with SST highlighted
      B) SST subtypes sorted by -log10(p), colored by layer
      C) Top contributing genes heatmap for SST types
      D) Layer depth vs enrichment scatter

    Parameters
    ----------
    enrichment_df : pd.DataFrame
        Enrichment results.
    layer_df : pd.DataFrame
        SST layer classification (with avg_depth, layer_group columns).
    gene_contributions : dict
        Mapping SST supertype -> contribution DataFrame.
    output_path : str or Path
    """
    set_figure_style()

    # Extract SST supertypes
    sst_enrichment = enrichment_df[
        enrichment_df["supertype"].str.startswith("Sst")
    ].copy()
    sst_enrichment["-log10p"] = -np.log10(sst_enrichment["p_value"].clip(lower=1e-20))

    # Merge with layer info
    if "supertype" in layer_df.columns:
        sst_merged = sst_enrichment.merge(
            layer_df[["supertype", "avg_depth", "layer_group"]],
            on="supertype", how="left"
        )
    else:
        sst_merged = sst_enrichment.copy()
        sst_merged["avg_depth"] = np.nan
        sst_merged["layer_group"] = "Unknown"

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Panel A: SST subtypes by -log10(p)
    ax = axes[0, 0]
    sst_sorted = sst_merged.sort_values("-log10p", ascending=True)
    layer_colors = {"Upper": "#E74C3C", "Middle": "#F39C12", "Deep": "#3498DB", "Unknown": "#95A5A6"}
    colors = [layer_colors.get(lg, "#95A5A6") for lg in sst_sorted["layer_group"]]
    ax.barh(range(len(sst_sorted)), sst_sorted["-log10p"], color=colors)
    ax.set_yticks(range(len(sst_sorted)))
    ax.set_yticklabels(sst_sorted["supertype"], fontsize=10)
    ax.set_xlabel("-log10(p-value)", fontsize=LABEL_FONTSIZE)
    ax.set_title("A) SST Subtype SCZ Enrichment", fontsize=TITLE_FONTSIZE - 2)
    # FDR line
    fdr_types = sst_enrichment[sst_enrichment["p_fdr"] < FDR_THRESHOLD]
    if len(fdr_types) > 0:
        fdr_line = -np.log10(fdr_types["p_value"].max())
        ax.axvline(fdr_line, color="red", linestyle="--", linewidth=1, label="FDR < 0.05")
    ax.legend(fontsize=10)
    sns.despine(ax=ax)

    # Panel B: Layer depth vs enrichment
    ax = axes[0, 1]
    valid = sst_merged.dropna(subset=["avg_depth"])
    for lg, color in layer_colors.items():
        mask = valid["layer_group"] == lg
        if mask.any():
            ax.scatter(valid.loc[mask, "avg_depth"], valid.loc[mask, "-log10p"],
                      c=color, s=100, label=lg, edgecolors="k", linewidth=0.5, zorder=3)
            for _, row in valid[mask].iterrows():
                ax.annotate(row["supertype"], (row["avg_depth"], row["-log10p"]),
                           fontsize=8, ha="left", va="bottom",
                           xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("Cortical Depth (0=pial, 1=WM)", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("-log10(p-value)", fontsize=LABEL_FONTSIZE)
    ax.set_title("B) Depth vs SCZ Enrichment", fontsize=TITLE_FONTSIZE - 2)
    ax.legend(fontsize=10)
    sns.despine(ax=ax)

    # Panel C: Top genes heatmap
    ax = axes[1, 0]
    # Get top 15 genes across significant SST types
    sig_sst = sst_enrichment[sst_enrichment["p_fdr"] < FDR_THRESHOLD]["supertype"].tolist()
    if gene_contributions and sig_sst:
        all_top_genes = []
        for st in sig_sst[:6]:  # Top 6 SST types
            if st in gene_contributions:
                top = gene_contributions[st].head(10)["gene"].tolist()
                all_top_genes.extend(top)
        top_genes = list(dict.fromkeys(all_top_genes))[:20]  # Unique, top 20

        if top_genes and sig_sst:
            heatmap_data = pd.DataFrame(index=top_genes, columns=sig_sst[:6])
            for st in sig_sst[:6]:
                if st in gene_contributions:
                    gc = gene_contributions[st].set_index("gene")
                    for g in top_genes:
                        if g in gc.index:
                            heatmap_data.loc[g, st] = gc.loc[g, "contribution"]
                        else:
                            heatmap_data.loc[g, st] = 0
            heatmap_data = heatmap_data.astype(float)
            sns.heatmap(heatmap_data, ax=ax, cmap="YlOrRd", annot=False,
                       xticklabels=True, yticklabels=True, cbar_kws={"shrink": 0.5})
            ax.set_title("C) Top Contributing Genes", fontsize=TITLE_FONTSIZE - 2)
        else:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
    else:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("C) Top Contributing Genes", fontsize=TITLE_FONTSIZE - 2)

    # Panel D: Summary statistics
    ax = axes[1, 1]
    n_bonf = (sst_enrichment["p_bonferroni"] < 0.05).sum()
    n_fdr = (sst_enrichment["p_fdr"] < FDR_THRESHOLD).sum()
    summary_text = (
        f"SST Enrichment Summary\n\n"
        f"Total SST supertypes: {len(sst_enrichment)}\n"
        f"Bonferroni significant: {n_bonf}\n"
        f"FDR significant: {n_fdr}\n\n"
    )
    if len(valid) > 0:
        upper = valid[valid["layer_group"] == "Upper"]
        deep = valid[valid["layer_group"] == "Deep"]
        summary_text += (
            f"Upper-layer types: {len(upper)}\n"
            f"  Mean -log10(p): {upper['-log10p'].mean():.1f}\n"
            f"Deep-layer types: {len(deep)}\n"
            f"  Mean -log10(p): {deep['-log10p'].mean():.1f}\n"
        )
    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes,
            fontsize=14, verticalalignment="top", fontfamily="monospace")
    ax.axis("off")
    ax.set_title("D) Summary", fontsize=TITLE_FONTSIZE - 2)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ============================================================================
# Figure 3: SST subtypes by p-value, colored by layer
# ============================================================================

def plot_sst_pvalues_by_layer(enrichment_df, layer_df, output_path):
    """
    Horizontal bar chart: SST supertypes sorted by -log10(p), colored by layer.

    Parameters
    ----------
    enrichment_df : pd.DataFrame
        Enrichment results.
    layer_df : pd.DataFrame
        SST layer info with supertype, avg_depth, layer_group.
    output_path : str or Path
    """
    set_figure_style()

    sst = enrichment_df[enrichment_df["supertype"].str.startswith("Sst")].copy()
    sst["-log10p"] = -np.log10(sst["p_value"].clip(lower=1e-20))

    sst = sst.merge(layer_df[["supertype", "avg_depth", "layer_group"]],
                     on="supertype", how="left")
    sst = sst.sort_values("-log10p", ascending=True)

    layer_colors = {"Upper": "#E74C3C", "Middle": "#F39C12", "Deep": "#3498DB"}
    colors = [layer_colors.get(lg, "#95A5A6") for lg in sst["layer_group"]]

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.barh(range(len(sst)), sst["-log10p"], color=colors)
    ax.set_yticks(range(len(sst)))
    ax.set_yticklabels(
        [f"{row['supertype']} (d={row['avg_depth']:.2f})"
         for _, row in sst.iterrows()],
        fontsize=11
    )

    # FDR line
    fdr_types = sst[sst["p_fdr"] < FDR_THRESHOLD]
    if len(fdr_types) > 0:
        fdr_line = -np.log10(fdr_types["p_value"].max())
        ax.axvline(fdr_line, color="black", linestyle="--", linewidth=1.5,
                   label="FDR < 0.05")

    ax.set_xlabel("-log10(p-value)", fontsize=LABEL_FONTSIZE)
    ax.set_title("SST Subtype SCZ Enrichment by Cortical Depth",
                 fontsize=TITLE_FONTSIZE)

    handles = [mpatches.Patch(color=c, label=l) for l, c in layer_colors.items()]
    ax.legend(handles=handles, loc="lower right", fontsize=LEGEND_FONTSIZE)

    sns.despine(ax=ax)
    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ============================================================================
# Figure 4: Upper-layer enrichment analysis
# ============================================================================

def plot_upper_layer_enrichment(spatial_df, output_path):
    """
    Scatter plot: upper_fraction vs GWAS Z, highlighting upper-biased SCZ risk genes.

    Parameters
    ----------
    spatial_df : pd.DataFrame
        From spatial.compute_upper_layer_gene_scores.
    output_path : str or Path
    """
    set_figure_style()

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Panel A: upper_fraction vs GWAS Z for SCZ genes
    ax = axes[0]
    if "is_scz_gwas" in spatial_df.columns:
        scz = spatial_df[spatial_df["is_scz_gwas"]].copy()
    else:
        scz = spatial_df.nlargest(100, "gwas_zstat").copy()

    scz = scz[scz["total_sst_spec"] > 0.001]  # Filter to genes with SST expression

    ax.scatter(scz["upper_fraction"], scz["gwas_zstat"],
              alpha=0.6, s=40, c="#3498DB", edgecolors="k", linewidth=0.3)

    # Label top genes
    top = scz.nlargest(12, "upper_gwas_score")
    for _, row in top.iterrows():
        ax.annotate(row["gene"],
                   (row["upper_fraction"], row["gwas_zstat"]),
                   fontsize=9, ha="left", va="bottom",
                   xytext=(3, 3), textcoords="offset points")

    ax.axvline(0.5, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("Upper-Layer Fraction of SST Expression", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("GWAS Z-statistic", fontsize=LABEL_FONTSIZE)
    ax.set_title("A) SCZ Risk Genes: Upper vs Deep SST", fontsize=TITLE_FONTSIZE - 2)
    sns.despine(ax=ax)

    # Panel B: Top genes by upper_gwas_score
    ax = axes[1]
    top_genes = spatial_df[spatial_df["total_sst_spec"] > 0.001].nlargest(20, "upper_gwas_score")
    colors = ["#E74C3C" if f > 0.5 else "#3498DB"
              for f in top_genes["upper_fraction"]]

    ax.barh(range(len(top_genes)), top_genes["upper_gwas_score"],
            color=colors)
    ax.set_yticks(range(len(top_genes)))
    labels = [
        f"{row['gene']} ({row['upper_fraction']:.0%} upper, Z={row['gwas_zstat']:.1f})"
        for _, row in top_genes.iterrows()
    ]
    ax.set_yticklabels(labels, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel("Upper SST Spec \u00d7 GWAS Z", fontsize=LABEL_FONTSIZE)
    ax.set_title("B) Top Upper-Layer SST SCZ Genes", fontsize=TITLE_FONTSIZE - 2)
    sns.despine(ax=ax)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ============================================================================
# Figure 9: SST layer vs ephys relationships
# ============================================================================

def plot_sst_layer_ephys(layer_df, output_path):
    """
    SST supertype layer position vs SAG and TAU.

    Parameters
    ----------
    layer_df : pd.DataFrame
        SST supertype info with avg_depth, mean_sag, mean_tau, layer_group.
    output_path : str or Path
    """
    set_figure_style()

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    layer_colors = {"Upper": "#E74C3C", "Middle": "#F39C12", "Deep": "#3498DB"}

    for ax, (feat, label) in zip(axes, [("mean_sag", "SAG"), ("mean_tau", "TAU (ms)")]):
        if feat not in layer_df.columns:
            ax.text(0.5, 0.5, f"No {feat} data", ha="center", va="center",
                    transform=ax.transAxes)
            continue

        for lg, color in layer_colors.items():
            mask = layer_df["layer_group"] == lg
            sub = layer_df[mask]
            if len(sub) > 0:
                ax.scatter(sub["avg_depth"], sub[feat], c=color, s=100,
                          label=lg, edgecolors="k", linewidth=0.5, zorder=3)
                for _, row in sub.iterrows():
                    ax.annotate(row["supertype"], (row["avg_depth"], row[feat]),
                               fontsize=8, ha="left", va="bottom",
                               xytext=(3, 3), textcoords="offset points")

        ax.set_xlabel("Cortical Depth (0=pial, 1=WM)", fontsize=LABEL_FONTSIZE)
        ax.set_ylabel(label, fontsize=LABEL_FONTSIZE)
        ax.set_title(f"Cortical Depth vs {label}", fontsize=TITLE_FONTSIZE - 2)
        ax.legend(fontsize=LEGEND_FONTSIZE)
        sns.despine(ax=ax)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ============================================================================
# Figure 10: SST depth volcano plot
# ============================================================================

def plot_depth_volcano(
    depth_corr_df,
    output_path,
    title="Gene-Depth Correlations Across SST Supertypes",
    n_label=15,
    highlight_genes=None,
    magma_threshold=3.0,
):
    """
    Volcano plot: gene specificity/expression correlation with cortical depth.

    X-axis: Spearman rho (negative = upper-layer enriched, positive = deep-layer)
    Y-axis: -log10(p-value)

    Points colored categorically: red = risk (signed MAGMA Z > 0),
    blue = protective (signed MAGMA Z < 0), gray = no direction data.
    Point size scales with |MAGMA Z| so high-signal genes stand out.
    SCZ GWAS loci genes get a black edge.

    Parameters
    ----------
    depth_corr_df : pd.DataFrame
        Must have: gene, spearman_rho, pval, pval_fdr.
        Optional: signed_magma_z (float), magma_zstat (float),
                  is_scz_gwas (bool), lead_snp_beta (float).
    output_path : str or Path
        Path to save figure.
    title : str
        Figure title.
    n_label : int
        Number of top genes to label on each side.
    highlight_genes : list of str, optional
        Specific genes to always label if present.
    magma_threshold : float
        |ZSTAT| threshold for labeling high-signal genes (default: 3.0).
    """
    set_figure_style()

    df = depth_corr_df.copy()
    df["-log10p"] = -np.log10(df["pval"].clip(lower=1e-300))

    has_gwas = "is_scz_gwas" in df.columns
    has_signed = "signed_magma_z" in df.columns
    has_magma = "magma_zstat" in df.columns

    # Colors
    RISK_COLOR = "#D32F2F"       # red
    PROTECTIVE_COLOR = "#1565C0" # blue
    NEUTRAL_COLOR = "#BDBDBD"    # light gray for no direction data

    fig, ax = plt.subplots(figsize=(13, 8))

    if has_signed:
        # --- Categorical risk/protective coloring ---

        # Assign category and color per gene
        df["direction"] = "neutral"
        df.loc[df["signed_magma_z"] > 0, "direction"] = "risk"
        df.loc[df["signed_magma_z"] < 0, "direction"] = "protective"
        df.loc[df["signed_magma_z"].isna(), "direction"] = "neutral"

        # Point size scales with |MAGMA Z| — small baseline, larger for high signal
        base_size = 8.0
        df["pt_size"] = base_size
        has_z = df["magma_zstat"].notna()
        # Scale: 8 for Z~0, up to ~60 for Z~8+
        df.loc[has_z, "pt_size"] = base_size + (
            df.loc[has_z, "magma_zstat"].abs().clip(upper=8) * 6.5
        )

        # Alpha: higher for genes with stronger signal
        df["pt_alpha"] = 0.15
        df.loc[has_z, "pt_alpha"] = np.clip(
            0.15 + df.loc[has_z, "magma_zstat"].abs() * 0.1, 0.15, 0.9
        )

        color_map = {
            "risk": RISK_COLOR,
            "protective": PROTECTIVE_COLOR,
            "neutral": NEUTRAL_COLOR,
        }

        # Separate GWAS from non-GWAS
        if has_gwas:
            non_gwas = df[~df["is_scz_gwas"]]
            gwas_df = df[df["is_scz_gwas"]]
        else:
            non_gwas = df
            gwas_df = df.iloc[0:0]

        # Layer 1: Non-GWAS genes, draw each direction group
        for direction in ["neutral", "protective", "risk"]:
            sub = non_gwas[non_gwas["direction"] == direction]
            if len(sub) == 0:
                continue
            ax.scatter(
                sub["spearman_rho"], sub["-log10p"],
                s=sub["pt_size"], c=color_map[direction],
                alpha=sub["pt_alpha"].values,
                edgecolors="none", zorder=2 if direction != "neutral" else 1,
                rasterized=True,
            )

        # Layer 2: SCZ GWAS loci — same categorical colors, black edge, larger
        for direction in ["neutral", "protective", "risk"]:
            sub = gwas_df[gwas_df["direction"] == direction]
            if len(sub) == 0:
                continue
            # GWAS genes get extra size boost
            gwas_size = sub["pt_size"] * 2.5
            gwas_size = gwas_size.clip(lower=40)
            ax.scatter(
                sub["spearman_rho"], sub["-log10p"],
                s=gwas_size, c=color_map[direction],
                alpha=0.85, edgecolors="black", linewidth=0.7,
                zorder=4,
            )

        # Legend
        legend_handles = [
            mpatches.Patch(color=RISK_COLOR, label="Risk (upregulation increases SCZ risk)"),
            mpatches.Patch(color=PROTECTIVE_COLOR, label="Protective (upregulation decreases SCZ risk)"),
            mpatches.Patch(color=NEUTRAL_COLOR, label="No direction data"),
        ]
        if has_gwas and len(gwas_df) > 0:
            legend_handles.append(
                plt.scatter([], [], s=60, c="white", edgecolors="black",
                            linewidth=1, label=f"SCZ GWAS loci ({len(gwas_df)})")
            )
        ax.legend(handles=legend_handles, loc="upper center",
                  fontsize=LEGEND_FONTSIZE - 1, framealpha=0.85, ncol=2)

    elif has_magma:
        # Fallback: color by unsigned MAGMA Z
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable

        cmap = plt.cm.RdYlBu_r
        valid = df["magma_zstat"].dropna()
        vmin = max(valid.quantile(0.02), -2)
        vmax = min(valid.quantile(0.98), 10)
        norm = Normalize(vmin=vmin, vmax=vmax, clip=True)

        no_z = df[df["magma_zstat"].isna()]
        ax.scatter(no_z["spearman_rho"], no_z["-log10p"],
                   s=5, c=NEUTRAL_COLOR, alpha=0.3, zorder=1, rasterized=True)

        with_z = df[df["magma_zstat"].notna()]
        ax.scatter(
            with_z["spearman_rho"], with_z["-log10p"],
            s=12, c=with_z["magma_zstat"], cmap=cmap, norm=norm,
            alpha=0.5, edgecolors="none", zorder=2, rasterized=True,
        )
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.02, aspect=25)
        cbar.set_label("MAGMA Z-stat", fontsize=LABEL_FONTSIZE - 2)
    else:
        ax.scatter(df["spearman_rho"], df["-log10p"],
                   s=8, c="#888888", alpha=0.4, zorder=1, rasterized=True)

    # FDR threshold line
    fdr_genes = df[df["pval_fdr"] < 0.05]
    if len(fdr_genes) > 0:
        fdr_line = -np.log10(fdr_genes["pval"].max())
        ax.axhline(fdr_line, color="gray", linestyle="--", linewidth=1,
                   alpha=0.6, zorder=0)
        ax.text(ax.get_xlim()[1] * 0.92, fdr_line + 0.15, "FDR = 0.05",
                ha="right", fontsize=10, color="gray")

    # Vertical center line
    ax.axvline(0, color="gray", linestyle="-", linewidth=0.5, alpha=0.3, zorder=0)

    # --- Labeling ---
    genes_to_label = set()

    # Top upper-layer genes by significance
    upper = df[df["spearman_rho"] < 0].nlargest(n_label, "-log10p")
    genes_to_label.update(upper["gene"].tolist())

    # Top deep-layer genes by significance
    deep = df[df["spearman_rho"] > 0].nlargest(n_label, "-log10p")
    genes_to_label.update(deep["gene"].tolist())

    # Always label highlighted genes
    if highlight_genes:
        for g in highlight_genes:
            if g in df["gene"].values:
                genes_to_label.add(g)

    # Label SCZ GWAS genes that are FDR-significant for depth
    if has_gwas:
        gwas_sig = df[(df["is_scz_gwas"]) & (df["pval_fdr"] < 0.05)]
        genes_to_label.update(gwas_sig["gene"].tolist())

    # Label high-signal genes that are FDR-significant for depth
    z_col = "signed_magma_z" if has_signed else ("magma_zstat" if has_magma else None)
    if z_col:
        magma_sig = df[(df[z_col].abs() > magma_threshold) & (df["pval_fdr"] < 0.05)]
        magma_sig_top = magma_sig.nlargest(10, "-log10p")
        genes_to_label.update(magma_sig_top["gene"].tolist())

    # Draw labels
    label_df = df[df["gene"].isin(genes_to_label)]

    def _label_color(row):
        """Color label text by risk/protective direction."""
        if has_signed and pd.notna(row.get("signed_magma_z")):
            z = row["signed_magma_z"]
            if z > 3:
                return "#B71C1C"    # dark red — strong risk
            elif z > 0:
                return "#E53935"    # lighter red — risk
            elif z < -3:
                return "#0D47A1"    # dark blue — strong protective
            elif z < 0:
                return "#1E88E5"    # lighter blue — protective
        elif has_gwas and row.get("is_scz_gwas", False):
            return "#E65100"
        return "#455A64"

    try:
        from adjustText import adjust_text
        texts = []
        for _, row in label_df.iterrows():
            color = _label_color(row)
            t = ax.text(row["spearman_rho"], row["-log10p"], row["gene"],
                       fontsize=8, color=color, fontweight="bold",
                       ha="center", va="center")
            texts.append(t)
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.5, alpha=0.4),
                    expand=(1.3, 1.5))
    except ImportError:
        for _, row in label_df.iterrows():
            color = _label_color(row)
            ax.annotate(row["gene"],
                       (row["spearman_rho"], row["-log10p"]),
                       fontsize=7, color=color, fontweight="bold",
                       ha="left", va="bottom",
                       xytext=(3, 3), textcoords="offset points")

    # Direction annotations
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    y_arrow = ylim[1] * 0.97
    ax.annotate("Upper Layer (pial) \u2190",
                xy=(xlim[0] * 0.85, y_arrow),
                fontsize=12, fontweight="bold", color="#8E44AD",
                ha="left", va="top")
    ax.annotate("\u2192 Deep Layer (WM)",
                xy=(xlim[1] * 0.85, y_arrow),
                fontsize=12, fontweight="bold", color="#2980B9",
                ha="right", va="top")

    ax.set_xlabel("Spearman rho (gene vs cortical depth across SST supertypes)",
                  fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("-log10(p-value)", fontsize=LABEL_FONTSIZE)
    ax.set_title(title, fontsize=TITLE_FONTSIZE)

    sns.despine()
    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")
