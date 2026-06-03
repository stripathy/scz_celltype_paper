"""Figure generation for cell density analysis."""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from . import config


def plot_diagnosis_barplots(df_R, df_L, stats_df, save_path, title=None):
    """Create the 4×2 grid figure (cell types × layers) with diagnosis on x-axis.

    Parameters
    ----------
    df_R : DataFrame
        R-section data with SST, VIP, subject, layer, diagnosis columns.
    df_L : DataFrame
        L-section data with PV, PYR, subject, layer, diagnosis columns.
    stats_df : DataFrame
        Model statistics with columns: cell_type, term, coef, p.
    save_path : str or Path
        Where to save the figure.
    title : str, optional
        Overall figure title.
    """
    # Cell types in row order: SST, PV, VIP, PYR
    cell_configs = [
        ("SST", df_R, "SST", "SST"),
        ("PV", df_L, "PV", "PV (PVALB)"),
        ("VIP", df_R, "VIP", "VIP"),
        ("PYR", df_L, "PYR", "PYR (SLC17A7)"),
    ]

    fig, axes = plt.subplots(4, 2, figsize=(14, 20))

    for row_idx, (cell_type, df, col, display_name) in enumerate(cell_configs):
        for col_idx, (layer, layer_label) in enumerate(
            [("L23", "L2/3 (Upper)"), ("L56", "L5/6 (Deep)")]
        ):
            ax = axes[row_idx, col_idx]
            layer_df = df[df["layer"] == layer].copy()

            # Compute per-subject means
            subj_means = (
                layer_df.groupby(["subject", "diagnosis"])[col]
                .mean()
                .reset_index()
            )

            # Plot bars
            positions = []
            for i, diag in enumerate(config.DIAG_ORDER):
                diag_data = subj_means[subj_means["diagnosis"] == diag][col]
                n = len(diag_data)
                color = config.DIAG_COLORS[diag]

                bar = ax.bar(
                    i, diag_data.mean(), width=0.7,
                    color=color, alpha=0.7, edgecolor="black", linewidth=0.5,
                )

                # SEM error bar
                if n > 1:
                    sem = diag_data.std() / np.sqrt(n)
                    ax.errorbar(i, diag_data.mean(), yerr=sem,
                                fmt="none", color="black", capsize=4, linewidth=1.5)

                # Jitter individual points
                jitter = np.random.default_rng(42).uniform(-0.2, 0.2, n)
                ax.scatter(
                    np.full(n, i) + jitter, diag_data.values,
                    color=color, edgecolor="black", linewidth=0.5,
                    s=25, zorder=5, alpha=0.8,
                )

                positions.append(i)

            # Labels
            n_per_diag = [
                len(subj_means[subj_means["diagnosis"] == d]["subject"].unique())
                for d in config.DIAG_ORDER
            ]
            tick_labels = [
                f"{d}\n(n={n})" for d, n in zip(config.DIAG_ORDER, n_per_diag)
            ]
            ax.set_xticks(positions)
            ax.set_xticklabels(tick_labels, fontsize=9)
            ax.set_ylabel("Cells per site", fontsize=11)
            ax.set_title(f"{display_name} — {layer_label}", fontsize=13, fontweight="bold")

            # Add stats annotations
            _annotate_stats(ax, cell_type, layer, stats_df)

    # Overall title
    n_R = df_R["subject"].nunique()
    n_L = df_L["subject"].nunique()
    if title is None:
        title = (
            f"Cell Density by Diagnosis (cells/site)\n"
            f"Mixed model: diagnosis × layer + age + sex + (1|subject)\n"
            f"sgACC, Pitt Tetrad (SST/VIP N={n_R}, PV/PYR N={n_L})"
        )
    fig.suptitle(title, fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {save_path}")


def _annotate_stats(ax, cell_type, layer, stats_df):
    """Add statistics text box to a subplot."""
    ct_stats = stats_df[stats_df["cell_type"] == cell_type]
    if ct_stats.empty:
        return

    annotations = []

    # For L23 panels: main effects are the L23 effects
    # For L56 panels: need to note the interaction terms
    for _, row in ct_stats.iterrows():
        term = row["term"]
        p = row["p"]

        if p >= 0.05:
            continue

        # Determine which panel this annotation belongs to
        is_interaction = ":" in term and "L56" in term
        is_main_diag = term in ("BP", "MDD", "SCHIZ") and ":" not in term

        if layer == "L23":
            if is_main_diag:
                annotations.append(f"{term}: p={p:.3f}")
            elif is_interaction:
                # Show interaction on both panels
                annotations.append(f"{term.replace(':L56', '')}×layer: p={p:.3f}")
        elif layer == "L56":
            if is_main_diag:
                annotations.append(f"{term}: p={p:.3f}")
            elif is_interaction:
                annotations.append(f"{term.replace(':L56', '')}×layer: p={p:.3f}")

    if annotations:
        text = "\n".join(annotations)
        ax.text(
            0.02, 0.98, text, transform=ax.transAxes,
            fontsize=8, verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow",
                      edgecolor="gray", alpha=0.9),
        )


def plot_validation_scatter(computed_df, dwight_df, save_path):
    """Plot our computed per-subject totals vs Dwight's summary.

    Parameters
    ----------
    computed_df : DataFrame
        Our computed totals with columns: subject, PV, PYR_23, PYR_56, SST, VIP
    dwight_df : DataFrame
        Dwight's summary with same columns.
    save_path : str or Path
    """
    merged = computed_df.merge(dwight_df, on="subject", suffixes=("_ours", "_dwight"))

    cell_types = ["SST", "VIP", "PV", "PYR_23", "PYR_56"]
    fig, axes = plt.subplots(1, 5, figsize=(25, 5))

    for ax, ct in zip(axes, cell_types):
        x = merged[f"{ct}_dwight"]
        y = merged[f"{ct}_ours"]
        mask = x.notna() & y.notna()
        ax.scatter(x[mask], y[mask], alpha=0.7, edgecolor="black", linewidth=0.5)

        # Identity line
        lims = [min(x[mask].min(), y[mask].min()), max(x[mask].max(), y[mask].max())]
        ax.plot(lims, lims, "k--", alpha=0.5)

        # Correlation
        r = np.corrcoef(x[mask], y[mask])[0, 1]
        ax.set_title(f"{ct} (r={r:.3f}, n={mask.sum()})", fontsize=14)
        ax.set_xlabel("Dwight's values", fontsize=12)
        ax.set_ylabel("Our computed values", fontsize=12)

    fig.suptitle("Validation: Our Computed Values vs Dwight's Summary",
                 fontsize=16, fontweight="bold")
    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {save_path}")


def plot_area_by_diagnosis(df_R, df_L, save_path):
    """Plot section area by diagnosis group."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, (label, df) in zip(axes, [("R-section", df_R), ("L-section", df_L)]):
        # Get area per subject (from first row with non-NaN area)
        areas = df.dropna(subset=["area_um2"]).groupby(["subject", "diagnosis"])["area_um2"].first().reset_index()

        for i, diag in enumerate(config.DIAG_ORDER):
            vals = areas[areas["diagnosis"] == diag]["area_um2"] / 1e6  # convert to mm²
            n = len(vals)
            color = config.DIAG_COLORS[diag]
            ax.bar(i, vals.mean(), color=color, alpha=0.7, edgecolor="black", linewidth=0.5)
            if n > 1:
                sem = vals.std() / np.sqrt(n)
                ax.errorbar(i, vals.mean(), yerr=sem, fmt="none", color="black", capsize=4)

        ax.set_xticks(range(4))
        ax.set_xticklabels(config.DIAG_ORDER, fontsize=11)
        ax.set_ylabel("Section area (mm²)", fontsize=12)
        ax.set_title(f"{label} area by diagnosis", fontsize=13, fontweight="bold")

    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {save_path}")
