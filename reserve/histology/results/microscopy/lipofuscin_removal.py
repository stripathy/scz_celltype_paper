"""
Lipofuscin removal and representative image generation for RNAscope composites.

Pipeline:
    1. Load raw RGB composite TIFF (1024×1024) from Arbabi et al. dataset
       - Blue = DAPI (405 nm)
       - Green = SST (488 nm) for R-section  [VIP for 568 swap - see channel map]
       - Red = VIP (568 nm)
    2. Crop to 800×800 counting frame (per Dwight's stereology protocol)
    3. Per-channel percentile normalization
    4. Spectral subtraction to suppress lipofuscin autofluorescence
    5. Optional: overlay manual cell annotations for figure generation

Usage (as script):
    python lipofuscin_removal.py
        --generates the main figures (fig_L23_manual_stars_on_raw.* and
         fig_all_manual_stars_on_raw.*) from the SVG with manual annotations.

Usage (as module):
    from lipofuscin_removal import load_clean, plain_original
    composite, red_clean, green_clean, blue = load_clean('1047-R-5.tif')

See README.md in this directory for the rationale behind the subtraction approach.
"""

import os
import re
import io
import base64
import json
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from scipy import ndimage

# Repo paths (relative to results/microscopy/)
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
RAW_IMAGE_DIR = os.path.join(REPO, 'data', 'dwight microscopy data representative images')

# ── Lipofuscin removal parameters ────────────────────────────────────────
# Rationale:
#   Lipofuscin autofluorescence is broadband — it fluoresces in BOTH the
#   488 nm (SST) and 568 nm (VIP) channels simultaneously. True RNAscope
#   signal is channel-specific (SST puncta appear only in green; VIP puncta
#   only in red). We therefore estimate the lipofuscin component as the
#   pixelwise minimum of the two signal channels and subtract it from each.
#   See README.md for the full rationale and literature precedent.
LIPO_FACTOR     = 1.2   # multiply lipo estimate by this (>1 = slight over-subtract,
                        # ensures residual lipofuscin is removed even when its
                        # red/green intensities aren't perfectly equal)
CROSS_CH_FACTOR = 0.15  # additionally subtract 15% of opposite channel (further
                        # suppresses residual crosstalk/bleedthrough)
THRESHOLD       = 0.03  # zero out dim residual pixels below this intensity

# Image acquisition parameters (from protocol)
IMAGE_SIZE_PX   = 1024                 # original TIFF dimensions
COUNTING_FRAME  = (112, 912)           # slice [112:912] → 800×800 counting region
FOV_UM          = 333                  # µm per side for 800 px counting frame (actually 260 µm;
                                       # 333 µm = full 1024 FOV; counting frame is 800/1024 * 333)
COUNTING_FRAME_UM = 800 / 1024 * 333   # ≈ 260 µm

# Channel assignments (R-section)
CHANNEL_DAPI = 2   # blue
CHANNEL_SST  = 1   # green  (488 nm)
CHANNEL_VIP  = 0   # red    (568 nm)


# ─────────────────────────────────────────────────────────────────────────
# Core image processing
# ─────────────────────────────────────────────────────────────────────────

def plain_original(fname, image_dir=RAW_IMAGE_DIR):
    """Load raw composite with per-channel percentile normalization.

    No lipofuscin removal is applied — this is what the raw tissue looks
    like on a microscope screen, with all autofluorescence present.

    Parameters
    ----------
    fname : str
        Filename under image_dir (e.g., '1047-R-5.tif').

    Returns
    -------
    np.ndarray
        800×800×3 float32 array in [0, 1].
    """
    img = np.array(Image.open(os.path.join(image_dir, fname))).astype(np.float32)
    img = img[COUNTING_FRAME[0]:COUNTING_FRAME[1],
              COUNTING_FRAME[0]:COUNTING_FRAME[1]]
    norm = np.zeros_like(img)
    for c in range(3):
        p1 = np.percentile(img[:, :, c], 1)
        p99 = np.percentile(img[:, :, c], 99.5)
        norm[:, :, c] = np.clip((img[:, :, c] - p1) / (p99 - p1 + 1e-6), 0, 1)
    return norm


def load_clean(fname, image_dir=RAW_IMAGE_DIR,
               lipo_factor=LIPO_FACTOR, cross_ch_factor=CROSS_CH_FACTOR,
               threshold=THRESHOLD, display_gain=2.5):
    """Load composite with lipofuscin suppression.

    Returns both the display-ready RGB composite and the cleaned individual
    channels (for cell detection / further analysis).

    Spectral subtraction:
        lipo        = min(red, green)
        red_clean   = red   - lipo_factor * lipo - cross_ch_factor * green
        green_clean = green - lipo_factor * lipo - cross_ch_factor * red

    Parameters
    ----------
    fname : str
        Filename under image_dir.
    lipo_factor, cross_ch_factor, threshold : float
        Lipofuscin suppression parameters. Defaults from module-level constants.
    display_gain : float
        Multiplier applied to cleaned red/green channels for display (does NOT
        affect the returned red_c, green_c arrays).

    Returns
    -------
    composite : np.ndarray
        800×800×3 display-ready RGB image with lipofuscin suppressed.
    red_c, green_c, blue : np.ndarray
        Per-channel cleaned arrays (800×800, float32 in [0, 1]).
    """
    img = np.array(Image.open(os.path.join(image_dir, fname))).astype(np.float32)
    img = img[COUNTING_FRAME[0]:COUNTING_FRAME[1],
              COUNTING_FRAME[0]:COUNTING_FRAME[1]]

    # DAPI: percentile-stretch (nuclei should always fill the range)
    p1_b = np.percentile(img[:, :, CHANNEL_DAPI], 1)
    p99_b = np.percentile(img[:, :, CHANNEL_DAPI], 99.5)
    blue = np.clip((img[:, :, CHANNEL_DAPI] - p1_b) / (p99_b - p1_b + 1e-6), 0, 1)

    # SST/VIP: fixed normalization (divide by 100) for comparability across images
    # and to avoid amplifying noise in dim sections
    red = np.clip(img[:, :, CHANNEL_VIP] / 100.0, 0, 1)
    green = np.clip(img[:, :, CHANNEL_SST] / 100.0, 0, 1)

    # Lipofuscin estimate = pixels bright in BOTH SST and VIP channels
    lipo = np.minimum(red, green)

    # Subtract lipofuscin and cross-channel contamination
    red_c = np.clip(red - lipo_factor * lipo - cross_ch_factor * green, 0, 1)
    green_c = np.clip(green - lipo_factor * lipo - cross_ch_factor * red, 0, 1)

    # Floor: zero out dim residual pixels
    red_c = np.where(red_c < threshold, 0, red_c)
    green_c = np.where(green_c < threshold, 0, green_c)

    # Build display composite (with gain applied so cleaned signal is visible)
    red_disp = np.clip(red_c * display_gain, 0, 1)
    green_disp = np.clip(green_c * display_gain, 0, 1)
    composite = np.stack([red_disp, green_disp, blue], axis=-1)
    return composite, red_c, green_c, blue


# ─────────────────────────────────────────────────────────────────────────
# Reading manual annotations from Inkscape-edited SVGs
# ─────────────────────────────────────────────────────────────────────────

def detect_stars_from_rendered_svg(svg_path, panel_svg_bounds, output_width=1600,
                                    dedup_distance_px=40):
    """Detect star marker positions by rendering SVG to PNG and finding colored blobs.

    This empirical approach avoids any coordinate-math errors in converting
    SVG transforms back to image pixel coordinates. The rendered SVG's
    visual layout is taken as ground truth.

    Parameters
    ----------
    svg_path : str
        Path to Inkscape-edited SVG with star annotations.
    panel_svg_bounds : dict of {panel_id: (x1, y1, x2, y2)}
        Each panel's bounding box in SVG point coordinates (after all transforms).
        Markers falling within a panel's bbox are assigned to that panel.
    output_width : int
        Width to render SVG (height auto-scaled). Higher = better spatial accuracy.
    dedup_distance_px : int
        Merge markers of the same type within same panel within this many pixels
        (the hollow star shape tends to be detected as 2–3 fragments).

    Returns
    -------
    list of dicts with keys: type ('sst' or 'vip'), panel (int), pixel (x, y)
        Pixel coords are in the 800×800 counting-frame reference.
    """
    try:
        import cairosvg
    except ImportError:
        raise ImportError('cairosvg is required to render SVG. Install via: pip install cairosvg')

    png_bytes = cairosvg.svg2png(url=svg_path, output_width=output_width)
    rendered = np.array(Image.open(io.BytesIO(png_bytes)).convert('RGB'))

    r, g, b = rendered[:, :, 0], rendered[:, :, 1], rendered[:, :, 2]
    # Green stars (#00FF66) — target colour used in the annotation script
    green_mask = (g > 200) & (r < 100) & (b > 50) & (b < 180)
    # Pink/Red stars (#FF3355)
    red_mask = (r > 200) & (g < 120) & (b < 150)

    # Morphological cleanup
    green_mask = ndimage.binary_opening(green_mask, iterations=1)
    red_mask = ndimage.binary_opening(red_mask, iterations=1)

    green_labels, n_green = ndimage.label(green_mask)
    red_labels, n_red = ndimage.label(red_mask)

    def _extract(labels, n, min_size=50, max_size=500):
        pts = []
        for i in range(1, n + 1):
            ys, xs = np.where(labels == i)
            if min_size <= len(ys) <= max_size:
                pts.append((ys.mean(), xs.mean(), len(ys)))
        return pts

    green_stars = _extract(green_labels, n_green)
    red_stars = _extract(red_labels, n_red)

    # SVG→PNG scale (from SVG width=917.4 pt standard output → PNG output_width)
    # We'll look up the actual SVG width from the file
    with open(svg_path) as f:
        header = f.read(4096)
    svg_width_match = re.search(r'width="([\d.]+)(?:pt)?"', header)
    svg_width_pt = float(svg_width_match.group(1)) if svg_width_match else 917.4
    svg_to_png = output_width / svg_width_pt

    # Convert panel SVG bounds to PNG pixel bounds
    panel_png_bounds = {
        pid: (x1 * svg_to_png, y1 * svg_to_png, x2 * svg_to_png, y2 * svg_to_png)
        for pid, (x1, y1, x2, y2) in panel_svg_bounds.items()
    }

    def _classify(sy, sx):
        """Return (panel_id, 800-px-x, 800-px-y) or None."""
        for pid, (x1, y1, x2, y2) in panel_png_bounds.items():
            if x1 <= sx <= x2 and y1 <= sy <= y2:
                rel_x = (sx - x1) / (x2 - x1)
                rel_y = (sy - y1) / (y2 - y1)
                return pid, rel_x * 800, rel_y * 800
        return None

    markers = []
    for (sy, sx, size) in green_stars:
        r = _classify(sy, sx)
        if r: markers.append({'type': 'sst', 'panel': r[0], 'pixel': (r[1], r[2]), 'size': size})
    for (sy, sx, size) in red_stars:
        r = _classify(sy, sx)
        if r: markers.append({'type': 'vip', 'panel': r[0], 'pixel': (r[1], r[2]), 'size': size})

    # Deduplicate: cluster same-type markers within same panel within distance threshold
    deduped = []
    used = set()
    for i, m in enumerate(markers):
        if i in used:
            continue
        cluster = [m]
        used.add(i)
        for j in range(i + 1, len(markers)):
            if j in used:
                continue
            o = markers[j]
            if o['type'] == m['type'] and o['panel'] == m['panel']:
                dx = o['pixel'][0] - m['pixel'][0]
                dy = o['pixel'][1] - m['pixel'][1]
                if np.sqrt(dx * dx + dy * dy) < dedup_distance_px:
                    cluster.append(o)
                    used.add(j)
        ax = np.mean([c['pixel'][0] for c in cluster])
        ay = np.mean([c['pixel'][1] for c in cluster])
        deduped.append({'type': m['type'], 'panel': m['panel'], 'pixel': (ax, ay)})
    return deduped


# ─────────────────────────────────────────────────────────────────────────
# Plotting helpers
# ─────────────────────────────────────────────────────────────────────────

def _add_scale_bar(ax, bar_length_um=100, y_pos_px=760, x_start_px=40):
    """Add a 100-µm scale bar to a matplotlib image axes."""
    bar_length_px = bar_length_um / (COUNTING_FRAME_UM / 800)
    ax.plot([x_start_px, x_start_px + bar_length_px], [y_pos_px, y_pos_px],
            color='white', linewidth=4)
    ax.text(x_start_px + bar_length_px / 2, y_pos_px - 20,
            f'{bar_length_um} µm', color='white', ha='center',
            fontsize=12, fontweight='bold')


def _plot_markers(ax, markers, bbox_max=800):
    """Plot SST (green) and VIP (red) stars on axes."""
    for m in markers:
        px, py = m['pixel']
        if not (0 <= px <= bbox_max and 0 <= py <= bbox_max):
            continue
        if m['type'] == 'sst':
            ax.plot(px, py, marker='*', markersize=22, color='#00FF66',
                    markeredgecolor='white', markeredgewidth=1.5,
                    linestyle='none', alpha=0.95, zorder=5)
        else:
            ax.plot(px, py, marker='*', markersize=22, color='#FF3355',
                    markeredgecolor='white', markeredgewidth=1.5,
                    linestyle='none', alpha=0.95, zorder=5)


# ─────────────────────────────────────────────────────────────────────────
# Main figure generation
# ─────────────────────────────────────────────────────────────────────────

# The 4 representative panels used in the manuscript
# (subject, diagnosis, layer, site, filename, panel_id for SVG mapping)
PANELS = [
    ('1047-R-5.tif',  1047, 'Control', 'L2/3', 5,  0),
    ('1211-R-6.tif',  1211, 'SCHIZ',   'L2/3', 6,  1),
    ('1047-R-12.tif', 1047, 'Control', 'L5/6', 12, 2),
    ('1211-R-17.tif', 1211, 'SCHIZ',   'L5/6', 17, 3),
]

# SVG panel bounding boxes (after matplotlib's transforms, in SVG pt)
# These match the layout matplotlib produces with figsize=(14, 14), 2×2 grid,
# tight_layout. Used to classify stars detected in the rendered SVG.
PANEL_SVG_BOUNDS = {
    0: (9.54, 23.17, 415.62, 429.25),      # top-left (Control L2/3)
    1: (506.99, 23.17, 913.07, 429.25),    # top-right (SCHIZ L2/3)
    2: (9.54, 470.59, 415.62, 876.67),     # bottom-left (Control L5/6)
    3: (506.99, 470.59, 913.07, 876.67),   # bottom-right (SCHIZ L5/6)
}


def generate_original_vs_cleaned_figure(output_path='fig_original_vs_cleaned_L23.png'):
    """L2/3 comparison figure: original (top row) vs lipofuscin-cleaned (bottom row).

    Illustrates the effect of the lipofuscin removal strategy. Cells are NOT
    annotated here — this is purely about showing what the cleaning does.
    """
    l23_panels = [('1047-R-5.tif', 1047, 5, 'Control', 'L2/3'),
                  ('1211-R-6.tif', 1211, 6, 'SCHIZ',   'L2/3')]

    fig, axes = plt.subplots(2, 2, figsize=(14, 14))
    for col, (fname, subj, site, diag, layer) in enumerate(l23_panels):
        diag_color = '#4C72B0' if diag == 'Control' else '#C44E52'
        # Top: original
        ax = axes[0, col]
        ax.imshow(plain_original(fname))
        ax.set_title(f'{subj} ({diag}) — {layer}  Site {site}\nOriginal',
                     fontsize=13, fontweight='bold', color=diag_color)
        ax.axis('off')
        _add_scale_bar(ax)
        # Bottom: cleaned
        ax = axes[1, col]
        composite, _, _, _ = load_clean(fname)
        ax.imshow(composite)
        ax.set_title(f'{subj} ({diag}) — {layer}  Site {site}\nLipofuscin-suppressed',
                     fontsize=13, fontweight='bold', color=diag_color)
        ax.axis('off')
        _add_scale_bar(ax)

    fig.text(0.01, 0.75, 'Original', fontsize=20, fontweight='bold',
             rotation=90, va='center', color='#333333')
    fig.text(0.01, 0.28, 'Lipofuscin\nsuppressed', fontsize=18, fontweight='bold',
             rotation=90, va='center', ha='center', color='#333333')
    fig.suptitle('L2/3 comparison: original vs lipofuscin-suppressed composites\n'
                 'Blue = DAPI  |  Green = SST (488 nm)  |  Red = VIP (568 nm)',
                 fontsize=16, fontweight='bold', y=0.97)
    fig.tight_layout(rect=[0.02, 0.02, 1, 0.96])
    fig.savefig(os.path.join(HERE, output_path), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {output_path}')


def generate_manual_annotations_figure(manual_svg,
                                        output_basename='fig_L23_manual_stars_on_raw',
                                        all_panels_basename='fig_all_manual_stars_on_raw'):
    """Read manual annotations from an Inkscape-edited SVG, overlay on raw images.

    Parameters
    ----------
    manual_svg : str
        Path to user's SVG with manually-placed/edited stars.
    output_basename : str
        Basename (no extension) for the focused 2-panel L2/3 figure.
    all_panels_basename : str
        Basename for the full 4-panel figure.
    """
    markers = detect_stars_from_rendered_svg(manual_svg, PANEL_SVG_BOUNDS)

    # Save coordinates for reference
    with open(os.path.join(HERE, 'marker_coordinates.csv'), 'w') as f:
        f.write('panel,type,pixel_x,pixel_y\n')
        for m in markers:
            f.write(f'{m["panel"]},{m["type"]},{m["pixel"][0]:.1f},{m["pixel"][1]:.1f}\n')

    # ── Focused L2/3 figure (panels 0 and 1) ──
    fig, axes = plt.subplots(1, 2, figsize=(16, 8.5))
    for ax_idx, panel_id in enumerate([0, 1]):
        fname, subj, diag, layer, site, _ = PANELS[panel_id]
        ax = axes[ax_idx]
        ax.imshow(plain_original(fname))
        panel_markers = [m for m in markers if m['panel'] == panel_id]
        _plot_markers(ax, panel_markers)
        n_sst = sum(1 for m in panel_markers if m['type'] == 'sst')
        n_vip = sum(1 for m in panel_markers if m['type'] == 'vip')
        color = '#4C72B0' if diag == 'Control' else '#C44E52'
        ax.set_title(f'{subj} ({diag}) — {layer}  Site {site}\n'
                     f'Manual annotations: {n_sst} SST, {n_vip} VIP',
                     fontsize=14, fontweight='bold', color=color)
        ax.axis('off')
        _add_scale_bar(ax)

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='*', color='w', markerfacecolor='#00FF66',
               markersize=20, markeredgecolor='white', markeredgewidth=1.5,
               label='SST+ cell', linestyle='none'),
        Line2D([0], [0], marker='*', color='w', markerfacecolor='#FF3355',
               markersize=20, markeredgecolor='white', markeredgewidth=1.5,
               label='VIP+ cell', linestyle='none'),
    ]
    fig.legend(handles=legend_elements, loc='lower center',
               bbox_to_anchor=(0.5, 0.02), fontsize=13, ncol=2,
               frameon=True, facecolor='#f0f0f0')
    fig.suptitle('L2/3 — manual cell annotations overlaid on raw composites\n'
                 'Blue = DAPI  |  Green = SST (488 nm)  |  Red = VIP (568 nm)',
                 fontsize=15, fontweight='bold', y=0.98)
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    fig.savefig(os.path.join(HERE, output_basename + '.png'), dpi=150, bbox_inches='tight')
    fig.savefig(os.path.join(HERE, output_basename + '.svg'), bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {output_basename}.png and .svg')

    # ── Full 4-panel figure ──
    fig, axes = plt.subplots(2, 2, figsize=(14, 14))
    for panel_id in range(4):
        row, col = panel_id // 2, panel_id % 2
        fname, subj, diag, layer, site, _ = PANELS[panel_id]
        ax = axes[row, col]
        ax.imshow(plain_original(fname))
        panel_markers = [m for m in markers if m['panel'] == panel_id]
        _plot_markers(ax, panel_markers)
        n_sst = sum(1 for m in panel_markers if m['type'] == 'sst')
        n_vip = sum(1 for m in panel_markers if m['type'] == 'vip')
        color = '#4C72B0' if diag == 'Control' else '#C44E52'
        ax.set_title(f'{subj} ({diag}) — {layer}  Site {site}\n'
                     f'Manual annotations: {n_sst} SST, {n_vip} VIP',
                     fontsize=12, fontweight='bold', color=color)
        ax.axis('off')
        _add_scale_bar(ax, y_pos_px=760)

    fig.legend(handles=legend_elements, loc='lower center',
               bbox_to_anchor=(0.5, -0.005), fontsize=13, ncol=2,
               frameon=True, facecolor='#f0f0f0')
    fig.suptitle('Manual cell annotations overlaid on raw composites — all panels',
                 fontsize=15, fontweight='bold', y=0.97)
    fig.tight_layout(rect=[0, 0.02, 1, 0.96])
    fig.savefig(os.path.join(HERE, all_panels_basename + '.png'), dpi=150, bbox_inches='tight')
    fig.savefig(os.path.join(HERE, all_panels_basename + '.svg'), bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {all_panels_basename}.png and .svg')


def main():
    """Regenerate all lipofuscin-related figures."""
    generate_original_vs_cleaned_figure()
    manual_svg = os.path.join(HERE, 'RNAscope_fig_representative_inkscape_SJT.svg')
    if os.path.exists(manual_svg):
        generate_manual_annotations_figure(manual_svg)
    else:
        print(f'No manual annotation SVG found at {manual_svg}; skipping annotated figures.')


if __name__ == '__main__':
    main()
