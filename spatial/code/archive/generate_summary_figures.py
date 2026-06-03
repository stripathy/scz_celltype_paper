import anndata as ad
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from collections import Counter
from multiprocessing import Pool
import glob
import time
import os

def build_raster(x, y, colors, bin_size=20, bg=(0.07, 0.07, 0.15)):
    x_edges = np.arange(x.min()-bin_size, x.max()+2*bin_size, bin_size)
    y_edges = np.arange(y.min()-bin_size, y.max()+2*bin_size, bin_size)
    nx, ny = len(x_edges)-1, len(y_edges)-1
    img_r = np.zeros((ny,nx)); img_g = np.zeros((ny,nx)); img_b = np.zeros((ny,nx))
    counts = np.zeros((ny,nx))
    xi = np.clip(np.digitize(x, x_edges)-1, 0, nx-1)
    yi = np.clip(np.digitize(y, y_edges)-1, 0, ny-1)
    for i in range(len(x)):
        img_r[yi[i],xi[i]] += colors[i][0]
        img_g[yi[i],xi[i]] += colors[i][1]
        img_b[yi[i],xi[i]] += colors[i][2]
        counts[yi[i],xi[i]] += 1
    mask = counts > 0
    img_r[mask]/=counts[mask]; img_g[mask]/=counts[mask]; img_b[mask]/=counts[mask]
    img_r[~mask]=bg[0]; img_g[~mask]=bg[1]; img_b[~mask]=bg[2]
    return np.stack([img_r,img_g,img_b], axis=-1), [x_edges[0],x_edges[-1],y_edges[-1],y_edges[0]]

def make_color_map(labels, palette_arrays):
    unique = sorted(set(labels))
    cmap = {}
    for i, u in enumerate(unique):
        cmap[u] = tuple(palette_arrays[i % len(palette_arrays)][:3])
    return cmap, unique

def process_sample(h5ad_path):
    t0 = time.time()
    sample_id = os.path.basename(h5ad_path).replace('_annotated.h5ad', '')
    
    try:
        adata = ad.read_h5ad(h5ad_path)
        coords = adata.obsm['spatial']
        x, y = coords[:, 0], coords[:, 1]
        n_cells = adata.shape[0]
        
        subclass = adata.obs['subclass_label'].values.astype(str)
        cluster = adata.obs['cluster_label'].values.astype(str)
        class_label = adata.obs['class_label'].values.astype(str)
        pred_depth = adata.obs['predicted_norm_depth'].values.astype(float)
        
        # Color palettes
        palette_main = np.vstack([
            plt.cm.tab20(np.linspace(0, 1, 20)),
            plt.cm.Set3(np.linspace(0, 1, 12))
        ])
        palette_cluster = np.vstack([
            plt.cm.tab20(np.linspace(0, 1, 20)),
            plt.cm.tab20b(np.linspace(0, 1, 20)),
            plt.cm.tab20c(np.linspace(0, 1, 20)),
            plt.cm.Set1(np.linspace(0, 1, 9)),
            plt.cm.Set2(np.linspace(0, 1, 8)),
            plt.cm.Set3(np.linspace(0, 1, 12)),
        ])
        
        # Class colors (matching actual class_label values in data)
        class_cmap = {
            'Neuronal: Glutamatergic': (0.2, 0.6, 0.9),
            'Neuronal: GABAergic': (0.9, 0.3, 0.2),
            'Non-neuronal and Non-neural': (0.3, 0.8, 0.3),
        }
        # Short display names for legend
        class_display = {
            'Neuronal: Glutamatergic': 'Glutamatergic',
            'Neuronal: GABAergic': 'GABAergic',
            'Non-neuronal and Non-neural': 'Non-neuronal',
        }
        
        # Subclass color map
        sub_cmap, unique_sub = make_color_map(subclass, palette_main)
        # Cluster color map
        clust_cmap, unique_clust = make_color_map(cluster, palette_cluster)
        
        # Layer definitions
        layer_bins = {
            'L1': (-np.inf, 0.1), 'L2/3': (0.1, 0.3), 'L4': (0.3, 0.45),
            'L5': (0.45, 0.65), 'L6': (0.65, 0.85), 'WM': (0.85, np.inf),
        }
        layer_colors = {
            'L1': (0.9, 0.3, 0.3), 'L2/3': (0.3, 0.8, 0.3), 'L4': (0.3, 0.3, 0.9),
            'L5': (0.9, 0.6, 0.1), 'L6': (0.7, 0.3, 0.8), 'WM': (0.5, 0.5, 0.5),
        }
        
        # Compute shared extent for all spatial panels
        x_pad = (x.max() - x.min()) * 0.02
        y_pad = (y.max() - y.min()) * 0.02
        shared_xlim = (x.min() - x_pad, x.max() + x_pad)
        shared_ylim = (y.max() + y_pad, y.min() - y_pad)  # inverted for image coords
        
        # ============ BUILD FIGURE ============
        fig, axes = plt.subplots(2, 3, figsize=(54, 36), facecolor='white',
                                 constrained_layout=False)
        
        # ---------- Row 1: Cell Type Annotations ----------
        
        # Panel (0,0): Subclass
        colors_sub = np.array([sub_cmap.get(s, (0.5,0.5,0.5)) for s in subclass])
        img, ext = build_raster(x, y, colors_sub)
        axes[0,0].imshow(img, extent=ext, aspect='equal', interpolation='nearest', origin='upper')
        axes[0,0].set_xlim(shared_xlim); axes[0,0].set_ylim(shared_ylim)
        counts_sub = Counter(subclass)
        top_sub = sorted(counts_sub.keys(), key=lambda k: -counts_sub[k])[:15]
        handles = [Line2D([0],[0], marker='o', color='w', markerfacecolor=sub_cmap.get(s,(0.5,0.5,0.5)),
                          markersize=14, label=f'{s} ({counts_sub[s]:,})') for s in top_sub]
        axes[0,0].legend(handles=handles, loc='upper right', fontsize=14, framealpha=0.9, markerscale=1.2)
        axes[0,0].set_title('Subclass Labels', fontsize=20, fontweight='bold')
        axes[0,0].set_xlabel('x (µm)', fontsize=16); axes[0,0].set_ylabel('y (µm)', fontsize=16)
        axes[0,0].tick_params(labelsize=13)
        
        # Panel (0,1): Cluster
        colors_clust = np.array([clust_cmap.get(c, (0.5,0.5,0.5)) for c in cluster])
        img, ext = build_raster(x, y, colors_clust)
        axes[0,1].imshow(img, extent=ext, aspect='equal', interpolation='nearest', origin='upper')
        axes[0,1].set_xlim(shared_xlim); axes[0,1].set_ylim(shared_ylim)
        counts_clust = Counter(cluster)
        top_clust = sorted(counts_clust.keys(), key=lambda k: -counts_clust[k])[:15]
        handles = [Line2D([0],[0], marker='o', color='w', markerfacecolor=clust_cmap.get(c,(0.5,0.5,0.5)),
                          markersize=14, label=f'{c} ({counts_clust[c]:,})') for c in top_clust]
        axes[0,1].legend(handles=handles, loc='upper right', fontsize=12, framealpha=0.9, markerscale=1.2)
        axes[0,1].set_title('Cluster Labels', fontsize=20, fontweight='bold')
        axes[0,1].set_xlabel('x (µm)', fontsize=16); axes[0,1].set_ylabel('y (µm)', fontsize=16)
        axes[0,1].tick_params(labelsize=13)
        
        # Panel (0,2): Class
        colors_class = np.array([class_cmap.get(c, (0.5,0.5,0.5)) for c in class_label])
        img, ext = build_raster(x, y, colors_class)
        axes[0,2].imshow(img, extent=ext, aspect='equal', interpolation='nearest', origin='upper')
        axes[0,2].set_xlim(shared_xlim); axes[0,2].set_ylim(shared_ylim)
        counts_class = Counter(class_label)
        handles = [Line2D([0],[0], marker='o', color='w', markerfacecolor=class_cmap.get(c,(0.5,0.5,0.5)),
                          markersize=16, label=f'{class_display.get(c, c)} ({counts_class[c]:,})')
                   for c in ['Neuronal: Glutamatergic', 'Neuronal: GABAergic', 'Non-neuronal and Non-neural'] if c in counts_class]
        axes[0,2].legend(handles=handles, loc='upper right', fontsize=16, framealpha=0.9, markerscale=1.2)
        axes[0,2].set_title('Class Labels', fontsize=20, fontweight='bold')
        axes[0,2].set_xlabel('x (µm)', fontsize=16); axes[0,2].set_ylabel('y (µm)', fontsize=16)
        axes[0,2].tick_params(labelsize=13)
        
        # ---------- Row 2: Cortical Depth ----------
        
        # Panel (1,0): Continuous depth
        cmap_depth = plt.cm.viridis
        vmin, vmax = np.nanmin(pred_depth), np.nanmax(pred_depth)
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        colors_depth = cmap_depth(norm(pred_depth))[:, :3]
        img, ext = build_raster(x, y, colors_depth)
        axes[1,0].imshow(img, extent=ext, aspect='equal', interpolation='nearest', origin='upper')
        axes[1,0].set_xlim(shared_xlim); axes[1,0].set_ylim(shared_ylim)
        sm = plt.cm.ScalarMappable(cmap=cmap_depth, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=axes[1,0], shrink=0.6, pad=0.02)
        cbar.set_label('Predicted normalized depth\n(0 = pia, 1 = WM)', fontsize=16)
        cbar.ax.tick_params(labelsize=14)
        axes[1,0].set_title(f'Predicted Normalized Depth\n(range: {vmin:.3f} to {vmax:.3f})', fontsize=20, fontweight='bold')
        axes[1,0].set_xlabel('x (µm)', fontsize=16); axes[1,0].set_ylabel('y (µm)', fontsize=16)
        axes[1,0].tick_params(labelsize=13)
        
        # Panel (1,1): Discrete layers
        colors_layers = np.zeros((n_cells, 3))
        cell_layers = np.full(n_cells, '', dtype=object)
        for lname, (lo, hi) in layer_bins.items():
            mask = (pred_depth >= lo) & (pred_depth < hi)
            colors_layers[mask] = layer_colors[lname]
            cell_layers[mask] = lname
        
        img, ext = build_raster(x, y, colors_layers)
        axes[1,1].imshow(img, extent=ext, aspect='equal', interpolation='nearest', origin='upper')
        axes[1,1].set_xlim(shared_xlim); axes[1,1].set_ylim(shared_ylim)
        layer_counts = Counter(cell_layers)
        handles = [Line2D([0],[0], marker='o', color='w', markerfacecolor=layer_colors[l],
                          markersize=16, label=f'{l} ({layer_counts.get(l,0):,})') for l in layer_bins.keys()]
        axes[1,1].legend(handles=handles, loc='upper right', fontsize=16, framealpha=0.9, markerscale=1.2)
        axes[1,1].set_title('Discrete Layers\n(from predicted depth)', fontsize=20, fontweight='bold')
        axes[1,1].set_xlabel('x (µm)', fontsize=16); axes[1,1].set_ylabel('y (µm)', fontsize=16)
        axes[1,1].tick_params(labelsize=13)
        
        # Panel (1,2): Depth histogram with layer boundaries
        ax_hist = axes[1,2]
        ax_hist.set_facecolor('white')
        ax_hist.hist(pred_depth, bins=80, color='steelblue', edgecolor='none', alpha=0.8)
        
        # Layer boundary lines
        boundaries = [0.1, 0.3, 0.45, 0.65, 0.85]
        boundary_labels = ['L1|L2/3', 'L2/3|L4', 'L4|L5', 'L5|L6', 'L6|WM']
        ylim = ax_hist.get_ylim()
        for b, bl in zip(boundaries, boundary_labels):
            ax_hist.axvline(b, color='red', linestyle='--', linewidth=2, alpha=0.8)
            ax_hist.text(b, ylim[1]*0.95, bl,
                        rotation=90, ha='right', va='top', fontsize=14, color='red', fontweight='bold')
        
        # Add layer shading
        for lname, (lo, hi) in layer_bins.items():
            lo_clip = max(lo, pred_depth.min() - 0.05)
            hi_clip = min(hi, pred_depth.max() + 0.05)
            ax_hist.axvspan(lo_clip, hi_clip, alpha=0.15, color=layer_colors[lname])
            mid = (max(lo, pred_depth.min()- 0.02) + min(hi, pred_depth.max()+0.02)) / 2
            ax_hist.text(mid, ylim[1]*0.5, lname, ha='center', va='center',
                        fontsize=18, fontweight='bold', color=layer_colors[lname], alpha=0.7)
        
        ax_hist.set_xlabel('Predicted normalized depth', fontsize=18)
        ax_hist.set_ylabel('Number of cells', fontsize=18)
        ax_hist.set_title(f'Depth Distribution\n(n={n_cells:,} cells)', fontsize=20, fontweight='bold')
        ax_hist.tick_params(labelsize=14)
        
        # ============ SUPTITLE ============
        fig.suptitle(f'{sample_id} — Xenium Spatial Transcriptomics Summary\n'
                     f'n = {n_cells:,} cells | SEA-AD label transfer | MERFISH-trained depth model',
                     fontsize=24, fontweight='bold', y=0.995)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        out_path = f'{os.path.expanduser("~/Github/SCZ_Xenium")}/output/plots/{sample_id}_summary.png'
        plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        
        elapsed = time.time() - t0
        print(f"  {sample_id}: {n_cells:,} cells, {elapsed:.1f}s", flush=True)
        return sample_id, True
        
    except Exception as e:
        print(f"  ERROR {sample_id}: {e}", flush=True)
        import traceback
        traceback.print_exc()
        return sample_id, False

if __name__ == '__main__':
    h5ad_dir = os.path.join(os.path.expanduser("~/Github/SCZ_Xenium"), "output", "h5ad")
    h5ad_files = sorted(glob.glob(os.path.join(h5ad_dir, '*_annotated.h5ad')))
    
    out_dir = os.path.join(os.path.expanduser("~/Github/SCZ_Xenium"), "output", "plots")
    os.makedirs(out_dir, exist_ok=True)
    
    print(f"Found {len(h5ad_files)} samples")
    print(f"Generating combined summary figures with 4 workers...")
    
    t_start = time.time()
    with Pool(4) as pool:
        results = pool.map(process_sample, h5ad_files)
    
    t_total = time.time() - t_start
    n_ok = sum(1 for _, ok in results if ok)
    print(f"\nDone! {n_ok}/{len(results)} samples in {t_total:.0f}s")
    for sid, ok in results:
        status = 'OK' if ok else 'FAILED'
        print(f"  {sid}: {status}")
