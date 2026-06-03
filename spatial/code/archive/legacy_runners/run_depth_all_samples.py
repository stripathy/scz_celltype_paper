import anndata as ad
import numpy as np
import pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.neighbors import NearestNeighbors
from collections import Counter
from multiprocessing import Pool
import glob
import time
import os

_BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")

# ---- Globals set by initializer ----
_model = None
_sub_to_idx = None
_subclass_names = None
_n_sub = None
_K = None

def init_worker(model_path):
    global _model, _sub_to_idx, _subclass_names, _n_sub, _K
    with open(model_path, 'rb') as f:
        mdata = pickle.load(f)
    _model = mdata['model']
    _sub_to_idx = mdata['sub_to_idx']
    _subclass_names = mdata['subclass_names']
    _n_sub = mdata['n_sub']
    _K = mdata['K']

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

def process_sample(h5ad_path):
    t0 = time.time()
    sample_id = os.path.basename(h5ad_path).replace('_annotated.h5ad', '')
    
    try:
        # Load
        adata = ad.read_h5ad(h5ad_path)
        coords = adata.obsm['spatial']
        subclass = adata.obs['subclass_label'].values.astype(str)
        n_cells = adata.shape[0]
        
        # Map subclass to index
        sub_idx = np.array([_sub_to_idx.get(s, -1) for s in subclass])
        
        # Build KNN
        nn = NearestNeighbors(n_neighbors=_K+1, algorithm='ball_tree')
        nn.fit(coords)
        _, indices = nn.kneighbors(coords)
        neighbor_idx = indices[:, 1:]
        
        # Build features
        n_features = _n_sub * 2
        features = np.zeros((n_cells, n_features))
        for i in range(n_cells):
            neigh_sub = sub_idx[neighbor_idx[i]]
            for ns in neigh_sub:
                if ns >= 0:
                    features[i, ns] += 1
            features[i, :_n_sub] /= _K
            own = sub_idx[i]
            if own >= 0:
                features[i, _n_sub + own] = 1
        
        # Predict — NO CLAMPING
        pred_depth = _model.predict(features)
        
        # Save to h5ad
        adata.obs['predicted_norm_depth'] = pred_depth
        adata.write_h5ad(h5ad_path)
        
        # --- Generate 3-panel plot ---
        x = coords[:, 0]
        y = coords[:, 1]
        
        fig, axes = plt.subplots(1, 3, figsize=(48, 16), facecolor='white')
        
        # Panel 1: Continuous depth
        cmap = plt.cm.viridis
        vmin, vmax = pred_depth.min(), pred_depth.max()
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        colors1 = np.array([cmap(norm(d))[:3] for d in pred_depth])
        img1, ext1 = build_raster(x, y, colors1)
        axes[0].imshow(img1, extent=ext1, aspect='equal', interpolation='nearest', origin='upper')
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=axes[0], shrink=0.6, pad=0.02)
        cbar.set_label('Predicted normalized depth (0=pia, 1=WM)', fontsize=13)
        axes[0].set_title(f'Predicted Normalized Depth\n(range: {vmin:.3f} to {vmax:.3f})', fontsize=18, fontweight='bold')
        axes[0].set_xlabel('x (µm)'); axes[0].set_ylabel('y (µm)')
        
        # Panel 2: Discrete layers
        layer_bins = {
            'L1': (-np.inf, 0.1), 'L2/3': (0.1, 0.3), 'L4': (0.3, 0.45),
            'L5': (0.45, 0.65), 'L6': (0.65, 0.85), 'WM': (0.85, np.inf),
        }
        layer_colors = {
            'L1': (0.9, 0.3, 0.3), 'L2/3': (0.3, 0.8, 0.3), 'L4': (0.3, 0.3, 0.9),
            'L5': (0.9, 0.6, 0.1), 'L6': (0.7, 0.3, 0.8), 'WM': (0.5, 0.5, 0.5),
        }
        colors2 = np.zeros((n_cells, 3))
        cell_layers = []
        for lname, (lo, hi) in layer_bins.items():
            mask = (pred_depth >= lo) & (pred_depth < hi)
            colors2[mask] = layer_colors[lname]
            cell_layers.extend([lname] * mask.sum())
        
        img2, ext2 = build_raster(x, y, colors2)
        axes[1].imshow(img2, extent=ext2, aspect='equal', interpolation='nearest', origin='upper')
        layer_counts = Counter(cell_layers)
        handles2 = [Line2D([0],[0], marker='o', color='w', markerfacecolor=layer_colors[l],
                           markersize=12, label=f'{l} ({layer_counts.get(l,0):,})') for l in layer_bins.keys()]
        axes[1].legend(handles=handles2, loc='upper right', fontsize=12, framealpha=0.9)
        axes[1].set_title('Discrete Layers from Predicted Depth', fontsize=18, fontweight='bold')
        axes[1].set_xlabel('x (µm)'); axes[1].set_ylabel('y (µm)')
        
        # Panel 3: Subclass
        unique_sub = sorted(set(subclass))
        palette = np.vstack([plt.cm.tab20(np.linspace(0,1,20)), plt.cm.Set3(np.linspace(0,1,12))])
        sub_cmap = {s: palette[i%len(palette)][:3] for i, s in enumerate(unique_sub)}
        colors3 = np.array([sub_cmap.get(s, (0.5,0.5,0.5)) for s in subclass])
        img3, ext3 = build_raster(x, y, colors3)
        axes[2].imshow(img3, extent=ext3, aspect='equal', interpolation='nearest', origin='upper')
        counts_sub = Counter(subclass)
        top_sub = sorted(counts_sub.keys(), key=lambda k: -counts_sub[k])[:20]
        handles3 = [Line2D([0],[0], marker='o', color='w', markerfacecolor=sub_cmap.get(s,(0.5,0.5,0.5)),
                           markersize=10, label=f'{s} ({counts_sub[s]:,})') for s in top_sub]
        axes[2].legend(handles=handles3, loc='upper right', fontsize=8, framealpha=0.9)
        axes[2].set_title(f'Cell Types by Subclass (n={n_cells:,})', fontsize=18, fontweight='bold')
        axes[2].set_xlabel('x (µm)'); axes[2].set_ylabel('y (µm)')
        
        fig.suptitle(f'{sample_id} — Predicted Normalized Cortical Depth',
                     fontsize=22, fontweight='bold', y=1.02)
        plt.tight_layout()
        plot_path = f'{_BASE_DIR}/output/plots/{sample_id}_depth.png'
        plt.savefig(plot_path, dpi=200, bbox_inches='tight')
        plt.close()
        
        elapsed = time.time() - t0
        print(f"  {sample_id}: {n_cells:,} cells, depth [{pred_depth.min():.3f}, {pred_depth.max():.3f}], {elapsed:.1f}s", flush=True)
        return sample_id, n_cells, pred_depth.min(), pred_depth.max(), True
        
    except Exception as e:
        print(f"  ERROR {sample_id}: {e}", flush=True)
        import traceback
        traceback.print_exc()
        return sample_id, 0, 0, 0, False

if __name__ == '__main__':
    model_path = f'{_BASE_DIR}/output/depth_model_normalized.pkl'
    h5ad_dir = f'{_BASE_DIR}/output/h5ad/'
    h5ad_files = sorted(glob.glob(os.path.join(h5ad_dir, '*_annotated.h5ad')))
    
    print(f"Found {len(h5ad_files)} samples")
    print(f"Starting parallel processing with 4 workers...")
    
    t_start = time.time()
    with Pool(4, initializer=init_worker, initargs=(model_path,)) as pool:
        results = pool.map(process_sample, h5ad_files)
    
    t_total = time.time() - t_start
    print(f"\n{'='*60}")
    print(f"All {len(results)} samples processed in {t_total:.0f}s")
    
    total_cells = 0
    for sid, nc, dmin, dmax, ok in results:
        total_cells += nc
        status = "OK" if ok else "FAILED"
        print(f"  {sid}: {nc:,} cells, depth [{dmin:.3f}, {dmax:.3f}] — {status}")
    print(f"\nTotal cells: {total_cells:,}")
    
    # Now merge all h5ad files into combined file
    print(f"\nMerging all samples into combined h5ad...")
    t_merge = time.time()
    adatas = []
    for h5 in h5ad_files:
        a = ad.read_h5ad(h5)
        adatas.append(a)
    
    combined = ad.concat(adatas, join='outer')
    combined_path = f'{_BASE_DIR}/output/all_samples_annotated.h5ad'
    combined.write_h5ad(combined_path)
    
    print(f"Merge took {time.time()-t_merge:.0f}s")
    print(f"Combined shape: {combined.shape}")
    print(f"Columns in .obs: {list(combined.obs.columns)}")
    print(f"'predicted_norm_depth' present: {'predicted_norm_depth' in combined.obs.columns}")
    print(f"Depth stats: min={combined.obs['predicted_norm_depth'].min():.3f}, max={combined.obs['predicted_norm_depth'].max():.3f}, mean={combined.obs['predicted_norm_depth'].mean():.3f}")
    
    file_size = os.path.getsize(combined_path) / 1e9
    print(f"\nSaved: {combined_path} ({file_size:.2f} GB)")
    print("Done!")
