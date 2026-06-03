"""
Cortical depth prediction from local cell type neighborhood composition.

Trains a GradientBoostingRegressor on SEA-AD MERFISH data to predict
normalized depth from pia (0 = pial surface, 1 = white matter boundary)
using K-nearest-neighbor subclass composition features.

The model can then be applied to Xenium data to assign a predicted
cortical depth to every cell, enabling depth-stratified analyses that
control for variable tissue sampling geometry across sections.

Key insight: predictions are NOT clamped to [0,1] — cells in white matter
may receive depth > 1 and cells above pia may receive depth < 0.

Note: OOD (out-of-distribution) detection was previously implemented here
via 1-NN distance thresholds but was never used in the pipeline. Domain
classification (Extra-cortical, Vascular, WM) is now handled by
banksy_domains.py using BANKSY spatial clustering, which correctly
identifies L1 border cells and white matter.
"""

import numpy as np
import pickle
import time
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import r2_score, mean_absolute_error
from scipy.stats import pearsonr


# Default depth strata for downstream analyses
# Boundaries derived from excitatory neuron marker crossovers in SEA-AD
# MERFISH (341K cells) and validated against Xenium Control samples.
# Method: for each pair of adjacent layers, compute the depth where the
# smoothed pairwise fraction A/(A+B) = 0.5 using excitatory marker neurons
# (L2/3 IT, L4 IT, L5 IT/ET/NP, L6 CT/IT/Car3/b). L1/L2-3 uses onset of
# excitatory density (>25% of peak); L6/WM uses dropoff (<10% of peak).
# See code/analysis/derive_layer_boundaries.py for derivation details.
DEPTH_STRATA = {
    'L2/3': (0.1225, 0.4696),
    'L4': (0.4696, 0.5443),
    'L5': (0.5443, 0.7079),
    'L6': (0.7079, 0.9275),
}

# Discrete layer bins (L1 through WM)
LAYER_BINS = {
    'L1': (-np.inf, 0.1225),
    'L2/3': (0.1225, 0.4696),
    'L4': (0.4696, 0.5443),
    'L5': (0.5443, 0.7079),
    'L6': (0.7079, 0.9275),
    'WM': (0.9275, np.inf),
}

LAYER_COLORS = {
    'L1': (0.9, 0.3, 0.3),
    'L2/3': (0.3, 0.8, 0.3),
    'L4': (0.3, 0.3, 0.9),
    'L5': (0.9, 0.6, 0.1),
    'L6': (0.7, 0.3, 0.8),
    'WM': (0.5, 0.5, 0.5),
    'Extra-cortical': (0.8, 0.8, 0.2),
}


def _vectorized_neighbor_fractions(nn_idx_neighbors, sub_idx_all, n_sub):
    """
    Vectorized computation of neighbor subclass fractions.

    Parameters
    ----------
    nn_idx_neighbors : np.ndarray (n_cells, K)
        Indices of K nearest neighbors for each cell (self excluded).
    sub_idx_all : np.ndarray (n_total,)
        Integer subclass index for each cell (-1 for unknown).
    n_sub : int
        Number of subclass categories.

    Returns
    -------
    np.ndarray (n_cells, n_sub)
        Fraction of each subclass among neighbors.
    """
    n_cells, K = nn_idx_neighbors.shape
    # Look up subclass index for every neighbor
    neighbor_types = sub_idx_all[nn_idx_neighbors]  # (n_cells, K)
    # One-hot encode and sum across neighbors
    # Use np.eye trick: create indicator matrix then sum
    valid = neighbor_types >= 0  # mask out unknown types
    neighbor_types_safe = np.where(valid, neighbor_types, 0)

    # Vectorized bincount per row using broadcasting
    fractions = np.zeros((n_cells, n_sub), dtype=np.float32)
    for k in range(K):
        col = neighbor_types_safe[:, k]
        mask = valid[:, k]
        np.add.at(fractions, (np.arange(n_cells)[mask], col[mask]), 1)

    # Normalize by actual neighbor count
    n_valid = valid.sum(axis=1, keepdims=True).astype(np.float32)
    n_valid = np.maximum(n_valid, 1)  # avoid division by zero
    fractions /= n_valid
    return fractions


def build_neighborhood_features(coords, subclass_labels, subclass_names,
                                 K=50, sections=None):
    """
    Build K-nearest-neighbor subclass composition features.

    For each cell, finds its K spatial nearest neighbors (within the same
    section if sections is provided), computes the fraction of each subclass
    among neighbors, and adds a one-hot encoding of the cell's own type.

    Parameters
    ----------
    coords : np.ndarray
        Spatial coordinates (n_cells x 2).
    subclass_labels : np.ndarray of str
        Subclass label per cell.
    subclass_names : list of str
        Ordered list of all possible subclass names.
    K : int
        Number of nearest neighbors.
    sections : np.ndarray of str, optional
        Section ID per cell. If provided, neighbors are found within
        each section separately (required for MERFISH multi-section data).

    Returns
    -------
    np.ndarray
        Feature matrix (n_cells x 2*n_subclass).
        First n_sub columns: neighbor fractions.
        Last n_sub columns: own-type one-hot.
    """
    n_cells = len(subclass_labels)
    n_sub = len(subclass_names)
    sub_to_idx = {s: i for i, s in enumerate(subclass_names)}

    sub_idx = np.array([sub_to_idx.get(s, -1) for s in subclass_labels])
    features = np.zeros((n_cells, n_sub * 2), dtype=np.float32)

    if sections is not None:
        # Build per-section
        unique_sections = np.unique(sections)
        for sec_i, sec in enumerate(unique_sections):
            if sec_i % 20 == 0:
                print(f"  Section {sec_i+1}/{len(unique_sections)}")
            sec_mask = sections == sec
            sec_indices = np.where(sec_mask)[0]
            sec_coords = coords[sec_indices]

            k_actual = min(K + 1, len(sec_indices))
            nn = NearestNeighbors(n_neighbors=k_actual, algorithm='ball_tree')
            nn.fit(sec_coords)
            _, nn_idx_local = nn.kneighbors(sec_coords)

            # Map local indices to global
            nn_neighbors_local = nn_idx_local[:, 1:K+1]  # exclude self
            nn_neighbors_global = sec_indices[nn_neighbors_local]

            # Vectorized neighbor fractions
            fracs = _vectorized_neighbor_fractions(
                nn_neighbors_global, sub_idx, n_sub
            )
            features[sec_indices, :n_sub] = fracs

            # Own type one-hot (vectorized)
            own_types = sub_idx[sec_indices]
            valid_own = own_types >= 0
            features[sec_indices[valid_own], n_sub + own_types[valid_own]] = 1
    else:
        # Single section / sample
        nn = NearestNeighbors(n_neighbors=K+1, algorithm='ball_tree')
        nn.fit(coords)
        _, nn_idx = nn.kneighbors(coords)

        nn_neighbors = nn_idx[:, 1:]  # exclude self
        fracs = _vectorized_neighbor_fractions(nn_neighbors, sub_idx, n_sub)
        features[:, :n_sub] = fracs

        # Own type one-hot (vectorized)
        valid_own = sub_idx >= 0
        features[np.where(valid_own)[0], n_sub + sub_idx[valid_own]] = 1

    return features


def train_depth_model(merfish_adata, test_donors=None, K=50,
                      n_estimators=300, max_depth=5, learning_rate=0.1):
    """
    Train a depth prediction model on MERFISH data.

    Uses cells with 'Normalized depth from pia' annotations as training
    targets, with K-nearest-neighbor subclass composition as features.

    Parameters
    ----------
    merfish_adata : anndata.AnnData
        SEA-AD MERFISH dataset with 'Subclass', 'Normalized depth from pia',
        'Section', 'Donor ID' in .obs, and spatial coordinates in
        .obsm['X_spatial_raw'].
    test_donors : list of str, optional
        Donor IDs to hold out for testing. If None, uses the 3 donors
        with the fewest depth-annotated cells.
    K : int
        Number of nearest neighbors for features.
    n_estimators : int
        Number of boosting rounds.
    max_depth : int
        Max tree depth.
    learning_rate : float
        Boosting learning rate.

    Returns
    -------
    dict
        Model bundle containing: 'model', 'subclass_names', 'sub_to_idx',
        'K', 'n_sub', 'feature_names', 'train_donors', 'test_donors',
        'train_r2', 'test_r2', 'train_mae', 'test_mae', 'target'.
    """
    print("Preparing MERFISH depth model training data...")

    # Get depth annotations
    depth = merfish_adata.obs['Normalized depth from pia'].values.astype(float)
    has_depth = ~np.isnan(depth)
    print(f"  Cells with depth annotation: {has_depth.sum():,} / {len(depth):,}")

    # Subclass info
    subclass = merfish_adata.obs['Subclass'].values.astype(str)
    subclass_names = sorted(set(subclass))
    n_sub = len(subclass_names)
    sub_to_idx = {s: i for i, s in enumerate(subclass_names)}

    # Donor info
    donor_col = 'Donor ID' if 'Donor ID' in merfish_adata.obs.columns else 'Specimen Barcode'
    donors = merfish_adata.obs[donor_col].values.astype(str)
    sections = merfish_adata.obs['Section'].values.astype(str)

    # Select test donors
    if test_donors is None:
        donor_depth_counts = {}
        for d in np.unique(donors[has_depth]):
            donor_depth_counts[d] = (donors[has_depth] == d).sum()
        all_donors = sorted(donor_depth_counts.keys(),
                            key=lambda d: -donor_depth_counts[d])
        test_donors = all_donors[-3:]

    train_donors = sorted(set(np.unique(donors[has_depth])) - set(test_donors))
    print(f"  Train donors: {len(train_donors)}, Test donors: {len(test_donors)}")

    # Build features for depth-annotated cells
    print(f"  Building K={K} neighborhood features...")
    t0 = time.time()
    coords = merfish_adata.obsm['X_spatial_raw']

    features = build_neighborhood_features(
        coords[has_depth], subclass[has_depth], subclass_names,
        K=K, sections=sections[has_depth]
    )
    depths = depth[has_depth]
    feat_donors = donors[has_depth]
    print(f"  Feature building: {time.time()-t0:.0f}s")

    # Train/test split
    train_mask = np.isin(feat_donors, train_donors)
    test_mask = np.isin(feat_donors, test_donors)
    X_train, y_train = features[train_mask], depths[train_mask]
    X_test, y_test = features[test_mask], depths[test_mask]
    print(f"  Train: {X_train.shape[0]:,}, Test: {X_test.shape[0]:,}")

    # Train
    print(f"  Training GBR (n_estimators={n_estimators})...")
    t1 = time.time()
    model = GradientBoostingRegressor(
        n_estimators=n_estimators, max_depth=max_depth,
        learning_rate=learning_rate, subsample=0.8,
        random_state=42, min_samples_leaf=20
    )
    model.fit(X_train, y_train)
    print(f"  Training: {time.time()-t1:.0f}s")

    # Evaluate (no clamping)
    y_pred_train = model.predict(X_train)
    y_pred_test = model.predict(X_test)

    metrics = {
        'train_r2': r2_score(y_train, y_pred_train),
        'test_r2': r2_score(y_test, y_pred_test),
        'train_mae': mean_absolute_error(y_train, y_pred_train),
        'test_mae': mean_absolute_error(y_test, y_pred_test),
        'train_r': pearsonr(y_train, y_pred_train)[0],
        'test_r': pearsonr(y_test, y_pred_test)[0],
    }
    print(f"\n  Train: R²={metrics['train_r2']:.3f}, MAE={metrics['train_mae']:.4f}")
    print(f"  Test:  R²={metrics['test_r2']:.3f}, MAE={metrics['test_mae']:.4f}")

    feature_names = ([f'neigh_{s}' for s in subclass_names] +
                     [f'own_{s}' for s in subclass_names])

    return {
        'model': model,
        'subclass_names': subclass_names,
        'sub_to_idx': sub_to_idx,
        'K': K,
        'n_sub': n_sub,
        'feature_names': feature_names,
        'target': 'normalized_depth_from_pia',
        'train_donors': train_donors,
        'test_donors': test_donors,
        **metrics,
    }


def save_model(model_bundle, path):
    """Save model bundle to pickle file."""
    with open(path, 'wb') as f:
        pickle.dump(model_bundle, f)
    print(f"Saved model: {path}")


def load_model(path):
    """Load model bundle from pickle file."""
    with open(path, 'rb') as f:
        return pickle.load(f)


def predict_depth(adata, model_bundle, subclass_col='subclass_label'):
    """
    Predict normalized cortical depth for all cells in a Xenium sample.

    Builds K-nearest-neighbor features from spatial coordinates and
    subclass labels, then applies the trained model. Predictions are
    NOT clamped to [0,1].

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated Xenium sample with subclass labels in .obs and
        spatial coordinates in .obsm['spatial'].
    model_bundle : dict
        Trained model from train_depth_model() or load_model().
    subclass_col : str
        Column name for subclass labels in adata.obs.

    Returns
    -------
    np.ndarray
        Predicted normalized depth per cell.
    """
    model = model_bundle['model']
    subclass_names = model_bundle['subclass_names']
    K = model_bundle['K']

    coords = adata.obsm['spatial']
    subclass = adata.obs[subclass_col].values.astype(str)

    features = build_neighborhood_features(
        coords, subclass, subclass_names, K=K, sections=None
    )

    return model.predict(features)


def assign_discrete_layers(depths, layer_bins=None):
    """
    Assign discrete cortical layer labels from continuous depth values.

    Parameters
    ----------
    depths : np.ndarray
        Predicted normalized depth values.
    layer_bins : dict, optional
        {layer_name: (lower, upper)} boundaries. Defaults to LAYER_BINS.

    Returns
    -------
    np.ndarray of str
        Layer label per cell.
    """
    if layer_bins is None:
        layer_bins = LAYER_BINS

    layers = np.full(len(depths), '', dtype=object)
    for lname, (lo, hi) in layer_bins.items():
        mask = (depths >= lo) & (depths < hi)
        layers[mask] = lname
    return layers


# ── Spatial layer smoothing constants ────────────────────────────────

SMOOTH_K = 30                   # spatial neighbors for smoothing
SMOOTH_ROUNDS = 2               # majority-vote iterations

VASC_TRIM_CORTICAL_THRESH = 0.33   # trim Vascular if >33% neighbors in L2/3–L6
VASC_TRIM_WM_L1_THRESH = 0.66     # trim Vascular if >66% neighbors in any non-Vasc

L1_PROMOTE_DEPTH_THRESH = 0.20     # promote banksy_is_l1 cells to L1 if depth < 0.20
L1_PROMOTE_NBR_THRESH = 0.05       # … and ≥5% of neighbors are already L1
L1_ISOLATED_BANKSY_THRESH = 0.05   # banksy L1 cells need ≥5% L1 neighbors to stay
L1_ISOLATED_OTHER_THRESH = 0.20    # non-banksy L1 cells need ≥20% L1 neighbors

# Layer categories used for indexing in smoothing
_SMOOTH_LAYER_ORDER = ['L1', 'L2/3', 'L4', 'L5', 'L6', 'WM', 'Vascular']


def smooth_layers_spatial(coords, layers, domains, is_l1_banksy, depths,
                          k=None, n_rounds=None, verbose=True):
    """
    Spatially smooth layer assignments using 3-step pipeline:

    1. Within-domain majority vote — smooths cortical layer boundaries
       without crossing BANKSY domain borders.
    2. Vascular border trim — reassigns border Vascular cells to cortical
       layers when a sufficient fraction of spatial neighbors are in
       cortical layers (L2/3–L6).
    3. BANKSY-anchored L1 contiguity — promotes banksy_is_l1 cells with
       shallow depth to L1, and removes isolated non-BANKSY L1 cells.

    Parameters
    ----------
    coords : np.ndarray (n_cells, 2)
        Spatial coordinates.
    layers : np.ndarray of str (n_cells,)
        Initial layer assignments (from assign_discrete_layers + Vascular
        override). Values must be in _SMOOTH_LAYER_ORDER.
    domains : np.ndarray of str (n_cells,)
        BANKSY domain labels ('Cortical', 'Vascular', 'WM').
    is_l1_banksy : np.ndarray of bool (n_cells,)
        BANKSY L1 border flag per cell.
    depths : np.ndarray (n_cells,)
        Predicted normalized depth per cell.
    k : int, optional
        Number of spatial neighbors. Default: SMOOTH_K (30).
    n_rounds : int, optional
        Number of within-domain smoothing rounds. Default: SMOOTH_ROUNDS (2).
    verbose : bool
        Print summary statistics.

    Returns
    -------
    np.ndarray of str (n_cells,)
        Smoothed layer assignments.
    """
    if k is None:
        k = SMOOTH_K
    if n_rounds is None:
        n_rounds = SMOOTH_ROUNDS

    n_cells = len(layers)
    layer_to_idx = {l: i for i, l in enumerate(_SMOOTH_LAYER_ORDER)}
    idx_to_layer = {i: l for l, i in layer_to_idx.items()}
    n_cats = len(_SMOOTH_LAYER_ORDER)

    CORTICAL_IDXS = {layer_to_idx[l] for l in ['L2/3', 'L4', 'L5', 'L6']}
    VASC_IDX = layer_to_idx['Vascular']
    WM_IDX = layer_to_idx['WM']
    L1_IDX = layer_to_idx['L1']

    # Encode layers as integers
    layer_int = np.array([layer_to_idx.get(l, 0) for l in layers])
    original = layer_int.copy()

    # Build spatial k-NN graph
    t0 = time.time()
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm='ball_tree')
    nn.fit(coords)
    _, indices = nn.kneighbors(coords)
    if verbose:
        print(f"    Spatial k-NN (k={k}): {time.time()-t0:.1f}s")

    # ── Step 1: Within-domain majority vote ───────────────────────
    t1 = time.time()
    for _ in range(n_rounds):
        nbr_layers = layer_int[indices]  # (n_cells, k+1)
        # Mask out cross-domain neighbors
        domain_self = domains[:, None]
        domain_nbrs = domains[indices]
        nbr_layers = nbr_layers.copy()
        nbr_layers[domain_self != domain_nbrs] = -1

        new = np.empty(n_cells, dtype=int)
        for i in range(n_cells):
            valid = nbr_layers[i][nbr_layers[i] >= 0]
            if len(valid) == 0:
                new[i] = layer_int[i]
            else:
                counts = np.bincount(valid, minlength=n_cats)
                new[i] = np.argmax(counts)
        layer_int = new

    n_step1 = (layer_int != original).sum()
    if verbose:
        print(f"    Step 1 (within-domain smooth, {n_rounds} rounds): "
              f"{n_step1:,} changed ({time.time()-t1:.1f}s)")

    # ── Step 2: Vascular border trim ──────────────────────────────
    t2 = time.time()
    n_vasc_trimmed = 0
    for i in range(n_cells):
        if layer_int[i] != VASC_IDX:
            continue
        nbr = layer_int[indices[i, 1:]]
        n_cortical = sum(1 for l in nbr if l in CORTICAL_IDXS)
        n_wm_l1 = sum(1 for l in nbr if l in (WM_IDX, L1_IDX))

        do_trim = False
        if n_cortical / k >= VASC_TRIM_CORTICAL_THRESH:
            do_trim = True
        elif (n_cortical + n_wm_l1) / k >= VASC_TRIM_WM_L1_THRESH:
            do_trim = True

        if do_trim:
            non_vasc = nbr[nbr != VASC_IDX]
            if len(non_vasc) > 0:
                counts = np.bincount(non_vasc, minlength=n_cats)
                layer_int[i] = np.argmax(counts)
                n_vasc_trimmed += 1

    if verbose:
        print(f"    Step 2 (vascular trim): {n_vasc_trimmed:,} cells "
              f"({time.time()-t2:.1f}s)")

    # ── Step 3: BANKSY-anchored L1 contiguity ─────────────────────
    t3 = time.time()
    n_l1_promoted = 0
    n_l1_removed = 0

    # 3a: Promote banksy_is_l1 cells with shallow depth to L1
    for i in range(n_cells):
        if (is_l1_banksy[i]
                and layer_int[i] != L1_IDX
                and depths[i] < L1_PROMOTE_DEPTH_THRESH):
            nbr = layer_int[indices[i, 1:]]
            l1_frac = (nbr == L1_IDX).sum() / k
            if l1_frac >= L1_PROMOTE_NBR_THRESH:
                layer_int[i] = L1_IDX
                n_l1_promoted += 1

    # 3b: Remove isolated L1 cells
    for i in range(n_cells):
        if layer_int[i] != L1_IDX:
            continue
        nbr = layer_int[indices[i, 1:]]
        l1_frac = (nbr == L1_IDX).sum() / k
        thresh = (L1_ISOLATED_BANKSY_THRESH if is_l1_banksy[i]
                  else L1_ISOLATED_OTHER_THRESH)
        if l1_frac < thresh:
            non_l1 = nbr[nbr != L1_IDX]
            if len(non_l1) > 0:
                counts = np.bincount(non_l1, minlength=n_cats)
                layer_int[i] = np.argmax(counts)
                n_l1_removed += 1

    if verbose:
        print(f"    Step 3 (L1 contiguity): +{n_l1_promoted:,} promoted, "
              f"-{n_l1_removed:,} removed ({time.time()-t3:.1f}s)")

    # Convert back to string labels
    smoothed = np.array([idx_to_layer[i] for i in layer_int])

    n_total_changed = (smoothed != layers).sum()
    if verbose:
        print(f"    Total changed: {n_total_changed:,} / {n_cells:,} "
              f"({100*n_total_changed/n_cells:.1f}%)")

    return smoothed


