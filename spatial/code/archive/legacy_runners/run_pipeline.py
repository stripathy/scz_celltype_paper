"""
Sequential pipeline: label transfer, depth prediction, and plotting.

Processes Xenium samples one at a time (useful for debugging).
For parallel processing, use run_parallel.py instead.

Usage:
    python run_pipeline.py --data_dir /path/to/xenium_data \
                           --reference /path/to/seaad_reference.h5ad \
                           --depth_model /path/to/depth_model.pkl \
                           --output_dir /path/to/output
"""

import os
import sys
import argparse
import anndata as ad

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from loading import load_xenium_sample, discover_samples
from label_transfer import load_reference, annotate_sample, get_seaad_colors
from depth_model import load_model, predict_depth
from plotting import plot_summary


def process_sample(sample_info, ref, colors, depth_model, output_dir):
    """Process a single sample: load -> annotate -> depth -> plot -> save."""
    sid = sample_info["sample_id"]
    print(f"\n{'='*60}")
    print(f"Processing {sid}")
    print(f"{'='*60}")

    # Load
    adata = load_xenium_sample(
        sample_info["h5_path"], sample_info["boundaries_path"]
    )
    print(f"  {adata.shape[0]:,} cells x {adata.shape[1]} genes")

    # Label transfer
    annotated = annotate_sample(adata, ref)

    # Depth prediction
    if depth_model is not None:
        print("  Predicting cortical depth...")
        pred_depth = predict_depth(annotated, depth_model)
        annotated.obs['predicted_norm_depth'] = pred_depth
        print(f"  Depth range: [{pred_depth.min():.3f}, {pred_depth.max():.3f}]")

    # Save h5ad
    h5ad_dir = os.path.join(output_dir, "h5ad")
    os.makedirs(h5ad_dir, exist_ok=True)
    h5ad_path = os.path.join(h5ad_dir, f"{sid}_annotated.h5ad")
    annotated.write_h5ad(h5ad_path)
    print(f"  Saved: {h5ad_path}")

    # Summary figure
    plot_dir = os.path.join(output_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    plot_path = os.path.join(plot_dir, f"{sid}_summary.png")
    plot_summary(annotated, plot_path, model_bundle=depth_model)

    return annotated


def main():
    parser = argparse.ArgumentParser(
        description="Sequential Xenium pipeline: label transfer + depth"
    )
    parser.add_argument("--data_dir", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--depth_model", default=None,
                        help="Path to trained depth model .pkl")
    parser.add_argument("--output_dir", default=None)
    parser.add_argument("--samples", nargs="*", default=None)
    parser.add_argument("--save_combined", action="store_true")
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = os.path.join(args.data_dir, "output")
    os.makedirs(args.output_dir, exist_ok=True)

    # Load reference
    ref = load_reference(args.reference)
    colors = get_seaad_colors(ref)

    # Load depth model
    depth_model = None
    if args.depth_model and os.path.exists(args.depth_model):
        depth_model = load_model(args.depth_model)

    # Discover samples
    all_samples = discover_samples(args.data_dir)
    if args.samples:
        all_samples = [s for s in all_samples
                       if s["sample_id"] in args.samples]
    print(f"\nWill process {len(all_samples)} samples")

    # Process
    annotated_list = []
    for info in all_samples:
        annotated = process_sample(info, ref, colors, depth_model,
                                   args.output_dir)
        annotated_list.append(annotated)

    # Combined
    if args.save_combined and annotated_list:
        print("\nSaving combined h5ad...")
        combined = ad.concat(annotated_list, join="outer")
        combined_path = os.path.join(args.output_dir,
                                     "all_samples_annotated.h5ad")
        combined.write_h5ad(combined_path)
        print(f"Saved: {combined_path} ({combined.shape})")

    print(f"\nDone! Processed {len(annotated_list)} samples.")


if __name__ == "__main__":
    main()
