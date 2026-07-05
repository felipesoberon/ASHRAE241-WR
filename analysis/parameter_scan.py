# parameter_scan.py
#
# Univariate sensitivity scan: fix all parameters at their median
# except one, sweep that one across its range, and observe the
# effect on infection probability P.
#
# This gives a simple "how much does each parameter matter?" curve
# per category. It is a first-order sensitivity measure -- it does
# not capture interactions between parameters (for that, see the
# exceed_vs_nonexceed or correlation_analysis scripts).
#
# Usage:
#   python analysis/parameter_scan.py [inputs.bin]
#       [--category Classroom] [--save FIG]
#
# Default inputs file: probabilityECAi_inputs.bin
# Default category: Classroom

import argparse
import numpy as np
import os
import sys
from inputs_reader import load_inputs, FIELD_NAMES

# Add parent dir for model.py (ECAi values)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from model import occupancy_params


def main():
    parser = argparse.ArgumentParser(
        description="Univariate parameter scan: sweep one parameter "
                    "at a time while holding others at median."
    )
    parser.add_argument(
        "infile", nargs="?", default="probabilityECAi_inputs.bin",
        help="Inputs .bin file (default: probabilityECAi_inputs.bin)"
    )
    parser.add_argument(
        "--category", type=str, default="Classroom",
        help="Category to analyze (default: Classroom)"
    )
    parser.add_argument(
        "--save", type=str, default=None,
        help="Save the figure to this file"
    )
    args = parser.parse_args()

    data = load_inputs(args.infile)

    if args.category not in data:
        print(f"Error: '{args.category}' not in data")
        print(f"Available: {', '.join(data.keys())}")
        sys.exit(1)

    arr = data[args.category]
    P_IDX = 8
    P = arr[:, P_IDX]
    ecai_val = occupancy_params.get(args.category, {}).get("ECAi", 0)

    # Parameters to scan (continuous ones only, skip n_infected
    # which is discrete)
    SCAN_FIELDS = [
        0,   # PBR_sample
        1,   # lambda_bio
        2,   # gamma
        5,   # QER_sum
        9,   # qer_PBR_qer
        10,  # qer_C_drop
        11,  # qer_d
        12,  # qer_E
        14,  # qer_GVL_ml
        16,  # qer_VF
        17,  # qer_RTD
        18,  # qer_DK
    ]
    scan_names = [FIELD_NAMES[fi] for fi in SCAN_FIELDS]

    print(f"\nCategory: {args.category}  (ECAi={ecai_val} L/s/p)")
    print(f"N simulations: {len(P)}")

    # Compute medians for all 21 fields
    medians = np.median(arr, axis=0)

    # For each parameter, scan from 10th to 90th percentile of the
    # observed range. We bin P by the parameter value and compute
    # the median P in each bin. This is a data-driven univariate
    # scan: we don't re-run the model, we exploit the fact that with
    # 1M simulations the marginal effect is visible by conditioning
    # on one parameter.

    n_bins = 20
    scan_results = {}

    for fi in SCAN_FIELDS:
        name = FIELD_NAMES[fi]
        vals = arr[:, fi]

        # Skip if all values are the same
        if vals.std() < 1e-15:
            scan_results[name] = None
            continue

        # Bin edges: 10th to 90th percentile to avoid extreme tails
        lo = np.percentile(vals, 10)
        hi = np.percentile(vals, 90)
        if lo == hi:
            scan_results[name] = None
            continue

        bin_edges = np.linspace(lo, hi, n_bins + 1)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        bin_p_medians = np.zeros(n_bins)

        for j in range(n_bins):
            mask = (vals >= bin_edges[j]) & (vals < bin_edges[j + 1])
            if j == n_bins - 1:
                mask = (vals >= bin_edges[j]) & (vals <= bin_edges[j + 1])
            if mask.sum() > 0:
                bin_p_medians[j] = np.median(P[mask])
            else:
                bin_p_medians[j] = np.nan

        scan_results[name] = (bin_centers, bin_p_medians)

        # Report the dynamic range (ratio of max to min median P)
        valid = bin_p_medians[~np.isnan(bin_p_medians)]
        if len(valid) > 1 and valid.min() > 0:
            dynamic_range = valid.max() / valid.min()
            print(f"  {name:<20s}  P range: {valid.min():.6f} to "
                  f"{valid.max():.6f}  ratio: {dynamic_range:.1f}x")
        else:
            print(f"  {name:<20s}  (insufficient variation)")

    # Plot
    if args.save:
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        mpl.use('Agg')

        valid_params = [(name, scan_results[name])
                        for name in scan_names
                        if scan_results[name] is not None]
        n = len(valid_params)

        n_cols = 3
        n_rows = int(np.ceil(n / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols,
                                figsize=(15, 3 * n_rows),
                                squeeze=False)

        for i, (name, (centers, p_meds)) in enumerate(valid_params):
            ax = axes[i // n_cols][i % n_cols]
            ax.plot(centers, p_meds * 100, 'o-', markersize=3,
                    color='steelblue')
            ax.axhline(0.1, color='red', linestyle='--', linewidth=0.5)
            ax.set_title(name, fontsize=9)
            ax.set_ylabel('Median P (%)', fontsize=7)
            ax.tick_params(labelsize=7)
            ax.set_yscale('log')
            ax.set_ylim(max(p_meds.min() * 100 * 0.5, 1e-4),
                        max(p_meds.max() * 100 * 2, 1))

        # Hide unused subplots
        for i in range(n, n_rows * n_cols):
            axes[i // n_cols][i % n_cols].set_visible(False)

        fig.suptitle(f'Univariate parameter scan: {args.category} '
                     f'(ECAi={ecai_val} L/s/p)', fontsize=11)
        fig.tight_layout(pad=2.0)
        fig.savefig(args.save, dpi=150)
        print(f"\nFigure saved to {args.save}")


if __name__ == "__main__":
    main()