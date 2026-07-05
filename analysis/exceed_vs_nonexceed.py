# exceed_vs_nonexceed.py
#
# Driver analysis: compare input parameter distributions between
# exceed (P > 0.1%) and non-exceed (P <= 0.1%) simulation groups.
#
# For each category, splits the 1M simulations into two groups based
# on whether the infection probability exceeds the target (0.1%).
# For each input parameter, computes the mean in each group and the
# ratio (exceed/non-exceed). Parameters with the largest ratio shift
# are the primary tail drivers.
#
# Optionally produces a tornado plot showing the mean-ratio per
# parameter per category.
#
# Usage:
#   python analysis/exceed_vs_nonexceed.py [inputs.bin]
#       [--target 0.1] [--csv outfile.csv] [--save FIG]
#       [--category CATEGORY]  (if omitted, does all categories)
#
# Default inputs file: probabilityECAi_inputs.bin

import argparse
import numpy as np
import csv
import os
import sys
from inputs_reader import load_inputs, FIELD_NAMES

# Add parent dir for model.py (ECAi values)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from model import occupancy_params


def main():
    parser = argparse.ArgumentParser(
        description="Exceed vs non-exceed input parameter comparison."
    )
    parser.add_argument(
        "infile", nargs="?", default="probabilityECAi_inputs.bin",
        help="Inputs .bin file (default: probabilityECAi_inputs.bin)"
    )
    parser.add_argument(
        "--target", type=float, default=0.1,
        help="Target probability %% for split (default: 0.1)"
    )
    parser.add_argument(
        "--csv", type=str, default=None,
        help="Write results to CSV"
    )
    parser.add_argument(
        "--save", type=str, default=None,
        help="Save a tornado plot to this file"
    )
    parser.add_argument(
        "--category", type=str, default=None,
        help="Only analyze this category (default: all)"
    )
    args = parser.parse_args()

    data = load_inputs(args.infile)
    target = args.target / 100.0  # convert from percent to fraction

    # P field is index 8
    P_IDX = 8

    # Parameters to analyze (skip P itself and derived Q/mask_factor
    # which are deterministic from other params)
    ANALYZE_FIELDS = [
        0,   # PBR_sample
        1,   # lambda_bio
        2,   # gamma
        3,   # n_infected
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

    results = []
    categories = list(data.keys()) if args.category is None else [args.category]

    for cat in categories:
        if cat not in data:
            print(f"Warning: {cat} not in data, skipping")
            continue
        arr = data[cat]
        P = arr[:, P_IDX]
        exceed_mask = P > target
        nonexceed_mask = ~exceed_mask

        n_exceed = exceed_mask.sum()
        n_nonexceed = nonexceed_mask.sum()

        ecai_val = occupancy_params.get(cat, {}).get("ECAi", 0)

        print(f"\n{'='*70}")
        print(f"Category: {cat}  (ECAi={ecai_val} L/s/p)")
        print(f"  Exceed (P>{args.target}%):    {n_exceed:>7d} ({100*n_exceed/len(P):.2f}%)")
        print(f"  Non-exceed (P<={args.target}%): {n_nonexceed:>7d} ({100*n_nonexceed/len(P):.2f}%)")

        if n_exceed == 0 or n_nonexceed == 0:
            print(f"  Skipping (one group is empty)")
            continue

        print(f"\n  {'Parameter':<20s} {'Exceed mean':>14s} {'NonEx mean':>14s} {'Ratio':>10s} {'Log shift':>10s}")
        print(f"  {'-'*20} {'-'*14} {'-'*14} {'-'*10} {'-'*10}")

        cat_results = []
        for fi in ANALYZE_FIELDS:
            name = FIELD_NAMES[fi]
            vals = arr[:, fi]
            mean_exceed = vals[exceed_mask].mean()
            mean_nonexceed = vals[nonexceed_mask].mean()

            if mean_nonexceed > 0 and mean_exceed > 0:
                ratio = mean_exceed / mean_nonexceed
                log_shift = np.log10(ratio)
            else:
                ratio = float('nan')
                log_shift = float('nan')

            print(f"  {name:<20s} {mean_exceed:>14.4g} {mean_nonexceed:>14.4g} {ratio:>10.3f} {log_shift:>10.3f}")

            cat_results.append({
                "Category": cat,
                "Parameter": name,
                "Exceed_mean": mean_exceed,
                "NonExceed_mean": mean_nonexceed,
                "Ratio": ratio,
                "Log10_shift": log_shift,
                "N_exceed": n_exceed,
                "N_nonexceed": n_nonexceed,
            })

        results.extend(cat_results)

    # CSV output
    if args.csv:
        fieldnames = ["Category", "Parameter", "Exceed_mean",
                      "NonExceed_mean", "Ratio", "Log10_shift",
                      "N_exceed", "N_nonexceed"]
        with open(args.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in results:
                writer.writerow(row)
        print(f"\nResults saved to {args.csv}")

    # Tornado plot
    if args.save and results:
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        mpl.use('Agg')

        cats_with_results = list(set(r["Category"] for r in results))
        n_cats = len(cats_with_results)

        fig, axes = plt.subplots(n_cats, 1, figsize=(10, 3 * n_cats),
                                 squeeze=False)
        for ci, cat in enumerate(cats_with_results):
            ax = axes[ci, 0]
            cat_rows = [r for r in results if r["Category"] == cat]
            cat_rows.sort(key=lambda r: abs(r["Log10_shift"]),
                           reverse=True)
            params = [r["Parameter"] for r in cat_rows]
            shifts = [r["Log10_shift"] for r in cat_rows]
            y_pos = np.arange(len(params))
            colors = ['crimson' if s > 0 else 'steelblue' for s in shifts]
            ax.barh(y_pos, shifts, color=colors, alpha=0.7)
            ax.set_yticks(y_pos)
            ax.set_yticklabels(params, fontsize=8)
            ax.axvline(0, color='black', linewidth=0.5)
            ax.set_xlabel('log10(mean_exceed / mean_nonexceed)')
            ax.set_title(f'{cat}', fontsize=10)
            ax.invert_yaxis()

        fig.tight_layout(pad=2.0)
        fig.savefig(args.save, dpi=150)
        print(f"Tornado plot saved to {args.save}")


if __name__ == "__main__":
    main()