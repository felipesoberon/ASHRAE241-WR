# correlation_analysis.py
#
# Driver analysis: compute Spearman rank correlation between each
# input parameter and the infection probability P.
#
# For each category, computes the Spearman correlation of every input
# parameter with P. High-correlation parameters (positive) are the
# primary drivers of the high-probability tail.
#
# Optionally produces a heatmap of correlations (categories x
# parameters).
#
# Usage:
#   python analysis/correlation_analysis.py [inputs.bin]
#       [--csv outfile.csv] [--save FIG]
#
# Default inputs file: probabilityECAi_inputs.bin

import argparse
import numpy as np
import csv
import os
import sys
from scipy.stats import spearmanr
from inputs_reader import load_inputs, FIELD_NAMES

# Add parent dir for model.py (ECAi values)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from model import occupancy_params


def main():
    parser = argparse.ArgumentParser(
        description="Spearman correlation between input parameters "
                    "and infection probability."
    )
    parser.add_argument(
        "infile", nargs="?", default="probabilityECAi_inputs.bin",
        help="Inputs .bin file (default: probabilityECAi_inputs.bin)"
    )
    parser.add_argument(
        "--csv", type=str, default=None,
        help="Write correlation table to CSV"
    )
    parser.add_argument(
        "--save", type=str, default=None,
        help="Save a correlation heatmap to this file"
    )
    args = parser.parse_args()

    data = load_inputs(args.infile)

    # P field is index 8
    P_IDX = 8

    # Parameters to correlate with P (skip P itself, Q which is
    # deterministic from P, mask_factor which is constant per category,
    # n_infected which is discrete)
    CORR_FIELDS = [
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
    corr_names = [FIELD_NAMES[fi] for fi in CORR_FIELDS]

    results = []
    all_corrs = {}  # for heatmap

    for cat, arr in data.items():
        P = arr[:, P_IDX]
        ecai_val = occupancy_params.get(cat, {}).get("ECAi", 0)

        print(f"\n{'='*60}")
        print(f"Category: {cat}  (ECAi={ecai_val} L/s/p)")

        corrs = []
        for fi in CORR_FIELDS:
            name = FIELD_NAMES[fi]
            vals = arr[:, fi]
            rho, pval = spearmanr(vals, P)
            if np.isnan(rho):
                rho = 0.0
            corrs.append(rho)
            print(f"  {name:<20s}  rho={rho:+.4f}  (p={pval:.2e})")
            results.append({
                "Category": cat,
                "Parameter": name,
                "Spearman_rho": rho,
                "p_value": pval,
            })

        all_corrs[cat] = corrs

    # CSV output
    if args.csv:
        fieldnames = ["Category", "Parameter", "Spearman_rho", "p_value"]
        with open(args.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in results:
                writer.writerow(row)
        print(f"\nResults saved to {args.csv}")

    # Heatmap
    if args.save:
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        mpl.use('Agg')

        cats = list(all_corrs.keys())
        matrix = np.array([all_corrs[c] for c in cats])

        fig, ax = plt.subplots(figsize=(10, 10))
        im = ax.imshow(matrix, cmap='RdBu_r', vmin=-1, vmax=1,
                       aspect='auto')
        ax.set_xticks(np.arange(len(corr_names)))
        ax.set_xticklabels(corr_names, rotation=45, ha='right',
                           fontsize=8)
        ax.set_yticks(np.arange(len(cats)))
        ax.set_yticklabels(cats, fontsize=9)
        ax.set_title("Spearman correlation: input parameters vs "
                     "infection probability P")
        fig.colorbar(im, ax=ax, shrink=0.6, label="Spearman rho")

        # Annotate cells
        for i in range(len(cats)):
            for j in range(len(corr_names)):
                val = matrix[i, j]
                color = "white" if abs(val) > 0.5 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=6, color=color)

        fig.tight_layout(pad=2.0)
        fig.savefig(args.save, dpi=150)
        print(f"Correlation heatmap saved to {args.save}")


if __name__ == "__main__":
    main()