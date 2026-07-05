# percentile_heatmap.py
#
# Percentile heatmap from raw C++ binary data.
#
# A 25-row x N-column grid where rows are categories (sorted by
# 96th percentile, worst at bottom) and columns are percentiles
# (50, 60, 65, 70, 75, 80, 85, 90, 95, 96, 99). Color = infection
# probability (%). A contour or boundary line marks the 0.1% target
# so it is immediately visible where the standard achieves its
# target and where it does not.
#
# Usage:
#   python analysis/percentile_heatmap.py [infile.bin]
#       [--target 0.1] [--save FILE] [--show]
#
# Default infile: probabilityECAi_raw.bin

import argparse
import numpy as np
import matplotlib.pyplot as plt
from bin_reader import load_bin


def main():
    parser = argparse.ArgumentParser(
        description="Percentile heatmap of raw C++ binary output."
    )
    parser.add_argument(
        "infile", nargs="?", default="probabilityECAi_raw.bin",
        help="Raw .bin file to analyze (default: probabilityECAi_raw.bin)"
    )
    parser.add_argument(
        "--target", type=float, default=0.1,
        help="Target probability %% drawn as a contour line "
             "(default: 0.1)"
    )
    parser.add_argument(
        "--percentiles", type=int, nargs="+",
        default=[50, 60, 65, 70, 75, 80, 85, 90, 95, 96, 99],
        help="Percentiles for the columns (default: 50 60 65 70 75 "
             "80 85 90 95 96 99)"
    )
    parser.add_argument(
        "--save", type=str, default=None,
        help="Save the figure to this file instead of displaying"
    )
    parser.add_argument(
        "--show", action="store_true",
        help="Display the figure interactively (in addition to --save)"
    )
    args = parser.parse_args()

    data = load_bin(args.infile)
    percentiles = args.percentiles

    # Sort categories by 96th percentile (most compliant at top,
    # least at bottom)
    cats_sorted = sorted(data.keys(),
                         key=lambda c: np.percentile(data[c], 96))

    # Build the heatmap matrix
    n_cats = len(cats_sorted)
    n_pcts = len(percentiles)
    matrix = np.zeros((n_cats, n_pcts))
    for i, cat in enumerate(cats_sorted):
        vals = data[cat].astype(np.float64) * 100.0  # to percent
        for j, p in enumerate(percentiles):
            matrix[i, j] = np.percentile(vals, p)

    fig, ax = plt.subplots(figsize=(10, 10))

    # Log-scaled color map: values span orders of magnitude
    # Use LogNorm with a floor to avoid log(0)
    from matplotlib.colors import LogNorm
    vmin = max(matrix.min(), 0.001)  # floor at 0.001%
    vmax = max(matrix.max(), vmin * 10)
    im = ax.pcolormesh(matrix.T, cmap="RdYlGn_r",
                       norm=LogNorm(vmin=vmin, vmax=vmax),
                       shading="auto")

    # Contour at the target -- uses the same matrix.T orientation
    # as pcolormesh so the contour aligns with the cells
    cs = ax.contour(matrix.T, levels=[args.target],
                    colors="black", linewidths=2.0)
    ax.clabel(cs, fmt="%g%% target", fontsize=8,
              inline=True)

    # Annotate each cell with the value
    for i in range(n_cats):
        for j in range(n_pcts):
            val = matrix[i, j]
            # Choose text color based on cell color intensity:
            # dark red cells get white text, light green cells get black
            log_val = np.log10(max(val, vmin))
            log_mid = 0.5 * (np.log10(vmin) + np.log10(vmax))
            text_color = "white" if log_val > log_mid else "black"
            ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                    fontsize=6, color=text_color)

    ax.set_yticks(np.arange(n_cats))
    ax.set_yticklabels(cats_sorted, fontsize=9)
    ax.set_xticks(np.arange(n_pcts))
    pct_labels = [str(p) for p in percentiles]
    ax.set_xticklabels(pct_labels, fontsize=9)
    ax.set_xlabel("Percentile")
    ax.set_ylabel("Category (sorted by P96, best at top)")
    ax.set_title("Infection probability (%) by category and "
                 "percentile\nASHRAE 241 ECAi values, "
                 "require_infectors=True")

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("Infection Probability (%)  [log scale]")

    fig.tight_layout(pad=2.0)

    if args.save:
        fig.savefig(args.save, dpi=150)
        print(f"Figure saved to {args.save}")
    if args.show or not args.save:
        plt.show()


if __name__ == "__main__":
    main()