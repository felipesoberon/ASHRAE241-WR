# ccdf_plot.py
#
# Complementary CDF (exceedance curve) plot from raw C++ binary data.
#
# For each category, plots the fraction of simulations exceeding a
# given infection probability threshold. The 0.1% target line is
# drawn vertically so the reader can see exactly what fraction of
# scenarios exceed the ASHRAE 241 target per category.
#
# Usage:
#   python analysis/ccdf_plot.py [infile.bin] [--target 0.1]
#                                [--xmin 0.001] [--xmax 100]
#                                [--save FILE] [--show]
#
# Default infile: probabilityECAi_raw.bin
# All axis values are in percent (0.001% to 100%).

import argparse
import numpy as np
import matplotlib.pyplot as plt
from bin_reader import load_bin


def main():
    parser = argparse.ArgumentParser(
        description="Complementary CDF (exceedance) plot of raw "
                    "C++ binary output."
    )
    parser.add_argument(
        "infile", nargs="?", default="probabilityECAi_raw.bin",
        help="Raw .bin file to plot (default: probabilityECAi_raw.bin)"
    )
    parser.add_argument(
        "--target", type=float, default=0.1,
        help="Target probability %% drawn as a vertical reference "
             "line (default: 0.1)"
    )
    parser.add_argument(
        "--xmin", type=float, default=0.001,
        help="Left edge of the log x-axis in %% (default: 0.001)"
    )
    parser.add_argument(
        "--xmax", type=float, default=100.0,
        help="Right edge of the log x-axis in %% (default: 100)"
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

    fig, ax = plt.subplots(figsize=(10, 8))

    # Colormap: one color per category, sorted by 96th percentile
    # (worst at top of legend)
    cats_sorted = sorted(data.keys(),
                         key=lambda c: np.percentile(data[c], 96))
    n = len(cats_sorted)
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, n))

    for i, cat in enumerate(cats_sorted):
        vals = data[cat].astype(np.float64) * 100.0  # to percent
        vals = np.sort(vals)
        # CCDF: fraction of simulations with P >= threshold
        ccdf = 1.0 - np.arange(1, len(vals) + 1) / len(vals)
        # Plot only positive values (log x-axis)
        mask = vals > 0
        ax.plot(vals[mask], ccdf[mask] * 100.0, color=colors[i],
                linewidth=0.8, alpha=0.8, label=cat)

    # Target vertical line
    ax.axvline(args.target, color="red", linestyle="--", linewidth=1.5,
              label=f"{args.target}% target")

    ax.set_xscale("log")
    ax.set_xlim(args.xmin, args.xmax)
    lo = int(np.floor(np.log10(args.xmin)))
    hi = int(np.ceil(np.log10(args.xmax)))
    ticks = [10.0 ** k for k in range(lo, hi + 1)]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())

    ax.set_xlabel("Infection Probability (%)  [log scale]")
    ax.set_ylabel("Fraction of simulations exceeding (%)")
    ax.set_title("Exceedance curves: infection probability at "
                 "ASHRAE 241 ECAi values")
    ax.grid(True, which="both", axis="both", linestyle="--",
            linewidth=0.5, alpha=0.5)
    ax.legend(loc="upper right", fontsize=7, ncol=2)
    fig.tight_layout(pad=3.0)

    if args.save:
        fig.savefig(args.save, dpi=150)
        print(f"Figure saved to {args.save}")
    if args.show or not args.save:
        plt.show()


if __name__ == "__main__":
    main()