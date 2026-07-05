# boxplot_bin.py
#
# Horizontal box-and-whisker plot (log-scaled probability axis) from
# raw C++ binary data. One box per occupancy category, categories on
# the y-axis and infection probability (%) on a log x-axis. Boxes show
# quartiles with a black median line and 1.5x-IQR whiskers (outlier
# fliers suppressed for readability). Each category's 96th percentile
# is marked with a red tick, and a dashed line marks the target
# probability (default 0.1%).
#
# This mirrors analysis/boxplot.py but reads the C++ .bin format
# directly instead of .npz.
#
# Usage:
#   python analysis/boxplot_bin.py [infile.bin] [--target 0.1]
#                                 [--xmin 0.001] [--xmax 100]
#                                 [--save FILE] [--show]
#
# Default infile: probabilityECAi_raw.bin

import argparse
import numpy as np
import matplotlib.pyplot as plt
from bin_reader import load_bin


def main():
    parser = argparse.ArgumentParser(
        description="Horizontal box-and-whisker plot (log scale) of "
                    "raw C++ binary output."
    )
    parser.add_argument(
        "infile", nargs="?", default="probabilityECAi_raw.bin",
        help="Raw .bin file to plot (default: probabilityECAi_raw.bin)"
    )
    parser.add_argument(
        "--target", type=float, default=0.1,
        help="Target probability %% drawn as a reference line "
             "(default: 0.1)"
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
    categories = list(data.keys())

    # Convert to percent (drop non-positive values for the log axis)
    # and record each category's 96th percentile.
    data_pct = []
    p96 = []
    for c in categories:
        vals = data[c].astype(np.float64) * 100.0
        vals = vals[vals > 0]
        data_pct.append(vals)
        p96.append(np.percentile(vals, 96))

    positions = np.arange(len(categories))

    fig, ax = plt.subplots(figsize=(10, 12))
    bp = ax.boxplot(
        data_pct, positions=positions, orientation="horizontal",
        widths=0.6, patch_artist=True, showfliers=False,
        medianprops=dict(color="black"))
    for box in bp["boxes"]:
        box.set_facecolor("steelblue")
        box.set_alpha(0.6)

    # 96th percentile markers (the ASHRAE 241 metric)
    ax.scatter(p96, positions, color="crimson", marker="|", s=200,
               zorder=3, label="96th percentile")

    # Target reference line
    ax.axvline(args.target, color="red", linestyle="--", linewidth=1,
              label=f"{args.target}% target")

    # Log axis clipped to [--xmin, --xmax], labelled with plain values
    ax.set_xscale("log")
    lo = int(np.floor(np.log10(args.xmin)))
    hi = int(np.ceil(np.log10(args.xmax)))
    ticks = [10.0 ** k for k in range(lo, hi + 1)]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.set_xlim(args.xmin, args.xmax)

    ax.set_yticks(positions)
    ax.set_yticklabels(categories)
    ax.invert_yaxis()
    ax.set_xlabel("Infection Probability (%)  [log scale]")
    ax.set_title("Infection probability distribution by occupancy "
                 "(ASHRAE 241 ECAi)")
    ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.7)
    ax.legend(loc="lower right")
    fig.tight_layout(pad=3.0)
    fig.subplots_adjust(left=0.16, right=0.94, top=0.94, bottom=0.08)

    if args.save:
        fig.savefig(args.save, dpi=150)
        print(f"Figure saved to {args.save}")
    if args.show or not args.save:
        plt.show()


if __name__ == "__main__":
    main()