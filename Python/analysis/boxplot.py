# analysis/boxplot.py
#
# Horizontal box-and-whisker plot (log-scaled probability axis) of the raw
# infection-probability data saved by probabilityECAi.py --save-all.
#
# One box per occupancy category, categories on the y-axis and infection
# probability (%) on a log x-axis. Boxes show the quartiles (median line in
# black) with 1.5x-IQR whiskers; outlier fliers are suppressed to keep the
# plot readable. A dashed line marks the target probability (default 0.1%)
# and a red tick marks each category's 96th percentile (the ASHRAE 241
# metric of interest).

import argparse
import numpy as np
import matplotlib.pyplot as plt


def load(path):
    d = np.load(path)
    data = {c: d[c] for c in d.files}
    d.close()  # release the file handle (Windows keeps it open lazily)
    return data


def main():
    parser = argparse.ArgumentParser(
        description="Horizontal box-and-whisker plot (log scale) of raw probabilityECAi output (.npz)."
    )
    parser.add_argument("infile", nargs="?", default="probabilityECAi_raw.npz",
                        help="Raw .npz file to plot (default probabilityECAi_raw.npz)")
    parser.add_argument("--target", type=float, default=0.1,
                        help="Target probability (%%) drawn as a reference line (default 0.1)")
    parser.add_argument("--xmin", type=float, default=0.001,
                        help="Left edge of the log x-axis in %% (default 0.001)")
    parser.add_argument("--xmax", type=float, default=100,
                        help="Right edge of the log x-axis in %% (default 100)")
    parser.add_argument("--save", type=str, default=None,
                        help="Save the figure to this file instead of showing it")
    args = parser.parse_args()

    data = load(args.infile)
    categories = list(data.keys())

    # Convert to percent (drop any non-positive values for the log axis) and
    # record each category's 96th percentile.
    data_pct = []
    p96 = []
    for c in categories:
        vals = data[c].astype(np.float64) * 100.0
        vals = vals[vals > 0]
        data_pct.append(vals)
        p96.append(np.percentile(vals, 96))

    positions = np.arange(len(categories))

    fig, ax = plt.subplots(figsize=(10, 12))
    bp = ax.boxplot(data_pct, positions=positions, vert=False, widths=0.6,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color="black"))
    for box in bp["boxes"]:
        box.set_facecolor("steelblue")
        box.set_alpha(0.6)

    # 96th percentile markers (the ASHRAE 241 metric).
    ax.scatter(p96, positions, color="crimson", marker="|", s=200,
               zorder=3, label="96th percentile")

    # Target reference line.
    ax.axvline(args.target, color="red", linestyle="--", linewidth=1,
               label=f"{args.target}% target")

    # Log axis clipped to [--xmin, --xmax], labelled with plain values.
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
    ax.set_title("Infection probability distribution by occupancy (ASHRAE 241 ECAi)")
    ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.7)
    ax.legend(loc="lower right")
    fig.tight_layout(pad=3.0)
    fig.subplots_adjust(left=0.16, right=0.94, top=0.94, bottom=0.08)

    if args.save:
        fig.savefig(args.save, dpi=150)
        print(f"Figure saved to {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
