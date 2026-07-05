# analysis/percentiles.py
#
# Post-processes the raw infection-probability data saved by
# probabilityECAi.py --save-all (a compressed .npz with one float32
# array of raw probabilities [0-1] per occupancy category).
#
# Two reports:
#   1. table    - per-category percentile table (default 50, 75, 96).
#   2. threshold - highest percentile (swept 50..95 by 5) at which ALL
#                  categories stay below a target probability (default 0.1%).

import argparse
import numpy as np


def load(path):
    d = np.load(path)
    data = {c: d[c] for c in d.files}
    d.close()  # release the file handle (Windows keeps it open lazily)
    return data


def report_table(data, percentiles):
    cols = "".join(f"{f'{p}th (%)':>12}" for p in percentiles)
    print(f"{'Category':<15}{cols}")
    print("-" * (15 + 12 * len(percentiles)))
    for cat, arr in data.items():
        vals = np.percentile(arr, percentiles) * 100
        print(f"{cat:<15}" + "".join(f"{v:>12.3f}" for v in vals))


def worst_at(data, p, cache):
    """Worst-case (max across categories) p-th percentile probability, in %."""
    if p not in cache:
        vals = {c: np.percentile(a, p) * 100 for c, a in data.items()}
        worst_cat = max(vals, key=vals.get)
        cache[p] = (vals[worst_cat], worst_cat)
    return cache[p]


def report_threshold(data, target_pct):
    """Coarse-to-fine search for the highest percentile at which every
    category stays below target_pct. Scans 25..95 in steps of 20, then
    refines within the bracketing interval at steps 10, 5, and 1."""
    cache = {}
    lo, hi = 25, 95
    best = None
    print(f"Coarse-to-fine search for highest percentile with all categories < {target_pct}%\n")

    for step in (20, 10, 5, 1):
        ps = list(range(lo, hi + 1, step))
        if ps[-1] != hi:
            ps.append(hi)

        print(f"Scan step {step} ({lo}..{hi}):")
        print(f"{'Percentile':>10} {'Max P (%)':>12} {'Worst category':>18} {'All < target?':>14}")
        print("-" * 58)
        results = []
        for p in ps:
            worst, worst_cat = worst_at(data, p, cache)
            ok = worst < target_pct
            results.append((p, ok))
            print(f"{p:>10} {worst:>12.4f} {worst_cat:>18} {str(ok):>14}")
        print()

        # Find the adjacent pair bracketing the crossing (ok -> not ok)
        interval = None
        for (a, oka), (b, okb) in zip(results, results[1:]):
            if oka and not okb:
                interval = (a, b)
                break

        if interval is None:
            if all(ok for _, ok in results):
                best = ps[-1]  # crossing is at or above the top of the range
                print(f"All categories stay below {target_pct}% through the {best}th percentile "
                      f"(top of the {lo}..{hi} range); crossing not reached.")
            else:
                best = None
                print(f"Even the {lo}th percentile exceeds {target_pct}%; crossing is below {lo}.")
            return

        lo, hi = interval
        best = lo

    print(f"Highest percentile where ALL categories < {target_pct}%: {best}")


def main():
    parser = argparse.ArgumentParser(
        description="Percentile analysis of raw probabilityECAi output (.npz)."
    )
    parser.add_argument("infile", nargs="?", default="probabilityECAi_raw.npz",
                        help="Raw .npz file to analyze (default probabilityECAi_raw.npz)")
    parser.add_argument("--report", choices=["table", "threshold", "both"],
                        default="both", help="Which report to produce (default both)")
    parser.add_argument("--percentiles", type=int, nargs="+", default=[50, 75, 96],
                        help="Percentiles for the table report (default 50 75 96)")
    parser.add_argument("--target", type=float, default=0.1,
                        help="Target probability in %% for the threshold report (default 0.1)")
    args = parser.parse_args()

    data = load(args.infile)

    if args.report in ("table", "both"):
        report_table(data, args.percentiles)
    if args.report == "both":
        print()
    if args.report in ("threshold", "both"):
        report_threshold(data, args.target)


if __name__ == "__main__":
    main()
