# summary_table.py
#
# Comprehensive summary table from raw C++ binary data.
#
# One row per category with: ECAi, median, P75, P90, P96, P99,
# exceedance fraction >0.1%, and compliance group (A/B/C).
#
# Group A (compliant):      P96 < 0.1%
# Group B (tail risk):      Median < 0.1% but P96 >= 0.1%
# Group C (systematic):     Median >= 0.1%
#
# Usage:
#   python analysis/summary_table.py [infile.bin] [--target 0.1]
#                                    [--csv outfile.csv]
#
# Default infile: probabilityECAi_raw.bin

import argparse
import numpy as np
import csv
import os
import sys
from bin_reader import load_bin

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from model import occupancy_params


def classify_group(median_pct, p96_pct, target):
    if p96_pct < target:
        return "A"
    elif median_pct < target:
        return "B"
    else:
        return "C"


def main():
    parser = argparse.ArgumentParser(
        description="Comprehensive summary table from raw C++ "
                    "binary output."
    )
    parser.add_argument(
        "infile", nargs="?", default="probabilityECAi_raw.bin",
        help="Raw .bin file to analyze (default: probabilityECAi_raw.bin)"
    )
    parser.add_argument(
        "--target", type=float, default=0.1,
        help="Target probability %% for exceedance and grouping "
             "(default: 0.1)"
    )
    parser.add_argument(
        "--csv", type=str, default=None,
        help="Write results to this CSV file"
    )
    args = parser.parse_args()

    data = load_bin(args.infile)
    target = args.target

    # Header
    header = (
        f"{'Category':<20s} {'ECAi':>5s} {'Median':>9s} {'P75':>9s} "
        f"{'P90':>9s} {'P96':>9s} {'P99':>9s} "
        f"{'>0.1%':>8s} {'Group':>6s}"
    )
    print(header)
    print("-" * len(header))

    results = []

    for cat, arr in data.items():
        vals = arr.astype(np.float64) * 100.0  # to percent
        ecai_val = occupancy_params.get(cat, {}).get("ECAi", 0)
        median_pct = np.percentile(vals, 50)
        p75 = np.percentile(vals, 75)
        p90 = np.percentile(vals, 90)
        p96 = np.percentile(vals, 96)
        p99 = np.percentile(vals, 99)
        exceed = np.mean(vals > target) * 100.0
        group = classify_group(median_pct, p96, target)

        print(
            f"{cat:<20s} {ecai_val:>5.0f} {median_pct:>9.4f} {p75:>9.4f} "
            f"{p90:>9.4f} {p96:>9.4f} {p99:>9.4f} "
            f"{exceed:>7.2f}% {group:>6s}"
        )

        results.append({
            "Category": cat,
            "ECAi_Lps": ecai_val,
            "Median_pct": median_pct,
            "P75_pct": p75,
            "P90_pct": p90,
            "P96_pct": p96,
            "P99_pct": p99,
            "Exceed_0.1pct": exceed,
            "Group": group,
        })

    # Summary by group
    groups = {}
    for r in results:
        groups.setdefault(r["Group"], []).append(r["Category"])
    print()
    for g in sorted(groups):
        label = {
            "A": "A (compliant: P96 < 0.1%)",
            "B": "B (tail risk: median < 0.1%, P96 >= 0.1%)",
            "C": "C (systematic: median >= 0.1%)",
        }[g]
        print(f"  Group {label}: {len(groups[g])} categories")
        print(f"    {', '.join(groups[g])}")

    if args.csv:
        fieldnames = [
            "Category", "ECAi_Lps", "Median_pct", "P75_pct",
            "P90_pct", "P96_pct", "P99_pct", "Exceed_0.1pct",
            "Group",
        ]
        with open(args.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in results:
                writer.writerow(row)
        print(f"\nResults saved to {args.csv}")


if __name__ == "__main__":
    main()