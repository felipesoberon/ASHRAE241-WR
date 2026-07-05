# exceedance.py
#
# Compute exceedance fractions from raw C++ binary data.
#
# For each category, reports:
#   - Mean and median infection probability (%)
#   - 96th percentile (%)
#   - Fraction of simulations exceeding 0.1%, 0.5%, 1.0% targets
#
# Usage:
#   python analysis/exceedance.py [infile.bin] [--targets 0.1 0.5 1.0]
#                                [--csv outfile.csv]
#
# Default infile: probabilityECAi_raw.bin
# Default targets: 0.1 0.5 1.0 (percent)

import argparse
import numpy as np
import csv
import os
import sys
from bin_reader import load_bin

# Import occupancy params for ECAi values (add parent dir to path)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from model import occupancy_params


def main():
    parser = argparse.ArgumentParser(
        description="Exceedance fraction analysis of raw C++ binary output."
    )
    parser.add_argument(
        "infile", nargs="?", default="probabilityECAi_raw.bin",
        help="Raw .bin file to analyze (default: probabilityECAi_raw.bin)"
    )
    parser.add_argument(
        "--targets", type=float, nargs="+", default=[0.1, 0.5, 1.0],
        help="Target probabilities in %% for exceedance (default: 0.1 0.5 1.0)"
    )
    parser.add_argument(
        "--csv", type=str, default=None,
        help="Write results to this CSV file (default: stdout only)"
    )
    args = parser.parse_args()

    data = load_bin(args.infile)
    targets = args.targets  # in percent

    # Header
    n_targets = len(targets)
    target_cols = [f">{t}%" for t in targets]

    # Print table header
    header = f"{'Category':<20s} {'ECAi':>6s} {'Mean':>10s} {'Median':>10s} {'P96':>10s}"
    for tc in target_cols:
        header += f" {tc:>10s}"
    print(header)
    print("-" * len(header))

    results = []

    # Category order: iterate by insertion order (the .bin preserves it,
    # and Python dict since 3.7 keeps it). To match the model's category_order
    # exactly, use the key order from the file.
    for cat, arr in data.items():
        vals = arr.astype(np.float64) * 100.0  # convert to percent

        ecai_val = occupancy_params.get(cat, {}).get("ECAi", 0)
        mean_pct = np.mean(vals)
        median_pct = np.median(vals)
        p96 = np.percentile(vals, 96)

        exceed = []
        for t in targets:
            frac = np.mean(vals > t) * 100.0  # fraction exceeding, in percent
            exceed.append(frac)

        row_str = (
            f"{cat:<20s} {ecai_val:>6.1f} {mean_pct:>10.4f} {median_pct:>10.4f} "
            f"{p96:>10.4f}"
        )
        for e in exceed:
            row_str += f" {e:>10.3f}"
        print(row_str)

        results.append({
            "Category": cat,
            "ECAi": ecai_val,
            "Mean_pct": mean_pct,
            "Median_pct": median_pct,
            "P96_pct": p96,
            **{
                f"Exceed_{t}pct": e
                for t, e in zip(targets, exceed)
            }
        })

    # CSV output
    if args.csv:
        fieldnames = ["Category", "ECAi", "Mean_pct", "Median_pct", "P96_pct"]
        fieldnames += [f"Exceed_{t}pct" for t in targets]
        with open(args.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in results:
                writer.writerow(row)
        print(f"\nResults saved to {args.csv}")


if __name__ == "__main__":
    main()