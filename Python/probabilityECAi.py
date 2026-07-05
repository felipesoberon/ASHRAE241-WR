# simulate_infection_probability_ecai.py

import numpy as np
import csv
import argparse
from model import sample_parameters, infection_probability, occupancy_params
from randomManager import RandomNumberManager

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate the 96th percentile infection probability at the ASHRAE 241 ECAi for each category."
    )
    parser.add_argument("N", nargs="?", type=int, default=10000,
                        help="Number of simulations per category (default 10000)")
    parser.add_argument("--save-all", action="store_true",
                        help="Save every raw simulation probability to a compressed .npz file")
    parser.add_argument("--outfile", type=str, default="probabilityECAi_raw.npz",
                        help="Output file for raw data when --save-all is set (default probabilityECAi_raw.npz)")
    args = parser.parse_args()

    N = args.N
    if N < 1:
        print("Invalid value for N, using default of 10000.")
        N = 10000

    target_prob = 0.001  # 0.1%
    results = []
    raw_data = {} if args.save_all else None

    print(f"\nEvaluating the 96th percentile probability at ECAi from occupancy_params for each category (N={N}):\n")
    print("| {:^20} | {:^15} | {:^20} |".format("Category", "ECAi (L/s/p)", "96th per P (%)"))
    print("|" + "-"*22 + "|" + "-"*17 + "|" + "-"*22 + "|")

    rng = RandomNumberManager()

    for category, params in occupancy_params.items():
        ECAi = params.get("ECAi")
        if ECAi is None:
            print(f"| {category:<20} | {'No ECAi':>15} | {'Skipped':>20} |")
            results.append({
                "Category": category,
                "ECAi_Lps": "No ECAi",
                "P_96th_percentile": "Skipped"
            })
            continue

        probabilities = []
        for _ in range(N):
            par = sample_parameters(rng, category=category)
            prob, _ = infection_probability(ECAi, par, rng, category=category, require_infectors=True)
            probabilities.append(prob)

        percentile_96 = np.percentile(probabilities, 96) * 100  # As percent
        if raw_data is not None:
            raw_data[category] = np.asarray(probabilities, dtype=np.float32)
        print("| {:<20} | {:>15.2f} | {:>20.3f} |".format(category, ECAi, percentile_96))
        results.append({
            "Category": category,
            "ECAi_Lps": ECAi,
            "P_96th_percentile": percentile_96
        })

    # Write to CSV
    csv_filename = "ecai_ashrae241_96th_percentile.csv"
    with open(csv_filename, mode="w", newline='') as file:
        writer = csv.DictWriter(file, fieldnames=["Category", "ECAi_Lps", "P_96th_percentile"])
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"\nSummary table saved to {csv_filename}")

    # Save all raw simulation probabilities if requested
    if raw_data is not None:
        np.savez_compressed(args.outfile, **raw_data)
        total_values = sum(arr.size for arr in raw_data.values())
        print(f"Raw data ({total_values} values, raw probability 0-1) saved to {args.outfile}")

if __name__ == "__main__":
    main()
