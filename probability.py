# simulate_infection_probability.py

import numpy as np
import csv
import sys
from model import sample_parameters, infection_probability, occupancy_params

def main():
    # Get N from command line argument (default to 10000)
    if len(sys.argv) > 1:
        try:
            N = int(sys.argv[1])
            if N < 1:
                raise ValueError
        except Exception:
            print("Invalid value for N, using default of 10000.")
            N = 10000
    else:
        N = 10000

    ECAi_values = list(range(5, 105, 5))  # Extend range if needed
    target_prob = 0.001  # 0.1%

    results = []

    print(f"\nSearching for ECAi (L/s/person) where mean probability < {target_prob*100:.2f}% for each category (N={N}):\n")
    print("| {:^20} | {:^28} |".format("Category", "Min ECAi (L/s/person)"))
    print("|" + "-"*22 + "|" + "-"*30 + "|")

    for category in occupancy_params:
        found = False
        for ECAi in ECAi_values:
            probabilities = []
            for _ in range(N):
                par = sample_parameters(category=category)
                prob = infection_probability(ECAi, par, category=category)
                probabilities.append(prob)

            mean_prob = np.mean(probabilities)
            if mean_prob < target_prob:
                print("| {:<20} | {:>28.2f} |".format(category, ECAi))
                results.append({
                    "Category": category,
                    "ECAi_Lps_for_P_lt_0.1pct": ECAi
                })
                found = True
                break
        if not found:
            print("| {:<20} | {:>28} |".format(category, "Not found"))
            results.append({
                "Category": category,
                "ECAi_Lps_for_P_lt_0.1pct": "Not found"
            })

    # Write to CSV
    csv_filename = "ecai_min_for_p_lt_0.1pct.csv"
    with open(csv_filename, mode="w", newline='') as file:
        writer = csv.DictWriter(file, fieldnames=["Category", "ECAi_Lps_for_P_lt_0.1pct"])
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"\nSummary table saved to {csv_filename}")

if __name__ == "__main__":
    main()
