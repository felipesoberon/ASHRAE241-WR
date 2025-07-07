# ecai_simulation.py

import numpy as np
import sys
import csv
from model import occupancy_params, sample_parameters, compute_ECAi

LPS_TO_CFM = 2.11888  # 1 L/s = 2.11888 CFM

if __name__ == "__main__":    

    # Read number of simulations from command line, default 10000
    if len(sys.argv) > 1:
        try:
            N = int(sys.argv[1])
            if N < 1:
                raise ValueError
        except Exception:
            print("Invalid argument for number of simulations. Using default N = 10000.")
            N = 10000
    else:
        N = 10000

    target_P = 0.001
    print("\nCalculating ECAi from fixed infection probability using", N, "samples for all categories...")
    print(f"\nTarget Probability: {target_P*100:.3f}%")

    print("\n| {:^15} | {:^33} | {:^33} | {:^33} |".format(
        "Category",
        "96th Percentile ECAi (L/s/person)",
        "Rounded (L/s/person, mult. of 5)",
        "Rounded (CFM/person, mult. of 10)"
    ))
    print("|{:-^17}|{:-^35}|{:-^35}|{:-^35}|".format('', '', '', ''))

    results = []

    for category in occupancy_params:
        ECAi_list = []
        for _ in range(N):
            par = sample_parameters(category=category)
            ECAi_val = compute_ECAi(par, target_P, category=category)
            ECAi_list.append(ECAi_val)

        ECAi_array = np.array(ECAi_list)
        percentile_96 = np.percentile(ECAi_array, 96)
        rounded_lps = int(np.ceil(percentile_96 / 5.0)) * 5
        percentile_96_CFM = percentile_96 * LPS_TO_CFM
        rounded_cfm = int(np.ceil(percentile_96_CFM / 10.0)) * 10

        print("| {:<15} | {:>33.2f} | {:>33} | {:>33} |".format(
            category,
            percentile_96,
            rounded_lps,
            rounded_cfm
        ))

        # Save results for CSV
        results.append({
            "Category": category,
            "Percentile_96_Lps": f"{percentile_96:.2f}",
            "Rounded_Lps": rounded_lps,
            "Rounded_CFM": rounded_cfm
        })

    # Write results to CSV
    csv_filename = "ecai_results.csv"
    with open(csv_filename, mode="w", newline='') as file:
        writer = csv.DictWriter(
            file, 
            fieldnames=["Category", "Percentile_96_Lps", "Rounded_Lps", "Rounded_CFM"]
        )
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"\nResults saved to {csv_filename}")
