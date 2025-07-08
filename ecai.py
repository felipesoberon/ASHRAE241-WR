# ecai_simulation.py

import numpy as np
import sys
import csv
from model import occupancy_params, sample_parameters, compute_ECAi
from randomManager import RandomNumberManager

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
    print(f"\nTarget Probability: {target_P*100:.3f}%\n")

    print("Round (L/s/p) is the ceiling to nearest multiple of 5.")
    print("Round (CFM/p) is the ceiling to nearest multiple of 10.\n")

    header = (
        "| {:<15} | {:>13} | {:>13} | {:>15} | {:>13} |"
        .format("Category", "ECAi (L/s/p)", "Round (L/s/p)", "Round (CFM/p)", "Zero Inf. (%)")
    )
    line = "+-----------------+---------------+---------------+-----------------+---------------+"

    print(line)
    print(header)
    print(line)

    results = []
    grand_zero_infected = 0
    rng = RandomNumberManager()  # LHS Random number manager instance

    for category in occupancy_params:
        ECAi_list = []
        zero_infected_count = 0

        for _ in range(N):
            par = sample_parameters(rng, category=category)
            ECAi_val, infected_flag = compute_ECAi(par, target_P, rng, category=category)
            ECAi_list.append(ECAi_val)
            if infected_flag == 0:
                zero_infected_count += 1                
                
        grand_zero_infected += zero_infected_count

        ECAi_array = np.array(ECAi_list)
        percentile_96 = np.percentile(ECAi_array, 96)
        rounded_lps = int(np.ceil(percentile_96 / 5.0)) * 5
        percentile_96_CFM = percentile_96 * LPS_TO_CFM
        rounded_cfm = int(np.ceil(percentile_96_CFM / 10.0)) * 10
        zero_infected_percent = 100 * zero_infected_count / N

        print("| {:<15} | {:13.2f} | {:13} | {:15} | {:12.1f}  |".format(
            category,
            percentile_96,
            rounded_lps,
            rounded_cfm,
            zero_infected_percent
        ))

        results.append({
            "Category": category,
            "Percentile_96_Lps": f"{percentile_96:.2f}",
            "Rounded_Lps": rounded_lps,
            "Rounded_CFM": rounded_cfm,
            "ZeroInfectedPercent": zero_infected_percent
        })

    print(line)
    
    # Calculate grand total
    grand_total_zero_infected = sum(int(round(row["ZeroInfectedPercent"] * N / 100)) for row in results)
    grand_total_simulations = len(occupancy_params) * N
    grand_percent = 100 * grand_total_zero_infected / grand_total_simulations

    print("\nGrand Total Simulations with Zero Infected: {} of {} ({:.1f} %)".format(
        grand_total_zero_infected,
        grand_total_simulations,
        grand_percent
    ))


    # Write results to CSV
    csv_filename = "ecai_results.csv"
    with open(csv_filename, mode="w", newline='') as file:
        writer = csv.DictWriter(
            file, 
            fieldnames=[
                "Category",
                "Percentile_96_Lps",
                "Rounded_Lps",
                "Rounded_CFM",
                "ZeroInfectedPercent"
            ]
        )
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"\nResults saved to {csv_filename}")


