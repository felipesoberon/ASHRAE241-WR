# simulate_infection_probability.py

import sys
import numpy as np
from model import sample_parameters, infection_probability, occupancy_params

def main():
    if len(sys.argv) < 2:
        print("Usage: python simulate_infection_probability.py <Category>")
        print("Example categories:", ', '.join(list(occupancy_params.keys())[:5]), "...")
        sys.exit(1)

    category = sys.argv[1]
    if category not in occupancy_params:
        print(f"Error: '{category}' not in defined categories.")
        print("Available categories:", ', '.join(occupancy_params.keys()))
        sys.exit(1)

    N = 10000
    ECAi_values = list(range(5, 65, 5))  # 5, 10, ..., 

    # Print table header
    print(f"\nResults for category '{category}' over {N} simulations for each ECAi (L/s/person):\n")
    print("| {:^8} | {:^8} | {:^8} | {:^8} | {:^8} | {:^8} |".format(
        "ECAi", "Min %", "Max %", "Median %", "Mean %", "96th %"))
    print("|" + "-"*10 + "|" + "-"*10 + "|" + "-"*10 + "|" + "-"*10 + "|" + "-"*10 + "|" + "-"*10 + "|")

    for ECAi in ECAi_values:
        probabilities = []
        for _ in range(N):
            par = sample_parameters(category=category)
            prob = infection_probability(ECAi, par, category=category)
            probabilities.append(prob)

        probabilities = np.array(probabilities)
        min_val = np.min(probabilities) * 100
        max_val = np.max(probabilities) * 100
        median = np.median(probabilities) * 100
        mean = np.mean(probabilities) * 100
        perc_96 = np.percentile(probabilities, 96) * 100

        print("| {:>8.2f} | {:>8.2f} | {:>8.2f} | {:>8.2f} | {:>8.2f} | {:>8.2f} |".format(
            ECAi, min_val, max_val, median, mean, perc_96))

if __name__ == "__main__":
    main()
