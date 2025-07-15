# singleProbability.py

import numpy as np
import matplotlib.pyplot as plt
from model import sample_parameters, infection_probability, occupancy_params
from randomManager import RandomNumberManager


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description="Run single occupancy category probability simulation.")
    parser.add_argument('--N', type=int, default=10000, help="Number of simulations (default 10000)")
    parser.add_argument('--category', type=str, default='Classroom', choices=occupancy_params.keys(), help="Occupancy category")
    parser.add_argument('--community_rate', type=float, default=None, help="Community infection rate (0–1), optional")
    parser.add_argument('--no_zero_infectors', action='store_false', dest='allow_zero_infectors', help="Disallow zero infected simulations (default allows zero infectors)")
    parser.add_argument('--ecai', type=float, default=None, help="ECAi in L/s/person (default from occupancy params)")
    parser.add_argument('--show_plots', action='store_true', help="Show plots of probability distribution")

    return parser.parse_args()


def print_progress(iteration, total):
    bar_length = 100
    progress = int((iteration / total) * bar_length)
    bar = '*' * progress + ' ' * (bar_length - progress)
    print(f"\rProgress: |{bar}| {iteration}/{total}", end='', flush=True)


def main():
    args = parse_args()

    N = args.N
    category = args.category
    community_rate = args.community_rate
    allow_zero_infectors = args.allow_zero_infectors
    ECAi = args.ecai if args.ecai is not None else occupancy_params[category]["ECAi"]

    rng = RandomNumberManager()
    probabilities = []
    zero_infectors_count = 0

    for i in range(1, N + 1):
        par = sample_parameters(rng, category=category)
        comm_rate = 0 if community_rate is None else community_rate
        prob, infected_flag = infection_probability(
            ECAi,
            par,
            rng,
            category=category,
            require_infectors=not allow_zero_infectors,
            override_community_rate=comm_rate
        )
        probabilities.append(prob)
        if infected_flag == 0:
            zero_infectors_count += 1

        # Update progress
        if N >= 100 and i % (N // 100) == 0:
            print_progress(i, N)

    print_progress(N, N)
    print("\n")  # Newline after progress bar

    probabilities = np.array(probabilities) * 100  # Convert to %
    p96 = np.percentile(probabilities, 96)
    zero_percent = (zero_infectors_count / N) * 100

    print(f"\nSimulation results for category '{category}':")
    print(f"  Simulations: {N}")
    print(f"  ECAi: {ECAi} L/s/person")
    if community_rate is None:
        default_rate = occupancy_params[category]["community_rate"]
        print(f"  Community infection rate: Default ({default_rate})")
    else:
        print(f"  Community infection rate: {community_rate}")
    print(f"  Allow zero infectors: {allow_zero_infectors}")
    print(f"  96th percentile of infection probability: {p96:.3f}%")
    print(f"  Percentage of simulations with zero infectors: {zero_percent:.1f}%\n")

    if args.show_plots:
        fig, axs = plt.subplots(2, 1, figsize=(10, 10))

        # Plot: Linear scale
        axs[0].hist(probabilities, bins=1000, color='steelblue', edgecolor='black')
        axs[0].set_xlabel('Infection Probability (%)')
        axs[0].set_ylabel('Count')
        axs[0].set_title(f'Probability Distribution (Linear Scale) - {category}')
        axs[0].grid(True, axis='y', linestyle='--', linewidth=0.5)

        # Plot: Log10-transformed data on linear scale
        probabilities_nozero = probabilities[probabilities > 0.001]
        log_probs = np.log10(probabilities_nozero)

        axs[1].hist(log_probs, bins=100, color='steelblue', edgecolor='black')
        axs[1].set_xlabel('log₁₀(Infection Probability [%])')
        axs[1].set_ylabel('Count')
        axs[1].set_title(f'Probability Distribution (Log10-Transformed) - {category}')
        axs[1].grid(True, axis='y', linestyle='--', linewidth=0.5)

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()
