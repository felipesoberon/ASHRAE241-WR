# ecai_simulation.py

import numpy as np
from model import occupancy_params, sample_parameters, compute_ECAi

if __name__ == "__main__":    

    target_P = 0.001
    N = 10000    
    print("\nCalculating ECAi from fixed infection probability using ", N, " samples for all categories...")

    print("\n| {:^15} | {:^22} | {:^33} | {:^36} |".format("Category", "Target Probability (%)", "96th Percentile ECAi (L/s/person)", "Rounded Minimum ECAi (L/s/person)"))
    print("|{:-^17}|{:-^24}|{:-^35}|{:-^38}|".format('', '', '', ''))

    for category in occupancy_params:
        ECAi_list = []
        for _ in range(N):
            par = sample_parameters(category=category)
            ECAi_val = compute_ECAi(par, target_P, category=category)  # pass category here!
            ECAi_list.append(ECAi_val)

        ECAi_array = np.array(ECAi_list)
        percentile_96 = np.percentile(ECAi_array, 96)
        rounded_min_ECAi = int(np.ceil(percentile_96 / 5.0)) * 5

        print("| {:<15} | {:<22} | {:>33.2f} | {:>36} |".format(category, target_P*100, percentile_96, rounded_min_ECAi))
