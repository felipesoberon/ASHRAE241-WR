import sys
import numpy as np
import matplotlib.pyplot as plt
from model import sample_parameters, QER, occupancy_params, infection_probability
from randomManager import RandomNumberManager

# -------- CATEGORY SELECTION --------
if len(sys.argv) > 1:
    category = sys.argv[1]
    if category not in occupancy_params:
        print(f"Error: '{category}' is not a valid category.")
        print("Available categories:", ', '.join(sorted(occupancy_params.keys())))
        sys.exit(1)
else:
    category = "Classroom"

print(f"Using occupancy category: {category}")

N = 10000
ECAi = 10  # Example value in L/s/person

rng = RandomNumberManager()

# Storage for variables from sample_parameters
D_list = []
VOL_list = []
I0_list = []
lambda_bio_list = []
gamma_list = []
mask_efficiency_list = []
community_rate_list = []
PBR_list = []

# Storage for QER, log10 QER, Probability P, and log10 P
QER_list = []
log10_QER_list = []
P_list = []
log10_P_list = []

for _ in range(N):
    par = sample_parameters(rng, category=category)
    D_list.append(par['D'])
    VOL_list.append(par['VOL'])
    I0_list.append(par['I0'])
    lambda_bio_list.append(par['lambda_bio'])
    gamma_list.append(par['gamma'])
    mask_efficiency_list.append(par['mask_efficiency'])
    community_rate_list.append(par['community_rate'])
    PBR_list.append(par['PBR'])

    qer_val = QER(rng, category=category)
    QER_list.append(qer_val)
    log10_QER_list.append(np.log10(qer_val) if qer_val > 0 else np.nan)

    # Use infection_probability from model.py
    P, _ = infection_probability(ECAi, par, rng, category=category)
    P_pct = P * 100
    P_list.append(P_pct)
    if P_pct > 0:
        log10_P_list.append(np.log10(P_pct))

# Remove any nan values from log10_QER_list
log10_QER_list = [v for v in log10_QER_list if not np.isnan(v)]

# Prepare plotting
variables = [
    ("D", D_list, "Exposure duration (h)", False),
    ("VOL", VOL_list, "Room volume (m³)", False),
    ("I₀", I0_list, "No. of occupants", False),
    ("lambda_bio", lambda_bio_list, "Biological removal rate (1/h)", False),
    ("gamma", gamma_list, "Other removal rate (1/h)", False),
    ("mask_efficiency", mask_efficiency_list, "Mask efficiency (fraction)", False),
    ("community_rate", community_rate_list, "Community infection rate (fraction)", False),
    ("PBR", PBR_list, "Pulmonary breathing rate (m³/h)", False),
    ("QER", QER_list, "Quanta emission rate (quanta/h)", True),
    ("log₁₀ QER", log10_QER_list, "log₁₀ [Quanta emission rate] (quanta/h)", False),
    ("Infection Probability", P_list, "Probability of Infection (%)", False),
    ("log₁₀ Probability", log10_P_list, "log₁₀ [Probability of Infection] (%)", False),
]

n_plots = len(variables)
n_cols = 3
n_rows = int(np.ceil(n_plots / n_cols))

fig = plt.figure(figsize=(20, 6 * n_rows))
for i, (name, data, xlabel, logx) in enumerate(variables, 1):
    plt.subplot(n_rows, n_cols, i)
    if logx:
        plt.hist(data, bins=50, alpha=0.75, edgecolor='k', log=True)
        plt.xscale('log')
    else:
        plt.hist(data, bins=50, alpha=0.75, edgecolor='k')
    plt.title(name)
    plt.xlabel(xlabel)
    plt.ylabel("Count")
plt.tight_layout(pad=6.0)
plt.subplots_adjust(hspace=0.65)
plt.show()
