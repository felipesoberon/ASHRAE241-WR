import sys
import numpy as np
import matplotlib.pyplot as plt
from model import get_occupancy_parameters, occupancy_params, QER
from randomManager import (
    RandomNumberManager,
    random_lognormal_lhs,
    random_uniform_lhs,
    random_beta_lhs,
    random_log10normal_lhs
)

# --- CATEGORY SELECTION ---
if len(sys.argv) > 1:
    category = sys.argv[1]
    if category not in occupancy_params:
        print(f"Error: '{category}' is not a valid category.")
        print("Available categories:", ', '.join(sorted(occupancy_params.keys())))
        sys.exit(1)
else:
    category = "Classroom"

print(f"Sampling QER internals for category: {category}")

N = 10000
occ = get_occupancy_parameters(category)
rng = RandomNumberManager()

# Sample arrays
PBR_qer_list = []
C_drop_list = []
d_list = []
E_list = []
Vdrop_list = []
GVL_ml_list = []
GVL_m3_list = []
log10_GVL_ml_list = []
log10_GVL_m3_list = []
VF_list = []
RTD_list = []
DK_list = []
QER_list = []
log10_QER_list = []

for _ in range(N):
    PBR_qer = random_lognormal_lhs(rng, mean=np.log(occ["PBR_GM"]), sigma=np.log(occ["PBR_GSD"]))  # m³/h
    C_drop = random_lognormal_lhs(rng, mean=np.log(occ["Cdrop_GM"]), sigma=np.log(occ["Cdrop_GSD"]))  # particles/m³
    d = random_lognormal_lhs(rng, mean=np.log(occ["d_GM"]), sigma=np.log(occ["d_GSD"]))  # m
    E = random_beta_lhs(rng, a=5.0, b=2.0, loc=2.0, scale=3.0)  # in [2, 5]
    Vdrop = (np.pi / 6.0) * (d * E) ** 3 * C_drop
    GVL_ml = random_log10normal_lhs(rng, mu=7.0, sigma=1.4)
    GVL_m3 = GVL_ml * 1e6
    log10_GVL_ml = np.log10(GVL_ml)
    log10_GVL_m3 = np.log10(GVL_m3)
    VF = random_beta_lhs(rng, a=5.0, b=2.0, loc=1e-4, scale=(1e-2 - 1e-4))
    RTD = random_uniform_lhs(rng, 0.43, 0.65)
    DK = random_uniform_lhs(rng, 5, 15)
    RD = 1  # always 1

    VER = PBR_qer * Vdrop * GVL_m3 * VF * RD
    QER_val = RTD * VER / DK
    log10_QER = np.log10(QER_val) if QER_val > 0 else np.nan

    # Save
    PBR_qer_list.append(PBR_qer)
    C_drop_list.append(C_drop)
    d_list.append(d)
    E_list.append(E)
    Vdrop_list.append(Vdrop)
    GVL_ml_list.append(GVL_ml)
    GVL_m3_list.append(GVL_m3)
    log10_GVL_ml_list.append(log10_GVL_ml)
    log10_GVL_m3_list.append(log10_GVL_m3)
    VF_list.append(VF)
    RTD_list.append(RTD)
    DK_list.append(DK)
    QER_list.append(QER_val)
    log10_QER_list.append(log10_QER)

# Remove any nan values from log10_QER_list (in case of zero QER)
log10_QER_list = [v for v in log10_QER_list if not np.isnan(v)]

# Variables to plot: (name, list, xlabel, logx)
variables = [
    ("PBR_qer", PBR_qer_list, "Pulmonary breathing rate (m³/h)", False),
    ("C_drop", C_drop_list, "Exhaled droplet conc. (particles/m³)", True),
    ("d", d_list, "Droplet diameter (m)", False),
    ("E", E_list, "Elongation factor (dimensionless)", False),
    ("Vdrop", Vdrop_list, "Exhaled air fraction (m³/m³)", True),
    ("log₁₀ GVL_ml", log10_GVL_ml_list, "log₁₀ [Viral load] (RNA copies/ml)", False),
    ("log₁₀ GVL_m3", log10_GVL_m3_list, "log₁₀ [Viral load] (RNA copies/m³)", False),
    ("VF", VF_list, "Viable fraction (dimensionless)", False),
    ("RTD", RTD_list, "Respiratory tract deposition (dimensionless)", False),
    ("DK", DK_list, "Virions/quanta (dimensionless)", False),
    ("QER", QER_list, "Quanta emission rate (quanta/h)", True),
    ("log₁₀ QER", log10_QER_list, "log₁₀ [Quanta emission rate] (quanta/h)", False),
]

n_plots = len(variables)
n_cols = 3
n_rows = int(np.ceil(n_plots / n_cols))

fig = plt.figure(figsize=(18, 5 * n_rows))
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
plt.tight_layout(pad=4.0)
plt.subplots_adjust(hspace=0.55)
plt.show()
