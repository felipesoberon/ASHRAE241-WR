import sys
import numpy as np
import matplotlib.pyplot as plt

# --- Occupancy-dependent input ---

occupancy_params = {
    "Cell":           {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 14e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 30, "volume_m3": 320, "community_rate": 0.01},
    "Dayroom":        {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 19e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 30, "volume_m3": 600, "community_rate": 0.01},
    "Food":           {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 15e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 50, "volume_m3": 360, "community_rate": 0.01},
    "Gym":            {"PBR_GM": 0.62, "PBR_GSD": 1.3, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 180, "volume_m3": 1900, "community_rate": 0.01},
    "Office":         {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 14e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 50, "volume_m3": 2700, "community_rate": 0.01},
    "Retail":         {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 11e4,  "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 150, "volume_m3": 6000, "community_rate": 0.01},
    "Transportation": {"PBR_GM": 0.64, "PBR_GSD": 1.2, "Cdrop_GM": 15e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 100, "volume_m3": 10000, "community_rate": 0.01},
    "Classroom":      {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 17e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 30, "volume_m3": 320, "community_rate": 0.01},
    "Lecture":        {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 150, "volume_m3": 1800, "community_rate": 0.01},
    "Manufacturing":  {"PBR_GM": 0.73, "PBR_GSD": 1.3, "Cdrop_GM": 15e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 70, "volume_m3": 12000, "community_rate": 0.01},
    "Sorting":        {"PBR_GM": 0.77, "PBR_GSD": 1.3, "Cdrop_GM": 19e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 20, "volume_m3": 4800, "community_rate": 0.01},
    "Warehouse":      {"PBR_GM": 0.77, "PBR_GSD": 1.3, "Cdrop_GM": 15e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 20, "volume_m3": 1200, "community_rate": 0.01},
    "Exam":           {"PBR_GM": 0.62, "PBR_GSD": 1.3, "Cdrop_GM": 32e4,  "Cdrop_GSD": 1.1, "d_GM": 2.0e-6, "d_GSD": 1.1, "mask_eff": 0.3, "I0": 3, "volume_m3": 41, "community_rate": 0.03},
    "Group":          {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 12e4,  "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.3, "I0": 20, "volume_m3": 270, "community_rate": 0.03},
    "Patient":        {"PBR_GM": 0.62, "PBR_GSD": 1.3, "Cdrop_GM": 30e4,  "Cdrop_GSD": 1.1, "d_GM": 2.0e-6, "d_GSD": 1.1, "mask_eff": 0.3, "I0": 3, "volume_m3": 81, "community_rate": 0.03},
    "Resident":       {"PBR_GM": 0.62, "PBR_GSD": 1.3, "Cdrop_GM": 21e4,  "Cdrop_GSD": 1.1, "d_GM": 2.0e-6, "d_GSD": 1.1, "mask_eff": 0.3, "I0": 3, "volume_m3": 81, "community_rate": 0.03},
    "Waiting":        {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 11e4,  "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.3, "I0": 30, "volume_m3": 270, "community_rate": 0.03},
    "Auditorium":     {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 150, "volume_m3": 4600, "community_rate": 0.01},
    "Place":          {"PBR_GM": 0.50, "PBR_GSD": 1.2, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 180, "volume_m3": 2400, "community_rate": 0.01},
    "Museum":         {"PBR_GM": 0.49, "PBR_GSD": 1.2, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 400, "volume_m3": 10000, "community_rate": 0.01},
    "Convention":     {"PBR_GM": 0.49, "PBR_GSD": 1.2, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 400, "volume_m3": 10000, "community_rate": 0.01},
    "Spectator":      {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 15e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 100, "volume_m3": 6300, "community_rate": 0.01},
    "Lobbies":        {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 150, "volume_m3": 4600, "community_rate": 0.01},
    "Common":         {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 150, "volume_m3": 4600, "community_rate": 0.01},
    "Dwelling":       {"PBR_GM": 0.051,"PBR_GSD": 4.3, "Cdrop_GM": 40e4, "Cdrop_GSD": 1.1, "d_GM": 2.2e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 6, "volume_m3": 540, "community_rate": 0.01}
}

def get_occupancy_parameters(category):
    default = {
        "PBR_GM": 0.55,
        "PBR_GSD": 1.2,
        "Cdrop_GM": 17e4,
        "Cdrop_GSD": 1.1,
        "d_GM": 1.9e-6,
        "d_GSD": 1.1,
        "mask_eff": 0.0,
        "I0": 30,
        "volume_m3": 320,
        "community_rate": 0.01
    }
    return occupancy_params.get(category, default)

def sample_parameters(category="Classroom"):
    occ = get_occupancy_parameters(category)
    PBR = np.random.lognormal(np.log(occ["PBR_GM"]), np.log(occ["PBR_GSD"]))  # m^3/h
    C_drop = np.random.lognormal(np.log(occ["Cdrop_GM"]), np.log(occ["Cdrop_GSD"]))  # particles/m^3
    d = np.random.lognormal(np.log(occ["d_GM"]), np.log(occ["d_GSD"]))  # m
    VOL = occ["volume_m3"]  # m³
    D = 1.0  # h
    I0 = occ["I0"]  # number of occupants
    E = np.random.beta(5.0, 2.0) * (5.0 - 2.0) + 2.0  # dimensionless
    hydrated_droplet_volume = (np.pi / 6.0) * (d * E)**3 * C_drop  # m^3/m^3 (exhaled air)
    GVL_ml = 10**np.random.normal(7.0, 1.4)  # RNA copies/ml
    GVL_m3 = GVL_ml * 1e6  # RNA copies/m^3
    VF = np.random.beta(5.0, 2.0) * (1e-2 - 1e-4) + 1e-4  # dimensionless
    RTD = np.random.uniform(0.43, 0.65)  # dimensionless
    RD = 1  # dimensionless
    DK = np.random.uniform(5, 15)  # virions/quanta
    lambda_bio = np.random.lognormal(np.log(0.52), np.log(1.9))  # 1/h
    gamma = np.random.uniform(0.42, 0.61)  # 1/h
    mask_efficiency = occ.get("mask_eff", 0.0)  # dimensionless (default to 0.0 if missing)
    community_rate = occ.get("community_rate", 0.01) # community infection rate 

    return dict(PBR=PBR, D=D, VOL=VOL, I0=I0,
            hydrated_droplet_volume=hydrated_droplet_volume,
            GVL=GVL_m3, VF=VF, RD=RD, RTD=RTD, DK=DK,
            lambda_bio=lambda_bio, gamma=gamma,
            mask_efficiency=mask_efficiency,
            community_rate=community_rate)
            

def infection_probability(ECAi, par):
    F_ec = ECAi * par['I0'] * 3.6  # m³/h total equivalent clean airflow
    phi = par['gamma'] + par['lambda_bio'] + F_ec / par['VOL']  # 1/h (total removal rate)
    liquid_per_hour = par['PBR'] * par['hydrated_droplet_volume']  # m³/h of hydrated droplets
    VER = liquid_per_hour * par['GVL'] * par['VF'] * par['RD']  # virions/h
    QER = par['RTD'] * VER / par['DK']  # quanta/h
    mask_factor = (1 - par['mask_efficiency']) ** 2  # dimensionless
    Q = (par['PBR'] * par['D'] * mask_factor / (phi * par['VOL'])) * (   par['I0'] * par['community_rate']    * QER) # dimensionless (quanta dose)
   #Q = (par['PBR'] * par['D'] * mask_factor / (phi * par['VOL'])) * (   max(1, round(par['I0'] * par['community_rate']))    * QER) # dimensionless (quanta dose)
   #Q = (par['PBR'] * par['D'] * mask_factor / (phi * par['VOL'])) * (   int(np.random.rand() < par['community_rate'])       * QER) # dimensionless (quanta dose)
     

    P = 1 - np.exp(-Q)  # dimensionless (probability of infection)
    return P

if __name__ == "__main__":
    N = 10000
    ECAi_values = list(range(5, 100, 5))
    category = sys.argv[1] if len(sys.argv) > 1 else "Classroom"

    print(f"\n{'ECAi':>6}  {'Mean':>9}  {'Median':>9}  {'Min':>9}  {'Max':>9}  {'96th pct':>9}")
    print('-' * 60)

    last_ECAi = None
    last_probs = None

    for ECAi in ECAi_values:
        probs = []
        for _ in range(N):
            par = sample_parameters(category=category)
            P = infection_probability(ECAi, par)
            probs.append(P)
        probs = np.array(probs)
        meanP = np.mean(probs)
        medianP = np.median(probs)
        minP = np.min(probs)
        maxP = np.max(probs)
        pct96 = np.percentile(probs, 96)
        print(f"{ECAi:6}  {meanP*100:9.4f}  {medianP*100:9.4f}  {minP*100:9.4f}  {maxP*100:9.4f}  {pct96*100:9.4f}")

        last_ECAi = ECAi
        last_probs = probs

        if pct96 < 0.001:  # 0.1%
            break

    # Plot only the last run
    plt.figure(figsize=(10, 6))
    bins = np.linspace(0, 1, 10000)  # Infection probability in fraction (0–1)
    hist, _ = np.histogram(last_probs, bins=bins, density=False)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    plt.plot(bin_centers * 100, hist, label=f"ECAi = {last_ECAi}")  # % scale

    plt.axvline(x=0.1, color='red', linestyle='--', linewidth=1.2)  # Vertical at 0.1%
    plt.yscale('log')
    plt.xlim(0, 1)
    plt.xlabel("Infection Probability (%)")
    plt.ylabel("Density (log scale)")
    plt.title("Infection Probability Distribution")
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.6)
    plt.tight_layout()
    plt.show()
