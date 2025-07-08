# model.py

from randomManager import (
    random_lognormal_lhs,
    random_uniform_lhs,
    random_beta_lhs,
    random_normal_lhs,
    random_binomial_lhs,
    random_log10normal_lhs
)

import numpy as np

# --- OCCUPANCY PARAMETERS ---

occupancy_params = {
    "Cell":           {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 14e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 30,  "volume_m3": 320, "community_rate": 0.01},
    "Dayroom":        {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 19e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 30,  "volume_m3": 600, "community_rate": 0.01},
    "Food":           {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 15e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 50,  "volume_m3": 360, "community_rate": 0.01},
    "Gym":            {"PBR_GM": 0.62, "PBR_GSD": 1.3, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 180, "volume_m3": 1900, "community_rate": 0.01},
    "Office":         {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 14e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 50,  "volume_m3": 2700, "community_rate": 0.01},
    "Retail":         {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 11e4,  "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 150, "volume_m3": 6000, "community_rate": 0.01},
    "Transportation": {"PBR_GM": 0.64, "PBR_GSD": 1.2, "Cdrop_GM": 15e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 100, "volume_m3": 10000, "community_rate": 0.01},
    "Classroom":      {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 17e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 30,  "volume_m3": 320, "community_rate": 0.01},
    "Lecture":        {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 150, "volume_m3": 1800, "community_rate": 0.01},
    "Manufacturing":  {"PBR_GM": 0.73, "PBR_GSD": 1.3, "Cdrop_GM": 15e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 70,  "volume_m3": 12000, "community_rate": 0.01},
    "Sorting":        {"PBR_GM": 0.77, "PBR_GSD": 1.3, "Cdrop_GM": 19e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 20,  "volume_m3": 4800, "community_rate": 0.01},
    "Warehouse":      {"PBR_GM": 0.77, "PBR_GSD": 1.3, "Cdrop_GM": 15e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 20,  "volume_m3": 1200, "community_rate": 0.01},
    "Exam":           {"PBR_GM": 0.62, "PBR_GSD": 1.3, "Cdrop_GM": 32e4,  "Cdrop_GSD": 1.1, "d_GM": 2.0e-6, "d_GSD": 1.1, "mask_eff": 0.3, "I0": 3,   "volume_m3": 41, "community_rate": 0.03},
    "Group":          {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 12e4,  "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.3, "I0": 20,  "volume_m3": 270, "community_rate": 0.03},
    "Patient":        {"PBR_GM": 0.62, "PBR_GSD": 1.3, "Cdrop_GM": 30e4,  "Cdrop_GSD": 1.1, "d_GM": 2.0e-6, "d_GSD": 1.1, "mask_eff": 0.3, "I0": 3,   "volume_m3": 81, "community_rate": 0.03},
    "Resident":       {"PBR_GM": 0.62, "PBR_GSD": 1.3, "Cdrop_GM": 21e4,  "Cdrop_GSD": 1.1, "d_GM": 2.0e-6, "d_GSD": 1.1, "mask_eff": 0.3, "I0": 3,   "volume_m3": 81, "community_rate": 0.03},
    "Waiting":        {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 11e4,  "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.3, "I0": 30,  "volume_m3": 270, "community_rate": 0.03},
    "Auditorium":     {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 150, "volume_m3": 4600, "community_rate": 0.01},
    "Place":          {"PBR_GM": 0.50, "PBR_GSD": 1.2, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 180, "volume_m3": 2400, "community_rate": 0.01},
    "Museum":         {"PBR_GM": 0.49, "PBR_GSD": 1.2, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 400, "volume_m3": 10000, "community_rate": 0.01},
    "Convention":     {"PBR_GM": 0.49, "PBR_GSD": 1.2, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 400, "volume_m3": 10000, "community_rate": 0.01},
    "Spectator":      {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 15e4,  "Cdrop_GSD": 1.1, "d_GM": 1.9e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 100, "volume_m3": 6300,  "community_rate": 0.01},
    "Lobbies":        {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 150, "volume_m3": 4600,  "community_rate": 0.01},
    "Common":         {"PBR_GM": 0.55, "PBR_GSD": 1.2, "Cdrop_GM": 9.8e4, "Cdrop_GSD": 1.1, "d_GM": 1.8e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 150, "volume_m3": 4600,  "community_rate": 0.01},
    "Dwelling":       {"PBR_GM": 0.51, "PBR_GSD": 4.3, "Cdrop_GM": 40e4,  "Cdrop_GSD": 1.1, "d_GM": 2.2e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 6,   "volume_m3": 540,   "community_rate": 0.01}
}

# --- FUNCTION DEFINITIONS ---

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

def sample_parameters(rng, category="Classroom"):
    occ = get_occupancy_parameters(category)
    VOL = occ["volume_m3"]                 # m³
    D = 1.0                                # h
    I0 = occ["I0"]                         # number of occupants
    mask_efficiency = occ.get("mask_eff", 0.0)
    community_rate = occ.get("community_rate", 0.01)
    lambda_bio = random_lognormal_lhs(rng, mean=np.log(0.52), sigma=np.log(1.9))  # 1/h
    gamma = random_uniform_lhs(rng, 0.42, 0.61)  # 1/h
    PBR = random_lognormal_lhs(rng, mean=np.log(occ["PBR_GM"]), sigma=np.log(occ["PBR_GSD"]))  # m³/h

    return dict(
        D=D,
        VOL=VOL,
        I0=I0,
        lambda_bio=lambda_bio,
        gamma=gamma,
        mask_efficiency=mask_efficiency,
        community_rate=community_rate,
        PBR=PBR
    )

def QER(rng, category="Classroom"):
    occ = get_occupancy_parameters(category)
    PBR = random_lognormal_lhs(rng, np.log(occ["PBR_GM"]), np.log(occ["PBR_GSD"]))
    C_drop = random_lognormal_lhs(rng, np.log(occ["Cdrop_GM"]), np.log(occ["Cdrop_GSD"]))
    d = random_lognormal_lhs(rng, np.log(occ["d_GM"]), np.log(occ["d_GSD"]))
    E = random_beta_lhs(rng, a=5.0, b=2.0, loc=2.0, scale=3.0)  # in [2, 5]
    Vdrop = (np.pi / 6.0) * (d * E) ** 3 * C_drop
    GVL_ml = random_log10normal_lhs(rng, mu=7.0, sigma=1.4)
    GVL_m3 = GVL_ml * 1e6
    VF = random_beta_lhs(rng, a=5.0, b=2.0, loc=1e-4, scale=(1e-2 - 1e-4))  # in [1e-4, 1e-2]
    RD = 1
    RTD = random_uniform_lhs(rng, 0.43, 0.65)
    DK = random_uniform_lhs(rng, 5, 15)

    VER = PBR * Vdrop * GVL_m3 * VF * RD
    QER_val = RTD * VER / DK
    return QER_val

def infection_probability(ECAi, par, rng, category="Classroom"):
    TECAi = ECAi * par['I0'] * 3.6
    phi = par['gamma'] + par['lambda_bio'] + TECAi / par['VOL']

    I0 = par['I0']
    community_rate = par['community_rate']
    n_infected = random_binomial_lhs(rng, I0, community_rate)

    infected_flag = 1 if n_infected > 0 else 0

    if n_infected > 0:
        QER_sum = sum(QER(rng, category) for _ in range(n_infected))
        mask_factor = (1 - par['mask_efficiency']) ** 2
        Q = (par['PBR'] * par['D'] * mask_factor / (phi * par['VOL'])) * QER_sum
        P = 1 - np.exp(-Q)
    else:
        P = 0

    return P, infected_flag

def compute_ECAi(par, target_P, rng, category="Classroom"):
    I0 = par['I0']
    community_rate = par['community_rate']
    n_infected = random_binomial_lhs(rng, I0, community_rate)

    infected_flag = 1 if n_infected > 0 else 0

    if n_infected > 0:
        Q = -np.log(1 - target_P)
        QER_sum = sum(QER(rng, category) for _ in range(n_infected))

        PBR = par['PBR']
        D = par['D']
        VOL = par['VOL']
        mu = par['mask_efficiency']
        mask_factor_squared = (1 - mu) ** 2
        gamma = par['gamma']
        lambda_bio = par['lambda_bio']

        term1 = (PBR * D * mask_factor_squared / Q) * QER_sum
        term2 = VOL * (gamma + lambda_bio)
        ECAi = (1 / (3.6 * I0)) * (term1 - term2)
        if ECAi < 0:
            ECAi = 0
    else:
        ECAi = 0

    return ECAi, infected_flag
