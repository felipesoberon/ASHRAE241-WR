# model.py

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
    "Dwelling":       {"PBR_GM": 0.051,"PBR_GSD": 4.3, "Cdrop_GM": 40e4,  "Cdrop_GSD": 1.1, "d_GM": 2.2e-6, "d_GSD": 1.1, "mask_eff": 0.0, "I0": 6,   "volume_m3": 540,   "community_rate": 0.01}
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

def sample_parameters(category="Classroom"):
    occ = get_occupancy_parameters(category)
    VOL = occ["volume_m3"]                 # m³
    D = 1.0                                # h
    I0 = occ["I0"]                         # number of occupants
    mask_efficiency = occ.get("mask_eff", 0.0)     # dimensionless
    community_rate = occ.get("community_rate", 0.01)  # dimensionless
    lambda_bio = np.random.lognormal(np.log(0.52), np.log(1.9))  # 1/h
    gamma = np.random.uniform(0.42, 0.61)  # 1/h
    PBR = np.random.lognormal(np.log(occ["PBR_GM"]), np.log(occ["PBR_GSD"]))  # m³/h

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
            
def QER(category="Classroom"):
    """
    Computes the quanta emission rate per person (QER) for a given occupancy category.
    Implements equations (3), (4), and (5) from the referenced publication.
    """
    occ = get_occupancy_parameters(category)
    PBR = np.random.lognormal(np.log(occ["PBR_GM"]), np.log(occ["PBR_GSD"]))        # m³/h (pulmonary breathing rate)
    C_drop = np.random.lognormal(np.log(occ["Cdrop_GM"]), np.log(occ["Cdrop_GSD"])) # particles/m³ (exhaled droplet conc.)
    d = np.random.lognormal(np.log(occ["d_GM"]), np.log(occ["d_GSD"]))              # m (diameter)
    E = np.random.beta(5.0, 2.0) * (5.0 - 2.0) + 2.0                                # dimensionless (elongation factor)
    Vdrop = (np.pi / 6.0) * (d * E) ** 3 * C_drop                                   # m³/m³ (exhaled air) -- Eq (5)
    GVL_ml = 10 ** np.random.normal(7.0, 1.4)                                       # RNA copies/ml
    GVL_m3 = GVL_ml * 1e6                                                           # RNA copies/m³
    VF = np.random.beta(5.0, 2.0) * (1e-2 - 1e-4) + 1e-4                            # dimensionless (viable fraction)
    RD = 1                                                                          # dimensionless (respiratory droplet)
    RTD = np.random.uniform(0.43, 0.65)                                             # dimensionless (respiratory tract deposition)
    DK = np.random.uniform(5, 15)                                                   # virions/quanta

    # Eq (4): Virion emission rate (VER)
    VER = PBR * Vdrop * GVL_m3 * VF * RD                                            # virions/h

    # Eq (3): Quanta emission rate (QER)
    QER_val = RTD * VER / DK                                                        # quanta/h

    return QER_val


def infection_probability(ECAi, par, category="Classroom"):
    TECAi = ECAi * par['I0'] * 3.6  # m³/h total equivalent clean airflow ... (7)
    phi = par['gamma'] + par['lambda_bio'] + TECAi / par['VOL']  # 1/h (total removal rate) ...(6)

    # Number of infected people
    I0 = par['I0']
    community_rate = par['community_rate']
    n_infected = np.random.binomial(I0, community_rate)
    if n_infected == 0:
        n_infected = 1  # risk-averse: minimum 1 infected

    # Sum QER for all infected
    QER_sum = sum(QER(category) for _ in range(n_infected))

    mask_factor = (1 - par['mask_efficiency']) ** 2  # dimensionless

    # Quanta dose (updated, matches equation 2 logic)
    Q = (par['PBR'] * par['D'] * mask_factor / (phi * par['VOL'])) * QER_sum

    P = 1 - np.exp(-Q)  # dimensionless (probability of infection) ...(1)
    return P


def compute_ECAi(par, target_P, category="Classroom"):
    Q = -np.log(1 - target_P)  # (1)

    # Sample the number of infected people (I) from a binomial distribution
    I0 = par['I0']
    community_rate = par['community_rate']
    n_infected = np.random.binomial(I0, community_rate)
    if n_infected == 0:
        n_infected = 1  # Risk-averse: minimum 1 infected

    # Sum QER values for all infected individuals
    QER_sum = sum(QER(category) for _ in range(n_infected))
    
    # parameters from par
    PBR = par['PBR']
    D = par['D']
    VOL = par['VOL']
    mu = par['mask_efficiency']
    mask_factor_squared = (1 - mu) ** 2
    gamma = par['gamma']
    lambda_bio = par['lambda_bio']

    # Equation (2)
    term1 = (PBR * D * mask_factor_squared / Q) * QER_sum
    term2 = VOL * (gamma + lambda_bio)
    ECAi = (1 / (3.6 * I0)) * (term1 - term2)
    if ECAi < 0:
        ECAi = 0
    return ECAi

