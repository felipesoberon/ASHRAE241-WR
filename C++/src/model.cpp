// model.cpp — ASHRAE 241 risk model: occupancy params, QER, infection probability, ECAi

#include "model.h"

#include <cmath>
#include <stdexcept>

namespace {
constexpr double PI = 3.14159265358979323846;
}

// --- OCCUPANCY PARAMETERS ---
// Ported exactly from model.py lines 17-41. Do not change any value.

const std::map<std::string, OccupancyParams> occupancy_params = {
    {"Cell",          {0.55, 1.2, 14e4, 1.1, 1.9e-6, 1.1, 0.0,  30,   320,  0.01, 15}},
    {"Dayroom",       {0.55, 1.2, 19e4, 1.1, 1.9e-6, 1.1, 0.0,  30,   600,  0.01, 20}},
    {"Food",          {0.55, 1.2, 15e4, 1.1, 1.9e-6, 1.1, 0.0,  50,   360,  0.01, 30}},
    {"Gym",           {0.62, 1.3, 9.8e4,1.1, 1.8e-6, 1.1, 0.0, 180,  1900,  0.01, 40}},
    {"Office",        {0.55, 1.2, 14e4, 1.1, 1.9e-6, 1.1, 0.0,  50,  2700,  0.01, 15}},
    {"Retail",        {0.55, 1.2, 11e4, 1.1, 1.8e-6, 1.1, 0.0, 150,  6000,  0.01, 20}},
    {"Transportation",{0.64, 1.2, 15e4, 1.1, 1.9e-6, 1.1, 0.0, 100, 10000,  0.01, 30}},
    {"Classroom",     {0.55, 1.2, 17e4, 1.1, 1.9e-6, 1.1, 0.0,  30,   320,  0.01, 20}},
    {"Lecture",       {0.55, 1.2, 9.8e4,1.1, 1.8e-6, 1.1, 0.0, 150,  1800,  0.01, 25}},
    {"Manufacturing", {0.73, 1.3, 15e4, 1.1, 1.9e-6, 1.1, 0.0,  70, 12000,  0.01, 25}},
    {"Sorting",       {0.77, 1.3, 19e4, 1.1, 1.9e-6, 1.1, 0.0,  20,  4800,  0.01, 10}},
    {"Warehouse",     {0.77, 1.3, 15e4, 1.1, 1.9e-6, 1.1, 0.0,  20,  1200,  0.01, 10}},
    {"Exam",          {0.62, 1.3, 32e4, 1.1, 2.0e-6, 1.1, 0.3,   3,    41,  0.03, 20}},
    {"Group",         {0.55, 1.2, 12e4, 1.1, 1.8e-6, 1.1, 0.3,  20,   270,  0.03, 35}},
    {"Patient",       {0.62, 1.3, 30e4, 1.1, 2.0e-6, 1.1, 0.3,   3,    81,  0.03, 35}},
    {"Resident",      {0.62, 1.3, 21e4, 1.1, 2.0e-6, 1.1, 0.3,   3,    81,  0.03, 25}},
    {"Waiting",       {0.55, 1.2, 11e4, 1.1, 1.8e-6, 1.1, 0.3,  30,   270,  0.03, 45}},
    {"Auditorium",    {0.55, 1.2, 9.8e4,1.1, 1.8e-6, 1.1, 0.0, 150,  4600,  0.01, 25}},
    {"Place",         {0.50, 1.2, 9.8e4,1.1, 1.8e-6, 1.1, 0.0, 180,  2400,  0.01, 25}},
    {"Museum",        {0.49, 1.2, 9.8e4,1.1, 1.8e-6, 1.1, 0.0, 400, 10000,  0.01, 30}},
    {"Convention",    {0.49, 1.2, 9.8e4,1.1, 1.8e-6, 1.1, 0.0, 400, 10000,  0.01, 30}},
    {"Spectator",     {0.55, 1.2, 15e4, 1.1, 1.9e-6, 1.1, 0.0, 100,  6300,  0.01, 25}},
    {"Lobbies",       {0.55, 1.2, 9.8e4,1.1, 1.8e-6, 1.1, 0.0, 150,  4600,  0.01, 25}},
    {"Common",        {0.55, 1.2, 9.8e4,1.1, 1.8e-6, 1.1, 0.0, 150,  4600,  0.01, 25}},
    {"Dwelling",      {0.51, 4.3, 40e4, 1.1, 2.2e-6, 1.1, 0.0,   6,   540,  0.01, 15}},
};

static const OccupancyParams default_params = {
    0.55, 1.2, 17e4, 1.1, 1.9e-6, 1.1, 0.0, 30, 320, 0.01, 0
};

const OccupancyParams& get_occupancy_parameters(const std::string& category) {
    auto it = occupancy_params.find(category);
    if (it == occupancy_params.end())
        return default_params;
    return it->second;
}

// --- sample_parameters ---
SimParameters sample_parameters(RandomNumberManager& rng, const std::string& category) {
    const auto& occ = get_occupancy_parameters(category);
    SimParameters par;
    par.VOL = occ.volume_m3;
    par.D = 1.0;
    par.I0 = occ.I0;
    par.mask_efficiency = occ.mask_eff;
    par.community_rate = occ.community_rate;
    par.lambda_bio = random_lognormal_lhs(rng, std::log(0.52), std::log(1.9));
    par.gamma = random_uniform_lhs(rng, 0.42, 0.61);
    par.PBR = random_lognormal_lhs(rng, std::log(occ.PBR_GM), std::log(occ.PBR_GSD));
    return par;
}

// --- QER ---
double QER(RandomNumberManager& rng, const std::string& category) {
    const auto& occ = get_occupancy_parameters(category);
    double PBR    = random_lognormal_lhs(rng, std::log(occ.PBR_GM),  std::log(occ.PBR_GSD));
    double C_drop = random_lognormal_lhs(rng, std::log(occ.Cdrop_GM), std::log(occ.Cdrop_GSD));
    double d      = random_lognormal_lhs(rng, std::log(occ.d_GM),     std::log(occ.d_GSD));
    double E      = random_beta_lhs(rng, 5.0, 2.0, 2.0, 3.0);  // in [2, 5]
    double Vdrop  = (PI / 6.0) * std::pow(d * E, 3) * C_drop;
    double GVL_ml = random_log10normal_lhs(rng, 7.0, 1.4);
    double GVL_m3 = GVL_ml * 1e6;
    double VF     = random_beta_lhs(rng, 5.0, 2.0, 1e-4, (1e-2 - 1e-4));  // in [1e-4, 1e-2]
    double RD = 1;
    double RTD = random_uniform_lhs(rng, 0.43, 0.65);
    double DK  = random_uniform_lhs(rng, 5.0, 15.0);
    double VER = PBR * Vdrop * GVL_m3 * VF * RD;
    return RTD * VER / DK;
}

// --- infection_probability ---
std::pair<double, int> infection_probability(
    double ECAi, const SimParameters& par, RandomNumberManager& rng,
    const std::string& category, bool require_infectors,
    double override_community_rate)
{
    double TECAi = ECAi * par.I0 * 3.6;
    double phi = par.gamma + par.lambda_bio + TECAi / par.VOL;
    double comm_rate = (override_community_rate > 0) ? override_community_rate : par.community_rate;

    int n_infected = random_binomial_lhs(rng, par.I0, comm_rate);
    if (require_infectors) {
        while (n_infected == 0)
            n_infected = random_binomial_lhs(rng, par.I0, comm_rate);
    }

    int infected_flag = (n_infected > 0) ? 1 : 0;
    double P = 0.0;

    if (n_infected > 0) {
        double QER_sum = 0;
        for (int i = 0; i < n_infected; i++)
            QER_sum += QER(rng, category);
        double mask_factor = std::pow(1.0 - par.mask_efficiency, 2);
        double Q = (par.PBR * par.D * mask_factor / (phi * par.VOL)) * QER_sum;
        P = 1.0 - std::exp(-Q);
    }

    return {P, infected_flag};
}

// --- compute_ECAi ---
std::pair<double, int> compute_ECAi(
    const SimParameters& par, double target_P, RandomNumberManager& rng,
    const std::string& category, bool require_infectors,
    double override_community_rate)
{
    double comm_rate = (override_community_rate > 0) ? override_community_rate : par.community_rate;

    int n_infected = random_binomial_lhs(rng, par.I0, comm_rate);
    if (require_infectors) {
        while (n_infected == 0)
            n_infected = random_binomial_lhs(rng, par.I0, comm_rate);
    }

    int infected_flag = (n_infected > 0) ? 1 : 0;
    double ECAi_val = 0;

    if (n_infected > 0) {
        double Q = -std::log(1.0 - target_P);
        double QER_sum = 0;
        for (int i = 0; i < n_infected; i++)
            QER_sum += QER(rng, category);

        double mask_factor_sq = std::pow(1.0 - par.mask_efficiency, 2);
        double term1 = (par.PBR * par.D * mask_factor_sq / Q) * QER_sum;
        double term2 = par.VOL * (par.gamma + par.lambda_bio);
        ECAi_val = (1.0 / (3.6 * par.I0)) * (term1 - term2);
        if (ECAi_val < 0) ECAi_val = 0;
    }

    return {ECAi_val, infected_flag};
}