#pragma once
#include "random_manager.h"
#include <string>
#include <map>
#include <utility>

// Occupancy parameters for each ASHRAE 241 category
struct OccupancyParams {
    double PBR_GM;
    double PBR_GSD;
    double Cdrop_GM;
    double Cdrop_GSD;
    double d_GM;
    double d_GSD;
    double mask_eff;
    int    I0;
    double volume_m3;
    double community_rate;
    double ECAi;
};

extern const std::map<std::string, OccupancyParams> occupancy_params;

const OccupancyParams& get_occupancy_parameters(const std::string& category);

// Sampled simulation parameters
struct SimParameters {
    double D;
    double VOL;
    int    I0;
    double lambda_bio;
    double gamma;
    double mask_efficiency;
    double community_rate;
    double PBR;
};

SimParameters sample_parameters(RandomNumberManager& rng, const std::string& category);

double QER(RandomNumberManager& rng, const std::string& category);

// Returns {probability, infected_flag}
std::pair<double, int> infection_probability(
    double ECAi, const SimParameters& par, RandomNumberManager& rng,
    const std::string& category, bool require_infectors = false,
    double override_community_rate = 0);

// Returns {ECAi_value, infected_flag}
std::pair<double, int> compute_ECAi(
    const SimParameters& par, double target_P, RandomNumberManager& rng,
    const std::string& category, bool require_infectors = false,
    double override_community_rate = 0);