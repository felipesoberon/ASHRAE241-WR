#pragma once
#include "random_manager.h"
#include "qer_distribution.h"
#include <string>
#include <map>
#include <vector>
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

// Category names in Python insertion order (not alphabetical).
extern const std::vector<std::string> category_order;

const OccupancyParams& get_occupancy_parameters(const std::string& category);

// Healthcare spaces are those with the elevated default community infection
// rate (3%): Exam, Group, Patient, Resident, Waiting (per Jones et al. 2025).
// General spaces default to 1%. Used to apply group-wide rate overrides.
bool is_healthcare_category(const std::string& category);

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

// QER calculation. The default is the full Jones et al. method.
double QER(RandomNumberManager& rng, const std::string& category,
           QERDistribution distribution = {});
// Backward-compatible overload: true selects JonesFitted, false selects Jones.
double QER(RandomNumberManager& rng, const std::string& category,
           bool use_fitted);

// QER with inputs: returns {QER_value, sampled parameter values}
struct QERInputs {
    double PBR_qer;
    double C_drop;
    double d;
    double E;
    double Vdrop;
    double GVL_ml;
    double GVL_m3;
    double VF;
    double RTD;
    double DK;
    double VER;
    double QER_val;
};
std::pair<double, QERInputs> QER_with_inputs(
    RandomNumberManager& rng, const std::string& category,
    QERDistribution distribution = {});

// Returns {probability, infected_flag}
// override_community_rate: a negative value (default) means "unset" and
// falls back to the category's default community rate; any value in
// [0, 1], including an explicit 0, is honored as-is.
// Throws std::invalid_argument if require_infectors is true and the
// effective community rate is 0: no infector could ever be drawn.
std::pair<double, int> infection_probability(
    double ECAi, const SimParameters& par, RandomNumberManager& rng,
    const std::string& category, bool require_infectors = false,
    double override_community_rate = -1,
    QERDistribution distribution = {});

// Infection probability with inputs: returns {probability, flag,
// intermediate values for driver analysis}
struct SimInputs {
    int    n_infected;
    double phi;
    double QER_sum;
    double mask_factor;
    double Q;
    double P;
    // QER breakdown for the first infector only (to keep file size
    // manageable; QER_sum captures the total across all infectors)
    QERInputs qer_first;
};
struct InfectionResultWithInputs {
    double P;
    int    infected_flag;
    SimInputs inputs;
};
InfectionResultWithInputs infection_probability_with_inputs(
    double ECAi, const SimParameters& par, RandomNumberManager& rng,
    const std::string& category, bool require_infectors = false,
    double override_community_rate = -1,
    QERDistribution distribution = {});

// Returns {ECAi_value, infected_flag}
// See infection_probability() for override_community_rate semantics.
std::pair<double, int> compute_ECAi(
    const SimParameters& par, double target_P, RandomNumberManager& rng,
    const std::string& category, bool require_infectors = false,
    double override_community_rate = -1,
    QERDistribution distribution = {});