// qer_fitted.cpp — Fitted log-normal QER_i distributions per category
//
// Parameters derived from fit_qer_lognormal.py analysis of
// probECAi_1M_inputs.bin (1M simulations per category).
// All values are in log10 space, consistent with Mikszewski et al. (2022)
// and Jones et al. (2025).

#include "qer_fitted.h"
#include "random_manager.h"
#include <cmath>
#include <stdexcept>

// 25-category dictionary. Categories sharing the same QER input
// parameters (PBR, C_drop, d) have identical fitted values.
// Values from the combined group fits (output/qer_lognormal_fit_results.csv).
const std::map<std::string, QERFittedParams> qer_fitted_params = {
    {"Cell",           {-1.123527, 1.428904}},  // Group: Cell, Office
    {"Dayroom",        {-0.992113, 1.429109}},
    {"Food",           {-1.093306, 1.427494}},  // Group: Food, Spectator
    {"Gym",            {-1.298071, 1.432531}},
    {"Office",         {-1.123527, 1.428904}},  // Group: Cell, Office
    {"Retail",         {-1.299721, 1.429381}},  // Group: Retail, Waiting
    {"Transportation", {-1.026098, 1.429956}},
    {"Classroom",      {-1.039436, 1.427156}},
    {"Lecture",        {-1.349078, 1.428965}},  // Group: Lecture, Auditorium, Lobbies, Common
    {"Manufacturing",  {-0.970826, 1.430100}},
    {"Sorting",        {-0.844935, 1.431221}},
    {"Warehouse",      {-0.947311, 1.432394}},
    {"Exam",           {-0.645766, 1.430231}},
    {"Group",          {-1.260304, 1.428842}},
    {"Patient",        {-0.673900, 1.430657}},
    {"Resident",       {-0.827809, 1.430225}},
    {"Waiting",        {-1.299721, 1.429381}},  // Group: Retail, Waiting
    {"Auditorium",     {-1.349078, 1.428965}},  // Group: Lecture, Auditorium, Lobbies, Common
    {"Place",          {-1.389396, 1.427938}},
    {"Museum",         {-1.399522, 1.429133}},  // Group: Museum, Convention
    {"Convention",     {-1.399522, 1.429133}},  // Group: Museum, Convention
    {"Spectator",      {-1.093306, 1.427494}},  // Group: Food, Spectator
    {"Lobbies",        {-1.349078, 1.428965}},  // Group: Lecture, Auditorium, Lobbies, Common
    {"Common",         {-1.349078, 1.428965}},  // Group: Lecture, Auditorium, Lobbies, Common
    {"Dwelling",       {-0.509175, 1.561887}},
};

static const QERFittedParams default_fitted = {-1.0, 1.43};

const QERFittedParams& get_qer_fitted_params(const std::string& category) {
    auto it = qer_fitted_params.find(category);
    if (it == qer_fitted_params.end())
        return default_fitted;
    return it->second;
}

// Generate a random QER_i from the fitted log-normal distribution.
// Uses a dedicated LHS buffer ("qer_fitted") for stratified sampling,
// then transforms through normal_ppf to produce a log10-normal draw.
//
// QER_i = 10^(normal_ppf(u) * sigma10 + mu10)
//
// This is equivalent to scipy.stats.lognorm.ppf(u, s=sigma_ln,
// scale=exp(mu_ln)) where sigma_ln = sigma10 * ln(10) and
// mu_ln = mu10 * ln(10).
double QER_fitted(RandomNumberManager& rng, const std::string& category) {
    const auto& p = get_qer_fitted_params(category);
    double u = rng.get("qer_fitted");
    // Clamp to avoid infinities from normal_ppf at exact 0 or 1
    if (u <= 0.0) u = 1e-15;
    if (u >= 1.0) u = 1.0 - 1e-15;
    double z = normal_ppf(u);               // standard normal quantile
    double log10_qer = z * p.sigma_log10 + p.mu_log10;
    return std::pow(10.0, log10_qer);       // QER_i in quanta/h
}