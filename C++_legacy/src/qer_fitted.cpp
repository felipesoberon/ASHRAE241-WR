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

static double qer_from_params(RandomNumberManager& rng,
                              const QERFittedParams& p) {
    double u = rng.get("qer_fitted");
    // Clamp to avoid infinities from normal_ppf at exact 0 or 1.
    if (u <= 0.0) u = 1e-15;
    if (u >= 1.0) u = 1.0 - 1e-15;
    double z = normal_ppf(u);
    return std::pow(10.0, z * p.sigma_log10 + p.mu_log10);
}

double QER_fitted(RandomNumberManager& rng, const std::string& category) {
    return qer_from_params(rng, get_qer_fitted_params(category));
}

// Mikszewski et al. (2022), Table 2 and Table 1. Each profile uses
// mu_log10 = log10(Table 2 median) and the organism-specific sigma_log10
// from Table 1.
const std::map<std::string, QERFittedParams> mikszewski_params = {
    {"mikszewski-sars-cov-1-resting",              {std::log10(0.0084), 1.30}},
    {"mikszewski-sars-cov-1-standing-speaking",   {std::log10(0.042),  1.30}},
    {"mikszewski-sars-cov-1-light-speaking-loudly",{std::log10(0.71),  1.30}},
    {"mikszewski-mers-resting",                   {std::log10(0.011),  1.60}},
    {"mikszewski-mers-standing-speaking",        {std::log10(0.056),  1.60}},
    {"mikszewski-mers-light-speaking-loudly",     {std::log10(0.96),  1.60}},
    {"mikszewski-tb-on-treatment-resting",        {std::log10(0.020),  1.40}},
    {"mikszewski-tb-on-treatment-standing-speaking", {std::log10(0.098), 1.40}},
    {"mikszewski-tb-on-treatment-light-speaking-loudly", {std::log10(1.7), 1.40}},
    {"mikszewski-influenza-resting",              {std::log10(0.035), 0.84}},
    {"mikszewski-influenza-standing-speaking",    {std::log10(0.17),  0.84}},
    {"mikszewski-influenza-light-speaking-loudly",{std::log10(3.0),   0.84}},
    {"mikszewski-coxsackievirus-resting",         {std::log10(0.062), 1.10}},
    {"mikszewski-coxsackievirus-standing-speaking", {std::log10(0.31), 1.10}},
    {"mikszewski-coxsackievirus-light-speaking-loudly", {std::log10(5.2), 1.10}},
    {"mikszewski-rhinovirus-resting",             {std::log10(0.21),  0.83}},
    {"mikszewski-rhinovirus-standing-speaking",   {std::log10(1.0),   0.83}},
    {"mikszewski-rhinovirus-light-speaking-loudly",{std::log10(18.0), 0.83}},
    {"mikszewski-sars-cov-2-resting",             {std::log10(0.55),  1.20}},
    {"mikszewski-sars-cov-2-standing-speaking",   {std::log10(2.7),   1.20}},
    {"mikszewski-sars-cov-2-light-speaking-loudly",{std::log10(46.0), 1.20}},
    {"mikszewski-tb-untreated-resting",           {std::log10(0.62),  1.30}},
    {"mikszewski-tb-untreated-standing-speaking",{std::log10(3.1),   1.30}},
    {"mikszewski-tb-untreated-light-speaking-loudly", {std::log10(52.0), 1.30}},
    {"mikszewski-adenovirus-resting",             {std::log10(0.78),  0.95}},
    {"mikszewski-adenovirus-standing-speaking",   {std::log10(3.9),   0.95}},
    {"mikszewski-adenovirus-light-speaking-loudly",{std::log10(66.0), 0.95}},
    {"mikszewski-measles-resting",                {std::log10(3.1),   1.60}},
    {"mikszewski-measles-standing-speaking",      {std::log10(15.0),  1.60}},
    {"mikszewski-measles-light-speaking-loudly",  {std::log10(260.0), 1.60}},
};

const QERFittedParams mikszewski_qer_params = {
    std::log10(2.7), 1.2
};

double QER_mikszewski(RandomNumberManager& rng, const std::string& profile) {
    auto it = mikszewski_params.find(profile);
    if (it == mikszewski_params.end())
        throw std::invalid_argument("Unknown Mikszewski QER profile: " + profile);
    return qer_from_params(rng, it->second);
}

double QER_mikszewski(RandomNumberManager& rng) {
    return qer_from_params(rng, mikszewski_qer_params);
}