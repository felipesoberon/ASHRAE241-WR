// qer_fitted.h — Fitted log-normal QER_i distributions per category
//
// Contains the fitted log10-normal parameters (mu10, sigma10) for each
// of the 25 ASHRAE 241 occupancy categories, derived from 1,000,000
// simulations per category of the Jones et al. (2025) model.
//
// Categories with identical QER input parameters (PBR, C_drop, d) share
// the same fitted values, but each has its own entry so that future
// updates to individual categories do not require remapping.
//
// The QER_fitted function draws a random QER_i from the fitted
// distribution using the LHS random manager, producing values
// statistically equivalent to the full 8-parameter calculation.

#pragma once
#include "random_manager.h"
#include <string>
#include <map>

// Fitted log-normal parameters in log10 space
struct QERFittedParams {
    double mu_log10;
    double sigma_log10;
};

// 25-category dictionary of fitted QER_i parameters
extern const std::map<std::string, QERFittedParams> qer_fitted_params;

// Get fitted parameters for a category (throws if not found)
const QERFittedParams& get_qer_fitted_params(const std::string& category);

// Generate a random QER_i from the fitted log-normal distribution
// Uses an LHS-buffered uniform via the random manager, transformed
// through normal_ppf to produce a log-normal draw in log10 space.
double QER_fitted(RandomNumberManager& rng, const std::string& category);