#pragma once

#include <string>

// Selects the emission-rate distribution used by all simulation modes.
enum class QERDistribution {
    Jones,       // Full eight-parameter Jones et al. calculation (default)
    JonesFitted, // Category-specific fitted log10-normal Jones distributions
    Mikszewski   // SARS-CoV-2 standing/speaking distribution from Table 2
};

// Parse a command-line distribution name. Accepted names are:
// jones/full, fitted/jones-fitted, and mikszewski.
bool parse_qer_distribution(const std::string& name, QERDistribution& result);
const char* qer_distribution_name(QERDistribution distribution);
