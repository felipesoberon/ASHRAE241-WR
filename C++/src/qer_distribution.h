#pragma once

#include <string>
#include <vector>

// Selects the QER_i source used by all simulation modes.
enum class QERDistributionKind {
    Jones,
    JonesFitted,
    Mikszewski
};

struct QERDistribution {
    QERDistributionKind kind = QERDistributionKind::Jones;
    std::string profile;  // Used only for Mikszewski profiles.
};

// Parse: jones, jones-fitted, or one of the 30
// mikszewski-<organism>-<activity> profile names.
bool parse_qer_distribution(const std::string& name, QERDistribution& result);
std::string qer_distribution_name(const QERDistribution& distribution);
std::vector<std::string> mikszewski_distribution_names();
