#include "qer_distribution.h"

#include <algorithm>
#include <cctype>

namespace {
std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

const std::vector<std::string> profiles = {
    "mikszewski-sars-cov-1-resting",
    "mikszewski-sars-cov-1-standing-speaking",
    "mikszewski-sars-cov-1-light-speaking-loudly",
    "mikszewski-mers-resting",
    "mikszewski-mers-standing-speaking",
    "mikszewski-mers-light-speaking-loudly",
    "mikszewski-tb-on-treatment-resting",
    "mikszewski-tb-on-treatment-standing-speaking",
    "mikszewski-tb-on-treatment-light-speaking-loudly",
    "mikszewski-influenza-resting",
    "mikszewski-influenza-standing-speaking",
    "mikszewski-influenza-light-speaking-loudly",
    "mikszewski-coxsackievirus-resting",
    "mikszewski-coxsackievirus-standing-speaking",
    "mikszewski-coxsackievirus-light-speaking-loudly",
    "mikszewski-rhinovirus-resting",
    "mikszewski-rhinovirus-standing-speaking",
    "mikszewski-rhinovirus-light-speaking-loudly",
    "mikszewski-sars-cov-2-resting",
    "mikszewski-sars-cov-2-standing-speaking",
    "mikszewski-sars-cov-2-light-speaking-loudly",
    "mikszewski-tb-untreated-resting",
    "mikszewski-tb-untreated-standing-speaking",
    "mikszewski-tb-untreated-light-speaking-loudly",
    "mikszewski-adenovirus-resting",
    "mikszewski-adenovirus-standing-speaking",
    "mikszewski-adenovirus-light-speaking-loudly",
    "mikszewski-measles-resting",
    "mikszewski-measles-standing-speaking",
    "mikszewski-measles-light-speaking-loudly",
};
}

bool parse_qer_distribution(const std::string& name, QERDistribution& result) {
    const std::string value = lower(name);
    if (value == "jones" || value == "full" || value == "default") {
        result = {QERDistributionKind::Jones, {}};
        return true;
    }
    if (value == "fitted" || value == "jones-fitted" || value == "jones_fitted") {
        result = {QERDistributionKind::JonesFitted, {}};
        return true;
    }
    if (value == "mikszewski" || value == "mikszewski-sars-cov-2") {
        result = {QERDistributionKind::Mikszewski,
                  "mikszewski-sars-cov-2-standing-speaking"};
        return true;
    }
    if (std::find(profiles.begin(), profiles.end(), value) != profiles.end()) {
        result = {QERDistributionKind::Mikszewski, value};
        return true;
    }
    return false;
}

std::string qer_distribution_name(const QERDistribution& distribution) {
    if (distribution.kind == QERDistributionKind::Jones) return "jones";
    if (distribution.kind == QERDistributionKind::JonesFitted) return "jones-fitted";
    return distribution.profile;
}

std::vector<std::string> mikszewski_distribution_names() {
    return profiles;
}
