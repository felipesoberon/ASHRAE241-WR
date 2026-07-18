#include "qer_distribution.h"

#include <algorithm>
#include <cctype>

namespace {
std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}
}

bool parse_qer_distribution(const std::string& name, QERDistribution& result) {
    const std::string value = lower(name);
    if (value == "jones" || value == "full" || value == "default") {
        result = QERDistribution::Jones;
        return true;
    }
    if (value == "fitted" || value == "jones-fitted" || value == "jones_fitted") {
        result = QERDistribution::JonesFitted;
        return true;
    }
    if (value == "mikszewski") {
        result = QERDistribution::Mikszewski;
        return true;
    }
    return false;
}

const char* qer_distribution_name(QERDistribution distribution) {
    switch (distribution) {
    case QERDistribution::Jones: return "jones";
    case QERDistribution::JonesFitted: return "jones-fitted";
    case QERDistribution::Mikszewski: return "mikszewski";
    }
    return "jones";
}
