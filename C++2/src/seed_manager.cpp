// Stage 4: deterministic per-replicate/per-stream seed derivation.
//
// Source of truth: docs/MODEL_SPECIFICATION.md section 11.3.
// IMPLEMENTATION_PLAN.md Section 7, Stage 4.

#include "seed_manager.h"

#include <charconv>
#include <cmath>
#include <stdexcept>

namespace ecai_cpp2 {
namespace seed_manager {

std::string stream_id_name(StreamId id) {
    switch (id) {
        case StreamId::infector_count:
            return "infector_count";
        case StreamId::susceptible_pbr:
            return "susceptible_pbr";
        case StreamId::gamma:
            return "gamma";
        case StreamId::lambda:
            return "lambda";
        case StreamId::infector_pbr:
            return "infector_pbr";
        case StreamId::c_drop:
            return "c_drop";
        case StreamId::d_bar:
            return "d_bar";
        case StreamId::E:
            return "E";
        case StreamId::GVL:
            return "GVL";
        case StreamId::VF:
            return "VF";
        case StreamId::RTD:
            return "RTD";
        case StreamId::DK:
            return "DK";
    }
    throw std::invalid_argument("seed_manager::stream_id_name: unrecognized StreamId");
}

std::uint64_t fnv1a_64(const std::string& data) {
    std::uint64_t hash = kFnvOffsetBasis64;
    for (unsigned char byte : data) {
        hash ^= static_cast<std::uint64_t>(byte);
        hash *= kFnvPrime64;
    }
    return hash;
}

std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    std::uint64_t z = x;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z = z ^ (z >> 31);
    return z;
}

namespace {

// Formats one canonical-serialization double field: rejects non-finite
// values, normalizes negative zero to positive zero, then formats
// equivalent to C's `%.17g` via std::to_chars, which never consults the
// current locale (docs/SEED_POLICY.md).
std::string format_double_canonical(double value) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(
            "seed_manager: effective configuration double field must be finite");
    }
    if (value == 0.0) {
        value = 0.0;  // -0.0 and +0.0 both canonicalize to +0.0
    }

    char buffer[64];
    const auto result =
        std::to_chars(buffer, buffer + sizeof(buffer), value, std::chars_format::general, 17);
    if (result.ec != std::errc()) {
        throw std::invalid_argument(
            "seed_manager: failed to format effective configuration double field");
    }
    return std::string(buffer, result.ptr);
}

}  // namespace

std::string canonical_config_string(const EffectiveConfig& config) {
    if (config.category_name.find('|') != std::string::npos) {
        throw std::invalid_argument(
            "seed_manager: category_name must not contain the '|' delimiter");
    }

    std::string out;
    out.reserve(128);
    out += config.category_name;
    out += '|';
    out += std::to_string(config.i0);
    out += '|';
    out += format_double_canonical(config.volume_m3);
    out += '|';
    out += format_double_canonical(config.mask_efficiency);
    out += '|';
    out += format_double_canonical(config.community_infection_rate);
    out += '|';
    out += format_double_canonical(config.target_probability);
    out += '|';
    out += format_double_canonical(config.duration_hours);
    out += '|';
    out += format_double_canonical(config.percentile);
    out += '|';
    out += format_double_canonical(config.confidence_level);
    out += '|';
    out += std::to_string(config.simulations_n);
    out += '|';
    out += std::to_string(config.replicates_r);
    return out;
}

std::uint64_t config_fingerprint(const EffectiveConfig& config) {
    return fnv1a_64(canonical_config_string(config));
}

std::uint64_t fold_stream_seed(std::uint64_t base_seed, std::uint64_t design_version,
                                std::uint64_t fingerprint, std::uint64_t replicate_index,
                                std::uint64_t stream_id, std::uint64_t list_length) {
    std::uint64_t state = splitmix64(base_seed ^ kDomainTag);
    const std::uint64_t elements[5] = {design_version, fingerprint, replicate_index, stream_id,
                                        list_length};
    for (std::uint64_t element : elements) {
        state = splitmix64(state ^ splitmix64(element));
    }
    return state;
}

std::uint64_t SeedManager::derive_seed(const EffectiveConfig& config,
                                        std::uint64_t replicate_index, StreamId stream_id,
                                        std::uint64_t list_length) const {
    // Throws std::invalid_argument for any stream_id not in kAllStreamIds.
    stream_id_name(stream_id);
    const std::uint64_t fingerprint = config_fingerprint(config);
    return fold_stream_seed(base_seed_, design_version_, fingerprint, replicate_index,
                             static_cast<std::uint64_t>(stream_id), list_length);
}

}  // namespace seed_manager
}  // namespace ecai_cpp2
