// Stage 4: deterministic per-replicate/per-stream seed derivation.
//
// Source of truth: docs/MODEL_SPECIFICATION.md section 11.3 (seed
// architecture: SplitMix64 mixing + canonical configuration fingerprint) and
// docs/SEED_POLICY.md (this module's full policy document).
// IMPLEMENTATION_PLAN.md Section 7, Stage 4.

#ifndef ECAI_CPP2_SEED_MANAGER_H
#define ECAI_CPP2_SEED_MANAGER_H

#include <array>
#include <cstdint>
#include <string>

namespace ecai_cpp2 {
namespace seed_manager {

// One fixed identifier per named stochastic parameter stream
// (MODEL_SPECIFICATION.md section 11.3). Values are explicit and frozen:
// assigned once, never renumbered, never reused for a different parameter,
// even if a parameter is later removed.
enum class StreamId : std::uint64_t {
    infector_count = 0,
    susceptible_pbr = 1,
    gamma = 2,
    lambda = 3,
    infector_pbr = 4,
    c_drop = 5,
    d_bar = 6,
    E = 7,
    GVL = 8,
    VF = 9,
    RTD = 10,
    DK = 11,
};

// Exhaustive, ordered list of every currently defined StreamId. Used by
// completeness/uniqueness tests and by nothing else; this is not a
// dispatch table.
inline constexpr std::array<StreamId, 12> kAllStreamIds = {
    StreamId::infector_count, StreamId::susceptible_pbr, StreamId::gamma,
    StreamId::lambda,         StreamId::infector_pbr,    StreamId::c_drop,
    StreamId::d_bar,          StreamId::E,               StreamId::GVL,
    StreamId::VF,             StreamId::RTD,             StreamId::DK,
};
static_assert(kAllStreamIds.size() == static_cast<std::size_t>(StreamId::DK) + 1,
              "kAllStreamIds must remain exhaustive for contiguous StreamId values");

// Human-readable name for a stream id (its enumerator spelling). Throws
// std::invalid_argument for any value not in kAllStreamIds, including a
// value produced by casting an out-of-range integer to StreamId.
std::string stream_id_name(StreamId id);

// Typed effective-configuration fingerprint input
// (MODEL_SPECIFICATION.md section 11.3): the ordered, fixed set of fields
// that make two runs' derived seeds distinguishable. Field order here is
// normative and matches canonical_config_string's serialization order.
struct EffectiveConfig {
    std::string category_name;
    std::int64_t i0 = 0;
    double volume_m3 = 0.0;
    double mask_efficiency = 0.0;
    double community_infection_rate = 0.0;
    double target_probability = 0.0;
    double duration_hours = 0.0;
    double percentile = 0.0;
    double confidence_level = 0.0;
    std::uint64_t simulations_n = 0;
    std::uint64_t replicates_r = 0;
};

// Current C++2 sampling-design version. Part of the seed-derivation tuple so
// a future sampling-design change cannot silently collide with seeds
// produced by an earlier design.
inline constexpr std::uint64_t kDesignVersion = 1;

// Human-readable version of the complete seed derivation policy. Output
// metadata consumes this owner directly; serializers must not duplicate it.
inline constexpr const char* kSeedPolicyVersion = "CPP2-SEED-v1";

// ASCII domain tag "CPP2SEED", XORed into the base seed before the fold
// (MODEL_SPECIFICATION.md section 11.3).
inline constexpr std::uint64_t kDomainTag = 0x4350503253454544ULL;

// 64-bit FNV-1a constants (MODEL_SPECIFICATION.md section 11.3).
inline constexpr std::uint64_t kFnvOffsetBasis64 = 0xcbf29ce484222325ULL;
inline constexpr std::uint64_t kFnvPrime64 = 0x100000001b3ULL;

// 64-bit FNV-1a over the raw bytes of `data`. No std::hash anywhere in this
// module: std::hash<std::string>'s output is not guaranteed stable across
// standard library versions/builds (see docs/SEED_POLICY.md).
std::uint64_t fnv1a_64(const std::string& data);

// The exact SplitMix64 value transformation (Vigna), including the initial
// addition of the golden-gamma constant, applied to `x` as a pure function
// of one input (docs/SEED_POLICY.md and MODEL_SPECIFICATION.md section
// 11.3). This is not a persistent stream generator: calling it twice with
// the same `x` returns the same `z`.
std::uint64_t splitmix64(std::uint64_t x);

// Canonical, fixed-order, ASCII, '|'-delimited serialization of `config`:
// category name, I0, VOL, mask efficiency, community rate, target
// probability, duration, percentile, confidence level, simulations N,
// replicates R. Integers are formatted decimal; doubles are formatted
// equivalent to C's `%.17g` in the C locale via std::to_chars (locale-
// independent by construction), with negative zero normalized to positive
// zero before formatting. Throws std::invalid_argument if any double field
// is not finite, or if category_name contains the '|' delimiter. Full
// policy: docs/SEED_POLICY.md.
std::string canonical_config_string(const EffectiveConfig& config);

// 64-bit FNV-1a hash of canonical_config_string(config). Throws under the
// same conditions as canonical_config_string.
std::uint64_t config_fingerprint(const EffectiveConfig& config);

// The normative ordered SplitMix64 fold (MODEL_SPECIFICATION.md section
// 11.3): state = splitmix64(base_seed ^ kDomainTag); then for each element
// of (design_version, fingerprint, replicate_index, stream_id, list_length),
// in that exact order, state = splitmix64(state ^ splitmix64(element));
// the final state is the returned stream seed. Exposed directly (in
// addition to SeedManager::derive_seed) so the fold has independent,
// config-free known-answer test coverage. No alternate "xor or additive"
// choice is permitted; this fold is normative.
std::uint64_t fold_stream_seed(std::uint64_t base_seed, std::uint64_t design_version,
                                std::uint64_t fingerprint, std::uint64_t replicate_index,
                                std::uint64_t stream_id, std::uint64_t list_length);

// Derives one stream seed per (base_seed, design_version,
// effective_configuration_fingerprint, replicate_index, stream_id,
// list_length) tuple (IMPLEMENTATION_PLAN.md section 6). Centralizes
// derivation so no component may improvise `seed + index` independently.
class SeedManager {
public:
    explicit SeedManager(std::uint64_t base_seed, std::uint64_t design_version = kDesignVersion)
        : base_seed_(base_seed), design_version_(design_version) {}

    // Throws std::invalid_argument if `stream_id` is not in kAllStreamIds,
    // or if `config` is invalid per canonical_config_string.
    std::uint64_t derive_seed(const EffectiveConfig& config, std::uint64_t replicate_index,
                               StreamId stream_id, std::uint64_t list_length) const;

private:
    std::uint64_t base_seed_;
    std::uint64_t design_version_;
};

}  // namespace seed_manager
}  // namespace ecai_cpp2

#endif  // ECAI_CPP2_SEED_MANAGER_H
