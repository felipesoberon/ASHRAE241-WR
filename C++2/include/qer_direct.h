// Work package 1: typed direct-QER family/profile selection and the
// immutable Jones-fitted and Mikszewski parameter tables.
//
// Source of truth: C++2/docs/QER_DISTRIBUTION_IMPLEMENTATION_PLAN.md Phase 2,
// and the legacy `qer_fitted_params`/`mikszewski_params` dictionaries in
// C++/src/qer_fitted.cpp (read-only reference; this module does not include
// or link against the legacy C++ tree).
//
// This module is intentionally self-contained: it defines the QER family and
// Mikszewski profile types, the two immutable parameter tables, their
// explicit legacy-short-to-canonical mappings, and the pure direct-QER
// log10-normal transform (built on distributions::log10normal_ppf). It does
// not touch the CLI, sampling plan, simulation loop, output schema, or seed
// policy -- those are later work packages per the plan.

#ifndef ECAI_CPP2_QER_DIRECT_H
#define ECAI_CPP2_QER_DIRECT_H

#include <array>
#include <string>

namespace ecai_cpp2 {
namespace qer_direct {

// The three QER calculation families from the plan. Jones (full mechanistic
// Jones et al. calculation) is the zero/default value; JonesFitted and
// Mikszewski are the two direct-QER families this module adds tables for.
enum class QerFamily {
    Jones = 0,
    JonesFitted,
    Mikszewski,
};

// Canonical lowercase-hyphenated rendering: "jones", "jones-fitted",
// "mikszewski".
std::string to_string(QerFamily family);

// Strict parse of the canonical rendering above. No aliases, no case
// folding. Throws std::invalid_argument for anything else, including the
// empty string.
QerFamily parse_qer_family(const std::string& name);

// (mu_log10, sigma_log10) of a direct-QER log10-normal distribution:
// QER_i = 10 ** (mu_log10 + sigma_log10 * normal_ppf(u)).
struct DirectQerParams {
    double mu_log10;
    double sigma_log10;
};

// ------------------------------------------------------------------
// Jones-fitted table: exactly 25 records, one per Jones et al. occupancy
// category (category_data.h publication order), copied exactly from
// C++/src/qer_fitted.cpp's `qer_fitted_params` map.
// ------------------------------------------------------------------

inline constexpr std::size_t kJonesFittedRecordCount = 25;

struct JonesFittedRecord {
    std::string canonical_category;  // Matches category_data.h canonical_name.
    DirectQerParams params;
};

// Identifies the fitted-parameter source and derivation method. Immutable
// provenance for output/audit use; not parsed or matched against anything.
extern const char* const kJonesFittedTableVersion;

// The 25 Jones-fitted records, in category_data.h publication order.
const std::array<JonesFittedRecord, kJonesFittedRecordCount>& jones_fitted_table();

// Strict lookup by C++2 canonical category name (category_data.h
// canonical_name, exact match, no normalization). Throws
// std::invalid_argument if `canonical_category_name` is not one of the 25
// canonical categories in jones_fitted_table().
const DirectQerParams& jones_fitted_params_for_canonical_category(
    const std::string& canonical_category_name);

// One legacy short name (C++/src/qer_fitted.cpp `qer_fitted_params` map key)
// mapped to its C++2 canonical category name.
struct LegacyJonesFittedAlias {
    std::string legacy_short_name;
    std::string canonical_category;
};

// The complete, explicit legacy-short-to-canonical mapping: exactly 25
// entries, one per canonical category in jones_fitted_table() (some legacy
// short names, e.g. "Cell", already equal their canonical category; others,
// e.g. "Food" -> "Food and beverage facilities", do not). No entry is
// omitted and no canonical category is left unmapped.
const std::array<LegacyJonesFittedAlias, kJonesFittedRecordCount>& jones_fitted_legacy_aliases();

// Strict lookup by legacy short name only (exact match against
// jones_fitted_legacy_aliases(), not against canonical names). Throws
// std::invalid_argument for any name not present in that mapping --
// including a canonical name that is not itself also a legacy short name --
// so there is no unknown-category or canonical-name fallback.
const DirectQerParams& jones_fitted_params_for_legacy_short_name(
    const std::string& legacy_short_name);

// ------------------------------------------------------------------
// Mikszewski table: exactly 30 profiles (10 organisms x 3 activities),
// copied exactly from C++/src/qer_fitted.cpp's `mikszewski_params` map.
// mu_log10 is computed as log10(median) at initialization, without
// rounding; sigma_log10 is the organism-specific value from that map.
// ------------------------------------------------------------------

enum class MikszewskiOrganism {
    SarsCov1,
    Mers,
    TbOnTreatment,
    Influenza,
    Coxsackievirus,
    Rhinovirus,
    SarsCov2,
    TbUntreated,
    Adenovirus,
    Measles,
};

enum class MikszewskiActivity {
    Resting,
    StandingSpeaking,
    LightSpeakingLoudly,
};

struct MikszewskiProfile {
    MikszewskiOrganism organism;
    MikszewskiActivity activity;
};

inline constexpr std::size_t kMikszewskiProfileCount = 30;

struct MikszewskiProfileRecord {
    MikszewskiProfile profile;
    double median_quanta_per_hour;  // Literal median, quanta/h, as printed.
    DirectQerParams params;         // mu_log10 = log10(median_quanta_per_hour).
};

extern const char* const kMikszewskiTableVersion;

// The 30 Mikszewski profile records, in
// C++/src/qer_fitted.cpp `mikszewski_params` map order (organism-major,
// resting/standing-speaking/light-speaking-loudly within each organism).
const std::array<MikszewskiProfileRecord, kMikszewskiProfileCount>& mikszewski_profile_table();

// Canonical name: "mikszewski-<organism>-<activity>", matching the legacy
// `mikszewski_params` map keys exactly (e.g.
// "mikszewski-sars-cov-2-standing-speaking").
std::string canonical_mikszewski_profile_name(MikszewskiProfile profile);

// Strict parse of the canonical name above, over all 30 profiles. No bare
// "mikszewski" default resolution (an explicit plan decision left to a later
// work package) and no partial/prefix match. Throws std::invalid_argument
// for anything not exactly one of the 30 canonical names.
MikszewskiProfile parse_mikszewski_profile(const std::string& canonical_name);

// Strict lookup by profile. Every MikszewskiProfile value constructible from
// the enums above and present in mikszewski_profile_table() is accepted;
// throws std::invalid_argument for an organism/activity combination not in
// the table (not reachable via parse_mikszewski_profile, which only ever
// returns table entries).
const DirectQerParams& mikszewski_params_for_profile(MikszewskiProfile profile);

// Convenience: parse_mikszewski_profile followed by mikszewski_params_for_profile.
const DirectQerParams& mikszewski_params_for_canonical_name(const std::string& canonical_name);

// ------------------------------------------------------------------
// Pure direct-QER log10-normal transform.
// ------------------------------------------------------------------

// QER_i = 10 ** (mu_log10 + sigma_log10 * normal_ppf(u)), computed via the
// existing distributions::log10normal_ppf. Pure function of (u, params):
// no RNG state, no side effects. Throws std::invalid_argument under exactly
// the same conditions as distributions::log10normal_ppf (u not strictly
// inside (0, 1); nonfinite mu_log10; nonfinite or non-positive sigma_log10).
double sample_direct_qer_log10normal(double u, const DirectQerParams& params);

}  // namespace qer_direct
}  // namespace ecai_cpp2

#endif  // ECAI_CPP2_QER_DIRECT_H
