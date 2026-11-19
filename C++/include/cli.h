// Stage 12: strict command-line parsing and run resolution.
//
// Source of truth: IMPLEMENTATION_PLAN.md Section 7, Stage 12 and sections
// 5.1 (frozen defaults) and 5.2 (the exact flag list and parsing rules), and
// docs/MODEL_SPECIFICATION.md D7 (category naming, normalization and legacy
// aliases; no fuzzy matching).
//
// The module is split in two deliberately:
//
//   parse_cli   -- a pure function of the argument vector. It performs every
//                  syntactic and domain check that can be made without
//                  knowing the category, and runs no lookup, no allocation
//                  and no simulation.
//   resolve_run -- performs the strict category lookup and builds the
//                  simulation::EffectiveConfiguration and
//                  simulation::EffectiveOverrides a run needs, after
//                  re-deriving the parser's provenance and requiring it to
//                  agree.
//
// Dependency direction: this header includes simulation.h (and, through it,
// category_data.h); nothing under simulation/statistics includes this
// header, so no cycle exists.

#ifndef ECAI_CPP2_CLI_H
#define ECAI_CPP2_CLI_H

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "category_data.h"
#include "qer_direct.h"
#include "simulation.h"

namespace ecai_cpp2 {
namespace cli {

// Frozen defaults applied when the corresponding flag is omitted
// (IMPLEMENTATION_PLAN.md section 5.1). The four category scientific inputs
// (occupants, volume, mask efficiency, community rate) have no default here:
// they come from the selected occupancy category record.
inline constexpr std::uint64_t kDefaultSimulations = 10000;
inline constexpr std::uint64_t kDefaultReplicates = 500;
inline constexpr double kDefaultPercentile = 0.96;
inline constexpr double kDefaultConfidenceLevel = 0.95;
inline constexpr double kDefaultTargetProbability = 0.001;
inline constexpr double kDefaultDurationHours = 1.0;
inline constexpr std::uint64_t kDefaultBaseSeed = 485;

// The output-mode default (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md
// section 4.4). There is no default constant for --input-ecai: in
// probability mode it resolves from the selected category's published Table
// 4 ECAi (simulation::resolve_scenario), not from a frozen literal here.
inline constexpr simulation::OutputMode kDefaultMode = simulation::OutputMode::ecai;

// Every accepted long option, in the order of IMPLEMENTATION_PLAN.md section
// 5.2. Exposed so the help text and the parser cannot drift apart.
const std::vector<std::string>& accepted_flags();

// The complete current help text. `--help` prints exactly this.
std::string help_text();

// One fully parsed command line. Values are already domain-checked against
// everything knowable without the category record; provenance flags record
// which values were supplied explicitly, so Stage 13 can report source
// defaults, effective values, and custom-scenario status separately.
struct CliOptions {
    // True only for the sole-argument `--help` invocation. When set, every
    // other field is at its default and no category is required.
    bool help_requested = false;

    // The category string exactly as supplied; never normalized here.
    // Lookup happens in resolve_run.
    std::string category_name;
    bool category_supplied = false;

    // Run controls. Not scientific inputs: supplying them does not make a
    // run a custom scenario.
    std::uint64_t base_seed = kDefaultBaseSeed;
    std::uint64_t simulations_n = kDefaultSimulations;
    std::uint64_t replicates_r = kDefaultReplicates;
    double confidence_level = kDefaultConfidenceLevel;

    // Scientific run inputs with a frozen default (no category source).
    double percentile = kDefaultPercentile;
    double target_probability = kDefaultTargetProbability;
    // Nonempty only for paired multi-target output, in seed-significant order.
    std::vector<double> target_probabilities;
    std::optional<std::string> raw_output_directory;
    double duration_hours = kDefaultDurationHours;

    // The four category scientific inputs. An unset optional means "use the
    // occupancy-category default"; the category record itself is never
    // copied over.
    std::optional<long long> occupants;
    std::optional<double> volume_m3;
    std::optional<double> mask_efficiency;
    std::optional<double> community_infection_rate;

    // Output mode (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md section 4.4).
    // `mode` alone never marks a run custom scenario; an explicit
    // --input-ecai does. An unset input_ecai_l_s_person in probability mode
    // resolves to the selected category's published Table 4 ECAi.
    simulation::OutputMode mode = kDefaultMode;
    bool mode_explicit = false;
    std::optional<double> input_ecai_l_s_person;

    // QER family/profile selector (docs/QER_DISTRIBUTION_IMPLEMENTATION_PLAN.md
    // Phase 5, work package 2). One strict `--qer-distribution` selector
    // covers all three families: "jones" (default, mechanistic Jones et al.),
    // "jones-fitted", or a full canonical Mikszewski profile name (e.g.
    // "mikszewski-sars-cov-2-standing-speaking"). Bare "mikszewski" is
    // rejected: it has no default-profile resolution. An explicit selector,
    // including an explicit "jones", always marks the run a custom scenario,
    // because it changes the QER calculation mechanism itself. This field is
    // parsed and resolved to an immutable descriptor before family-aware
    // sampling, simulation, and output-schema provenance integration.
    qer_direct::QerFamily qer_family = qer_direct::QerFamily::Jones;
    bool qer_family_explicit = false;
    std::optional<qer_direct::MikszewskiProfile> qer_mikszewski_profile;

    // Infector policy (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md Phase 3):
    // "community" (default) or "risk-averse", case sensitive, no aliases.
    // Independent of `mode`; never marks the run a custom scenario by
    // itself (see EffectiveScenario::infector_policy in simulation.h).
    simulation::InfectorPolicy infector_policy = simulation::InfectorPolicy::community;
    bool infector_policy_explicit = false;

    // Output controls. Stage 12 parses and records these only: no directory
    // or file is created, opened, written or checked.
    std::optional<std::string> output_directory;
    std::optional<std::string> summary_file;
    std::optional<std::string> replicates_file;
    bool write_replicates = false;

    // Per-field explicit-value provenance. True iff the flag appeared on the
    // command line, even when the supplied value equals the default.
    bool base_seed_explicit = false;
    bool simulations_explicit = false;
    bool replicates_explicit = false;
    bool confidence_explicit = false;
    bool percentile_explicit = false;
    bool target_probability_explicit = false;
    bool duration_explicit = false;

    // True iff at least one explicit scientific override was supplied:
    // target probability, duration, percentile, community rate, occupants,
    // volume, mask efficiency, input ECAi, or QER family/profile. Explicitly
    // supplying a value equal to the default still counts. Simulation count,
    // replicate count, base seed, confidence level, output mode, infector
    // policy and the output controls never set this on their own.
    bool custom_scenario = false;

    // True iff the effective duration is exactly 1.0 h. The target
    // probability always applies over the effective duration; only an
    // exactly-1 h duration may additionally be described as hourly, and an
    // explicitly supplied 1.0 h remains hourly while still being custom.
    bool target_is_hourly = true;
};

// A parsed command line; the name used where the whole invocation, rather
// than the option values alone, is the subject.
using ParsedCommand = CliOptions;

// The immutable resolved QER family/profile descriptor for one run, built
// once in resolve_run from the parsed selector and (for Jones-fitted) the
// resolved category, entirely through the work package 1 qer_direct API.
// Every field is set once at construction and consumed by the family-aware
// sampling, simulation, and output provenance paths.
//
// `direct_params` is populated for JonesFitted and Mikszewski only; for the
// default Jones family it is left at its zero-initialized value and must not
// be read, since full-Jones mechanistic sampling does not use a direct
// log10-normal distribution at all.
//
// The direct parameters are used by family-aware sampling and simulation;
// output builders emit their immutable table/profile provenance once per
// summary and replicate metadata record.
struct ResolvedQerSelection {
    qer_direct::QerFamily family = qer_direct::QerFamily::Jones;
    std::optional<qer_direct::MikszewskiProfile> mikszewski_profile;

    // Canonical rendering: "jones", "jones-fitted", or the full Mikszewski
    // profile name (qer_direct::canonical_mikszewski_profile_name).
    std::string canonical_name;

    // Meaningful only when family != Jones; see the struct comment above.
    qer_direct::DirectQerParams direct_params{};
};

// One resolved run: the selected category record, retained rather than
// copied over, alongside the effective configuration/overrides/scenario
// derived from it and the parsed options.
struct ResolvedRun {
    // Points into the static all_occupancy_categories() table, so the
    // published source defaults stay available next to the effective values.
    const OccupancyCategory* category = nullptr;

    simulation::EffectiveConfiguration config;
    simulation::EffectiveOverrides overrides;
    simulation::EffectiveScenario scenario;
    ResolvedQerSelection qer_selection;

    std::uint64_t base_seed = kDefaultBaseSeed;
    bool custom_scenario = false;
    bool target_is_hourly = true;
};

// Parses `args` (the argument vector *without* the program name). Pure: no
// lookup, no allocation beyond the returned value, no simulation, no
// filesystem access.
//
// Throws std::invalid_argument, with a single-line message, for: an unknown
// flag; a positional argument; a repeated flag (including a repeated boolean
// flag); a missing value; a value that begins with "--" where a value is
// required; `--help` combined with any other argument; an empty category or
// path string; any numeric value that is not a strict, locale independent,
// full-string parse in range for its field; an --infector-policy value other
// than "community" or "risk-averse"; --infector-policy risk-averse
// combined with an explicit --community-rate <= 0 (the definitive check,
// covering the category-default rate too, runs later in resolve_run); and a
// --qer-distribution value other than "jones", "jones-fitted", or a full
// canonical Mikszewski profile name (bare "mikszewski" is rejected).
CliOptions parse_cli(const std::vector<std::string>& args);

// Performs the strict category lookup (docs/MODEL_SPECIFICATION.md D7) and
// builds the effective configuration, overrides and scenario for `options`.
//
// Explicit CLI overrides win directly over the selected category defaults.
// The parser's provenance is re-derived here and must agree with the
// resulting configuration exactly; a disagreement is a defect and throws
// rather than running.
//
// Throws std::invalid_argument if `options` requests help, if no category
// was supplied, if the category does not match a canonical name or an
// enumerated legacy alias, or if any resolved value fails the Stage 10/11
// domain checks.
ResolvedRun resolve_run(const CliOptions& options);

}  // namespace cli
}  // namespace ecai_cpp2

#endif  // ECAI_CPP2_CLI_H
