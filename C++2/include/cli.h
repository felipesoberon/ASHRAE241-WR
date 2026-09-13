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
    // volume, mask efficiency, or input ECAi. Explicitly supplying a value
    // equal to the default still counts. Simulation count, replicate count,
    // base seed, confidence level, output mode and the output controls never
    // set this on their own.
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
// path string; and any numeric value that is not a strict, locale
// independent, full-string parse in range for its field.
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
