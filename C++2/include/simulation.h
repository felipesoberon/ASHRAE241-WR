// Stage 10: one replicate, end to end.
//
// Source of truth: docs/MODEL_SPECIFICATION.md section 3 (QER/ECAi), section
// 11.1 (person-specific sampling and the zero-infector shortcut, D3/D4),
// section 11.4 (resource budget), D8 and section 12 (the frozen linear
// percentile-interpolation convention). IMPLEMENTATION_PLAN.md Section 7,
// Stage 10 and section 4.3/4.5.
//
// This module is orchestration only: it owns no equation and no
// distribution. It combines Stage 6 (build_sampling_plan), Stage 7
// (build_parameter_banks + assign_replicate), Stage 8 (compute_qer), Stage 9
// (compute_ecai) and the frozen percentile helper below into one replicate,
// emitting exactly N ECAi values in simulation-index order.

#ifndef ECAI_CPP2_SIMULATION_H
#define ECAI_CPP2_SIMULATION_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "category_data.h"
#include "sampling_plan.h"
#include "seed_manager.h"
#include "statistics.h"

namespace ecai_cpp2 {
namespace simulation {

// The frozen finite-sample percentile convention (D8 / section 12): linear
// interpolation matching numpy.percentile's default "linear" (Hyndman-Fan
// Type 7) method. rank = p * (n - 1); the result interpolates between the
// order statistics at floor(rank) and ceil(rank).
//
// `values` is taken by value-independent read: the function sorts an
// internal copy and never reorders or otherwise mutates the caller's
// index-ordered vector. p is accepted over the closed interval [0, 1]
// (p = 0 yields the minimum, p = 1 the maximum). n == 1 yields that single
// value for every accepted p. No rounding or ceiling is applied here; Stage
// 11/12 own the ceiling-to-5 comparator.
double percentile_linear(const std::vector<double>& values, double p);

// The output mode a run resolves to: the frozen ECAi inversion (default), or
// the forward infection-probability evaluation at a fixed input ECAi
// (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md section 4.2). `ecai` is
// enumerator 0, so a value-initialized EffectiveScenario/EffectiveOverrides
// defaults to ECAi mode without any explicit initialization.
enum class OutputMode { ecai = 0, probability = 1 };
static_assert(static_cast<int>(OutputMode::ecai) == 0,
              "OutputMode::ecai must be enumerator 0 so value-initialization defaults to it");

// The infector-handling policy a run resolves to
// (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md): the current unconditional
// binomial infector count (default, including zero-infector rows), or a
// risk-averse policy that conditions every simulation on at least one
// infector. Independent of OutputMode: either policy may be combined with
// either output mode. `community` is enumerator 0, so a value-initialized
// EffectiveScenario/EffectiveOverrides defaults to the unchanged community
// behavior without any explicit initialization. Phase 3 (this header/CLI)
// resolves and validates the policy; the conditional sampling and simulation
// integration that consumes it are added in later phases.
enum class InfectorPolicy { community = 0, risk_averse = 1 };
static_assert(static_cast<int>(InfectorPolicy::community) == 0,
              "InfectorPolicy::community must be enumerator 0 so value-initialization defaults "
              "to it");

// The fully resolved per-run scenario used by one replicate: category
// defaults after any explicit override, plus the run-level scientific
// inputs. Units follow model::ECAiInputs (volume m^3, duration h, everything
// else dimensionless or a count).
struct EffectiveScenario {
    std::string category_name;
    long long occupants = 0;  // I0
    double volume_m3 = 0.0;
    double mask_efficiency = 0.0;
    double community_infection_rate = 0.0;  // binomial p for the Stage 6 plan
    double target_probability = 0.0;
    double duration_hours = 0.0;
    double percentile = 0.0;
    std::size_t simulations_n = 0;

    // Per-field provenance: true when the value above came from an explicit
    // override rather than the occupancy-category record.
    bool occupants_overridden = false;
    bool volume_overridden = false;
    bool mask_efficiency_overridden = false;
    bool community_infection_rate_overridden = false;
    // True iff at least one of the four flags above is set
    // (IMPLEMENTATION_PLAN.md section 5.2: any override of a category
    // scientific input makes the run a custom scenario). Deliberately not
    // ORed with input_ecai_overridden below: this keeps its ECAi-mode
    // meaning ("one of the four category scientific inputs was overridden")
    // exactly as it was before probability mode existed.
    bool custom_scenario = false;

    // Probability-mode resolution (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md
    // section 4.2/7). In ECAi mode both stay at their defaults regardless of
    // any --input-ecai value (the CLI rejects that combination before
    // resolution). In probability mode input_ecai_l_s_person is either the
    // explicit override or the selected category's published Table 4 ECAi.
    OutputMode mode = OutputMode::ecai;
    double input_ecai_l_s_person = 0.0;
    bool input_ecai_overridden = false;

    // Infector policy (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md Phase 3),
    // resolved once here and never re-derived from a raw CLI flag by any
    // downstream consumer. `infector_policy_explicit` records whether
    // --infector-policy appeared on the command line; it is provenance only
    // and, like `mode`, is deliberately not ORed into `custom_scenario`
    // above -- the policy is an estimand change, not a four-field scientific
    // override (Phase 5).
    InfectorPolicy infector_policy = InfectorPolicy::community;
    bool infector_policy_explicit = false;
};

// The four category scientific inputs a run may override. An unset optional
// means "use the occupancy-category default". Target probability, duration
// and percentile have no category default and are always taken from the
// effective configuration.
struct EffectiveOverrides {
    std::optional<long long> occupants;
    std::optional<double> volume_m3;
    std::optional<double> mask_efficiency;
    std::optional<double> community_infection_rate;

    // Probability-mode overrides. `mode` selects the output mode (default
    // ECAi); `input_ecai_l_s_person` is only consulted in probability mode,
    // where an unset optional resolves to the category's published ECAi.
    OutputMode mode = OutputMode::ecai;
    std::optional<double> input_ecai_l_s_person;

    // Infector policy override (default community) plus its explicit-value
    // provenance; see EffectiveScenario::infector_policy above.
    InfectorPolicy infector_policy = InfectorPolicy::community;
    bool infector_policy_explicit = false;
};

// Resolves the category record plus explicit overrides into one scenario,
// and verifies that `config` -- the fingerprint the seeds are derived from --
// agrees with the result exactly. Runs before any allocation or sampling.
//
// Throws std::invalid_argument if:
//   - config.category_name is not exactly category.canonical_name;
//   - any of I0 / VOL / mask efficiency / community rate in `config` differs
//     from the resolved (category default or overridden) value;
//   - config.simulations_n is zero, or exceeds std::size_t;
//   - I0 < 1, VOL not finite and positive, mask efficiency outside [0, 1],
//     community rate outside [0, 1], target probability outside the open
//     interval (0, 1), duration not finite and positive, or percentile not
//     finite and inside [0, 1];
//   - overrides.mode is probability and the resolved input ECAi (the
//     override, or else category.published_ecai_Lsp) is not finite or is
//     negative;
//   - overrides.infector_policy is not a declared InfectorPolicy enumerator,
//     or it is risk_averse and the resolved community_infection_rate is not
//     strictly positive (I0 < 1 is already rejected unconditionally above).
//     This authoritative check runs here, before any LHS generation or
//     allocation, and applies identically whether the caller went through
//     the CLI or constructed `overrides` directly.
// Mode, the resolved input ECAi and the infector policy are deliberately not
// part of `config`: they are not fingerprint fields
// (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md section 6;
// docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md Phase 4/7), so this function reads
// them only from `overrides`.
EffectiveScenario resolve_scenario(const OccupancyCategory& category,
                                    const seed_manager::EffectiveConfig& config,
                                    const EffectiveOverrides& overrides = EffectiveOverrides{});

// The order in which already-assigned active rows are handed to the model.
// This is a computational schedule only: assignments are materialized in
// stable simulation index order first, so no order below can change a
// person-to-value mapping, an output position, a counter, or the percentile.
// Infectors within a row are never reordered, so a row's QER sum is
// bit-identical under every schedule.
enum class EvaluationOrder {
    forward = 0,      // active rows in ascending simulation index
    reverse = 1,      // active rows in descending simulation index
    interleaved = 2,  // first, last, second, second-to-last, ...
};

// Per-replicate options. `policy` is forwarded unchanged to the Stage 6/7
// pre-allocation checks; `collect_audit` turns on the (otherwise unbuilt)
// consumption audit below.
struct ReplicateOptions {
    sampling_plan::ResourcePolicy policy{};
    EvaluationOrder evaluation_order = EvaluationOrder::forward;
    bool collect_audit = false;
};

// Public run configuration for multi-replicate orchestration. It is the same
// canonical seed fingerprint input used by Stage 10, with the Stage 11
// structural defaults applied without executing a batch.
struct EffectiveConfiguration : seed_manager::EffectiveConfig {
    EffectiveConfiguration() {
        confidence_level = 0.95;
        replicates_r = 500;
    }
};

// Batch orchestration forwards the existing per-replicate options unchanged.
using SimulationOptions = ReplicateOptions;

// Explicit, separately tracked counts for one replicate. zero_infector_count
// and negative_clamp_count are deliberately distinct: output_zero_count may
// include both (and, in principle, an active row whose raw ECAi is exactly
// zero without clamping).
struct ReplicateCounters {
    std::size_t n = 0;                    // emitted ECAi values
    std::size_t zero_infector_count = 0;  // rows short-circuited to +0.0
    std::size_t active_count = 0;         // rows with at least one infector
    std::size_t k_total = 0;              // realized infectors in the replicate
    std::size_t negative_clamp_count = 0;  // active rows whose raw ECAi was negative
    std::size_t output_zero_count = 0;     // outputs exactly equal to zero
    std::size_t qer_evaluations = 0;       // one per realized infector
    std::size_t ecai_evaluations = 0;      // one per active row
    std::size_t susceptible_values_consumed = 0;  // 3 * active_count
    std::size_t infector_values_consumed = 0;     // 8 * k_total

    // Probability-mode-only counters (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md
    // section 4.2). Both stay zero in ECAi mode, exactly as negative_clamp_count
    // and ecai_evaluations stay zero in probability mode.
    std::size_t probability_evaluations = 0;  // one per active row, probability mode only
    std::size_t saturated_count = 0;          // active rows whose probability == 1.0
};

// One simulation's audited identity, recorded from the materialized
// assignment and therefore independent of the evaluation schedule.
struct RowAudit {
    std::size_t simulation_index = 0;
    long long infector_count = 0;
    bool active = false;
    double qer_sum_per_h = 0.0;
    double ecai_raw_l_s_person = 0.0;
    double ecai_l_s_person = 0.0;
    bool clamped = false;
    // Probability-mode-only fields, populated instead of the two above when
    // the replicate ran in probability mode; left at their defaults (0.0 /
    // false) in ECAi mode.
    double probability = 0.0;
    bool saturated = false;
    // {susceptible_pbr, gamma, lambda}; all zero for an inactive row, which
    // consumes nothing.
    std::array<std::uint64_t, 3> susceptible_ids{};
    // One entry per realized infector, in infector-list order; each holds
    // that infector's own eight ids in
    // {infector_pbr, c_drop, d_bar, E, GVL, VF, RTD, DK} order.
    std::vector<std::array<std::uint64_t, 8>> infector_ids;
};

// Consumption audit: every value handed out by the Stage 7 cursor, in
// consumption order, identified by sampling_plan::consumption_id, plus the
// per-simulation identity rows and the schedule actually used. Built only
// when ReplicateOptions::collect_audit is true.
struct ReplicateAudit {
    std::vector<std::uint64_t> consumption_ids;
    // Exactly N entries, in stable simulation index order.
    std::vector<RowAudit> rows;
    // The active rows' simulation indices in the order they were evaluated.
    std::vector<std::size_t> evaluation_order;
};

// One replicate's complete result. Deliberately carries no ParameterBanks
// and no cursor: once the traversal is finished the banks are destroyed, so
// no consumer can re-read, refill, or re-traverse them.
struct ReplicateResult {
    // Stable top-level batch identity. The injected from-banks fixture path
    // intentionally defaults to zero because it has no replicate-index input.
    std::uint64_t replicate_index = 0;
    // Exactly N values, in original simulation index order (never sorted).
    // Authoritative Stage 6 realized infector counts in the same stable order.
    std::vector<long long> infector_counts;
    // The resolved infector policy this replicate ran under
    // (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md Phase 5), copied verbatim from
    // the scenario that produced this replicate -- not re-derived from the
    // plan -- so a caller can reconcile scenario/plan/result agreement (and,
    // in a batch, detect a mixed-policy run) without re-reading the plan.
    // `community` (enumerator 0) is the default, matching InfectorPolicy's
    // own zero-valued default.
    InfectorPolicy infector_policy = InfectorPolicy::community;
    // The resolved mode this replicate ran in (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md
    // section 4.2). Determines which of the two pairs below is populated;
    // the other stays at its default (empty vector, 0.0 percentile).
    OutputMode mode = OutputMode::ecai;
    std::vector<double> ecai_l_s_person;
    double percentile = 0.0;                // the p that was applied
    double percentile_l_s_person = 0.0;     // raw, unrounded; ECAi mode only
    std::vector<double> probability;        // probability mode only
    double percentile_probability = 0.0;    // raw, unrounded; probability mode only
    ReplicateCounters counters;
    // Indexed by BankId ordinal.
    std::array<std::size_t, 11> planned_bank_lengths{};
    std::array<std::size_t, 11> final_cursors{};
    // Whether each bank was actually generated (a zero-length bank derives no
    // seed and draws no LHS list). Every entry is false for an all-zero
    // replicate: no non-binomial value is generated at all.
    std::array<bool, 11> bank_generated{};
    bool cursors_exhausted = false;
    bool audit_collected = false;
    ReplicateAudit audit;
};

// Complete, index-ordered replicate results and the summary computed from
// their raw (unrounded and unadjusted) percentile estimates.
struct ReplicateBatch {
    std::vector<ReplicateResult> replicates;
    statistics::ReplicateSummary summary;
};

// Runs one replicate over an already-built plan and already-built parameter
// banks (the injectable path: both are supplied explicitly, so a fixture can
// pin every sampled value).
//
// Order of operations is normative:
//   1. every zero-infector row's output is set to exactly +0.0 and counted,
//      before any susceptible/infector assignment or model input is
//      requested (section 11.1, D3);
//   2. sampling_plan::assign_replicate materializes the whole
//      simulation/person -> value mapping, in stable simulation index order,
//      before any equation runs;
//   3. every realized infector's QER is computed, summed in infector-list
//      order, and each active row's ECAi is written back to its own stable
//      simulation index;
//   4. the frozen percentile is taken from a sorted copy.
// A successful return proves every cursor was exactly exhausted.
ReplicateResult run_replicate_from_banks(const sampling_plan::SamplingPlan& plan,
                                          const sampling_plan::ParameterBanks& banks,
                                          const EffectiveScenario& scenario,
                                          const ReplicateOptions& options = ReplicateOptions{});

// Runs one replicate over an explicitly supplied plan (an injected `ninf`
// vector), building this replicate's parameter banks from it and delegating
// to run_replicate_from_banks. The banks are local to this call and are
// destroyed before it returns, so nothing can re-read or refill them.
ReplicateResult run_replicate_from_plan(const seed_manager::SeedManager& seeds,
                                         const seed_manager::EffectiveConfig& config,
                                         std::uint64_t replicate_index,
                                         const sampling_plan::SamplingPlan& plan,
                                         const OccupancyCategory& category,
                                         const EffectiveScenario& scenario,
                                         const ReplicateOptions& options = ReplicateOptions{});

// The full Stage 10 entry point: one replicate from an explicit base seed,
// effective configuration, occupancy category, effective overrides and
// resource policy.
//
//   resolve_scenario -> build_sampling_plan -> build_parameter_banks ->
//   assign_replicate -> compute_qer -> compute_ecai -> N index-ordered
//   outputs -> percentile_linear
//
// Seeds are derived only through seed_manager::SeedManager(base_seed); no
// component improvises a seed. Nothing persists between calls, so two calls
// with the same arguments are bit-identical and no state can leak from an
// intervening call with different arguments.
ReplicateResult run_replicate(std::uint64_t base_seed,
                               const seed_manager::EffectiveConfig& config,
                               std::uint64_t replicate_index, const OccupancyCategory& category,
                               const EffectiveOverrides& overrides = EffectiveOverrides{},
                               const ReplicateOptions& options = ReplicateOptions{});

// Runs replicates serially at indices [0, R), preserving each complete Stage
// 10 result and appending each result's raw (unrounded, unadjusted)
// percentile estimate to raw_percentiles_out, in replicate-index order. R and
// confidence are preflighted before any replicate runs or batch allocation.
// ResourcePolicy remains scoped to each run_replicate call.
//
// This is the Stage 10-only half of run_replicates below: it calls no Stage
// 11 statistics function, so it links against the Stage 3-10 sources alone
// (see ecai_cpp2_simulation_tests in CMakeLists.txt, which -- unlike
// ecai_cpp2_statistics_tests -- does not link src/statistics.cpp). Only
// run_replicates itself, which turns raw_percentiles_out into a summary via
// statistics::summarize_replicates, requires that Stage 11 link.
inline ReplicateBatch run_replicate_batch(std::uint64_t base_seed, EffectiveConfiguration config,
                                           const OccupancyCategory& category,
                                           const EffectiveOverrides& overrides,
                                           const SimulationOptions& options,
                                           std::vector<double>& raw_percentiles_out) {
    // Batch-only inputs are validated before scenario resolution, allocation,
    // or stochastic work. The ordered comparisons also reject NaN.
    if (!(config.confidence_level > 0.0) || !(config.confidence_level < 1.0)) {
        throw std::invalid_argument(
            "simulation: confidence_level must lie strictly inside (0, 1)");
    }
    if (config.replicates_r == 0) {
        throw std::invalid_argument("simulation: replicates R must be at least 1");
    }
    if (config.replicates_r > statistics::kMaxSupportedReplicateCount) {
        throw std::invalid_argument(
            "simulation: replicates R exceeds the supported Student-t range");
    }
    if constexpr (std::numeric_limits<std::size_t>::max() <
                  std::numeric_limits<std::uint64_t>::max()) {
        if (config.replicates_r >
            static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
            throw std::invalid_argument("simulation: replicates R exceeds size_t");
        }
    }
    const std::size_t replicate_count = static_cast<std::size_t>(config.replicates_r);

    ReplicateBatch batch;
    if (replicate_count > batch.replicates.max_size() ||
        replicate_count > raw_percentiles_out.max_size()) {
        throw std::invalid_argument("simulation: replicates R exceeds vector max_size");
    }
    try {
        batch.replicates.reserve(replicate_count);
        raw_percentiles_out.reserve(replicate_count);
    } catch (const std::bad_alloc&) {
        throw std::runtime_error("simulation: unable to allocate replicate batch storage");
    }

    // All replicates in one batch must resolve to the same mode and the same
    // infector policy: neither is a fingerprint field, so nothing else
    // enforces this, and this is the first point that observes every
    // replicate's resolved mode/policy together
    // (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md section 4.2;
    // docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md Phase 4's "batch-consistency
    // check rejecting mixed policies").
    std::optional<OutputMode> batch_mode;
    std::optional<InfectorPolicy> batch_infector_policy;
    for (std::size_t index = 0; index < replicate_count; ++index) {
        ReplicateResult result = run_replicate(base_seed, config,
                                               static_cast<std::uint64_t>(index), category,
                                               overrides, options);
        if (!batch_mode.has_value()) {
            batch_mode = result.mode;
        } else if (*batch_mode != result.mode) {
            throw std::runtime_error(
                "simulation: replicates within one batch resolved to different output modes");
        }
        if (!batch_infector_policy.has_value()) {
            batch_infector_policy = result.infector_policy;
        } else if (*batch_infector_policy != result.infector_policy) {
            throw std::runtime_error(
                "simulation: replicates within one batch resolved to different infector "
                "policies");
        }
        raw_percentiles_out.push_back(result.mode == OutputMode::probability
                                           ? result.percentile_probability
                                           : result.percentile_l_s_person);
        batch.replicates.push_back(std::move(result));
    }
    return batch;
}

// Runs replicates serially at indices [0, R), preserving each complete Stage
// 10 result, then summarizes their raw percentile estimates via
// statistics::summarize_replicates. Thin wrapper over run_replicate_batch
// above; see that function for the preflighting and per-replicate work this
// delegates.
inline ReplicateBatch run_replicates(
    std::uint64_t base_seed, EffectiveConfiguration config, const OccupancyCategory& category,
    const EffectiveOverrides& overrides = EffectiveOverrides{},
    const SimulationOptions& options = SimulationOptions{}) {
    std::vector<double> raw_percentiles;
    ReplicateBatch batch =
        run_replicate_batch(base_seed, config, category, overrides, options, raw_percentiles);
    batch.summary = statistics::summarize_replicates(raw_percentiles, config.confidence_level);
    return batch;
}

}  // namespace simulation
}  // namespace ecai_cpp2

#endif  // ECAI_CPP2_SIMULATION_H
