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
    // scientific input makes the run a custom scenario).
    bool custom_scenario = false;
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
//     finite and inside [0, 1].
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
    std::vector<double> ecai_l_s_person;
    double percentile = 0.0;                // the p that was applied
    double percentile_l_s_person = 0.0;     // raw, unrounded
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
// 10 result. R and confidence are preflighted before any replicate runs or
// batch allocation. ResourcePolicy remains scoped to each run_replicate call.
inline ReplicateBatch run_replicates(
    std::uint64_t base_seed, EffectiveConfiguration config, const OccupancyCategory& category,
    const EffectiveOverrides& overrides = EffectiveOverrides{},
    const SimulationOptions& options = SimulationOptions{}) {
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
    std::vector<double> raw_percentiles;
    if (replicate_count > batch.replicates.max_size() ||
        replicate_count > raw_percentiles.max_size()) {
        throw std::invalid_argument("simulation: replicates R exceeds vector max_size");
    }
    try {
        batch.replicates.reserve(replicate_count);
        raw_percentiles.reserve(replicate_count);
    } catch (const std::bad_alloc&) {
        throw std::runtime_error("simulation: unable to allocate replicate batch storage");
    }

    for (std::size_t index = 0; index < replicate_count; ++index) {
        ReplicateResult result = run_replicate(base_seed, config,
                                               static_cast<std::uint64_t>(index), category,
                                               overrides, options);
        raw_percentiles.push_back(result.percentile_l_s_person);
        batch.replicates.push_back(std::move(result));
    }
    batch.summary = statistics::summarize_replicates(raw_percentiles, config.confidence_level);
    return batch;
}

}  // namespace simulation
}  // namespace ecai_cpp2

#endif  // ECAI_CPP2_SIMULATION_H
