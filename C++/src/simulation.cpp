// Stage 10: one replicate, end to end.
// See include/simulation.h for the full contract.

#include "simulation.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "model.h"

namespace ecai_cpp2 {
namespace simulation {

namespace {

void require_finite(double value, const std::string& name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument("simulation: " + name + " must be finite");
    }
}

void require_positive(double value, const std::string& name) {
    require_finite(value, name);
    if (!(value > 0.0)) {
        throw std::invalid_argument("simulation: " + name + " must be strictly positive");
    }
}

void require_in_closed_unit(double value, const std::string& name) {
    require_finite(value, name);
    if (!(value >= 0.0) || !(value <= 1.0)) {
        throw std::invalid_argument("simulation: " + name + " must lie in [0, 1]");
    }
}

void require_matching(double configured, double resolved, const std::string& name,
                      const std::string& context = "resolve_scenario") {
    if (!(configured == resolved)) {
        throw std::invalid_argument(context + ": effective config " + name +
                                     " disagrees with the resolved category/override value");
    }
}

// Rejects an EvaluationOrder outside the declared enum before any output
// allocation or assignment (defense in depth: build_schedule below also
// rejects an unknown order, but only after the output vector and audit rows
// already exist).
void validate_evaluation_order(EvaluationOrder order) {
    switch (order) {
        case EvaluationOrder::forward:
        case EvaluationOrder::reverse:
        case EvaluationOrder::interleaved:
            return;
    }
    throw std::invalid_argument("simulation: unknown EvaluationOrder");
}

// Rejects an OutputMode outside the declared enum before any scenario field
// is trusted (same defensive pattern as validate_evaluation_order above; no
// `default:` label, so a future third enumerator trips -Wswitch).
void validate_output_mode(OutputMode mode) {
    switch (mode) {
        case OutputMode::ecai:
        case OutputMode::probability:
            return;
    }
    throw std::invalid_argument("simulation: unknown OutputMode");
}

// Rejects an InfectorPolicy outside the declared enum before any scenario
// field is trusted (same defensive pattern as validate_output_mode above).
void validate_infector_policy(InfectorPolicy policy) {
    switch (policy) {
        case InfectorPolicy::community:
        case InfectorPolicy::risk_averse:
            return;
    }
    throw std::invalid_argument("simulation: unknown InfectorPolicy");
}

// Rejects a QerFamily outside the declared enum before any scenario field is
// trusted (same defensive pattern as validate_infector_policy above).
void validate_qer_family(qer_direct::QerFamily family) {
    switch (family) {
        case qer_direct::QerFamily::Jones:
        case qer_direct::QerFamily::JonesFitted:
        case qer_direct::QerFamily::Mikszewski:
            return;
    }
    throw std::invalid_argument("simulation: unknown QerFamily");
}

// Cross-checks the four inputs to run_replicate_from_plan before any bank is
// built: config and scenario must name the same canonical category as the
// supplied record, every config scientific/run input must equal the
// resolved scenario value exactly, and the simulation count must agree
// across config, scenario and plan.
void validate_plan_context(const seed_manager::EffectiveConfig& config,
                            const OccupancyCategory& category, const EffectiveScenario& scenario,
                            const sampling_plan::SamplingPlan& plan) {
    static const std::string kContext = "run_replicate_from_plan";
    if (config.category_name != category.canonical_name) {
        throw std::invalid_argument(kContext +
                                     ": effective config names a different category than the "
                                     "supplied record");
    }
    if (scenario.category_name != category.canonical_name) {
        throw std::invalid_argument(
            kContext + ": scenario names a different category than the supplied record");
    }
    if (config.i0 != scenario.occupants) {
        throw std::invalid_argument(kContext + ": effective config I0 disagrees with the scenario");
    }
    require_matching(config.volume_m3, scenario.volume_m3, "VOL", kContext);
    require_matching(config.mask_efficiency, scenario.mask_efficiency, "mask efficiency", kContext);
    require_matching(config.community_infection_rate, scenario.community_infection_rate,
                     "community infection rate", kContext);
    require_matching(config.target_probability, scenario.target_probability, "target probability",
                     kContext);
    require_matching(config.duration_hours, scenario.duration_hours, "duration", kContext);
    require_matching(config.percentile, scenario.percentile, "percentile", kContext);
    if (config.simulations_n != plan.n || scenario.simulations_n != plan.n) {
        throw std::invalid_argument(
            kContext + ": simulations N disagrees among config, scenario and plan");
    }
    // Policy agreement (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md Phase 5): the
    // enumerator values are frozen to agree between simulation::InfectorPolicy
    // and sampling_plan::InfectorPolicy, so a plain static_cast is the
    // authoritative conversion (see sampling_plan.h). A mismatch here is a
    // tampered or stale plan, not a user error.
    if (static_cast<sampling_plan::InfectorPolicy>(scenario.infector_policy) != plan.policy) {
        throw std::invalid_argument(
            kContext + ": scenario infector policy disagrees with the plan's policy");
    }
}

// Every scenario domain check, run before any allocation or sampling.
void validate_scenario(const EffectiveScenario& scenario) {
    if (scenario.occupants < 1) {
        throw std::invalid_argument("simulation: occupants (I0) must be at least 1");
    }
    require_positive(scenario.volume_m3, "volume_m3");
    require_in_closed_unit(scenario.mask_efficiency, "mask_efficiency");
    require_in_closed_unit(scenario.community_infection_rate, "community_infection_rate");
    require_finite(scenario.target_probability, "target_probability");
    if (!(scenario.target_probability > 0.0) || !(scenario.target_probability < 1.0)) {
        throw std::invalid_argument(
            "simulation: target_probability must lie strictly inside (0, 1)");
    }
    require_positive(scenario.duration_hours, "duration_hours");
    require_in_closed_unit(scenario.percentile, "percentile");

    validate_output_mode(scenario.mode);
    if (scenario.mode == OutputMode::probability) {
        require_finite(scenario.input_ecai_l_s_person, "input_ecai_l_s_person");
        if (!(scenario.input_ecai_l_s_person >= 0.0)) {
            throw std::invalid_argument("simulation: input_ecai_l_s_person must be >= 0");
        }
    } else {
        // An ECAi scenario carrying a probability-mode input is a defect in
        // the caller, not a user error: reject rather than silently ignore.
        if (scenario.input_ecai_overridden || scenario.input_ecai_l_s_person != 0.0) {
            throw std::invalid_argument(
                "simulation: an ECAi-mode scenario must not carry a probability-mode input ECAi");
        }
    }

    // Risk-averse infector policy (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md
    // Phase 3): authoritative, applies to both output modes, and runs before
    // any LHS generation or allocation (this function is called from
    // resolve_scenario and from every replicate entry point, all of which
    // run before Stage 6 sampling). I0 < 1 is already rejected unconditionally
    // above; the community_infection_rate <= 0 check below is the one
    // domain restriction that differs between the two policies (community
    // mode accepts a rate of exactly 0).
    validate_infector_policy(scenario.infector_policy);
    if (scenario.infector_policy == InfectorPolicy::risk_averse) {
        if (scenario.occupants < 1) {
            throw std::invalid_argument(
                "simulation: risk-averse infector policy requires occupants (I0) >= 1");
        }
        if (!(scenario.community_infection_rate > 0.0)) {
            throw std::invalid_argument(
                "simulation: risk-averse infector policy requires community_infection_rate > 0");
        }
    }

    // QER family/profile (docs/QER_DISTRIBUTION_IMPLEMENTATION_PLAN.md Phase
    // 5): authoritative, applies to both output modes and both infector
    // policies, and runs before any LHS generation or allocation. For the
    // default Jones family, qer_mikszewski_profile must stay unset and
    // qer_direct_params is not consulted (see EffectiveScenario::qer_family).
    validate_qer_family(scenario.qer_family);
    if (scenario.qer_family == qer_direct::QerFamily::Jones) {
        if (scenario.qer_mikszewski_profile.has_value()) {
            throw std::invalid_argument(
                "simulation: a Jones-family scenario must not carry a resolved Mikszewski profile");
        }
    } else {
        if (!std::isfinite(scenario.qer_direct_params.mu_log10) ||
            !std::isfinite(scenario.qer_direct_params.sigma_log10) ||
            !(scenario.qer_direct_params.sigma_log10 > 0.0)) {
            throw std::invalid_argument(
                "simulation: a direct-QER-family scenario requires a finite mu_log10 and a finite, "
                "positive sigma_log10");
        }
        if (scenario.qer_family == qer_direct::QerFamily::Mikszewski &&
            !scenario.qer_mikszewski_profile.has_value()) {
            throw std::invalid_argument(
                "simulation: a Mikszewski-family scenario requires a resolved profile");
        }
        if (scenario.qer_family == qer_direct::QerFamily::JonesFitted &&
            scenario.qer_mikszewski_profile.has_value()) {
            throw std::invalid_argument(
                "simulation: a Jones-fitted scenario must not carry a resolved Mikszewski profile");
        }
    }
}

// Rejects a malformed plan before any work: the same accounting contract
// Stage 7 enforces, re-checked here so nothing is allocated first.
void validate_plan(const sampling_plan::SamplingPlan& plan) {
    if (plan.counts.size() != plan.n) {
        throw std::invalid_argument("simulation: plan.n disagrees with plan.counts.size()");
    }
    const sampling_plan::CountSummary summary = sampling_plan::summarize_counts(plan.counts);
    if (summary.n_zero != plan.n_zero || summary.n_active != plan.n_active ||
        summary.k_total != plan.k_total || plan.susceptible_draws != summary.n_active ||
        plan.infector_draws != summary.k_total) {
        throw std::invalid_argument("simulation: plan accounting disagrees with its own counts");
    }
}

// Rejects banks that were not sized for this plan before any work.
void validate_banks(const sampling_plan::SamplingPlan& plan,
                    const sampling_plan::ParameterBanks& banks) {
    if (banks.n != plan.n || banks.n_active != plan.n_active || banks.k_total != plan.k_total) {
        throw std::invalid_argument(
            "simulation: parameter banks were not sized for this sampling plan");
    }
    for (const sampling_plan::BankId id : sampling_plan::kAllBankIds) {
        const std::size_t expected =
            sampling_plan::is_infector_bank(id) ? plan.k_total : plan.n_active;
        if (banks.bank(id).values.size() != expected) {
            throw std::invalid_argument("simulation: bank " + sampling_plan::bank_id_name(id) +
                                         " has the wrong planned length");
        }
    }
}

// Rejects family-aware banks that were not sized for this plan, or that name
// a different family than the scenario, before any work (package 3b's
// counterpart of validate_banks above).
void validate_family_banks(const sampling_plan::SamplingPlan& plan,
                            const sampling_plan::FamilyParameterBanks& banks,
                            qer_direct::QerFamily expected_family) {
    if (banks.family != expected_family) {
        throw std::invalid_argument(
            "simulation: parameter banks were built for a different QER family than the scenario "
            "names");
    }
    if (banks.n != plan.n || banks.n_active != plan.n_active || banks.k_total != plan.k_total) {
        throw std::invalid_argument(
            "simulation: family parameter banks were not sized for this sampling plan");
    }
    for (const sampling_plan::ParameterBank& bank : banks.susceptible_banks) {
        if (bank.values.size() != plan.n_active) {
            throw std::invalid_argument(
                "simulation: a susceptible bank has the wrong planned length");
        }
    }
    if (banks.direct_qer_bank.values.size() != plan.k_total) {
        throw std::invalid_argument("simulation: the direct-QER bank has the wrong planned length");
    }
}

// The evaluation schedule over already-assigned active rows. Purely a
// computational order: the mapping it iterates was fixed before it was built.
std::vector<std::size_t> build_schedule(const std::vector<std::size_t>& active_indices,
                                         EvaluationOrder order) {
    switch (order) {
        case EvaluationOrder::forward:
            return active_indices;
        case EvaluationOrder::reverse: {
            std::vector<std::size_t> schedule(active_indices.rbegin(), active_indices.rend());
            return schedule;
        }
        case EvaluationOrder::interleaved: {
            std::vector<std::size_t> schedule;
            schedule.reserve(active_indices.size());
            std::size_t front = 0;
            std::size_t back = active_indices.size();
            while (front < back) {
                schedule.push_back(active_indices[front]);
                ++front;
                if (front < back) {
                    --back;
                    schedule.push_back(active_indices[back]);
                }
            }
            return schedule;
        }
    }
    throw std::invalid_argument("simulation: unknown EvaluationOrder");
}

// One active row's ECAi evaluation (Jones/full-QER family). Factored out of
// run_replicate_from_banks' ECAi loop so the schedule-based path (reverse,
// interleaved, or collect_audit=true) and the fused forward/no-audit path
// below execute the exact same statements in the exact same order for a
// given row -- the two paths differ only in how a row is discovered, never
// in what happens once it is.
void evaluate_ecai_row(std::size_t row, const sampling_plan::SimulationAssignment& simulation,
                        model::ECAiInputs& inputs, const ReplicateOptions& options,
                        ReplicateResult& result) {
    double qer_sum_per_h = 0.0;
    for (const sampling_plan::InfectorParameters& infector : simulation.infectors) {
        const double qer = model::compute_qer(infector).qer_per_h;
        model::require_qer_value(qer, "ECAi QER_i");
        qer_sum_per_h += qer;
        ++result.counters.qer_evaluations;
    }
    if (!std::isfinite(qer_sum_per_h)) {
        throw std::runtime_error("ECAi QER_sum overflowed while summing the supplied QER values");
    }

    inputs.susceptible_pbr_m3h = simulation.susceptible.pbr_m3h;
    inputs.gamma_per_h = simulation.susceptible.gamma_per_h;
    inputs.lambda_per_h = simulation.susceptible.lambda_per_h;

    const model::ECAiResult ecai = model::compute_ecai_from_qer_sum(inputs, qer_sum_per_h);
    ++result.counters.ecai_evaluations;
    if (ecai.clamped) {
        ++result.counters.negative_clamp_count;
    }
    if (!std::isfinite(ecai.ecai_l_s_person)) {
        throw std::runtime_error("simulation: a nonfinite ECAi was produced at simulation " +
                                  std::to_string(row));
    }
    result.ecai_l_s_person[row] = ecai.ecai_l_s_person;
    const double original_target = inputs.target_probability;
    for (std::size_t target = 0; target < options.paired_target_probabilities.size(); ++target) {
        inputs.target_probability = options.paired_target_probabilities[target];
        result.paired_ecai_l_s_person[target][row] =
            model::compute_ecai_from_qer_sum(inputs, qer_sum_per_h).ecai_l_s_person;
    }
    inputs.target_probability = original_target;
    if (options.collect_audit) {
        RowAudit& audit_row = result.audit.rows[row];
        audit_row.susceptible_pbr_m3h = inputs.susceptible_pbr_m3h;
        audit_row.gamma_per_h = inputs.gamma_per_h;
        audit_row.lambda_per_h = inputs.lambda_per_h;
        audit_row.qer_sum_per_h = ecai.qer_sum_per_h;
        audit_row.ecai_raw_l_s_person = ecai.ecai_raw_l_s_person;
        audit_row.ecai_l_s_person = ecai.ecai_l_s_person;
        audit_row.clamped = ecai.clamped;
    }
}

// Package 3b's direct-QER counterpart of evaluate_ecai_row above: identical
// except that each infector's QER comes directly from
// FamilyInfectorAssignment::direct_qer_quanta_per_h rather than from
// model::compute_qer(InfectorParameters).
void evaluate_family_ecai_row(std::size_t row,
                               const sampling_plan::FamilySimulationAssignment& simulation,
                               model::ECAiInputs& inputs, const ReplicateOptions& options,
                               ReplicateResult& result) {
    double qer_sum_per_h = 0.0;
    for (const sampling_plan::FamilyInfectorAssignment& infector : simulation.infectors) {
        const double qer = infector.direct_qer_quanta_per_h;
        model::require_qer_value(qer, "ECAi QER_i");
        qer_sum_per_h += qer;
        ++result.counters.qer_evaluations;
    }
    if (!std::isfinite(qer_sum_per_h)) {
        throw std::runtime_error("ECAi QER_sum overflowed while summing the supplied QER values");
    }

    inputs.susceptible_pbr_m3h = simulation.susceptible.pbr_m3h;
    inputs.gamma_per_h = simulation.susceptible.gamma_per_h;
    inputs.lambda_per_h = simulation.susceptible.lambda_per_h;

    const model::ECAiResult ecai = model::compute_ecai_from_qer_sum(inputs, qer_sum_per_h);
    ++result.counters.ecai_evaluations;
    if (ecai.clamped) {
        ++result.counters.negative_clamp_count;
    }
    if (!std::isfinite(ecai.ecai_l_s_person)) {
        throw std::runtime_error("simulation: a nonfinite ECAi was produced at simulation " +
                                  std::to_string(row));
    }
    result.ecai_l_s_person[row] = ecai.ecai_l_s_person;
    const double original_target = inputs.target_probability;
    for (std::size_t target = 0; target < options.paired_target_probabilities.size(); ++target) {
        inputs.target_probability = options.paired_target_probabilities[target];
        result.paired_ecai_l_s_person[target][row] =
            model::compute_ecai_from_qer_sum(inputs, qer_sum_per_h).ecai_l_s_person;
    }
    inputs.target_probability = original_target;
    if (options.collect_audit) {
        RowAudit& audit_row = result.audit.rows[row];
        audit_row.susceptible_pbr_m3h = inputs.susceptible_pbr_m3h;
        audit_row.gamma_per_h = inputs.gamma_per_h;
        audit_row.lambda_per_h = inputs.lambda_per_h;
        audit_row.qer_sum_per_h = ecai.qer_sum_per_h;
        audit_row.ecai_raw_l_s_person = ecai.ecai_raw_l_s_person;
        audit_row.ecai_l_s_person = ecai.ecai_l_s_person;
        audit_row.clamped = ecai.clamped;
    }
}

// Step 3's probability-mode counterpart of evaluate_ecai_row above: same
// per-row shape (QER accumulation, then one model call), reused by both the
// schedule-based loop and the fused forward/no-audit loop
// (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md section 4.3).
void evaluate_probability_row(std::size_t row, const sampling_plan::SimulationAssignment& simulation,
                               model::ProbabilityInputs& inputs, const ReplicateOptions& options,
                               ReplicateResult& result) {
    double qer_sum_per_h = 0.0;
    for (const sampling_plan::InfectorParameters& infector : simulation.infectors) {
        const double qer = model::compute_qer(infector).qer_per_h;
        model::require_qer_value(qer, "Probability QER_i");
        qer_sum_per_h += qer;
        ++result.counters.qer_evaluations;
    }
    if (!std::isfinite(qer_sum_per_h)) {
        throw std::runtime_error(
            "Probability QER_sum overflowed while summing the supplied QER values");
    }

    inputs.susceptible_pbr_m3h = simulation.susceptible.pbr_m3h;
    inputs.gamma_per_h = simulation.susceptible.gamma_per_h;
    inputs.lambda_per_h = simulation.susceptible.lambda_per_h;

    const model::ProbabilityResult probability =
        model::compute_probability_from_qer_sum(inputs, qer_sum_per_h);
    ++result.counters.probability_evaluations;
    if (probability.saturated) {
        ++result.counters.saturated_count;
    }
    if (!std::isfinite(probability.probability)) {
        throw std::runtime_error(
            "simulation: a nonfinite probability was produced at simulation " +
            std::to_string(row));
    }
    result.probability[row] = probability.probability;
    if (options.collect_audit) {
        RowAudit& audit_row = result.audit.rows[row];
        audit_row.qer_sum_per_h = probability.qer_sum_per_h;
        audit_row.probability = probability.probability;
        audit_row.saturated = probability.saturated;
    }
}

void evaluate_probability_schedule(const std::vector<std::size_t>& schedule,
                                    const sampling_plan::ReplicateAssignments& assignments,
                                    const EffectiveScenario& scenario,
                                    const ReplicateOptions& options, ReplicateResult& result) {
    model::ProbabilityInputs inputs;
    inputs.volume_m3 = scenario.volume_m3;
    inputs.occupants = scenario.occupants;
    inputs.mask_efficiency = scenario.mask_efficiency;
    inputs.duration_h = scenario.duration_hours;
    inputs.input_ecai_l_s_person = scenario.input_ecai_l_s_person;
    for (const std::size_t row : schedule) {
        evaluate_probability_row(row, assignments.simulations[row], inputs, options, result);
    }
}

// Package 3b's direct-QER counterpart of evaluate_probability_row above: the
// row loop and QER accumulation are identical; the only difference is that
// each infector's QER comes directly from
// FamilyInfectorAssignment::direct_qer_quanta_per_h -- already computed once
// when the direct-QER bank was generated -- rather than from
// model::compute_qer(InfectorParameters).
void evaluate_family_probability_row(std::size_t row,
                                      const sampling_plan::FamilySimulationAssignment& simulation,
                                      model::ProbabilityInputs& inputs, const ReplicateOptions& options,
                                      ReplicateResult& result) {
    double qer_sum_per_h = 0.0;
    for (const sampling_plan::FamilyInfectorAssignment& infector : simulation.infectors) {
        const double qer = infector.direct_qer_quanta_per_h;
        model::require_qer_value(qer, "Probability QER_i");
        qer_sum_per_h += qer;
        ++result.counters.qer_evaluations;
    }
    if (!std::isfinite(qer_sum_per_h)) {
        throw std::runtime_error(
            "Probability QER_sum overflowed while summing the supplied QER values");
    }

    inputs.susceptible_pbr_m3h = simulation.susceptible.pbr_m3h;
    inputs.gamma_per_h = simulation.susceptible.gamma_per_h;
    inputs.lambda_per_h = simulation.susceptible.lambda_per_h;

    const model::ProbabilityResult probability =
        model::compute_probability_from_qer_sum(inputs, qer_sum_per_h);
    ++result.counters.probability_evaluations;
    if (probability.saturated) {
        ++result.counters.saturated_count;
    }
    if (!std::isfinite(probability.probability)) {
        throw std::runtime_error(
            "simulation: a nonfinite probability was produced at simulation " +
            std::to_string(row));
    }
    result.probability[row] = probability.probability;
    if (options.collect_audit) {
        RowAudit& audit_row = result.audit.rows[row];
        audit_row.qer_sum_per_h = probability.qer_sum_per_h;
        audit_row.probability = probability.probability;
        audit_row.saturated = probability.saturated;
    }
}

void evaluate_family_probability_schedule(const std::vector<std::size_t>& schedule,
                                           const sampling_plan::FamilyReplicateAssignments& assignments,
                                           const EffectiveScenario& scenario,
                                           const ReplicateOptions& options, ReplicateResult& result) {
    model::ProbabilityInputs inputs;
    inputs.volume_m3 = scenario.volume_m3;
    inputs.occupants = scenario.occupants;
    inputs.mask_efficiency = scenario.mask_efficiency;
    inputs.duration_h = scenario.duration_hours;
    inputs.input_ecai_l_s_person = scenario.input_ecai_l_s_person;
    for (const std::size_t row : schedule) {
        evaluate_family_probability_row(row, assignments.simulations[row], inputs, options, result);
    }
}

}  // namespace

double percentile_linear(const std::vector<double>& values, double p) {
    if (values.empty()) {
        throw std::invalid_argument("percentile_linear: the value list must not be empty");
    }
    for (const double value : values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("percentile_linear: every value must be finite");
        }
    }
    if (!std::isfinite(p)) {
        throw std::invalid_argument("percentile_linear: p must be finite");
    }
    if (!(p >= 0.0) || !(p <= 1.0)) {
        throw std::invalid_argument("percentile_linear: p must lie in [0, 1]");
    }

    // A working copy: the caller's index-ordered vector is never reordered.
    std::vector<double> working(values);

    // Frozen convention (D8 / section 12): rank = p * (n - 1), interpolated
    // linearly between the order statistics at floor(rank) and ceil(rank).
    const double rank = p * static_cast<double>(working.size() - 1);
    const double floor_rank = std::floor(rank);
    const double ceil_rank = std::ceil(rank);
    if (!(floor_rank >= 0.0) || !(ceil_rank <= static_cast<double>(working.size() - 1))) {
        throw std::runtime_error("percentile_linear: interpolation rank fell outside the list");
    }
    const std::size_t lo = static_cast<std::size_t>(floor_rank);
    const std::size_t hi = static_cast<std::size_t>(ceil_rank);
    const double frac = rank - floor_rank;

    // Only the two order statistics at lo and hi (hi is lo or lo + 1) are
    // needed, so select them in O(n) instead of sorting the whole vector.
    // nth_element guarantees working[lo] equals the lo-th order statistic and
    // partitions the remainder so its minimum equals the (lo + 1)-th.
    const auto lo_it = working.begin() + static_cast<std::ptrdiff_t>(lo);
    std::nth_element(working.begin(), lo_it, working.end());
    const double lo_value = *lo_it;
    const double hi_value =
        (hi == lo) ? lo_value : *std::min_element(lo_it + 1, working.end());
    return lo_value * (1.0 - frac) + hi_value * frac;
}

EffectiveScenario resolve_scenario(const OccupancyCategory& category,
                                    const seed_manager::EffectiveConfig& config,
                                    const EffectiveOverrides& overrides) {
    if (config.category_name != category.canonical_name) {
        throw std::invalid_argument("resolve_scenario: effective config names category \"" +
                                     config.category_name + "\" but the supplied record is \"" +
                                     category.canonical_name + "\"");
    }

    EffectiveScenario scenario;
    scenario.category_name = category.canonical_name;
    scenario.occupants = overrides.occupants ? *overrides.occupants : category.i0;
    scenario.volume_m3 = overrides.volume_m3 ? *overrides.volume_m3 : category.volume_m3;
    scenario.mask_efficiency =
        overrides.mask_efficiency ? *overrides.mask_efficiency : category.mask_efficiency;
    scenario.community_infection_rate = overrides.community_infection_rate
                                            ? *overrides.community_infection_rate
                                            : category.community_infection_rate;
    scenario.target_probability = config.target_probability;
    scenario.duration_hours = config.duration_hours;
    scenario.percentile = config.percentile;

    // Branch from the resolved scenario mode (overrides.mode), never from a
    // raw override value alone: this is the single source of truth every
    // downstream consumer (simulation, and eventually output) dispatches on.
    validate_output_mode(overrides.mode);
    scenario.mode = overrides.mode;
    if (scenario.mode == OutputMode::probability) {
        scenario.input_ecai_l_s_person = overrides.input_ecai_l_s_person
                                              ? *overrides.input_ecai_l_s_person
                                              : category.published_ecai_Lsp;
        scenario.input_ecai_overridden = overrides.input_ecai_l_s_person.has_value();
        if (!overrides.input_ecai_l_s_person) {
            // Defensive: the published Table 4 value is never expected to be
            // out of domain, but resolution still checks it before it is
            // trusted as the probability-mode default.
            if (!std::isfinite(category.published_ecai_Lsp) || category.published_ecai_Lsp < 0.0) {
                throw std::invalid_argument(
                    "resolve_scenario: category published_ecai_Lsp is not a finite, "
                    "non-negative value");
            }
        }
    }
    // In ECAi mode, input_ecai_l_s_person and input_ecai_overridden stay at
    // their EffectiveScenario defaults (0.0 / false) regardless of
    // overrides.input_ecai_l_s_person: the CLI has already rejected supplying
    // --input-ecai without --mode probability.

    // Resolved once here, and forwarded unchanged to every downstream
    // consumer; nothing branches on a raw CLI flag (Phase 3).
    scenario.infector_policy = overrides.infector_policy;
    scenario.infector_policy_explicit = overrides.infector_policy_explicit;

    // Resolved once here, and forwarded unchanged to every downstream
    // consumer, mirroring infector_policy above (docs/QER_DISTRIBUTION_
    // IMPLEMENTATION_PLAN.md Phase 5). qer_direct_params/qer_mikszewski_
    // profile are copied through as supplied -- resolve_scenario performs no
    // qer_direct table lookup itself; that already happened once, in the CLI
    // layer, through the work package 1 API.
    scenario.qer_family = overrides.qer_family;
    scenario.qer_family_explicit = overrides.qer_family_explicit;
    scenario.qer_mikszewski_profile = overrides.qer_mikszewski_profile;
    scenario.qer_direct_params = overrides.qer_direct_params;

    scenario.occupants_overridden = overrides.occupants.has_value();
    scenario.volume_overridden = overrides.volume_m3.has_value();
    scenario.mask_efficiency_overridden = overrides.mask_efficiency.has_value();
    scenario.community_infection_rate_overridden = overrides.community_infection_rate.has_value();
    scenario.custom_scenario = scenario.occupants_overridden || scenario.volume_overridden ||
                                scenario.mask_efficiency_overridden ||
                                scenario.community_infection_rate_overridden;

    // The effective configuration is what the seeds are fingerprinted from,
    // so it must agree with the resolved values exactly.
    if (config.i0 != scenario.occupants) {
        throw std::invalid_argument(
            "resolve_scenario: effective config I0 disagrees with the resolved "
            "category/override value");
    }
    require_matching(config.volume_m3, scenario.volume_m3, "VOL");
    require_matching(config.mask_efficiency, scenario.mask_efficiency, "mask efficiency");
    require_matching(config.community_infection_rate, scenario.community_infection_rate,
                     "community infection rate");

    if (config.simulations_n == 0) {
        throw std::invalid_argument("resolve_scenario: simulations N must be at least 1");
    }
    if (config.simulations_n > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::invalid_argument("resolve_scenario: simulations N exceeds size_t");
    }
    scenario.simulations_n = static_cast<std::size_t>(config.simulations_n);

    validate_scenario(scenario);
    return scenario;
}

ReplicateResult run_replicate_from_banks(const sampling_plan::SamplingPlan& plan,
                                          const sampling_plan::ParameterBanks& banks,
                                          const EffectiveScenario& scenario,
                                          const ReplicateOptions& options) {
    // --- Every rejection happens before the output vector is allocated. ---
    validate_scenario(scenario);
    validate_plan(plan);
    // Cross-checks scenario.infector_policy against plan.policy and, for a
    // risk-averse plan, re-asserts the no-zero-infector invariant
    // (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md Phase 5: "tampered
    // scenario/plan/output policy disagreement" must be rejected here, before
    // any output vector is allocated).
    sampling_plan::validate_plan_policy_context(
        plan, static_cast<sampling_plan::InfectorPolicy>(scenario.infector_policy),
        options.policy);
    validate_banks(plan, banks);
    if (scenario.simulations_n != plan.n) {
        throw std::invalid_argument(
            "simulation: the scenario's simulation count disagrees with the plan's N");
    }
    validate_evaluation_order(options.evaluation_order);

    ReplicateResult result;
    result.mode = scenario.mode;
    result.infector_policy = scenario.infector_policy;
    result.percentile = scenario.percentile;
    result.counters.n = plan.n;
    result.infector_counts = plan.counts;
    for (double target : options.paired_target_probabilities) {
        if (!(target > 0.0 && target < 1.0) || !std::isfinite(target) ||
            scenario.mode != OutputMode::ecai) {
            throw std::invalid_argument("simulation: invalid paired ECAi target");
        }
        result.paired_ecai_l_s_person.emplace_back(plan.n, 0.0);
    }

    // Step 1: every zero-infector row becomes exactly +0.0 before any
    // assignment or model input is requested (MODEL_SPECIFICATION.md 11.1,
    // D3). The whole output vector is initialized to +0.0, so an active row
    // that is later clamped also lands on +0.0. `outputs` names whichever of
    // the two mode-specific vectors this replicate populates; the other
    // stays empty.
    std::vector<double>& outputs =
        (scenario.mode == OutputMode::probability) ? result.probability : result.ecai_l_s_person;
    outputs.assign(plan.n, 0.0);
    for (std::size_t i = 0; i < plan.counts.size(); ++i) {
        if (plan.counts[i] == 0) {
            ++result.counters.zero_infector_count;
        }
    }

    // Step 2: materialize the full simulation/person -> value mapping in
    // stable simulation index order, before any equation runs.
    const sampling_plan::ReplicateAssignments assignments =
        sampling_plan::assign_replicate(banks, plan);
    if (assignments.simulations.size() != plan.n) {
        throw std::runtime_error("simulation: the assignment did not cover every simulation");
    }
    result.counters.active_count = assignments.n_active;
    result.counters.k_total = assignments.k_total;

    // Step 2b: record identity (consumption ids, per-row people) from the
    // materialized mapping, in stable simulation index order. Nothing here
    // depends on the evaluation schedule chosen below.
    result.audit_collected = options.collect_audit;
    if (options.collect_audit) {
        result.audit.rows.resize(plan.n);
    }

    // Step 3 (forward, collect_audit=false only): the schedule is exactly
    // the active rows in ascending simulation index -- the same order this
    // pass already discovers them in -- and there is no audit.evaluation_
    // order to populate, so each active row is evaluated as soon as it is
    // found, without ever materializing an active_indices vector or the
    // schedule copy build_schedule(active_indices, forward) would otherwise
    // return. Each active row's evaluation (counters, per-row model call,
    // output write) is identical to how the schedule-based path in the else
    // branch evaluates that same row; only the interleaving with the
    // simulation-index validation differs between the two branches, not any
    // individual row's result.
    if (!options.collect_audit && options.evaluation_order == EvaluationOrder::forward) {
        if (scenario.mode == OutputMode::ecai) {
            model::ECAiInputs inputs;
            inputs.volume_m3 = scenario.volume_m3;
            inputs.occupants = scenario.occupants;
            inputs.mask_efficiency = scenario.mask_efficiency;
            inputs.duration_h = scenario.duration_hours;
            inputs.target_probability = scenario.target_probability;
            for (std::size_t row = 0; row < assignments.simulations.size(); ++row) {
                const sampling_plan::SimulationAssignment& simulation = assignments.simulations[row];
                if (simulation.simulation_index != row) {
                    throw std::runtime_error("simulation: assignments are not in simulation index order");
                }
                if (!simulation.active) {
                    continue;
                }
                result.counters.susceptible_values_consumed += simulation.susceptible.ids.size();
                result.counters.infector_values_consumed +=
                    simulation.infectors.size() * sampling_plan::InfectorParameters{}.ids.size();
                evaluate_ecai_row(row, simulation, inputs, options, result);
            }
        } else {
            model::ProbabilityInputs inputs;
            inputs.volume_m3 = scenario.volume_m3;
            inputs.occupants = scenario.occupants;
            inputs.mask_efficiency = scenario.mask_efficiency;
            inputs.duration_h = scenario.duration_hours;
            inputs.input_ecai_l_s_person = scenario.input_ecai_l_s_person;
            for (std::size_t row = 0; row < assignments.simulations.size(); ++row) {
                const sampling_plan::SimulationAssignment& simulation = assignments.simulations[row];
                if (simulation.simulation_index != row) {
                    throw std::runtime_error("simulation: assignments are not in simulation index order");
                }
                if (!simulation.active) {
                    continue;
                }
                result.counters.susceptible_values_consumed += simulation.susceptible.ids.size();
                result.counters.infector_values_consumed +=
                    simulation.infectors.size() * sampling_plan::InfectorParameters{}.ids.size();
                evaluate_probability_row(row, simulation, inputs, options, result);
            }
        }
    } else {
        std::vector<std::size_t> active_indices;
        active_indices.reserve(assignments.n_active);
        for (std::size_t row = 0; row < assignments.simulations.size(); ++row) {
            const sampling_plan::SimulationAssignment& simulation = assignments.simulations[row];
            if (simulation.simulation_index != row) {
                throw std::runtime_error("simulation: assignments are not in simulation index order");
            }
            if (options.collect_audit) {
                RowAudit& audit_row = result.audit.rows[row];
                audit_row.simulation_index = row;
                audit_row.infector_count = simulation.infector_count;
                audit_row.active = simulation.active;
            }
            if (!simulation.active) {
                continue;
            }
            active_indices.push_back(row);
            result.counters.susceptible_values_consumed += simulation.susceptible.ids.size();
            result.counters.infector_values_consumed +=
                simulation.infectors.size() * sampling_plan::InfectorParameters{}.ids.size();
            if (options.collect_audit) {
                RowAudit& audit_row = result.audit.rows[row];
                audit_row.susceptible_ids = simulation.susceptible.ids;
                audit_row.infector_ids.reserve(simulation.infectors.size());
                for (const sampling_plan::InfectorParameters& infector : simulation.infectors) {
                    audit_row.infector_ids.push_back(infector.ids);
                }
                for (const std::uint64_t id : simulation.susceptible.ids) {
                    result.audit.consumption_ids.push_back(id);
                }
                for (const sampling_plan::InfectorParameters& infector : simulation.infectors) {
                    for (const std::uint64_t id : infector.ids) {
                        result.audit.consumption_ids.push_back(id);
                    }
                }
            }
        }

        // Step 3: evaluate the active rows on the requested schedule,
        // writing each result back to its own stable simulation index.
        // Infectors within a row are always summed in list order, so no
        // schedule can change a row's QER sum.
        const std::vector<std::size_t> schedule =
            build_schedule(active_indices, options.evaluation_order);
        if (options.collect_audit) {
            result.audit.evaluation_order = schedule;
        }
        if (scenario.mode == OutputMode::ecai) {
            model::ECAiInputs inputs;
            inputs.volume_m3 = scenario.volume_m3;
            inputs.occupants = scenario.occupants;
            inputs.mask_efficiency = scenario.mask_efficiency;
            inputs.duration_h = scenario.duration_hours;
            inputs.target_probability = scenario.target_probability;
            for (const std::size_t row : schedule) {
                evaluate_ecai_row(row, assignments.simulations[row], inputs, options, result);
            }
        } else {
            evaluate_probability_schedule(schedule, assignments, scenario, options, result);
        }
    }

    for (const double value : outputs) {
        if (value == 0.0) {
            ++result.counters.output_zero_count;
        }
    }

    // Cursor reconciliation: assign_replicate already called finalize(), so
    // this records the proof rather than establishing it.
    result.cursors_exhausted = true;
    for (const sampling_plan::BankId id : sampling_plan::kAllBankIds) {
        const std::size_t ordinal = static_cast<std::size_t>(id);
        result.planned_bank_lengths[ordinal] = banks.bank(id).values.size();
        result.final_cursors[ordinal] = assignments.final_cursor(id);
        result.bank_generated[ordinal] = banks.bank(id).generated;
        if (result.final_cursors[ordinal] != result.planned_bank_lengths[ordinal]) {
            result.cursors_exhausted = false;
        }
    }

    // Step 4: the frozen percentile, from a sorted copy.
    if (scenario.mode == OutputMode::probability) {
        result.percentile_probability = percentile_linear(outputs, scenario.percentile);
    } else {
        result.percentile_l_s_person = percentile_linear(outputs, scenario.percentile);
    }
    return result;
}

ReplicateResult run_replicate_from_family_banks(const sampling_plan::SamplingPlan& plan,
                                                 const sampling_plan::FamilyParameterBanks& banks,
                                                 const EffectiveScenario& scenario,
                                                 const ReplicateOptions& options) {
    // --- Every rejection happens before the output vector is allocated. ---
    validate_scenario(scenario);
    validate_plan(plan);
    sampling_plan::validate_plan_policy_context(
        plan, static_cast<sampling_plan::InfectorPolicy>(scenario.infector_policy),
        options.policy);
    validate_family_banks(plan, banks, scenario.qer_family);
    if (scenario.simulations_n != plan.n) {
        throw std::invalid_argument(
            "simulation: the scenario's simulation count disagrees with the plan's N");
    }
    validate_evaluation_order(options.evaluation_order);

    ReplicateResult result;
    result.mode = scenario.mode;
    result.infector_policy = scenario.infector_policy;
    result.percentile = scenario.percentile;
    result.counters.n = plan.n;
    result.infector_counts = plan.counts;
    for (double target : options.paired_target_probabilities) {
        if (!(target > 0.0 && target < 1.0) || !std::isfinite(target) ||
            scenario.mode != OutputMode::ecai) {
            throw std::invalid_argument("simulation: invalid paired ECAi target");
        }
        result.paired_ecai_l_s_person.emplace_back(plan.n, 0.0);
    }

    // Step 1: every zero-infector row becomes exactly +0.0 before any
    // assignment or model input is requested (MODEL_SPECIFICATION.md 11.1,
    // D3), identically to run_replicate_from_banks.
    std::vector<double>& outputs =
        (scenario.mode == OutputMode::probability) ? result.probability : result.ecai_l_s_person;
    outputs.assign(plan.n, 0.0);
    for (std::size_t i = 0; i < plan.counts.size(); ++i) {
        if (plan.counts[i] == 0) {
            ++result.counters.zero_infector_count;
        }
    }

    // Step 2: materialize the full simulation/person -> value mapping in
    // stable simulation index order, before any equation runs. Every
    // infector here is direct-sourced: banks.family != Jones, so
    // assign_family_replicate never constructs an InfectorParameters or
    // calls model::compute_qer.
    const sampling_plan::FamilyReplicateAssignments assignments =
        sampling_plan::assign_family_replicate(banks, plan);
    if (assignments.simulations.size() != plan.n) {
        throw std::runtime_error("simulation: the assignment did not cover every simulation");
    }
    result.counters.active_count = assignments.n_active;
    result.counters.k_total = assignments.k_total;

    result.audit_collected = options.collect_audit;
    if (options.collect_audit) {
        result.audit.rows.resize(plan.n);
    }

    // Step 3 (forward, collect_audit=false only): mirrors the fused forward
    // path in run_replicate_from_banks. Each active row's evaluation
    // (counters, per-row model call, output write) is identical to how the
    // scheduled path in the else branch below evaluates that same row; only
    // the interleaving with the simulation-index validation differs between
    // the two branches, not any individual row's result.
    if (!options.collect_audit && options.evaluation_order == EvaluationOrder::forward) {
        if (scenario.mode == OutputMode::ecai) {
            model::ECAiInputs inputs;
            inputs.volume_m3 = scenario.volume_m3;
            inputs.occupants = scenario.occupants;
            inputs.mask_efficiency = scenario.mask_efficiency;
            inputs.duration_h = scenario.duration_hours;
            inputs.target_probability = scenario.target_probability;
            for (std::size_t row = 0; row < assignments.simulations.size(); ++row) {
                const sampling_plan::FamilySimulationAssignment& simulation = assignments.simulations[row];
                if (simulation.simulation_index != row) {
                    throw std::runtime_error("simulation: assignments are not in simulation index order");
                }
                if (!simulation.active) {
                    continue;
                }
                result.counters.susceptible_values_consumed += simulation.susceptible.ids.size();
                result.counters.infector_values_consumed += simulation.infectors.size();
                evaluate_family_ecai_row(row, simulation, inputs, options, result);
            }
        } else {
            model::ProbabilityInputs inputs;
            inputs.volume_m3 = scenario.volume_m3;
            inputs.occupants = scenario.occupants;
            inputs.mask_efficiency = scenario.mask_efficiency;
            inputs.duration_h = scenario.duration_hours;
            inputs.input_ecai_l_s_person = scenario.input_ecai_l_s_person;
            for (std::size_t row = 0; row < assignments.simulations.size(); ++row) {
                const sampling_plan::FamilySimulationAssignment& simulation = assignments.simulations[row];
                if (simulation.simulation_index != row) {
                    throw std::runtime_error("simulation: assignments are not in simulation index order");
                }
                if (!simulation.active) {
                    continue;
                }
                result.counters.susceptible_values_consumed += simulation.susceptible.ids.size();
                result.counters.infector_values_consumed += simulation.infectors.size();
                evaluate_family_probability_row(row, simulation, inputs, options, result);
            }
        }
    } else {
        std::vector<std::size_t> active_indices;
        active_indices.reserve(assignments.n_active);
        for (std::size_t row = 0; row < assignments.simulations.size(); ++row) {
            const sampling_plan::FamilySimulationAssignment& simulation = assignments.simulations[row];
            if (simulation.simulation_index != row) {
                throw std::runtime_error("simulation: assignments are not in simulation index order");
            }
            if (options.collect_audit) {
                RowAudit& audit_row = result.audit.rows[row];
                audit_row.simulation_index = row;
                audit_row.infector_count = simulation.infector_count;
                audit_row.active = simulation.active;
            }
            if (!simulation.active) {
                continue;
            }
            active_indices.push_back(row);
            result.counters.susceptible_values_consumed += simulation.susceptible.ids.size();
            // One direct-QER value per infector, not eight: the direct-family
            // persistent-bank formula is N + 3*N_active + K_total, not
            // N + 3*N_active + 8*K_total (docs/QER_DISTRIBUTION_IMPLEMENTATION_
            // PLAN.md Phase 3).
            result.counters.infector_values_consumed += simulation.infectors.size();
            if (options.collect_audit) {
                RowAudit& audit_row = result.audit.rows[row];
                audit_row.susceptible_ids = simulation.susceptible.ids;
                for (const std::uint64_t id : simulation.susceptible.ids) {
                    result.audit.consumption_ids.push_back(id);
                }
                for (const sampling_plan::FamilyInfectorAssignment& infector : simulation.infectors) {
                    result.audit.consumption_ids.push_back(infector.id);
                }
            }
        }

        // Step 3: evaluate the active rows on the requested schedule, writing
        // each result back to its own stable simulation index.
        const std::vector<std::size_t> schedule =
            build_schedule(active_indices, options.evaluation_order);
        if (options.collect_audit) {
            result.audit.evaluation_order = schedule;
        }
        if (scenario.mode == OutputMode::ecai) {
            model::ECAiInputs inputs;
            inputs.volume_m3 = scenario.volume_m3;
            inputs.occupants = scenario.occupants;
            inputs.mask_efficiency = scenario.mask_efficiency;
            inputs.duration_h = scenario.duration_hours;
            inputs.target_probability = scenario.target_probability;
            for (const std::size_t row : schedule) {
                evaluate_family_ecai_row(row, assignments.simulations[row], inputs, options, result);
            }
        } else {
            evaluate_family_probability_schedule(schedule, assignments, scenario, options, result);
        }
    }

    for (const double value : outputs) {
        if (value == 0.0) {
            ++result.counters.output_zero_count;
        }
    }

    // Cursor reconciliation: assign_family_replicate already called
    // finalize(), so this records the proof rather than establishing it. The
    // 11-entry Jones-BankId-indexed arrays below do not describe a direct-
    // family bank layout (there is no per-BankId direct-QER bank), so they
    // are deliberately left at their default (all-zero/false) values rather
    // than force-fit a non-Jones layout into them.
    result.cursors_exhausted = true;

    // Step 4: the frozen percentile, from a sorted copy.
    if (scenario.mode == OutputMode::probability) {
        result.percentile_probability = percentile_linear(outputs, scenario.percentile);
    } else {
        result.percentile_l_s_person = percentile_linear(outputs, scenario.percentile);
    }
    return result;
}

ReplicateResult run_replicate_from_plan(const seed_manager::SeedManager& seeds,
                                         const seed_manager::EffectiveConfig& config,
                                         std::uint64_t replicate_index,
                                         const sampling_plan::SamplingPlan& plan,
                                         const OccupancyCategory& category,
                                         const EffectiveScenario& scenario,
                                         const ReplicateOptions& options) {
    validate_scenario(scenario);
    validate_plan(plan);
    // Reject a policy mismatch (or a risk-averse plan carrying a zero-infector
    // row) before build_parameter_banks allocates or generates any LHS list
    // for this plan (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md Phase 5).
    // validate_plan_context below repeats the plain policy-equality half of
    // this check for its own aggregate diagnostic; run_replicate_from_banks
    // repeats the full check again once banks exist, so a caller that skips
    // straight to it is still protected.
    sampling_plan::validate_plan_policy_context(
        plan, static_cast<sampling_plan::InfectorPolicy>(scenario.infector_policy),
        options.policy);
    validate_plan_context(config, category, scenario, plan);
    // Branch once, after plan construction (docs/QER_DISTRIBUTION_
    // IMPLEMENTATION_PLAN.md Phase 5): full Jones keeps its existing,
    // byte-identical build_parameter_banks/run_replicate_from_banks path.
    // Jones-fitted and Mikszewski dispatch through the family-aware bank
    // builder and run_replicate_from_family_banks instead, and never
    // construct an InfectorParameters or call the full-Jones compute_qer().
    // The banks live only for this call: they are destroyed on return, so no
    // caller can re-read, re-extend, or re-traverse them.
    if (scenario.qer_family == qer_direct::QerFamily::Jones) {
        const sampling_plan::ParameterBanks banks = sampling_plan::build_parameter_banks(
            seeds, config, replicate_index, plan, category, options.policy);
        ReplicateResult result = run_replicate_from_banks(plan, banks, scenario, options);
        result.replicate_index = replicate_index;
        return result;
    }
    const sampling_plan::FamilyParameterBanks banks = sampling_plan::build_family_parameter_banks(
        seeds, config, replicate_index, plan, category, scenario.qer_family,
        scenario.qer_direct_params, options.policy);
    ReplicateResult result = run_replicate_from_family_banks(plan, banks, scenario, options);
    result.replicate_index = replicate_index;
    return result;
}

ReplicateResult run_replicate(std::uint64_t base_seed,
                               const seed_manager::EffectiveConfig& config,
                               std::uint64_t replicate_index, const OccupancyCategory& category,
                               const EffectiveOverrides& overrides,
                               const ReplicateOptions& options) {
    const EffectiveScenario scenario = resolve_scenario(category, config, overrides);
    const seed_manager::SeedManager seeds(base_seed);
    // Forwards the resolved policy unchanged (docs/RISK_AVERSE_
    // IMPLEMENTATION_PLAN.md Phase 5): build_sampling_plan's own seed
    // derivation is identical regardless of `infector_policy` (Phase 4), so
    // this branch changes only how each u is transformed into a count, never
    // which u values are generated or consumed.
    const sampling_plan::SamplingPlan plan = sampling_plan::build_sampling_plan(
        seeds, config, replicate_index, scenario.simulations_n, scenario.occupants,
        scenario.community_infection_rate, options.policy,
        static_cast<sampling_plan::InfectorPolicy>(scenario.infector_policy));
    return run_replicate_from_plan(seeds, config, replicate_index, plan, category, scenario,
                                    options);
}

}  // namespace simulation
}  // namespace ecai_cpp2
