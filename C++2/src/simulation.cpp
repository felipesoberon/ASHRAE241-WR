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

    // A sorted copy: the caller's index-ordered vector is never reordered.
    std::vector<double> sorted(values);
    std::sort(sorted.begin(), sorted.end());

    // Frozen convention (D8 / section 12): rank = p * (n - 1), interpolated
    // linearly between the order statistics at floor(rank) and ceil(rank).
    const double rank = p * static_cast<double>(sorted.size() - 1);
    const double floor_rank = std::floor(rank);
    const double ceil_rank = std::ceil(rank);
    if (!(floor_rank >= 0.0) || !(ceil_rank <= static_cast<double>(sorted.size() - 1))) {
        throw std::runtime_error("percentile_linear: interpolation rank fell outside the list");
    }
    const std::size_t lo = static_cast<std::size_t>(floor_rank);
    const std::size_t hi = static_cast<std::size_t>(ceil_rank);
    const double frac = rank - floor_rank;
    return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
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
    validate_banks(plan, banks);
    if (scenario.simulations_n != plan.n) {
        throw std::invalid_argument(
            "simulation: the scenario's simulation count disagrees with the plan's N");
    }
    validate_evaluation_order(options.evaluation_order);

    ReplicateResult result;
    result.percentile = scenario.percentile;
    result.counters.n = plan.n;
    result.infector_counts = plan.counts;

    // Step 1: every zero-infector row becomes exactly +0.0 before any
    // assignment or model input is requested (MODEL_SPECIFICATION.md 11.1,
    // D3). The whole output vector is initialized to +0.0, so an active row
    // that is later clamped also lands on +0.0.
    result.ecai_l_s_person.assign(plan.n, 0.0);
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

    // Step 3: evaluate the active rows on the requested schedule, writing
    // each result back to its own stable simulation index. Infectors within
    // a row are always summed in list order, so no schedule can change a
    // row's QER sum.
    const std::vector<std::size_t> schedule =
        build_schedule(active_indices, options.evaluation_order);
    if (options.collect_audit) {
        result.audit.evaluation_order = schedule;
    }
    // One scratch QER list for the whole schedule: it is emptied and then
    // filled in list order for each row, so its contents when compute_ecai
    // runs are exactly that row's per-infector QERs in the same order as
    // before. The room-level inputs are the same for every row and are set
    // once; only the susceptible's three values change per row.
    std::vector<double> qer_per_h;
    model::ECAiInputs inputs;
    inputs.volume_m3 = scenario.volume_m3;
    inputs.occupants = scenario.occupants;
    inputs.mask_efficiency = scenario.mask_efficiency;
    inputs.duration_h = scenario.duration_hours;
    inputs.target_probability = scenario.target_probability;
    for (const std::size_t row : schedule) {
        const sampling_plan::SimulationAssignment& simulation = assignments.simulations[row];
        qer_per_h.clear();
        qer_per_h.reserve(simulation.infectors.size());
        for (const sampling_plan::InfectorParameters& infector : simulation.infectors) {
            qer_per_h.push_back(model::compute_qer(infector).qer_per_h);
            ++result.counters.qer_evaluations;
        }

        inputs.susceptible_pbr_m3h = simulation.susceptible.pbr_m3h;
        inputs.gamma_per_h = simulation.susceptible.gamma_per_h;
        inputs.lambda_per_h = simulation.susceptible.lambda_per_h;

        const model::ECAiResult ecai = model::compute_ecai(inputs, qer_per_h);
        ++result.counters.ecai_evaluations;
        if (ecai.clamped) {
            ++result.counters.negative_clamp_count;
        }
        if (!std::isfinite(ecai.ecai_l_s_person)) {
            throw std::runtime_error("simulation: a nonfinite ECAi was produced at simulation " +
                                      std::to_string(row));
        }
        result.ecai_l_s_person[row] = ecai.ecai_l_s_person;
        if (options.collect_audit) {
            RowAudit& audit_row = result.audit.rows[row];
            audit_row.qer_sum_per_h = ecai.qer_sum_per_h;
            audit_row.ecai_raw_l_s_person = ecai.ecai_raw_l_s_person;
            audit_row.ecai_l_s_person = ecai.ecai_l_s_person;
            audit_row.clamped = ecai.clamped;
        }
    }

    for (const double value : result.ecai_l_s_person) {
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
    result.percentile_l_s_person = percentile_linear(result.ecai_l_s_person, scenario.percentile);
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
    validate_plan_context(config, category, scenario, plan);
    // The banks live only for this call: they are destroyed on return, so no
    // caller can re-read, re-extend, or re-traverse them.
    const sampling_plan::ParameterBanks banks = sampling_plan::build_parameter_banks(
        seeds, config, replicate_index, plan, category, options.policy);
    ReplicateResult result = run_replicate_from_banks(plan, banks, scenario, options);
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
    const sampling_plan::SamplingPlan plan =
        sampling_plan::build_sampling_plan(seeds, config, replicate_index, scenario.simulations_n,
                                           scenario.occupants, scenario.community_infection_rate,
                                           options.policy);
    return run_replicate_from_plan(seeds, config, replicate_index, plan, category, scenario,
                                    options);
}

}  // namespace simulation
}  // namespace ecai_cpp2
