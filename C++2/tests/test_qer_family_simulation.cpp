// Package 3b tests: direct-QER family integration into the Stage 10
// simulation ECAi and Probability loops, per
// docs/QER_DISTRIBUTION_IMPLEMENTATION_PLAN.md Phase 5.
//
// Covers: all three families (Jones, JonesFitted, Mikszewski) in both output
// modes; community mixed/all-zero infector rows; risk-averse (N_active == N,
// no zero rows); one direct draw vs eight Jones draws per infector; direct
// cursor exhaustion/mismatch rejection; finite/nonnegative ECAi and [0, 1]
// probability; and that a direct-family replicate never touches the
// full-Jones infector banks or model::compute_qer (proven both by the
// full-Jones banks being left ungenerated -- so any accidental use would
// throw -- and by bit-exact reproduction of the ECAi/probability outputs
// directly from the precomputed direct-QER bank values).

#include "simulation.h"

#include "model.h"
#include "qer_direct.h"
#include "sampling_plan.h"
#include "seed_manager.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

void check(bool condition, const std::string& description) {
    if (!condition) {
        std::cerr << "FAIL: " << description << "\n";
        ++g_failures;
    } else {
        std::cout << "PASS: " << description << "\n";
    }
}

template <typename Fn>
bool throws_runtime_error(Fn&& fn) {
    try {
        fn();
    } catch (const std::runtime_error&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

template <typename Fn>
bool throws_invalid_argument(Fn&& fn) {
    try {
        fn();
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

bool is_positive_zero(double v) { return v == 0.0 && !std::signbit(v); }

using ecai_cpp2::OccupancyCategory;
using ecai_cpp2::qer_direct::DirectQerParams;
using ecai_cpp2::qer_direct::QerFamily;
using ecai_cpp2::sampling_plan::build_family_parameter_banks;
using ecai_cpp2::sampling_plan::build_parameter_banks;
using ecai_cpp2::sampling_plan::CountSummary;
using ecai_cpp2::sampling_plan::FamilyParameterBanks;
using ecai_cpp2::sampling_plan::InfectorPolicy;
using ecai_cpp2::sampling_plan::ParameterBanks;
using ecai_cpp2::sampling_plan::SamplingPlan;
using ecai_cpp2::sampling_plan::summarize_counts;
using ecai_cpp2::seed_manager::EffectiveConfig;
using ecai_cpp2::seed_manager::SeedManager;
using ecai_cpp2::simulation::EffectiveOverrides;
using ecai_cpp2::simulation::EffectiveScenario;
using ecai_cpp2::simulation::OutputMode;
using ecai_cpp2::simulation::ReplicateOptions;
using ecai_cpp2::simulation::ReplicateResult;

constexpr double kPi = 3.14159265358979323846;

// Delegates to the production equation helpers (model.cpp
// compute_ecai_from_qer_sum / compute_probability_from_qer_sum) instead of a
// hand-expanded restatement, so this file's oracle is bit-identical to the
// production operation order by construction. The independent checks in this
// file are that the direct bank values were consumed in infector order and
// summed correctly (qer_sum below), not the equation arithmetic itself.
double expected_ecai(double pbr_s, double gamma, double lambda, double volume, long long i0,
                     double mu, double duration_h, double target_p, double qer_sum) {
    ecai_cpp2::model::ECAiInputs inputs;
    inputs.susceptible_pbr_m3h = pbr_s;
    inputs.gamma_per_h = gamma;
    inputs.lambda_per_h = lambda;
    inputs.volume_m3 = volume;
    inputs.occupants = i0;
    inputs.mask_efficiency = mu;
    inputs.duration_h = duration_h;
    inputs.target_probability = target_p;
    return ecai_cpp2::model::compute_ecai_from_qer_sum(inputs, qer_sum).ecai_l_s_person;
}

double expected_probability(double pbr_s, double gamma, double lambda, double volume, long long i0,
                            double mu, double duration_h, double input_ecai, double qer_sum) {
    ecai_cpp2::model::ProbabilityInputs inputs;
    inputs.susceptible_pbr_m3h = pbr_s;
    inputs.gamma_per_h = gamma;
    inputs.lambda_per_h = lambda;
    inputs.volume_m3 = volume;
    inputs.occupants = i0;
    inputs.mask_efficiency = mu;
    inputs.duration_h = duration_h;
    inputs.input_ecai_l_s_person = input_ecai;
    return ecai_cpp2::model::compute_probability_from_qer_sum(inputs, qer_sum).probability;
}

EffectiveConfig make_config(double community_rate = 0.5) {
    EffectiveConfig config;
    config.category_name = "Patient room";
    config.i0 = 3;
    config.volume_m3 = 81.0;
    config.mask_efficiency = 0.3;
    config.community_infection_rate = community_rate;
    config.target_probability = 0.001;
    config.duration_hours = 1.0;
    config.percentile = 0.95;
    config.confidence_level = 0.95;
    config.simulations_n = 20;
    config.replicates_r = 5;
    return config;
}

OccupancyCategory make_test_category() {
    OccupancyCategory category{};
    category.canonical_name = "Patient room";
    category.group = "Health Care";
    category.is_health_care_group = true;
    category.i0 = 3;
    category.volume_m3 = 81.0;
    category.probability_susceptible_present = 0.057;
    category.people_outdoor_airflow_Lsp = 15.0;
    category.area_outdoor_airflow_Lsm2 = 0.0;
    category.baseline_outdoor_airflow_Lsp = 15.0;
    category.community_infection_rate = 0.5;
    category.mask_efficiency = 0.30;
    category.pbr_gm_m3h = 0.62;
    category.pbr_gsd = 1.3;
    category.activity_breathing_pct = 100.0;
    category.activity_talking_pct = 0.0;
    category.activity_vocalizing_pct = 0.0;
    category.cdrop_gm_e4_per_m3 = 30.0;
    category.cdrop_gsd = 1.1;
    category.d_gm_e6_m = 2.0;
    category.d_gsd = 1.1;
    category.published_ecai_Lsp = 35.0;
    category.published_rr = 0.53;
    return category;
}

SamplingPlan make_plan(const std::vector<long long>& counts,
                        InfectorPolicy policy = InfectorPolicy::community) {
    const CountSummary summary = summarize_counts(counts);
    SamplingPlan plan;
    plan.counts = counts;
    plan.n = counts.size();
    plan.n_zero = summary.n_zero;
    plan.n_active = summary.n_active;
    plan.k_total = summary.k_total;
    plan.susceptible_draws = summary.n_active;
    plan.infector_draws = summary.k_total;
    plan.policy = policy;
    return plan;
}

const DirectQerParams kDirectParams{-1.0, 1.4};  // arbitrary, finite, sigma > 0

EffectiveScenario make_scenario(const SamplingPlan& plan, const OccupancyCategory& category,
                                 QerFamily family, OutputMode mode,
                                 InfectorPolicy policy = InfectorPolicy::community) {
    EffectiveScenario scenario;
    scenario.category_name = category.canonical_name;
    scenario.occupants = category.i0;
    scenario.volume_m3 = category.volume_m3;
    scenario.mask_efficiency = category.mask_efficiency;
    scenario.community_infection_rate = category.community_infection_rate;
    scenario.target_probability = 0.001;
    scenario.duration_hours = 1.0;
    scenario.percentile = 0.95;
    scenario.simulations_n = plan.n;
    scenario.mode = mode;
    if (mode == OutputMode::probability) {
        scenario.input_ecai_l_s_person = category.published_ecai_Lsp;
    }
    scenario.infector_policy = static_cast<ecai_cpp2::simulation::InfectorPolicy>(policy);
    scenario.qer_family = family;
    scenario.qer_family_explicit = (family != QerFamily::Jones);
    if (family == QerFamily::Mikszewski) {
        scenario.qer_mikszewski_profile = ecai_cpp2::qer_direct::MikszewskiProfile{
            ecai_cpp2::qer_direct::MikszewskiOrganism::SarsCov2,
            ecai_cpp2::qer_direct::MikszewskiActivity::StandingSpeaking};
    }
    if (family != QerFamily::Jones) {
        scenario.qer_direct_params = kDirectParams;
    }
    return scenario;
}

}  // namespace

int main() {
    using ecai_cpp2::simulation::run_replicate_from_banks;
    using ecai_cpp2::simulation::run_replicate_from_family_banks;

    const SeedManager seeds(42);
    const EffectiveConfig config = make_config();
    const OccupancyCategory category = make_test_category();

    // ================================================================
    // Bit-exact reproduction: a JonesFitted ECAi replicate's outputs are
    // exactly expected_ecai() evaluated at the precomputed direct-QER bank
    // value(s) for that row -- never a Jones-mechanistic recomputation.
    // ================================================================
    {
        const SamplingPlan plan = make_plan({1, 2, 0, 1});  // n_active=3, k_total=4
        const FamilyParameterBanks banks = build_family_parameter_banks(
            seeds, config, 0, plan, category, QerFamily::JonesFitted, kDirectParams);
        const EffectiveScenario scenario =
            make_scenario(plan, category, QerFamily::JonesFitted, OutputMode::ecai);
        const ReplicateResult result =
            run_replicate_from_family_banks(plan, banks, scenario, ReplicateOptions{});

        check(result.ecai_l_s_person.size() == plan.n, "ecai output has N entries");
        check(is_positive_zero(result.ecai_l_s_person[2]), "zero-infector row is exactly +0.0");

        // Reconstruct susceptible/direct-QER values independently from the
        // banks (already produced by build_family_parameter_banks, whose own
        // correctness is package 3a's contract) and re-derive the expected
        // ECAi via the frozen equation restated in this file.
        std::size_t susceptible_index = 0;
        std::size_t direct_index = 0;
        for (std::size_t row = 0; row < plan.counts.size(); ++row) {
            if (plan.counts[row] == 0) continue;
            const double pbr_s = banks.susceptible_banks[0].values[susceptible_index];
            const double gamma = banks.susceptible_banks[1].values[susceptible_index];
            const double lambda = banks.susceptible_banks[2].values[susceptible_index];
            ++susceptible_index;
            double qer_sum = 0.0;
            for (long long k = 0; k < plan.counts[row]; ++k) {
                qer_sum += banks.direct_qer_bank.values[direct_index];
                ++direct_index;
            }
            const double expected = expected_ecai(pbr_s, gamma, lambda, category.volume_m3,
                                                   category.i0, category.mask_efficiency, 1.0,
                                                   0.001, qer_sum);
            check(std::memcmp(&expected, &result.ecai_l_s_person[row], sizeof(double)) == 0,
                  "row " + std::to_string(row) + " ECAi matches the direct-QER-sourced equation exactly");
        }
        check(direct_index == plan.k_total, "every direct-QER value was consumed exactly once");
        check(result.cursors_exhausted, "direct-family cursors report exhausted");
        check(result.counters.qer_evaluations == plan.k_total,
              "qer_evaluations counts one per realized infector");
    }

    // ================================================================
    // Same bit-exact reproduction in probability mode.
    // ================================================================
    {
        const SamplingPlan plan = make_plan({2, 1});  // n_active=2, k_total=3
        const FamilyParameterBanks banks = build_family_parameter_banks(
            seeds, config, 1, plan, category, QerFamily::Mikszewski, kDirectParams);
        const EffectiveScenario scenario =
            make_scenario(plan, category, QerFamily::Mikszewski, OutputMode::probability);
        const ReplicateResult result =
            run_replicate_from_family_banks(plan, banks, scenario, ReplicateOptions{});

        check(result.probability.size() == plan.n, "probability output has N entries");
        std::size_t susceptible_index = 0;
        std::size_t direct_index = 0;
        for (std::size_t row = 0; row < plan.counts.size(); ++row) {
            const double pbr_s = banks.susceptible_banks[0].values[susceptible_index];
            const double gamma = banks.susceptible_banks[1].values[susceptible_index];
            const double lambda = banks.susceptible_banks[2].values[susceptible_index];
            ++susceptible_index;
            double qer_sum = 0.0;
            for (long long k = 0; k < plan.counts[row]; ++k) {
                qer_sum += banks.direct_qer_bank.values[direct_index];
                ++direct_index;
            }
            const double expected =
                expected_probability(pbr_s, gamma, lambda, category.volume_m3, category.i0,
                                      category.mask_efficiency, 1.0, category.published_ecai_Lsp,
                                      qer_sum);
            check(std::memcmp(&expected, &result.probability[row], sizeof(double)) == 0,
                  "row " + std::to_string(row) + " probability matches the direct-QER equation exactly");
            check(result.probability[row] >= 0.0 && result.probability[row] <= 1.0,
                  "row " + std::to_string(row) + " probability lies in [0, 1]");
        }
    }

    // ================================================================
    // All three families, both modes: finite/in-domain outputs, and one
    // direct draw versus eight Jones draws per infector (counters).
    // ================================================================
    for (const QerFamily family : {QerFamily::Jones, QerFamily::JonesFitted, QerFamily::Mikszewski}) {
        for (const OutputMode mode : {OutputMode::ecai, OutputMode::probability}) {
            const SamplingPlan plan = make_plan({3, 0, 2, 1, 0});  // n_active=3, k_total=6
            const EffectiveScenario scenario = make_scenario(plan, category, family, mode);
            ReplicateResult result;
            if (family == QerFamily::Jones) {
                const ParameterBanks banks =
                    build_parameter_banks(seeds, config, 2, plan, category);
                result = run_replicate_from_banks(plan, banks, scenario, ReplicateOptions{});
                check(result.counters.infector_values_consumed == 8 * plan.k_total,
                      "Jones family consumes 8 values per infector (family=" +
                          ecai_cpp2::qer_direct::to_string(family) + ")");
            } else {
                const FamilyParameterBanks banks = build_family_parameter_banks(
                    seeds, config, 2, plan, category, family, kDirectParams);
                result = run_replicate_from_family_banks(plan, banks, scenario, ReplicateOptions{});
                check(result.counters.infector_values_consumed == plan.k_total,
                      "direct family consumes exactly 1 value per infector (family=" +
                          ecai_cpp2::qer_direct::to_string(family) + ")");
            }
            check(result.counters.qer_evaluations == plan.k_total,
                  "qer_evaluations == K_total regardless of family");
            const std::vector<double>& outputs =
                (mode == OutputMode::probability) ? result.probability : result.ecai_l_s_person;
            for (std::size_t row = 0; row < outputs.size(); ++row) {
                check(std::isfinite(outputs[row]), "output is finite");
                if (mode == OutputMode::ecai) {
                    check(outputs[row] >= 0.0, "ECAi output is nonnegative");
                } else {
                    check(outputs[row] >= 0.0 && outputs[row] <= 1.0, "probability output lies in [0, 1]");
                }
            }
            check(is_positive_zero(outputs[1]) && is_positive_zero(outputs[4]),
                  "zero-infector rows remain exactly +0.0");
        }
    }

    // ================================================================
    // Community, all-zero: nothing is drawn, every output is exactly +0.0.
    // ================================================================
    {
        const SamplingPlan plan = make_plan({0, 0, 0, 0});
        const FamilyParameterBanks banks = build_family_parameter_banks(
            seeds, config, 3, plan, category, QerFamily::JonesFitted, kDirectParams);
        check(!banks.direct_qer_bank.generated, "an all-zero plan generates no direct-QER bank");
        const EffectiveScenario scenario =
            make_scenario(plan, category, QerFamily::JonesFitted, OutputMode::ecai);
        const ReplicateResult result =
            run_replicate_from_family_banks(plan, banks, scenario, ReplicateOptions{});
        for (const double value : result.ecai_l_s_person) {
            check(is_positive_zero(value), "all-zero community plan output is exactly +0.0");
        }
        check(result.counters.zero_infector_count == plan.n, "every row is counted zero-infector");
        check(result.counters.qer_evaluations == 0, "no QER evaluation for an all-zero plan");
    }

    // ================================================================
    // Risk-averse: N_active == N, no zero rows, direct family.
    // ================================================================
    {
        const SamplingPlan plan = make_plan({2, 1, 3, 1}, InfectorPolicy::risk_averse);
        const FamilyParameterBanks banks = build_family_parameter_banks(
            seeds, config, 4, plan, category, QerFamily::Mikszewski, kDirectParams);
        const EffectiveScenario scenario = make_scenario(
            plan, category, QerFamily::Mikszewski, OutputMode::ecai, InfectorPolicy::risk_averse);
        const ReplicateResult result =
            run_replicate_from_family_banks(plan, banks, scenario, ReplicateOptions{});
        check(result.counters.active_count == plan.n, "risk-averse: every row is active");
        check(result.counters.zero_infector_count == 0, "risk-averse: no zero-infector rows");
        check(result.counters.k_total == plan.k_total, "risk-averse: K_total matches the plan");
        for (const double value : result.ecai_l_s_person) {
            check(std::isfinite(value) && value >= 0.0, "risk-averse ECAi output is finite and nonnegative");
        }
    }

    // ================================================================
    // Direct cursor exhaustion / mismatch: banks built for a different plan
    // must be rejected rather than silently under/over-consumed.
    // ================================================================
    {
        const SamplingPlan plan_small = make_plan({1, 1});
        const SamplingPlan plan_big = make_plan({1, 1, 1});
        const FamilyParameterBanks banks_for_small = build_family_parameter_banks(
            seeds, config, 5, plan_small, category, QerFamily::JonesFitted, kDirectParams);
        const EffectiveScenario scenario =
            make_scenario(plan_big, category, QerFamily::JonesFitted, OutputMode::ecai);
        check(throws_invalid_argument([&]() {
                  run_replicate_from_family_banks(plan_big, banks_for_small, scenario,
                                                   ReplicateOptions{});
              }),
              "banks sized for a smaller plan are rejected before traversal");
    }

    // ================================================================
    // Reproducibility: two traversals of the same banks/plan/scenario, via
    // fresh cursors, agree bit-for-bit.
    // ================================================================
    {
        const SamplingPlan plan = make_plan({1, 2, 1});
        const FamilyParameterBanks banks = build_family_parameter_banks(
            seeds, config, 6, plan, category, QerFamily::Mikszewski, kDirectParams);
        const EffectiveScenario scenario =
            make_scenario(plan, category, QerFamily::Mikszewski, OutputMode::ecai);
        const ReplicateResult first =
            run_replicate_from_family_banks(plan, banks, scenario, ReplicateOptions{});
        const ReplicateResult second =
            run_replicate_from_family_banks(plan, banks, scenario, ReplicateOptions{});
        check(first.ecai_l_s_person.size() == second.ecai_l_s_person.size() &&
                  std::memcmp(first.ecai_l_s_person.data(), second.ecai_l_s_person.data(),
                              first.ecai_l_s_person.size() * sizeof(double)) == 0,
              "two traversals of the same direct-family banks agree bit-for-bit");
    }

    // ================================================================
    // Family/scenario mismatch is rejected before traversal.
    // ================================================================
    {
        const SamplingPlan plan = make_plan({1, 1});
        const FamilyParameterBanks banks = build_family_parameter_banks(
            seeds, config, 7, plan, category, QerFamily::JonesFitted, kDirectParams);
        const EffectiveScenario mismatched =
            make_scenario(plan, category, QerFamily::Mikszewski, OutputMode::ecai);
        check(throws_invalid_argument([&]() {
                  run_replicate_from_family_banks(plan, banks, mismatched, ReplicateOptions{});
              }),
              "a scenario naming a different QER family than the banks is rejected");
    }

    // ================================================================
    // Top-level run_replicate integration: EffectiveOverrides selects the
    // family, resolve_scenario/run_replicate dispatch correctly end to end.
    // ================================================================
    {
        ecai_cpp2::seed_manager::EffectiveConfig top_config;
        top_config.category_name = category.canonical_name;
        top_config.i0 = category.i0;
        top_config.volume_m3 = category.volume_m3;
        top_config.mask_efficiency = category.mask_efficiency;
        top_config.community_infection_rate = category.community_infection_rate;
        top_config.target_probability = 0.001;
        top_config.duration_hours = 1.0;
        top_config.percentile = 0.95;
        top_config.confidence_level = 0.95;
        top_config.simulations_n = 12;
        top_config.replicates_r = 3;

        for (const QerFamily family : {QerFamily::Jones, QerFamily::JonesFitted, QerFamily::Mikszewski}) {
            EffectiveOverrides overrides;
            overrides.qer_family = family;
            overrides.qer_family_explicit = (family != QerFamily::Jones);
            if (family == QerFamily::Mikszewski) {
                overrides.qer_mikszewski_profile = ecai_cpp2::qer_direct::MikszewskiProfile{
                    ecai_cpp2::qer_direct::MikszewskiOrganism::Influenza,
                    ecai_cpp2::qer_direct::MikszewskiActivity::Resting};
            }
            if (family != QerFamily::Jones) {
                overrides.qer_direct_params = kDirectParams;
            }
            const ReplicateResult result =
                ecai_cpp2::simulation::run_replicate(99, top_config, 0, category, overrides);
            check(result.ecai_l_s_person.size() == top_config.simulations_n,
                  "run_replicate produces N outputs for family=" +
                      ecai_cpp2::qer_direct::to_string(family));
            for (const double value : result.ecai_l_s_person) {
                check(std::isfinite(value) && value >= 0.0,
                      "run_replicate output is finite and nonnegative for family=" +
                          ecai_cpp2::qer_direct::to_string(family));
            }
        }
    }

    // ================================================================
    // An out-of-band QerFamily on the scenario is rejected before any
    // allocation (mirrors the existing InfectorPolicy/OutputMode pattern).
    // ================================================================
    {
        const SamplingPlan plan = make_plan({1, 1});
        EffectiveScenario scenario = make_scenario(plan, category, QerFamily::Jones, OutputMode::ecai);
        scenario.qer_family = static_cast<QerFamily>(99);
        const ParameterBanks banks = build_parameter_banks(seeds, config, 8, plan, category);
        check(throws_invalid_argument([&]() {
                  run_replicate_from_banks(plan, banks, scenario, ReplicateOptions{});
              }),
              "an out-of-band QerFamily on the scenario is rejected");
    }

    std::cout << (g_failures == 0 ? "All qer_family_simulation checks passed.\n"
                                   : std::to_string(g_failures) + " qer_family_simulation check(s) failed.\n");
    return g_failures == 0 ? 0 : 1;
}
