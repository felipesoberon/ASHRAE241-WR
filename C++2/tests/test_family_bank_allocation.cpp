// Package 3a tests: family-aware Stage 7 bank allocation and accounting.
//
// Covers docs/QER_DISTRIBUTION_IMPLEMENTATION_PLAN.md Phase 3's family-aware
// sampling-plan requirements: Jones-fitted/Mikszewski allocate only the 3
// susceptible banks plus one direct-QER bank of length K_total (no unused
// full-Jones infector banks/InfectorParameters); zero-infector rows consume
// nothing; mixed-infector-count realizations consume exactly K_total direct
// values in order; risk-averse N_active == N and K_total >= N; the two
// family-specific persistent-bank formulas checked once before allocation;
// cursor/consumption reconciliation; rejection of invalid/mixed-family
// requests before allocation; same-seed reproducibility; and full
// preservation of the existing Jones sampling_plan behavior.

#include "sampling_plan.h"

#include "category_data.h"
#include "qer_direct.h"
#include "seed_manager.h"

#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
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

using ecai_cpp2::OccupancyCategory;
using ecai_cpp2::qer_direct::DirectQerParams;
using ecai_cpp2::qer_direct::QerFamily;
using ecai_cpp2::sampling_plan::build_sampling_plan;
using ecai_cpp2::sampling_plan::CountSummary;
using ecai_cpp2::sampling_plan::FamilyParameterBanks;
using ecai_cpp2::sampling_plan::FamilyParameterCursor;
using ecai_cpp2::sampling_plan::FamilyReplicateAssignments;
using ecai_cpp2::sampling_plan::InfectorPolicy;
using ecai_cpp2::sampling_plan::ResourcePolicy;
using ecai_cpp2::sampling_plan::SamplingPlan;
using ecai_cpp2::sampling_plan::summarize_counts;
using ecai_cpp2::sampling_plan::validate_plan_policy_context;
using ecai_cpp2::seed_manager::EffectiveConfig;
using ecai_cpp2::seed_manager::SeedManager;

EffectiveConfig make_config() {
    EffectiveConfig config;
    config.category_name = "Patient room";
    config.i0 = 3;
    config.volume_m3 = 81.0;
    config.mask_efficiency = 0.3;
    config.community_infection_rate = 0.03;
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
    category.community_infection_rate = 0.03;
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

const DirectQerParams kDirectParams{-1.0, 1.4};

}  // namespace

int main() {
    using namespace ecai_cpp2::sampling_plan;

    const SeedManager seeds(42);
    const EffectiveConfig config = make_config();
    const OccupancyCategory category = make_test_category();

    // ================================================================
    // Bank shapes for each family.
    // ================================================================
    {
        const SamplingPlan plan = make_plan({2, 0, 1, 3});  // n_active=3, k_total=6
        for (QerFamily family : {QerFamily::Jones, QerFamily::JonesFitted, QerFamily::Mikszewski}) {
            const FamilyParameterBanks banks =
                build_family_parameter_banks(seeds, config, 0, plan, category, family, kDirectParams);
            check(banks.family == family, "banks.family matches requested family");
            for (const auto& bank : banks.susceptible_banks) {
                check(bank.values.size() == plan.n_active,
                      "susceptible bank has length N_active for family " +
                          ecai_cpp2::qer_direct::to_string(family));
            }
            if (family == QerFamily::Jones) {
                for (const auto& bank : banks.infector_banks) {
                    check(bank.values.size() == plan.k_total,
                          "Jones infector bank has length K_total");
                }
                check(!banks.direct_qer_bank.generated,
                      "Jones plan leaves the direct-QER bank ungenerated");
                check(banks.direct_qer_bank.values.empty(),
                      "Jones plan leaves the direct-QER bank empty");
            } else {
                check(banks.direct_qer_bank.generated, "direct family generates the direct-QER bank");
                check(banks.direct_qer_bank.values.size() == plan.k_total,
                      "direct-QER bank length equals K_total exactly");
                for (const auto& bank : banks.infector_banks) {
                    check(!bank.generated,
                          "direct family leaves every full-Jones infector bank ungenerated (" +
                              ecai_cpp2::qer_direct::to_string(family) + ")");
                    check(bank.values.empty(),
                          "direct family leaves every full-Jones infector bank empty (" +
                              ecai_cpp2::qer_direct::to_string(family) + ")");
                }
            }
        }
    }

    // ================================================================
    // Direct bank length equals K_total exactly, across several K_total.
    // ================================================================
    {
        for (const auto& counts : std::vector<std::vector<long long>>{
                 {0, 0, 0}, {1}, {3, 0, 2, 0, 1}, {0, 5, 0, 0}}) {
            const SamplingPlan plan = make_plan(counts);
            const FamilyParameterBanks banks = build_family_parameter_banks(
                seeds, config, 0, plan, category, QerFamily::JonesFitted, kDirectParams);
            check(banks.direct_qer_bank.values.size() == plan.k_total,
                  "direct bank length == K_total=" + std::to_string(plan.k_total));
            check(banks.direct_qer_bank.generated == (plan.k_total > 0),
                  "direct bank generated flag matches K_total > 0");
        }
    }

    // ================================================================
    // Zero-infector-count rows consume zero direct-bank entries; mixed
    // infector counts consume exactly K_total total, in order.
    // ================================================================
    {
        const std::vector<long long> counts = {0, 2, 0, 1, 3, 0, 0, 4};
        const SamplingPlan plan = make_plan(counts);
        const FamilyParameterBanks banks = build_family_parameter_banks(
            seeds, config, 0, plan, category, QerFamily::Mikszewski, kDirectParams);
        const FamilyReplicateAssignments assignments = assign_family_replicate(banks, plan);

        std::size_t consumed = 0;
        for (const auto& row : assignments.simulations) {
            if (row.infector_count == 0) {
                check(!row.active && row.infectors.empty(),
                      "zero-infector row consumes zero direct-bank entries");
                continue;
            }
            check(row.active && row.infectors.size() == static_cast<std::size_t>(row.infector_count),
                  "active row consumes exactly infector_count direct values");
            for (std::size_t k = 0; k < row.infectors.size(); ++k) {
                check(row.infectors[k].is_direct, "every direct-family infector is a direct value");
                check(row.infectors[k].bank_index == consumed,
                      "direct-bank entries are consumed in ascending order");
                check(row.infectors[k].direct_qer_quanta_per_h ==
                          banks.direct_qer_bank.values[consumed],
                      "consumed direct value matches the bank at its index");
                ++consumed;
            }
        }
        check(consumed == plan.k_total, "total direct-bank entries consumed == K_total exactly");
        check(assignments.final_direct_cursor == plan.k_total,
              "final direct cursor lands exactly on K_total");
        check(assignments.final_susceptible_cursor == plan.n_active,
              "final susceptible cursor lands exactly on N_active");
        check(assignments.final_infector_cursor == 0,
              "direct family never advances the (unused) Jones infector cursor");
    }

    // ================================================================
    // No unused Jones infector banks / InfectorParameters materialization
    // for direct families: every FamilyInfectorAssignment is a direct value,
    // never a Jones InfectorParameters, and the infector_banks array stays
    // default (checked above); here we also confirm the assignment-level
    // absence of any Jones-shaped payload.
    // ================================================================
    {
        const SamplingPlan plan = make_plan({2, 1});
        const FamilyParameterBanks banks = build_family_parameter_banks(
            seeds, config, 0, plan, category, QerFamily::JonesFitted, kDirectParams);
        const FamilyReplicateAssignments assignments = assign_family_replicate(banks, plan);
        for (const auto& row : assignments.simulations) {
            for (const auto& infector : row.infectors) {
                check(infector.is_direct, "direct-family assignment never carries a Jones payload");
                check(infector.jones.pbr_m3h == 0.0 && infector.jones.c_drop == 0.0,
                      "unused Jones InfectorParameters field stays default-valued");
            }
        }
    }

    // ================================================================
    // Full-Jones family still produces InfectorParameters (not direct
    // values) and matches build_parameter_banks' own layout exactly.
    // ================================================================
    {
        const SamplingPlan plan = make_plan({1, 0, 2});
        const FamilyParameterBanks family_banks = build_family_parameter_banks(
            seeds, config, 0, plan, category, QerFamily::Jones, kDirectParams);
        const ParameterBanks plain_banks = build_parameter_banks(seeds, config, 0, plan, category);

        check(family_banks.susceptible_banks[0].values == plain_banks.bank(BankId::susceptible_pbr).values &&
                  family_banks.susceptible_banks[1].values == plain_banks.bank(BankId::gamma).values &&
                  family_banks.susceptible_banks[2].values == plain_banks.bank(BankId::lambda).values,
              "Jones-family susceptible banks match build_parameter_banks exactly");

        static constexpr std::array<BankId, 8> kInfectorOrder = {
            BankId::infector_pbr, BankId::c_drop, BankId::d_bar, BankId::E,
            BankId::GVL,          BankId::VF,     BankId::RTD,   BankId::DK};
        bool infector_banks_match = true;
        for (std::size_t i = 0; i < kInfectorOrder.size(); ++i) {
            if (family_banks.infector_banks[i].values != plain_banks.bank(kInfectorOrder[i]).values) {
                infector_banks_match = false;
            }
        }
        check(infector_banks_match, "Jones-family infector banks match build_parameter_banks exactly");

        const FamilyReplicateAssignments assignments = assign_family_replicate(family_banks, plan);
        for (const auto& row : assignments.simulations) {
            for (const auto& infector : row.infectors) {
                check(!infector.is_direct, "Jones family produces Jones InfectorParameters, not direct values");
            }
        }
    }

    // ================================================================
    // Risk-averse policy: N_active == N, K_total >= N enforced.
    // ================================================================
    {
        const SamplingPlan plan =
            build_sampling_plan(seeds, config, 0, 6, 3, 0.4, ResourcePolicy{}, InfectorPolicy::risk_averse);
        check(plan.n_active == plan.n, "risk-averse plan has N_active == N");
        check(plan.k_total >= plan.n, "risk-averse plan has K_total >= N");
        // validate_plan_policy_context is the authoritative risk-averse
        // invariant check (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md); it must
        // accept this plan for both families' downstream allocation.
        validate_plan_policy_context(plan, InfectorPolicy::risk_averse);
        const FamilyParameterBanks banks = build_family_parameter_banks(
            seeds, config, 0, plan, category, QerFamily::Mikszewski, kDirectParams);
        check(banks.direct_qer_bank.values.size() == plan.k_total,
              "risk-averse direct bank length == K_total");
        check(banks.susceptible_banks[0].values.size() == plan.n,
              "risk-averse N_active == N, so susceptible banks are length N");

        // A tampered plan claiming N_active < N under risk-averse policy is
        // rejected by validate_plan_policy_context before any allocation.
        SamplingPlan tampered = plan;
        tampered.n_zero = 1;
        tampered.n_active = plan.n - 1;
        check(throws_invalid_argument([&] { validate_plan_policy_context(tampered, InfectorPolicy::risk_averse); }),
              "a risk-averse plan with N_active != N is rejected");
    }

    // ================================================================
    // Persistent memory/resource formulas, checked once before allocation.
    // ================================================================
    {
        const std::size_t n = 10, n_active = 7, k_total = 20;
        const BankBudget jones_budget =
            compute_family_bank_budget(QerFamily::Jones, n, n_active, k_total, ResourcePolicy{});
        check(jones_budget.total_elements == n + 3 * n_active + 8 * k_total,
              "Jones family budget == N + 3*N_active + 8*K_total");

        const BankBudget direct_budget =
            compute_family_bank_budget(QerFamily::JonesFitted, n, n_active, k_total, ResourcePolicy{});
        check(direct_budget.total_elements == n + 3 * n_active + k_total,
              "direct family budget == N + 3*N_active + K_total");
        check(direct_budget.total_elements < jones_budget.total_elements,
              "direct family requires strictly less persistent capacity than Jones for the same plan");

        // The transient LHS/transform/assignment allowance is included in a
        // single ceiling check: a K_total large enough that the persistent
        // formula fits max_bytes but the persistent-plus-transient total does
        // not must still be rejected.
        ResourcePolicy tight;
        tight.max_bytes = sizeof(double) * (n + 3 * n_active + k_total);  // no room for the transient buffer
        check(throws_runtime_error([&] {
                  compute_family_bank_budget(QerFamily::JonesFitted, n, n_active, k_total, tight);
              }),
              "direct-family budget check rejects a plan whose persistent formula fits but whose "
              "transient allowance overflows the byte ceiling");
    }

    // ================================================================
    // Overflow/limit behavior.
    // ================================================================
    {
        ResourcePolicy tiny;
        tiny.max_total_infectors = 5;
        check(throws_runtime_error([&] {
                  compute_family_bank_budget(QerFamily::JonesFitted, 10, 3, 6, tiny);
              }),
              "K_total exceeding max_total_infectors is rejected for a direct family");

        const SamplingPlan plan = make_plan({1, 2, 3});
        ResourcePolicy tiny_bytes;
        tiny_bytes.max_bytes = 4;
        check(throws_runtime_error([&] {
                  build_family_parameter_banks(seeds, config, 0, plan, category, QerFamily::JonesFitted,
                                                kDirectParams, tiny_bytes);
              }),
              "build_family_parameter_banks rejects a byte ceiling too small to allocate");

        check(throws_invalid_argument([&] {
                  compute_family_bank_budget(QerFamily::Jones, 3, 5, 5, ResourcePolicy{});
              }),
              "N_active exceeding N is rejected for the Jones formula");
        check(throws_invalid_argument([&] {
                  compute_family_bank_budget(QerFamily::JonesFitted, 3, 5, 5, ResourcePolicy{});
              }),
              "N_active exceeding N is rejected for the direct formula");
    }

    // ================================================================
    // Cursor/bank-ID tracking, consumption counters, and final
    // reconciliation: total consumed == allocated for every relevant cursor.
    // ================================================================
    {
        const SamplingPlan plan = make_plan({2, 0, 3, 1, 0});
        const FamilyParameterBanks banks = build_family_parameter_banks(
            seeds, config, 0, plan, category, QerFamily::Mikszewski, kDirectParams);
        FamilyParameterCursor cursor(banks);
        check(!cursor.exhausted(), "a fresh cursor over a nonempty plan is not exhausted");
        const FamilyReplicateAssignments assignments = assign_family_replicate(cursor, plan);
        check(cursor.exhausted(), "cursor is exhausted after a full traversal");
        check(assignments.final_susceptible_cursor == plan.n_active &&
                  assignments.final_direct_cursor == plan.k_total,
              "reconciliation: total consumed == allocated for both susceptible and direct banks");

        // A second traversal of the same cursor is rejected (no refill path).
        check(throws_runtime_error([&] { assign_family_replicate(cursor, plan); }),
              "re-traversing an already-traversed cursor is rejected");

        // Distinct, stable ids across every direct-bank entry.
        FamilyParameterBanks fresh_banks = build_family_parameter_banks(
            seeds, config, 0, plan, category, QerFamily::Mikszewski, kDirectParams);
        const FamilyReplicateAssignments fresh = assign_family_replicate(fresh_banks, plan);
        std::vector<std::uint64_t> ids;
        for (const auto& row : fresh.simulations) {
            for (const auto& infector : row.infectors) ids.push_back(infector.id);
        }
        bool all_unique = true;
        for (std::size_t i = 0; i < ids.size() && all_unique; ++i) {
            for (std::size_t j = i + 1; j < ids.size(); ++j) {
                if (ids[i] == ids[j]) all_unique = false;
            }
        }
        check(all_unique && ids.size() == plan.k_total,
              "every consumed direct-bank id is unique and there are exactly K_total of them");
    }

    // ================================================================
    // Rejection of invalid/mixed-family requests BEFORE any allocation.
    // ================================================================
    {
        const SamplingPlan plan = make_plan({1, 2});
        check(throws_invalid_argument([&] {
                  build_family_parameter_banks(seeds, config, 0, plan, category,
                                                static_cast<QerFamily>(99), kDirectParams);
              }),
              "an out-of-band QerFamily value is rejected before any allocation");
        check(throws_invalid_argument([&] { compute_family_bank_budget(static_cast<QerFamily>(99), 1, 1, 1); }),
              "compute_family_bank_budget rejects an out-of-band QerFamily value");

        // A direct family with a non-positive or nonfinite sigma is rejected
        // up front, before any bank is generated.
        check(throws_invalid_argument([&] {
                  build_family_parameter_banks(seeds, config, 0, plan, category, QerFamily::JonesFitted,
                                                DirectQerParams{-1.0, 0.0});
              }),
              "a direct family with non-positive sigma_log10 is rejected before allocation");
        check(throws_invalid_argument([&] {
                  build_family_parameter_banks(seeds, config, 0, plan, category, QerFamily::Mikszewski,
                                                DirectQerParams{std::numeric_limits<double>::quiet_NaN(), 1.0});
              }),
              "a direct family with nonfinite mu_log10 is rejected before allocation");

        // Rejection happens before allocation: no bank should be populated
        // even transiently. We cannot observe intermediate state directly,
        // but we can confirm the call throws deterministically and does not
        // partially succeed by calling it twice and getting the same error
        // type both times (no side effect that would change behavior).
        check(throws_invalid_argument([&] {
                  build_family_parameter_banks(seeds, config, 0, plan, category, QerFamily::JonesFitted,
                                                DirectQerParams{-1.0, -2.0});
              }) && throws_invalid_argument([&] {
                  build_family_parameter_banks(seeds, config, 0, plan, category, QerFamily::JonesFitted,
                                                DirectQerParams{-1.0, -2.0});
              }),
              "an invalid direct-family request is rejected deterministically and repeatedly");
    }

    // ================================================================
    // Reproducibility: same seed + same plan => identical bank contents and
    // identical consumption order.
    // ================================================================
    {
        const SamplingPlan plan = make_plan({0, 2, 1, 4, 0});
        for (QerFamily family : {QerFamily::Jones, QerFamily::JonesFitted, QerFamily::Mikszewski}) {
            const FamilyParameterBanks first =
                build_family_parameter_banks(seeds, config, 7, plan, category, family, kDirectParams);
            const FamilyParameterBanks second =
                build_family_parameter_banks(seeds, config, 7, plan, category, family, kDirectParams);
            for (std::size_t i = 0; i < 3; ++i) {
                check(first.susceptible_banks[i].values == second.susceptible_banks[i].values &&
                          first.susceptible_banks[i].seed == second.susceptible_banks[i].seed,
                      "reproducibility: susceptible bank " + std::to_string(i) + " identical for " +
                          ecai_cpp2::qer_direct::to_string(family));
            }
            if (family == QerFamily::Jones) {
                for (std::size_t i = 0; i < 8; ++i) {
                    check(first.infector_banks[i].values == second.infector_banks[i].values,
                          "reproducibility: Jones infector bank " + std::to_string(i) + " identical");
                }
            } else {
                check(first.direct_qer_bank.values == second.direct_qer_bank.values &&
                          first.direct_qer_bank.seed == second.direct_qer_bank.seed,
                      "reproducibility: direct-QER bank identical for " +
                          ecai_cpp2::qer_direct::to_string(family));
            }

            const FamilyReplicateAssignments a1 = assign_family_replicate(first, plan);
            const FamilyReplicateAssignments a2 = assign_family_replicate(second, plan);
            check(a1.simulations.size() == a2.simulations.size(), "reproducibility: same row count");
            for (std::size_t i = 0; i < a1.simulations.size(); ++i) {
                check(a1.simulations[i].infectors.size() == a2.simulations[i].infectors.size(),
                      "reproducibility: identical consumption order (row " + std::to_string(i) + ")");
                for (std::size_t k = 0; k < a1.simulations[i].infectors.size(); ++k) {
                    check(a1.simulations[i].infectors[k].id == a2.simulations[i].infectors[k].id,
                          "reproducibility: identical consumption id");
                }
            }
        }
    }

    // ================================================================
    // Full preservation of existing Jones layout/order/streams/
    // reproducibility/default behavior: build_parameter_banks/
    // ParameterCursor/assign_replicate remain usable and unaffected by the
    // new family-aware path.
    // ================================================================
    {
        const SamplingPlan plan = make_plan({2, 0, 1, 3});
        const ParameterBanks banks = build_parameter_banks(seeds, config, 0, plan, category);
        check(banks.bank(BankId::susceptible_pbr).values.size() == plan.n_active,
              "existing build_parameter_banks: susceptible bank length unaffected");
        check(banks.bank(BankId::infector_pbr).values.size() == plan.k_total,
              "existing build_parameter_banks: infector bank length unaffected");
        const ReplicateAssignments assignments = assign_replicate(banks, plan);
        check(assignments.final_cursor(BankId::DK) == plan.k_total,
              "existing assign_replicate: final cursor reconciliation unaffected");
    }

    if (g_failures > 0) {
        std::cerr << g_failures << " test(s) failed.\n";
        return 1;
    }
    std::cout << "All family bank allocation tests passed.\n";
    return 0;
}
