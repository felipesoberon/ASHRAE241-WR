// Stage 6 tests for sampling_plan: binomial pre-sampling and exact
// draw-count plan, per docs/MODEL_SPECIFICATION.md section 11.4 (resource /
// memory budget policy) and IMPLEMENTATION_PLAN.md Section 7, Stage 6.
//
// Every fixture below either (a) hand-reconciles a small counts vector
// through the pure summarize_counts() helper, independent of any RNG/LHS
// pipeline, or (b) independently re-derives the seed and re-runs
// lhs::generate_lhs_1d + distributions::binomial_ppf directly (the same
// already-verified Stage 3/4/5 primitives sampling_plan.cpp itself calls)
// to confirm build_sampling_plan wires them together with the documented
// StreamId::infector_count / list_length = n convention, per Stage 6's
// "Independent verification: Confirm binomial transform convention and
// manually reconcile several small fixtures."

#include "sampling_plan.h"

#include "category_data.h"
#include "distributions.h"
#include "lhs.h"
#include "seed_manager.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
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

using ecai_cpp2::sampling_plan::build_sampling_plan;
using ecai_cpp2::sampling_plan::CountSummary;
using ecai_cpp2::sampling_plan::ResourcePolicy;
using ecai_cpp2::sampling_plan::SamplingPlan;
using ecai_cpp2::sampling_plan::summarize_counts;
using ecai_cpp2::seed_manager::EffectiveConfig;
using ecai_cpp2::seed_manager::SeedManager;
using ecai_cpp2::seed_manager::StreamId;

EffectiveConfig make_config() {
    EffectiveConfig config;
    config.category_name = "Patient room";
    config.i0 = 400;
    config.volume_m3 = 3000.0;
    config.mask_efficiency = 0.3;
    config.community_infection_rate = 0.03;
    config.target_probability = 0.001;
    config.duration_hours = 1.0;
    config.percentile = 0.96;
    config.confidence_level = 0.95;
    config.simulations_n = 10000;
    config.replicates_r = 500;
    return config;
}

// A hand-written OccupancyCategory fixture holding the Jones et al. Table 2 /
// Table 3 values for "Patient room" (docs/MODEL_SPECIFICATION.md sections 6
// and 7). Written out literally rather than looked up so these Stage 7 tests
// stay independent of the Stage 2 lookup implementation.
ecai_cpp2::OccupancyCategory make_test_category() {
    ecai_cpp2::OccupancyCategory category{};
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

// Builds a fully consistent SamplingPlan directly from a hand-written counts
// vector, without running the Stage 6 RNG/LHS/binomial pipeline, so Stage 7
// fixtures can pin exact N_active / K_total values.
SamplingPlan make_plan(const std::vector<long long>& counts) {
    const CountSummary summary = summarize_counts(counts);
    SamplingPlan plan;
    plan.counts = counts;
    plan.n = counts.size();
    plan.n_zero = summary.n_zero;
    plan.n_active = summary.n_active;
    plan.k_total = summary.k_total;
    plan.susceptible_draws = summary.n_active;
    plan.infector_draws = summary.k_total;
    return plan;
}

// The frozen Stage 5 inverse transform for each bank, written out here
// independently of src/sampling_plan.cpp so a wiring change in the
// implementation cannot silently redefine what the test expects.
double reference_transform(ecai_cpp2::sampling_plan::BankId bank, double u,
                            const ecai_cpp2::OccupancyCategory& category) {
    namespace d = ecai_cpp2::distributions;
    using ecai_cpp2::sampling_plan::BankId;
    switch (bank) {
        case BankId::susceptible_pbr:
        case BankId::infector_pbr:
            return d::lognormal_gm_gsd_ppf(u, category.pbr_gm_m3h, category.pbr_gsd);
        case BankId::gamma:
            return d::uniform_ppf(u, 0.42, 0.61);
        case BankId::lambda:
            return d::lognormal_gm_gsd_ppf(u, 0.52, 1.9);
        case BankId::c_drop:
            return d::lognormal_gm_gsd_ppf(u, category.cdrop_gm_e4_per_m3, category.cdrop_gsd);
        case BankId::d_bar:
            return d::lognormal_gm_gsd_ppf(u, category.d_gm_e6_m, category.d_gsd);
        case BankId::E:
            return d::beta_ppf(u, 5.0, 2.0, 2.0, 5.0);
        case BankId::GVL:
            return d::log10normal_ppf(u, 7.0, 1.4);
        case BankId::VF:
            return d::beta_ppf(u, 5.0, 2.0, 1e-4, 1e-2);
        case BankId::RTD:
            return d::uniform_ppf(u, 0.43, 0.65);
        case BankId::DK:
            return d::uniform_ppf(u, 5.0, 15.0);
    }
    throw std::invalid_argument("reference_transform: unknown BankId");
}

// The frozen hard support bounds of each bank's distribution
// (docs/MODEL_SPECIFICATION.md sections 5 and 7). Unbounded-above lognormal
// families are only required to be strictly positive.
bool in_frozen_domain(ecai_cpp2::sampling_plan::BankId bank, double v) {
    using ecai_cpp2::sampling_plan::BankId;
    switch (bank) {
        case BankId::susceptible_pbr:
        case BankId::infector_pbr:
        case BankId::lambda:
        case BankId::c_drop:
        case BankId::d_bar:
        case BankId::GVL:
            return v > 0.0;
        case BankId::gamma:
            return v >= 0.42 && v <= 0.61;
        case BankId::E:
            return v >= 2.0 && v <= 5.0;
        case BankId::VF:
            return v >= 1e-4 && v <= 1e-2;
        case BankId::RTD:
            return v >= 0.43 && v <= 0.65;
        case BankId::DK:
            return v >= 5.0 && v <= 15.0;
    }
    return false;
}

// Rank of each element within its own vector (0 = smallest), used to compare
// two same-size banks for rank-sequence equality.
std::vector<std::size_t> rank_sequence(const std::vector<double>& values) {
    std::vector<std::size_t> order(values.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(),
              [&values](std::size_t a, std::size_t b) { return values[a] < values[b]; });
    std::vector<std::size_t> ranks(values.size());
    for (std::size_t r = 0; r < order.size(); ++r) ranks[order[r]] = r;
    return ranks;
}

// Every consumption id that a completed traversal actually handed out, in
// assignment order.
std::vector<std::uint64_t> collect_consumed_ids(
    const ecai_cpp2::sampling_plan::ReplicateAssignments& assignments) {
    std::vector<std::uint64_t> ids;
    for (const auto& row : assignments.simulations) {
        if (!row.active) continue;
        for (std::uint64_t id : row.susceptible.ids) ids.push_back(id);
        for (const auto& infector : row.infectors) {
            for (std::uint64_t id : infector.ids) ids.push_back(id);
        }
    }
    return ids;
}

// Every consumption id that exists across all 11 generated banks.
std::vector<std::uint64_t> collect_generated_ids(
    const ecai_cpp2::sampling_plan::ParameterBanks& banks) {
    std::vector<std::uint64_t> ids;
    for (ecai_cpp2::sampling_plan::BankId id : ecai_cpp2::sampling_plan::kAllBankIds) {
        const auto& bank = banks.bank(id);
        for (std::size_t i = 0; i < bank.values.size(); ++i) ids.push_back(bank.value_id(i));
    }
    return ids;
}

std::size_t sum_of_positive_counts(const std::vector<long long>& counts) {
    std::size_t total = 0;
    for (long long c : counts) {
        if (c > 0) total += static_cast<std::size_t>(c);
    }
    return total;
}

}  // namespace

int main() {
    const EffectiveConfig config = make_config();
    const SeedManager seeds(4242);

    // ---------------------------------------------------------------
    // summarize_counts: pure helper, hand-reconciled mixed fixture.
    // counts = {0, 0, 3, 1, 0, 2, 5}: 3 zeros, 4 active, K_total = 11.
    // ---------------------------------------------------------------
    {
        const std::vector<long long> counts = {0, 0, 3, 1, 0, 2, 5};
        const CountSummary summary = summarize_counts(counts);
        check(summary.n_zero == 3, "summarize_counts: hand-reconciled mixed fixture n_zero == 3");
        check(summary.n_active == 4, "summarize_counts: hand-reconciled mixed fixture n_active == 4");
        check(summary.k_total == 11, "summarize_counts: hand-reconciled mixed fixture k_total == 11");
        check(summary.n_zero + summary.n_active == counts.size(),
              "summarize_counts: hand-reconciled mixed fixture n_zero + n_active == counts.size()");
    }

    // ---------------------------------------------------------------
    // summarize_counts: all-zero fixture.
    // ---------------------------------------------------------------
    {
        const std::vector<long long> counts(9, 0);
        const CountSummary summary = summarize_counts(counts);
        check(summary.n_zero == 9, "summarize_counts: all-zero fixture n_zero == 9");
        check(summary.n_active == 0, "summarize_counts: all-zero fixture n_active == 0");
        check(summary.k_total == 0, "summarize_counts: all-zero fixture k_total == 0");
    }

    // ---------------------------------------------------------------
    // summarize_counts: all-active fixture, every count equal to a fixed
    // I0-like value.
    // ---------------------------------------------------------------
    {
        const std::vector<long long> counts(6, 7);
        const CountSummary summary = summarize_counts(counts);
        check(summary.n_zero == 0, "summarize_counts: all-active fixture n_zero == 0");
        check(summary.n_active == 6, "summarize_counts: all-active fixture n_active == 6");
        check(summary.k_total == 42, "summarize_counts: all-active fixture k_total == 6*7 == 42");
    }

    // ---------------------------------------------------------------
    // summarize_counts: rejects a negative count.
    // ---------------------------------------------------------------
    {
        bool threw = false;
        try {
            summarize_counts({1, -1, 2});
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "summarize_counts: rejects a negative count with std::invalid_argument");
    }

    // ---------------------------------------------------------------
    // summarize_counts: checked K_total accumulation overflow. Three
    // LLONG_MAX-valued counts sum past std::size_t's maximum on a 64-bit
    // platform (2 * LLONG_MAX < SIZE_MAX < 3 * LLONG_MAX).
    // ---------------------------------------------------------------
    {
        const long long huge = std::numeric_limits<long long>::max();
        const std::vector<long long> counts = {huge, huge, huge};
        bool threw = false;
        try {
            summarize_counts(counts);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw,
              "summarize_counts: checked K_total accumulation overflow throws std::runtime_error");
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: exactly n counts, every count in [0, I0], and
    // the N_zero + N_active == N / K_total == sum(counts) identities, for a
    // small n so results can be sanity-checked by eye alongside the
    // programmatic checks.
    // ---------------------------------------------------------------
    {
        const std::size_t n = 25;
        const long long i0 = 10;
        const double p = 0.4;
        const SamplingPlan plan = build_sampling_plan(seeds, config, 0, n, i0, p);

        check(plan.counts.size() == n, "build_sampling_plan: returns exactly n counts");
        check(plan.n == n, "build_sampling_plan: plan.n == n");

        bool all_in_range = true;
        for (long long c : plan.counts) {
            if (c < 0 || c > i0) all_in_range = false;
        }
        check(all_in_range, "build_sampling_plan: every count is in [0, I0]");

        check(plan.n_zero + plan.n_active == n,
              "build_sampling_plan: N_zero + N_active == N");
        check(plan.k_total == sum_of_positive_counts(plan.counts),
              "build_sampling_plan: K_total == sum(ninf)");
        check(plan.susceptible_draws == plan.n_active,
              "build_sampling_plan: susceptible_draws == N_active");
        check(plan.infector_draws == plan.k_total,
              "build_sampling_plan: infector_draws == K_total");

        // Cross-check against the pure summarize_counts helper.
        const CountSummary summary = summarize_counts(plan.counts);
        check(summary.n_zero == plan.n_zero && summary.n_active == plan.n_active &&
                  summary.k_total == plan.k_total,
              "build_sampling_plan: plan accounting matches summarize_counts(plan.counts)");
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: all-zero synthetic case (p == 0). Every count is
    // 0 regardless of I0; both susceptible and infector draw counts are
    // zero, and zero-infector rows imply no other draws.
    // ---------------------------------------------------------------
    {
        const SamplingPlan plan = build_sampling_plan(seeds, config, 0, 40, 250, 0.0);
        bool all_zero = true;
        for (long long c : plan.counts) {
            if (c != 0) all_zero = false;
        }
        check(all_zero, "build_sampling_plan: p=0 yields every count == 0");
        check(plan.n_zero == 40, "build_sampling_plan: p=0 yields N_zero == N");
        check(plan.n_active == 0, "build_sampling_plan: p=0 yields N_active == 0");
        check(plan.k_total == 0, "build_sampling_plan: p=0 yields K_total == 0");
        check(plan.susceptible_draws == 0,
              "build_sampling_plan: p=0 yields susceptible_draws == 0");
        check(plan.infector_draws == 0, "build_sampling_plan: p=0 yields infector_draws == 0");
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: all-infected case (p == 1). Every count equals
    // I0; N_zero == 0, N_active == N, K_total == N * I0.
    // ---------------------------------------------------------------
    {
        const std::size_t n = 17;
        const long long i0 = 9;
        const SamplingPlan plan = build_sampling_plan(seeds, config, 0, n, i0, 1.0);
        bool all_full = true;
        for (long long c : plan.counts) {
            if (c != i0) all_full = false;
        }
        check(all_full, "build_sampling_plan: p=1 yields every count == I0");
        check(plan.n_zero == 0, "build_sampling_plan: p=1 yields N_zero == 0");
        check(plan.n_active == n, "build_sampling_plan: p=1 yields N_active == N");
        check(plan.k_total == n * static_cast<std::size_t>(i0),
              "build_sampling_plan: p=1 yields K_total == N * I0");
        check(plan.susceptible_draws == n,
              "build_sampling_plan: p=1 yields one susceptible draw per active simulation");
        check(plan.infector_draws == n * static_cast<std::size_t>(i0),
              "build_sampling_plan: p=1 yields one infector draw per realized infector");
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: N=10000, I0=400, p=0.01 (the publication default
    // scale). Structural invariants plus a loose, non-flaky sanity bound on
    // the mean count against the analytic binomial mean I0*p = 4.
    // ---------------------------------------------------------------
    {
        const std::size_t n = 10000;
        const long long i0 = 400;
        const double p = 0.01;
        const SamplingPlan plan = build_sampling_plan(seeds, config, 0, n, i0, p);

        check(plan.counts.size() == n, "build_sampling_plan(N=10000,I0=400,p=.01): returns N counts");
        bool all_in_range = true;
        for (long long c : plan.counts) {
            if (c < 0 || c > i0) all_in_range = false;
        }
        check(all_in_range, "build_sampling_plan(N=10000,I0=400,p=.01): every count is in [0, I0]");
        check(plan.n_zero + plan.n_active == n,
              "build_sampling_plan(N=10000,I0=400,p=.01): N_zero + N_active == N");
        check(plan.k_total == sum_of_positive_counts(plan.counts),
              "build_sampling_plan(N=10000,I0=400,p=.01): K_total == sum(ninf)");
        check(plan.susceptible_draws == plan.n_active,
              "build_sampling_plan(N=10000,I0=400,p=.01): susceptible draws equal N_active");
        check(plan.infector_draws == plan.k_total,
              "build_sampling_plan(N=10000,I0=400,p=.01): infector draws equal K_total");

        const double mean = static_cast<double>(plan.k_total) / static_cast<double>(n);
        check(mean > 1.0 && mean < 10.0,
              "build_sampling_plan(N=10000,I0=400,p=.01): mean count is a plausible, "
              "non-flaky bound around the analytic mean I0*p=4");
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: stream wiring / binomial transform convention,
    // independent reconciliation. Manually re-derive the seed via
    // StreamId::infector_count with list_length == n, regenerate the LHS
    // list, and transform through distributions::binomial_ppf directly,
    // then compare bit-for-bit against build_sampling_plan's output.
    // ---------------------------------------------------------------
    {
        const std::size_t n = 13;
        const long long i0 = 6;
        const double p = 0.35;
        const SamplingPlan plan = build_sampling_plan(seeds, config, 2, n, i0, p);

        const std::vector<long long> pinned_counts = {2, 2, 2, 2, 6, 3, 1, 0, 1, 4, 3, 2, 1};
        check(plan.counts == pinned_counts,
              "build_sampling_plan: fixed seed/config/replicate has the pinned count vector");

        const std::uint64_t expected_seed =
            seeds.derive_seed(config, 2, StreamId::infector_count, n);
        const std::vector<double> expected_uniforms =
            ecai_cpp2::lhs::generate_lhs_1d(n, expected_seed);
        std::vector<long long> expected_counts;
        expected_counts.reserve(n);
        for (double u : expected_uniforms) {
            expected_counts.push_back(ecai_cpp2::distributions::binomial_ppf(u, i0, p));
        }
        check(plan.counts == expected_counts,
              "build_sampling_plan: matches an independent StreamId::infector_count / "
              "list_length=n reconciliation, bit-for-bit");
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: determinism. Same (config, replicate, n, I0, p)
    // tuple reproduces an identical plan.
    // ---------------------------------------------------------------
    {
        const SamplingPlan a = build_sampling_plan(seeds, config, 3, 30, 50, 0.2);
        const SamplingPlan b = build_sampling_plan(seeds, config, 3, 30, 50, 0.2);
        check(a.counts == b.counts,
              "build_sampling_plan: identical tuples produce bit-for-bit identical counts");
        check(a.n_zero == b.n_zero && a.n_active == b.n_active && a.k_total == b.k_total,
              "build_sampling_plan: identical tuples produce identical accounting");
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: cross-replicate difference. Changing only
    // replicate_index changes the resulting counts.
    // ---------------------------------------------------------------
    {
        const SamplingPlan a = build_sampling_plan(seeds, config, 4, 30, 50, 0.2);
        const SamplingPlan b = build_sampling_plan(seeds, config, 5, 30, 50, 0.2);
        check(a.counts != b.counts,
              "build_sampling_plan: different replicate_index yields different counts");
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: invalid inputs.
    // ---------------------------------------------------------------
    {
        bool threw = false;
        try {
            build_sampling_plan(seeds, config, 0, 0, 10, 0.5);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "build_sampling_plan: n=0 throws std::invalid_argument");
    }
    {
        bool threw = false;
        try {
            build_sampling_plan(seeds, config, 0, 10, -1, 0.5);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "build_sampling_plan: negative I0 throws std::invalid_argument");
    }
    {
        bool threw = false;
        try {
            build_sampling_plan(seeds, config, 0, 10, 10, std::numeric_limits<double>::quiet_NaN());
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "build_sampling_plan: p=NaN throws std::invalid_argument");
    }
    {
        bool threw = false;
        try {
            build_sampling_plan(seeds, config, 0, 10, 10, -0.1);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "build_sampling_plan: p<0 throws std::invalid_argument");
    }
    {
        bool threw = false;
        try {
            build_sampling_plan(seeds, config, 0, 10, 10, 1.1);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "build_sampling_plan: p>1 throws std::invalid_argument");
    }
    {
        bool threw = false;
        try {
            build_sampling_plan(seeds, config, 0, 10, 10, std::numeric_limits<double>::infinity());
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "build_sampling_plan: p=+inf throws std::invalid_argument");
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: resource-policy rejection, N exceeds
    // max_simulations, before any allocation.
    // ---------------------------------------------------------------
    {
        ResourcePolicy policy;
        policy.max_simulations = 5;
        bool threw = false;
        try {
            build_sampling_plan(seeds, config, 0, 6, 10, 0.5, policy);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw,
              "build_sampling_plan: n exceeding ResourcePolicy::max_simulations throws "
              "std::runtime_error");
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: resource-policy rejection, N*I0 exceeds
    // max_total_infectors (default policy), before any allocation.
    // ---------------------------------------------------------------
    {
        bool threw = false;
        try {
            build_sampling_plan(seeds, config, 0, 1000, 1000000, 0.5);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw,
              "build_sampling_plan: N*I0 exceeding the default ResourcePolicy::"
              "max_total_infectors throws std::runtime_error");
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: checked N*I0 multiplication overflow of
    // std::size_t itself (independent of any policy ceiling — the policy is
    // widened to SIZE_MAX so only the overflow guard can fire).
    // ---------------------------------------------------------------
    {
        ResourcePolicy policy;
        policy.max_simulations = std::numeric_limits<std::size_t>::max();
        policy.max_total_infectors = std::numeric_limits<std::size_t>::max();
        policy.max_bytes = std::numeric_limits<std::size_t>::max();
        const std::size_t huge_n = std::numeric_limits<std::size_t>::max() / 2 + 10;
        bool threw = false;
        try {
            build_sampling_plan(seeds, config, 0, huge_n, 3, 0.5, policy);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw, "build_sampling_plan: N*I0 overflowing std::size_t throws std::runtime_error "
                     "before any allocation");
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: resource-policy rejection, required byte count
    // exceeds a tightened max_bytes, before any allocation.
    // ---------------------------------------------------------------
    {
        ResourcePolicy policy;
        policy.max_bytes = 1000;
        bool threw = false;
        try {
            build_sampling_plan(seeds, config, 0, 1000000, 1, 0.5, policy);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw,
              "build_sampling_plan: required bytes exceeding a tightened "
              "ResourcePolicy::max_bytes throws std::runtime_error");
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: max_bytes preflight must account for the
    // simultaneously-live binomial PMF buffer ((I0+1) doubles), not just
    // the N-length uniform/count vectors. N=1, I0=100, p=.5 needs only 16
    // bytes for the N-length vectors but 808 additional bytes per binomial
    // buffer; with max_simulations/max_total_infectors widened so only the
    // byte check can fire, and max_bytes=100 (below the true peak but above
    // the under-counted 16-byte figure), this must throw
    // std::runtime_error before any allocation rather than silently
    // succeeding.
    // ---------------------------------------------------------------
    {
        ResourcePolicy policy;
        policy.max_simulations = std::numeric_limits<std::size_t>::max();
        policy.max_total_infectors = std::numeric_limits<std::size_t>::max();
        policy.max_bytes = 100;
        bool threw = false;
        try {
            build_sampling_plan(seeds, config, 0, 1, 100, 0.5, policy);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw,
              "build_sampling_plan: max_bytes preflight accounts for the live binomial PMF "
              "buffer (N=1,I0=100,p=.5,max_bytes=100) and throws std::runtime_error before any "
              "allocation");
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: max_bytes preflight must also account for the
    // BinomialPpfTable's retained CDF buffer, not just the PMF buffer.
    // N=1, I0=100 needs 16 bytes for the N-length vectors plus 2*808 bytes
    // for the PMF and CDF buffers (1632-byte peak); max_bytes=1000 is above
    // the PMF-only 824-byte figure but below the true peak, so this must
    // throw std::runtime_error before any allocation.
    // ---------------------------------------------------------------
    {
        ResourcePolicy policy;
        policy.max_simulations = std::numeric_limits<std::size_t>::max();
        policy.max_total_infectors = std::numeric_limits<std::size_t>::max();
        policy.max_bytes = 1000;
        bool threw = false;
        try {
            build_sampling_plan(seeds, config, 0, 1, 100, 0.5, policy);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw,
              "build_sampling_plan: max_bytes preflight accounts for the binomial CDF buffer "
              "(N=1,I0=100,p=.5,max_bytes=1000) and throws std::runtime_error before any "
              "allocation");
        policy.max_bytes = 1632;
        const SamplingPlan plan = build_sampling_plan(seeds, config, 0, 1, 100, 0.5, policy);
        check(plan.n == 1 && plan.counts.size() == 1,
              "build_sampling_plan: max_bytes exactly equal to the 16N + 16(I0+1) peak succeeds");
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: per-replicate BinomialPpfTable equivalence at
    // production scale, including a p whose pmf(0) underflows. Every count
    // must match the single-shot binomial_ppf transform of the same LHS
    // list, bit-for-bit.
    // ---------------------------------------------------------------
    {
        struct Case { std::size_t n; long long i0; double p; };
        const Case cases[] = {{5000, 400, 0.01}, {5000, 400, 0.9}, {2000, 400, 0.99}};
        for (const Case& c : cases) {
            const SamplingPlan plan = build_sampling_plan(seeds, config, 7, c.n, c.i0, c.p);
            const std::uint64_t seed = seeds.derive_seed(config, 7, StreamId::infector_count, c.n);
            const std::vector<double> uniforms = ecai_cpp2::lhs::generate_lhs_1d(c.n, seed);
            bool equal = plan.counts.size() == uniforms.size();
            for (std::size_t i = 0; equal && i < uniforms.size(); ++i) {
                equal = plan.counts[i] ==
                        ecai_cpp2::distributions::binomial_ppf(uniforms[i], c.i0, c.p);
            }
            check(equal, "build_sampling_plan(N=" + std::to_string(c.n) + ",I0=" +
                             std::to_string(c.i0) + ",p=" + std::to_string(c.p) +
                             "): table-driven counts equal single-shot binomial_ppf bit-for-bit");
        }
    }

    // ---------------------------------------------------------------
    // build_sampling_plan: I0+1 (the binomial_ppf PMF-buffer feasibility
    // check) exceeding std::vector<double>::max_size(), with the N*I0 and
    // byte-budget policy ceilings widened so only this specific guard can
    // fire, before any allocation.
    // ---------------------------------------------------------------
    {
        ResourcePolicy policy;
        policy.max_simulations = std::numeric_limits<std::size_t>::max();
        policy.max_total_infectors = std::numeric_limits<std::size_t>::max();
        policy.max_bytes = std::numeric_limits<std::size_t>::max();
        const long long huge_i0 = std::numeric_limits<long long>::max();
        bool threw = false;
        try {
            build_sampling_plan(seeds, config, 0, 1, huge_i0, 0.5, policy);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw,
              "build_sampling_plan: I0 exceeding the binomial_ppf PMF-buffer feasibility bound "
              "throws std::runtime_error before any allocation");
    }

    // ===============================================================
    // Stage 7 (IMPLEMENTATION_PLAN.md Section 7, Stage 7): exact-size
    // parameter banks and cursor accounting.
    //
    // Bank identity: the 11 named physical parameter banks of
    // MODEL_SPECIFICATION.md section 11.1 (3 susceptible-sized, 8
    // infector-sized), each bound to its own frozen Stage 4 StreamId.
    // ===============================================================
    {
        using ecai_cpp2::sampling_plan::bank_id_name;
        using ecai_cpp2::sampling_plan::bank_stream_id;
        using ecai_cpp2::sampling_plan::BankId;
        using ecai_cpp2::sampling_plan::is_infector_bank;
        using ecai_cpp2::sampling_plan::kAllBankIds;

        check(kAllBankIds.size() == 11,
              "Stage 7 banks: exactly 11 named parameter banks (3 susceptible + 8 infector)");

        const std::vector<std::pair<BankId, StreamId>> expected = {
            {BankId::susceptible_pbr, StreamId::susceptible_pbr},
            {BankId::gamma, StreamId::gamma},
            {BankId::lambda, StreamId::lambda},
            {BankId::infector_pbr, StreamId::infector_pbr},
            {BankId::c_drop, StreamId::c_drop},
            {BankId::d_bar, StreamId::d_bar},
            {BankId::E, StreamId::E},
            {BankId::GVL, StreamId::GVL},
            {BankId::VF, StreamId::VF},
            {BankId::RTD, StreamId::RTD},
            {BankId::DK, StreamId::DK},
        };
        bool mapping_ok = expected.size() == kAllBankIds.size();
        for (std::size_t i = 0; mapping_ok && i < expected.size(); ++i) {
            mapping_ok = kAllBankIds[i] == expected[i].first &&
                          bank_stream_id(expected[i].first) == expected[i].second;
        }
        check(mapping_ok,
              "Stage 7 banks: kAllBankIds order and per-bank StreamId mapping match the frozen "
              "Stage 4 stream ids");

        // No bank may reuse StreamId::infector_count (Stage 6's own stream),
        // and no two banks may share a stream id.
        std::vector<std::uint64_t> stream_values;
        bool none_is_infector_count = true;
        for (BankId id : kAllBankIds) {
            const StreamId sid = bank_stream_id(id);
            none_is_infector_count =
                none_is_infector_count && sid != StreamId::infector_count;
            stream_values.push_back(static_cast<std::uint64_t>(sid));
        }
        check(none_is_infector_count,
              "Stage 7 banks: no parameter bank reuses StreamId::infector_count");
        std::sort(stream_values.begin(), stream_values.end());
        check(std::adjacent_find(stream_values.begin(), stream_values.end()) ==
                  stream_values.end(),
              "Stage 7 banks: all 11 bank StreamIds are distinct");

        // Susceptible/infector classification drives bank length (N_active
        // vs K_total) in MODEL_SPECIFICATION.md section 11.1.
        std::size_t infector_banks = 0;
        for (BankId id : kAllBankIds) {
            if (is_infector_bank(id)) ++infector_banks;
        }
        check(infector_banks == 8,
              "Stage 7 banks: exactly 8 banks are infector-sized, the other 3 susceptible-sized");
        check(!is_infector_bank(BankId::susceptible_pbr) && !is_infector_bank(BankId::gamma) &&
                  !is_infector_bank(BankId::lambda) && is_infector_bank(BankId::infector_pbr) &&
                  is_infector_bank(BankId::DK),
              "Stage 7 banks: susceptible PBR/gamma/lambda are susceptible-sized; infector PBR "
              "and DK are infector-sized");

        // Names are stable, distinct, human-readable identifiers.
        std::vector<std::string> names;
        for (BankId id : kAllBankIds) names.push_back(bank_id_name(id));
        check(names.front() == "susceptible_pbr" && names.back() == "DK",
              "Stage 7 banks: bank_id_name returns the frozen enumerator spelling");
        std::sort(names.begin(), names.end());
        check(std::adjacent_find(names.begin(), names.end()) == names.end(),
              "Stage 7 banks: all 11 bank names are distinct");

        bool threw = false;
        try {
            bank_stream_id(static_cast<BankId>(99));
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "Stage 7 banks: bank_stream_id rejects an out-of-range BankId");
        threw = false;
        try {
            bank_id_name(static_cast<BankId>(99));
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "Stage 7 banks: bank_id_name rejects an out-of-range BankId");
    }

    // ===============================================================
    // Stage 7: combined-bank budget (MODEL_SPECIFICATION.md section 11.4).
    // total_elements = N + 3*N_active + 8*K_total, total_bytes = 8 *
    // total_elements, computed with checked arithmetic and compared against
    // the ResourcePolicy ceilings before any bank exists.
    // ===============================================================
    {
        using ecai_cpp2::sampling_plan::BankBudget;
        using ecai_cpp2::sampling_plan::compute_bank_budget;

        // Hand-reconciled: 7 + 3*4 + 8*11 = 7 + 12 + 88 = 107 elements.
        const BankBudget budget = compute_bank_budget(7, 4, 11);
        check(budget.total_elements == 107,
              "bank budget: hand-reconciled N=7,N_active=4,K_total=11 gives 107 elements");
        check(budget.total_bytes == 856,
              "bank budget: hand-reconciled N=7,N_active=4,K_total=11 gives 856 bytes");

        // Degenerate all-zero-infector replicate: only the N binomial values.
        const BankBudget zero_budget = compute_bank_budget(9, 0, 0);
        check(zero_budget.total_elements == 9 && zero_budget.total_bytes == 72,
              "bank budget: all-zero-infector plan budgets only the N binomial elements");

        // Randomized synthetic triples: the formula must match an
        // independently computed reference for every case.
        std::mt19937_64 rng(20250909ULL);
        bool formula_ok = true;
        for (int trial = 0; trial < 400; ++trial) {
            const std::size_t n = 1 + static_cast<std::size_t>(rng() % 5000);
            const std::size_t n_active = static_cast<std::size_t>(rng() % (n + 1));
            const std::size_t k_total =
                n_active == 0 ? 0 : n_active + static_cast<std::size_t>(rng() % 4000);
            const BankBudget b = compute_bank_budget(n, n_active, k_total);
            const std::size_t reference = n + 3 * n_active + 8 * k_total;
            formula_ok = formula_ok && b.total_elements == reference &&
                          b.total_bytes == 8 * reference;
        }
        check(formula_ok,
              "bank budget: N + 3*N_active + 8*K_total matches an independent reference over "
              "400 randomized synthetic triples");

        // Consumption IDs reserve 56 bits for an index, so a bank length may
        // not exceed 2^56 even when caller-overridden memory ceilings allow it.
        {
            ResourcePolicy widened;
            widened.max_simulations = std::numeric_limits<std::size_t>::max();
            widened.max_total_infectors = std::numeric_limits<std::size_t>::max();
            widened.max_bytes = std::numeric_limits<std::size_t>::max();
            bool threw = false;
            try {
                compute_bank_budget(1, 1, (std::size_t{1} << 56) + 1, widened);
            } catch (const std::runtime_error&) {
                threw = true;
            }
            check(threw,
                  "bank budget: rejects a bank longer than the 56-bit consumption-ID capacity");
            bool exact_capacity_ok = false;
            try {
                const BankBudget edge =
                    compute_bank_budget(1, 1, std::size_t{1} << 56, widened);
                exact_capacity_ok = edge.total_elements ==
                                    1 + 3 + 8 * (std::size_t{1} << 56);
            } catch (const std::exception&) {
                exact_capacity_ok = false;
            }
            check(exact_capacity_ok,
                  "bank budget: accepts exactly 2^56 values whose maximum ID index is 2^56-1");
        }

        // The frozen 4 GiB default: 2^29 elements is exactly at budget.
        bool at_budget_ok = false;
        try {
            const BankBudget edge = compute_bank_budget(536870912ULL, 0, 0);
            at_budget_ok = edge.total_bytes == 4ULL * 1024 * 1024 * 1024;
        } catch (const std::exception&) {
            at_budget_ok = false;
        }
        check(at_budget_ok,
              "bank budget: exactly 2^29 elements (4 GiB) is accepted under the frozen default");

        bool threw = false;
        try {
            compute_bank_budget(536870913ULL, 0, 0);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw,
              "bank budget: 2^29 + 1 elements exceeds the frozen 4 GiB default and is rejected");

        // Checked arithmetic, with every policy ceiling widened to size_t
        // max so only the arithmetic/feasibility guards can fire.
        ResourcePolicy wide;
        wide.max_simulations = std::numeric_limits<std::size_t>::max();
        wide.max_total_infectors = std::numeric_limits<std::size_t>::max();
        wide.max_bytes = std::numeric_limits<std::size_t>::max();
        const std::size_t max_bank = std::vector<double>().max_size();
        threw = false;
        try {
            compute_bank_budget(1, 1, max_bank + 1, wide);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw,
              "bank budget: a K_total bank length exceeding std::vector<double>::max_size() is "
              "rejected");

        // The products 3*N_active and 8*K_total fit, but adding them to an
        // N of SIZE_MAX cannot: the checked combined sum must catch it before
        // container-feasibility guards.
        threw = false;
        try {
            compute_bank_budget(std::numeric_limits<std::size_t>::max(), 1, 1, wide);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw,
              "bank budget: N + 3*N_active + 8*K_total overflowing std::size_t is rejected");

        // Policy K_total ceiling is enforced independently of the byte budget.
        ResourcePolicy tight_k;
        tight_k.max_total_infectors = 10;
        threw = false;
        try {
            compute_bank_budget(4, 2, 11, tight_k);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw,
              "bank budget: K_total exceeding ResourcePolicy::max_total_infectors is rejected");

        // N_active can never exceed N: an inconsistent accounting triple is
        // a caller error, not a resource-limit failure.
        threw = false;
        try {
            compute_bank_budget(3, 4, 4);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "bank budget: N_active > N is rejected as std::invalid_argument");

        // K_total < N_active is likewise impossible (each active simulation
        // has at least one infector).
        threw = false;
        try {
            compute_bank_budget(9, 4, 3);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "bank budget: K_total < N_active is rejected as std::invalid_argument");
    }

    // ===============================================================
    // Stage 7: exact-size bank construction. One vector per named physical
    // parameter, sized N_active (susceptible) or K_total (infector) per
    // MODEL_SPECIFICATION.md section 11.1, with the combined-bank budget
    // enforced before any bank is allocated.
    // ===============================================================
    {
        using ecai_cpp2::sampling_plan::BankId;
        using ecai_cpp2::sampling_plan::build_parameter_banks;
        using ecai_cpp2::sampling_plan::is_infector_bank;
        using ecai_cpp2::sampling_plan::kAllBankIds;
        using ecai_cpp2::sampling_plan::ParameterBanks;

        const ecai_cpp2::OccupancyCategory category = make_test_category();

        // Mixed plan: counts {0, 2, 0, 3, 1} -> N=5, N_active=3, K_total=6.
        const ParameterBanks banks =
            build_parameter_banks(seeds, config, 0, make_plan({0, 2, 0, 3, 1}), category);
        check(banks.n == 5 && banks.n_active == 3 && banks.k_total == 6,
              "parameter banks: mixed plan carries N=5, N_active=3, K_total=6");
        check(banks.budget.total_elements == 5 + 3 * 3 + 8 * 6,
              "parameter banks: mixed plan budget is N + 3*N_active + 8*K_total");

        bool lengths_ok = true;
        bool generated_ok = true;
        bool identity_ok = true;
        for (BankId id : kAllBankIds) {
            const std::size_t expected_length =
                is_infector_bank(id) ? banks.k_total : banks.n_active;
            const auto& bank = banks.bank(id);
            lengths_ok = lengths_ok && bank.values.size() == expected_length;
            generated_ok = generated_ok && bank.generated;
            identity_ok =
                identity_ok && bank.id == id &&
                bank.stream == ecai_cpp2::sampling_plan::bank_stream_id(id) &&
                bank.seed == seeds.derive_seed(config, 0, bank.stream,
                                                static_cast<std::uint64_t>(expected_length));
        }
        check(lengths_ok,
              "parameter banks: each susceptible bank has exactly N_active values and each "
              "infector bank exactly K_total values");
        check(generated_ok, "parameter banks: every nonzero-length bank is marked generated");
        check(identity_ok,
              "parameter banks: each bank carries its own BankId/StreamId and the seed derived "
              "from its own exact list length");

        // All-active plan.
        const ParameterBanks all_active =
            build_parameter_banks(seeds, config, 0, make_plan({1, 2, 3, 4}), category);
        bool all_active_ok = all_active.n_active == 4 && all_active.k_total == 10;
        for (BankId id : kAllBankIds) {
            all_active_ok =
                all_active_ok &&
                all_active.bank(id).values.size() == (is_infector_bank(id) ? 10u : 4u);
        }
        check(all_active_ok,
              "parameter banks: all-active plan sizes susceptible banks 4 and infector banks 10");

        // All-zero plan: every bank has length zero, and no bank is
        // generated (no LHS list is created for a zero-length bank).
        const ParameterBanks all_zero =
            build_parameter_banks(seeds, config, 0, make_plan({0, 0, 0}), category);
        bool all_zero_ok = all_zero.n == 3 && all_zero.n_active == 0 && all_zero.k_total == 0 &&
                            all_zero.budget.total_elements == 3;
        for (BankId id : kAllBankIds) {
            const auto& bank = all_zero.bank(id);
            all_zero_ok = all_zero_ok && bank.values.empty() && !bank.generated && bank.seed == 0;
        }
        check(all_zero_ok,
              "parameter banks: all-zero-infector plan yields 11 empty, ungenerated banks with "
              "no seed derived");

        // Single-infector plan: susceptible banks 1, infector banks 1.
        const ParameterBanks single =
            build_parameter_banks(seeds, config, 0, make_plan({0, 1}), category);
        bool single_ok = single.n_active == 1 && single.k_total == 1;
        for (BankId id : kAllBankIds) {
            single_ok = single_ok && single.bank(id).values.size() == 1;
        }
        check(single_ok, "parameter banks: single-infector plan yields length-1 banks");

        // Randomized synthetic count vectors: bank lengths must always equal
        // the accounting carried by the plan itself.
        std::mt19937_64 rng(777777ULL);
        bool randomized_ok = true;
        for (int trial = 0; trial < 40; ++trial) {
            std::vector<long long> counts;
            const std::size_t rows = 1 + static_cast<std::size_t>(rng() % 12);
            for (std::size_t i = 0; i < rows; ++i) {
                counts.push_back(static_cast<long long>(rng() % 5));
            }
            const SamplingPlan p = make_plan(counts);
            const ParameterBanks b = build_parameter_banks(seeds, config, 3, p, category);
            for (BankId id : kAllBankIds) {
                const std::size_t expected = is_infector_bank(id) ? p.k_total : p.n_active;
                randomized_ok = randomized_ok && b.bank(id).values.size() == expected &&
                                 b.bank(id).generated == (expected > 0);
            }
        }
        check(randomized_ok,
              "parameter banks: bank lengths match N_active/K_total over 40 randomized "
              "synthetic count vectors");

        // Budget rejection happens before any bank is produced: no partial
        // output is observable, only the thrown diagnostic.
        ResourcePolicy tiny;
        tiny.max_bytes = 64;  // 8 elements; the mixed plan needs 62.
        bool threw = false;
        bool no_partial_output = true;
        try {
            const ParameterBanks rejected = build_parameter_banks(
                seeds, config, 0, make_plan({0, 2, 0, 3, 1}), category, tiny);
            no_partial_output = rejected.bank(BankId::gamma).values.empty();
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw && no_partial_output,
              "parameter banks: a combined-bank budget overrun is rejected before any bank is "
              "allocated, with no partial output");

        // The Stage 6 working-set ceiling and the Stage 7 combined-bank
        // ceiling are separate checks: a plan that Stage 6 accepted can
        // still be rejected here.
        ResourcePolicy tight_k;
        tight_k.max_total_infectors = 5;
        threw = false;
        try {
            build_parameter_banks(seeds, config, 0, make_plan({0, 2, 0, 3, 1}), category, tight_k);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw,
              "parameter banks: the Stage 7 budget re-checks "
              "ResourcePolicy::max_total_infectors independently of Stage 6");

        // An internally inconsistent plan is a caller error.
        threw = false;
        try {
            SamplingPlan bad = make_plan({0, 2});
            bad.k_total = 99;
            build_parameter_banks(seeds, config, 0, bad, category);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw,
              "parameter banks: a plan whose accounting contradicts its counts is rejected");

        threw = false;
        try {
            SamplingPlan bad = make_plan({0, 2});
            bad.n = 7;
            build_parameter_banks(seeds, config, 0, bad, category);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "parameter banks: a plan whose n contradicts counts.size() is rejected");

        threw = false;
        try {
            banks.bank(static_cast<BankId>(42));
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "parameter banks: ParameterBanks::bank rejects an out-of-range BankId");
    }

    // ===============================================================
    // Stage 7: per-bank transforms, seeds, and independence. Each bank is
    // independently regenerated here from the same already-verified Stage
    // 3/4/5 primitives (derive_seed -> generate_lhs_1d -> the frozen
    // inverse transform) and compared against what build_parameter_banks
    // produced.
    // ===============================================================
    {
        using ecai_cpp2::sampling_plan::BankId;
        using ecai_cpp2::sampling_plan::build_parameter_banks;
        using ecai_cpp2::sampling_plan::is_infector_bank;
        using ecai_cpp2::sampling_plan::kAllBankIds;
        using ecai_cpp2::sampling_plan::ParameterBanks;

        const ecai_cpp2::OccupancyCategory category = make_test_category();

        // Equal-length banks (N_active == K_total == 8) so all 11 banks are
        // directly comparable to one another.
        const SamplingPlan plan = make_plan({1, 1, 1, 1, 1, 1, 1, 1});
        const ParameterBanks banks = build_parameter_banks(seeds, config, 11, plan, category);

        // Independent regeneration of every bank, transform included.
        bool transforms_ok = true;
        bool strata_ok = true;
        for (BankId id : kAllBankIds) {
            const auto& bank = banks.bank(id);
            const std::uint64_t expected_seed = seeds.derive_seed(
                config, 11, ecai_cpp2::sampling_plan::bank_stream_id(id),
                static_cast<std::uint64_t>(bank.values.size()));
            const std::vector<double> uniforms =
                ecai_cpp2::lhs::generate_lhs_1d(bank.values.size(), expected_seed);
            std::vector<double> expected;
            for (double u : uniforms) {
                expected.push_back(reference_transform(id, u, category));
            }
            transforms_ok = transforms_ok && bank.seed == expected_seed &&
                             bank.values.size() == expected.size() && bank.values == expected;

            // Exactly one LHS uniform per stratum, checked on the
            // regenerated uniforms before any transform is applied.
            std::vector<int> occupancy(uniforms.size(), 0);
            for (double u : uniforms) {
                const std::size_t stratum = static_cast<std::size_t>(
                    static_cast<double>(uniforms.size()) * u);
                if (stratum >= uniforms.size()) {
                    strata_ok = false;
                } else {
                    ++occupancy[stratum];
                }
            }
            for (int count : occupancy) {
                strata_ok = strata_ok && count == 1;
            }
        }
        check(transforms_ok,
              "bank transforms: every bank equals an independently regenerated "
              "derive_seed -> generate_lhs_1d -> frozen inverse transform sequence");
        check(strata_ok,
              "bank transforms: each bank draws exactly one LHS uniform per stratum before the "
              "transform is applied");

        // All 11 derived seeds are distinct even though all 11 banks share
        // the same list length (only the stream id differs).
        std::vector<std::uint64_t> bank_seeds;
        for (BankId id : kAllBankIds) bank_seeds.push_back(banks.bank(id).seed);
        std::sort(bank_seeds.begin(), bank_seeds.end());
        check(std::adjacent_find(bank_seeds.begin(), bank_seeds.end()) == bank_seeds.end(),
              "bank transforms: all 11 derived seeds are distinct at a common list length");

        // No two same-size banks share an exact value vector or an exact
        // rank sequence.
        bool no_vector_equality = true;
        bool no_rank_equality = true;
        for (std::size_t a = 0; a < kAllBankIds.size(); ++a) {
            for (std::size_t b = a + 1; b < kAllBankIds.size(); ++b) {
                const auto& va = banks.bank(kAllBankIds[a]).values;
                const auto& vb = banks.bank(kAllBankIds[b]).values;
                if (va.size() != vb.size()) continue;
                no_vector_equality = no_vector_equality && va != vb;
                no_rank_equality = no_rank_equality && rank_sequence(va) != rank_sequence(vb);
            }
        }
        check(no_vector_equality,
              "bank transforms: no two same-size banks have an identical value vector");
        check(no_rank_equality,
              "bank transforms: no two same-size banks have an identical rank sequence");

        // MODEL_SPECIFICATION.md section 11.1 / D4: no realized PBR value is
        // ever shared between the susceptible bank and the infector bank.
        bool no_shared_pbr = true;
        for (double s : banks.bank(BankId::susceptible_pbr).values) {
            for (double i : banks.bank(BankId::infector_pbr).values) {
                no_shared_pbr = no_shared_pbr && s != i;
            }
        }
        check(no_shared_pbr,
              "bank transforms: susceptible PBR and infector PBR share no realized value");

        // Every value is finite and inside its frozen physical domain.
        bool domains_ok = true;
        for (BankId id : kAllBankIds) {
            for (double v : banks.bank(id).values) {
                domains_ok = domains_ok && std::isfinite(v) && in_frozen_domain(id, v);
            }
        }
        check(domains_ok,
              "bank transforms: every generated value is finite and inside its frozen physical "
              "domain");

        // Determinism and per-replicate difference.
        const ParameterBanks again = build_parameter_banks(seeds, config, 11, plan, category);
        const ParameterBanks other_replicate =
            build_parameter_banks(seeds, config, 12, plan, category);
        bool deterministic = true;
        bool replicate_differs = true;
        for (BankId id : kAllBankIds) {
            deterministic = deterministic && banks.bank(id).values == again.bank(id).values;
            replicate_differs =
                replicate_differs && banks.bank(id).values != other_replicate.bank(id).values;
        }
        check(deterministic,
              "bank transforms: repeating the same (config, replicate, plan) reproduces every "
              "bank exactly");
        check(replicate_differs,
              "bank transforms: a different replicate index changes every bank");

        // Consumption ids are stable, bank-scoped, and decodable.
        bool ids_ok = true;
        std::vector<std::uint64_t> all_ids;
        for (BankId id : kAllBankIds) {
            const auto& bank = banks.bank(id);
            for (std::size_t i = 0; i < bank.values.size(); ++i) {
                const std::uint64_t cid = bank.value_id(i);
                ids_ok = ids_ok && cid == ecai_cpp2::sampling_plan::consumption_id(id, i) &&
                          ecai_cpp2::sampling_plan::consumption_id_bank(cid) == id &&
                          ecai_cpp2::sampling_plan::consumption_id_index(cid) == i;
                all_ids.push_back(cid);
            }
        }
        std::sort(all_ids.begin(), all_ids.end());
        check(ids_ok && all_ids.size() == 88,
              "bank transforms: every value carries a decodable bank/index consumption id");
        check(std::adjacent_find(all_ids.begin(), all_ids.end()) == all_ids.end(),
              "bank transforms: consumption ids are unique across all 11 banks");

        bool threw = false;
        try {
            banks.bank(BankId::gamma).value_id(8);
        } catch (const std::out_of_range&) {
            threw = true;
        }
        check(threw, "bank transforms: value_id rejects an index past the end of its bank");
    }

    // ===============================================================
    // Stage 7: cursor accounting and typed person assignments. One
    // traversal of the plan in simulation order consumes every generated
    // value exactly once: a zero-infector row consumes nothing, an active
    // row consumes one element from each susceptible bank and k elements
    // from every infector bank (MODEL_SPECIFICATION.md section 11.1, D3).
    // ===============================================================
    {
        using ecai_cpp2::sampling_plan::assign_replicate;
        using ecai_cpp2::sampling_plan::BankId;
        using ecai_cpp2::sampling_plan::build_parameter_banks;
        using ecai_cpp2::sampling_plan::InfectorParameters;
        using ecai_cpp2::sampling_plan::is_infector_bank;
        using ecai_cpp2::sampling_plan::kAllBankIds;
        using ecai_cpp2::sampling_plan::ParameterBanks;
        using ecai_cpp2::sampling_plan::ParameterCursor;
        using ecai_cpp2::sampling_plan::ReplicateAssignments;
        using ecai_cpp2::sampling_plan::SimulationAssignment;
        using ecai_cpp2::sampling_plan::SusceptibleParameters;

        const ecai_cpp2::OccupancyCategory category = make_test_category();

        // Mixed plan: {0, 2, 0, 3, 1} -> N=5, N_active=3, K_total=6.
        const SamplingPlan plan = make_plan({0, 2, 0, 3, 1});
        const ParameterBanks banks = build_parameter_banks(seeds, config, 5, plan, category);
        const ReplicateAssignments assignments = assign_replicate(banks, plan);

        // Direct coverage of the public banks() and exhausted() APIs. Use
        // literal bank contents to prove banks() exposes the source object
        // rather than a copied or reconstructed view.
        ParameterBanks accessor_fixture;
        for (BankId id : kAllBankIds) {
            const std::size_t ordinal = static_cast<std::size_t>(id);
            accessor_fixture.banks[ordinal].id = id;
            accessor_fixture.banks[ordinal].values = {
                static_cast<double>(ordinal) + 0.25,
                static_cast<double>(ordinal) + 0.75,
            };
        }
        ParameterCursor accessor_cursor(accessor_fixture);
        bool bank_reference_ok = &accessor_cursor.banks() == &accessor_fixture;
        for (BankId id : kAllBankIds) {
            const std::size_t ordinal = static_cast<std::size_t>(id);
            const auto& observed = accessor_cursor.banks().bank(id);
            bank_reference_ok = bank_reference_ok &&
                                &observed == &accessor_fixture.banks[ordinal] &&
                                &observed.values == &accessor_fixture.banks[ordinal].values &&
                                observed.values.size() == 2 &&
                                observed.values[0] == static_cast<double>(ordinal) + 0.25 &&
                                observed.values[1] == static_cast<double>(ordinal) + 0.75;
        }
        check(bank_reference_ok,
              "cursor banks(): returns the identical source banks with exact literal contents");

        ParameterCursor directly_observed(banks);
        check(!directly_observed.exhausted(),
              "cursor exhausted(): a fresh nonempty cursor is not exhausted");
        directly_observed.take_susceptible();
        directly_observed.take_infector();
        check(!directly_observed.exhausted(),
              "cursor exhausted(): partial susceptible and infector consumption remains false");
        for (int i = 0; i < 2; ++i) directly_observed.take_susceptible();
        for (int i = 0; i < 5; ++i) directly_observed.take_infector();
        check(directly_observed.exhausted(),
              "cursor exhausted(): exact full consumption of all 11 banks returns true");

        check(assignments.simulations.size() == 5,
              "cursor: one assignment row per simulation, in simulation order");
        check(assignments.n_active == 3 && assignments.k_total == 6,
              "cursor: assignments carry the plan N_active / K_total accounting");

        // Zero-infector rows are inert: no susceptible tuple, no infectors.
        bool zero_rows_ok = true;
        for (std::size_t i : {std::size_t{0}, std::size_t{2}}) {
            const SimulationAssignment& row = assignments.simulations[i];
            zero_rows_ok = zero_rows_ok && !row.active && row.infector_count == 0 &&
                            row.infectors.empty() && row.simulation_index == i;
        }
        check(zero_rows_ok,
              "cursor: a zero-infector row emits no susceptible tuple and no infector tuples");

        // Active rows carry exactly k infector tuples.
        bool active_rows_ok = assignments.simulations[1].active &&
                               assignments.simulations[1].infectors.size() == 2 &&
                               assignments.simulations[3].active &&
                               assignments.simulations[3].infectors.size() == 3 &&
                               assignments.simulations[4].active &&
                               assignments.simulations[4].infectors.size() == 1;
        check(active_rows_ok,
              "cursor: each active row carries exactly its own infector_count infector tuples");

        // Susceptible tuples take consecutive positions in every susceptible
        // bank, in active-row order; infector tuples take consecutive
        // positions in every infector bank, in (row, infector) order.
        bool positions_ok = true;
        std::size_t next_susceptible = 0;
        std::size_t next_infector = 0;
        for (const SimulationAssignment& row : assignments.simulations) {
            if (!row.active) continue;
            const SusceptibleParameters& s = row.susceptible;
            positions_ok = positions_ok && s.bank_index == next_susceptible &&
                            s.pbr_m3h == banks.bank(BankId::susceptible_pbr).values[next_susceptible] &&
                            s.gamma_per_h == banks.bank(BankId::gamma).values[next_susceptible] &&
                            s.lambda_per_h == banks.bank(BankId::lambda).values[next_susceptible];
            ++next_susceptible;
            for (const InfectorParameters& inf : row.infectors) {
                positions_ok =
                    positions_ok && inf.bank_index == next_infector &&
                    inf.pbr_m3h == banks.bank(BankId::infector_pbr).values[next_infector] &&
                    inf.c_drop == banks.bank(BankId::c_drop).values[next_infector] &&
                    inf.d_bar == banks.bank(BankId::d_bar).values[next_infector] &&
                    inf.evaporation_factor == banks.bank(BankId::E).values[next_infector] &&
                    inf.gvl_copies_per_ml == banks.bank(BankId::GVL).values[next_infector] &&
                    inf.viable_fraction == banks.bank(BankId::VF).values[next_infector] &&
                    inf.rtd == banks.bank(BankId::RTD).values[next_infector] &&
                    inf.dose_constant == banks.bank(BankId::DK).values[next_infector];
                ++next_infector;
            }
        }
        check(positions_ok,
              "cursor: each assigned value is the next unconsumed element of its own bank");
        check(next_susceptible == 3 && next_infector == 6,
              "cursor: an active row advances susceptible cursors once and infector cursors k "
              "times");

        // Direct, per-row cursor observation of the zero-infector
        // short-circuit (MODEL_SPECIFICATION.md D3): drive the same mixed
        // plan through a fresh single-use cursor one row at a time and
        // snapshot all 11 cursor positions before and after each row. A
        // zero-infector row must leave every cursor exactly where it was;
        // an active row must advance each susceptible cursor by exactly 1
        // and each infector cursor by exactly its own infector_count. The
        // checks above establish this only indirectly (through bank sizing
        // and the consumed-id set); this establishes it by observation.
        {
            ParameterCursor stepwise(banks);
            bool zero_rows_advance_nothing = true;
            bool active_rows_advance_exactly = true;
            for (std::size_t i = 0; i < plan.counts.size(); ++i) {
                std::array<std::size_t, 11> before{};
                for (std::size_t b = 0; b < kAllBankIds.size(); ++b) {
                    before[b] = stepwise.cursor(kAllBankIds[b]);
                }

                const long long count = plan.counts[i];
                if (count > 0) {
                    stepwise.take_susceptible();
                    for (long long k = 0; k < count; ++k) stepwise.take_infector();
                }

                for (std::size_t b = 0; b < kAllBankIds.size(); ++b) {
                    const BankId id = kAllBankIds[b];
                    const std::size_t moved = stepwise.cursor(id) - before[b];
                    if (count == 0) {
                        zero_rows_advance_nothing = zero_rows_advance_nothing && moved == 0;
                    } else {
                        const std::size_t expected =
                            is_infector_bank(id) ? static_cast<std::size_t>(count) : 1;
                        active_rows_advance_exactly = active_rows_advance_exactly &&
                                                       moved == expected;
                    }
                }
            }
            check(zero_rows_advance_nothing,
                  "cursor: in a mixed plan, a zero-infector row advances none of the 11 cursors "
                  "(observed directly, per row)");
            check(active_rows_advance_exactly,
                  "cursor: in the same mixed plan, each active row advances every susceptible "
                  "cursor by exactly 1 and every infector cursor by exactly its infector_count");
            stepwise.finalize();
        }

        // Every infector receives its own complete 8-value tuple; no two
        // infectors and no susceptible/infector pair share a PBR id.
        std::vector<std::uint64_t> infector_pbr_ids;
        std::vector<std::uint64_t> susceptible_pbr_ids;
        bool tuple_ids_ok = true;
        for (const SimulationAssignment& row : assignments.simulations) {
            if (!row.active) continue;
            susceptible_pbr_ids.push_back(row.susceptible.ids[0]);
            tuple_ids_ok = tuple_ids_ok && row.susceptible.ids.size() == 3;
            for (const InfectorParameters& inf : row.infectors) {
                tuple_ids_ok = tuple_ids_ok && inf.ids.size() == 8;
                infector_pbr_ids.push_back(inf.ids[0]);
            }
        }
        bool pbr_never_shared = true;
        for (std::uint64_t s : susceptible_pbr_ids) {
            for (std::uint64_t i : infector_pbr_ids) {
                pbr_never_shared = pbr_never_shared && s != i;
            }
        }
        std::sort(infector_pbr_ids.begin(), infector_pbr_ids.end());
        check(tuple_ids_ok && infector_pbr_ids.size() == 6,
              "cursor: every realized infector receives a complete, distinct 8-value tuple");
        check(std::adjacent_find(infector_pbr_ids.begin(), infector_pbr_ids.end()) ==
                  infector_pbr_ids.end(),
              "cursor: no two infectors reuse the same infector PBR id");
        check(pbr_never_shared,
              "cursor: a susceptible PBR id can never equal an infector PBR id");

        // Consumption ids: every generated value used exactly once, none
        // left unused.
        std::vector<std::uint64_t> consumed = collect_consumed_ids(assignments);
        std::vector<std::uint64_t> generated = collect_generated_ids(banks);
        std::sort(consumed.begin(), consumed.end());
        std::sort(generated.begin(), generated.end());
        check(std::adjacent_find(consumed.begin(), consumed.end()) == consumed.end(),
              "cursor: no consumption id is used twice");
        check(consumed == generated,
              "cursor: the consumed id set is exactly the generated id set (none unused)");

        // Cursors are exhausted exactly at replicate completion.
        bool cursors_exhausted = true;
        for (BankId id : kAllBankIds) {
            const std::size_t expected = is_infector_bank(id) ? plan.k_total : plan.n_active;
            cursors_exhausted = cursors_exhausted && assignments.final_cursor(id) == expected;
        }
        check(cursors_exhausted,
              "cursor: after all rows, every one of the 11 cursors equals its planned length");

        // All-zero plan: nothing is consumed and the traversal still
        // finalizes cleanly.
        const SamplingPlan zero_plan = make_plan({0, 0, 0});
        const ParameterBanks zero_banks =
            build_parameter_banks(seeds, config, 5, zero_plan, category);
        ParameterCursor fresh_zero_cursor(zero_banks);
        check(fresh_zero_cursor.exhausted(),
              "cursor exhausted(): an all-zero bank set is exhausted at construction");
        const ReplicateAssignments zero_assignments = assign_replicate(zero_banks, zero_plan);
        bool zero_ok = zero_assignments.simulations.size() == 3;
        for (const SimulationAssignment& row : zero_assignments.simulations) {
            zero_ok = zero_ok && !row.active && row.infectors.empty();
        }
        for (BankId id : kAllBankIds) {
            zero_ok = zero_ok && zero_assignments.final_cursor(id) == 0;
        }
        check(zero_ok,
              "cursor: an all-zero-infector plan consumes nothing and still finalizes exhausted");
        {
            ParameterCursor zero_cursor(zero_banks);
            assign_replicate(zero_cursor, zero_plan);
            bool threw = false;
            try {
                assign_replicate(zero_cursor, zero_plan);
            } catch (const std::runtime_error&) {
                threw = true;
            }
            check(threw,
                  "cursor: an all-zero plan still rejects a second traversal of one cursor");
        }

        // Cursor overrun: consuming past the end of a bank is detected.
        {
            ParameterCursor cursor(banks);
            for (int i = 0; i < 3; ++i) cursor.take_susceptible();
            bool threw = false;
            try {
                cursor.take_susceptible();
            } catch (const std::runtime_error&) {
                threw = true;
            }
            check(threw, "cursor: taking a susceptible tuple past the end of its banks throws");
        }
        {
            ParameterCursor cursor(banks);
            for (int i = 0; i < 6; ++i) cursor.take_infector();
            bool threw = false;
            try {
                cursor.take_infector();
            } catch (const std::runtime_error&) {
                threw = true;
            }
            check(threw, "cursor: taking an infector tuple past the end of its banks throws");
        }
        {
            // A zero-length bank set has nothing to give: there is no refill
            // path that could produce a value on demand.
            ParameterCursor cursor(zero_banks);
            bool threw_s = false;
            bool threw_i = false;
            try {
                cursor.take_susceptible();
            } catch (const std::runtime_error&) {
                threw_s = true;
            }
            try {
                cursor.take_infector();
            } catch (const std::runtime_error&) {
                threw_i = true;
            }
            check(threw_s && threw_i,
                  "cursor: an empty bank set cannot be drawn from -- there is no refill path");
        }

        // Double traversal of a single cursor is detected as an overrun.
        {
            ParameterCursor cursor(banks);
            assign_replicate(cursor, plan);
            bool threw = false;
            try {
                assign_replicate(cursor, plan);
            } catch (const std::runtime_error&) {
                threw = true;
            }
            check(threw, "cursor: a second traversal of an already-exhausted cursor throws");
        }

        // Finalizing before exhaustion is detected (cursor underflow).
        {
            ParameterCursor cursor(banks);
            cursor.take_susceptible();
            bool threw = false;
            try {
                cursor.finalize();
            } catch (const std::runtime_error&) {
                threw = true;
            }
            check(threw, "cursor: finalize() before every cursor is exhausted throws");
            check(cursor.cursor(BankId::susceptible_pbr) == 1 &&
                      cursor.cursor(BankId::infector_pbr) == 0,
                  "cursor: per-bank cursor positions are observable during traversal");
        }

        // A plan that does not match the banks it is traversed against is a
        // caller error, not a silent mis-assignment.
        {
            bool threw = false;
            try {
                assign_replicate(banks, make_plan({1, 1, 1}));
            } catch (const std::invalid_argument&) {
                threw = true;
            }
            check(threw, "cursor: a plan whose accounting does not match the banks is rejected");
        }
        {
            SamplingPlan malformed = plan;
            malformed.counts[1] = -1;
            bool threw = false;
            try {
                assign_replicate(banks, malformed);
            } catch (const std::invalid_argument&) {
                threw = true;
            }
            check(threw, "cursor: a malformed plan with a negative count is rejected");
        }
        {
            bool threw = false;
            try {
                assignments.final_cursor(static_cast<BankId>(77));
            } catch (const std::invalid_argument&) {
                threw = true;
            }
            check(threw, "cursor: final_cursor rejects an out-of-range BankId");
        }

        // Randomized synthetic count vectors: for every plan shape, one
        // traversal consumes each generated id exactly once, leaves none
        // unused, and ends with all 11 cursors exhausted.
        std::mt19937_64 rng(31337ULL);
        bool property_ok = true;
        for (int trial = 0; trial < 60 && property_ok; ++trial) {
            std::vector<long long> counts;
            const std::size_t rows = 1 + static_cast<std::size_t>(rng() % 15);
            for (std::size_t i = 0; i < rows; ++i) {
                counts.push_back(static_cast<long long>(rng() % 4));
            }
            const SamplingPlan p = make_plan(counts);
            const ParameterBanks b = build_parameter_banks(seeds, config, 9, p, category);
            const ReplicateAssignments a = assign_replicate(b, p);
            std::vector<std::uint64_t> used = collect_consumed_ids(a);
            std::vector<std::uint64_t> made = collect_generated_ids(b);
            std::sort(used.begin(), used.end());
            std::sort(made.begin(), made.end());
            property_ok = property_ok && a.simulations.size() == rows && used == made &&
                           std::adjacent_find(used.begin(), used.end()) == used.end();
            for (std::size_t i = 0; i < rows; ++i) {
                const SimulationAssignment& row = a.simulations[i];
                property_ok = property_ok && row.infector_count == counts[i] &&
                               row.active == (counts[i] > 0) &&
                               row.infectors.size() == static_cast<std::size_t>(counts[i]);
            }
            for (BankId id : kAllBankIds) {
                property_ok = property_ok &&
                               a.final_cursor(id) == (is_infector_bank(id) ? p.k_total : p.n_active);
            }
        }
        check(property_ok,
              "cursor: over 60 randomized synthetic count vectors, every id is consumed exactly "
              "once, none is unused, and all 11 cursors end exhausted");

        // Direct-indexed cursor path: every take_* must read the value at
        // the bank position it reports, advance exactly the banks of its
        // own sizing class by one, and leave the public cursor() accessor
        // (which still range-checks its BankId) in agreement with the
        // consumption ids it handed out.
        {
            ParameterCursor stepper(banks);
            bool threw = false;
            try {
                stepper.cursor(static_cast<BankId>(55));
            } catch (const std::invalid_argument&) {
                threw = true;
            }
            check(threw, "cursor: ParameterCursor::cursor rejects an out-of-range BankId");

            bool consistent = true;
            for (std::size_t draw = 0; draw < plan.n_active; ++draw) {
                const SusceptibleParameters person = stepper.take_susceptible();
                consistent = consistent && person.bank_index == draw &&
                             person.pbr_m3h == banks.bank(BankId::susceptible_pbr).values[draw] &&
                             person.gamma_per_h == banks.bank(BankId::gamma).values[draw] &&
                             person.lambda_per_h == banks.bank(BankId::lambda).values[draw] &&
                             person.ids[0] == banks.bank(BankId::susceptible_pbr).value_id(draw) &&
                             person.ids[1] == banks.bank(BankId::gamma).value_id(draw) &&
                             person.ids[2] == banks.bank(BankId::lambda).value_id(draw);
                for (BankId id : kAllBankIds) {
                    const std::size_t expected = is_infector_bank(id) ? 0 : draw + 1;
                    consistent = consistent && stepper.cursor(id) == expected;
                }
            }
            for (std::size_t draw = 0; draw < plan.k_total; ++draw) {
                const InfectorParameters person = stepper.take_infector();
                consistent = consistent && person.bank_index == draw &&
                             person.pbr_m3h == banks.bank(BankId::infector_pbr).values[draw] &&
                             person.c_drop == banks.bank(BankId::c_drop).values[draw] &&
                             person.d_bar == banks.bank(BankId::d_bar).values[draw] &&
                             person.evaporation_factor == banks.bank(BankId::E).values[draw] &&
                             person.gvl_copies_per_ml == banks.bank(BankId::GVL).values[draw] &&
                             person.viable_fraction == banks.bank(BankId::VF).values[draw] &&
                             person.rtd == banks.bank(BankId::RTD).values[draw] &&
                             person.dose_constant == banks.bank(BankId::DK).values[draw];
                const std::array<BankId, 8> infector_order = {
                    BankId::infector_pbr, BankId::c_drop, BankId::d_bar, BankId::E,
                    BankId::GVL,          BankId::VF,     BankId::RTD,   BankId::DK};
                for (std::size_t i = 0; i < infector_order.size(); ++i) {
                    consistent = consistent &&
                                 person.ids[i] == banks.bank(infector_order[i]).value_id(draw);
                }
                for (BankId id : kAllBankIds) {
                    const std::size_t expected = is_infector_bank(id) ? draw + 1 : plan.n_active;
                    consistent = consistent && stepper.cursor(id) == expected;
                }
            }
            check(consistent && stepper.exhausted(),
                  "cursor: every take_* returns the value and id at its reported bank position "
                  "and advances only its own sizing class by exactly one");

            threw = false;
            try {
                stepper.take_susceptible();
            } catch (const std::runtime_error&) {
                threw = true;
            }
            check(threw, "cursor: take_susceptible on an exhausted cursor still overruns");
            threw = false;
            try {
                stepper.take_infector();
            } catch (const std::runtime_error&) {
                threw = true;
            }
            check(threw, "cursor: take_infector on an exhausted cursor still overruns");
        }
    }

    if (g_failures > 0) {
        std::cerr << g_failures << " test(s) failed.\n";
        return 1;
    }
    std::cout << "All tests passed.\n";
    return 0;
}
