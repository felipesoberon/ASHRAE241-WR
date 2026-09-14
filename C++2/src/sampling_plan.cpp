// Stage 6: binomial pre-sampling and exact draw-count plan.
// See include/sampling_plan.h for the full contract.

#include "sampling_plan.h"

#include "distributions.h"
#include "lhs.h"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace ecai_cpp2 {
namespace sampling_plan {

namespace {

// Checked size_t multiplication: returns false (leaving *out* unspecified)
// if a * b would overflow std::size_t.
bool checked_mul(std::size_t a, std::size_t b, std::size_t& out) {
    if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a) {
        return false;
    }
    out = a * b;
    return true;
}

// Checked size_t addition: returns false (leaving *out* unspecified) if
// a + b would overflow std::size_t.
bool checked_add(std::size_t a, std::size_t b, std::size_t& out) {
    if (a > std::numeric_limits<std::size_t>::max() - b) {
        return false;
    }
    out = a + b;
    return true;
}

// Rejects an InfectorPolicy outside the declared enum before it is trusted
// (the same defensive pattern simulation.cpp's validate_infector_policy
// uses). No `default:` label, so a future third enumerator trips -Wswitch.
void require_known_infector_policy(InfectorPolicy policy, const char* context) {
    switch (policy) {
        case InfectorPolicy::community:
        case InfectorPolicy::risk_averse:
            return;
    }
    throw std::invalid_argument(std::string(context) + ": unknown InfectorPolicy");
}

}  // namespace

CountSummary summarize_counts(const std::vector<long long>& counts) {
    CountSummary summary;
    for (long long c : counts) {
        if (c < 0) {
            throw std::invalid_argument("summarize_counts: infector count must be non-negative");
        }
        if (c == 0) {
            ++summary.n_zero;
            continue;
        }
        if (static_cast<unsigned long long>(c) >
            static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max())) {
            throw std::runtime_error("summarize_counts: infector count exceeds size_t");
        }
        ++summary.n_active;
        std::size_t new_k_total = 0;
        if (!checked_add(summary.k_total, static_cast<std::size_t>(c), new_k_total)) {
            throw std::runtime_error("summarize_counts: K_total accumulation overflows size_t");
        }
        summary.k_total = new_k_total;
    }
    return summary;
}

SamplingPlan build_sampling_plan(const seed_manager::SeedManager& seeds,
                                  const seed_manager::EffectiveConfig& config,
                                  std::uint64_t replicate_index, std::size_t n, long long i0,
                                  double p, const ResourcePolicy& policy,
                                  InfectorPolicy infector_policy) {
    if (n == 0) {
        throw std::invalid_argument("build_sampling_plan: n must be nonzero");
    }
    if (i0 < 0) {
        throw std::invalid_argument("build_sampling_plan: i0 must be non-negative");
    }
    if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
        throw std::invalid_argument("build_sampling_plan: p must be finite and in [0, 1]");
    }
    require_known_infector_policy(infector_policy, "build_sampling_plan");
    // Risk-averse's own domain: I0 >= 1 and p strictly positive, checked
    // here as the authoritative direct-call guard regardless of whether a
    // CLI/scenario layer already performed the same check
    // (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md Phase 4).
    if (infector_policy == InfectorPolicy::risk_averse) {
        if (i0 < 1) {
            throw std::invalid_argument(
                "build_sampling_plan: risk-averse infector policy requires i0 >= 1");
        }
        if (!(p > 0.0)) {
            throw std::invalid_argument(
                "build_sampling_plan: risk-averse infector policy requires p > 0");
        }
    }

    // --- Every remaining check runs before any allocation. ---

    if (n > policy.max_simulations) {
        throw std::runtime_error("build_sampling_plan: n=" + std::to_string(n) +
                                  " exceeds ResourcePolicy::max_simulations=" +
                                  std::to_string(policy.max_simulations));
    }

    if (static_cast<unsigned long long>(i0) >
        static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max())) {
        throw std::runtime_error("build_sampling_plan: i0 exceeds size_t");
    }
    const std::size_t i0_unsigned = static_cast<std::size_t>(i0);

    // K_total's checked upper bound: N * I0.
    std::size_t max_possible_k_total = 0;
    if (!checked_mul(n, i0_unsigned, max_possible_k_total)) {
        throw std::runtime_error("build_sampling_plan: n * i0 overflows size_t");
    }
    if (max_possible_k_total > policy.max_total_infectors) {
        throw std::runtime_error(
            "build_sampling_plan: n * i0 = " + std::to_string(max_possible_k_total) +
            " exceeds ResourcePolicy::max_total_infectors=" +
            std::to_string(policy.max_total_infectors));
    }

    // Peak Stage-6 working set: the n-length uniform and count vectors are
    // conservatively counted as live alongside the BinomialPpfTable's peak of
    // three (i0+1)-double buffers (the transient PMF buffer, the retained
    // CDF buffer, and the retained conditional-tail buffer, see
    // distributions.h). The conditional-tail buffer is built unconditionally
    // by the table's constructor for every policy -- not only for
    // risk_averse -- so it is charged against every build_sampling_plan
    // call regardless of `infector_policy`.
    constexpr std::size_t kBytesPerSimulation = sizeof(double) + sizeof(long long);
    constexpr std::size_t kBinomialTableBuffers = 3;
    std::size_t vectors_bytes = 0;
    if (!checked_mul(n, kBytesPerSimulation, vectors_bytes)) {
        throw std::runtime_error("build_sampling_plan: vector byte-count computation overflows size_t");
    }

    if (n > std::vector<double>().max_size() || n > std::vector<long long>().max_size()) {
        throw std::runtime_error(
            "build_sampling_plan: n exceeds the maximum size supported for the "
            "uniform/count vectors");
    }

    // Mirrors distributions::BinomialPpfTable's (n+1)-length PMF/CDF buffers.
    std::size_t i0_plus_one = 0;
    if (!checked_add(i0_unsigned, static_cast<std::size_t>(1), i0_plus_one) ||
        i0_plus_one > std::vector<double>().max_size()) {
        throw std::runtime_error(
            "build_sampling_plan: i0=" + std::to_string(i0) +
            " exceeds the maximum size supported for the binomial PMF/CDF buffers");
    }
    std::size_t table_bytes = 0;
    std::size_t peak_bytes = 0;
    if (!checked_mul(i0_plus_one, kBinomialTableBuffers * sizeof(double), table_bytes) ||
        !checked_add(vectors_bytes, table_bytes, peak_bytes)) {
        throw std::runtime_error("build_sampling_plan: peak byte-count computation overflows size_t");
    }
    if (peak_bytes > policy.max_bytes) {
        throw std::runtime_error("build_sampling_plan: required peak " +
                                  std::to_string(peak_bytes) +
                                  " bytes exceeds ResourcePolicy::max_bytes=" +
                                  std::to_string(policy.max_bytes));
    }

    // --- All checks passed; allocation and sampling proceed. ---

    const std::uint64_t seed = seeds.derive_seed(config, replicate_index,
                                                  seed_manager::StreamId::infector_count,
                                                  static_cast<std::uint64_t>(n));
    const std::vector<double> uniforms = lhs::generate_lhs_1d(n, seed);

    // One table per replicate: the constant-(i0, p) PMF/CDF/conditional-tail
    // is built once here, and each u is then resolved against it -- through
    // ppf() (community) or conditional_ppf_at_least_one() (risk_averse) --
    // with exactly the same result the corresponding distributions.h entry
    // point would give for that u.
    const distributions::BinomialPpfTable table(i0, p);

    InfectorPolicyProvenance provenance;
    provenance.policy = infector_policy;
    if (infector_policy == InfectorPolicy::risk_averse) {
        provenance.zero_mass = table.zero_mass();
        provenance.zero_mass_rounded_to_zero = (provenance.zero_mass == 0.0);
        provenance.zero_mass_rounded_to_one = (provenance.zero_mass >= 1.0);
        provenance.conditional_ppf_algorithm_version = kConditionalPpfAlgorithmVersion;
        provenance.deterministic_branch_used = (p >= 1.0) || provenance.zero_mass_rounded_to_one;
    }

    std::vector<long long> counts;
    counts.reserve(n);
    if (infector_policy == InfectorPolicy::risk_averse) {
        for (double u : uniforms) {
            const long long count = table.conditional_ppf_at_least_one(u);
            if (count < 1 || count > i0) {
                throw std::runtime_error(
                    "build_sampling_plan: risk-averse conditional ppf returned " +
                    std::to_string(count) + ", outside [1, " + std::to_string(i0) + "]");
            }
            counts.push_back(count);
        }
    } else {
        for (double u : uniforms) {
            const long long count = table.ppf(u);
            if (count < 0 || count > i0) {
                throw std::runtime_error("build_sampling_plan: binomial ppf returned " +
                                          std::to_string(count) + ", outside [0, " +
                                          std::to_string(i0) + "]");
            }
            counts.push_back(count);
        }
    }

    const CountSummary summary = summarize_counts(counts);

    SamplingPlan plan;
    plan.counts = std::move(counts);
    plan.n = n;
    plan.n_zero = summary.n_zero;
    plan.n_active = summary.n_active;
    plan.k_total = summary.k_total;
    plan.susceptible_draws = summary.n_active;
    plan.infector_draws = summary.k_total;
    plan.policy = infector_policy;
    plan.infector_provenance = provenance;
    plan.lhs_generated = uniforms.size();
    plan.lhs_consumed = n;
    return plan;
}

void validate_plan_policy_context(const SamplingPlan& plan, InfectorPolicy expected_policy,
                                   const ResourcePolicy& policy) {
    require_known_infector_policy(expected_policy, "validate_plan_policy_context");
    require_known_infector_policy(plan.policy, "validate_plan_policy_context");
    if (plan.policy != expected_policy) {
        throw std::invalid_argument(
            "validate_plan_policy_context: plan policy disagrees with the expected policy");
    }

    if (plan.counts.size() != plan.n) {
        throw std::invalid_argument(
            "validate_plan_policy_context: plan.n disagrees with plan.counts.size()");
    }
    const CountSummary summary = summarize_counts(plan.counts);
    if (summary.n_zero != plan.n_zero || summary.n_active != plan.n_active ||
        summary.k_total != plan.k_total || plan.susceptible_draws != summary.n_active ||
        plan.infector_draws != summary.k_total) {
        throw std::invalid_argument(
            "validate_plan_policy_context: plan accounting disagrees with an independent "
            "summarize_counts(plan.counts)");
    }

    // Re-checks n_active <= n, k_total >= n_active, and every ResourcePolicy
    // ceiling; a plan that could never be banked is rejected here too,
    // before any bank exists.
    (void)compute_bank_budget(plan.n, plan.n_active, plan.k_total, policy);

    if (plan.policy == InfectorPolicy::risk_averse) {
        if (plan.n_zero != 0 || plan.n_active != plan.n) {
            throw std::invalid_argument(
                "validate_plan_policy_context: a risk-averse plan must have no zero-infector "
                "rows (n_active == n)");
        }
        if (plan.k_total < plan.n) {
            throw std::invalid_argument(
                "validate_plan_policy_context: a risk-averse plan must have k_total >= n");
        }
        for (long long c : plan.counts) {
            if (c < 1) {
                throw std::invalid_argument(
                    "validate_plan_policy_context: a risk-averse plan must have every count "
                    ">= 1");
            }
        }
    }
}


// ===================================================================
// Stage 7: exact-size parameter banks and cursor accounting.
// See include/sampling_plan.h for the full contract.
// ===================================================================

namespace {

// The single frozen table binding each Stage 7 bank to its stream id, its
// susceptible/infector sizing class, and its stable name. Defined once so
// bank_stream_id/is_infector_bank/bank_id_name can never disagree.
struct BankTraits {
    BankId bank;
    seed_manager::StreamId stream;
    bool infector_sized;
    const char* name;
};

constexpr std::array<BankTraits, 11> kBankTraits = {{
    {BankId::susceptible_pbr, seed_manager::StreamId::susceptible_pbr, false, "susceptible_pbr"},
    {BankId::gamma, seed_manager::StreamId::gamma, false, "gamma"},
    {BankId::lambda, seed_manager::StreamId::lambda, false, "lambda"},
    {BankId::infector_pbr, seed_manager::StreamId::infector_pbr, true, "infector_pbr"},
    {BankId::c_drop, seed_manager::StreamId::c_drop, true, "c_drop"},
    {BankId::d_bar, seed_manager::StreamId::d_bar, true, "d_bar"},
    {BankId::E, seed_manager::StreamId::E, true, "E"},
    {BankId::GVL, seed_manager::StreamId::GVL, true, "GVL"},
    {BankId::VF, seed_manager::StreamId::VF, true, "VF"},
    {BankId::RTD, seed_manager::StreamId::RTD, true, "RTD"},
    {BankId::DK, seed_manager::StreamId::DK, true, "DK"},
}};

// Bank-count constants of MODEL_SPECIFICATION.md section 11.1: 3
// susceptible-sized banks, 8 infector-sized banks.
constexpr std::size_t kSusceptibleBankCount = 3;
constexpr std::size_t kInfectorBankCount = 8;
static_assert(kSusceptibleBankCount + kInfectorBankCount == kBankTraits.size(),
              "bank-count constants must stay consistent with kBankTraits");

// kBankTraits is indexed directly by BankId ordinal below, so row i must
// describe BankId(i). Proven here at compile time rather than searched for
// at run time.
constexpr bool bank_traits_are_ordinal_indexed() {
    for (std::size_t i = 0; i < kBankTraits.size(); ++i) {
        if (static_cast<std::size_t>(kBankTraits[i].bank) != i) return false;
    }
    return true;
}
static_assert(bank_traits_are_ordinal_indexed(),
              "kBankTraits row i must describe BankId ordinal i");
static_assert(kBankTraits.size() == kAllBankIds.size(),
              "kBankTraits must cover exactly the banks in kAllBankIds");

// The ordinal of a BankId already known to be valid (a member of kAllBankIds
// or kBankTraits): no range check, so only for compile-time-fixed ids.
constexpr std::size_t ordinal_of(BankId bank) { return static_cast<std::size_t>(bank); }

// Range-checked, O(1) lookup: a BankId outside kAllBankIds is rejected with
// std::invalid_argument exactly as before.
const BankTraits& traits_for(BankId bank) {
    if (static_cast<std::uint64_t>(bank) >= kBankTraits.size()) {
        throw std::invalid_argument("sampling_plan: unknown BankId value " +
                                     std::to_string(static_cast<std::uint64_t>(bank)));
    }
    return kBankTraits[static_cast<std::size_t>(bank)];
}

}  // namespace


BankBudget compute_bank_budget(std::size_t n, std::size_t n_active, std::size_t k_total,
                                const ResourcePolicy& policy) {
    if (n_active > n) {
        throw std::invalid_argument("compute_bank_budget: N_active=" + std::to_string(n_active) +
                                     " exceeds N=" + std::to_string(n));
    }
    if (k_total < n_active) {
        throw std::invalid_argument("compute_bank_budget: K_total=" + std::to_string(k_total) +
                                     " is below N_active=" + std::to_string(n_active) +
                                     "; every active simulation has at least one infector");
    }
    if (n_active == 0 && k_total != 0) {
        throw std::invalid_argument(
            "compute_bank_budget: K_total is nonzero with N_active == 0");
    }
    if (k_total > policy.max_total_infectors) {
        throw std::runtime_error("compute_bank_budget: K_total=" + std::to_string(k_total) +
                                  " exceeds ResourcePolicy::max_total_infectors=" +
                                  std::to_string(policy.max_total_infectors));
    }

    std::size_t susceptible_elements = 0;
    std::size_t infector_elements = 0;
    std::size_t total_elements = 0;
    if (!checked_mul(n_active, kSusceptibleBankCount, susceptible_elements) ||
        !checked_mul(k_total, kInfectorBankCount, infector_elements) ||
        !checked_add(n, susceptible_elements, total_elements) ||
        !checked_add(total_elements, infector_elements, total_elements)) {
        throw std::runtime_error(
            "compute_bank_budget: N + 3*N_active + 8*K_total overflows std::size_t");
    }

    const std::size_t max_bank = std::vector<double>().max_size();
    if (n > max_bank || n_active > max_bank || k_total > max_bank) {
        throw std::runtime_error(
            "compute_bank_budget: a required bank length exceeds "
            "std::vector<double>::max_size()");
    }

    constexpr std::uint64_t kConsumptionIdCapacity = std::uint64_t{1} << 56;
    if (sizeof(std::size_t) >= sizeof(std::uint64_t) &&
        (static_cast<std::uint64_t>(n_active) > kConsumptionIdCapacity ||
         static_cast<std::uint64_t>(k_total) > kConsumptionIdCapacity)) {
        throw std::runtime_error(
            "compute_bank_budget: a required bank length exceeds the 56-bit "
            "consumption-ID index capacity");
    }

    std::size_t total_bytes = 0;
    if (!checked_mul(total_elements, sizeof(double), total_bytes)) {
        throw std::runtime_error(
            "compute_bank_budget: combined-bank byte count overflows std::size_t");
    }
    if (total_bytes > policy.max_bytes) {
        throw std::runtime_error("compute_bank_budget: combined banks require " +
                                  std::to_string(total_bytes) +
                                  " bytes, exceeding ResourcePolicy::max_bytes=" +
                                  std::to_string(policy.max_bytes));
    }

    BankBudget budget;
    budget.total_elements = total_elements;
    budget.total_bytes = total_bytes;
    return budget;
}

seed_manager::StreamId bank_stream_id(BankId bank) { return traits_for(bank).stream; }

bool is_infector_bank(BankId bank) { return traits_for(bank).infector_sized; }

std::string bank_id_name(BankId bank) { return traits_for(bank).name; }

namespace {

// consumption_id packs a bank ordinal into the high 8 bits and the index into
// the low 56, so ids are unique across banks and stable for a fixed plan.
constexpr int kConsumptionIndexBits = 56;
constexpr std::uint64_t kConsumptionIndexMask = (std::uint64_t{1} << kConsumptionIndexBits) - 1;

// Applies the frozen Stage 5 inverse transform for `bank` to one open-unit
// LHS value. Category-dependent banks take their GM/GSD from `category`.
double transform_value(BankId bank, double u, const OccupancyCategory& category) {
    switch (bank) {
        case BankId::susceptible_pbr:
        case BankId::infector_pbr:
            return distributions::lognormal_gm_gsd_ppf(u, category.pbr_gm_m3h, category.pbr_gsd);
        case BankId::gamma:
            return distributions::uniform_ppf(u, kGammaMinPerH, kGammaMaxPerH);
        case BankId::lambda:
            return distributions::lognormal_gm_gsd_ppf(u, kLambdaGmPerH, kLambdaGsd);
        case BankId::c_drop:
            return distributions::lognormal_gm_gsd_ppf(u, category.cdrop_gm_e4_per_m3,
                                                        category.cdrop_gsd);
        case BankId::d_bar:
            return distributions::lognormal_gm_gsd_ppf(u, category.d_gm_e6_m, category.d_gsd);
        case BankId::E:
            return distributions::beta_ppf(u, kEvaporationAlpha, kEvaporationBeta,
                                            kEvaporationMin, kEvaporationMax);
        case BankId::GVL:
            return distributions::log10normal_ppf(u, kGvlLog10Mu, kGvlLog10Sigma);
        case BankId::VF:
            return distributions::beta_ppf(u, kViableFractionAlpha, kViableFractionBeta,
                                            kViableFractionMin, kViableFractionMax);
        case BankId::RTD:
            return distributions::uniform_ppf(u, kRtdMin, kRtdMax);
        case BankId::DK:
            return distributions::uniform_ppf(u, kDoseConstantMin, kDoseConstantMax);
    }
    throw std::invalid_argument("transform_value: unknown BankId value " +
                                 std::to_string(static_cast<std::uint64_t>(bank)));
}

// Re-summarizes plan.counts independently and rejects any field that
// disagrees, so a hand-assembled or mutated plan cannot silently mis-size a
// bank. Runs before any bank allocation.
void validate_plan(const SamplingPlan& plan) {
    if (plan.counts.size() != plan.n) {
        throw std::invalid_argument("build_parameter_banks: plan.n=" + std::to_string(plan.n) +
                                     " disagrees with plan.counts.size()=" +
                                     std::to_string(plan.counts.size()));
    }
    const CountSummary summary = summarize_counts(plan.counts);
    if (summary.n_zero != plan.n_zero || summary.n_active != plan.n_active ||
        summary.k_total != plan.k_total || plan.susceptible_draws != summary.n_active ||
        plan.infector_draws != summary.k_total) {
        throw std::invalid_argument(
            "build_parameter_banks: plan accounting disagrees with an independent "
            "summarize_counts(plan.counts)");
    }
}

}  // namespace

std::uint64_t consumption_id(BankId bank, std::size_t index) {
    const std::uint64_t ordinal = static_cast<std::uint64_t>(traits_for(bank).bank);
    if (index > kConsumptionIndexMask) {
        throw std::invalid_argument("consumption_id: index " + std::to_string(index) +
                                     " does not fit in the low 56 bits");
    }
    return (ordinal << kConsumptionIndexBits) | static_cast<std::uint64_t>(index);
}

BankId consumption_id_bank(std::uint64_t id) {
    const BankId bank = static_cast<BankId>(id >> kConsumptionIndexBits);
    (void)traits_for(bank);  // Throws std::invalid_argument for an unknown bank.
    return bank;
}

std::size_t consumption_id_index(std::uint64_t id) {
    (void)consumption_id_bank(id);
    return static_cast<std::size_t>(id & kConsumptionIndexMask);
}

std::uint64_t ParameterBank::value_id(std::size_t index) const {
    if (index >= values.size()) {
        throw std::out_of_range("ParameterBank::value_id: index " + std::to_string(index) +
                                 " is outside bank " + bank_id_name(id) + " of length " +
                                 std::to_string(values.size()));
    }
    return consumption_id(id, index);
}

const ParameterBank& ParameterBanks::bank(BankId id) const {
    return banks[ordinal_of(traits_for(id).bank)];
}

ParameterBanks build_parameter_banks(const seed_manager::SeedManager& seeds,
                                      const seed_manager::EffectiveConfig& config,
                                      std::uint64_t replicate_index, const SamplingPlan& plan,
                                      const OccupancyCategory& category,
                                      const ResourcePolicy& policy) {
    // --- Every check below runs before the first bank allocation. ---
    validate_plan(plan);
    const BankBudget budget = compute_bank_budget(plan.n, plan.n_active, plan.k_total, policy);

    // --- All checks passed; bank generation proceeds. ---
    ParameterBanks result;
    result.n = plan.n;
    result.n_active = plan.n_active;
    result.k_total = plan.k_total;
    result.budget = budget;

    for (const BankTraits& traits : kBankTraits) {
        const std::size_t length = traits.infector_sized ? plan.k_total : plan.n_active;
        ParameterBank bank;
        bank.id = traits.bank;
        bank.stream = traits.stream;
        if (length == 0) {
            // No seed derived and no LHS list generated for an empty bank.
            result.banks[static_cast<std::size_t>(traits.bank)] = std::move(bank);
            continue;
        }
        bank.seed = seeds.derive_seed(config, replicate_index, traits.stream,
                                       static_cast<std::uint64_t>(length));
        bank.generated = true;
        const std::vector<double> uniforms = lhs::generate_lhs_1d(length, bank.seed);
        bank.values.reserve(length);
        for (std::size_t i = 0; i < uniforms.size(); ++i) {
            const double value = transform_value(traits.bank, uniforms[i], category);
            if (!std::isfinite(value)) {
                throw std::runtime_error("build_parameter_banks: bank " +
                                          std::string(traits.name) + " index " +
                                          std::to_string(i) + " transformed to a nonfinite value");
            }
            bank.values.push_back(value);
        }
        result.banks[static_cast<std::size_t>(traits.bank)] = std::move(bank);
    }
    return result;
}

namespace {

// The 3 susceptible banks and the 8 infector banks, in the fixed order used
// by SusceptibleParameters::ids and InfectorParameters::ids.
constexpr std::array<BankId, 3> kSusceptibleBanks = {BankId::susceptible_pbr, BankId::gamma,
                                                      BankId::lambda};
constexpr std::array<BankId, 8> kInfectorBanks = {
    BankId::infector_pbr, BankId::c_drop, BankId::d_bar, BankId::E,
    BankId::GVL,          BankId::VF,     BankId::RTD,   BankId::DK};

// Every id the cursor indexes banks/cursors_ with directly is a member of
// kBankTraits with the sizing class its list name claims.
template <std::size_t N>
constexpr bool bank_list_is_valid(const std::array<BankId, N>& ids, bool infector_sized) {
    for (BankId id : ids) {
        const std::size_t ordinal = static_cast<std::size_t>(id);
        if (ordinal >= kBankTraits.size()) return false;
        if (kBankTraits[ordinal].infector_sized != infector_sized) return false;
    }
    return true;
}
static_assert(bank_list_is_valid(kSusceptibleBanks, false),
              "kSusceptibleBanks must list exactly susceptible-sized kBankTraits rows");
static_assert(bank_list_is_valid(kInfectorBanks, true),
              "kInfectorBanks must list exactly infector-sized kBankTraits rows");
static_assert(kSusceptibleBanks.size() == kSusceptibleBankCount &&
                  kInfectorBanks.size() == kInfectorBankCount,
              "cursor bank lists must match the section 11.1 bank counts");

}  // namespace

std::size_t ReplicateAssignments::final_cursor(BankId id) const {
    return final_cursors[ordinal_of(traits_for(id).bank)];
}

std::size_t ParameterCursor::cursor(BankId id) const {
    return cursors_[ordinal_of(traits_for(id).bank)];
}

// The cursor paths below index banks_->banks and cursors_ directly by the
// ordinals of the compile-time-verified kSusceptibleBanks / kInfectorBanks
// lists instead of re-validating each id through bank()/cursor() on every
// draw. The per-bank exhaustion check, the per-value bounds-checked
// value_id, the bank ids and the forward-only cursor behavior are unchanged.
SusceptibleParameters ParameterCursor::take_susceptible() {
    const std::array<ParameterBank, 11>& banks = banks_->banks;
    for (BankId id : kSusceptibleBanks) {
        const std::size_t ordinal = ordinal_of(id);
        if (cursors_[ordinal] >= banks[ordinal].values.size()) {
            throw std::runtime_error("ParameterCursor::take_susceptible: bank " +
                                      bank_id_name(id) + " is exhausted at index " +
                                      std::to_string(cursors_[ordinal]) +
                                      "; there is no refill path");
        }
    }
    const auto value_at = [&](BankId id) {
        return banks[ordinal_of(id)].values[cursors_[ordinal_of(id)]];
    };
    SusceptibleParameters person;
    person.bank_index = cursors_[ordinal_of(BankId::susceptible_pbr)];
    person.pbr_m3h = value_at(BankId::susceptible_pbr);
    person.gamma_per_h = value_at(BankId::gamma);
    person.lambda_per_h = value_at(BankId::lambda);
    for (std::size_t i = 0; i < kSusceptibleBanks.size(); ++i) {
        const std::size_t ordinal = ordinal_of(kSusceptibleBanks[i]);
        person.ids[i] = banks[ordinal].value_id(cursors_[ordinal]);
        ++cursors_[ordinal];
    }
    return person;
}

InfectorParameters ParameterCursor::take_infector() {
    const std::array<ParameterBank, 11>& banks = banks_->banks;
    for (BankId id : kInfectorBanks) {
        const std::size_t ordinal = ordinal_of(id);
        if (cursors_[ordinal] >= banks[ordinal].values.size()) {
            throw std::runtime_error("ParameterCursor::take_infector: bank " + bank_id_name(id) +
                                      " is exhausted at index " +
                                      std::to_string(cursors_[ordinal]) +
                                      "; there is no refill path");
        }
    }
    const auto value_at = [&](BankId id) {
        return banks[ordinal_of(id)].values[cursors_[ordinal_of(id)]];
    };
    InfectorParameters person;
    person.bank_index = cursors_[ordinal_of(BankId::infector_pbr)];
    person.pbr_m3h = value_at(BankId::infector_pbr);
    person.c_drop = value_at(BankId::c_drop);
    person.d_bar = value_at(BankId::d_bar);
    person.evaporation_factor = value_at(BankId::E);
    person.gvl_copies_per_ml = value_at(BankId::GVL);
    person.viable_fraction = value_at(BankId::VF);
    person.rtd = value_at(BankId::RTD);
    person.dose_constant = value_at(BankId::DK);
    for (std::size_t i = 0; i < kInfectorBanks.size(); ++i) {
        const std::size_t ordinal = ordinal_of(kInfectorBanks[i]);
        person.ids[i] = banks[ordinal].value_id(cursors_[ordinal]);
        ++cursors_[ordinal];
    }
    return person;
}

bool ParameterCursor::exhausted() const {
    for (std::size_t ordinal = 0; ordinal < kBankTraits.size(); ++ordinal) {
        if (cursors_[ordinal] != banks_->banks[ordinal].values.size()) {
            return false;
        }
    }
    return true;
}

void ParameterCursor::finalize() const {
    for (std::size_t ordinal = 0; ordinal < kBankTraits.size(); ++ordinal) {
        const BankTraits& traits = kBankTraits[ordinal];
        const std::size_t position = cursors_[ordinal];
        const std::size_t length = banks_->banks[ordinal].values.size();
        if (position != length) {
            throw std::runtime_error("ParameterCursor::finalize: bank " +
                                      std::string(traits.name) + " consumed " +
                                      std::to_string(position) + " of " +
                                      std::to_string(length) + " planned values");
        }
    }
}

ReplicateAssignments assign_replicate(ParameterCursor& cursor, const SamplingPlan& plan) {
    validate_plan(plan);
    const ParameterBanks& banks = cursor.banks();
    if (banks.n != plan.n || banks.n_active != plan.n_active || banks.k_total != plan.k_total) {
        throw std::invalid_argument(
            "assign_replicate: plan accounting does not match the banks it is traversed against");
    }
    if (cursor.traversal_started_) {
        throw std::runtime_error("assign_replicate: this cursor has already been traversed");
    }
    cursor.traversal_started_ = true;

    ReplicateAssignments assignments;
    assignments.n_active = plan.n_active;
    assignments.k_total = plan.k_total;
    assignments.simulations.reserve(plan.n);
    for (std::size_t i = 0; i < plan.counts.size(); ++i) {
        SimulationAssignment row;
        row.simulation_index = i;
        row.infector_count = plan.counts[i];
        if (plan.counts[i] == 0) {
            // D3: a zero-infector row advances no cursor at all.
            assignments.simulations.push_back(std::move(row));
            continue;
        }
        row.active = true;
        row.susceptible = cursor.take_susceptible();
        row.infectors.reserve(static_cast<std::size_t>(plan.counts[i]));
        for (long long k = 0; k < plan.counts[i]; ++k) {
            row.infectors.push_back(cursor.take_infector());
        }
        assignments.simulations.push_back(std::move(row));
    }

    cursor.finalize();
    // cursors_ and final_cursors share the BankId-ordinal layout, so the
    // finalized positions are copied across in one step.
    assignments.final_cursors = cursor.cursors_;
    return assignments;
}

ReplicateAssignments assign_replicate(const ParameterBanks& banks, const SamplingPlan& plan) {
    ParameterCursor cursor(banks);
    return assign_replicate(cursor, plan);
}
}  // namespace sampling_plan
}  // namespace ecai_cpp2
