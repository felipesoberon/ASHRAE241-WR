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
                                  double p, const ResourcePolicy& policy) {
    if (n == 0) {
        throw std::invalid_argument("build_sampling_plan: n must be nonzero");
    }
    if (i0 < 0) {
        throw std::invalid_argument("build_sampling_plan: i0 must be non-negative");
    }
    if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
        throw std::invalid_argument("build_sampling_plan: p must be finite and in [0, 1]");
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

    // Peak Stage-6 working set: the n-length uniform and count vectors remain
    // live while binomial_ppf temporarily owns one (i0+1)-double PMF buffer.
    constexpr std::size_t kBytesPerSimulation = sizeof(double) + sizeof(long long);
    std::size_t vectors_bytes = 0;
    if (!checked_mul(n, kBytesPerSimulation, vectors_bytes)) {
        throw std::runtime_error("build_sampling_plan: vector byte-count computation overflows size_t");
    }

    if (n > std::vector<double>().max_size() || n > std::vector<long long>().max_size()) {
        throw std::runtime_error(
            "build_sampling_plan: n exceeds the maximum size supported for the "
            "uniform/count vectors");
    }

    // Mirrors distributions::binomial_ppf's internal (n+1)-length PMF buffer.
    std::size_t i0_plus_one = 0;
    if (!checked_add(i0_unsigned, static_cast<std::size_t>(1), i0_plus_one) ||
        i0_plus_one > std::vector<double>().max_size()) {
        throw std::runtime_error(
            "build_sampling_plan: i0=" + std::to_string(i0) +
            " exceeds the maximum size supported for the binomial_ppf PMF buffer");
    }
    std::size_t pmf_bytes = 0;
    std::size_t peak_bytes = 0;
    if (!checked_mul(i0_plus_one, sizeof(double), pmf_bytes) ||
        !checked_add(vectors_bytes, pmf_bytes, peak_bytes)) {
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

    std::vector<long long> counts;
    counts.reserve(n);
    for (double u : uniforms) {
        const long long count = distributions::binomial_ppf(u, i0, p);
        if (count < 0 || count > i0) {
            throw std::runtime_error("build_sampling_plan: binomial_ppf returned " +
                                      std::to_string(count) + ", outside [0, " +
                                      std::to_string(i0) + "]");
        }
        counts.push_back(count);
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
    return plan;
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

const BankTraits& traits_for(BankId bank) {
    for (const BankTraits& t : kBankTraits) {
        if (t.bank == bank) return t;
    }
    throw std::invalid_argument("sampling_plan: unknown BankId value " +
                                 std::to_string(static_cast<std::uint64_t>(bank)));
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
    const BankTraits& traits = traits_for(id);
    return banks[static_cast<std::size_t>(traits.bank)];
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

}  // namespace

std::size_t ReplicateAssignments::final_cursor(BankId id) const {
    return final_cursors[static_cast<std::size_t>(traits_for(id).bank)];
}

std::size_t ParameterCursor::cursor(BankId id) const {
    return cursors_[static_cast<std::size_t>(traits_for(id).bank)];
}

SusceptibleParameters ParameterCursor::take_susceptible() {
    for (BankId id : kSusceptibleBanks) {
        const std::size_t ordinal = static_cast<std::size_t>(id);
        if (cursors_[ordinal] >= banks_->bank(id).values.size()) {
            throw std::runtime_error("ParameterCursor::take_susceptible: bank " +
                                      bank_id_name(id) + " is exhausted at index " +
                                      std::to_string(cursors_[ordinal]) +
                                      "; there is no refill path");
        }
    }
    const std::size_t index = cursors_[static_cast<std::size_t>(BankId::susceptible_pbr)];
    SusceptibleParameters person;
    person.bank_index = index;
    person.pbr_m3h = banks_->bank(BankId::susceptible_pbr).values[index];
    person.gamma_per_h = banks_->bank(BankId::gamma).values[cursor(BankId::gamma)];
    person.lambda_per_h = banks_->bank(BankId::lambda).values[cursor(BankId::lambda)];
    for (std::size_t i = 0; i < kSusceptibleBanks.size(); ++i) {
        const BankId id = kSusceptibleBanks[i];
        person.ids[i] = banks_->bank(id).value_id(cursor(id));
        ++cursors_[static_cast<std::size_t>(id)];
    }
    return person;
}

InfectorParameters ParameterCursor::take_infector() {
    for (BankId id : kInfectorBanks) {
        const std::size_t ordinal = static_cast<std::size_t>(id);
        if (cursors_[ordinal] >= banks_->bank(id).values.size()) {
            throw std::runtime_error("ParameterCursor::take_infector: bank " + bank_id_name(id) +
                                      " is exhausted at index " +
                                      std::to_string(cursors_[ordinal]) +
                                      "; there is no refill path");
        }
    }
    InfectorParameters person;
    person.bank_index = cursor(BankId::infector_pbr);
    person.pbr_m3h = banks_->bank(BankId::infector_pbr).values[cursor(BankId::infector_pbr)];
    person.c_drop = banks_->bank(BankId::c_drop).values[cursor(BankId::c_drop)];
    person.d_bar = banks_->bank(BankId::d_bar).values[cursor(BankId::d_bar)];
    person.evaporation_factor = banks_->bank(BankId::E).values[cursor(BankId::E)];
    person.gvl_copies_per_ml = banks_->bank(BankId::GVL).values[cursor(BankId::GVL)];
    person.viable_fraction = banks_->bank(BankId::VF).values[cursor(BankId::VF)];
    person.rtd = banks_->bank(BankId::RTD).values[cursor(BankId::RTD)];
    person.dose_constant = banks_->bank(BankId::DK).values[cursor(BankId::DK)];
    for (std::size_t i = 0; i < kInfectorBanks.size(); ++i) {
        const BankId id = kInfectorBanks[i];
        person.ids[i] = banks_->bank(id).value_id(cursor(id));
        ++cursors_[static_cast<std::size_t>(id)];
    }
    return person;
}

bool ParameterCursor::exhausted() const {
    for (const BankTraits& traits : kBankTraits) {
        if (cursors_[static_cast<std::size_t>(traits.bank)] !=
            banks_->bank(traits.bank).values.size()) {
            return false;
        }
    }
    return true;
}

void ParameterCursor::finalize() const {
    for (const BankTraits& traits : kBankTraits) {
        const std::size_t position = cursors_[static_cast<std::size_t>(traits.bank)];
        const std::size_t length = banks_->bank(traits.bank).values.size();
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
    for (const BankTraits& traits : kBankTraits) {
        assignments.final_cursors[static_cast<std::size_t>(traits.bank)] = cursor.cursor(traits.bank);
    }
    return assignments;
}

ReplicateAssignments assign_replicate(const ParameterBanks& banks, const SamplingPlan& plan) {
    ParameterCursor cursor(banks);
    return assign_replicate(cursor, plan);
}
}  // namespace sampling_plan
}  // namespace ecai_cpp2
