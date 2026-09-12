// Stage 6: binomial pre-sampling and exact draw-count plan.
//
// Source of truth: docs/MODEL_SPECIFICATION.md section 11.4 (resource /
// memory budget policy) and IMPLEMENTATION_PLAN.md Section 7, Stage 6.
//
// Objective: generate all `ninf` (infector count) values for one
// (configuration, replicate) pair first, from a single exact-size 1D LHS
// list transformed through the binomial inverse CDF (distributions.h), then
// compute the exact downstream draw-list sizes (one susceptible-parameter
// draw per active simulation, one infector-parameter draw per infector
// across every active simulation) without allocating any Stage 7 parameter
// bank. Zero-infector simulations are accounted for but imply no other
// per-simulation draw.

#ifndef ECAI_CPP2_SAMPLING_PLAN_H
#define ECAI_CPP2_SAMPLING_PLAN_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "category_data.h"
#include "seed_manager.h"

namespace ecai_cpp2 {
namespace sampling_plan {

// Pre-allocation resource limits for a single build_sampling_plan call.
// Every field is an explicit, caller-overridable ceiling checked before any
// allocation. max_bytes is a Stage-6 peak-working-set safeguard. Stage 7 must
// separately enforce MODEL_SPECIFICATION.md section 11.4's combined-bank
// formula N + 3*N_active + 8*K_total before creating parameter banks. These
// limits are engineering defaults, not Jones et al. values.
struct ResourcePolicy {
    // Ceiling on N, the requested LHS list length / simulation count.
    std::size_t max_simulations = 536870912ULL;  // 2^29
    // Ceiling on K_total (the sum of every simulation's infector count) and
    // on its checked upper bound N * I0.
    std::size_t max_total_infectors = 536870912ULL;  // 2^29
    // Ceiling on peak bytes: N uniforms + N counts + the two simultaneously-
    // live (I0+1)-double binomial buffers (PMF and CDF) of the per-replicate
    // distributions::BinomialPpfTable.
    std::size_t max_bytes = 4ULL * 1024 * 1024 * 1024;  // 4 GiB
};

// Pure accounting summary of a fixed infector-count vector: no seeding, no
// sampling, no allocation beyond its own return value. Exposed
// independently of build_sampling_plan so hand-reconciled fixtures can
// exercise the N_zero/N_active/K_total accounting in isolation from the
// RNG/LHS/binomial pipeline.
struct CountSummary {
    std::size_t n_zero = 0;
    std::size_t n_active = 0;
    std::size_t k_total = 0;
};

// Summarizes `counts` into N_zero (simulations with count == 0), N_active
// (simulations with count > 0), and K_total (checked sum of every positive
// count). Throws std::invalid_argument if any count is negative. Throws
// std::runtime_error if the checked K_total accumulation would overflow
// std::size_t.
CountSummary summarize_counts(const std::vector<long long>& counts);

// The full Stage 6 sampling plan for one (configuration, replicate) pair:
// exactly `n` infector counts, in simulation order, plus their accounting.
struct SamplingPlan {
    // counts[i] is simulation i's infector count (ninf), in the same order
    // as the underlying LHS list; simulation order is preserved.
    std::vector<long long> counts;
    std::size_t n = 0;
    std::size_t n_zero = 0;
    std::size_t n_active = 0;
    std::size_t k_total = 0;
    std::size_t susceptible_draws = 0;  // == n_active
    std::size_t infector_draws = 0;     // == k_total
};

// Generates exactly `n` open-unit 1D LHS values (seeded from
// seeds.derive_seed(config, replicate_index, StreamId::infector_count, n)),
// transforms each through one per-call distributions::BinomialPpfTable(i0,
// p) -- whose ppf(u) equals distributions::binomial_ppf(u, i0, p) exactly,
// with the constant (i0, p) PMF/CDF built once rather than once per u -- to
// obtain one infector count per simulation, and returns the resulting plan.
//
// Every check below runs, in order, before any allocation:
//   - n == 0: std::invalid_argument.
//   - i0 < 0: std::invalid_argument.
//   - p nonfinite or outside [0, 1]: std::invalid_argument.
//   - n > policy.max_simulations: std::runtime_error.
//   - checked n * static_cast<size_t>(i0) (K_total's maximum possible
//     value) overflowing std::size_t, or exceeding
//     policy.max_total_infectors: std::runtime_error.
//   - checked peak byte cost of the n-length uniform/count vectors plus the
//     two `(i0+1)`-double binomial buffers (PMF and CDF), overflowing
//     std::size_t or exceeding policy.max_bytes, or n exceeding
//     std::vector<double>::max_size() / std::vector<long long>::max_size():
//     std::runtime_error.
//   - checked i0 + 1 (distributions::BinomialPpfTable's PMF/CDF buffer
//     size) overflowing std::size_t or exceeding
//     std::vector<double>::max_size(): std::runtime_error.
//
// Zero-infector simulations (counts[i] == 0) are accounted for but imply no
// other per-simulation draw; susceptible_draws counts only active
// simulations and infector_draws is the checked sum of every positive
// count (see summarize_counts).
SamplingPlan build_sampling_plan(const seed_manager::SeedManager& seeds,
                                  const seed_manager::EffectiveConfig& config,
                                  std::uint64_t replicate_index, std::size_t n, long long i0,
                                  double p, const ResourcePolicy& policy = ResourcePolicy{});


// ===================================================================
// Stage 7: exact-size parameter banks and cursor accounting.
//
// Source of truth: docs/MODEL_SPECIFICATION.md section 11.1
// (person-specific sampling: 3 susceptible banks sized N_active, 8 infector
// banks sized K_total, nothing shared between a susceptible and an infector
// or between two infectors), section 5 / section 7 (the frozen Table 1 and
// Table 3 distribution families and parameters), section 11.3 (one fixed
// StreamId per named parameter), and section 11.4 (the combined-bank
// N + 3*N_active + 8*K_total budget, checked before any allocation).
// IMPLEMENTATION_PLAN.md Section 7, Stage 7.
// ===================================================================


// The exact combined size of every Stage 7 parameter bank in one replicate,
// per MODEL_SPECIFICATION.md section 11.4.
struct BankBudget {
    // N (the Stage 6 binomial list) + 3*N_active (susceptible banks) +
    // 8*K_total (infector banks).
    std::size_t total_elements = 0;
    // sizeof(double) * total_elements.
    std::size_t total_bytes = 0;
};

// Computes the combined-bank budget with fully checked arithmetic and
// enforces every pre-allocation ceiling, before any bank is created. Runs,
// in order:
//   - n_active > n, or (k_total < n_active), or (n_active == 0 &&
//     k_total != 0): std::invalid_argument (an inconsistent accounting
//     triple; every active simulation has at least one infector).
//   - k_total > policy.max_total_infectors: std::runtime_error.
//   - n, n_active, or k_total exceeding std::vector<double>::max_size():
//     std::runtime_error.
//   - checked 3*n_active, 8*k_total, and their sum with n overflowing
//     std::size_t: std::runtime_error.
//   - checked sizeof(double) * total_elements overflowing std::size_t, or
//     exceeding policy.max_bytes (4 GiB default): std::runtime_error.
// This check is separate from, and additional to, build_sampling_plan's
// Stage 6 peak-working-set check.
BankBudget compute_bank_budget(std::size_t n, std::size_t n_active, std::size_t k_total,
                                const ResourcePolicy& policy = ResourcePolicy{});

// One identifier per named Stage 7 physical parameter bank. Deterministic
// constants (R, D) and shared per-room values (VOL, mask efficiency) are
// not banks: they are not sampled and consume no LHS values.
enum class BankId : std::uint64_t {
    // Susceptible-sized banks (length N_active).
    susceptible_pbr = 0,
    gamma = 1,
    lambda = 2,
    // Infector-sized banks (length K_total).
    infector_pbr = 3,
    c_drop = 4,
    d_bar = 5,
    E = 6,
    GVL = 7,
    VF = 8,
    RTD = 9,
    DK = 10,
};

// Exhaustive, ordered list of every Stage 7 bank: the 3 susceptible-sized
// banks first, then the 8 infector-sized banks. Order is frozen and is the
// bank-construction order used by build_parameter_banks.
inline constexpr std::array<BankId, 11> kAllBankIds = {
    BankId::susceptible_pbr, BankId::gamma, BankId::lambda, BankId::infector_pbr,
    BankId::c_drop,          BankId::d_bar, BankId::E,      BankId::GVL,
    BankId::VF,              BankId::RTD,   BankId::DK,
};
static_assert(kAllBankIds.size() == static_cast<std::size_t>(BankId::DK) + 1,
              "kAllBankIds must remain exhaustive for contiguous BankId values");

// The frozen Stage 4 stream id owned by `bank`. Each bank has its own
// stream; none reuses StreamId::infector_count (Stage 6's binomial stream).
// Throws std::invalid_argument for a BankId not in kAllBankIds.
seed_manager::StreamId bank_stream_id(BankId bank);

// True for the 8 infector-sized banks (length K_total), false for the 3
// susceptible-sized banks (length N_active). Throws std::invalid_argument
// for a BankId not in kAllBankIds.
bool is_infector_bank(BankId bank);

// Human-readable name for a bank (its enumerator spelling). Throws
// std::invalid_argument for a BankId not in kAllBankIds.
std::string bank_id_name(BankId bank);

// Frozen Jones et al. Table 1 distribution parameters (MODEL_SPECIFICATION.md
// section 5) for the category-independent banks. The category-dependent banks
// (susceptible/infector PBR, C_drop, d_bar) take their GM/GSD from the
// OccupancyCategory record (Table 3, MODEL_SPECIFICATION.md section 7).
inline constexpr double kGammaMinPerH = 0.42;  // Surface deposition rate, uniform.
inline constexpr double kGammaMaxPerH = 0.61;
inline constexpr double kLambdaGmPerH = 0.52;  // Biological decay rate, lognormal GM/GSD.
inline constexpr double kLambdaGsd = 1.9;
inline constexpr double kEvaporationAlpha = 5.0;  // Evaporation factor E, scaled beta.
inline constexpr double kEvaporationBeta = 2.0;
inline constexpr double kEvaporationMin = 2.0;
inline constexpr double kEvaporationMax = 5.0;
inline constexpr double kGvlLog10Mu = 7.0;  // Genomic viral load, base-10 lognormal.
inline constexpr double kGvlLog10Sigma = 1.4;
inline constexpr double kViableFractionAlpha = 5.0;  // Viable fraction VF, scaled beta.
inline constexpr double kViableFractionBeta = 2.0;
inline constexpr double kViableFractionMin = 1e-4;
inline constexpr double kViableFractionMax = 1e-2;
inline constexpr double kRtdMin = 0.43;  // Respiratory tract deposition, uniform.
inline constexpr double kRtdMax = 0.65;
inline constexpr double kDoseConstantMin = 5.0;  // Dose constant DK, uniform.
inline constexpr double kDoseConstantMax = 15.0;

// Stable, globally unique identifier for one generated value: its bank and
// its index within that bank. Encoded as (bank ordinal << 56) | index, so two
// different banks can never produce the same id and an id is stable across
// runs for a fixed plan. Throws std::invalid_argument for a BankId not in
// kAllBankIds, or an index that does not fit in the low 56 bits.
std::uint64_t consumption_id(BankId bank, std::size_t index);

// Inverse of consumption_id: recovers the bank ordinal and index. Throws
// std::invalid_argument if `id` does not decode to a bank in kAllBankIds.
BankId consumption_id_bank(std::uint64_t id);
std::size_t consumption_id_index(std::uint64_t id);

// One exact-size bank of transformed physical values for a single named
// parameter, plus the provenance needed to reproduce it.
//
// Units are exactly those of the frozen distribution parameters the bank was
// drawn from, with no unit conversion applied: PBR in m^3/h, gamma and lambda
// in 1/h, C_drop in the Table 3 printed x10^4 m^-3 units, d_bar in the Table 3
// printed x10^-6 m units, GVL in RNA copies/mL (Table 1 log10 units
// exponentiated), and E/VF/RTD/DK dimensionless. Stage 8 applies each unit
// conversion exactly once; Stage 7 performs no model arithmetic.
struct ParameterBank {
    BankId id = BankId::susceptible_pbr;
    seed_manager::StreamId stream = seed_manager::StreamId::susceptible_pbr;
    // Seed derived from this bank's own stream id and its own exact list
    // length. Zero and meaningless when `generated` is false.
    std::uint64_t seed = 0;
    // False for a zero-length bank: no seed is derived and no LHS list is
    // generated for it (there is nothing to draw).
    bool generated = false;
    std::vector<double> values;

    // consumption_id(id, index), with index bounds-checked against values.
    // Throws std::out_of_range if index >= values.size().
    std::uint64_t value_id(std::size_t index) const;
};

// The complete set of 11 Stage 7 banks for one (configuration, replicate)
// pair, together with the plan accounting they were sized from.
struct ParameterBanks {
    // Indexed by BankId ordinal, in kAllBankIds order.
    std::array<ParameterBank, 11> banks{};
    std::size_t n = 0;
    std::size_t n_active = 0;
    std::size_t k_total = 0;
    BankBudget budget;

    // Throws std::invalid_argument for a BankId not in kAllBankIds.
    const ParameterBank& bank(BankId id) const;
};

// Generates every Stage 7 parameter bank for one completed Stage 6 plan.
//
// Sizing (MODEL_SPECIFICATION.md section 11.1): the 3 susceptible banks are
// each exactly plan.n_active long; the 8 infector banks are each exactly
// plan.k_total long (aggregate infector demand across the whole replicate).
//
// Per bank, in kAllBankIds order and only when its length is nonzero:
//   1. seed = seeds.derive_seed(config, replicate_index, bank_stream_id(bank),
//      that bank's exact list length);
//   2. one lhs::generate_lhs_1d(length, seed) call -- exactly one LHS list per
//      bank, never appended to, never regenerated, with no refill path;
//   3. the frozen Stage 5 inverse transform for that parameter, applied
//      element-wise in list order:
//        susceptible_pbr / infector_pbr : lognormal_gm_gsd_ppf(u, category
//                                          PBR GM, category PBR GSD)
//        gamma                          : uniform_ppf(u, 0.42, 0.61)
//        lambda                         : lognormal_gm_gsd_ppf(u, 0.52, 1.9)
//        c_drop                         : lognormal_gm_gsd_ppf(u, category
//                                          C_drop GM, category C_drop GSD)
//        d_bar                          : lognormal_gm_gsd_ppf(u, category
//                                          d_bar GM, category d_bar GSD)
//        E                              : beta_ppf(u, 5, 2, 2, 5)
//        GVL                            : log10normal_ppf(u, 7.0, 1.4)
//        VF                             : beta_ppf(u, 5, 2, 1e-4, 1e-2)
//        RTD                            : uniform_ppf(u, 0.43, 0.65)
//        DK                             : uniform_ppf(u, 5, 15)
// A zero-length bank is left empty and ungenerated: no seed is derived and
// generate_lhs_1d is never called for it.
//
// Susceptible and infector PBR are drawn from two separate banks with two
// separate stream ids, so no realized PBR is ever shared between a
// susceptible and an infector, or between two infectors (section 11.1, D4).
//
// Every check runs before any bank is allocated:
//   - plan.counts.size() != plan.n, any negative count, or plan accounting
//     (n_zero, n_active, k_total, susceptible_draws, infector_draws)
//     disagreeing with a fresh summarize_counts(plan.counts):
//     std::invalid_argument.
//   - compute_bank_budget(plan.n, plan.n_active, plan.k_total, policy):
//     std::runtime_error on any overflow or ceiling overrun, leaving no
//     partial bank output.
// Category GM/GSD values outside their transform domains are rejected with
// std::invalid_argument by the Stage 5 transform contract. A transform that
// would produce a nonfinite value throws std::runtime_error naming the
// offending bank and index.
ParameterBanks build_parameter_banks(const seed_manager::SeedManager& seeds,
                                      const seed_manager::EffectiveConfig& config,
                                      std::uint64_t replicate_index, const SamplingPlan& plan,
                                      const OccupancyCategory& category,
                                      const ResourcePolicy& policy = ResourcePolicy{});

// One active simulation's susceptible person, as a typed tuple of already-
// sampled values. Carries no model equations: Stage 8/9 consume these.
// Units follow ParameterBank: PBR in m^3/h, gamma and lambda in 1/h.
struct SusceptibleParameters {
    double pbr_m3h = 0.0;
    double gamma_per_h = 0.0;
    double lambda_per_h = 0.0;
    // Position consumed in every susceptible bank (all three susceptible
    // cursors advance together, so one index describes all of them).
    std::size_t bank_index = 0;
    // Stable consumption ids, in {susceptible_pbr, gamma, lambda} order.
    std::array<std::uint64_t, 3> ids{};
};

// One realized infector's complete parameter set. Every infector gets its
// own eight values; nothing is shared with the susceptible or with another
// infector (MODEL_SPECIFICATION.md section 11.1, D4). Units follow
// ParameterBank: PBR in m^3/h, C_drop in printed x10^4 m^-3 units, d_bar in
// printed x10^-6 m units, GVL in RNA copies/mL, E/VF/RTD/DK dimensionless.
struct InfectorParameters {
    double pbr_m3h = 0.0;
    double c_drop = 0.0;
    double d_bar = 0.0;
    double evaporation_factor = 0.0;
    double gvl_copies_per_ml = 0.0;
    double viable_fraction = 0.0;
    double rtd = 0.0;
    double dose_constant = 0.0;
    // Position consumed in every infector bank (all eight infector cursors
    // advance together).
    std::size_t bank_index = 0;
    // Stable consumption ids, in {infector_pbr, c_drop, d_bar, E, GVL, VF,
    // RTD, DK} order.
    std::array<std::uint64_t, 8> ids{};
};

// One simulation row of the replicate, in simulation order.
struct SimulationAssignment {
    std::size_t simulation_index = 0;
    long long infector_count = 0;
    // False for a zero-infector row: `susceptible` is then default-valued
    // and meaningless, and `infectors` is empty (D3: a zero-infector
    // simulation consumes nothing from any bank).
    bool active = false;
    SusceptibleParameters susceptible;
    std::vector<InfectorParameters> infectors;  // size == infector_count
};

// A completed, fully accounted traversal of one replicate.
struct ReplicateAssignments {
    std::vector<SimulationAssignment> simulations;  // one row per simulation
    std::size_t n_active = 0;
    std::size_t k_total = 0;
    // Final cursor position per bank, indexed by BankId ordinal. After a
    // successful traversal each equals that bank's planned length.
    std::array<std::size_t, 11> final_cursors{};

    // Throws std::invalid_argument for a BankId not in kAllBankIds.
    std::size_t final_cursor(BankId id) const;
};

// A single-use, forward-only read cursor over one ParameterBanks. There is
// deliberately no reset, rewind, or refill operation: once a value has been
// handed out it can never be handed out again, and a bank that has run out
// cannot be extended.
class ParameterCursor {
public:
    // `banks` must outlive the cursor; it is referenced, not copied.
    explicit ParameterCursor(const ParameterBanks& banks) : banks_(&banks) {}

    const ParameterBanks& banks() const { return *banks_; }

    // Consumes the next element of each of the 3 susceptible banks.
    // Throws std::runtime_error naming the offending bank if any susceptible
    // cursor is already at the end of its bank (overrun, including a second
    // traversal of an exhausted cursor).
    SusceptibleParameters take_susceptible();

    // Consumes the next element of each of the 8 infector banks. Throws
    // std::runtime_error naming the offending bank on overrun.
    InfectorParameters take_infector();

    // Current position of one bank's cursor. Throws std::invalid_argument
    // for a BankId not in kAllBankIds.
    std::size_t cursor(BankId id) const;

    // True when every one of the 11 cursors equals its bank's length.
    bool exhausted() const;

    // Throws std::runtime_error naming the first under-consumed bank unless
    // every cursor exactly equals its bank's length.
    void finalize() const;

private:
    friend ReplicateAssignments assign_replicate(ParameterCursor& cursor,
                                                   const SamplingPlan& plan);
    const ParameterBanks* banks_;
    std::array<std::size_t, 11> cursors_{};
    bool traversal_started_ = false;
};

// Traverses `plan` in simulation order using `cursor`, emitting one row per
// simulation: a zero-infector row consumes nothing, and an active row
// consumes one element from each susceptible bank followed by
// counts[i] elements from each infector bank. Calls cursor.finalize()
// afterwards, so a successful return proves every cursor is exactly
// exhausted.
//
// Throws std::invalid_argument if `plan` is internally inconsistent (see
// build_parameter_banks) or if its n / n_active / k_total do not match the
// cursor's banks. Throws std::runtime_error on any cursor overrun (including
// a second traversal of the same cursor) or if the traversal ends with any
// cursor short of its bank length.
ReplicateAssignments assign_replicate(ParameterCursor& cursor, const SamplingPlan& plan);

// Convenience overload: constructs a fresh single-use cursor over `banks`
// and traverses `plan` once.
ReplicateAssignments assign_replicate(const ParameterBanks& banks, const SamplingPlan& plan);
}  // namespace sampling_plan
}  // namespace ecai_cpp2

#endif  // ECAI_CPP2_SAMPLING_PLAN_H
