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
#include "qer_direct.h"
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

// Infector-handling policy for one Stage 6 plan
// (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md). Defined independently here
// rather than reused from simulation.h, which includes this header (a
// reverse include would be circular). The enumerator values are frozen to
// agree with simulation::InfectorPolicy exactly, so a caller may convert
// between the two with a plain static_cast rather than a lookup table.
// `community` is enumerator 0, so a value-initialized SamplingPlan defaults
// to the unchanged current binomial behavior without any explicit
// initialization.
enum class InfectorPolicy { community = 0, risk_averse = 1 };
static_assert(static_cast<int>(InfectorPolicy::community) == 0,
              "InfectorPolicy::community must be enumerator 0 so value-initialization defaults "
              "to it");

// Version tag for the conditional-PPF algorithm/table formulation used to
// realize a risk-averse plan's counts: distributions::BinomialPpfTable's
// retained conditional-tail cumulative-mass table
// (conditional_ppf_at_least_one), not a redraw-until-positive or affine
// endpoint construction. Recorded verbatim in InfectorPolicyProvenance so a
// later reader can tell which algorithm produced a stored zero_mass/cutoff
// without re-deriving it from the code.
inline constexpr const char* kConditionalPpfAlgorithmVersion = "conditional-tail-v1";

// Risk-averse-specific numerical provenance for one Stage 6 plan
// (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md Phase 4, section 8's
// auditability requirement). Meaningful only when the owning SamplingPlan's
// `policy` is risk_averse; every field stays at its default for a community
// plan, matching InfectorPolicy::community's own zero-valued default.
struct InfectorPolicyProvenance {
    InfectorPolicy policy = InfectorPolicy::community;
    // distributions::BinomialPpfTable::zero_mass() for this plan's (I0, p)
    // -- the conditional-inversion cutoff -- captured before any LHS value
    // is consumed. Stays 0.0 (unused) for a community plan.
    double zero_mass = 0.0;
    // True iff zero_mass rounded to exactly 0.0 in binary64 (the p = 1
    // deterministic branch, or an underflowed pmf[0] at p < 1).
    bool zero_mass_rounded_to_zero = false;
    // True iff zero_mass rounded to exactly 1.0 in binary64 for a
    // representable positive p, triggering conditional_ppf_at_least_one's
    // explicit limiting branch (returns 1 for every u) rather than an
    // unresolvable-tail table lookup.
    bool zero_mass_rounded_to_one = false;
    // kConditionalPpfAlgorithmVersion, recorded only for a risk-averse plan;
    // empty for a community plan.
    std::string conditional_ppf_algorithm_version;
    // True iff every count in this plan was produced by a deterministic or
    // limiting branch of conditional_ppf_at_least_one rather than an
    // ordinary conditional-tail table lookup: p == 1 (point mass at I0) or
    // zero_mass rounding to exactly 1.0. Always false for a community plan.
    bool deterministic_branch_used = false;
};

// The full Stage 6 sampling plan for one (configuration, replicate) pair:
// exactly `n` infector counts, in simulation order, plus their accounting.
struct SamplingPlan {
    // counts[i] is simulation i's infector count (ninf), in the same order
    // as the underlying LHS list; simulation order is preserved. Every
    // count is in [0, I0] under InfectorPolicy::community and in [1, I0]
    // under InfectorPolicy::risk_averse.
    std::vector<long long> counts;
    std::size_t n = 0;
    std::size_t n_zero = 0;
    std::size_t n_active = 0;
    std::size_t k_total = 0;
    std::size_t susceptible_draws = 0;  // == n_active
    std::size_t infector_draws = 0;     // == k_total

    // The infector policy this plan was built under (docs/RISK_AVERSE_
    // IMPLEMENTATION_PLAN.md Phase 4). `community` (enumerator 0) is the
    // default, so a value-initialized SamplingPlan is a community plan.
    InfectorPolicy policy = InfectorPolicy::community;
    // Risk-averse numerical provenance; left at its all-default value for a
    // community plan.
    InfectorPolicyProvenance infector_provenance;

    // Explicit count-stream invariants, tracked independently of
    // counts.size(): the number of count-stream uniforms actually generated
    // by the single StreamId::infector_count / list_length = n LHS call,
    // and the number actually consumed to produce `counts`. Both equal `n`
    // exactly for every policy -- there is no retry/rejection loop that
    // could generate or consume a variable number of count-stream uniforms
    // -- and are tracked separately from parameter-bank/QER consumption
    // counters (susceptible_draws/infector_draws) below.
    std::size_t lhs_generated = 0;
    std::size_t lhs_consumed = 0;
};

// Generates exactly `n` open-unit 1D LHS values (seeded from
// seeds.derive_seed(config, replicate_index, StreamId::infector_count, n)
// -- the same seed regardless of `infector_policy`, so the count-stream is
// matched across policies), transforms each through one per-call
// distributions::BinomialPpfTable(i0, p) to obtain one infector count per
// simulation, and returns the resulting plan.
//
// Under the default InfectorPolicy::community, each u is transformed by
// table.ppf(u) -- equal to distributions::binomial_ppf(u, i0, p) exactly,
// with the constant (i0, p) PMF/CDF built once rather than once per u --
// preserving the current community behavior, including zero-infector rows,
// exactly and unconditionally.
//
// Under InfectorPolicy::risk_averse, each u is instead transformed by
// table.conditional_ppf_at_least_one(u), producing exactly `n` counts each
// in [1, i0] with no zero-infector row and no redraw/rejection loop --
// exactly `n` count-stream uniforms are generated and exactly `n` are
// consumed, identically to the community case. Because the seed and LHS
// list are unchanged by policy, a risk-averse count is never lower than the
// community count that the same u would have produced. The plan's
// `infector_provenance` records the cutoff (zero_mass), whether it rounded
// to 0.0/1.0, the conditional-PPF algorithm version, and whether a
// deterministic/limiting branch was used; these all stay at their default,
// unused values for a community plan.
//
// Every check below runs, in order, before any allocation:
//   - n == 0: std::invalid_argument.
//   - i0 < 0: std::invalid_argument.
//   - p nonfinite or outside [0, 1]: std::invalid_argument.
//   - infector_policy is not a declared InfectorPolicy enumerator:
//     std::invalid_argument.
//   - infector_policy is risk_averse and i0 < 1, or p is not strictly
//     positive: std::invalid_argument (the authoritative direct-call check;
//     it applies identically whether or not a CLI/scenario layer already
//     performed the same check).
//   - n > policy.max_simulations: std::runtime_error.
//   - checked n * static_cast<size_t>(i0) (K_total's maximum possible
//     value) overflowing std::size_t, or exceeding
//     policy.max_total_infectors: std::runtime_error.
//   - checked peak byte cost of the n-length uniform/count vectors plus the
//     three `(i0+1)`-double binomial buffers retained/transiently owned by
//     distributions::BinomialPpfTable (PMF, CDF, and conditional-tail --
//     the conditional-tail buffer is built unconditionally by the table's
//     constructor for every policy, not only for risk_averse, so it is
//     charged against every call), overflowing std::size_t or exceeding
//     policy.max_bytes, or n exceeding std::vector<double>::max_size() /
//     std::vector<long long>::max_size(): std::runtime_error.
//   - checked i0 + 1 (distributions::BinomialPpfTable's PMF/CDF/
//     conditional-tail buffer size) overflowing std::size_t or exceeding
//     std::vector<double>::max_size(): std::runtime_error.
//
// Zero-infector simulations (counts[i] == 0; only possible under
// InfectorPolicy::community) are accounted for but imply no other
// per-simulation draw; susceptible_draws counts only active simulations and
// infector_draws is the checked sum of every positive count (see
// summarize_counts).
SamplingPlan build_sampling_plan(const seed_manager::SeedManager& seeds,
                                  const seed_manager::EffectiveConfig& config,
                                  std::uint64_t replicate_index, std::size_t n, long long i0,
                                  double p, const ResourcePolicy& policy = ResourcePolicy{},
                                  InfectorPolicy infector_policy = InfectorPolicy::community);

// Validates a completed plan's policy-specific invariants before it is
// trusted by a caller expecting `expected_policy` (docs/RISK_AVERSE_
// IMPLEMENTATION_PLAN.md Phase 4: "policy is stored in the plan and checked
// against the resolved scenario during plan validation and replicate
// execution"). Runs, in order:
//   - expected_policy or plan.policy is not a declared InfectorPolicy
//     enumerator: std::invalid_argument.
//   - plan.policy != expected_policy: std::invalid_argument (a tampered or
//     mismatched scenario/plan policy pair).
//   - plan.counts.size() != plan.n, or plan's accounting disagrees with an
//     independent summarize_counts(plan.counts): std::invalid_argument.
//   - compute_bank_budget(plan.n, plan.n_active, plan.k_total, policy)
//     rejects the plan (an inconsistent accounting triple or a resource
//     ceiling overrun): propagates its std::invalid_argument or
//     std::runtime_error.
//   - plan.policy == risk_averse and (plan.n_zero != 0, plan.n_active !=
//     plan.n, plan.k_total < plan.n, or any element of plan.counts < 1):
//     std::invalid_argument (a risk-averse plan must contain no
//     zero-infector rows, so every simulation is active and K_total is at
//     least N).
// A community plan is exempt from the last group of checks: it may
// legitimately contain zero-infector rows.
void validate_plan_policy_context(const SamplingPlan& plan, InfectorPolicy expected_policy,
                                   const ResourcePolicy& policy = ResourcePolicy{});


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


// ===================================================================
// Package 3a: family-aware Stage 7 bank allocation and accounting.
//
// Extends Stage 7 so Jones-fitted and Mikszewski QER families (docs/
// QER_DISTRIBUTION_IMPLEMENTATION_PLAN.md Phase 3) allocate only the 3
// susceptible-sized banks plus one direct-QER bank of length K_total,
// without allocating, seeding, or iterating the 8 full-Jones infector banks.
// The full-Jones family continues to use build_parameter_banks/
// ParameterBanks/ParameterCursor/assign_replicate above, byte-identically.
// ===================================================================

// The family-specific persistent-bank formula (Phase 3):
//   Jones:                     N + 3*N_active + 8*K_total (delegates to
//                               compute_bank_budget, unchanged).
//   JonesFitted / Mikszewski:  N + 3*N_active + K_total, plus the transient
//                               K_total-length LHS uniform buffer alive
//                               alongside the direct-QER bank during its one
//                               generation pass.
// Checked once, here, before any bank is allocated; never re-checked in a
// hot loop. Runs, in order:
//   - `family` is not a declared qer_direct::QerFamily enumerator:
//     std::invalid_argument.
//   - n_active > n, k_total < n_active, or (n_active == 0 && k_total != 0):
//     std::invalid_argument.
//   - k_total > policy.max_total_infectors: std::runtime_error.
//   - n, n_active, or k_total exceeding std::vector<double>::max_size():
//     std::runtime_error.
//   - checked combined-element and combined-byte (including the transient
//     buffer noted above) arithmetic overflowing std::size_t, or the byte
//     total exceeding policy.max_bytes: std::runtime_error.
BankBudget compute_family_bank_budget(qer_direct::QerFamily family, std::size_t n,
                                       std::size_t n_active, std::size_t k_total,
                                       const ResourcePolicy& policy = ResourcePolicy{});

// One exact-size direct-QER bank: length K_total, one LHS draw and one
// qer_direct::sample_direct_qer_log10normal transform per realized infector,
// generated exactly once. Meaningful only when the owning
// FamilyParameterBanks' family is JonesFitted or Mikszewski; left default
// (ungenerated, empty) for a Jones plan.
struct DirectQerBank {
    seed_manager::StreamId stream = seed_manager::StreamId::direct_qer;
    std::uint64_t seed = 0;
    bool generated = false;
    std::vector<double> values;  // quanta/h, length K_total when generated.

    // Stable id for one direct-bank value: its index (there is only ever one
    // direct bank per replicate, so no bank-ordinal packing is needed).
    // Throws std::out_of_range if index >= values.size().
    std::uint64_t value_id(std::size_t index) const;
};

// The complete family-aware Stage 7 bank set for one (configuration,
// replicate) pair. Exactly one of {infector_banks, direct_qer_bank} is ever
// generated -- the other is left entirely default/ungenerated -- so no
// unused Jones infector bank or direct-QER bank is ever materialized for the
// wrong family.
struct FamilyParameterBanks {
    qer_direct::QerFamily family = qer_direct::QerFamily::Jones;
    // The 3 susceptible-sized banks, in {susceptible_pbr, gamma, lambda}
    // order, generated identically for every family (same seeds/order/values
    // build_parameter_banks would produce for the same plan/category).
    std::array<ParameterBank, 3> susceptible_banks{};
    // The 8 full-Jones infector-sized banks, in {infector_pbr, c_drop, d_bar,
    // E, GVL, VF, RTD, DK} order. Generated only when family == Jones; left
    // default (ungenerated, empty banks) otherwise.
    std::array<ParameterBank, 8> infector_banks{};
    // The single direct-QER bank. Generated only when family != Jones; left
    // default (ungenerated, empty) for a Jones plan.
    DirectQerBank direct_qer_bank;
    std::size_t n = 0;
    std::size_t n_active = 0;
    std::size_t k_total = 0;
    BankBudget budget;
};

// Generates the family-aware Stage 7 banks for one completed Stage 6 plan.
//
// `family` selects the layout: Jones builds the 3 susceptible banks and the
// 8 full-Jones infector banks exactly as build_parameter_banks does (same
// seeds, order, and values for the same plan/category/replicate); JonesFitted
// and Mikszewski build the 3 susceptible banks and exactly one direct-QER
// bank of length K_total from `direct_params`, seeded from
// seeds.derive_seed(config, replicate_index, StreamId::direct_qer, k_total),
// and never construct, seed, or generate any of the 8 full-Jones infector
// banks.
//
// Every check below runs before any bank is allocated:
//   - `family` is not a declared qer_direct::QerFamily enumerator:
//     std::invalid_argument (rejects an invalid/mixed-family request before
//     any allocation).
//   - family != Jones and direct_params has a nonfinite mu_log10 or a
//     nonfinite/non-positive sigma_log10: std::invalid_argument.
//   - plan.counts.size() != plan.n, or plan accounting disagreeing with a
//     fresh summarize_counts(plan.counts): std::invalid_argument.
//   - compute_family_bank_budget(family, plan.n, plan.n_active, plan.k_total,
//     policy): propagates std::runtime_error on any overflow or ceiling
//     overrun, leaving no partial bank output.
// A zero-length bank (K_total == 0 for the direct bank, or N_active == 0 for
// the susceptible banks) is left empty and ungenerated: no seed is derived
// and no LHS list is generated for it.
FamilyParameterBanks build_family_parameter_banks(
    const seed_manager::SeedManager& seeds, const seed_manager::EffectiveConfig& config,
    std::uint64_t replicate_index, const SamplingPlan& plan, const OccupancyCategory& category,
    qer_direct::QerFamily family, const qer_direct::DirectQerParams& direct_params,
    const ResourcePolicy& policy = ResourcePolicy{});

// One realized infector's family-aware value: either its full eight Jones
// parameters or its single direct QER value, never both.
struct FamilyInfectorAssignment {
    bool is_direct = false;
    InfectorParameters jones;              // meaningful iff !is_direct
    double direct_qer_quanta_per_h = 0.0;  // meaningful iff is_direct
    std::size_t bank_index = 0;            // position consumed in the relevant bank(s)
    std::uint64_t id = 0;                  // stable consumption id
};

// One simulation row under family-aware assignment. A zero-infector row
// (infector_count == 0) has `active == false`, a default-valued
// `susceptible`, and an empty `infectors`: it consumes nothing from any bank.
struct FamilySimulationAssignment {
    std::size_t simulation_index = 0;
    long long infector_count = 0;
    bool active = false;
    SusceptibleParameters susceptible;
    std::vector<FamilyInfectorAssignment> infectors;  // size == infector_count
};

// A completed, fully accounted family-aware traversal of one replicate.
struct FamilyReplicateAssignments {
    qer_direct::QerFamily family = qer_direct::QerFamily::Jones;
    std::vector<FamilySimulationAssignment> simulations;
    std::size_t n_active = 0;
    std::size_t k_total = 0;
    // Final cursor positions after a successful traversal: the susceptible
    // cursor (all 3 susceptible banks advance together) and, depending on
    // family, either the infector cursor (the 8 Jones banks, advancing
    // together) or the direct cursor (the single direct-QER bank). The
    // cursor not used by this plan's family stays 0.
    std::size_t final_susceptible_cursor = 0;
    std::size_t final_infector_cursor = 0;
    std::size_t final_direct_cursor = 0;
};

// A single-use, forward-only, family-aware read cursor over one
// FamilyParameterBanks. As with ParameterCursor there is no reset, rewind, or
// refill: once handed out, a value can never be handed out again.
class FamilyParameterCursor {
public:
    // `banks` must outlive the cursor; it is referenced, not copied.
    explicit FamilyParameterCursor(const FamilyParameterBanks& banks) : banks_(&banks) {}

    const FamilyParameterBanks& banks() const { return *banks_; }

    // Consumes the next element of each of the 3 susceptible banks. Throws
    // std::runtime_error on overrun.
    SusceptibleParameters take_susceptible();

    // Consumes the next element of each of the 8 Jones infector banks.
    // Throws std::runtime_error if banks().family != Jones, or on overrun.
    InfectorParameters take_jones_infector();

    // Consumes the next element of the direct-QER bank. Throws
    // std::runtime_error if banks().family == Jones, or on overrun.
    double take_direct_infector();

    // True when the susceptible cursor and the one family-relevant infector/
    // direct cursor each equal their bank's length.
    bool exhausted() const;

    // Throws std::runtime_error naming the under-consumed bank group unless
    // every family-relevant cursor exactly equals its bank's length -- the
    // total-consumed-equals-allocated reconciliation.
    void finalize() const;

private:
    friend FamilyReplicateAssignments assign_family_replicate(FamilyParameterCursor& cursor,
                                                                const SamplingPlan& plan);
    const FamilyParameterBanks* banks_;
    std::size_t susceptible_cursor_ = 0;
    std::size_t infector_cursor_ = 0;
    std::size_t direct_cursor_ = 0;
    bool traversal_started_ = false;
};

// Traverses `plan` in simulation order using `cursor`: a zero-infector row
// consumes nothing; an active row consumes one susceptible draw followed by
// counts[i] infector draws, dispatched once by banks().family (all Jones or
// all direct within one traversal -- never mixed). Calls cursor.finalize()
// afterwards, so a successful return proves every family-relevant cursor is
// exactly exhausted.
//
// Throws std::invalid_argument if `plan` is internally inconsistent, or if
// its n / n_active / k_total do not match the cursor's banks. Throws
// std::runtime_error on any cursor overrun (including a second traversal of
// the same cursor) or if the traversal ends with any relevant cursor short of
// its bank length.
FamilyReplicateAssignments assign_family_replicate(FamilyParameterCursor& cursor,
                                                     const SamplingPlan& plan);

// Convenience overload: constructs a fresh single-use cursor over `banks`
// and traverses `plan` once.
FamilyReplicateAssignments assign_family_replicate(const FamilyParameterBanks& banks,
                                                     const SamplingPlan& plan);
}  // namespace sampling_plan
}  // namespace ecai_cpp2

#endif  // ECAI_CPP2_SAMPLING_PLAN_H
