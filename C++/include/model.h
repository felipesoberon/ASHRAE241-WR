// Stage 8: Jones QER calculation -- a pure, deterministic function of
// already-sampled infector values.
//
// Source of truth: docs/MODEL_SPECIFICATION.md section 3 (Eq. 3 QER_i, Eq. 4
// VER_i, Eq. 5 V*_drop), section 5 / Table 1 (constants and hard bounds),
// section 7 / Table 3 (printed x1e4 C_drop and x1e-6 d_bar scaling), and
// section 10 (dimensional analysis). IMPLEMENTATION_PLAN.md Section 7,
// Stage 8.
//
// This file computes QER_i only, from a caller-supplied
// sampling_plan::InfectorParameters. It performs no random-number
// generation, no LHS draws, no category lookup, and holds no global mutable
// state: every input is a value the caller already sampled (Stage 5/6/7).

#ifndef ECAI_CPP2_MODEL_H
#define ECAI_CPP2_MODEL_H

#include <vector>

#include "sampling_plan.h"

namespace ecai_cpp2 {
namespace model {

// Fractional reduction of viability after emission ("R" in Eq. 4; Table 1,
// p.15 fixes this at 1.0). Frozen constant, never sampled. See
// MODEL_SPECIFICATION.md Discrepancy D1 for the R-vs-RD typesetting note:
// this is the prose/Table-1 symbol R, not RTD (Eq. 3's separate,
// independently-sampled respiratory tract deposition/absorption fraction).
inline constexpr double kViabilityAfterEmissionR = 1.0;

// Every frozen intermediate IMPLEMENTATION_PLAN.md Stage 8 requires exposed,
// plus the final quanta emission rate. Units:
//   vdrop             : dimensionless (Eq. 5; m^3 particle volume per m^3
//                       air, after applying Table 3's printed x1e4 C_drop
//                       and x1e-6 d_bar scalings)
//   gvl_copies_per_m3 : RNA copies m^-3 (Eq. 4's GVL_i, after the Table 1
//                       footnote's x1e6 mL^-1 -> m^-3 conversion)
//   ver_virions_per_h : virions h^-1 per person (Eq. 4)
//   qer_per_h         : quanta h^-1 per person (Eq. 3)
struct QERResult {
    double vdrop = 0.0;
    double gvl_copies_per_m3 = 0.0;
    double ver_virions_per_h = 0.0;
    double qer_per_h = 0.0;
};

// Computes one infector's QER_i (Eq. 3) from already-sampled parameters,
// applying Eq. 5 (V*_drop) and Eq. 4 (VER_i) as intermediates. Pure
// function: no randomness, no allocation beyond the return value, no
// reference to any occupancy category table.
//
// Applies, exactly once each:
//   - Table 3's printed C_drop scaling (x1e4 -> # m^-3) and d_bar scaling
//     (x1e-6 -> m), before cubing (d_bar * E)^3 in Eq. 5;
//   - Table 1's GVL_i footnote scaling (x1e6, RNA copies mL^-1 -> m^-3);
//   - the frozen viability-after-emission constant R = 1.0 in Eq. 4;
//   - RTD and DK exactly as sampled, per Eq. 3.
//
// Throws std::invalid_argument if any field of `infector` is nonfinite, or
// falls outside the physical domain MODEL_SPECIFICATION.md section 5 /
// Table 1 documents:
//   - pbr_m3h, c_drop, d_bar, gvl_copies_per_ml: unbounded (lognormal /
//     normal) distributions in Table 1, so only finiteness and strict
//     positivity are physically required;
//   - evaporation_factor: Table 1's Beta hard bounds [2, 5];
//   - viable_fraction: Table 1's Beta hard bounds [1e-4, 1e-2];
//   - rtd: Table 1's Uniform hard bounds [0.43, 0.65];
//   - dose_constant: Table 1's Uniform hard bounds [5, 15].
// No input is ever clamped into range; an out-of-domain value is rejected.
QERResult compute_qer(const sampling_plan::InfectorParameters& infector);

// ===================================================================
// Stage 9: ECAi calculation (frozen inversion) and the zero-QER case.
//
// Source of truth: docs/MODEL_SPECIFICATION.md section 3 ("Frozen ECAi
// inversion used by C++2", the algebraic rearrangement of Eqs. 1, 2, 6 and
// 7), section 10 (dimensional analysis: both inversion terms in m^3/h, and
// division by 3.6 * I0 yielding L/s/person), and section 12 (frozen
// defaults). IMPLEMENTATION_PLAN.md Section 7, Stage 9 and section 4.4.
//
// This is pure algebra over caller-supplied values. It owns no parameter
// bank and no cursor, and therefore cannot and does not implement the
// zero-infector consumption shortcut: orchestration (sampling_plan.cpp /
// Stage 10 simulation) must short-circuit a zero-infector row before
// requesting any susceptible or infector parameter. What is implemented
// here is only the algebraic half of that rule -- an explicitly supplied
// QER_sum of zero yields exactly 0 L/s/person.
// ===================================================================

// Converts m^3/h to L/s (1 m^3/h = 1/3.6 L/s), applied together with the
// per-person division by I0 in Eq. 7. Frozen constant, never sampled.
inline constexpr double kM3PerHourPerLitrePerSecond = 3.6;

// One active simulation's non-QER inputs to the frozen inversion. Units:
//   susceptible_pbr_m3h : m^3/h   (PBR_s, Eq. 2)
//   gamma_per_h         : 1/h     (surface deposition rate, Eq. 6)
//   lambda_per_h        : 1/h     (biological decay, Eq. 6)
//   volume_m3           : m^3     (VOL)
//   occupants           : count   (I0, Eq. 7)
//   mask_efficiency     : -       (mu, Eq. 2)
//   duration_h          : h       (D, Eq. 2)
//   target_probability  : -       (P(I), Eq. 1)
struct ECAiInputs {
    double susceptible_pbr_m3h = 0.0;
    double gamma_per_h = 0.0;
    double lambda_per_h = 0.0;
    double volume_m3 = 0.0;
    long long occupants = 0;
    double mask_efficiency = 0.0;
    double duration_h = 0.0;
    double target_probability = 0.0;
};

// Every intermediate of the frozen inversion, exposed so the result can be
// audited term by term against the publication. Units:
//   q_target             : dimensionless (Eq. 1 inverted)
//   mask_factor          : dimensionless
//   qer_sum_per_h        : quanta/h (sum over the realized infectors)
//   term1_m3_per_h       : m^3/h
//   term2_m3_per_h       : m^3/h
//   t_ecai_raw_m3_per_h  : m^3/h  (T_ECAi before clamping)
//   ecai_raw_l_s_person  : L/s/person, unclamped (may be negative)
//   ecai_l_s_person      : L/s/person, clamped at exactly 0
//   clamped              : true iff the raw value was negative
struct ECAiResult {
    double q_target = 0.0;
    double mask_factor = 0.0;
    double qer_sum_per_h = 0.0;
    double term1_m3_per_h = 0.0;
    double term2_m3_per_h = 0.0;
    double t_ecai_raw_m3_per_h = 0.0;
    double ecai_raw_l_s_person = 0.0;
    double ecai_l_s_person = 0.0;
    bool clamped = false;
};

// Computes one simulation's ECAi from a pre-summed QER. Pure function: no
// randomness, no global mutable state, no bank or cursor access.
ECAiResult compute_ecai_from_qer_sum(const ECAiInputs& inputs, double qer_sum_per_h);

// Convenience overload: sums `qer_per_h` in the supplied order (so the sum
// is deterministic for a fixed infector ordering) and delegates to
// compute_ecai_from_qer_sum. An empty vector is the explicit zero-QER case
// and yields exactly 0 L/s/person.
ECAiResult compute_ecai(const ECAiInputs& inputs, const std::vector<double>& qer_per_h);

// ===================================================================
// Probability mode: forward Wells-Riley infection-probability evaluation at
// a fixed input ECAi.
//
// Source of truth: docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md section 3.2
// ("New forward evaluation (probability mode)"). This is the algebraic
// mirror of compute_ecai_from_qer_sum/compute_ecai above: instead of
// inverting Eq. 1 for the ECAi that reproduces a target probability, it
// evaluates Eq. 1 forward for the probability produced by a caller-supplied
// ECAi. It is deliberately decoupled from any target probability: it takes
// no target, validates none, and reads none. Pure function: no randomness,
// no global mutable state, no bank or cursor access.
// ===================================================================

// One active simulation's non-QER inputs to the forward evaluation. Same
// seven physical fields as ECAiInputs (everything except
// target_probability), plus the fixed input ECAi the forward equation
// consumes instead. Units:
//   susceptible_pbr_m3h   : m^3/h   (PBR_s, Eq. 2)
//   gamma_per_h           : 1/h     (surface deposition rate, Eq. 6)
//   lambda_per_h          : 1/h     (biological decay, Eq. 6)
//   volume_m3             : m^3     (VOL)
//   occupants             : count   (I0, Eq. 7)
//   mask_efficiency       : -       (mu, Eq. 2)
//   duration_h            : h       (D, Eq. 2)
//   input_ecai_l_s_person : L/s/person, finite and >= 0 (the fixed input
//                           ECAi; Eq. 7 forward)
struct ProbabilityInputs {
    double susceptible_pbr_m3h = 0.0;
    double gamma_per_h = 0.0;
    double lambda_per_h = 0.0;
    double volume_m3 = 0.0;
    long long occupants = 0;
    double mask_efficiency = 0.0;
    double duration_h = 0.0;
    double input_ecai_l_s_person = 0.0;
};

// Every intermediate of the forward evaluation, exposed so the result can be
// audited term by term. Units:
//   eca_total_m3_per_h : m^3/h    (Eq. 7 forward: 3.6 * I0 * input ECAi)
//   removal_m3_per_h   : m^3/h    (total removal, ECA_total + VOL*(gamma+lambda))
//   mask_factor        : dimensionless ((1 - mu)^2)
//   qer_sum_per_h      : quanta/h (sum over the realized infectors)
//   dose_numerator     : m^3 quanta / h (PBR_s * D * mask_factor * QER_sum)
//   quanta_inhaled     : quanta, dimensionless (Q = dose_numerator / removal)
//   probability        : dimensionless, in the closed range [0, 1], never
//                        -0.0 (P = -expm1(-Q), sign-normalized to +0.0)
//   saturated          : true iff probability == 1.0 exactly, i.e. Q is
//                        large enough that -expm1(-Q) rounds to 1.0 in
//                        double precision (numerical saturation, not an
//                        error; see the plan's saturation policy)
struct ProbabilityResult {
    double eca_total_m3_per_h = 0.0;
    double removal_m3_per_h = 0.0;
    double mask_factor = 0.0;
    double qer_sum_per_h = 0.0;
    double dose_numerator = 0.0;
    double quanta_inhaled = 0.0;
    double probability = 0.0;
    bool saturated = false;
};

// Computes one simulation's forward infection probability from a pre-summed
// QER, at the caller-supplied fixed input ECAi. Pure function: no
// randomness, no global mutable state, no bank or cursor access.
ProbabilityResult compute_probability_from_qer_sum(const ProbabilityInputs& inputs,
                                                    double qer_sum_per_h);

// Convenience overload: sums `qer_per_h` in the supplied order (matching
// compute_ecai's summation policy) and delegates to
// compute_probability_from_qer_sum. An empty vector is the explicit
// zero-QER case and yields exactly +0.0 probability.
ProbabilityResult compute_probability(const ProbabilityInputs& inputs,
                                       const std::vector<double>& qer_per_h);

}  // namespace model
}  // namespace ecai_cpp2

#endif  // ECAI_CPP2_MODEL_H
