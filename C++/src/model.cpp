#include "model.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace ecai_cpp2 {
namespace model {

namespace {

constexpr double kPi = 3.14159265358979323846;

void require_finite(double value, const char* field_name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string("InfectorParameters.") + field_name +
                                     " must be finite");
    }
}

void require_positive(double value, const char* field_name) {
    require_finite(value, field_name);
    if (!(value > 0.0)) {
        throw std::invalid_argument(std::string("InfectorParameters.") + field_name +
                                     " must be strictly positive");
    }
}

void require_in_range(double value, double lo, double hi, const char* field_name) {
    require_finite(value, field_name);
    if (value < lo || value > hi) {
        throw std::invalid_argument(std::string("InfectorParameters.") + field_name +
                                     " must lie within its Table 1 hard bounds [" +
                                     std::to_string(lo) + ", " + std::to_string(hi) + "]");
    }
}


// Stage 9 validation helpers. Separate from the Stage 8 ones only so the
// thrown message names ECAiInputs rather than InfectorParameters; the
// policy is identical (finiteness first, then the physical domain, and no
// input is ever clamped into range).
void require_ecai_finite(double value, const char* field_name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string("ECAiInputs.") + field_name +
                                     " must be finite");
    }
}

void require_ecai_positive(double value, const char* field_name) {
    require_ecai_finite(value, field_name);
    if (!(value > 0.0)) {
        throw std::invalid_argument(std::string("ECAiInputs.") + field_name +
                                     " must be strictly positive");
    }
}

void require_ecai_in_range(double value, double lo, double hi, const char* field_name) {
    require_ecai_finite(value, field_name);
    if (value < lo || value > hi) {
        throw std::invalid_argument(std::string("ECAiInputs.") + field_name +
                                     " must lie within [" + std::to_string(lo) + ", " +
                                     std::to_string(hi) + "]");
    }
}

}  // namespace

// One supplied quanta emission rate (or pre-summed QER_sum): finite and
// nonnegative; a -0.0 entry is a legitimate zero. Declared in model.h (not
// file-local) so simulation.cpp's per-row loops can apply the same check
// while accumulating QER_i values without a temporary std::vector<double>.
void require_qer_value(double value, const char* what) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(what) + " must be finite");
    }
    if (!(value >= 0.0)) {
        throw std::invalid_argument(std::string(what) + " must be nonnegative");
    }
}

namespace {

// Every computed intermediate must stay finite: an overflow to +/-Inf (or a
// NaN) is reported rather than returned as an ECAi value.
void require_finite_intermediate(double value, const char* what) {
    if (!std::isfinite(value)) {
        throw std::runtime_error(std::string("ECAi ") + what +
                                  " is not finite (overflow or invalid intermediate)");
    }
}

// Full domain check for the eight ECAi inputs, run before any arithmetic.
// Bounds: PBR_s / lambda / VOL strictly positive (unbounded lognormal
// families, MODEL_SPECIFICATION.md section 5 / Table 1, so only positivity
// is physically required); gamma inside its Table 1 Uniform hard bounds
// [0.42, 0.61]; I0 a positive count; mu in [0, 1]; D strictly positive;
// P(I) strictly inside the open interval (0, 1)
// (IMPLEMENTATION_PLAN.md section 5.2).
void validate_ecai_inputs(const ECAiInputs& inputs) {
    require_ecai_positive(inputs.susceptible_pbr_m3h, "susceptible_pbr_m3h");
    require_ecai_in_range(inputs.gamma_per_h, 0.42, 0.61, "gamma_per_h");
    require_ecai_positive(inputs.lambda_per_h, "lambda_per_h");
    require_ecai_positive(inputs.volume_m3, "volume_m3");
    if (inputs.occupants <= 0) {
        throw std::invalid_argument("ECAiInputs.occupants must be a positive count");
    }
    require_ecai_in_range(inputs.mask_efficiency, 0.0, 1.0, "mask_efficiency");
    require_ecai_positive(inputs.duration_h, "duration_h");
    require_ecai_finite(inputs.target_probability, "target_probability");
    if (!(inputs.target_probability > 0.0) || !(inputs.target_probability < 1.0)) {
        throw std::invalid_argument(
            "ECAiInputs.target_probability must lie strictly inside the open interval (0, 1)");
    }
}

}  // namespace

QERResult compute_qer(const sampling_plan::InfectorParameters& infector) {
    require_positive(infector.pbr_m3h, "pbr_m3h");
    require_positive(infector.c_drop, "c_drop");
    require_positive(infector.d_bar, "d_bar");
    require_in_range(infector.evaporation_factor, 2.0, 5.0, "evaporation_factor");
    require_positive(infector.gvl_copies_per_ml, "gvl_copies_per_ml");
    require_in_range(infector.viable_fraction, 1e-4, 1e-2, "viable_fraction");
    require_in_range(infector.rtd, 0.43, 0.65, "rtd");
    require_in_range(infector.dose_constant, 5.0, 15.0, "dose_constant");

    // Table 3 printed-unit scalings (applied exactly once): C_drop is
    // sampled/stored in units of x1e4 m^-3; d_bar is sampled/stored in units
    // of x1e-6 m.
    const double c_drop_per_m3 = infector.c_drop * 1e4;
    const double d_bar_m = infector.d_bar * 1e-6;

    // Eq. 5: V*_drop = (pi/6) * (d_bar * E)^3 * C_drop.
    const double hydrated_diameter_m = d_bar_m * infector.evaporation_factor;
    const double vdrop =
        (kPi / 6.0) * hydrated_diameter_m * hydrated_diameter_m * hydrated_diameter_m *
        c_drop_per_m3;

    // Table 1 footnote: GVL_i is converted to RNA copies m^-3 by x1e6.
    const double gvl_copies_per_m3 = infector.gvl_copies_per_ml * 1e6;

    // Eq. 4: VER_i = PBR * V*_drop * GVL_i * VF * R.
    const double ver_virions_per_h = infector.pbr_m3h * vdrop * gvl_copies_per_m3 *
                                      infector.viable_fraction * kViabilityAfterEmissionR;

    // Eq. 3: QER_i = RTD * VER_i / DK.
    const double qer_per_h = infector.rtd * ver_virions_per_h / infector.dose_constant;

    QERResult result;
    result.vdrop = vdrop;
    result.gvl_copies_per_m3 = gvl_copies_per_m3;
    result.ver_virions_per_h = ver_virions_per_h;
    result.qer_per_h = qer_per_h;
    return result;
}


ECAiResult compute_ecai_from_qer_sum(const ECAiInputs& inputs, double qer_sum_per_h) {
    validate_ecai_inputs(inputs);
    require_qer_value(qer_sum_per_h, "ECAi QER_sum");

    // Q_target = -log(1 - target_probability), evaluated as -log1p(-P): the
    // two are algebraically identical, but log1p keeps full relative
    // precision for the small target probabilities the standard uses (the
    // frozen default is 0.001), where forming 1 - P first would discard
    // low-order bits and, below about 1e-16, collapse to exactly 1.0 and
    // make the divisor zero. Validation above guarantees 0 < P < 1, so
    // Q_target is strictly positive and the division below can never be by
    // zero.
    const double q_target = -std::log1p(-inputs.target_probability);
    require_finite_intermediate(q_target, "Q_target");

    // mask_factor = (1 - mu)^2 (Eq. 2), formed by one multiplication so the
    // factor is applied exactly once.
    const double one_minus_mu = 1.0 - inputs.mask_efficiency;
    const double mask_factor = one_minus_mu * one_minus_mu;

    const double term1_m3_per_h =
        inputs.susceptible_pbr_m3h * inputs.duration_h * mask_factor * qer_sum_per_h / q_target;
    require_finite_intermediate(term1_m3_per_h, "term1");

    const double term2_m3_per_h = inputs.volume_m3 * (inputs.gamma_per_h + inputs.lambda_per_h);
    require_finite_intermediate(term2_m3_per_h, "term2");

    const double t_ecai_raw_m3_per_h = term1_m3_per_h - term2_m3_per_h;
    require_finite_intermediate(t_ecai_raw_m3_per_h, "T_ECAi_raw");

    // Eq. 7: m^3/h -> L/s/person. occupants is a validated positive count,
    // so the divisor is strictly positive.
    const double ecai_raw_l_s_person =
        t_ecai_raw_m3_per_h /
        (kM3PerHourPerLitrePerSecond * static_cast<double>(inputs.occupants));
    require_finite_intermediate(ecai_raw_l_s_person, "ECAi_raw");

    ECAiResult result;
    result.q_target = q_target;
    result.mask_factor = mask_factor;
    result.qer_sum_per_h = qer_sum_per_h;
    result.term1_m3_per_h = term1_m3_per_h;
    result.term2_m3_per_h = term2_m3_per_h;
    result.t_ecai_raw_m3_per_h = t_ecai_raw_m3_per_h;
    result.ecai_raw_l_s_person = ecai_raw_l_s_person;
    result.clamped = ecai_raw_l_s_person < 0.0;
    // Only a strictly positive raw value passes through; everything else
    // becomes exactly +0.0, so a clamped (or -0.0) result can never carry a
    // negative sign bit into the percentile stage.
    result.ecai_l_s_person = (ecai_raw_l_s_person > 0.0) ? ecai_raw_l_s_person : 0.0;
    return result;
}

ECAiResult compute_ecai(const ECAiInputs& inputs, const std::vector<double>& qer_per_h) {
    // Plain left-to-right accumulation in the supplied order: the sum is a
    // deterministic function of the infector ordering the caller fixed, and
    // no compensated-summation reordering is introduced that the frozen
    // specification does not call for. Every term is validated nonnegative
    // and finite first, so the running sum is nondecreasing and can only
    // leave the finite range by overflowing to +Inf, which is checked.
    double qer_sum_per_h = 0.0;
    for (const double qer : qer_per_h) {
        require_qer_value(qer, "ECAi QER_i");
        qer_sum_per_h += qer;
    }
    if (!std::isfinite(qer_sum_per_h)) {
        throw std::runtime_error("ECAi QER_sum overflowed while summing the supplied QER values");
    }
    return compute_ecai_from_qer_sum(inputs, qer_sum_per_h);
}

// ===================================================================
// Probability mode: forward evaluation (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md
// section 3.2, 4.1). Deliberately decoupled from target_probability: the
// validation and arithmetic below never read or reject a target.
// ===================================================================

namespace {

// Stage "probability mode" validation helpers. Separate from both the
// Stage 8 (InfectorParameters) and Stage 9 (ECAiInputs) helpers only so the
// thrown message names ProbabilityInputs; the policy is identical
// (finiteness first, then the physical domain, and no input is ever
// clamped into range).
void require_probability_finite(double value, const char* field_name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string("ProbabilityInputs.") + field_name +
                                     " must be finite");
    }
}

void require_probability_positive(double value, const char* field_name) {
    require_probability_finite(value, field_name);
    if (!(value > 0.0)) {
        throw std::invalid_argument(std::string("ProbabilityInputs.") + field_name +
                                     " must be strictly positive");
    }
}

void require_probability_in_range(double value, double lo, double hi, const char* field_name) {
    require_probability_finite(value, field_name);
    if (value < lo || value > hi) {
        throw std::invalid_argument(std::string("ProbabilityInputs.") + field_name +
                                     " must lie within [" + std::to_string(lo) + ", " +
                                     std::to_string(hi) + "]");
    }
}

void require_probability_nonnegative(double value, const char* field_name) {
    require_probability_finite(value, field_name);
    if (!(value >= 0.0)) {
        throw std::invalid_argument(std::string("ProbabilityInputs.") + field_name +
                                     " must be nonnegative");
    }
}

// Every computed intermediate must stay finite: an overflow to +/-Inf (or a
// NaN) is reported rather than returned as a probability.
void require_finite_probability_intermediate(double value, const char* what) {
    if (!std::isfinite(value)) {
        throw std::runtime_error(std::string("probability ") + what +
                                  " is not finite (overflow or invalid intermediate)");
    }
}

// Full domain check for the seven physical inputs plus the fixed input
// ECAi, run before any arithmetic. Bounds mirror validate_ecai_inputs
// exactly (same seven physical fields; PROBABILITY_MODE_IMPLEMENTATION_PLAN.md
// section 4.1), minus target_probability (there is none) plus
// input_ecai_l_s_person (finite, >= 0).
void validate_probability_inputs(const ProbabilityInputs& inputs) {
    require_probability_positive(inputs.susceptible_pbr_m3h, "susceptible_pbr_m3h");
    require_probability_in_range(inputs.gamma_per_h, 0.42, 0.61, "gamma_per_h");
    require_probability_positive(inputs.lambda_per_h, "lambda_per_h");
    require_probability_positive(inputs.volume_m3, "volume_m3");
    if (inputs.occupants <= 0) {
        throw std::invalid_argument("ProbabilityInputs.occupants must be a positive count");
    }
    require_probability_in_range(inputs.mask_efficiency, 0.0, 1.0, "mask_efficiency");
    require_probability_positive(inputs.duration_h, "duration_h");
    require_probability_nonnegative(inputs.input_ecai_l_s_person, "input_ecai_l_s_person");
}

}  // namespace

ProbabilityResult compute_probability_from_qer_sum(const ProbabilityInputs& inputs,
                                                    double qer_sum_per_h) {
    validate_probability_inputs(inputs);
    require_qer_value(qer_sum_per_h, "Probability QER_sum");

    // Eq. 7 forward: input ECAi (L/s/person) -> total effective clean air
    // delivery (m^3/h) at the effective occupancy.
    const double eca_total_m3_per_h = kM3PerHourPerLitrePerSecond *
                                       static_cast<double>(inputs.occupants) *
                                       inputs.input_ecai_l_s_person;
    require_finite_probability_intermediate(eca_total_m3_per_h, "ECA_total");

    // Total removal: the fixed ECA plus the sampled surface/decay loss term.
    // Strictly positive given validated inputs (ECA_total >= 0, VOL > 0,
    // gamma + lambda > 0), so the division for Q below is never by zero.
    const double removal_m3_per_h =
        eca_total_m3_per_h + inputs.volume_m3 * (inputs.gamma_per_h + inputs.lambda_per_h);
    require_finite_probability_intermediate(removal_m3_per_h, "Lambda");

    // mask_factor = (1 - mu)^2 (Eq. 2), formed by one multiplication so the
    // factor is applied exactly once.
    const double one_minus_mu = 1.0 - inputs.mask_efficiency;
    const double mask_factor = one_minus_mu * one_minus_mu;

    const double dose_numerator =
        inputs.susceptible_pbr_m3h * inputs.duration_h * mask_factor * qer_sum_per_h;
    require_finite_probability_intermediate(dose_numerator, "dose_numerator");

    const double quanta_inhaled = dose_numerator / removal_m3_per_h;
    require_finite_probability_intermediate(quanta_inhaled, "Q");
    if (!(quanta_inhaled >= 0.0)) {
        // Impossible after the validation above (every factor of
        // dose_numerator and removal_m3_per_h is nonnegative/positive); a
        // negative Q here would indicate a defect elsewhere, not a
        // reachable input state, so it is reported rather than clamped.
        throw std::runtime_error("probability Q is negative (impossible after validation; "
                                  "indicates a defect)");
    }

    // Eq. 1 forward: P = 1 - exp(-Q), evaluated as -expm1(-Q) to retain
    // relative precision for the small Q the default target implies.
    const double p_raw = -std::expm1(-quanta_inhaled);
    require_finite_probability_intermediate(p_raw, "P_raw");

    ProbabilityResult result;
    result.eca_total_m3_per_h = eca_total_m3_per_h;
    result.removal_m3_per_h = removal_m3_per_h;
    result.mask_factor = mask_factor;
    result.qer_sum_per_h = qer_sum_per_h;
    result.dose_numerator = dose_numerator;
    result.quanta_inhaled = quanta_inhaled;
    // Sign normalization, mirroring compute_ecai_from_qer_sum's clamp: only
    // a strictly positive raw value passes through, so a zero-dose or
    // signed-zero result can never carry a negative sign bit forward.
    result.probability = (p_raw > 0.0) ? p_raw : 0.0;
    // P == 1.0 exactly is reachable once exp(-Q) underflows below 2^-53;
    // this is numerical saturation of -expm1(-Q) in double precision, not
    // an error, and is recorded rather than rejected.
    result.saturated = (result.probability == 1.0);
    return result;
}

ProbabilityResult compute_probability(const ProbabilityInputs& inputs,
                                       const std::vector<double>& qer_per_h) {
    // Same left-to-right, order-preserving accumulation policy as
    // compute_ecai above.
    double qer_sum_per_h = 0.0;
    for (const double qer : qer_per_h) {
        require_qer_value(qer, "Probability QER_i");
        qer_sum_per_h += qer;
    }
    if (!std::isfinite(qer_sum_per_h)) {
        throw std::runtime_error(
            "Probability QER_sum overflowed while summing the supplied QER values");
    }
    return compute_probability_from_qer_sum(inputs, qer_sum_per_h);
}

}  // namespace model
}  // namespace ecai_cpp2
