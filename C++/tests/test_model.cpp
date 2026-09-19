// Stage 8 tests for model.h/model.cpp: Jones QER (Eq. 3, 4, 5) as a pure
// deterministic function of already-sampled InfectorParameters, per
// docs/MODEL_SPECIFICATION.md sections 3, 5, 7, 10, 11.1 and
// IMPLEMENTATION_PLAN.md Section 7, Stage 8.
//
// Every expected literal below is computed independently in this file (by
// hand, from the frozen equations), never by re-deriving or copying
// model.cpp's own arithmetic, so this test genuinely exercises the
// production code against the specification rather than against itself.

#include "model.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <string>
#include <vector>

#ifndef ECAI_CPP2_MODEL_SOURCE_PATH
#error "ECAI_CPP2_MODEL_SOURCE_PATH must be defined by the build system"
#endif

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

bool approx(double actual, double expected, double abs_tol = 1e-12, double rel_tol = 1e-10) {
    if (!std::isfinite(actual) || !std::isfinite(expected)) return false;
    const double diff = std::fabs(actual - expected);
    if (diff <= abs_tol) return true;
    return diff <= rel_tol * std::fabs(expected);
}

using ecai_cpp2::model::compute_qer;
using ecai_cpp2::model::QERResult;
using ecai_cpp2::sampling_plan::InfectorParameters;

// Builds a valid, hand-pickable InfectorParameters with round numbers so the
// expected Vdrop/GVL_m3/VER/QER literals below can be derived by hand.
// Units follow sampling_plan.h's InfectorParameters contract:
//   pbr_m3h              : m^3/h
//   c_drop               : Table 3 printed x1e4 m^-3 units (e.g. 10 -> 1e5 m^-3)
//   d_bar                : Table 3 printed x1e-6 m units (e.g. 2.0 -> 2e-6 m)
//   evaporation_factor E : dimensionless, Table 1 hard bounds [2, 5]
//   gvl_copies_per_ml    : RNA copies mL^-1 (linear, already exponentiated)
//   viable_fraction VF   : dimensionless, Table 1 hard bounds [1e-4, 1e-2]
//   rtd                  : dimensionless, Table 1 hard bounds [0.43, 0.65]
//   dose_constant DK     : dimensionless, Table 1 hard bounds [5, 15]
InfectorParameters make_infector(double pbr_m3h, double c_drop, double d_bar,
                                  double evaporation_factor, double gvl_copies_per_ml,
                                  double viable_fraction, double rtd, double dose_constant) {
    InfectorParameters p;
    p.pbr_m3h = pbr_m3h;
    p.c_drop = c_drop;
    p.d_bar = d_bar;
    p.evaporation_factor = evaporation_factor;
    p.gvl_copies_per_ml = gvl_copies_per_ml;
    p.viable_fraction = viable_fraction;
    p.rtd = rtd;
    p.dose_constant = dose_constant;
    return p;
}

// A fully valid infector (every field inside its Table 1 domain), used as
// the shared base case for scaling and domain-boundary tests below.
InfectorParameters valid_baseline_infector() {
    return make_infector(/*pbr_m3h=*/0.55, /*c_drop=*/14.0, /*d_bar=*/1.9,
                          /*evaporation_factor=*/3.5, /*gvl_copies_per_ml=*/10.0,
                          /*viable_fraction=*/0.005, /*rtd=*/0.5, /*dose_constant=*/10.0);
}

bool throws_invalid_argument(const InfectorParameters& infector) {
    try {
        compute_qer(infector);
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

// One named field of InfectorParameters, addressed via a setter closure, so
// the NaN/+Inf/-Inf sweep below can iterate all 8 fields without repeating
// eight near-identical check() call sites.
struct FieldSpec {
    const char* name;
    std::function<void(InfectorParameters&, double)> set;
};

const std::vector<FieldSpec>& all_infector_fields() {
    static const std::vector<FieldSpec> fields = {
        {"pbr_m3h", [](InfectorParameters& p, double v) { p.pbr_m3h = v; }},
        {"c_drop", [](InfectorParameters& p, double v) { p.c_drop = v; }},
        {"d_bar", [](InfectorParameters& p, double v) { p.d_bar = v; }},
        {"evaporation_factor", [](InfectorParameters& p, double v) { p.evaporation_factor = v; }},
        {"gvl_copies_per_ml", [](InfectorParameters& p, double v) { p.gvl_copies_per_ml = v; }},
        {"viable_fraction", [](InfectorParameters& p, double v) { p.viable_fraction = v; }},
        {"rtd", [](InfectorParameters& p, double v) { p.rtd = v; }},
        {"dose_constant", [](InfectorParameters& p, double v) { p.dose_constant = v; }},
    };
    return fields;
}

// Strips // line comments and /* */ block comments from C++ source text so
// the compliance scan below cannot false-positive on a forbidden word
// appearing only in prose (e.g. a comment explaining what the file does
// NOT do), matching the intent (not the literal implementation) of the
// analogous scan in test_lhs.cpp.
std::string strip_cpp_comments(const std::string& source) {
    std::string out;
    out.reserve(source.size());
    enum class State { kCode, kLineComment, kBlockComment, kString, kChar };
    State state = State::kCode;
    for (std::size_t i = 0; i < source.size(); ++i) {
        const char c = source[i];
        const char next = (i + 1 < source.size()) ? source[i + 1] : '\0';
        switch (state) {
            case State::kCode:
                if (c == '/' && next == '/') {
                    state = State::kLineComment;
                    ++i;
                } else if (c == '/' && next == '*') {
                    state = State::kBlockComment;
                    ++i;
                } else if (c == '"') {
                    state = State::kString;
                    out.push_back(c);
                } else if (c == '\'') {
                    state = State::kChar;
                    out.push_back(c);
                } else {
                    out.push_back(c);
                }
                break;
            case State::kLineComment:
                if (c == '\n') {
                    state = State::kCode;
                    out.push_back(c);
                }
                break;
            case State::kBlockComment:
                if (c == '*' && next == '/') {
                    state = State::kCode;
                    ++i;
                }
                break;
            case State::kString:
                out.push_back(c);
                if (c == '\\' && i + 1 < source.size()) {
                    out.push_back(next);
                    ++i;
                } else if (c == '"') {
                    state = State::kCode;
                }
                break;
            case State::kChar:
                out.push_back(c);
                if (c == '\\' && i + 1 < source.size()) {
                    out.push_back(next);
                    ++i;
                } else if (c == '\'') {
                    state = State::kCode;
                }
                break;
        }
    }
    return out;
}

using ecai_cpp2::model::compute_ecai;
using ecai_cpp2::model::compute_ecai_from_qer_sum;
using ecai_cpp2::model::ECAiInputs;
using ecai_cpp2::model::ECAiResult;
using ecai_cpp2::model::compute_probability;
using ecai_cpp2::model::compute_probability_from_qer_sum;
using ecai_cpp2::model::ProbabilityInputs;
using ecai_cpp2::model::ProbabilityResult;

// The shared, fully valid Stage 9 fixture used by the blocks below:
//   PBR_s = 0.6 m^3/h, gamma = lambda = 0.5 1/h, VOL = 50 m^3, I0 = 10,
//   mu = 0, D = 1.0 h, target probability = 0.5 (so Q_target = ln 2).
// Every field sits inside the domain MODEL_SPECIFICATION.md section 5 /
// section 12 and IMPLEMENTATION_PLAN.md section 5.2 permit.
ECAiInputs baseline_ecai_inputs() {
    ECAiInputs in;
    in.susceptible_pbr_m3h = 0.6;
    in.gamma_per_h = 0.5;
    in.lambda_per_h = 0.5;
    in.volume_m3 = 50.0;
    in.occupants = 10;
    in.mask_efficiency = 0.0;
    in.duration_h = 1.0;
    in.target_probability = 0.5;
    return in;
}

bool ecai_throws_invalid_argument(const ECAiInputs& in, double qer_sum) {
    try {
        compute_ecai_from_qer_sum(in, qer_sum);
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

bool ecai_vector_throws_invalid_argument(const ECAiInputs& in, const std::vector<double>& qer) {
    try {
        compute_ecai(in, qer);
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

bool ecai_throws_runtime_error(const ECAiInputs& in, double qer_sum) {
    try {
        compute_ecai_from_qer_sum(in, qer_sum);
    } catch (const std::runtime_error&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

bool ecai_vector_throws_runtime_error(const ECAiInputs& in, const std::vector<double>& qer) {
    try {
        compute_ecai(in, qer);
    } catch (const std::runtime_error&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

bool probability_throws_invalid_argument(const ProbabilityInputs& in, double qer_sum) {
    try {
        compute_probability_from_qer_sum(in, qer_sum);
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

bool probability_vector_throws_invalid_argument(const ProbabilityInputs& in,
                                                 const std::vector<double>& qer) {
    try {
        compute_probability(in, qer);
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

bool probability_throws_runtime_error(const ProbabilityInputs& in, double qer_sum) {
    try {
        compute_probability_from_qer_sum(in, qer_sum);
    } catch (const std::runtime_error&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

bool probability_vector_throws_runtime_error(const ProbabilityInputs& in,
                                              const std::vector<double>& qer) {
    try {
        compute_probability(in, qer);
    } catch (const std::runtime_error&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

// Every double-valued ECAiInputs field, addressed via a setter closure, so
// the NaN/+Inf/-Inf sweep can cover all seven of them without seven
// near-identical call sites. The eighth input, `occupants`, is an integral
// count that cannot represent NaN or an infinity at all; it is swept
// separately for its own nonpositive rejection.
struct ECAiFieldSpec {
    const char* name;
    std::function<void(ECAiInputs&, double)> set;
};

const std::vector<ECAiFieldSpec>& all_ecai_double_fields() {
    static const std::vector<ECAiFieldSpec> fields = {
        {"susceptible_pbr_m3h", [](ECAiInputs& p, double v) { p.susceptible_pbr_m3h = v; }},
        {"gamma_per_h", [](ECAiInputs& p, double v) { p.gamma_per_h = v; }},
        {"lambda_per_h", [](ECAiInputs& p, double v) { p.lambda_per_h = v; }},
        {"volume_m3", [](ECAiInputs& p, double v) { p.volume_m3 = v; }},
        {"mask_efficiency", [](ECAiInputs& p, double v) { p.mask_efficiency = v; }},
        {"duration_h", [](ECAiInputs& p, double v) { p.duration_h = v; }},
        {"target_probability", [](ECAiInputs& p, double v) { p.target_probability = v; }},
    };
    return fields;
}

// True when every reported intermediate and both ECAi values are finite.
bool all_ecai_fields_finite(const ECAiResult& r) {
    return std::isfinite(r.q_target) && std::isfinite(r.mask_factor) &&
           std::isfinite(r.qer_sum_per_h) && std::isfinite(r.term1_m3_per_h) &&
           std::isfinite(r.term2_m3_per_h) && std::isfinite(r.t_ecai_raw_m3_per_h) &&
           std::isfinite(r.ecai_raw_l_s_person) && std::isfinite(r.ecai_l_s_person);
}

// Exactly +0.0 (never -0.0): the frozen clamp must produce a positive zero
// so downstream percentile/zero-counting cannot distinguish clamped rows by
// a sign bit.
bool is_exact_positive_zero(double v) { return v == 0.0 && !std::signbit(v); }

// The shared, fully valid probability-mode fixture: the same physical
// values as baseline_ecai_inputs() (PBR_s = 0.6 m^3/h, gamma = lambda = 0.5
// 1/h, VOL = 50 m^3, I0 = 10, mu = 0, D = 1.0 h), with a fixed input ECAi
// of 0 L/s/person in place of a target probability.
ProbabilityInputs baseline_probability_inputs() {
    ProbabilityInputs in;
    in.susceptible_pbr_m3h = 0.6;
    in.gamma_per_h = 0.5;
    in.lambda_per_h = 0.5;
    in.volume_m3 = 50.0;
    in.occupants = 10;
    in.mask_efficiency = 0.0;
    in.duration_h = 1.0;
    in.input_ecai_l_s_person = 0.0;
    return in;
}

// Builds a ProbabilityInputs sharing the same seven physical fields as an
// ECAiInputs fixture (everything but target_probability), at a
// caller-supplied fixed input ECAi.
ProbabilityInputs probability_inputs_from_ecai(const ECAiInputs& in, double input_ecai_l_s_person) {
    ProbabilityInputs p;
    p.susceptible_pbr_m3h = in.susceptible_pbr_m3h;
    p.gamma_per_h = in.gamma_per_h;
    p.lambda_per_h = in.lambda_per_h;
    p.volume_m3 = in.volume_m3;
    p.occupants = in.occupants;
    p.mask_efficiency = in.mask_efficiency;
    p.duration_h = in.duration_h;
    p.input_ecai_l_s_person = input_ecai_l_s_person;
    return p;
}

// Every double-valued ProbabilityInputs field, addressed via a setter
// closure, so the NaN/+Inf/-Inf sweep can cover all seven of them (the same
// six physical double fields as ECAiInputs, minus target_probability, plus
// input_ecai_l_s_person) without seven near-identical call sites.
// `occupants` is swept separately, like ECAiInputs.occupants.
struct ProbabilityFieldSpec {
    const char* name;
    std::function<void(ProbabilityInputs&, double)> set;
};

const std::vector<ProbabilityFieldSpec>& all_probability_double_fields() {
    static const std::vector<ProbabilityFieldSpec> fields = {
        {"susceptible_pbr_m3h",
         [](ProbabilityInputs& p, double v) { p.susceptible_pbr_m3h = v; }},
        {"gamma_per_h", [](ProbabilityInputs& p, double v) { p.gamma_per_h = v; }},
        {"lambda_per_h", [](ProbabilityInputs& p, double v) { p.lambda_per_h = v; }},
        {"volume_m3", [](ProbabilityInputs& p, double v) { p.volume_m3 = v; }},
        {"mask_efficiency", [](ProbabilityInputs& p, double v) { p.mask_efficiency = v; }},
        {"duration_h", [](ProbabilityInputs& p, double v) { p.duration_h = v; }},
        {"input_ecai_l_s_person",
         [](ProbabilityInputs& p, double v) { p.input_ecai_l_s_person = v; }},
    };
    return fields;
}

}  // namespace

int main() {
    const double kPi = 3.14159265358979323846;

    // ================================================================
    // (1) Known hand-calculated QER fixture (Eq. 3, 4, 5).
    //
    // Inputs (all round numbers, hand-traceable):
    //   PBR = 2.0 m^3/h
    //   C_drop (printed) = 10        -> real C_drop = 10 * 1e4 = 1e5 m^-3
    //   d_bar (printed)  = 2.0       -> real d_bar  = 2.0 * 1e-6 = 2e-6 m
    //   E = 3.0
    //   GVL = 1.0 copies/mL          -> GVL_m3 = 1.0 * 1e6 = 1e6 copies/m^3
    //   VF = 0.01
    //   RTD = 0.5
    //   DK = 10.0
    //
    // Hand algebra:
    //   hydrated_diameter = d_bar_m * E = 2e-6 * 3.0 = 6e-6 m
    //   Vdrop = (pi/6) * (6e-6)^3 * 1e5
    //         = (pi/6) * 2.16e-16 * 1e5
    //         = (pi/6) * 2.16e-11
    //   GVL_m3 = 1e6
    //   VER = PBR * Vdrop * GVL_m3 * VF * R(=1.0)
    //       = 2.0 * Vdrop * 1e6 * 0.01
    //       = 2.0e4 * Vdrop
    //   QER = RTD * VER / DK = 0.5 * VER / 10.0 = 0.05 * VER
    // ================================================================
    {
        const InfectorParameters infector =
            make_infector(2.0, 10.0, 2.0, 3.0, 1.0, 0.01, 0.5, 10.0);

        const double expected_hydrated_diameter_m = (2.0 * 1e-6) * 3.0;  // 6e-6 m
        const double expected_vdrop =
            (kPi / 6.0) * expected_hydrated_diameter_m * expected_hydrated_diameter_m *
            expected_hydrated_diameter_m * (10.0 * 1e4);
        const double expected_gvl_m3 = 1.0 * 1e6;
        const double expected_ver = 2.0 * expected_vdrop * expected_gvl_m3 * 0.01 * 1.0;
        const double expected_qer = 0.5 * expected_ver / 10.0;

        const QERResult result = compute_qer(infector);

        check(approx(result.vdrop, expected_vdrop), "hand fixture: Vdrop matches");
        check(approx(result.gvl_copies_per_m3, expected_gvl_m3), "hand fixture: GVL_m3 matches");
        check(approx(result.ver_virions_per_h, expected_ver), "hand fixture: VER matches");
        check(approx(result.qer_per_h, expected_qer), "hand fixture: QER_i matches");
        check(result.qer_per_h > 0.0 && std::isfinite(result.qer_per_h),
              "hand fixture: QER_i is finite and positive");
    }

    // ================================================================
    // (2) Independently pinned numeric literals for the same fixture as
    // (1). These four constants were computed once, externally (double-
    // precision arithmetic evaluated outside this file), and hard-coded
    // here as decimal literals -- they are NOT re-derived from the
    // (kPi/6.0)*... formula written above, so a bug that corrupts both the
    // production formula and the test's restated formula identically would
    // still be caught here.
    //   Vdrop    = 1.1309733552923255e-11
    //   GVL_m3   = 1000000.0
    //   VER      = 2.261946710584651e-07
    //   QER      = 1.1309733552923255e-08
    // ================================================================
    {
        const InfectorParameters infector =
            make_infector(2.0, 10.0, 2.0, 3.0, 1.0, 0.01, 0.5, 10.0);
        const QERResult result = compute_qer(infector);

        check(approx(result.vdrop, 1.1309733552923255e-11),
              "pinned literal: Vdrop matches independently computed constant");
        check(approx(result.gvl_copies_per_m3, 1000000.0),
              "pinned literal: GVL_m3 matches independently computed constant");
        check(approx(result.ver_virions_per_h, 2.261946710584651e-07),
              "pinned literal: VER matches independently computed constant");
        check(approx(result.qer_per_h, 1.1309733552923255e-08),
              "pinned literal: QER_i matches independently computed constant");
    }

    // ================================================================
    // (3) Scaling: doubling PBR doubles VER and QER; Vdrop is unaffected
    // (PBR does not appear in Eq. 5).
    // ================================================================
    {
        const InfectorParameters base = valid_baseline_infector();
        InfectorParameters doubled = base;
        doubled.pbr_m3h *= 2.0;

        const QERResult r_base = compute_qer(base);
        const QERResult r_doubled = compute_qer(doubled);

        check(approx(r_doubled.vdrop, r_base.vdrop), "doubling PBR: Vdrop unchanged");
        check(approx(r_doubled.ver_virions_per_h, 2.0 * r_base.ver_virions_per_h),
              "doubling PBR: VER doubles");
        check(approx(r_doubled.qer_per_h, 2.0 * r_base.qer_per_h), "doubling PBR: QER doubles");
    }

    // ================================================================
    // (4) Scaling: doubling C_drop doubles Vdrop, VER, and QER (C_drop is a
    // direct linear multiplicative factor in Eq. 5).
    // ================================================================
    {
        const InfectorParameters base = valid_baseline_infector();
        InfectorParameters doubled = base;
        doubled.c_drop *= 2.0;

        const QERResult r_base = compute_qer(base);
        const QERResult r_doubled = compute_qer(doubled);

        check(approx(r_doubled.vdrop, 2.0 * r_base.vdrop), "doubling C_drop: Vdrop doubles");
        check(approx(r_doubled.ver_virions_per_h, 2.0 * r_base.ver_virions_per_h),
              "doubling C_drop: VER doubles");
        check(approx(r_doubled.qer_per_h, 2.0 * r_base.qer_per_h),
              "doubling C_drop: QER doubles");
    }

    // ================================================================
    // (5) Scaling: doubling GVL doubles VER and QER; Vdrop is unaffected
    // (GVL does not appear in Eq. 5).
    // ================================================================
    {
        const InfectorParameters base = valid_baseline_infector();
        InfectorParameters doubled = base;
        doubled.gvl_copies_per_ml *= 2.0;

        const QERResult r_base = compute_qer(base);
        const QERResult r_doubled = compute_qer(doubled);

        check(approx(r_doubled.vdrop, r_base.vdrop), "doubling GVL: Vdrop unchanged");
        check(approx(r_doubled.gvl_copies_per_m3, 2.0 * r_base.gvl_copies_per_m3),
              "doubling GVL: GVL_m3 doubles");
        check(approx(r_doubled.ver_virions_per_h, 2.0 * r_base.ver_virions_per_h),
              "doubling GVL: VER doubles");
        check(approx(r_doubled.qer_per_h, 2.0 * r_base.qer_per_h), "doubling GVL: QER doubles");
    }

    // ================================================================
    // (6) Scaling: doubling VF doubles VER and QER; Vdrop is unaffected
    // (VF does not appear in Eq. 5). 0.005 -> 0.01 stays within Table 1's
    // Beta hard bounds [1e-4, 1e-2] (0.01 is the exact upper bound).
    // ================================================================
    {
        const InfectorParameters base = valid_baseline_infector();
        InfectorParameters doubled = base;
        doubled.viable_fraction *= 2.0;

        const QERResult r_base = compute_qer(base);
        const QERResult r_doubled = compute_qer(doubled);

        check(approx(r_doubled.vdrop, r_base.vdrop), "doubling VF: Vdrop unchanged");
        check(approx(r_doubled.ver_virions_per_h, 2.0 * r_base.ver_virions_per_h),
              "doubling VF: VER doubles");
        check(approx(r_doubled.qer_per_h, 2.0 * r_base.qer_per_h), "doubling VF: QER doubles");
    }

    // ================================================================
    // (7) Scaling: RTD is a direct linear multiplicative factor in Eq. 3
    // (QER = RTD * VER / DK), so QER scales exactly with RTD while VER/Vdrop
    // (which do not depend on RTD) are unaffected. Table 1's Uniform hard
    // bounds for RTD are [0.43, 0.65] -- a ratio of only ~1.512, too narrow
    // to admit an exact-doubling pair (2 * 0.43 = 0.86 > 0.65) -- so this
    // test uses the widest domain-valid ratio, the two hard-bound endpoints
    // themselves, instead of a literal factor of 2.
    // ================================================================
    {
        InfectorParameters low_rtd = valid_baseline_infector();
        low_rtd.rtd = 0.43;
        InfectorParameters high_rtd = valid_baseline_infector();
        high_rtd.rtd = 0.65;
        const double ratio = 0.65 / 0.43;

        const QERResult r_low = compute_qer(low_rtd);
        const QERResult r_high = compute_qer(high_rtd);

        check(approx(r_high.vdrop, r_low.vdrop), "RTD endpoints: Vdrop unchanged");
        check(approx(r_high.ver_virions_per_h, r_low.ver_virions_per_h),
              "RTD endpoints: VER unchanged");
        check(approx(r_high.qer_per_h, ratio * r_low.qer_per_h),
              "RTD endpoints: QER scales by RTD's hi/lo ratio (0.65/0.43)");
    }

    // ================================================================
    // (8) Scaling: doubling d_bar scales Vdrop, VER, and QER by 8, since
    // Eq. 5 cubes (d_bar * E). d_bar has no Table 1 upper hard bound
    // (lognormal support), so 1.9 -> 3.8 is unconstrained.
    // ================================================================
    {
        const InfectorParameters base = valid_baseline_infector();
        InfectorParameters doubled = base;
        doubled.d_bar *= 2.0;

        const QERResult r_base = compute_qer(base);
        const QERResult r_doubled = compute_qer(doubled);

        check(approx(r_doubled.vdrop, 8.0 * r_base.vdrop), "doubling d_bar: Vdrop scales by 8");
        check(approx(r_doubled.ver_virions_per_h, 8.0 * r_base.ver_virions_per_h),
              "doubling d_bar: VER scales by 8");
        check(approx(r_doubled.qer_per_h, 8.0 * r_base.qer_per_h),
              "doubling d_bar: QER scales by 8");
    }

    // ================================================================
    // (9) Scaling: doubling E scales Vdrop, VER, and QER by 8, symmetric
    // with (8) since Eq. 5 cubes (d_bar * E). Uses a dedicated base (E=2.0)
    // so the doubled value (4.0) stays within Table 1's Beta hard bounds
    // [2, 5] (doubling the shared baseline's E=3.5 would give 7.0, out of
    // range).
    // ================================================================
    {
        InfectorParameters base = valid_baseline_infector();
        base.evaporation_factor = 2.0;
        InfectorParameters doubled = base;
        doubled.evaporation_factor = 4.0;

        const QERResult r_base = compute_qer(base);
        const QERResult r_doubled = compute_qer(doubled);

        check(approx(r_doubled.vdrop, 8.0 * r_base.vdrop), "doubling E: Vdrop scales by 8");
        check(approx(r_doubled.ver_virions_per_h, 8.0 * r_base.ver_virions_per_h),
              "doubling E: VER scales by 8");
        check(approx(r_doubled.qer_per_h, 8.0 * r_base.qer_per_h), "doubling E: QER scales by 8");
    }

    // ================================================================
    // (10) Scaling: doubling DK halves QER only; Vdrop and VER do not
    // depend on DK (Eq. 3: QER = RTD * VER / DK). Uses a dedicated base
    // (DK=5.0) so the doubled value (10.0) stays within Table 1's Uniform
    // hard bounds [5, 15] (doubling the shared baseline's DK=10.0 would
    // give 20.0, out of range).
    // ================================================================
    {
        InfectorParameters base = valid_baseline_infector();
        base.dose_constant = 5.0;
        InfectorParameters doubled = base;
        doubled.dose_constant = 10.0;

        const QERResult r_base = compute_qer(base);
        const QERResult r_doubled = compute_qer(doubled);

        check(approx(r_doubled.vdrop, r_base.vdrop), "doubling DK: Vdrop unchanged");
        check(approx(r_doubled.ver_virions_per_h, r_base.ver_virions_per_h),
              "doubling DK: VER unchanged");
        check(approx(r_doubled.qer_per_h, 0.5 * r_base.qer_per_h), "doubling DK: QER halves");
    }

    // ================================================================
    // (11) Domain: the exact Table 1 hard-bound endpoints for E, VF, RTD,
    // DK are accepted (inclusive bounds), and every returned intermediate
    // is finite and positive at each endpoint.
    // ================================================================
    {
        struct BoundEndpoint {
            const char* field_name;
            std::function<InfectorParameters(double)> build;
            double lo;
            double hi;
        };
        const std::vector<BoundEndpoint> endpoints = {
            {"evaporation_factor",
             [](double v) {
                 InfectorParameters p = valid_baseline_infector();
                 p.evaporation_factor = v;
                 return p;
             },
             2.0, 5.0},
            {"viable_fraction",
             [](double v) {
                 InfectorParameters p = valid_baseline_infector();
                 p.viable_fraction = v;
                 return p;
             },
             1e-4, 1e-2},
            {"rtd",
             [](double v) {
                 InfectorParameters p = valid_baseline_infector();
                 p.rtd = v;
                 return p;
             },
             0.43, 0.65},
            {"dose_constant",
             [](double v) {
                 InfectorParameters p = valid_baseline_infector();
                 p.dose_constant = v;
                 return p;
             },
             5.0, 15.0},
        };

        for (const auto& e : endpoints) {
            for (double bound : {e.lo, e.hi}) {
                const InfectorParameters infector = e.build(bound);
                bool threw = false;
                QERResult result;
                try {
                    result = compute_qer(infector);
                } catch (const std::invalid_argument&) {
                    threw = true;
                }
                check(!threw, std::string(e.field_name) + " = " + std::to_string(bound) +
                                   " (exact hard bound): accepted, does not throw");
                if (!threw) {
                    check(std::isfinite(result.vdrop) && result.vdrop > 0.0 &&
                              std::isfinite(result.gvl_copies_per_m3) &&
                              result.gvl_copies_per_m3 > 0.0 &&
                              std::isfinite(result.ver_virions_per_h) &&
                              result.ver_virions_per_h > 0.0 && std::isfinite(result.qer_per_h) &&
                              result.qer_per_h > 0.0,
                          std::string(e.field_name) + " = " + std::to_string(bound) +
                              ": all intermediates finite and positive");
                }
            }
        }
    }

    // ================================================================
    // (12) Domain: values just below/above each Table 1 hard bound (E, VF,
    // RTD, DK) are rejected with std::invalid_argument.
    // ================================================================
    {
        struct RangeField {
            const char* name;
            double lo;
            double hi;
            std::function<void(InfectorParameters&, double)> set;
        };
        const std::vector<RangeField> ranged_fields = {
            {"evaporation_factor", 2.0, 5.0,
             [](InfectorParameters& p, double v) { p.evaporation_factor = v; }},
            {"viable_fraction", 1e-4, 1e-2,
             [](InfectorParameters& p, double v) { p.viable_fraction = v; }},
            {"rtd", 0.43, 0.65, [](InfectorParameters& p, double v) { p.rtd = v; }},
            {"dose_constant", 5.0, 15.0,
             [](InfectorParameters& p, double v) { p.dose_constant = v; }},
        };

        for (const auto& f : ranged_fields) {
            InfectorParameters below = valid_baseline_infector();
            f.set(below, f.lo - 1e-9);
            check(throws_invalid_argument(below),
                  std::string(f.name) + " just below its hard bound (" +
                      std::to_string(f.lo) + "): rejected");

            InfectorParameters above = valid_baseline_infector();
            f.set(above, f.hi + 1e-9);
            check(throws_invalid_argument(above),
                  std::string(f.name) + " just above its hard bound (" +
                      std::to_string(f.hi) + "): rejected");
        }
    }

    // ================================================================
    // (13) Domain: zero and negative values are rejected for every field
    // that is physically required to be strictly positive (the four
    // lognormal-derived, otherwise-unbounded fields: pbr_m3h, c_drop,
    // d_bar, gvl_copies_per_ml).
    // ================================================================
    {
        struct PositiveField {
            const char* name;
            std::function<void(InfectorParameters&, double)> set;
        };
        const std::vector<PositiveField> positive_fields = {
            {"pbr_m3h", [](InfectorParameters& p, double v) { p.pbr_m3h = v; }},
            {"c_drop", [](InfectorParameters& p, double v) { p.c_drop = v; }},
            {"d_bar", [](InfectorParameters& p, double v) { p.d_bar = v; }},
            {"gvl_copies_per_ml", [](InfectorParameters& p, double v) { p.gvl_copies_per_ml = v; }},
        };

        for (const auto& f : positive_fields) {
            InfectorParameters zero = valid_baseline_infector();
            f.set(zero, 0.0);
            check(throws_invalid_argument(zero), std::string(f.name) + " = 0.0: rejected");

            InfectorParameters negative = valid_baseline_infector();
            f.set(negative, -1.0);
            check(throws_invalid_argument(negative), std::string(f.name) + " = -1.0: rejected");
        }
    }

    // ================================================================
    // (14) Domain: NaN, +Infinity, and -Infinity are rejected for every one
    // of the 8 InfectorParameters fields.
    // ================================================================
    {
        const double nan_value = std::numeric_limits<double>::quiet_NaN();
        const double pos_inf = std::numeric_limits<double>::infinity();
        const double neg_inf = -std::numeric_limits<double>::infinity();
        const std::vector<std::pair<const char*, double>> nonfinite_values = {
            {"NaN", nan_value}, {"+Inf", pos_inf}, {"-Inf", neg_inf}};

        for (const auto& field : all_infector_fields()) {
            for (const auto& nv : nonfinite_values) {
                InfectorParameters infector = valid_baseline_infector();
                field.set(infector, nv.second);
                check(throws_invalid_argument(infector),
                      std::string(field.name) + " = " + nv.first + ": rejected");
            }
        }
    }

    // ================================================================
    // (15) Valid extreme boundary fixtures: every bounded field pinned to
    // one Table 1 hard-bound corner simultaneously (all-max, all-min),
    // combined with safely-scaled (no overflow/underflow) values for the
    // unbounded fields, still yields finite, positive intermediates.
    // ================================================================
    {
        const InfectorParameters all_max =
            make_infector(/*pbr_m3h=*/5.0, /*c_drop=*/1000.0, /*d_bar=*/100.0,
                          /*evaporation_factor=*/5.0, /*gvl_copies_per_ml=*/1e8,
                          /*viable_fraction=*/1e-2, /*rtd=*/0.65, /*dose_constant=*/15.0);
        const QERResult r_max = compute_qer(all_max);
        check(std::isfinite(r_max.vdrop) && r_max.vdrop > 0.0 &&
                  std::isfinite(r_max.gvl_copies_per_m3) && r_max.gvl_copies_per_m3 > 0.0 &&
                  std::isfinite(r_max.ver_virions_per_h) && r_max.ver_virions_per_h > 0.0 &&
                  std::isfinite(r_max.qer_per_h) && r_max.qer_per_h > 0.0,
              "extreme boundary fixture (all bounded fields at Table 1 max): all intermediates "
              "finite and positive");

        const InfectorParameters all_min =
            make_infector(/*pbr_m3h=*/0.001, /*c_drop=*/0.001, /*d_bar=*/0.001,
                          /*evaporation_factor=*/2.0, /*gvl_copies_per_ml=*/0.001,
                          /*viable_fraction=*/1e-4, /*rtd=*/0.43, /*dose_constant=*/5.0);
        const QERResult r_min = compute_qer(all_min);
        check(std::isfinite(r_min.vdrop) && r_min.vdrop > 0.0 &&
                  std::isfinite(r_min.gvl_copies_per_m3) && r_min.gvl_copies_per_m3 > 0.0 &&
                  std::isfinite(r_min.ver_virions_per_h) && r_min.ver_virions_per_h > 0.0 &&
                  std::isfinite(r_min.qer_per_h) && r_min.qer_per_h > 0.0,
              "extreme boundary fixture (all bounded fields at Table 1 min): all intermediates "
              "finite and positive");
    }

    // ================================================================
    // (16) Stage 9 tracer: one-infector hand-calculated ECAi fixture.
    //
    // Frozen inversion (MODEL_SPECIFICATION.md section 3, "Frozen ECAi
    // inversion used by C++2"; IMPLEMENTATION_PLAN.md section 4.4):
    //   Q_target    = -log(1 - target_probability)
    //   mask_factor = (1 - mu)^2
    //   term1_m3_h  = PBR_s * D * mask_factor * QER_sum / Q_target
    //   term2_m3_h  = VOL * (gamma + lambda)
    //   T_ECAi_raw  = term1_m3_h - term2_m3_h            [m^3/h]
    //   ECAi_raw    = T_ECAi_raw / (3.6 * I0)            [L/s/person]
    //   ECAi        = max(ECAi_raw, 0)
    //
    // Inputs (chosen so every intermediate is hand-traceable):
    //   PBR_s = 0.6 m^3/h, D = 1.0 h, mu = 0.0, VOL = 50 m^3,
    //   gamma = 0.5 1/h, lambda = 0.5 1/h, I0 = 10,
    //   target_probability = 0.5, one infector with QER = 100 quanta/h.
    //
    // Every expected value below is pinned as an independently computed
    // decimal literal (40-digit arbitrary-precision evaluation of the
    // equations above, printed to 30 significant digits), never as a
    // re-execution of model.cpp own arithmetic:
    //   Q_target = -ln(1 - 0.5) = ln 2 = 0.693147180559945309417232121458
    //   mask_factor = (1 - 0)^2 = 1
    //   QER_sum  = 100
    //   term1    = 0.6 * 1.0 * 1 * 100 / ln 2
    //            = 60 / ln 2 = 86.5617024533378044415954808601
    //   term2    = 50 * (0.5 + 0.5) = 50
    //   T_raw    = 36.5617024533378044415954808601
    //   ECAi_raw = T_raw / (3.6 * 10) = T_raw / 36
    //            = 1.01560284592605012337765224611  (positive: unclamped)
    // ================================================================
    {
        ecai_cpp2::model::ECAiInputs inputs;
        inputs.susceptible_pbr_m3h = 0.6;
        inputs.gamma_per_h = 0.5;
        inputs.lambda_per_h = 0.5;
        inputs.volume_m3 = 50.0;
        inputs.occupants = 10;
        inputs.mask_efficiency = 0.0;
        inputs.duration_h = 1.0;
        inputs.target_probability = 0.5;

        const ecai_cpp2::model::ECAiResult r =
            ecai_cpp2::model::compute_ecai(inputs, std::vector<double>{100.0});

        check(approx(r.q_target, 0.693147180559945309417232121458),
              "ECAi one-infector fixture: Q_target == ln 2");
        check(approx(r.mask_factor, 1.0),
              "ECAi one-infector fixture: mask_factor == 1 for mu == 0");
        check(approx(r.qer_sum_per_h, 100.0),
              "ECAi one-infector fixture: QER_sum == the single supplied QER");
        check(approx(r.term1_m3_per_h, 86.5617024533378044415954808601),
              "ECAi one-infector fixture: term1 == 60 / ln 2 m^3/h");
        check(approx(r.term2_m3_per_h, 50.0),
              "ECAi one-infector fixture: term2 == VOL * (gamma + lambda) == 50 m^3/h");
        check(approx(r.t_ecai_raw_m3_per_h, 36.5617024533378044415954808601),
              "ECAi one-infector fixture: T_ECAi_raw == term1 - term2 m^3/h");
        check(approx(r.ecai_raw_l_s_person, 1.01560284592605012337765224611),
              "ECAi one-infector fixture: unclamped ECAi_raw == T_raw / (3.6 * I0) L/s/person");
        check(approx(r.ecai_l_s_person, 1.01560284592605012337765224611),
              "ECAi one-infector fixture: positive ECAi is returned unchanged");
        check(!r.clamped, "ECAi one-infector fixture: a positive raw result is not clamped");
    }

    // ================================================================
    // (17) Stage 9: an explicitly supplied zero QER returns exactly +0.0.
    //
    // The pure algebra owns no parameter bank and no cursor, so it cannot
    // and does not prove that nothing was consumed -- that is
    // orchestration's job (Stage 7 assign_replicate, covered in
    // tests/test_sampling_plan.cpp). What is proven here is only the
    // algebraic half of the frozen rule: with QER_sum == 0, term1 == 0, so
    // T_ECAi_raw == -VOL * (gamma + lambda) <= 0 and the clamp gives
    // exactly +0.0 L/s/person.
    //
    // Pinned independently for the shared fixture: term2 = 50 * (0.5 + 0.5)
    // = 50 m^3/h, so the auditable raw value is
    //   ECAi_raw = -50 / (3.6 * 10) = -50/36
    //            = -1.38888888888888888888888888889 L/s/person.
    // ================================================================
    {
        const ECAiInputs in = baseline_ecai_inputs();

        const ECAiResult empty = compute_ecai(in, std::vector<double>{});
        const ECAiResult zeros = compute_ecai(in, std::vector<double>{0.0, 0.0, 0.0});
        const ECAiResult presummed = compute_ecai_from_qer_sum(in, 0.0);

        check(is_exact_positive_zero(empty.ecai_l_s_person),
              "ECAi zero-QER: an empty QER vector yields exactly +0.0 L/s/person");
        check(is_exact_positive_zero(zeros.ecai_l_s_person),
              "ECAi zero-QER: an explicit {0,0,0} QER vector yields exactly +0.0 L/s/person");
        check(is_exact_positive_zero(presummed.ecai_l_s_person),
              "ECAi zero-QER: a pre-summed QER_sum of 0 yields exactly +0.0 L/s/person");
        check(empty.qer_sum_per_h == 0.0 && zeros.qer_sum_per_h == 0.0,
              "ECAi zero-QER: QER_sum is exactly 0 for both the empty and the all-zero vector");
        check(empty.term1_m3_per_h == 0.0,
              "ECAi zero-QER: term1 is exactly 0 when QER_sum is 0");
        check(approx(empty.term2_m3_per_h, 50.0),
              "ECAi zero-QER: term2 is still the full VOL * (gamma + lambda) loss term");
        check(approx(empty.ecai_raw_l_s_person, -1.38888888888888888888888888889),
              "ECAi zero-QER: the auditable unclamped raw value is -50/36 L/s/person");
        check(empty.clamped, "ECAi zero-QER: the result is reported as clamped");

        // Property: for any admissible parameter set, an explicitly zero
        // QER_sum yields exactly +0.0. Five independent parameter sets, all
        // inside their frozen domains.
        bool zero_everywhere = true;
        const double gammas[5] = {0.42, 0.45, 0.50, 0.55, 0.61};
        const double lambdas[5] = {0.05, 0.52, 1.0, 3.0, 9.5};
        const double volumes[5] = {1.0, 33.0, 250.0, 1000.0, 12500.0};
        const double masks[5] = {0.0, 0.3, 0.5, 0.9, 1.0};
        const double targets[5] = {1e-9, 0.001, 0.01, 0.5, 0.99};
        const long long occupant_counts[5] = {1, 2, 10, 400, 5000};
        const double durations[5] = {0.25, 1.0, 2.0, 8.0, 24.0};
        for (int i = 0; i < 5; ++i) {
            ECAiInputs v;
            v.susceptible_pbr_m3h = 0.1 + 0.2 * i;
            v.gamma_per_h = gammas[i];
            v.lambda_per_h = lambdas[i];
            v.volume_m3 = volumes[i];
            v.occupants = occupant_counts[i];
            v.mask_efficiency = masks[i];
            v.duration_h = durations[i];
            v.target_probability = targets[i];
            const ECAiResult r = compute_ecai_from_qer_sum(v, 0.0);
            zero_everywhere = zero_everywhere && is_exact_positive_zero(r.ecai_l_s_person) &&
                               r.term1_m3_per_h == 0.0 && r.ecai_raw_l_s_person <= 0.0 &&
                               all_ecai_fields_finite(r);
        }
        check(zero_everywhere,
              "ECAi zero-QER: exactly +0.0 with a nonpositive auditable raw value for every "
              "admissible parameter set (5 sweeps)");
    }

    // ================================================================
    // (18) Stage 9: multiple-infector QER summation, in supplied order.
    //
    // Same shared fixture, three infectors with QER 40 + 60 + 100 = 200
    // quanta/h. Pinned independently (40-digit arbitrary-precision
    // evaluation of the frozen equations, printed to 30 significant
    // digits), not by restating model.cpp's arithmetic:
    //   term1 = 0.6 * 1.0 * 1 * 200 / ln 2 = 120 / ln 2
    //         = 173.12340490667560247633795409 m^3/h
    //   term2 = 50 m^3/h
    //   T_raw = 123.12340490667560247633795409 m^3/h
    //   ECAi  = T_raw / (3.6 * 10) = 3.42009458074098895767605428029
    // ================================================================
    {
        const ECAiInputs in = baseline_ecai_inputs();
        const ECAiResult r = compute_ecai(in, std::vector<double>{40.0, 60.0, 100.0});

        check(approx(r.qer_sum_per_h, 200.0),
              "ECAi multi-infector: QER_sum is the sum of the supplied per-infector values");
        check(approx(r.term1_m3_per_h, 173.12340490667560247633795409),
              "ECAi multi-infector: term1 == 120 / ln 2 m^3/h");
        check(approx(r.t_ecai_raw_m3_per_h, 123.12340490667560247633795409),
              "ECAi multi-infector: T_ECAi_raw == 123.12340490667560247634 m^3/h");
        check(approx(r.ecai_l_s_person, 3.42009458074098895767605428029),
              "ECAi multi-infector: ECAi == 3.42009458074098895768 L/s/person");
        check(!r.clamped && r.ecai_l_s_person == r.ecai_raw_l_s_person,
              "ECAi multi-infector: a positive raw ECAi is returned bit-for-bit unclamped");

        const ECAiResult presummed = compute_ecai_from_qer_sum(in, 200.0);
        check(presummed.ecai_l_s_person == r.ecai_l_s_person &&
                  presummed.term1_m3_per_h == r.term1_m3_per_h &&
                  presummed.t_ecai_raw_m3_per_h == r.t_ecai_raw_m3_per_h,
              "ECAi multi-infector: the vector overload matches the pre-summed overload "
              "bit-for-bit");

        const ECAiResult single = compute_ecai(in, std::vector<double>{200.0});
        check(single.ecai_l_s_person == r.ecai_l_s_person,
              "ECAi multi-infector: one infector carrying the whole sum gives the identical "
              "result (only the sum enters the equation)");

        const ECAiResult with_zero =
            compute_ecai(in, std::vector<double>{40.0, 0.0, 60.0, 100.0, 0.0});
        check(with_zero.ecai_l_s_person == r.ecai_l_s_person,
              "ECAi multi-infector: appending zero-QER entries leaves the result unchanged");

        // Supplied-order determinism: the sum is a plain left-to-right
        // accumulation, so a fixture whose reordering is not exactly
        // representable must reproduce its own supplied order bit-for-bit
        // on every call, and the same values summed left-to-right by the
        // test independently must match.
        const std::vector<double> ordered = {1e16, 1.0, -0.0 + 2.0, 3.0};
        double independent_sum = 0.0;
        for (double q : ordered) independent_sum += q;
        const ECAiResult ord_a = compute_ecai(in, ordered);
        const ECAiResult ord_b = compute_ecai(in, ordered);
        check(ord_a.qer_sum_per_h == independent_sum && ord_a.qer_sum_per_h == ord_b.qer_sum_per_h,
              "ECAi multi-infector: QER summation is a deterministic left-to-right sum in the "
              "supplied order");

        // Summation overflow is detected, never silently returned as an
        // infinite QER_sum.
        const double kMax = std::numeric_limits<double>::max();
        check(ecai_vector_throws_runtime_error(in, std::vector<double>{kMax, kMax}),
              "ECAi multi-infector: a QER vector whose sum overflows to infinity throws "
              "std::runtime_error");
    }

    // ================================================================
    // (19) Stage 9: negative raw ECAi clamps to exactly +0.0; a positive
    // raw value is returned unchanged (IMPLEMENTATION_PLAN.md section 4.4,
    // "ECAi = max(ECAi, 0)"; clamping is not optional in version 1, per
    // section 5.1).
    //
    // Shared fixture with QER_sum = 12.5 + 30.25 + 7.25 = 50, pinned:
    //   term1 = 0.6 * 50 / ln 2 = 30 / ln 2
    //         = 43.2808512266689006190844885226 m^3/h
    //   T_raw = -6.7191487733310993809155114774 m^3/h  (term2 = 50)
    //   ECAi_raw = T_raw / 36 = -0.186643021481419427247653096595
    //   ECAi     = exactly +0.0
    // ================================================================
    {
        const ECAiInputs in = baseline_ecai_inputs();
        const ECAiResult r = compute_ecai(in, std::vector<double>{12.5, 30.25, 7.25});

        check(approx(r.term1_m3_per_h, 43.2808512266689006190844885226),
              "ECAi clamp: term1 == 30 / ln 2 m^3/h for the negative fixture");
        check(approx(r.t_ecai_raw_m3_per_h, -6.7191487733310993809155114774),
              "ECAi clamp: T_ECAi_raw is negative when the loss term exceeds term1");
        check(approx(r.ecai_raw_l_s_person, -0.186643021481419427247653096595),
              "ECAi clamp: the unclamped raw ECAi is preserved for audit");
        check(is_exact_positive_zero(r.ecai_l_s_person),
              "ECAi clamp: a negative raw ECAi clamps to exactly +0.0, never -0.0");
        check(r.clamped, "ECAi clamp: the clamped flag is set for a negative raw ECAi");

        const ECAiResult positive = compute_ecai_from_qer_sum(in, 200.0);
        check(!positive.clamped && positive.ecai_l_s_person == positive.ecai_raw_l_s_person,
              "ECAi clamp: a positive raw ECAi is returned bit-for-bit unchanged and not "
              "flagged");

        // term1 == term2 exactly: the result is zero, nonnegative, and not
        // reported as clamped. QER_sum = VOL * (gamma + lambda) * Q_target
        // / (PBR_s * D) = 100 * 1.0 * ln 2 / 1.0 = 100 * ln 2.
        ECAiInputs exact = baseline_ecai_inputs();
        exact.susceptible_pbr_m3h = 1.0;
        exact.volume_m3 = 100.0;
        exact.gamma_per_h = 0.42;
        exact.lambda_per_h = 0.58;  // term2 = 100 * 1.0 = 100 m^3/h
        const ECAiResult at_zero =
            compute_ecai_from_qer_sum(exact, 69.3147180559945309417232121458);
        check(std::fabs(at_zero.ecai_l_s_person) < 1e-12 && at_zero.ecai_l_s_person >= 0.0,
              "ECAi clamp: term1 == term2 gives an ECAi of zero to within 1e-12, never "
              "negative");
    }

    // ================================================================
    // (20) Stage 9: mask factor (1 - mu)^2, applied exactly once (Eq. 2;
    // mu = 0.30 is the frozen publication Health Care value,
    // MODEL_SPECIFICATION.md section 12).
    //
    // Shared fixture, QER_sum = 200, mu = 0.30. Pinned independently:
    //   mask_factor = (1 - 0.3)^2 = 0.49
    //   term1 = 0.49 * 120 / ln 2 = 84.8304684042710479042838607088 m^3/h
    //   T_raw = 34.8304684042710479042838607088 m^3/h
    //   ECAi  = 0.967513011229751330674551686356 L/s/person
    // ================================================================
    {
        ECAiInputs in = baseline_ecai_inputs();
        in.mask_efficiency = 0.30;
        const ECAiResult r = compute_ecai_from_qer_sum(in, 200.0);

        check(approx(r.mask_factor, 0.49), "ECAi mask: mask_factor == (1 - 0.30)^2 == 0.49");
        check(approx(r.term1_m3_per_h, 84.8304684042710479042838607088),
              "ECAi mask: term1 == 0.49 * 120 / ln 2 m^3/h at the publication mu = 0.30");
        check(approx(r.ecai_l_s_person, 0.967513011229751330674551686356),
              "ECAi mask: ECAi == 0.96751301122975133067 L/s/person at mu = 0.30");

        // The factor is squared, not linear, and is applied exactly once:
        // the ratio to the unmasked term1 is exactly (1 - mu)^2.
        const ECAiResult unmasked = compute_ecai_from_qer_sum(baseline_ecai_inputs(), 200.0);
        check(approx(r.term1_m3_per_h / unmasked.term1_m3_per_h, 0.49),
              "ECAi mask: term1 scales by exactly (1 - mu)^2, not (1 - mu), and only once");
        check(approx(unmasked.mask_factor, 1.0),
              "ECAi mask: mu = 0 (the frozen general-category value) gives mask_factor 1");

        ECAiInputs half = baseline_ecai_inputs();
        half.mask_efficiency = 0.5;
        check(approx(compute_ecai_from_qer_sum(half, 200.0).mask_factor, 0.25),
              "ECAi mask: mu = 0.5 gives mask_factor 0.25");

        // mu = 1 (a perfect mask) zeroes term1, so the clamp gives +0.0.
        ECAiInputs perfect = baseline_ecai_inputs();
        perfect.mask_efficiency = 1.0;
        const ECAiResult r_perfect = compute_ecai_from_qer_sum(perfect, 200.0);
        check(r_perfect.mask_factor == 0.0 && r_perfect.term1_m3_per_h == 0.0 &&
                  is_exact_positive_zero(r_perfect.ecai_l_s_person),
              "ECAi mask: mu = 1.0 zeroes term1 and yields exactly +0.0 L/s/person");

        // Both endpoints of [0, 1] are admissible; anything outside is not
        // (IMPLEMENTATION_PLAN.md section 5.2).
        ECAiInputs below = baseline_ecai_inputs();
        below.mask_efficiency = std::nextafter(0.0, -1.0);
        ECAiInputs above = baseline_ecai_inputs();
        above.mask_efficiency = std::nextafter(1.0, 2.0);
        check(ecai_throws_invalid_argument(below, 200.0) &&
                  ecai_throws_invalid_argument(above, 200.0),
              "ECAi mask: the first double below 0 and the first above 1 are both rejected");
    }

    // ================================================================
    // (21) Stage 9: target-probability boundary validation and the
    // numerically stable Q_target = -log1p(-P(I)) evaluation.
    //
    // IMPLEMENTATION_PLAN.md section 5.2 rejects target probabilities
    // outside the open interval (0, 1). Inside it, Q_target must stay
    // accurate for the very small probabilities the standard uses (the
    // frozen default is 0.001): -log1p(-P) keeps full relative precision
    // where forming (1 - P) first would not.
    // ================================================================
    {
        const double kNaN = std::numeric_limits<double>::quiet_NaN();
        const double kInf = std::numeric_limits<double>::infinity();

        bool rejected = true;
        const double bad_targets[9] = {0.0,  -0.0, 1.0,  -1e-300, 1.0 + 1e-12,
                                       2.0,  kNaN, kInf, -kInf};
        for (double t : bad_targets) {
            ECAiInputs in = baseline_ecai_inputs();
            in.target_probability = t;
            rejected = rejected && ecai_throws_invalid_argument(in, 200.0);
        }
        check(rejected,
              "ECAi target probability: 0, -0, 1, negative, >1, NaN and +/-Inf are all "
              "rejected");

        // The frozen Jones default, pinned independently:
        //   -ln(1 - 0.001) = 0.00100050033358353352098050148502
        ECAiInputs jones = baseline_ecai_inputs();
        jones.target_probability = 0.001;
        check(approx(compute_ecai_from_qer_sum(jones, 200.0).q_target,
                     0.00100050033358353352098050148502),
              "ECAi target probability: Q_target at the frozen 0.001 default == "
              "0.0010005003335835335");

        // log1p discriminator: at P = 1e-17 the quantity (1 - P) rounds to
        // exactly 1.0 in double precision, so a naive -log(1 - P) would
        // return exactly 0 and then be divided by. -log1p(-P) returns
        // 1.0000000000000000e-17. Checked with a relative-only tolerance so
        // a zero Q_target cannot pass on the absolute tolerance.
        ECAiInputs tiny = baseline_ecai_inputs();
        tiny.target_probability = 1e-17;
        const ECAiResult r_tiny = compute_ecai_from_qer_sum(tiny, 200.0);
        check(approx(r_tiny.q_target, 1e-17, /*abs_tol=*/0.0, /*rel_tol=*/1e-10),
              "ECAi target probability: Q_target == 1e-17 at P = 1e-17 (log1p, not log(1-P))");
        check(all_ecai_fields_finite(r_tiny) && r_tiny.ecai_l_s_person > 0.0,
              "ECAi target probability: P = 1e-17 still gives a finite, positive ECAi");

        // Largest representable target probability strictly below 1:
        // 1 - 2^-53, for which Q_target = 53 * ln 2 =
        // 36.7368005696771013991133024373.
        ECAiInputs near_one = baseline_ecai_inputs();
        near_one.target_probability = std::nextafter(1.0, 0.0);
        const ECAiResult r_near_one = compute_ecai_from_qer_sum(near_one, 200.0);
        check(approx(r_near_one.q_target, 36.7368005696771013991133024373),
              "ECAi target probability: Q_target == 53 * ln 2 at the largest double below 1");
        check(all_ecai_fields_finite(r_near_one),
              "ECAi target probability: the largest double below 1 keeps every intermediate "
              "finite");

        // Smallest representable positive target probability (the denormal
        // 5e-324, nextafter(0, 1)): admissible, and with QER_sum == 0 there
        // is nothing to divide by it, so the result is still exactly +0.0
        // and never 0/0 -> NaN.
        ECAiInputs denormal = baseline_ecai_inputs();
        denormal.target_probability = std::nextafter(0.0, 1.0);
        const ECAiResult r_denormal = compute_ecai_from_qer_sum(denormal, 0.0);
        check(is_exact_positive_zero(r_denormal.ecai_l_s_person) &&
                  all_ecai_fields_finite(r_denormal),
              "ECAi target probability: the smallest positive double with QER_sum == 0 still "
              "yields exactly +0.0 and no NaN");

        // The same denormal target with a positive QER_sum overflows term1
        // to infinity; that must be reported, never returned as a nonfinite
        // ECAi.
        check(ecai_throws_runtime_error(denormal, 200.0),
              "ECAi overflow: a denormal target probability that overflows term1 throws "
              "std::runtime_error rather than returning a nonfinite ECAi");
        // A QER_sum of DBL_MAX alone does not overflow term1 for the
        // shared fixture (0.6 * DBL_MAX / ln 2 is still finite); it does
        // once the 24 h duration multiplier is applied.
        ECAiInputs long_day = baseline_ecai_inputs();
        long_day.duration_h = 24.0;
        check(ecai_throws_runtime_error(long_day, std::numeric_limits<double>::max()),
              "ECAi overflow: a QER_sum that overflows term1 throws std::runtime_error");
    }

    // ================================================================
    // (22) Stage 9: exposure duration scaling and the m^3/h -> L/s/person
    // unit conversion (Eq. 7's 3.6 * I0 divisor; MODEL_SPECIFICATION.md
    // section 10).
    //
    // Duration fixture: shared inputs, D = 2.0 h, QER_sum = 200. Pinned
    // independently:
    //   term1 = 0.6 * 2.0 * 200 / ln 2 = 240 / ln 2
    //         = 346.246809813351204952675908181 m^3/h
    //   T_raw = 296.246809813351204952675908181 m^3/h
    //   ECAi  = 8.22907805037086680424099744947 L/s/person
    //
    // Conversion fixture (chosen so T_ECAi_raw is exactly 36 m^3/h and the
    // answer is exactly 1 L/s/person at I0 = 10):
    //   P(I) = 0.6321205588285577 (= 1 - e^-1, so Q_target = 1),
    //   PBR_s = 0.5 m^3/h, D = 2.0 h, mu = 0, QER_sum = 86 quanta/h,
    //   VOL = 100 m^3, gamma = 0.45, lambda = 0.05 1/h, I0 = 10:
    //   term1 = 86, term2 = 50, T_raw = 36 m^3/h, ECAi = 36/36 = 1.
    // ================================================================
    {
        ECAiInputs two_hours = baseline_ecai_inputs();
        two_hours.duration_h = 2.0;
        const ECAiResult r2 = compute_ecai_from_qer_sum(two_hours, 200.0);
        const ECAiResult r1 = compute_ecai_from_qer_sum(baseline_ecai_inputs(), 200.0);

        check(approx(r2.term1_m3_per_h, 346.246809813351204952675908181),
              "ECAi duration: term1 == 240 / ln 2 m^3/h at D = 2 h");
        check(approx(r2.ecai_l_s_person, 8.22907805037086680424099744947),
              "ECAi duration: ECAi == 8.2290780503708668042 L/s/person at D = 2 h");
        check(approx(r2.term1_m3_per_h, 2.0 * r1.term1_m3_per_h),
              "ECAi duration: doubling D exactly doubles term1 (D enters term1 linearly)");
        check(approx(r2.term2_m3_per_h, r1.term2_m3_per_h),
              "ECAi duration: D does not enter the VOL * (gamma + lambda) loss term");

        bool duration_rejected = true;
        for (double bad : {0.0, -0.0, -1.0}) {
            ECAiInputs in = baseline_ecai_inputs();
            in.duration_h = bad;
            duration_rejected = duration_rejected && ecai_throws_invalid_argument(in, 200.0);
        }
        check(duration_rejected,
              "ECAi duration: a nonpositive duration is rejected; only a positive finite "
              "duration is admissible");

        ECAiInputs conv;
        conv.susceptible_pbr_m3h = 0.5;
        conv.gamma_per_h = 0.45;
        conv.lambda_per_h = 0.05;
        conv.volume_m3 = 100.0;
        conv.occupants = 10;
        conv.mask_efficiency = 0.0;
        conv.duration_h = 2.0;
        conv.target_probability = 0.6321205588285577;
        const ECAiResult rc = compute_ecai_from_qer_sum(conv, 86.0);
        check(approx(rc.q_target, 1.0), "ECAi conversion: Q_target == 1 at P(I) = 1 - e^-1");
        check(approx(rc.t_ecai_raw_m3_per_h, 36.0),
              "ECAi conversion: T_ECAi_raw == 36 m^3/h for the conversion fixture");
        check(approx(rc.ecai_l_s_person, 1.0),
              "ECAi conversion: 36 m^3/h shared by 10 people == 1.0 L/s/person");

        // 1 m^3/h == 1/3.6 L/s: with I0 = 1 the whole 36 m^3/h goes to one
        // person, giving 10 L/s/person.
        ECAiInputs solo = conv;
        solo.occupants = 1;
        check(approx(compute_ecai_from_qer_sum(solo, 86.0).ecai_l_s_person, 10.0),
              "ECAi conversion: the same 36 m^3/h at I0 = 1 is 10 L/s/person (1 m^3/h == "
              "1/3.6 L/s)");

        check(approx(r1.ecai_raw_l_s_person * 3.6 * 10.0, r1.t_ecai_raw_m3_per_h),
              "ECAi conversion: ECAi_raw * 3.6 * I0 == T_ECAi_raw");

        // I0 enters only the conversion: doubling it halves ECAi_raw while
        // leaving both m^3/h terms untouched.
        ECAiInputs twenty = baseline_ecai_inputs();
        twenty.occupants = 20;
        const ECAiResult r20 = compute_ecai_from_qer_sum(twenty, 200.0);
        check(approx(r20.t_ecai_raw_m3_per_h, r1.t_ecai_raw_m3_per_h) &&
                  approx(r20.ecai_raw_l_s_person, r1.ecai_raw_l_s_person / 2.0),
              "ECAi conversion: doubling I0 halves ECAi while leaving T_ECAi_raw unchanged");
    }

    // ================================================================
    // (23) Stage 9: finite-output requirement and invalid / nonfinite
    // input rejection over every physical domain
    // (IMPLEMENTATION_PLAN.md section 5.2 for VOL / I0 / mu / D / P(I);
    // MODEL_SPECIFICATION.md section 5 / Table 1 for the sampled loss
    // rates, using the same hard-bounds policy Stage 8 applies to
    // E/VF/RTD/DK).
    // ================================================================
    {
        const double kNaN = std::numeric_limits<double>::quiet_NaN();
        const double kInf = std::numeric_limits<double>::infinity();

        // Seven double fields x {NaN, +Inf, -Inf}. `occupants`, the eighth
        // input, is an integral count and cannot hold a nonfinite value at
        // all; it is covered by its own nonpositive sweep below.
        bool nonfinite_rejected = true;
        for (const auto& field : all_ecai_double_fields()) {
            for (double bad : {kNaN, kInf, -kInf}) {
                ECAiInputs in = baseline_ecai_inputs();
                field.set(in, bad);
                nonfinite_rejected =
                    nonfinite_rejected && ecai_throws_invalid_argument(in, 200.0);
            }
        }
        check(nonfinite_rejected,
              "ECAi validation: a NaN or +/-Inf in any of the seven double inputs is rejected");
        static_assert(std::is_integral<decltype(ECAiInputs{}.occupants)>::value,
                      "occupants must stay an integral count, which cannot represent NaN/Inf");

        // Strictly positive, unbounded-distribution inputs.
        bool positivity_rejected = true;
        for (double bad : {0.0, -0.0, -1e-300, -1.0}) {
            ECAiInputs pbr = baseline_ecai_inputs();
            pbr.susceptible_pbr_m3h = bad;
            ECAiInputs lam = baseline_ecai_inputs();
            lam.lambda_per_h = bad;
            ECAiInputs vol = baseline_ecai_inputs();
            vol.volume_m3 = bad;
            positivity_rejected =
                positivity_rejected && ecai_throws_invalid_argument(pbr, 200.0) &&
                ecai_throws_invalid_argument(lam, 200.0) &&
                ecai_throws_invalid_argument(vol, 200.0);
        }
        check(positivity_rejected,
              "ECAi validation: a nonpositive PBR_s, lambda, or VOL is rejected");

        // gamma keeps its Table 1 uniform hard bounds [0.42, 0.61].
        bool gamma_bounds_ok = true;
        for (double bad : {0.41, 0.0, 0.62, 5.0}) {
            ECAiInputs in = baseline_ecai_inputs();
            in.gamma_per_h = bad;
            gamma_bounds_ok = gamma_bounds_ok && ecai_throws_invalid_argument(in, 200.0);
        }
        for (double good : {0.42, 0.61}) {
            ECAiInputs in = baseline_ecai_inputs();
            in.gamma_per_h = good;
            gamma_bounds_ok = gamma_bounds_ok && !ecai_throws_invalid_argument(in, 200.0);
        }
        check(gamma_bounds_ok,
              "ECAi validation: gamma outside its Table 1 [0.42, 0.61] hard bounds is "
              "rejected, both endpoints accepted");

        // Occupancy must be a positive count.
        bool occupants_rejected = true;
        for (long long bad : {0LL, -1LL, -1000LL}) {
            ECAiInputs in = baseline_ecai_inputs();
            in.occupants = bad;
            occupants_rejected = occupants_rejected && ecai_throws_invalid_argument(in, 200.0);
        }
        check(occupants_rejected, "ECAi validation: a nonpositive occupant count is rejected");

        // Supplied QER values must be finite and nonnegative -- a negative
        // quanta emission rate is unphysical and must never silently
        // propagate into a smaller ECAi.
        bool qer_rejected =
            ecai_throws_invalid_argument(baseline_ecai_inputs(), -1e-12) &&
            ecai_throws_invalid_argument(baseline_ecai_inputs(), -200.0) &&
            ecai_throws_invalid_argument(baseline_ecai_inputs(), kNaN) &&
            ecai_throws_invalid_argument(baseline_ecai_inputs(), kInf) &&
            ecai_throws_invalid_argument(baseline_ecai_inputs(), -kInf) &&
            ecai_vector_throws_invalid_argument(baseline_ecai_inputs(),
                                                 std::vector<double>{40.0, -1.0, 60.0}) &&
            ecai_vector_throws_invalid_argument(baseline_ecai_inputs(),
                                                 std::vector<double>{40.0, kNaN}) &&
            ecai_vector_throws_invalid_argument(baseline_ecai_inputs(),
                                                 std::vector<double>{kInf}) &&
            ecai_vector_throws_invalid_argument(baseline_ecai_inputs(),
                                                 std::vector<double>{-kInf});
        check(qer_rejected,
              "ECAi validation: a negative or nonfinite QER value (per-infector or pre-summed) "
              "is rejected");

        // A -0.0 QER entry is a legitimate zero and must be accepted, and
        // must not turn the exact +0.0 clamp into -0.0.
        const ECAiResult negative_zero =
            compute_ecai(baseline_ecai_inputs(), std::vector<double>{-0.0, -0.0});
        check(is_exact_positive_zero(negative_zero.ecai_l_s_person),
              "ECAi validation: -0.0 QER entries are accepted and still clamp to exactly +0.0");

        // Finite output requirement across an admissible corner sweep: two
        // values per input dimension, all inside their frozen domains.
        const double pbrs[2] = {0.02, 5.0};
        const double gammas[2] = {0.42, 0.61};
        const double lambdas[2] = {1e-3, 20.0};
        const double vols[2] = {1.0, 1e5};
        const long long i0s[2] = {1, 100000};
        const double mus[2] = {0.0, 1.0};
        const double ds[2] = {1e-3, 24.0};
        const double ps[2] = {1e-6, 0.999};
        const double sums[2] = {0.0, 1e6};
        bool all_finite = true;
        bool never_negative = true;
        for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b)
                for (int c = 0; c < 2; ++c)
                    for (int e = 0; e < 2; ++e)
                        for (int f = 0; f < 2; ++f)
                            for (int g = 0; g < 2; ++g)
                                for (int h = 0; h < 2; ++h)
                                    for (int k = 0; k < 2; ++k)
                                        for (int m = 0; m < 2; ++m) {
                                            ECAiInputs in;
                                            in.susceptible_pbr_m3h = pbrs[a];
                                            in.gamma_per_h = gammas[b];
                                            in.lambda_per_h = lambdas[c];
                                            in.volume_m3 = vols[e];
                                            in.occupants = i0s[f];
                                            in.mask_efficiency = mus[g];
                                            in.duration_h = ds[h];
                                            in.target_probability = ps[k];
                                            const ECAiResult r =
                                                compute_ecai_from_qer_sum(in, sums[m]);
                                            all_finite = all_finite && all_ecai_fields_finite(r);
                                            never_negative =
                                                never_negative && r.ecai_l_s_person >= 0.0 &&
                                                !std::signbit(r.ecai_l_s_person);
                                        }
        check(all_finite,
              "ECAi finite output: every reported intermediate is finite across the 512-point "
              "admissible corner sweep");
        check(never_negative,
              "ECAi finite output: the clamped ECAi is never negative and never -0.0 across "
              "the same sweep");
    }

    // ================================================================
    // (24) Source-level compliance: model.cpp contains no RNG, no LHS
    // generation, no seed derivation, and no occupancy-category lookup --
    // it is a pure function of its InfectorParameters argument only (Stage
    // 8 objective; IMPLEMENTATION_PLAN.md Section 7, Stage 8). The source
    // is scanned with comments stripped first, so a forbidden word
    // appearing only in prose cannot cause a false positive.
    // ================================================================
    {
        std::ifstream source_file(ECAI_CPP2_MODEL_SOURCE_PATH);
        check(source_file.good(), "model.cpp source file is readable for compliance scan");
        std::stringstream buffer;
        buffer << source_file.rdbuf();
        const std::string code_only = strip_cpp_comments(buffer.str());

        const std::vector<std::string> forbidden_tokens = {
            "#include <random>",
            "#include \"lhs.h\"",
            "#include \"seed_manager.h\"",
            "#include \"category_data.h\"",
            "mt19937",
            "uniform_real_distribution",
            "std::shuffle",
            "generate_lhs",
            "SeedManager",
            "OccupancyCategory",
            "lookup_occupancy_category",
            "rand(",
            "srand(",
        };
        for (const auto& token : forbidden_tokens) {
            check(code_only.find(token) == std::string::npos,
                  "model.cpp does not contain forbidden token \"" + token + "\"");
        }

        // No global mutable state: no non-const, non-static-local variable
        // declared at namespace scope (a `static` local inside a function,
        // e.g. this test's own helpers, is not global state and is not
        // flagged here since it is scoped to strip_cpp_comments()'s own
        // translation unit, not model.cpp).
        check(code_only.find("\nstatic ") == std::string::npos &&
                  code_only.find("\nstatic\t") == std::string::npos,
              "model.cpp declares no namespace-scope `static` (mutable global) state");
    }

    // ================================================================
    // (25) Probability mode: round-trip identity
    // (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md section 3.3-1). For an
    // active row with T_ECAi > 0 (unclamped) and P_t in (0, 1), evaluating
    // the forward probability at ECAi_in = ECAi_raw (unclamped) must
    // return P == P_t to relative 1e-12. Exercised over several fixtures
    // already pinned above (mask factor, duration, unit conversion, the
    // frozen default target), so the round trip is checked against
    // independently-derived ECAi_raw values, not against a value this test
    // computes itself.
    // ================================================================
    {
        struct RoundTripCase {
            std::string label;
            ECAiInputs ecai_in;
            double qer_sum;
        };
        ECAiInputs masked = baseline_ecai_inputs();
        masked.mask_efficiency = 0.30;
        ECAiInputs two_hours = baseline_ecai_inputs();
        two_hours.duration_h = 2.0;
        ECAiInputs frozen_default = baseline_ecai_inputs();
        frozen_default.target_probability = 0.001;

        const std::vector<RoundTripCase> cases = {
            {"baseline (P_t = 0.5)", baseline_ecai_inputs(), 200.0},
            {"masked (mu = 0.30)", masked, 200.0},
            {"two-hour duration", two_hours, 200.0},
            {"frozen default target (P_t = 0.001)", frozen_default, 200.0},
            {"single infector", baseline_ecai_inputs(), 100.0},
        };

        for (const auto& c : cases) {
            const ECAiResult ecai_result = compute_ecai_from_qer_sum(c.ecai_in, c.qer_sum);
            check(!ecai_result.clamped && ecai_result.ecai_raw_l_s_person > 0.0,
                  "round trip [" + c.label + "]: precondition ECAi_raw > 0 (unclamped)");

            const ProbabilityInputs prob_in =
                probability_inputs_from_ecai(c.ecai_in, ecai_result.ecai_raw_l_s_person);
            const ProbabilityResult prob_result =
                compute_probability_from_qer_sum(prob_in, c.qer_sum);

            check(approx(prob_result.probability, c.ecai_in.target_probability,
                         /*abs_tol=*/0.0, /*rel_tol=*/1e-12),
                  "round trip [" + c.label + "]: P(ECAi_raw) == target_probability to relative "
                  "1e-12");
            check(!prob_result.saturated, "round trip [" + c.label + "]: not saturated");
        }
    }

    // ================================================================
    // (26) Probability mode: indicator identity (section 3.3-2). For any
    // X >= 0 and any row, P_i(X) > P_t  <=>  ECAi_raw_i > X, where
    // ECAi_raw_i is the row's own unclamped ECAi_raw at the same target
    // P_t. Six rows (distinct QER_sum, hence distinct ECAi_raw) at the
    // shared fixture's physical values and the frozen default target
    // P_t = 0.001; three threshold values X spanning below/within/above the
    // spread of ECAi_raw. Numerically the two sides may disagree by one
    // row when |ECAi_raw_i - X| <= 1e-9 * max(1, |X|); those near-tie rows
    // are excluded from the identity check, as the plan specifies.
    // ================================================================
    {
        ECAiInputs ecai_fixture = baseline_ecai_inputs();
        ecai_fixture.target_probability = 0.001;
        const double target_p = ecai_fixture.target_probability;

        const std::vector<double> qer_sums = {1.0, 5.0, 10.0, 50.0, 100.0, 500.0};
        std::vector<double> ecai_raws;
        for (double qer : qer_sums) {
            const ECAiResult r = compute_ecai_from_qer_sum(ecai_fixture, qer);
            check(!r.clamped, "indicator identity: fixture ECAi_raw is unclamped for QER_sum = " +
                                   std::to_string(qer));
            ecai_raws.push_back(r.ecai_raw_l_s_person);
        }

        const double min_raw = *std::min_element(ecai_raws.begin(), ecai_raws.end());
        const double max_raw = *std::max_element(ecai_raws.begin(), ecai_raws.end());
        const std::vector<double> thresholds = {std::max(0.0, min_raw - 1.0),
                                                  (min_raw + max_raw) / 2.0, max_raw + 1.0};

        for (double x : thresholds) {
            for (std::size_t i = 0; i < qer_sums.size(); ++i) {
                const double tol = 1e-9 * std::max(1.0, std::fabs(x));
                if (std::fabs(ecai_raws[i] - x) <= tol) continue;  // excluded near-tie band

                const ProbabilityInputs prob_in = probability_inputs_from_ecai(ecai_fixture, x);
                const ProbabilityResult prob_result =
                    compute_probability_from_qer_sum(prob_in, qer_sums[i]);

                const bool p_above_target = prob_result.probability > target_p;
                const bool ecai_above_x = ecai_raws[i] > x;
                check(p_above_target == ecai_above_x,
                      "indicator identity: P_i(X) > P_t iff ECAi_raw_i > X (X = " +
                          std::to_string(x) + ", row QER_sum = " + std::to_string(qer_sums[i]) +
                          ")");
            }
        }
    }

    // ================================================================
    // (27) Probability mode: order-statistic identity (section 3.3-3),
    // computed without depending on simulation.cpp's percentile_linear
    // (out of scope for this phase): N = 26 rows (distinct QER_sum, hence
    // distinct unclamped ECAi_raw with no tie), p = 0.96, so
    // rank = p * (N - 1) = 24.0 exactly. X is the row-wise ECAi_raw order
    // statistic at that exact rank. In probability mode with ECAi_in = X,
    // the same-rank order statistic of the row-wise probabilities equals
    // P_t: indicator identity (26) makes exactly 24 rows have
    // ECAi_raw < X (hence P < P_t), exactly 1 row have ECAi_raw > X (hence
    // P > P_t), and the row at X itself have P == P_t -- not because the
    // sort order of probabilities matches the sort order of ECAi_raw in
    // general (it need not, when PBR_s/gamma/lambda vary per row; the plan
    // requires the test not assert full order preservation).
    // ================================================================
    {
        constexpr int kN = 26;
        constexpr int kRank = 24;  // p * (N - 1) = 0.96 * 25 = 24.0 exactly
        ECAiInputs ecai_fixture = baseline_ecai_inputs();
        ecai_fixture.target_probability = 0.001;
        const double target_p = ecai_fixture.target_probability;

        std::vector<double> qer_sums;
        std::vector<double> ecai_raws;
        for (int i = 0; i < kN; ++i) {
            const double qer = 1.0 + static_cast<double>(i) * 3.0;  // 1, 4, 7, ..., distinct
            qer_sums.push_back(qer);
            const ECAiResult r = compute_ecai_from_qer_sum(ecai_fixture, qer);
            check(!r.clamped,
                  "order statistic: row " + std::to_string(i) + " ECAi_raw is unclamped");
            ecai_raws.push_back(r.ecai_raw_l_s_person);
        }

        std::vector<double> sorted_raws = ecai_raws;
        std::sort(sorted_raws.begin(), sorted_raws.end());
        check(sorted_raws[kRank - 1] != sorted_raws[kRank] &&
                  sorted_raws[kRank] != sorted_raws[kRank + 1],
              "order statistic: the rank-24 ECAi_raw order statistic is unique (no tie)");
        const double x = sorted_raws[kRank];

        std::vector<double> probabilities;
        for (int i = 0; i < kN; ++i) {
            const ProbabilityInputs prob_in = probability_inputs_from_ecai(ecai_fixture, x);
            const ProbabilityResult prob_result =
                compute_probability_from_qer_sum(prob_in, qer_sums[i]);
            probabilities.push_back(prob_result.probability);
        }
        std::sort(probabilities.begin(), probabilities.end());

        check(approx(probabilities[kRank], target_p, /*abs_tol=*/0.0, /*rel_tol=*/1e-10),
              "order statistic: sorted_P[24] == P_t to relative 1e-10 at ECAi_in = the rank-24 "
              "ECAi_raw order statistic");
    }

    // ================================================================
    // (28) Probability mode: saturation (section 3.2's numerical
    // saturation policy). Fixture engineered so quanta_inhaled Q equals
    // qer_sum exactly: input_ecai = 0 (so ECA_total = 0), PBR_s = D =
    // mask_factor = 1, and VOL * (gamma + lambda) = 1.
    // ================================================================
    {
        ProbabilityInputs sat = baseline_probability_inputs();
        sat.susceptible_pbr_m3h = 1.0;
        sat.gamma_per_h = 0.5;
        sat.lambda_per_h = 0.5;
        sat.volume_m3 = 1.0;
        sat.duration_h = 1.0;
        sat.mask_efficiency = 0.0;
        sat.input_ecai_l_s_person = 0.0;

        const ProbabilityResult q_1000 = compute_probability_from_qer_sum(sat, 1000.0);
        check(approx(q_1000.quanta_inhaled, 1000.0),
              "saturation: Q == qer_sum for the engineered fixture (Q = 1000)");
        check(q_1000.probability == 1.0 && q_1000.saturated,
              "saturation: Q = 1e3 gives P == 1.0 exactly and saturated == true");

        const ProbabilityResult q_50 = compute_probability_from_qer_sum(sat, 50.0);
        check(q_50.probability == 1.0 && q_50.saturated,
              "saturation: Q = 50 gives P == 1.0 exactly and saturated == true");

        const ProbabilityResult q_20 = compute_probability_from_qer_sum(sat, 20.0);
        check(q_20.probability < 1.0 && !q_20.saturated,
              "saturation: Q = 20 gives P < 1.0 and saturated == false");
    }

    // ================================================================
    // (29) Probability mode: overflow rejection (section 3.2, "Finiteness
    // (every intermediate, in order)"). A pathological input_ecai overflows
    // ECA_total; a pathological duration/QER_sum overflows dose_numerator;
    // a QER vector whose sum overflows is rejected before either fixed
    // point is reached. Every case throws std::runtime_error naming the
    // offending intermediate, never a nonfinite probability.
    // ================================================================
    {
        ProbabilityInputs overflow_eca = baseline_probability_inputs();
        overflow_eca.occupants = 1000;
        overflow_eca.input_ecai_l_s_person = 1e308;
        check(probability_throws_runtime_error(overflow_eca, 1.0),
              "overflow: ECA_total overflow (input_ecai = 1e308, I0 = 1000) throws "
              "std::runtime_error");
        try {
            compute_probability_from_qer_sum(overflow_eca, 1.0);
            check(false, "overflow: expected an exception for the ECA_total overflow fixture");
        } catch (const std::runtime_error& e) {
            const std::string message = e.what();
            check(message.find("ECA_total") != std::string::npos,
                  "overflow: the ECA_total overflow message names the offending intermediate");
        }

        ProbabilityInputs overflow_dose = baseline_probability_inputs();
        overflow_dose.duration_h = 24.0;
        check(probability_throws_runtime_error(overflow_dose,
                                                std::numeric_limits<double>::max()),
              "overflow: dose_numerator overflow (duration = 24h, QER_sum = DBL_MAX) throws "
              "std::runtime_error");
        try {
            compute_probability_from_qer_sum(overflow_dose, std::numeric_limits<double>::max());
            check(false, "overflow: expected an exception for the dose_numerator overflow "
                          "fixture");
        } catch (const std::runtime_error& e) {
            const std::string message = e.what();
            check(message.find("dose_numerator") != std::string::npos,
                  "overflow: the dose_numerator overflow message names the offending "
                  "intermediate");
        }

        check(probability_vector_throws_runtime_error(
                  baseline_probability_inputs(),
                  std::vector<double>{std::numeric_limits<double>::max(),
                                       std::numeric_limits<double>::max()}),
              "overflow: a QER vector whose sum overflows to infinity throws "
              "std::runtime_error");
    }

    // ================================================================
    // (30) Probability mode: invalid-input rejection (section 4.1, "the
    // seven checks of validate_ecai_inputs ... with 'ProbabilityInputs.'
    // message prefix", plus input_ecai_l_s_person finite and >= 0).
    // ProbabilityInputs has no target field to reject.
    // ================================================================
    {
        const double kNaN = std::numeric_limits<double>::quiet_NaN();
        const double kInf = std::numeric_limits<double>::infinity();

        bool nonfinite_rejected = true;
        for (const auto& field : all_probability_double_fields()) {
            for (double bad : {kNaN, kInf, -kInf}) {
                ProbabilityInputs in = baseline_probability_inputs();
                field.set(in, bad);
                nonfinite_rejected =
                    nonfinite_rejected && probability_throws_invalid_argument(in, 200.0);
            }
        }
        check(nonfinite_rejected,
              "probability validation: NaN or +/-Inf in any of the seven double inputs (incl. "
              "input_ecai_l_s_person) is rejected");

        bool positivity_rejected = true;
        for (double bad : {0.0, -0.0, -1e-300, -1.0}) {
            ProbabilityInputs pbr = baseline_probability_inputs();
            pbr.susceptible_pbr_m3h = bad;
            ProbabilityInputs lam = baseline_probability_inputs();
            lam.lambda_per_h = bad;
            ProbabilityInputs vol = baseline_probability_inputs();
            vol.volume_m3 = bad;
            positivity_rejected = positivity_rejected &&
                                   probability_throws_invalid_argument(pbr, 200.0) &&
                                   probability_throws_invalid_argument(lam, 200.0) &&
                                   probability_throws_invalid_argument(vol, 200.0);
        }
        check(positivity_rejected,
              "probability validation: a nonpositive PBR_s, lambda, or VOL is rejected");

        bool gamma_bounds_ok = true;
        for (double bad : {0.41, 0.0, 0.62, 5.0}) {
            ProbabilityInputs in = baseline_probability_inputs();
            in.gamma_per_h = bad;
            gamma_bounds_ok = gamma_bounds_ok && probability_throws_invalid_argument(in, 200.0);
        }
        for (double good : {0.42, 0.61}) {
            ProbabilityInputs in = baseline_probability_inputs();
            in.gamma_per_h = good;
            gamma_bounds_ok = gamma_bounds_ok && !probability_throws_invalid_argument(in, 200.0);
        }
        check(gamma_bounds_ok,
              "probability validation: gamma outside its Table 1 [0.42, 0.61] hard bounds is "
              "rejected, both endpoints accepted");

        bool occupants_rejected = true;
        for (long long bad : {0LL, -1LL, -1000LL}) {
            ProbabilityInputs in = baseline_probability_inputs();
            in.occupants = bad;
            occupants_rejected =
                occupants_rejected && probability_throws_invalid_argument(in, 200.0);
        }
        check(occupants_rejected, "probability validation: a nonpositive occupant count is "
                                   "rejected");

        bool input_ecai_rejected = true;
        for (double bad : {-1e-300, -1.0, -1000.0}) {
            ProbabilityInputs in = baseline_probability_inputs();
            in.input_ecai_l_s_person = bad;
            input_ecai_rejected =
                input_ecai_rejected && probability_throws_invalid_argument(in, 200.0);
        }
        check(input_ecai_rejected,
              "probability validation: a negative input_ecai_l_s_person is rejected");

        ProbabilityInputs zero_ecai = baseline_probability_inputs();
        zero_ecai.input_ecai_l_s_person = 0.0;
        check(!probability_throws_invalid_argument(zero_ecai, 200.0),
              "probability validation: input_ecai_l_s_person == 0.0 is accepted (boundary)");

        ProbabilityInputs mask_below = baseline_probability_inputs();
        mask_below.mask_efficiency = std::nextafter(0.0, -1.0);
        ProbabilityInputs mask_above = baseline_probability_inputs();
        mask_above.mask_efficiency = std::nextafter(1.0, 2.0);
        check(probability_throws_invalid_argument(mask_below, 200.0) &&
                  probability_throws_invalid_argument(mask_above, 200.0),
              "probability validation: the first double below 0 and the first above 1 for "
              "mask_efficiency are both rejected");

        bool qer_rejected =
            probability_throws_invalid_argument(baseline_probability_inputs(), -1e-12) &&
            probability_throws_invalid_argument(baseline_probability_inputs(), -200.0) &&
            probability_throws_invalid_argument(baseline_probability_inputs(), kNaN) &&
            probability_throws_invalid_argument(baseline_probability_inputs(), kInf) &&
            probability_throws_invalid_argument(baseline_probability_inputs(), -kInf) &&
            probability_vector_throws_invalid_argument(baseline_probability_inputs(),
                                                        std::vector<double>{40.0, -1.0, 60.0}) &&
            probability_vector_throws_invalid_argument(baseline_probability_inputs(),
                                                        std::vector<double>{40.0, kNaN});
        check(qer_rejected,
              "probability validation: a negative or nonfinite QER value (per-infector or "
              "pre-summed) is rejected");

        try {
            ProbabilityInputs bad_pbr = baseline_probability_inputs();
            bad_pbr.susceptible_pbr_m3h = -1.0;
            compute_probability_from_qer_sum(bad_pbr, 200.0);
            check(false, "probability validation: expected an exception for a negative PBR_s");
        } catch (const std::invalid_argument& e) {
            const std::string message = e.what();
            check(message.rfind("ProbabilityInputs.", 0) == 0,
                  "probability validation: the rejection message is prefixed "
                  "\"ProbabilityInputs.\"");
        }
    }

    // ================================================================
    // (31) Probability mode: zero dose (section 3.2, "P == 0 for an active
    // row occurs only when dose_num == 0"). An empty QER vector, an
    // explicit {0, 0, 0} QER vector, mu == 1 (perfect mask) with a positive
    // QER, and a pre-summed QER_sum of exactly 0 all yield P == 0.0
    // exactly, never -0.0.
    // ================================================================
    {
        ProbabilityInputs in = baseline_probability_inputs();
        in.input_ecai_l_s_person = 5.0;  // nonzero fixed ECA; removal/Q are otherwise unremarkable

        const ProbabilityResult empty = compute_probability(in, std::vector<double>{});
        const ProbabilityResult zeros = compute_probability(in, std::vector<double>{0.0, 0.0, 0.0});
        const ProbabilityResult presummed = compute_probability_from_qer_sum(in, 0.0);

        check(empty.probability == 0.0 && !std::signbit(empty.probability),
              "probability zero dose: an empty QER vector yields exactly +0.0");
        check(zeros.probability == 0.0 && !std::signbit(zeros.probability),
              "probability zero dose: an explicit {0,0,0} QER vector yields exactly +0.0");
        check(presummed.probability == 0.0 && !std::signbit(presummed.probability),
              "probability zero dose: a pre-summed QER_sum of 0 yields exactly +0.0");
        check(!empty.saturated && !zeros.saturated && !presummed.saturated,
              "probability zero dose: none of the zero-dose results are saturated");

        ProbabilityInputs perfect_mask = baseline_probability_inputs();
        perfect_mask.mask_efficiency = 1.0;
        const ProbabilityResult masked = compute_probability_from_qer_sum(perfect_mask, 200.0);
        check(masked.mask_factor == 0.0 && masked.dose_numerator == 0.0 &&
                  masked.probability == 0.0 && !std::signbit(masked.probability),
              "probability zero dose: mu = 1.0 (perfect mask) zeroes dose_numerator and yields "
              "exactly +0.0 with a positive QER");
    }

    // ================================================================
    // (32) Probability mode: signed zero. A -0.0 QER entry is a legitimate
    // zero (accepted, not rejected) and must not turn the exact +0.0
    // zero-dose result into -0.0; likewise Q == 0 must give P == +0.0,
    // mirroring compute_ecai_from_qer_sum's own sign-normalization clamp.
    // ================================================================
    {
        ProbabilityInputs in = baseline_probability_inputs();
        in.input_ecai_l_s_person = 5.0;

        const ProbabilityResult negative_zero_qer =
            compute_probability(in, std::vector<double>{-0.0, -0.0});
        check(negative_zero_qer.probability == 0.0 && !std::signbit(negative_zero_qer.probability),
              "probability signed zero: -0.0 QER entries are accepted and still yield exactly "
              "+0.0 probability");

        const ProbabilityResult zero_dose = compute_probability_from_qer_sum(in, 0.0);
        check(zero_dose.quanta_inhaled == 0.0 && !std::signbit(zero_dose.quanta_inhaled) &&
                  zero_dose.probability == 0.0 && !std::signbit(zero_dose.probability),
              "probability signed zero: Q == 0 gives P == +0.0 exactly (sign-normalized), never "
              "-0.0");
    }

    if (g_failures > 0) {
        std::cerr << g_failures << " test(s) failed.\n";
        return 1;
    }
    std::cout << "All tests passed.\n";
    return 0;
}
