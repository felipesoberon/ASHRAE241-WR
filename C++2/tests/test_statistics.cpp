// Stage 11 tests for statistics.h/statistics.cpp and the multi-replicate
// orchestration added to simulation.h/simulation.cpp, per
// IMPLEMENTATION_PLAN.md Section 7, Stage 11 and section 4.5, and
// docs/MODEL_SPECIFICATION.md section 11.5 (tolerances) and section 12
// (frozen rounding/percentile conventions).
//
// Every expected value below is either hand-derivable from the frozen
// formulas or an independently sourced literal (SciPy / published Student-t
// tables), never a value read back from this implementation.

#include "statistics.h"

#include "simulation.h"

#include <cmath>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef ECAI_CPP2_STATISTICS_HEADER_PATH
#error "ECAI_CPP2_STATISTICS_HEADER_PATH must be defined by the build system"
#endif
#ifndef ECAI_CPP2_STATISTICS_SOURCE_PATH
#error "ECAI_CPP2_STATISTICS_SOURCE_PATH must be defined by the build system"
#endif
#ifndef ECAI_CPP2_SIMULATION_HEADER_PATH
#error "ECAI_CPP2_SIMULATION_HEADER_PATH must be defined by the build system"
#endif

namespace {

int g_checks = 0;
int g_failures = 0;

void check(bool condition, const std::string& description) {
    ++g_checks;
    if (!condition) {
        std::cerr << "FAIL: " << description << "\n";
        ++g_failures;
    } else {
        std::cout << "PASS: " << description << "\n";
    }
}

// docs/MODEL_SPECIFICATION.md section 11.5: "Percentile interpolation / mean
// / sample SD / SE known-answer vectors | absolute 1e-9".
bool approx(double actual, double expected, double abs_tol = 1e-9) {
    if (!std::isfinite(actual) || !std::isfinite(expected)) return false;
    return std::fabs(actual - expected) <= abs_tol;
}

bool approx_relative(double actual, double expected, double rel_tol) {
    if (!std::isfinite(actual) || !std::isfinite(expected)) return false;
    return std::fabs(actual - expected) <= rel_tol * std::fabs(expected);
}

std::uint64_t positive_ulp_distance(double a, double b) {
    std::uint64_t a_bits = 0;
    std::uint64_t b_bits = 0;
    std::memcpy(&a_bits, &a, sizeof(double));
    std::memcpy(&b_bits, &b, sizeof(double));
    return a_bits > b_bits ? a_bits - b_bits : b_bits - a_bits;
}

// Bit-for-bit comparison, so a -0.0/+0.0 or last-ulp difference is caught.
bool bitwise_equal(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

bool bitwise_equal(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!bitwise_equal(a[i], b[i])) return false;
    }
    return true;
}

bool bitwise_equal(const std::optional<double>& a, const std::optional<double>& b) {
    if (a.has_value() != b.has_value()) return false;
    return !a || bitwise_equal(*a, *b);
}

template <typename Callable>
bool throws_invalid_argument(Callable&& fn) {
    try {
        fn();
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

template <typename Callable>
bool throws_runtime_error(Callable&& fn) {
    try {
        fn();
    } catch (const std::runtime_error&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

std::string read_file(const char* path) {
    std::ifstream in(path);
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

using ecai_cpp2::statistics::ceiling_to_5;
using ecai_cpp2::statistics::DispersionStatus;
using ecai_cpp2::statistics::kMaxSupportedStudentTDegreesOfFreedom;
using ecai_cpp2::statistics::ReplicateSummary;
using ecai_cpp2::statistics::student_t_two_sided_critical;
using ecai_cpp2::statistics::summarize_replicates;
using ecai_cpp2::OccupancyCategory;
using ecai_cpp2::lookup_occupancy_category;
using ecai_cpp2::simulation::EffectiveConfiguration;
using ecai_cpp2::simulation::EffectiveOverrides;
using ecai_cpp2::simulation::ReplicateBatch;
using ecai_cpp2::simulation::SimulationOptions;

EffectiveConfiguration make_batch_config(const OccupancyCategory& category, std::uint64_t n,
                                         std::uint64_t r,
                                         const EffectiveOverrides& overrides = {}) {
    EffectiveConfiguration config;
    config.category_name = category.canonical_name;
    config.i0 = overrides.occupants ? *overrides.occupants : category.i0;
    config.volume_m3 = overrides.volume_m3 ? *overrides.volume_m3 : category.volume_m3;
    config.mask_efficiency =
        overrides.mask_efficiency ? *overrides.mask_efficiency : category.mask_efficiency;
    config.community_infection_rate = overrides.community_infection_rate
                                          ? *overrides.community_infection_rate
                                          : category.community_infection_rate;
    config.target_probability = 0.001;
    config.duration_hours = 1.0;
    config.percentile = 0.96;
    config.confidence_level = 0.95;
    config.simulations_n = n;
    config.replicates_r = r;
    return config;
}

bool summaries_bitwise_equal(const ReplicateSummary& a, const ReplicateSummary& b) {
    return a.replicate_count == b.replicate_count &&
           bitwise_equal(a.replicate_estimates_l_s_person,
                         b.replicate_estimates_l_s_person) &&
           bitwise_equal(a.confidence_level, b.confidence_level) &&
           bitwise_equal(a.mean_l_s_person, b.mean_l_s_person) &&
           bitwise_equal(a.min_l_s_person, b.min_l_s_person) &&
           bitwise_equal(a.max_l_s_person, b.max_l_s_person) &&
           a.dispersion_status == b.dispersion_status &&
           bitwise_equal(a.degrees_of_freedom, b.degrees_of_freedom) &&
           bitwise_equal(a.sample_sd_l_s_person, b.sample_sd_l_s_person) &&
           bitwise_equal(a.standard_error_l_s_person, b.standard_error_l_s_person) &&
           bitwise_equal(a.t_critical, b.t_critical) &&
           bitwise_equal(a.ci_lower_l_s_person, b.ci_lower_l_s_person) &&
           bitwise_equal(a.ci_upper_l_s_person, b.ci_upper_l_s_person) &&
           bitwise_equal(a.adjusted_ceiling_to_5_l_s_person,
                         b.adjusted_ceiling_to_5_l_s_person);
}

// Known-answer vector, hand-derivable in exact arithmetic:
//   x        = {12, 15, 11, 18, 14}
//   sum      = 70, R = 5, mean = 14 exactly
//   devs     = {-2, 1, -3, 4, 0}, squares {4, 1, 9, 16, 0}, sum 30
//   variance = 30 / (5 - 1) = 7.5, SD = sqrt(7.5)
//   SE       = sqrt(7.5) / sqrt(5) = sqrt(1.5)
// Cross-checked against numpy (mean, std(ddof=1)).
void test_mean_sd_se_known_vector() {
    const std::vector<double> estimates = {12.0, 15.0, 11.0, 18.0, 14.0};
    const ReplicateSummary s = summarize_replicates(estimates, 0.95);

    check(s.replicate_count == 5, "replicate count is the supplied estimate count");
    check(approx(s.mean_l_s_person, 14.0), "raw mean matches the hand-derived value 14");
    check(s.sample_sd_l_s_person.has_value(), "sample SD is available for R = 5");
    check(s.standard_error_l_s_person.has_value(), "standard error is available for R = 5");
    check(s.sample_sd_l_s_person && approx(*s.sample_sd_l_s_person, 2.7386127875258306),
          "sample SD matches sqrt(30/4) = 2.7386127875258306");
    check(s.standard_error_l_s_person && approx(*s.standard_error_l_s_person, 1.224744871391589),
          "standard error matches sqrt(7.5)/sqrt(5) = 1.224744871391589");
}

using ecai_cpp2::statistics::student_t_two_sided_critical;

// Independently sourced two-sided Student-t critical values t_{1 - a/2, df}.
// Every entry below was produced by SciPy (scipy.stats.t.ppf(1 - (1 - conf)/2,
// df), SciPy 1.17.1 / NumPy 2.4.3) and the classical confidence levels are
// additionally corroborated by published Student-t tables (e.g. df = 1 / 0.95
// -> 12.706, df = 30 / 0.95 -> 2.042, df = 4 / 0.99 -> 4.604). They are
// literals here: nothing in this table is read back from the implementation.
// Tolerance: docs/MODEL_SPECIFICATION.md section 11.5, absolute 1e-6.
struct TCriticalCase {
    double confidence;
    double df;
    double expected;
};

const TCriticalCase kTCriticalCases[] = {
    {0.80, 1.0, 3.0776835371752544},   {0.80, 2.0, 1.8856180831641272},
    {0.80, 4.0, 1.533206274058944},    {0.80, 9.0, 1.3830287383966329},
    {0.80, 30.0, 1.3104150253913958},  {0.80, 499.0, 1.2832504230098853},
    {0.90, 1.0, 6.313751514675037},    {0.90, 2.0, 2.9199855803537242},
    {0.90, 4.0, 2.1318467863266495},   {0.90, 9.0, 1.8331129326562368},
    {0.90, 30.0, 1.6972608865939571},  {0.90, 499.0, 1.6479129840597124},
    {0.95, 1.0, 12.706204736174694},   {0.95, 2.0, 4.302652729749462},
    {0.95, 3.0, 3.1824463052837078},   {0.95, 4.0, 2.7764451051977934},
    {0.95, 9.0, 2.262157162798205},    {0.95, 10.0, 2.2281388519862744},
    {0.95, 30.0, 2.0422724563012378},  {0.95, 99.0, 1.9842169515864172},
    {0.95, 499.0, 1.9647293909876886}, {0.95, 4999.0, 1.9604386466615247},
    {0.98, 1.0, 31.820515953773935},   {0.98, 4.0, 3.746947387979196},
    {0.98, 9.0, 2.821437925025809},    {0.98, 30.0, 2.4572615424005915},
    {0.98, 499.0, 2.3338439936573887}, {0.99, 1.0, 63.656741162871526},
    {0.99, 2.0, 9.924843200918287},    {0.99, 4.0, 4.604094871349992},
    {0.99, 9.0, 3.2498355415921254},   {0.99, 30.0, 2.7499956535672254},
    {0.99, 499.0, 2.585717683111754},  {0.999, 1.0, 636.6192487687897},
    {0.999, 2.0, 31.599054576445365},  {0.999, 4.0, 8.610301581379522},
    {0.999, 9.0, 4.780912585931217},   {0.999, 30.0, 3.6459586350420623},
    {0.999, 499.0, 3.3101305703978303},
};

void test_student_t_critical_values() {
    double worst_absolute_error = 0.0;
    for (const TCriticalCase& c : kTCriticalCases) {
        const double actual = student_t_two_sided_critical(c.confidence, c.df);
        const double error = std::fabs(actual - c.expected);
        if (std::isfinite(error) && error > worst_absolute_error) worst_absolute_error = error;
        check(approx(actual, c.expected, 1e-6),
              "Student-t two-sided critical at confidence " + std::to_string(c.confidence) +
                  ", df " + std::to_string(static_cast<long long>(c.df)) + " matches " +
                  std::to_string(c.expected));
    }
    std::cout << "worst Student-t absolute error over "
              << (sizeof(kTCriticalCases) / sizeof(kTCriticalCases[0]))
              << " cases: " << worst_absolute_error << "\n";
    check(worst_absolute_error <= 1e-6,
          "worst Student-t absolute error over the whole grid is within the frozen 1e-6");

    // Monotonicity properties, checked independently of the table.
    check(student_t_two_sided_critical(0.95, 10.0) < student_t_two_sided_critical(0.99, 10.0),
          "a higher confidence level yields a larger critical value at fixed df");
    check(student_t_two_sided_critical(0.95, 10.0) > student_t_two_sided_critical(0.95, 100.0),
          "a larger df yields a smaller critical value at fixed confidence");
}

// For df = 1 the Student-t distribution is Cauchy, so the exact central
// probability inverse is t = tan(pi * confidence / 2). These literals were
// evaluated independently with Python mpmath at 100 decimal digits and then
// rounded to binary64; unlike (1 + confidence) / 2, this formulation does not
// erase a very small confidence level by addition to 1.
void test_student_t_low_confidence_df1() {
    const double at_1e20 = student_t_two_sided_critical(1e-20, 1.0);
    const double at_1e12 = student_t_two_sided_critical(1e-12, 1.0);
    check(std::isfinite(at_1e20) &&
              approx_relative(at_1e20, 1.5707963267948965e-20, 2e-15),
          "df = 1 at confidence 1e-20 converges to tan(pi*confidence/2)");
    check(std::isfinite(at_1e12) &&
              approx_relative(at_1e12, 1.5707963267948965e-12, 2e-15),
          "df = 1 at confidence 1e-12 converges without subtracting from 0.5");

    const double at_1e16 = student_t_two_sided_critical(1e-16, 1.0);
    const double at_1e8 = student_t_two_sided_critical(1e-8, 1.0);
    check(at_1e20 > 0.0 && at_1e20 < at_1e16 && at_1e16 < at_1e12 &&
              at_1e12 < at_1e8,
          "low-confidence critical values are positive and strictly monotone");
    std::ostringstream evidence;
    evidence.precision(17);
    evidence << "low-confidence df=1 actuals/errors: 1e-20=" << at_1e20 << " (rel "
             << std::fabs(at_1e20 - 1.5707963267948965e-20) /
                    1.5707963267948965e-20
             << "), 1e-12=" << at_1e12 << " (rel "
             << std::fabs(at_1e12 - 1.5707963267948965e-12) /
                    1.5707963267948965e-12
             << ")";
    std::cout << evidence.str() << "\n";
}

// References are C/(2*f(0)) evaluated with mpmath at 200 decimal digits and
// rounded to binary64. At these confidence levels the cubic inverse term is
// far below binary64 resolution. Subnormal results use an explicit ULP check;
// normal results use relative error so scaling does not hide a zero result.
void test_student_t_tiny_confidence_domain() {
    struct TinyCase {
        double confidence;
        double df;
        double expected;
    };
    const double denorm = std::numeric_limits<double>::denorm_min();
    const double minimum_normal = std::numeric_limits<double>::min();
    const TinyCase cases[] = {
        {denorm, 1.0, 0x0.0000000000002p-1022},
        {minimum_normal, 1.0, 0x1.921fb54442d18p-1022},
        {1e-300, 1.0, 0x1.0d4cab14b6bc0p-996},
        {1e-100, 1.0, 0x1.5fd17ae3bf80dp-332},
        {denorm, 2.0, 0x0.0000000000001p-1022},
        {minimum_normal, 2.0, 0x1.6a09e667f3bcdp-1022},
        {1e-300, 2.0, 0x1.e4e8d12762225p-997},
        {1e-100, 2.0, 0x1.3cbf6bd3ca8b4p-332},
        {denorm, 499.0, 0x0.0000000000001p-1022},
        {minimum_normal, 499.0, 0x1.41025b429aee4p-1022},
        {1e-300, 499.0, 0x1.adf48537255c5p-997},
        {1e-100, 499.0, 0x1.18d9de2f92242p-332},
        {denorm, 1e7, 0x0.0000000000001p-1022},
        {minimum_normal, 1e7, 0x1.40d93285f546ap-1022},
        {1e-300, 1e7, 0x1.adbd6475c6819p-997},
        {1e-100, 1e7, 0x1.18b5db97768e8p-332},
    };

    for (const TinyCase& c : cases) {
        const double actual = student_t_two_sided_critical(c.confidence, c.df);
        const bool subnormal_reference = c.expected < std::numeric_limits<double>::min();
        const bool accurate = subnormal_reference
                                  ? positive_ulp_distance(actual, c.expected) <= 1
                                  : approx_relative(actual, c.expected, 5e-15);
        check(std::isfinite(actual) && actual >= 0.0 && accurate,
              "tiny finite confidence returns the high-precision Student-t critical at df " +
                  std::to_string(static_cast<long long>(c.df)));
        std::cout.precision(17);
        std::cout << "tiny-confidence actual/reference: C=" << c.confidence
                  << ", df=" << c.df << ", actual=" << actual
                  << ", reference=" << c.expected << "\n";
    }
}

void test_student_t_supported_df_boundary() {
    const double at_lower_bound = student_t_two_sided_critical(0.95, 1.0);
    check(approx(at_lower_bound, 12.706204736174694, 1e-6),
          "the minimum supported df = 1 is accepted");
    check(throws_invalid_argument(
              [&] { student_t_two_sided_critical(0.95, std::nextafter(1.0, 0.0)); }),
          "one representable df below 1 is rejected before evaluation");
    check(throws_invalid_argument([&] { student_t_two_sided_critical(0.95, 0.5); }),
          "fractional df = 0.5 is outside the supported replicate-df domain");
    check(throws_invalid_argument([&] { student_t_two_sided_critical(0.95, 1.5); }),
          "non-integer df is outside the supported replicate-df domain");
    check(throws_invalid_argument([&] {
              student_t_two_sided_critical(0.95,
                                           std::numeric_limits<double>::denorm_min());
          }),
          "the smallest positive df is outside the supported replicate-df domain");

    // SciPy 1.17.1: scipy.stats.t.isf(0.025, 1e7).
    const double at_bound =
        student_t_two_sided_critical(0.95, kMaxSupportedStudentTDegreesOfFreedom);
    check(approx(at_bound, 1.9599642217672055, 1e-6),
          "the maximum supported df = 1e7 matches SciPy at 95% confidence");
    std::cout.precision(17);
    std::cout << "df=1e7 actual/absolute error: " << at_bound << ", "
              << std::fabs(at_bound - 1.9599642217672055) << "\n";
    check(throws_invalid_argument([&] {
              student_t_two_sided_critical(
                  0.95, std::nextafter(kMaxSupportedStudentTDegreesOfFreedom,
                                       std::numeric_limits<double>::infinity()));
          }),
          "one representable df above 1e7 is rejected before evaluation");
    check(throws_invalid_argument([&] { student_t_two_sided_critical(0.95, 1e20); }),
          "an inaccurate huge finite df is rejected rather than silently evaluated");
}

void test_student_t_branch_boundaries_are_monotone() {
    const double confidence_boundaries[] = {1e-8, 0.01, 0.1, 0.5};
    const double dfs[] = {1.0, 2.0, 499.0, 99999.0, 100000.0, 1e7};
    for (const double df : dfs) {
        for (const double boundary : confidence_boundaries) {
            const double below = student_t_two_sided_critical(
                std::nextafter(boundary, 0.0), df);
            const double at = student_t_two_sided_critical(boundary, df);
            const double above = student_t_two_sided_critical(
                std::nextafter(boundary, 1.0), df);
            check(below <= at && at <= above,
                  "Student-t critical is nondecreasing across adjacent confidence values at " +
                      std::to_string(boundary) + ", df " + std::to_string(df));
        }
    }

    const double df_pairs[][2] = {{99.0, 100.0}, {100.0, 101.0}, {2999.0, 3000.0},
                                  {9999.0, 10000.0}, {99999.0, 100000.0},
                                  {1000000.0, 10000000.0}};
    const double df_confidences[] = {
        std::nextafter(0.1, 1.0), 0.5, 0.8, 0.95};
    for (const auto& pair : df_pairs) {
        for (const double confidence : df_confidences) {
            const double lower_df = student_t_two_sided_critical(confidence, pair[0]);
            const double higher_df = student_t_two_sided_critical(confidence, pair[1]);
            check(lower_df >= higher_df,
                  "Student-t critical is nonincreasing across supported integer df values");
        }
    }
}

// Extreme-confidence behaviour: the critical value must stay finite and
// accurate, or fail loudly; it must never return a nonfinite or fabricated
// number. SciPy references use isf(0.5 * (1 - confidence), 1), matching the
// API's binary64 tail calculation exactly: 6366197.727026646 and
// 637129015469184.4. These are far outside the frozen
// 1e-6 absolute row of section 11.5, so they are checked relatively.
void test_student_t_extreme_confidence_and_validation() {
    const double very_high = student_t_two_sided_critical(0.9999999, 1.0);
    check(std::isfinite(very_high), "an extreme confidence level still yields a finite critical");
    check(approx_relative(very_high, 6366197.727026646, 1e-9),
          "t(1 - 1e-7 two-sided, df = 1) matches SciPy to a relative 1e-9");

    const double extreme = student_t_two_sided_critical(1.0 - 1e-15, 1.0);
    check(std::isfinite(extreme), "a 1 - 1e-15 confidence level still yields a finite critical");
    check(approx_relative(extreme, 637129015469184.4, 1e-9),
          "t(1 - 5e-16 two-sided, df = 1) matches SciPy to a relative 1e-9");

    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    const double inf_value = std::numeric_limits<double>::infinity();

    check(throws_invalid_argument([&] { student_t_two_sided_critical(nan_value, 4.0); }),
          "a nonfinite confidence level is rejected");
    check(throws_invalid_argument([&] { student_t_two_sided_critical(0.0, 4.0); }),
          "a confidence level of exactly 0 is rejected");
    check(throws_invalid_argument([&] { student_t_two_sided_critical(1.0, 4.0); }),
          "a confidence level of exactly 1 is rejected");
    check(throws_invalid_argument([&] { student_t_two_sided_critical(-0.5, 4.0); }),
          "a negative confidence level is rejected");
    check(throws_invalid_argument([&] { student_t_two_sided_critical(1.5, 4.0); }),
          "a confidence level above 1 is rejected");
    check(throws_invalid_argument([&] { student_t_two_sided_critical(0.95, 0.0); }),
          "df of exactly zero is rejected");
    check(throws_invalid_argument([&] { student_t_two_sided_critical(0.95, -4.0); }),
          "negative df is rejected");
    check(throws_invalid_argument([&] { student_t_two_sided_critical(0.95, nan_value); }),
          "nonfinite df is rejected");
    check(throws_invalid_argument([&] { student_t_two_sided_critical(0.95, inf_value); }),
          "infinite df is rejected");
}

// Full summary known-answer fixture: the same {12, 15, 11, 18, 14} vector,
// with CI bounds independently computed as mean +/- t * SE using the SciPy
// critical values above.
void test_confidence_interval_known_answers() {
    const std::vector<double> estimates = {12.0, 15.0, 11.0, 18.0, 14.0};

    const ReplicateSummary at95 = summarize_replicates(estimates, 0.95);
    check(at95.confidence_level == 0.95, "the summary echoes the confidence level it was given");
    check(at95.degrees_of_freedom && *at95.degrees_of_freedom == 4.0,
          "degrees of freedom are R - 1 = 4");
    check(at95.t_critical && approx(*at95.t_critical, 2.7764451051977934, 1e-6),
          "the summary's t critical at 95% / df 4 matches SciPy");
    check(at95.ci_lower_l_s_person && approx(*at95.ci_lower_l_s_person, 10.59956309670872),
          "95% CI lower bound matches mean - t * SE = 10.59956309670872");
    check(at95.ci_upper_l_s_person && approx(*at95.ci_upper_l_s_person, 17.40043690329128),
          "95% CI upper bound matches mean + t * SE = 17.40043690329128");
    check(at95.dispersion_status == DispersionStatus::available,
          "dispersion is reported available for R = 5");

    const ReplicateSummary at99 = summarize_replicates(estimates, 0.99);
    check(at99.t_critical && approx(*at99.t_critical, 4.604094871349992, 1e-6),
          "the summary's t critical at 99% / df 4 matches SciPy");
    check(at99.ci_lower_l_s_person && approx(*at99.ci_lower_l_s_person, 8.36115841891378),
          "99% CI lower bound matches mean - t * SE = 8.36115841891378");
    check(at99.ci_upper_l_s_person && approx(*at99.ci_upper_l_s_person, 19.63884158108622),
          "99% CI upper bound matches mean + t * SE = 19.63884158108622");
    check(at99.mean_l_s_person == at95.mean_l_s_person &&
              at99.sample_sd_l_s_person == at95.sample_sd_l_s_person,
          "changing only the confidence level leaves the mean and SD untouched");

    // Raw min/max, taken from the unsorted supplied vector.
    check(at95.min_l_s_person == 11.0, "raw minimum is the smallest supplied estimate");
    check(at95.max_l_s_person == 18.0, "raw maximum is the largest supplied estimate");

    // Diagnostics: the raw estimates are preserved verbatim, in the supplied
    // replicate-index order, so replicate-estimate shape can be assessed
    // without changing the estimator.
    check(bitwise_equal(at95.replicate_estimates_l_s_person, estimates),
          "the raw replicate estimates are preserved bit-for-bit in supplied order");

    // Deterministic supplied-order semantics: the same input reproduces a
    // bit-identical summary.
    const ReplicateSummary repeat = summarize_replicates(estimates, 0.95);
    check(bitwise_equal(repeat.mean_l_s_person, at95.mean_l_s_person) &&
              bitwise_equal(*repeat.sample_sd_l_s_person, *at95.sample_sd_l_s_person) &&
              bitwise_equal(*repeat.ci_lower_l_s_person, *at95.ci_lower_l_s_person),
          "summarize_replicates is a deterministic pure function of the supplied order");
}

// A zero-variance vector: SD and SE are exactly zero, and the CI collapses
// to the mean exactly (not approximately).
void test_zero_variance_vector() {
    const std::vector<double> estimates = {7.0, 7.0, 7.0, 7.0};
    const ReplicateSummary s = summarize_replicates(estimates, 0.95);
    check(s.sample_sd_l_s_person && *s.sample_sd_l_s_person == 0.0,
          "an identical-valued vector has exactly zero sample SD");
    check(s.standard_error_l_s_person && *s.standard_error_l_s_person == 0.0,
          "an identical-valued vector has exactly zero standard error");
    check(s.ci_lower_l_s_person && *s.ci_lower_l_s_person == 7.0,
          "a zero-variance CI lower bound is exactly the mean");
    check(s.ci_upper_l_s_person && *s.ci_upper_l_s_person == 7.0,
          "a zero-variance CI upper bound is exactly the mean");
    check(s.mean_l_s_person == 7.0 && s.min_l_s_person == 7.0 && s.max_l_s_person == 7.0,
          "a zero-variance vector's mean, min and max are all exactly the common value");
}

// R = 1: mean, min, max, the raw vector and the ceiling-to-5 comparator are
// available; SD, SE, t and the CI are typed-unavailable, never zero, NaN or
// fabricated. R = 0 is rejected outright.
void test_single_replicate_policy() {
    const std::vector<double> one = {13.25};
    const ReplicateSummary s = summarize_replicates(one, 0.95);

    check(s.replicate_count == 1, "R = 1 is accepted and reported as a count of 1");
    check(s.mean_l_s_person == 13.25, "R = 1 mean is exactly the single estimate");
    check(s.min_l_s_person == 13.25 && s.max_l_s_person == 13.25,
          "R = 1 min and max are exactly the single estimate");
    check(bitwise_equal(s.replicate_estimates_l_s_person, one),
          "R = 1 preserves the raw estimate vector");
    check(s.adjusted_ceiling_to_5_l_s_person == 15.0,
          "R = 1 still produces the ceiling-to-5 comparator 5 * ceil(13.25 / 5) = 15");
    check(s.dispersion_status == DispersionStatus::unavailable_single_replicate,
          "R = 1 reports the single-replicate dispersion status explicitly");
    check(!s.sample_sd_l_s_person.has_value(), "R = 1 sample SD is unavailable, not zero");
    check(!s.standard_error_l_s_person.has_value(), "R = 1 standard error is unavailable");
    check(!s.t_critical.has_value(), "R = 1 t critical is unavailable");
    check(!s.ci_lower_l_s_person.has_value(), "R = 1 CI lower bound is unavailable");
    check(!s.ci_upper_l_s_person.has_value(), "R = 1 CI upper bound is unavailable");
    check(!s.degrees_of_freedom.has_value(), "R = 1 degrees of freedom are unavailable");

    // Confidence is still validated at R = 1, even though no CI is produced.
    check(throws_invalid_argument([&] { summarize_replicates(one, 0.0); }),
          "R = 1 still rejects a confidence level of exactly 0");
    check(throws_invalid_argument([&] { summarize_replicates(one, 1.0); }),
          "R = 1 still rejects a confidence level of exactly 1");
    check(throws_invalid_argument(
              [&] { summarize_replicates(one, std::numeric_limits<double>::quiet_NaN()); }),
          "R = 1 still rejects a nonfinite confidence level");

    check(throws_invalid_argument([&] { summarize_replicates(std::vector<double>{}, 0.95); }),
          "R = 0 is rejected outright");
}

void test_invalid_inputs_rejected() {
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    const double inf_value = std::numeric_limits<double>::infinity();

    check(throws_invalid_argument([&] { summarize_replicates({1.0, nan_value, 3.0}, 0.95); }),
          "a NaN replicate estimate is rejected");
    check(throws_invalid_argument([&] { summarize_replicates({1.0, inf_value, 3.0}, 0.95); }),
          "an infinite replicate estimate is rejected");
    check(throws_invalid_argument([&] { summarize_replicates({1.0, -inf_value}, 0.95); }),
          "a negative-infinite replicate estimate is rejected");
    check(throws_invalid_argument([&] { summarize_replicates({1.0, -0.5, 3.0}, 0.95); }),
          "a negative replicate estimate is rejected (ECAi is clamped at zero)");
    check(throws_invalid_argument([&] { summarize_replicates({1.0, 2.0}, nan_value); }),
          "a nonfinite confidence level is rejected");
    check(throws_invalid_argument([&] { summarize_replicates({1.0, 2.0}, -0.1); }),
          "a negative confidence level is rejected");
    check(throws_invalid_argument([&] { summarize_replicates({1.0, 2.0}, 1.0); }),
          "a confidence level of exactly 1 is rejected");

    // Exactly zero is a legal ECAi estimate (an all-zero-infector replicate).
    const ReplicateSummary zeros = summarize_replicates({0.0, 0.0, 0.0}, 0.95);
    check(zeros.mean_l_s_person == 0.0 && zeros.adjusted_ceiling_to_5_l_s_person == 0.0,
          "an all-zero estimate vector is legal and its ceiling-to-5 comparator is 0");
}

// The literals below are first rounded to their actual binary64 values. An
// independent Python fractions.Fraction calculation over those exact doubles
// gives mean 1.000000000000002015984...e150, sample SD
// 1.580978652745069862034...e135, and SE 7.070351477028021408...e134.
void test_high_offset_statistics() {
    const std::vector<double> estimates = {
        1e150, 1e150 + 1e135, 1e150 + 2e135, 1e150 + 3e135, 1e150 + 4e135};
    const ReplicateSummary s = summarize_replicates(estimates, 0.95);
    check(bitwise_equal(s.mean_l_s_person, 1.000000000000002e150),
          "high-offset mean rounds from the exact representable inputs correctly");
    check(s.sample_sd_l_s_person &&
              approx_relative(*s.sample_sd_l_s_person, 1.5809786527450697e135, 2e-15),
          "high-offset sample SD matches exact-rational input arithmetic");
    check(s.standard_error_l_s_person &&
              approx_relative(*s.standard_error_l_s_person, 7.070351477028021e134, 2e-15),
          "high-offset standard error matches exact-rational input arithmetic");
    if (s.sample_sd_l_s_person && s.standard_error_l_s_person) {
        std::ostringstream evidence;
        evidence.precision(17);
        evidence << "high-offset actual mean/SD/SE: " << s.mean_l_s_person << ", "
                 << *s.sample_sd_l_s_person << ", " << *s.standard_error_l_s_person
                 << "; SD relative error: "
                 << std::fabs(*s.sample_sd_l_s_person - 1.5809786527450697e135) /
                        1.5809786527450697e135;
        std::cout << evidence.str() << "\n";
    }
}

// Overflow protection. A finite input whose squared deviations overflow must
// be rejected loudly rather than reported as an infinite or NaN SD.
void test_overflow_rejected() {
    const double huge = 1e308;
    check(throws_runtime_error([&] { summarize_replicates({0.0, huge}, 0.95); }),
          "a variance that overflows to infinity is rejected with a runtime error");

    // Welford keeps the mean of large equal values finite where a naive
    // running sum would overflow: this must succeed, not throw.
    const ReplicateSummary large = summarize_replicates({huge, huge, huge}, 0.95);
    check(large.mean_l_s_person == huge,
          "a stable mean of three 1e308 values is exactly 1e308 (no summation overflow)");
    check(large.sample_sd_l_s_person && *large.sample_sd_l_s_person == 0.0,
          "the stable variance of three identical large values is exactly zero");
}

void test_confidence_interval_overflow_rejected_directly() {
    const double maximum = std::numeric_limits<double>::max();
    check(throws_runtime_error([&] {
              ecai_cpp2::statistics::detail::checked_confidence_interval(0.0, maximum, 2.0);
          }),
          "a directly overflowing confidence-interval half-width is rejected loudly");
    check(throws_runtime_error([&] {
              ecai_cpp2::statistics::detail::checked_confidence_interval(maximum,
                                                                          maximum / 4.0, 1.0);
          }),
          "a finite half-width whose upper confidence bound overflows is rejected loudly");
}

// The frozen ceiling-to-5 comparator, D6 / section 12: 5 * ceil(raw / 5).
void test_ceiling_to_5_boundaries() {
    check(ceiling_to_5(25.0) == 25.0, "an exact multiple of 5 is left unchanged");
    check(ceiling_to_5(std::nextafter(25.0, std::numeric_limits<double>::infinity())) == 30.0,
          "one ulp above an exact multiple of 5 rounds up to the next multiple");
    check(ceiling_to_5(std::nextafter(25.0, 0.0)) == 25.0,
          "one ulp below an exact multiple of 5 rounds up to that multiple");
    check(ceiling_to_5(0.0) == 0.0, "zero maps to zero, not to 5");
    check(ceiling_to_5(std::numeric_limits<double>::denorm_min()) == 5.0,
          "the smallest positive value maps up to 5");
    check(ceiling_to_5(20.000001) == 25.0, "a value just above 20 maps to 25");
    check(ceiling_to_5(19.999999) == 20.0, "a value just below 20 maps to 20");
    const double dbl_max = std::numeric_limits<double>::max();
    check(std::isfinite(ceiling_to_5(dbl_max)) && ceiling_to_5(dbl_max) == dbl_max,
          "the largest finite double survives the ceiling adjustment without overflowing");
    check(throws_invalid_argument([&] { ceiling_to_5(std::numeric_limits<double>::infinity()); }),
          "a nonfinite raw value is rejected by the ceiling adjustment");

    // Boundary cases routed through a single-estimate summary, so the mean is
    // exactly the boundary value.
    check(summarize_replicates({25.0}, 0.95).adjusted_ceiling_to_5_l_s_person == 25.0,
          "a mean of exactly 25 yields an adjusted 25");
    check(summarize_replicates({std::nextafter(25.0, 1e300)}, 0.95)
                  .adjusted_ceiling_to_5_l_s_person == 30.0,
          "a mean one ulp above 25 yields an adjusted 30");
    check(summarize_replicates({std::nextafter(25.0, 0.0)}, 0.95)
                  .adjusted_ceiling_to_5_l_s_person == 25.0,
          "a mean one ulp below 25 yields an adjusted 25");
}

// No adjusted or rounded value may enter the mean, SD, SE, min, max or CI.
// The fixture is chosen so that the adjusted values differ from the raw ones
// in both mean and SD: raw {1, 2, 3} has mean 2 and SD 1, while the
// ceiling-adjusted values {5, 5, 5} would have mean 5 and SD 0.
void test_no_rounding_before_statistics() {
    const std::vector<double> estimates = {1.0, 2.0, 3.0};
    const ReplicateSummary s = summarize_replicates(estimates, 0.95);

    check(s.mean_l_s_person == 2.0, "the mean is the raw mean 2, not the adjusted mean 5");
    check(s.sample_sd_l_s_person && *s.sample_sd_l_s_person == 1.0,
          "the sample SD is the raw SD 1, not the adjusted SD 0");
    check(s.standard_error_l_s_person && approx(*s.standard_error_l_s_person, 0.5773502691896258),
          "the standard error is computed from the raw SD");
    check(s.ci_lower_l_s_person && approx(*s.ci_lower_l_s_person, -0.48413771175033027),
          "the CI lower bound follows from the raw mean and raw SE");
    check(s.ci_upper_l_s_person && approx(*s.ci_upper_l_s_person, 4.48413771175033),
          "the CI upper bound follows from the raw mean and raw SE");
    check(s.min_l_s_person == 1.0 && s.max_l_s_person == 3.0,
          "min and max are raw, not adjusted");
    check(s.adjusted_ceiling_to_5_l_s_person == 5.0,
          "the adjusted comparator is 5 * ceil(2 / 5) = 5, reported separately");
    check(bitwise_equal(s.replicate_estimates_l_s_person, estimates),
          "the preserved estimate vector is raw and unadjusted");

    // Explicit counter-fixture: summarizing the ceiling-adjusted values gives
    // a demonstrably different answer, proving the adjusted values were not
    // the ones used above.
    const ReplicateSummary adjusted_instead = summarize_replicates({5.0, 5.0, 5.0}, 0.95);
    check(adjusted_instead.mean_l_s_person != s.mean_l_s_person &&
              *adjusted_instead.sample_sd_l_s_person != *s.sample_sd_l_s_person,
          "summarizing the adjusted values would give a different mean and SD");

    // 5-grid boundary fixture: every estimate straddles a multiple of 5.
    const std::vector<double> grid = {4.999999, 5.0, 5.000001};
    const ReplicateSummary g = summarize_replicates(grid, 0.95);
    check(approx(g.mean_l_s_person, 5.0), "the 5-grid fixture's raw mean is 5");
    check(g.sample_sd_l_s_person && *g.sample_sd_l_s_person > 0.0,
          "the 5-grid fixture has a strictly positive raw SD, so no value was snapped to 5");
    check(g.adjusted_ceiling_to_5_l_s_person == 5.0,
          "the 5-grid fixture's adjusted comparator is 5");
    check(g.min_l_s_person == 4.999999 && g.max_l_s_person == 5.000001,
          "the 5-grid fixture's raw min and max are the unsnapped extremes");
}

// The frozen percentile convention is owned by Stage 10
// (simulation::percentile_linear, D8 / section 12). Stage 11 exercises it
// through that helper rather than reimplementing it, so ownership stays with
// one module. Expected values from numpy.percentile(..., method="linear").
void test_percentile_convention_is_stage_10_owned() {
    using ecai_cpp2::simulation::percentile_linear;

    const std::vector<double> odd = {1.0, 2.0, 4.0, 8.0, 16.0};
    const std::vector<double> even = {1.0, 2.0, 4.0, 8.0};
    const std::vector<double> duplicates = {3.0, 3.0, 3.0, 7.0, 7.0, 9.0};

    check(approx(percentile_linear(odd, 0.0), 1.0), "odd-length p = 0 is the minimum");
    check(approx(percentile_linear(odd, 0.25), 2.0), "odd-length p = 0.25 matches numpy");
    check(approx(percentile_linear(odd, 0.5), 4.0), "odd-length median matches numpy");
    check(approx(percentile_linear(odd, 0.96), 14.719999999999999),
          "odd-length p = 0.96 matches numpy");
    check(approx(percentile_linear(odd, 1.0), 16.0), "odd-length p = 1 is the maximum");

    check(approx(percentile_linear(even, 0.25), 1.75), "even-length p = 0.25 matches numpy");
    check(approx(percentile_linear(even, 0.5), 3.0), "even-length median interpolates to 3");
    check(approx(percentile_linear(even, 0.96), 7.52), "even-length p = 0.96 matches numpy");

    check(approx(percentile_linear(duplicates, 0.25), 3.0),
          "duplicate-valued p = 0.25 matches numpy");
    check(approx(percentile_linear(duplicates, 0.5), 5.0),
          "duplicate-valued median matches numpy");
    check(approx(percentile_linear(duplicates, 0.96), 8.6),
          "duplicate-valued p = 0.96 matches numpy");

    // The statistics module must not own or duplicate the percentile
    // convention: it consumes already-computed per-replicate estimates.
    const std::string source = read_file(ECAI_CPP2_STATISTICS_SOURCE_PATH);
    check(!source.empty(), "statistics.cpp is readable for the ownership check");
    check(source.find("percentile") == std::string::npos ||
              source.find("std::sort") == std::string::npos,
          "statistics.cpp does not reimplement a sorting percentile interpolation");
    check(source.find("std::sort") == std::string::npos,
          "statistics.cpp sorts nothing: it preserves the supplied replicate order");
}

// The documented estimator semantics and the absence of any silent data
// treatment are part of the contract, per IMPLEMENTATION_PLAN.md section 4.5.
void test_documented_assumptions_and_scope() {
    const std::string header = read_file(ECAI_CPP2_STATISTICS_HEADER_PATH);
    check(!header.empty(), "statistics.h is readable for the wording check");

    const char* required[] = {
        "approximate randomized-LHS Monte Carlo",
        "mean of the replicate percentile estimator",
        "not a confidence interval for the population percentile",
        "individual ECAi",
        "model-input uncertainty",
        "approximate-normality",
        "R - 1",
        "5 * ceil(raw_mean / 5)",
    };
    for (const char* token : required) {
        check(header.find(token) != std::string::npos,
              std::string("statistics.h documents \"") + token + "\"");
    }

    // No trimming, outlier removal or normality "correction" anywhere.
    const std::string source = read_file(ECAI_CPP2_STATISTICS_SOURCE_PATH);
    for (const std::string& token :
         {std::string("trim"), std::string("outlier"), std::string("winsor"),
          std::string("bootstrap"), std::string("shapiro")}) {
        check(source.find(token) == std::string::npos,
              "statistics.cpp contains no \"" + token + "\" operation");
        check(header.find(token) == std::string::npos,
              "statistics.h declares no \"" + token + "\" operation");
    }
}

void test_multi_replicate_orchestration() {
    const OccupancyCategory classroom = lookup_occupancy_category("Classroom");
    const EffectiveConfiguration config = make_batch_config(classroom, 64, 5);
    SimulationOptions options;
    options.collect_audit = true;

    const ReplicateBatch batch =
        ecai_cpp2::simulation::run_replicates(485, config, classroom, {}, options);
    check(batch.replicates.size() == 5,
          "run_replicates preserves exactly R = 5 complete results in index order");

    std::vector<double> raw;
    raw.reserve(batch.replicates.size());
    bool every_result_reconciles = true;
    bool every_result_keeps_audit = true;
    bool exact_replicate_index_order = true;
    bool at_least_one_raw_is_not_its_adjusted_value = false;
    for (std::size_t index = 0; index < batch.replicates.size(); ++index) {
        const auto& result = batch.replicates[index];
        raw.push_back(result.percentile_l_s_person);
        const auto direct =
            ecai_cpp2::simulation::run_replicate(485, config, index, classroom, {}, options);
        exact_replicate_index_order =
            exact_replicate_index_order &&
            bitwise_equal(result.ecai_l_s_person, direct.ecai_l_s_person) &&
            bitwise_equal(result.percentile_l_s_person, direct.percentile_l_s_person);
        every_result_reconciles =
            every_result_reconciles && result.ecai_l_s_person.size() == 64 &&
            result.counters.n == 64 &&
            result.counters.zero_infector_count + result.counters.active_count == 64 &&
            result.counters.qer_evaluations == result.counters.k_total &&
            result.counters.ecai_evaluations == result.counters.active_count &&
            result.counters.susceptible_values_consumed == 3 * result.counters.active_count &&
            result.counters.infector_values_consumed == 8 * result.counters.k_total &&
            result.cursors_exhausted;
        every_result_keeps_audit = every_result_keeps_audit && result.audit_collected &&
                                   result.audit.rows.size() == 64;
        at_least_one_raw_is_not_its_adjusted_value =
            at_least_one_raw_is_not_its_adjusted_value ||
            result.percentile_l_s_person != ceiling_to_5(result.percentile_l_s_person);
    }
    check(every_result_reconciles,
          "every batch result preserves N outputs and all Stage 10 counters reconcile");
    check(every_result_keeps_audit,
          "SimulationOptions are forwarded unchanged and every complete audit is retained");
    check(exact_replicate_index_order,
          "batch slot i bitwise equals a direct run_replicate at replicate index i");
    check(bitwise_equal(batch.summary.replicate_estimates_l_s_person, raw),
          "the summary raw vector bitwise equals result percentiles in replicate-index order");
    check(at_least_one_raw_is_not_its_adjusted_value,
          "the fixture distinguishes raw replicate percentiles from ceiling-adjusted values");

    const ReplicateSummary independently_summarized = summarize_replicates(raw, 0.95);
    check(batch.summary.replicate_count == 5 && batch.summary.confidence_level == 0.95,
          "the batch summary reports exact count and configured confidence");
    check(bitwise_equal(batch.summary.mean_l_s_person,
                        independently_summarized.mean_l_s_person) &&
              bitwise_equal(batch.summary.min_l_s_person,
                            independently_summarized.min_l_s_person) &&
              bitwise_equal(batch.summary.max_l_s_person,
                            independently_summarized.max_l_s_person),
          "batch mean, minimum and maximum reconcile from raw percentiles only");
    check(summaries_bitwise_equal(batch.summary, independently_summarized),
          "the batch summary is exactly statistics::summarize_replicates(raw, confidence)");

    const ReplicateBatch repeat =
        ecai_cpp2::simulation::run_replicates(485, config, classroom, {}, options);
    bool same_results = repeat.replicates.size() == batch.replicates.size();
    for (std::size_t i = 0; same_results && i < batch.replicates.size(); ++i) {
        same_results = bitwise_equal(batch.replicates[i].ecai_l_s_person,
                                     repeat.replicates[i].ecai_l_s_person) &&
                       bitwise_equal(batch.replicates[i].percentile_l_s_person,
                                     repeat.replicates[i].percentile_l_s_person);
    }
    check(same_results && summaries_bitwise_equal(batch.summary, repeat.summary),
          "the same base seed reproduces all results and the complete summary bitwise");

    const ReplicateBatch changed =
        ecai_cpp2::simulation::run_replicates(486, config, classroom, {}, options);
    check(!bitwise_equal(batch.summary.replicate_estimates_l_s_person,
                         changed.summary.replicate_estimates_l_s_person),
          "changing the base seed changes at least one raw replicate estimate");
    const ReplicateBatch after_changed =
        ecai_cpp2::simulation::run_replicates(485, config, classroom, {}, options);
    check(summaries_bitwise_equal(batch.summary, after_changed.summary),
          "rerunning the original seed after another batch has no leaked state");
}

void test_multi_replicate_edge_policies() {
    const OccupancyCategory classroom = lookup_occupancy_category("Classroom");

    EffectiveOverrides zero_overrides;
    zero_overrides.community_infection_rate = 0.0;
    const EffectiveConfiguration zero_config =
        make_batch_config(classroom, 16, 3, zero_overrides);
    const ReplicateBatch zeros = ecai_cpp2::simulation::run_replicates(
        485, zero_config, classroom, zero_overrides);
    bool all_zero_replicates = zeros.replicates.size() == 3;
    for (const auto& result : zeros.replicates) {
        all_zero_replicates = all_zero_replicates && result.counters.zero_infector_count == 16 &&
                              result.counters.active_count == 0 && result.counters.k_total == 0 &&
                              result.counters.qer_evaluations == 0 &&
                              result.counters.ecai_evaluations == 0 &&
                              result.percentile_l_s_person == 0.0;
    }
    check(all_zero_replicates,
          "p = 0 produces all-zero percentiles and zero-infector/no-model-work counters");
    check(bitwise_equal(zeros.summary.replicate_estimates_l_s_person,
                        std::vector<double>{0.0, 0.0, 0.0}) &&
              zeros.summary.mean_l_s_person == 0.0 &&
              zeros.summary.adjusted_ceiling_to_5_l_s_person == 0.0,
          "the p = 0 batch summary preserves three raw positive-zero estimates");

    const ReplicateBatch one = ecai_cpp2::simulation::run_replicates(
        485, make_batch_config(classroom, 16, 1), classroom);
    check(one.replicates.size() == 1 && one.summary.replicate_count == 1,
          "R = 1 produces one typed batch result");
    check(one.summary.dispersion_status == DispersionStatus::unavailable_single_replicate &&
              !one.summary.sample_sd_l_s_person && !one.summary.standard_error_l_s_person &&
              !one.summary.ci_lower_l_s_person && !one.summary.ci_upper_l_s_person,
          "R = 1 batch dispersion is typed unavailable rather than fabricated");

    EffectiveConfiguration zero_r = make_batch_config(classroom, 16, 0);
    check(throws_invalid_argument(
              [&] { ecai_cpp2::simulation::run_replicates(485, zero_r, classroom); }),
          "R = 0 is rejected before replicate execution");

    SimulationOptions impossible_work;
    impossible_work.policy.max_simulations = 0;
    const double invalid_confidences[] = {
        0.0,
        -0.1,
        1.0,
        1.1,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity(),
    };
    const char* invalid_confidence_names[] = {
        "zero", "negative", "one", "above one", "NaN", "infinity",
    };
    for (std::size_t i = 0; i < sizeof(invalid_confidences) / sizeof(invalid_confidences[0]);
         ++i) {
        EffectiveConfiguration invalid_confidence = make_batch_config(classroom, 16, 5);
        invalid_confidence.confidence_level = invalid_confidences[i];
        check(throws_invalid_argument([&] {
                  ecai_cpp2::simulation::run_replicates(485, invalid_confidence, classroom, {},
                                                         impossible_work);
              }),
              std::string("confidence ") + invalid_confidence_names[i] +
                  " is rejected in batch preflight before a replicate resource check");
    }

    EffectiveConfiguration impossible_r = make_batch_config(classroom, 16, 5);
    impossible_r.replicates_r = std::numeric_limits<std::uint64_t>::max();
    check(throws_invalid_argument(
              [&] { ecai_cpp2::simulation::run_replicates(485, impossible_r, classroom); }),
          "a replicate count exceeding vector max_size is rejected before reserve");

    EffectiveConfiguration unsupported_r = make_batch_config(classroom, 16, 5);
    unsupported_r.replicates_r =
        ecai_cpp2::statistics::kMaxSupportedReplicateCount + 1ULL;
    check(throws_invalid_argument(
              [&] { ecai_cpp2::simulation::run_replicates(485, unsupported_r, classroom); }),
          "R above 10,000,001 is rejected in batch preflight before allocation or work");
}

void test_multi_replicate_structural_defaults() {
    const EffectiveConfiguration defaults;
    check(defaults.replicates_r == 500,
          "EffectiveConfiguration structurally defaults to R = 500 without running it");
    check(defaults.confidence_level == 0.95,
          "EffectiveConfiguration structurally defaults to confidence = 0.95 without running it");
}

}  // namespace

int main() {
    test_mean_sd_se_known_vector();
    test_student_t_critical_values();
    test_student_t_low_confidence_df1();
    test_student_t_tiny_confidence_domain();
    test_student_t_supported_df_boundary();
    test_student_t_branch_boundaries_are_monotone();
    test_student_t_extreme_confidence_and_validation();
    test_confidence_interval_known_answers();
    test_zero_variance_vector();
    test_single_replicate_policy();
    test_invalid_inputs_rejected();
    test_high_offset_statistics();
    test_overflow_rejected();
    test_confidence_interval_overflow_rejected_directly();
    test_ceiling_to_5_boundaries();
    test_no_rounding_before_statistics();
    test_percentile_convention_is_stage_10_owned();
    test_documented_assumptions_and_scope();
    test_multi_replicate_orchestration();
    test_multi_replicate_edge_policies();
    test_multi_replicate_structural_defaults();

    std::cout << g_checks << " checks executed.\n";
    if (g_checks == 0) {
        std::cerr << "FAIL: no checks executed.\n";
        return 1;
    }
    if (g_failures > 0) {
        std::cerr << g_failures << " test(s) failed.\n";
        return 1;
    }
    std::cout << "All tests passed.\n";
    return 0;
}
