// Stage 11: across-replicate statistics. See include/statistics.h for the
// full contract, the frozen formulas, and the scope of the reported interval.

#include "statistics.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "distributions.h"

namespace ecai_cpp2 {
namespace statistics {

namespace {

void require_finite(double value, const std::string& name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument("statistics: " + name + " must be finite");
    }
}

// Open interval (0, 1), both endpoints excluded.
void require_open_unit(double value, const std::string& name) {
    require_finite(value, name);
    if (!(value > 0.0) || !(value < 1.0)) {
        throw std::invalid_argument("statistics: " + name +
                                     " must lie strictly inside the open interval (0, 1)");
    }
}

void require_finite_result(double value, const std::string& name) {
    if (!std::isfinite(value)) {
        throw std::runtime_error("statistics: " + name + " overflowed to a nonfinite value");
    }
}

double checked_double_result(long double value, const std::string& name) {
    if (!std::isfinite(value) ||
        std::fabs(value) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::runtime_error("statistics: " + name +
                                 " overflowed the finite double result range");
    }
    const double narrowed = static_cast<double>(value);
    require_finite_result(narrowed, name);
    return narrowed;
}

void compensated_add(long double value, long double& sum, long double& correction) {
    const long double next = sum + value;
    if (std::fabs(sum) >= std::fabs(value)) {
        correction += (sum - next) + value;
    } else {
        correction += (value - next) + sum;
    }
    sum = next;
}

// Student-t survival function P(T > t) for t >= 0 and df > 0, in the
// standard incomplete-beta form
//   P(T > t) = 0.5 * I_{df/(df + t^2)}(df/2, 1/2).
// t * t may overflow to infinity for a very large bracket candidate; the
// quotient then evaluates to exactly 0, I_0 = 0, and the survival value is
// 0, which is the mathematically correct limit and keeps the bracket search
// finite rather than producing a NaN.
double student_t_survival(double t, double df) {
    const double t_squared = t * t;
    const double x = std::isfinite(t_squared) ? (df / (df + t_squared)) : 0.0;
    return 0.5 * distributions::regularized_incomplete_beta(x, 0.5 * df, 0.5);
}

// Density at zero for Student-t(df). At large df, a gamma-ratio expansion
// avoids subtracting two large, nearly equal lgamma values.
long double student_t_density_at_zero(double df) {
    const long double nu = static_cast<long double>(df);
    if (df >= 100.0) {
        const long double inv_nu = 1.0L / nu;
        const long double inv_nu_2 = inv_nu * inv_nu;
        const long double ratio =
            1.0L - inv_nu / 4.0L + inv_nu_2 / 32.0L +
            5.0L * inv_nu_2 * inv_nu / 128.0L -
            21.0L * inv_nu_2 * inv_nu_2 / 2048.0L;
        return ratio / std::sqrt(2.0L * std::acos(-1.0L));
    }
    const long double half_nu = 0.5L * nu;
    const long double log_density =
        std::lgamma(half_nu + 0.5L) - std::lgamma(half_nu) -
        0.5L * (std::log(nu) + std::log(std::acos(-1.0L)));
    return std::exp(log_density);
}

// Inverse Student-t expansion about the normal quantile. At df >= 10000,
// terms through df^-3 are well below the frozen 1e-6 absolute tolerance over
// the supported binary64 confidence domain and avoid incomplete-beta noise as
// df grows. z is obtained from the upper-tail probability without forming
// (1 + confidence)/2, which would round to 1 near the endpoint.
double student_t_large_df_critical(double confidence_level, double df) {
    const double upper_tail = 0.5 * (1.0 - confidence_level);
    const long double z =
        -static_cast<long double>(distributions::normal_ppf(upper_tail));
    const long double z2 = z * z;
    const long double z3 = z2 * z;
    const long double z5 = z3 * z2;
    const long double z7 = z5 * z2;
    const long double nu = static_cast<long double>(df);
    const long double correction_1 = (z3 + z) / (4.0L * nu);
    const long double correction_2 =
        (5.0L * z5 + 16.0L * z3 + 3.0L * z) / (96.0L * nu * nu);
    const long double correction_3 =
        (3.0L * z7 + 19.0L * z5 + 17.0L * z3 - 15.0L * z) /
        (384.0L * nu * nu * nu);
    return checked_double_result(z + correction_1 + correction_2 + correction_3,
                                 "Student-t critical value");
}

}  // namespace

std::pair<double, double> detail::checked_confidence_interval(double mean,
                                                              double standard_error,
                                                              double t_critical) {
    const double half_width = t_critical * standard_error;
    require_finite_result(half_width, "confidence interval half-width");
    const double lower = mean - half_width;
    const double upper = mean + half_width;
    require_finite_result(lower, "confidence interval lower bound");
    require_finite_result(upper, "confidence interval upper bound");
    return {lower, upper};
}

double ceiling_to_5(double raw) {
    require_finite(raw, "raw value for the ceiling-to-5 adjustment");
    const double adjusted = 5.0 * std::ceil(raw / 5.0);
    // Preserve the mathematical ceiling for positive subnormals when raw / 5
    // underflows to zero in binary64.
    if (raw > 0.0 && adjusted == 0.0) return 5.0;
    require_finite_result(adjusted, "ceiling-to-5 adjustment");
    return adjusted;
}

double student_t_two_sided_critical(double confidence_level, double df) {
    require_open_unit(confidence_level, "confidence level");
    require_finite(df, "degrees of freedom");
    if (!(df >= 1.0)) {
        throw std::invalid_argument(
            "statistics: degrees of freedom must be at least 1 (R - 1 replicate df)");
    }
    if (std::floor(df) != df) {
        throw std::invalid_argument(
            "statistics: degrees of freedom must be an integer R - 1");
    }
    if (df > kMaxSupportedStudentTDegreesOfFreedom) {
        throw std::invalid_argument(
            "statistics: degrees of freedom exceed the supported maximum of 1e7");
    }

    // For small central probability C, invert the density expansion through
    // fifth order. With y=C/(2*f(0)), t=y+A*y^3+(3*A^2-B)*y^5+O(y^7),
    // A=(df+1)/(6*df), B=(df+1)(df+3)/(40*df^2). At C<=0.1 the omitted
    // term remains below the frozen 1e-6 absolute tolerance even at df=1.
    // Long-double evaluation preserves all positive binary64 confidence values,
    // including subnormals.
    if (confidence_level <= 0.1) {
        const long double nu = static_cast<long double>(df);
        const long double density = student_t_density_at_zero(df);
        const long double y = static_cast<long double>(confidence_level) / (2.0L * density);
        const long double y2 = y * y;
        const long double a = (nu + 1.0L) / (6.0L * nu);
        const long double b = (nu + 1.0L) * (nu + 3.0L) / (40.0L * nu * nu);
        const long double critical_long =
            y * (1.0L + a * y2 + (3.0L * a * a - b) * y2 * y2);
        const double critical =
            checked_double_result(critical_long, "Student-t critical value");
        if (!(critical > 0.0)) {
            throw std::runtime_error(
                "statistics: positive confidence produced no representable positive critical value");
        }
        return critical;
    }

    if (df >= 10000.0) {
        return student_t_large_df_critical(confidence_level, df);
    }

    // Above the local-series range, alpha/2 is separated far enough from 0.5
    // to be represented accurately. One survival inversion covers the rest of
    // the confidence domain and avoids a second method transition at 0.5.
    const double target = 0.5 * (1.0 - confidence_level);

    // Bracket [lo, hi] with survival(lo) > target >= survival(hi). survival(0)
    // is exactly 0.5 > target, so lo = 0 is always a valid lower end.
    double lo = 0.0;
    double hi = 1.0;
    int doublings = 0;
    while (student_t_survival(hi, df) > target) {
        lo = hi;
        hi *= 2.0;
        if (++doublings > 200) {
            throw std::runtime_error(
                "statistics: no Student-t bracket found within 200 doublings for the requested "
                "confidence level");
        }
    }

    bool converged = false;
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (student_t_survival(mid, df) > target) {
            lo = mid;
        } else {
            hi = mid;
        }
        if (hi - lo <= 1e-15 * hi) {
            converged = true;
            break;
        }
    }
    if (!converged) {
        throw std::runtime_error(
            "statistics: Student-t survival inversion did not converge within 200 bisections");
    }

    const double critical = 0.5 * (lo + hi);
    require_finite_result(critical, "Student-t critical value");
    return critical;
}

ReplicateSummary summarize_replicates(const std::vector<double>& estimates,
                                       double confidence_level) {
    // Confidence is validated first and unconditionally, so an R = 1 summary
    // (which reports no interval) still rejects an invalid level.
    require_open_unit(confidence_level, "confidence level");

    if (estimates.empty()) {
        throw std::invalid_argument(
            "statistics: at least one replicate estimate is required (R = 0 is rejected)");
    }
    if (estimates.size() > static_cast<std::size_t>(kMaxSupportedReplicateCount)) {
        throw std::invalid_argument(
            "statistics: replicate count exceeds the supported Student-t range");
    }

    for (std::size_t i = 0; i < estimates.size(); ++i) {
        const double value = estimates[i];
        if (!std::isfinite(value)) {
            throw std::invalid_argument("statistics: replicate estimate at index " +
                                         std::to_string(i) + " is not finite");
        }
        if (value < 0.0) {
            throw std::invalid_argument("statistics: replicate estimate at index " +
                                         std::to_string(i) +
                                         " is negative (ECAi is clamped at zero)");
        }
    }

    ReplicateSummary summary;
    summary.replicate_count = estimates.size();
    summary.replicate_estimates_l_s_person = estimates;  // raw, supplied order
    summary.confidence_level = confidence_level;

    // Shift by the first input, then use compensated long-double accumulation
    // in supplied order. The shift retains small differences atop a large
    // common offset; the second pass evaluates deviations from the improved
    // mean rather than inheriting binary64 Welford rounding at every step.
    const long double reference = static_cast<long double>(estimates[0]);
    long double offset_sum = 0.0L;
    long double offset_correction = 0.0L;
    double minimum = estimates[0];
    double maximum = estimates[0];
    for (const double value : estimates) {
        compensated_add(static_cast<long double>(value) - reference, offset_sum,
                        offset_correction);
        if (value < minimum) minimum = value;
        if (value > maximum) maximum = value;
    }

    const long double r_long_double = static_cast<long double>(estimates.size());
    const long double mean_offset = (offset_sum + offset_correction) / r_long_double;
    const long double mean_long_double = reference + mean_offset;
    const double mean = checked_double_result(mean_long_double, "replicate mean");

    long double sum_squared_deviations = 0.0L;
    long double squared_correction = 0.0L;
    for (const double value : estimates) {
        const long double deviation =
            (static_cast<long double>(value) - reference) - mean_offset;
        compensated_add(deviation * deviation, sum_squared_deviations, squared_correction);
    }
    sum_squared_deviations += squared_correction;
    if (!std::isfinite(sum_squared_deviations)) {
        throw std::runtime_error(
            "statistics: sum of squared deviations overflowed to a nonfinite value");
    }

    summary.mean_l_s_person = mean;
    summary.min_l_s_person = minimum;
    summary.max_l_s_person = maximum;

    const double r = static_cast<double>(estimates.size());
    if (estimates.size() >= 2) {
        const double degrees_of_freedom = r - 1.0;
        const long double variance_long_double =
            sum_squared_deviations / static_cast<long double>(estimates.size() - 1);
        // Preserve the existing public policy: a variance outside binary64's
        // finite range is rejected even though SD itself might fit in double.
        (void)checked_double_result(variance_long_double, "sample variance");
        const double sd =
            checked_double_result(std::sqrt(variance_long_double), "sample standard deviation");
        const double se = sd / std::sqrt(r);
        require_finite_result(se, "standard error");
        const double t_critical = student_t_two_sided_critical(confidence_level,
                                                                degrees_of_freedom);
        const std::pair<double, double> interval =
            detail::checked_confidence_interval(mean, se, t_critical);

        summary.dispersion_status = DispersionStatus::available;
        summary.degrees_of_freedom = degrees_of_freedom;
        summary.sample_sd_l_s_person = sd;
        summary.standard_error_l_s_person = se;
        summary.t_critical = t_critical;
        summary.ci_lower_l_s_person = interval.first;
        summary.ci_upper_l_s_person = interval.second;
    } else {
        // R == 1: no unbiased sample dispersion exists. Every dispersion
        // field stays empty; nothing is fabricated.
        summary.dispersion_status = DispersionStatus::unavailable_single_replicate;
    }

    // Reported separately, and never an input to any statistic above.
    summary.adjusted_ceiling_to_5_l_s_person = ceiling_to_5(mean);
    return summary;
}

}  // namespace statistics
}  // namespace ecai_cpp2
