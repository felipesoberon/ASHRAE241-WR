// Stage 11: across-replicate statistics.
//
// Source of truth: IMPLEMENTATION_PLAN.md Section 7, Stage 11 and section 4.5
// (across-R mean, sample SD with denominator R - 1, SE, two-sided Student-t
// CI, min/max, exact replicate count, finite-value checks, and the
// separately reported ceiling-to-5 comparator), together with
// docs/MODEL_SPECIFICATION.md section 11.5 (frozen tolerances) and section
// 12 (frozen rounding direction D6 and the frozen percentile convention D8).
//
// This module consumes already-computed per-replicate percentile estimates.
// It does not own, compute or duplicate the percentile convention itself:
// that belongs to Stage 10 (simulation::percentile_linear), and this header
// deliberately does not include simulation.h, so the dependency runs one way
// only (simulation.h may include this header; this header never includes
// simulation.h).
//
// ---------------------------------------------------------------------------
// What the reported interval is, and what it is not
// ---------------------------------------------------------------------------
//
// The interval produced here is an approximate randomized-LHS Monte Carlo
// Student-t confidence interval for the mean of the replicate percentile estimator
// -- that is, for the expected value of the per-replicate
// percentile estimate under repeated randomized-LHS replication of the same
// effective configuration. It is explicitly:
//
//   - not a confidence interval for the population percentile of the
//     underlying ECAi distribution;
//   - not an interval for individual ECAi values or for any individual
//     simulation within a replicate;
//   - not a characterization of model-input uncertainty (distribution
//     choices, Jones et al. parameter values, category inputs);
//   - not protection against bias in the percentile estimator itself.
//
// The interval additionally rests on an approximate-normality assumption:
// the Student-t construction assumes the R replicate percentile estimates
// behave approximately as independent draws from a normal distribution.
// That assumption is documented here rather than enforced, and it is
// assessed diagnostically by the caller: ReplicateSummary preserves the raw
// per-replicate estimates in supplied order together with the raw minimum
// and maximum, which is sufficient to inspect the replicate-estimate shape
// without changing the reported estimator.
//
// No data treatment of any kind is applied: the estimates are used exactly
// as supplied, in the order supplied. Nothing is discarded, re-weighted,
// re-ordered, or adjusted for non-normality, and no such option exists.
//
// ---------------------------------------------------------------------------
// Frozen formulas
// ---------------------------------------------------------------------------
//
//   R      = number of supplied replicate estimates (>= 1)
//   mean   = sum(x_i) / R                  (shifted/compensated, supplied order)
//   SD     = sqrt( sum((x_i - mean)^2) / (R - 1) )   -- unbiased, R - 1
//   SE     = SD / sqrt(R)
//   CI     = mean -/+ t_{1 - alpha/2, R - 1} * SE,  alpha = 1 - confidence
//   min    = min(x_i),  max = max(x_i)           -- raw, unsorted scan
//   adjusted_ceiling_to_5 = 5 * ceil(raw_mean / 5)   (D6, ceiling only)
//
// Every raw statistic above is computed from the raw estimates only. The
// ceiling-adjusted value is a separately labelled comparator and never
// enters the mean, SD, SE, CI, minimum or maximum.

#ifndef ECAI_CPP2_STATISTICS_H
#define ECAI_CPP2_STATISTICS_H

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace ecai_cpp2 {
namespace statistics {

namespace detail {
// Checked primitive shared by summarize_replicates and direct overflow
// regressions. This is an implementation detail, not a stable public API.
std::pair<double, double> checked_confidence_interval(double mean, double standard_error,
                                                      double t_critical);
}  // namespace detail

// Whether the dispersion statistics (sample SD, SE, t critical, CI) exist
// for this summary. A single replicate has no unbiased sample dispersion, so
// the corresponding optionals are empty rather than zero, NaN or fabricated.
enum class DispersionStatus {
    // R >= 2: sample SD, SE, t critical and both CI bounds are populated.
    available = 0,
    // R == 1: the R - 1 denominator is zero, so no sample SD, SE, t critical
    // or CI is defined. Mean, minimum, maximum, the raw estimate vector and
    // the ceiling-to-5 comparator remain available.
    unavailable_single_replicate = 1,
};

// One across-replicate summary. Every "raw" field is in L/s per person and
// is computed from the unrounded per-replicate percentile estimates.
struct ReplicateSummary {
    // Exact replicate count R (the supplied estimate count).
    std::size_t replicate_count = 0;

    // Diagnostics: the raw per-replicate estimates, preserved bit-for-bit in
    // the supplied replicate-index order. Preserved so replicate-estimate
    // shape can be assessed without changing the estimator above.
    std::vector<double> replicate_estimates_l_s_person;

    // The confidence level the CI below was requested at, echoed back.
    double confidence_level = 0.0;

    // Raw statistics (never adjusted, never rounded).
    double mean_l_s_person = 0.0;
    double min_l_s_person = 0.0;
    double max_l_s_person = 0.0;

    DispersionStatus dispersion_status = DispersionStatus::unavailable_single_replicate;

    // Populated iff dispersion_status == available.
    std::optional<double> degrees_of_freedom;         // R - 1
    std::optional<double> sample_sd_l_s_person;       // denominator R - 1
    std::optional<double> standard_error_l_s_person;  // SD / sqrt(R)
    std::optional<double> t_critical;                 // t_{1 - alpha/2, R - 1}
    std::optional<double> ci_lower_l_s_person;        // mean - t * SE
    std::optional<double> ci_upper_l_s_person;        // mean + t * SE

    // The separately labelled ceiling-to-5 comparator, 5 * ceil(raw_mean / 5)
    // (D6). Reported alongside the raw statistics; never an input to them.
    double adjusted_ceiling_to_5_l_s_person = 0.0;
};

// The frozen multiple-of-5 adjustment, 5 * ceil(raw_mean / 5)
// (docs/MODEL_SPECIFICATION.md D6 / section 12; ceiling direction only).
// Throws std::invalid_argument for a nonfinite input, and
// std::runtime_error if the adjusted value overflows to a nonfinite result
// (defensive: unreachable for any finite input in IEEE-754 double, since
// 5 * ceil(DBL_MAX / 5) == DBL_MAX exactly).
double ceiling_to_5(double raw);

// Two-sided Student-t critical value t_{1 - alpha/2, df} with
// alpha = 1 - confidence_level, computed in pure C++17 with no external
// runtime dependency.
//
// Method: confidence levels through 0.1 use the local central-density inverse
// expansion through fifth order; its omitted term is below the frozen 1e-6
// absolute tolerance in the supported df domain. At df >= 10000, higher levels
// use the inverse Student-t expansion about the cancellation-free upper-tail
// normal quantile through df^-3. Lower df invert the exact survival identity
// P(T > t) = 0.5 * I_{df/(df + t^2)}(df/2, 1/2) using the Stage 5 regularized
// incomplete beta function and bounded bisection directly in t-space (at most
// 200 bracket doublings and 200 bisections, requiring a relative-width
// convergence test). Exhausting an iteration bound throws rather than returning
// an estimate.
// Inverting in t-space rather than in
// beta-quantile space keeps the achieved accuracy governed by the incomplete
// beta evaluation rather than by an absolute bracket width on [0, 1].
// Verified against SciPy across df in {1, 2, 3, 4, 9, 10, 30, 99, 499, 4999}
// and confidence levels {0.80, 0.90, 0.95, 0.98, 0.99, 0.999}; the frozen
// tolerance is absolute 1e-6 (docs/MODEL_SPECIFICATION.md section 11.5).
//
// Degrees of freedom exist here solely as R - 1 replicate degrees of freedom.
// The supported public domain is therefore integer 1 <= df <= 1e7. The upper limit avoids loss of the
// df/(df+t^2) decrement in the Stage 5 incomplete-beta evaluation, is verified
// against SciPy at the exact bound, and is far above the default df = 499.
inline constexpr double kMaxSupportedStudentTDegreesOfFreedom = 1e7;
inline constexpr std::uint64_t kMaxSupportedReplicateCount = 10000001ULL;

// Throws std::invalid_argument if confidence_level is nonfinite or outside
// the open interval (0, 1), or if df is nonfinite, non-integer, less than 1, or
// greater than kMaxSupportedStudentTDegreesOfFreedom.
// Throws std::runtime_error if no bracket can be established within the
// bound above, or if the result is not finite; it never returns a nonfinite
// or fabricated value.
double student_t_two_sided_critical(double confidence_level, double df);

// Summarizes raw per-replicate percentile estimates, supplied in
// replicate-index order (index 0 first). The order is preserved exactly and
// the computation is a deterministic pure function of that order: the same
// vector always yields a bit-identical summary.
//
// The mean is accumulated from long-double offsets relative to the first
// estimate, with compensation, in supplied order. A second compensated pass
// accumulates squared deviations from that shifted mean. This preserves small
// representable spreads atop a large common offset without a naive sum of the
// original values. Results are narrowed to double only after explicit finite
// range checks; sample variance outside the finite double range is rejected.
//
// Throws std::invalid_argument if:
//   - `estimates` is empty (R = 0 is rejected);
//   - R exceeds kMaxSupportedReplicateCount, so R - 1 would exceed the
//     verified Student-t degrees-of-freedom domain;
//   - any estimate is nonfinite or negative (ECAi is clamped at zero, so a
//     negative replicate percentile estimate is a caller defect);
//   - confidence_level is nonfinite or outside the open interval (0, 1) --
//     validated even at R = 1, where no CI is produced.
// Throws std::runtime_error if any intermediate overflows to a nonfinite
// value (mean, sum of squared deviations, variance, SD, SE, CI bound or the
// ceiling adjustment).
ReplicateSummary summarize_replicates(const std::vector<double>& estimates,
                                       double confidence_level);

}  // namespace statistics
}  // namespace ecai_cpp2

#endif  // ECAI_CPP2_STATISTICS_H
