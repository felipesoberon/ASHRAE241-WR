// Stage 5: inverse CDFs (percent-point functions). See include/distributions.h
// for the full parameter/endpoint-policy contract and its
// docs/MODEL_SPECIFICATION.md provenance.

#include "distributions.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace ecai_cpp2 {
namespace distributions {

namespace {

constexpr double kPi = 3.14159265358979323846;

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::invalid_argument(message);
    }
}

// Acklam's rational approximation for the standard normal ppf, followed by
// one Halley refinement step against erfc (max relative error before
// refinement ~1.15e-9; refinement tightens further). `u` must already be
// validated as finite and strictly inside (0, 1).
double acklam_normal_ppf(double u) {
    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,
                                -2.759285188461852e+02, 1.383577518672690e+02,
                                -3.066479806614716e+01, 2.506628277453239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02,
                                -1.556989798598866e+02, 6.680131188771972e+01,
                                -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
                                -2.400758277161838e+00, -2.549732539343734e+00,
                                4.374664141118049e+00,  2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
                                2.445134137142996e+00, 3.754408661907416e+00};
    const double plow = 0.02425, phigh = 1.0 - plow;
    double q, r, x;
    if (u < plow) {
        q = std::sqrt(-2.0 * std::log(u));
        x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    } else if (u <= phigh) {
        q = u - 0.5;
        r = q * q;
        x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
    } else {
        q = std::sqrt(-2.0 * std::log(1.0 - u));
        x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    // At the smallest positive doubles, the ordinary Halley derivative
    // overflows. Refine there by solving log(Phi(x)) = log(u) with Newton's
    // method and the negative-tail Mills-ratio expansion. Ten alternating
    // terms are ample at |x| > 37; two Newton steps reach double precision.
    if (u < 1e-250) {
        const double log_u = std::log(u);
        for (int iteration = 0; iteration < 2; ++iteration) {
            const double y = -x;
            const double inv_y2 = 1.0 / (y * y);
            double term = 1.0;
            double mills_series = 1.0;
            for (int k = 1; k < 10; ++k) {
                term *= -static_cast<double>(2 * k - 1) * inv_y2;
                mills_series += term;
            }
            const double log_cdf = -0.5 * y * y - 0.5 * std::log(2.0 * kPi) -
                                   std::log(y) + std::log(mills_series);
            const double log_cdf_derivative = y / mills_series;
            x -= (log_cdf - log_u) / log_cdf_derivative;
        }
        return x;
    }

    // One Halley refinement step against erfc, tightening the ~1.15e-9
    // rational-approximation error to comfortably clear the 1e-9 tolerance.
    //
    // e must equal Phi(x) - u without catastrophic cancellation. For u > 0.5
    // (upper tail), both Phi(x) and u are close to 1, so computing
    // 0.5*erfc(-x/sqrt2) - u loses precision by subtracting two near-1
    // values; reflecting through the lower tail -- (1-u) - (1-Phi(x)), where
    // 1-Phi(x) = 0.5*erfc(x/sqrt2) -- keeps both operands small instead.
    const double e = (u > 0.5) ? (1.0 - u) - 0.5 * std::erfc(x / std::sqrt(2.0))
                                : 0.5 * std::erfc(-x / std::sqrt(2.0)) - u;
    const double u1 = e * std::sqrt(2.0 * kPi) * std::exp(x * x / 2.0);
    const double denom = 1.0 + x * u1 / 2.0;
    const double refined = x - u1 / denom;
    // At extreme tails (e.g. u = denorm_min), exp(x*x/2) overflows to +inf,
    // poisoning the refinement with a nonfinite result. The unrefined
    // rational approximation is still finite and within its own ~1.15e-9
    // advertised accuracy, so prefer it over propagating a NaN/inf.
    if (std::isfinite(e) && std::isfinite(u1) && std::isfinite(denom) && std::isfinite(refined)) {
        return refined;
    }
    return x;
}

// The "<name>: <detail>" message is assembled only on the failure path so
// the per-draw validators below never allocate when the input is valid.
[[noreturn]] void throw_named_invalid_argument(const char* name, const char* detail) {
    throw std::invalid_argument(std::string(name) + detail);
}

void require_open_unit_u(double u, const char* name) {
    if (!std::isfinite(u)) throw_named_invalid_argument(name, ": u must be finite");
    if (!(u > 0.0 && u < 1.0)) {
        throw_named_invalid_argument(name, ": u must be strictly inside (0, 1)");
    }
}

void require_closed_unit_u(double u, const char* name) {
    if (!std::isfinite(u)) throw_named_invalid_argument(name, ": u must be finite");
    if (!(u >= 0.0 && u <= 1.0)) {
        throw_named_invalid_argument(name, ": u must be in the closed interval [0, 1]");
    }
}

// ------------------------------------------------------------------
// Regularized incomplete beta: bounded Lentz continued fraction (Numerical
// Recipes betacf form).
// ------------------------------------------------------------------
constexpr int kBetaCfMaxIterations = 200;
constexpr double kBetaCfEps = 3e-16;
constexpr double kBetaCfFpMin = 1e-300;

double incomplete_beta_continued_fraction(double a, double b, double x) {
    const double qab = a + b, qap = a + 1.0, qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::fabs(d) < kBetaCfFpMin) d = kBetaCfFpMin;
    d = 1.0 / d;
    double h = d;
    bool converged = false;
    for (int m = 1, m2 = 2; m <= kBetaCfMaxIterations; ++m, m2 += 2) {
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < kBetaCfFpMin) d = kBetaCfFpMin;
        c = 1.0 + aa / c;
        if (std::fabs(c) < kBetaCfFpMin) c = kBetaCfFpMin;
        d = 1.0 / d;
        h *= d * c;

        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < kBetaCfFpMin) d = kBetaCfFpMin;
        c = 1.0 + aa / c;
        if (std::fabs(c) < kBetaCfFpMin) c = kBetaCfFpMin;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;

        if (std::fabs(del - 1.0) < kBetaCfEps) {
            converged = true;
            break;
        }
    }
    if (!converged || !std::isfinite(h)) {
        throw std::runtime_error(
            "regularized_incomplete_beta: continued fraction failed to converge within " +
            std::to_string(kBetaCfMaxIterations) + " iterations (a=" + std::to_string(a) +
            ", b=" + std::to_string(b) + ", x=" + std::to_string(x) + ")");
    }
    return h;
}

}  // namespace

double normal_ppf(double u) {
    require_open_unit_u(u, "normal_ppf");
    return acklam_normal_ppf(u);
}

double normal_ppf(double u, double mu, double sigma) {
    require(std::isfinite(mu), "normal_ppf: mu must be finite");
    require(std::isfinite(sigma) && sigma > 0.0, "normal_ppf: sigma must be finite and positive");
    return mu + sigma * normal_ppf(u);
}

double lognormal_gm_gsd_ppf(double u, double gm, double gsd) {
    require(std::isfinite(gm) && gm > 0.0, "lognormal_gm_gsd_ppf: gm must be finite and positive");
    require(std::isfinite(gsd) && gsd > 1.0,
            "lognormal_gm_gsd_ppf: gsd must be finite and strictly greater than 1");
    const double mu = std::log(gm);
    const double sigma = std::log(gsd);
    return std::exp(normal_ppf(u) * sigma + mu);
}

double log10normal_ppf(double u, double mu, double sigma) {
    require(std::isfinite(mu), "log10normal_ppf: mu must be finite");
    require(std::isfinite(sigma) && sigma > 0.0,
            "log10normal_ppf: sigma must be finite and positive");
    return std::pow(10.0, normal_ppf(u) * sigma + mu);
}

double uniform_ppf(double u, double min, double max) {
    require_closed_unit_u(u, "uniform_ppf");
    require(std::isfinite(min) && std::isfinite(max), "uniform_ppf: min and max must be finite");
    require(min < max, "uniform_ppf: min must be strictly less than max");
    if (u <= 0.0) return min;
    if (u >= 1.0) return max;
    return min + (max - min) * u;
}

double regularized_incomplete_beta(double x, double a, double b) {
    require(std::isfinite(a) && a > 0.0, "regularized_incomplete_beta: a must be finite and positive");
    require(std::isfinite(b) && b > 0.0, "regularized_incomplete_beta: b must be finite and positive");
    require(std::isfinite(x), "regularized_incomplete_beta: x must be finite");
    require(x >= 0.0 && x <= 1.0, "regularized_incomplete_beta: x must be in the closed interval [0, 1]");
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    const double log_beta = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b);
    const double bt = std::exp(log_beta + a * std::log(x) + b * std::log1p(-x));
    if (x < (a + 1.0) / (a + b + 2.0)) {
        return bt * incomplete_beta_continued_fraction(a, b, x) / a;
    }
    return 1.0 - bt * incomplete_beta_continued_fraction(b, a, 1.0 - x) / b;
}

namespace {

// Bisection on regularized_incomplete_beta for the standard (unscaled,
// support [0,1]) beta ppf. `u` is already validated by the caller (beta_ppf)
// to be finite and inside the closed interval [0, 1]; endpoints are handled
// by the caller directly and never reach this function.
double beta_ppf_standard(double u, double a, double b) {
    double lo = 0.0, hi = 1.0;
    constexpr int kMaxIterations = 200;
    constexpr double kToleranceWidth = 1e-15;
    for (int i = 0; i < kMaxIterations && (hi - lo) > kToleranceWidth; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double cdf = regularized_incomplete_beta(mid, a, b);
        if (cdf < u) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

}  // namespace

double beta_ppf(double u, double a, double b, double min, double max) {
    require_closed_unit_u(u, "beta_ppf");
    require(std::isfinite(a) && a > 0.0, "beta_ppf: a must be finite and positive");
    require(std::isfinite(b) && b > 0.0, "beta_ppf: b must be finite and positive");
    require(std::isfinite(min) && std::isfinite(max), "beta_ppf: min and max must be finite");
    require(min < max, "beta_ppf: min must be strictly less than max");
    if (u <= 0.0) return min;
    if (u >= 1.0) return max;
    return min + (max - min) * beta_ppf_standard(u, a, b);
}

namespace {

// Neumaier's improved Kahan compensated summation accumulator, used so the
// binomial normalizer and cumulative scan below are not dominated by
// rounding error from summing up to n+1 PMF terms.
class CompensatedSum {
  public:
    void add(double value) {
        const double t = sum_ + value;
        if (std::fabs(sum_) >= std::fabs(value)) {
            compensation_ += (sum_ - t) + value;
        } else {
            compensation_ += (value - t) + sum_;
        }
        sum_ = t;
    }
    double total() const { return sum_ + compensation_; }

  private:
    double sum_ = 0.0;
    double compensation_ = 0.0;
};

// Shared (n, p) parameter contract of binomial_ppf and BinomialPpfTable.
// `name` prefixes every diagnostic so each entry point keeps its own
// message text.
void require_binomial_parameters(long long n, double p, const char* name) {
    if (n < 0) throw_named_invalid_argument(name, ": n must be non-negative");
    if (!(std::isfinite(p) && p >= 0.0 && p <= 1.0)) {
        throw_named_invalid_argument(name, ": p must be finite and in [0, 1]");
    }
}

// Guard the n -> size_t widening and the (n+1)-length allocation against
// overflow *before* touching any vector, per a fixed platform limit
// (vector<double>::max_size()) rather than discovering the overflow via UB
// in the cast or an allocation-time exception. Returns the buffer length
// n+1. `n` must already be validated as non-negative.
std::size_t checked_binomial_buffer_count(long long n, const char* name) {
    const unsigned long long un = static_cast<unsigned long long>(n);
    const unsigned long long max_supported =
        static_cast<unsigned long long>(std::vector<double>().max_size()) - 1ULL;
    if (un > max_supported) {
        throw std::runtime_error(std::string(name) + ": n=" + std::to_string(n) +
                                 " exceeds the maximum size supported for the PMF buffer");
    }
    return static_cast<std::size_t>(n) + 1;
}

// The mode-centered PMF recurrence shared, term for term and in the same
// evaluation order, by binomial_ppf and BinomialPpfTable. `pmf` must be
// n+1 zero-initialized doubles; on return it holds the unnormalized PMF and
// the compensated normalizer (the sum of every term) is returned. Requires
// n >= 1 and p strictly inside (0, 1): the callers short-circuit every
// degenerate case before reaching the recurrence.
//
// Computing pmf(0) = (1-p)^n and scanning forward underflows to exactly 0.0
// whenever (1-p)^n is smaller than the smallest representable double (e.g.
// n=400, p=0.9: 0.1^400), after which multiplying the zero PMF by finite
// ratios stays zero forever. Starting from the mode (where the PMF is
// largest, so it is the term least likely to underflow) and filling outward
// in both directions avoids that failure mode entirely.
double fill_binomial_pmf(long long n, double p, std::vector<double>& pmf) {
    long long mode =
        static_cast<long long>(std::floor((static_cast<double>(n) + 1.0) * p));
    mode = std::max<long long>(0, std::min<long long>(mode, n));

    const double log_p = std::log(p);
    const double log_1mp = std::log1p(-p);
    const double log_pmf_mode = std::lgamma(static_cast<double>(n) + 1.0) -
                                 std::lgamma(static_cast<double>(mode) + 1.0) -
                                 std::lgamma(static_cast<double>(n - mode) + 1.0) +
                                 static_cast<double>(mode) * log_p +
                                 static_cast<double>(n - mode) * log_1mp;

    pmf[static_cast<std::size_t>(mode)] = std::exp(log_pmf_mode);

    // Fill upward from the mode: pmf(k+1) = pmf(k) * (n-k)*p / ((k+1)*(1-p)).
    for (long long k = mode; k < n; ++k) {
        const double next = pmf[static_cast<std::size_t>(k)] *
                             (static_cast<double>(n - k) * p) /
                             (static_cast<double>(k + 1) * (1.0 - p));
        pmf[static_cast<std::size_t>(k + 1)] = next;
        if (next == 0.0) break;  // Remaining upper tail is negligible.
    }
    // Fill downward from the mode: pmf(k-1) = pmf(k) * k*(1-p) / ((n-k+1)*p).
    for (long long k = mode; k > 0; --k) {
        const double prev = pmf[static_cast<std::size_t>(k)] *
                             (static_cast<double>(k) * (1.0 - p)) /
                             (static_cast<double>(n - k + 1) * p);
        pmf[static_cast<std::size_t>(k - 1)] = prev;
        if (prev == 0.0) break;  // Remaining lower tail is negligible.
    }

    CompensatedSum normalizer_sum;
    for (double v : pmf) normalizer_sum.add(v);
    return normalizer_sum.total();
}

}  // namespace

long long binomial_ppf(double u, long long n, double p) {
    require(std::isfinite(u), "binomial_ppf: u must be finite");
    require(u > 0.0 && u <= 1.0, "binomial_ppf: u must be in the half-open interval (0, 1]");
    require_binomial_parameters(n, p, "binomial_ppf");

    if (n == 0) return 0;
    // Endpoint convention is intentionally aligned with scipy.stats.binom.ppf:
    // q == 1 returns the upper support endpoint n, including degenerate p=0.
    if (u >= 1.0) return n;
    if (p <= 0.0) return 0;
    if (p >= 1.0) return n;

    const std::size_t count = checked_binomial_buffer_count(n, "binomial_ppf");
    std::vector<double> pmf(count, 0.0);
    const double normalizer = fill_binomial_pmf(n, p, pmf);

    CompensatedSum cumulative;
    for (long long k = 0; k <= n; ++k) {
        cumulative.add(pmf[static_cast<std::size_t>(k)]);
        if (cumulative.total() / normalizer >= u) return k;
    }
    // Explicit, documented fallback: floating-point rounding of the
    // normalized cumulative sum left it short of u even at k = n (whose true
    // CDF is exactly 1). Bounded by the n+1 terms evaluated above; never
    // loops indefinitely.
    return n;
}

BinomialPpfTable::BinomialPpfTable(long long n, double p) : n_(n), p_(p) {
    require_binomial_parameters(n, p, "BinomialPpfTable");
    // The same degenerate cases binomial_ppf short-circuits before its
    // recurrence need no table at all; ppf() answers them directly.
    if (n == 0 || p <= 0.0 || p >= 1.0) return;

    const std::size_t count = checked_binomial_buffer_count(n, "BinomialPpfTable");
    std::vector<double> pmf(count, 0.0);
    const double normalizer = fill_binomial_pmf(n, p, pmf);

    // Store exactly the value binomial_ppf compares against u at each k --
    // the same compensated cumulative sum, divided by the same normalizer,
    // in the same order -- so ppf() reproduces its decision bit for bit.
    cdf_.assign(count, 0.0);
    CompensatedSum cumulative;
    for (std::size_t k = 0; k < count; ++k) {
        cumulative.add(pmf[k]);
        cdf_[k] = cumulative.total() / normalizer;
    }
}

long long BinomialPpfTable::ppf(double u) const {
    require(std::isfinite(u), "BinomialPpfTable::ppf: u must be finite");
    require(u > 0.0 && u <= 1.0,
            "BinomialPpfTable::ppf: u must be in the half-open interval (0, 1]");

    if (n_ == 0) return 0;
    if (u >= 1.0) return n_;
    if (p_ <= 0.0) return 0;
    if (p_ >= 1.0) return n_;

    // Exact linear scan over the stored comparison values: the first k whose
    // stored value is >= u is precisely the k at which binomial_ppf's
    // cumulative scan returns.
    for (std::size_t k = 0; k < cdf_.size(); ++k) {
        if (cdf_[k] >= u) return static_cast<long long>(k);
    }
    // Same documented fallback as binomial_ppf.
    return n_;
}

}  // namespace distributions
}  // namespace ecai_cpp2
