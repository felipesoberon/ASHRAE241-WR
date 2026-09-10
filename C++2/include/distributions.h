// Stage 5: inverse CDFs (percent-point functions) transforming a single
// open-unit-interval LHS draw into a Jones et al. parameter distribution.
//
// Source of truth: docs/MODEL_SPECIFICATION.md section 5 (Table 1 — every
// distribution family, shape parameter, and bound used by a Jones parameter)
// and section 11.5 (validation acceptance tolerances). IMPLEMENTATION_PLAN.md
// Section 7, Stage 5.
//
// Endpoint policy (frozen here, since MODEL_SPECIFICATION.md section 11.5
// documents tolerances but not the domain convention itself):
//   - Unbounded transforms (normal, lognormal-from-GM/GSD, log10normal) are
//     defined only for u strictly inside (0, 1): the transform maps u -> 0
//     to -infinity and u -> 1 to +infinity, so an endpoint has no finite,
//     meaningful result. u <= 0, u >= 1, or a nonfinite u is rejected with
//     std::invalid_argument. The exact-size 1D LHS engine (lhs.h) never
//     actually emits 0.0 or 1.0, so this is not a live path in production,
//     only a documented input contract.
//   - Bounded transforms (uniform, beta) are defined on the closed interval
//     [0, 1]: u = 0 and u = 1 have finite, well-defined results (the
//     distribution's min and max respectively), so both endpoints are
//     accepted. u < 0, u > 1, or a nonfinite u is rejected with
//     std::invalid_argument.
//   - The binomial PPF (a discrete step function) is defined on the
// half-open interval (0, 1]: u = 1 returns the upper support endpoint n,
// matching scipy.stats.binom.ppf even for degenerate p=0; u = 0 is excluded because every k already
//     satisfies CDF(k) >= 0, so u = 0 carries no distinguishing information.
//     u <= 0, u > 1, or a nonfinite u is rejected with std::invalid_argument.
//
// Every function additionally rejects nonfinite or out-of-domain
// distribution parameters (std::invalid_argument) so no inverse CDF can
// silently produce a nonfinite output anywhere in the Jones-relevant open
// tail range.

#ifndef ECAI_CPP2_DISTRIBUTIONS_H
#define ECAI_CPP2_DISTRIBUTIONS_H

namespace ecai_cpp2 {
namespace distributions {

// Standard normal inverse CDF (quantile function). u strictly in (0, 1);
// throws std::invalid_argument for nonfinite u or u outside (0, 1). Uses
// Acklam's rational approximation (max relative error ~1.15e-9) followed by
// one Halley refinement step against erfc (computed via whichever tail keeps
// both operands small, avoiding cancellation near u close to 0 or 1),
// matching the tolerance frozen in MODEL_SPECIFICATION.md section 11.5. If
// the refinement's intermediates are nonfinite (only reachable at extreme
// tails far outside the production LHS engine's range, e.g. denorm_min),
// the unrefined rational approximation -- itself finite and within its own
// ~1.15e-9 advertised accuracy -- is returned instead of a nonfinite value.
double normal_ppf(double u);

// Location-scale normal inverse CDF: mu + sigma * normal_ppf(u). Throws
// std::invalid_argument for nonfinite mu, nonfinite or non-positive sigma,
// or an out-of-domain u (see normal_ppf(double)).
double normal_ppf(double u, double mu, double sigma);

// Natural-lognormal inverse CDF parameterized directly by geometric mean and
// geometric standard deviation (Jones Table 1 style, e.g. lambda: GM=0.52,
// GSD=1.9): exp(normal_ppf(u) * ln(gsd) + ln(gm)). Throws
// std::invalid_argument for nonfinite or non-positive gm, nonfinite gsd or
// gsd <= 1 (a GSD of exactly 1 is a degenerate zero-spread distribution and
// is scientifically meaningless for every Jones parameter that uses this
// family), or an out-of-domain u (see normal_ppf(double)).
double lognormal_gm_gsd_ppf(double u, double gm, double gsd);

// Base-10 lognormal inverse CDF parameterized directly by log10-scale mu and
// sigma (Jones Table 1 GVL style, e.g. mu=7.0, sigma=1.4 -- already in
// log10 RNA-copies-per-mL units, not a geometric mean/GSD pair):
// pow(10, normal_ppf(u) * sigma + mu). Throws std::invalid_argument for
// nonfinite mu, nonfinite or non-positive sigma, or an out-of-domain u (see
// normal_ppf(double)).
double log10normal_ppf(double u, double mu, double sigma);

// Affine uniform inverse CDF: min + (max - min) * u. u in the closed
// interval [0, 1] (see the bounded-transform endpoint policy above). Throws
// std::invalid_argument for nonfinite or out-of-[0,1] u, nonfinite min/max,
// or min >= max.
double uniform_ppf(double u, double min, double max);

// Regularized incomplete beta function I_x(a, b), evaluated via a bounded
// Lentz continued fraction (Numerical Recipes betacf form) on whichever of
// I_x(a,b)/I_{1-x}(b,a) has the faster-converging series, selected by the
// standard x < (a+1)/(a+b+2) split. Throws std::invalid_argument for
// nonfinite or non-positive a/b, or nonfinite/out-of-[0,1] x. Throws
// std::runtime_error, naming the offending (x, a, b), if the continued
// fraction fails to converge within its fixed 200-iteration bound.
double regularized_incomplete_beta(double x, double a, double b);

// Scaled beta inverse CDF: min + (max - min) * beta_standard_ppf(u, a, b),
// where beta_standard_ppf is found by bisection (a fixed 200-iteration bound,
// stopping early once the bracket width is <= 1e-15) on
// regularized_incomplete_beta. u in the closed interval [0, 1] (see the
// bounded-transform endpoint policy above; u = 0 and u = 1 return min and
// max directly without consulting the continued fraction). Throws
// std::invalid_argument for nonfinite or out-of-[0,1] u, nonfinite
// non-positive a/b, nonfinite min/max, or min >= max. May propagate
// regularized_incomplete_beta's std::runtime_error for a pathological (a, b)
// that fails to converge.
double beta_ppf(double u, double a, double b, double min, double max);

// Binomial inverse CDF (percent-point function): the smallest integer k in
// [0, n] such that the binomial(n, p) CDF at k is >= u. u in the half-open
// interval (0, 1] (see the binomial endpoint policy above). n >= 0; p in the
// closed interval [0, 1]. Throws std::invalid_argument for nonfinite or
// out-of-(0,1] u, negative n, or nonfinite/out-of-[0,1] p.
//
// Implemented as a mode-centered PMF recurrence: the PMF is computed exactly
// once, at the mode floor((n+1)*p), via lgamma/log/log1p (so it never starts
// from a term that could underflow), then filled outward in both directions
// by the ratio recurrences pmf(k+1) = pmf(k) * (n-k) * p / ((k+1) * (1-p))
// and pmf(k-1) = pmf(k) * k * (1-p) / ((n-k+1) * p), stopping a direction
// early once a term underflows to exactly 0.0 (the remaining tail is then
// negligible). The n+1 terms are summed, and then scanned cumulatively, with
// Neumaier-compensated summation to control rounding error; the returned k
// is the smallest one whose normalized compensated cumulative sum is >= u.
// This is bounded by exactly n+1 evaluated terms -- an explicit, finite
// iteration bound, not an iterative approximation that could fail to
// converge -- and (unlike a forward-from-zero scan) does not silently
// return n for every u whenever pmf(0) = (1-p)^n underflows to 0.0 (e.g.
// n=400, p=0.9, where 0.1^400 underflows double precision). The n -> size_t
// widening and the (n+1)-element buffer allocation are checked against
// vector<double>::max_size() and rejected with std::runtime_error before any
// allocation if n is too large to support. If floating-point rounding of
// the normalized cumulative sum leaves it short of u even at k = n (whose
// true CDF is exactly 1), the function returns n as an explicit, documented
// fallback rather than continuing indefinitely or returning an out-of-range
// value; this is the only case where a returned k's CDF may differ from u
// by more than ordinary per-term rounding (MODEL_SPECIFICATION.md section
// 11.5 binomial tolerance row). p = 0 and p = 1 are short-circuited to
// return 0 and n exactly, without consulting the recurrence.
long long binomial_ppf(double u, long long n, double p);

}  // namespace distributions
}  // namespace ecai_cpp2

#endif  // ECAI_CPP2_DISTRIBUTIONS_H
