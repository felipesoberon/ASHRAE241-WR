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

#include <vector>

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

// Reusable binomial inverse CDF for one fixed (n, p): the mode-centered PMF
// recurrence and the Neumaier-compensated normalizer and cumulative scan of
// binomial_ppf are evaluated exactly once, in the constructor, and the
// normalized cumulative value binomial_ppf compares against u at each k is
// stored, so ppf(u) reduces to the same >= comparison over the stored
// values. For every valid u, ppf(u) == binomial_ppf(u, n, p) exactly,
// including the documented tolerated CDF-boundary ties and the k = n
// fallback. Intended for callers that transform many u values through the
// same (n, p) (Stage 6's per-replicate infector-count list), where
// recomputing the (n+1)-term recurrence per u is pure repeated work.
//
// The constructor validates (n, p) with binomial_ppf's parameter contract
// (std::invalid_argument for negative n or nonfinite/out-of-[0,1] p) and
// applies the same vector<double>::max_size() guard (std::runtime_error)
// before allocating. It transiently owns one (n+1)-double PMF buffer while
// building the retained (n+1)-double CDF buffer, so its peak working set is
// two (n+1)-double buffers; the PMF buffer is released when the constructor
// returns. n = 0, p = 0 and p = 1 allocate nothing and are answered by
// ppf() with binomial_ppf's short-circuit conventions.
//
// ppf(u) validates u exactly as binomial_ppf does (half-open (0, 1],
// std::invalid_argument otherwise) and is const and reentrant.
class BinomialPpfTable {
  public:
    BinomialPpfTable(long long n, double p);

    long long ppf(double u) const;

    long long n() const { return n_; }
    double p() const { return p_; }

    // Zero-infector mass P(X = 0) for this table's (n, p), exposed for
    // risk-averse infector-policy callers that must know it before consuming
    // any LHS draws. Degenerate cases are answered without touching cdf_
    // (which is empty for all of them): n = 0 always returns 1.0 (X is
    // identically 0 regardless of p), p <= 0 returns 1.0 (all mass at 0),
    // and p >= 1.0 (with n >= 1) returns 0.0 (all mass at n). Otherwise
    // returns the stored cdf_[0], identical to the value ppf() itself
    // compares against u at k = 0.
    double zero_mass() const;

    // Conditional binomial ppf given X >= 1 ("at least one infector" out of
    // n = I0 candidates each with success probability p): the smallest k in
    // [1, n] such that the retained conditional-tail cumulative mass
    // T[k] = sum(pmf[1..k]) / sum(pmf[1..n]) is >= u. T is built once, in
    // the constructor, by compensated-summing the same PMF terms cdf_ is
    // built from -- never by subtracting cdf_[k] - zero_mass, which would
    // cancel catastrophically whenever zero_mass is close to 1. Requires
    // n (I0) >= 1 and p strictly positive; throws std::invalid_argument
    // otherwise (this is a stricter domain than the table's own community
    // (n, p) contract, so a table built for community use with n = 0 or
    // p = 0 correctly rejects a conditional query). u must be finite and in
    // the half-open interval (0, 1], validated exactly as ppf()'s u; throws
    // std::invalid_argument otherwise.
    //
    // p = 1 short-circuits before any table access -- the conditional
    // distribution is a point mass at n for every valid u -- and returns n
    // unconditionally; zero_mass() reports 0.0 in that branch. If
    // zero_mass() rounds to exactly 1.0 in binary64 for this representable
    // positive p (the conditional tail is too small relative to pmf[0] to
    // resolve numerically), an explicit limiting branch returns 1 rather
    // than searching an unresolvable tail or relying on a coincidental
    // 0 >= 0 comparison to produce the same answer. Otherwise this is
    // bounded by exactly n evaluated terms, the same explicit finite
    // iteration bound as ppf().
    long long conditional_ppf_at_least_one(double u) const;

  private:
    long long n_;
    double p_;
    std::vector<double> cdf_;
    std::vector<double> conditional_tail_;
};

}  // namespace distributions
}  // namespace ecai_cpp2

#endif  // ECAI_CPP2_DISTRIBUTIONS_H
