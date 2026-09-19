// Stage 3: exact-size one-dimensional Latin Hypercube Sample (LHS) engine,
// independent of any physical distribution transform.
//
// Source of truth: docs/MODEL_SPECIFICATION.md section 11.2 (RNG, open-
// interval conversion, unbiased bounded-integer generation, Fisher-Yates
// shuffle) and section 11.5 (validation tolerances).
// IMPLEMENTATION_PLAN.md Section 7, Stage 3.

#ifndef ECAI_CPP2_LHS_H
#define ECAI_CPP2_LHS_H

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

namespace ecai_cpp2 {
namespace lhs {

// Frozen open-interval integer -> double conversion (MODEL_SPECIFICATION.md
// section 11.2): takes the high 53 bits of `raw` (m = raw >> 11, so
// 0 <= m < 2^53) and maps to u = (m + 0.5) / 2^53. Guarantees 0 < u < 1
// strictly for every possible `raw` value, using the full 53-bit mantissa
// precision of an IEEE-754 double.
double open_interval_from_raw(std::uint64_t raw);

// Unbiased bounded-integer generation in [0, bound) via threshold-rejection
// (MODEL_SPECIFICATION.md section 11.2): threshold = uint64_t(0 - bound) %
// bound using defined unsigned wraparound; repeatedly draw raw uint64_t
// values from `engine` until one is >= threshold; return that value % bound.
// `bound` must be nonzero; throws std::invalid_argument otherwise.
std::uint64_t bounded_uint64(std::mt19937_64& engine, std::uint64_t bound);

// Deterministic Fisher-Yates (Durstenfeld) shuffle, in place: for i from
// values.size()-1 down to 1, draw j = bounded_uint64(engine, i + 1)
// (uniformly in [0, i]) and swap values[i] with values[j]. Fully
// deterministic given `engine`'s state.
void fisher_yates_shuffle(std::vector<double>& values, std::mt19937_64& engine);

// Exact-size, one-dimensional LHS, independent of any distribution
// transform. For size `n`, generates one jittered point strictly inside
// every stratum (i/n, (i+1)/n) for i in [0, n), using the open-interval
// conversion above to draw the jitter, then applies the frozen deterministic
// Fisher-Yates shuffle (continuing to draw from the same engine) so the
// returned order is not stratum-ordered. Every returned value is finite and
// strictly in (0, 1); floating-point rounding that would otherwise place a
// jittered point on (or past) its stratum boundary — including the final
// stratum's upper bound of exactly 1.0 — is deterministically nudged back
// inside the stratum with std::nextafter rather than ever returning 0.0 or
// 1.0. `size` must be nonzero; throws std::invalid_argument otherwise.
std::vector<double> generate_lhs_1d(std::size_t size, std::uint64_t seed);

}  // namespace lhs
}  // namespace ecai_cpp2

#endif  // ECAI_CPP2_LHS_H
