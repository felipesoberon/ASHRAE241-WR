// Stage 3: exact-size one-dimensional LHS engine, independent of any
// physical distribution transform.
//
// Source of truth: docs/MODEL_SPECIFICATION.md section 11.2.
// IMPLEMENTATION_PLAN.md Section 7, Stage 3.

#include "lhs.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace ecai_cpp2 {
namespace lhs {

double open_interval_from_raw(std::uint64_t raw) {
    const std::uint64_t m = raw >> 11;  // high 53 bits, 0 <= m < 2^53
    constexpr double kTwoPow53 = 9007199254740992.0;  // 2^53
    double u = (static_cast<double>(m) + 0.5) / kTwoPow53;

    // At m = 2^53-1, `m + 0.5` (9007199254740991.5) is not exactly
    // representable in a double and round-half-to-even rounds it up to
    // exactly 2^53, which would otherwise make u == 1.0 for the top 2^11
    // raw values. Deterministically guard both endpoints rather than ever
    // returning a value outside the open interval (0, 1).
    if (!(u < 1.0)) {
        u = std::nextafter(1.0, 0.0);
    }
    if (!(u > 0.0)) {
        u = std::nextafter(0.0, 1.0);
    }
    return u;
}

std::uint64_t bounded_uint64(std::mt19937_64& engine, std::uint64_t bound) {
    if (bound == 0) {
        throw std::invalid_argument("lhs::bounded_uint64: bound must be nonzero");
    }
    const std::uint64_t threshold = (static_cast<std::uint64_t>(0) - bound) % bound;
    std::uint64_t x;
    do {
        x = engine();
    } while (x < threshold);
    return x % bound;
}

void fisher_yates_shuffle(std::vector<double>& values, std::mt19937_64& engine) {
    const std::size_t n = values.size();
    for (std::size_t i = n; i > 1; --i) {
        const std::size_t idx = i - 1;
        const std::uint64_t j = bounded_uint64(engine, static_cast<std::uint64_t>(idx) + 1);
        std::swap(values[idx], values[static_cast<std::size_t>(j)]);
    }
}

std::vector<double> generate_lhs_1d(std::size_t size, std::uint64_t seed) {
    if (size == 0) {
        throw std::invalid_argument("lhs::generate_lhs_1d: size must be nonzero");
    }

    std::mt19937_64 engine(seed);
    const double count = static_cast<double>(size);

    std::vector<double> points;
    points.reserve(size);
    for (std::size_t i = 0; i < size; ++i) {
        const double u = open_interval_from_raw(engine());
        const double lower = static_cast<double>(i) / count;
        const double upper = static_cast<double>(i + 1) / count;
        double point = lower + u * (upper - lower);

        // Combining the stratum bounds with u can round the result onto (or
        // past) the stratum boundary — most notably the final stratum's
        // upper bound of exactly 1.0 — because the true deviation from the
        // boundary can be smaller than a double's ULP at that magnitude.
        // Nudge deterministically back inside the stratum rather than ever
        // returning an endpoint.
        if (!(point > lower)) {
            point = std::nextafter(lower, upper);
        }
        if (!(point < upper)) {
            point = std::nextafter(upper, lower);
        }
        points.push_back(point);
    }

    fisher_yates_shuffle(points, engine);
    return points;
}

}  // namespace lhs
}  // namespace ecai_cpp2
