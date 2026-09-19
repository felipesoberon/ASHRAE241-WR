// Stage 3 tests for lhs: exact-size one-dimensional LHS engine, independent
// of any distribution transform, per docs/MODEL_SPECIFICATION.md section
// 11.2 (RNG, open-interval conversion, unbiased bounded-integer generation,
// Fisher-Yates shuffle) and IMPLEMENTATION_PLAN.md Section 7, Stage 3.
//
// All known-answer values below (raw mt19937_64 outputs, open-interval
// conversions, bounded-integer draws, shuffle results, and full
// generate_lhs_1d vectors) were independently derived by hand-assembling
// the exact frozen algorithm from section 11.2 in a throwaway reference
// program and running it against std::mt19937_64, then hard-coded here —
// not copied from lhs.cpp after the fact.

#include "lhs.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef ECAI_CPP2_LHS_SOURCE_PATH
#error "ECAI_CPP2_LHS_SOURCE_PATH must be defined by the build system"
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

bool approx(double a, double b, double tol = 1e-9) {
    return std::fabs(a - b) <= tol;
}

// Every stratum floor(value * n) visited exactly once, for values in
// arbitrary (post-shuffle) order.
bool covers_every_stratum_exactly_once(const std::vector<double>& values, std::size_t n) {
    if (values.size() != n) return false;
    std::set<std::size_t> seen;
    for (double v : values) {
        const double scaled = v * static_cast<double>(n);
        const auto stratum = static_cast<std::size_t>(std::floor(scaled));
        if (stratum >= n) return false;
        if (!seen.insert(stratum).second) return false;  // duplicate stratum
    }
    return seen.size() == n;
}

std::vector<std::size_t> argsort(const std::vector<double>& values) {
    std::vector<std::size_t> indices(values.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&values](std::size_t a, std::size_t b) { return values[a] < values[b]; });
    return indices;
}

bool all_distinct(const std::vector<double>& values) {
    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    for (std::size_t i = 1; i < sorted.size(); ++i) {
        if (sorted[i] == sorted[i - 1]) return false;
    }
    return true;
}

}  // namespace

int main() {
    // ---------------------------------------------------------------
    // open_interval_from_raw: known-answer conversion for min/max synthetic
    // uint64 values, and the smallest/largest RNG integer mappings.
    // ---------------------------------------------------------------
    {
        const double u_min = ecai_cpp2::lhs::open_interval_from_raw(0ULL);
        check(approx(u_min, 5.5511151231257827021e-17, 1e-30),
              "open_interval_from_raw(0) matches the known-answer minimum (0.5/2^53)");
        check(u_min > 0.0, "open_interval_from_raw(0) is strictly greater than 0.0");

        const double u_max =
            ecai_cpp2::lhs::open_interval_from_raw(std::numeric_limits<std::uint64_t>::max());
        // The literal section-11.2 formula (m+0.5)/2^53 with m = 2^53-1 hits
        // round-half-to-even in `m + 0.5`, producing exactly 2^53 and hence
        // u == 1.0 for this raw value unless deterministically guarded; this
        // known-answer value asserts the guard is in place.
        check(approx(u_max, 0.99999999999999988898, 1e-17),
              "open_interval_from_raw(UINT64_MAX) matches the known-answer, "
              "round-half-to-even-guarded maximum");
        check(u_max < 1.0,
              "open_interval_from_raw(UINT64_MAX) is strictly less than 1.0 "
              "(guards the round-half-to-even edge case at m = 2^53-1)");

        // A raw value one m-step below the maximum does not need guarding
        // and should match the unguarded formula directly.
        const std::uint64_t m_second_highest = (1ULL << 53) - 2;
        const std::uint64_t raw_for_m_second_highest = m_second_highest << 11;
        const double u_second_highest =
            ecai_cpp2::lhs::open_interval_from_raw(raw_for_m_second_highest);
        check(approx(u_second_highest, 0.99999999999999977796, 1e-17),
              "open_interval_from_raw at m=2^53-2 matches the known-answer unguarded value");
    }

    // ---------------------------------------------------------------
    // std::mt19937_64 sequence-derived vectors: reference raw outputs for
    // seed 42, independently reproduced here via a fresh engine, and the
    // conversions applied to them.
    // ---------------------------------------------------------------
    {
        std::mt19937_64 reference_engine(42);
        const std::uint64_t expected_raw[6] = {
            13930160852258120406ULL, 11788048577503494824ULL, 13874630024467741450ULL,
            2513787319205155662ULL,  16662371453428439381ULL, 1735254072534978428ULL,
        };
        const double expected_u[6] = {
            0.7551555329545389661,  0.63903139385469742706, 0.7521452007480267099,
            0.13627268363243710825, 0.9032689664283783948,  0.094068311762837086842,
        };
        for (int i = 0; i < 6; ++i) {
            const std::uint64_t raw = reference_engine();
            check(raw == expected_raw[i],
                  "mt19937_64(42) raw draw " + std::to_string(i) +
                      " matches the known-answer sequence");
            const double u = ecai_cpp2::lhs::open_interval_from_raw(raw);
            check(approx(u, expected_u[i], 1e-15),
                  "open_interval_from_raw of mt19937_64(42) draw " + std::to_string(i) +
                      " matches the known-answer conversion");
        }
    }

    // ---------------------------------------------------------------
    // bounded_uint64: known-answer sequence, plus an independent parallel
    // reference implementation of the exact threshold-rejection algorithm
    // from section 11.2, run against a fresh, identically-seeded engine.
    // ---------------------------------------------------------------
    {
        std::mt19937_64 engine(7);
        const std::uint64_t expected[8] = {0, 0, 3, 1, 1, 3, 4, 3};
        for (int i = 0; i < 8; ++i) {
            const std::uint64_t drawn = ecai_cpp2::lhs::bounded_uint64(engine, 5);
            check(drawn == expected[i],
                  "bounded_uint64(engine, 5) draw " + std::to_string(i) +
                      " matches the known-answer sequence");
            check(drawn < 5, "bounded_uint64(engine, 5) draw " + std::to_string(i) +
                                  " is within [0, bound)");
        }

        // Independent parallel reference: same seed, same bound, algorithm
        // re-derived from the frozen specification text rather than calling
        // the library function under test.
        std::mt19937_64 reference_engine(7);
        bool all_match = true;
        for (int i = 0; i < 8; ++i) {
            const std::uint64_t bound = 5;
            const std::uint64_t threshold = (static_cast<std::uint64_t>(0) - bound) % bound;
            std::uint64_t x;
            do {
                x = reference_engine();
            } while (x < threshold);
            const std::uint64_t reference_value = x % bound;
            if (reference_value != expected[i]) all_match = false;
        }
        check(all_match,
              "independent parallel reference implementation of threshold-rejection "
              "matches the known-answer sequence");
    }

    // ---------------------------------------------------------------
    // bounded_uint64: reject impossible bounds (bound == 0).
    // ---------------------------------------------------------------
    {
        std::mt19937_64 engine(1);
        bool threw = false;
        try {
            ecai_cpp2::lhs::bounded_uint64(engine, 0);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "bounded_uint64 rejects an impossible bound of 0 with std::invalid_argument");
    }

    // ---------------------------------------------------------------
    // fisher_yates_shuffle: frozen known-answer vector, seed 99, n=5.
    // ---------------------------------------------------------------
    {
        std::mt19937_64 engine(99);
        std::vector<double> values = {10.0, 20.0, 30.0, 40.0, 50.0};
        ecai_cpp2::lhs::fisher_yates_shuffle(values, engine);
        const std::vector<double> expected = {30.0, 10.0, 50.0, 40.0, 20.0};
        check(values == expected,
              "fisher_yates_shuffle(seed=99, n=5) matches the frozen known-answer vector");
    }
    {
        // n=1 and n=0 must be no-ops (nothing to swap), never crash or loop.
        std::mt19937_64 engine(1);
        std::vector<double> one = {7.0};
        ecai_cpp2::lhs::fisher_yates_shuffle(one, engine);
        check(one.size() == 1 && one[0] == 7.0, "fisher_yates_shuffle on a size-1 vector is a no-op");

        std::vector<double> empty;
        ecai_cpp2::lhs::fisher_yates_shuffle(empty, engine);
        check(empty.empty(), "fisher_yates_shuffle on an empty vector is a no-op");
    }

    // ---------------------------------------------------------------
    // generate_lhs_1d: reject size zero.
    // ---------------------------------------------------------------
    {
        bool threw = false;
        try {
            ecai_cpp2::lhs::generate_lhs_1d(0, 123);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "generate_lhs_1d(0, seed) throws std::invalid_argument");
    }

    // ---------------------------------------------------------------
    // generate_lhs_1d: known-answer full vectors for small sizes.
    // ---------------------------------------------------------------
    {
        const std::vector<double> result = ecai_cpp2::lhs::generate_lhs_1d(3, 123);
        const std::vector<double> expected = {
            0.97942836939258959283,
            0.10440005955949024918,
            0.5186597064649527411,
        };
        check(result.size() == expected.size(), "generate_lhs_1d(3, 123) returns exactly 3 values");
        bool all_close = result.size() == expected.size();
        for (std::size_t i = 0; all_close && i < expected.size(); ++i) {
            all_close = approx(result[i], expected[i], 1e-12);
        }
        check(all_close, "generate_lhs_1d(3, 123) matches the known-answer vector exactly");
    }
    {
        const std::vector<double> result = ecai_cpp2::lhs::generate_lhs_1d(2, 7);
        const std::vector<double> expected = {0.97465060144632209571, 0.37719265207642904691};
        bool all_close = result.size() == expected.size();
        for (std::size_t i = 0; all_close && i < expected.size(); ++i) {
            all_close = approx(result[i], expected[i], 1e-12);
        }
        check(all_close, "generate_lhs_1d(2, 7) matches the known-answer vector");
    }
    {
        const std::vector<double> result = ecai_cpp2::lhs::generate_lhs_1d(1, 555);
        check(result.size() == 1, "generate_lhs_1d(1, 555) returns exactly 1 value");
        check(approx(result[0], 0.85485072223962199622, 1e-12),
              "generate_lhs_1d(1, 555) matches the known-answer single value");
    }

    // ---------------------------------------------------------------
    // generate_lhs_1d: sizes 1, 2, 3, odd sizes, and 10,000 — structural
    // invariants (length, finiteness, strict-(0,1), exact stratum coverage,
    // no duplicates) without any flaky statistical/correlation claim.
    // ---------------------------------------------------------------
    {
        const std::vector<std::size_t> sizes = {1, 2, 3, 5, 7, 9999, 10000};
        for (std::size_t n : sizes) {
            const std::vector<double> result = ecai_cpp2::lhs::generate_lhs_1d(n, 2000 + n);
            const std::string tag = "generate_lhs_1d(n=" + std::to_string(n) + ")";

            check(result.size() == n, tag + " returns exactly n values");

            bool all_finite_strict = true;
            for (double v : result) {
                if (!std::isfinite(v) || !(v > 0.0) || !(v < 1.0)) {
                    all_finite_strict = false;
                    break;
                }
            }
            check(all_finite_strict, tag + " every value is finite and strictly in (0, 1)");

            check(covers_every_stratum_exactly_once(result, n),
                  tag + " floor(value * n) covers every stratum exactly once");

            check(all_distinct(result), tag + " contains no duplicate values");
        }
    }

    // ---------------------------------------------------------------
    // generate_lhs_1d: endpoint exclusion — no value is exactly 0.0 or 1.0,
    // explicitly, across the full 10,000-size list (covers the boundary
    // strata i=0 and i=n-1 where rounding risk is highest).
    // ---------------------------------------------------------------
    {
        const std::vector<double> result = ecai_cpp2::lhs::generate_lhs_1d(10000, 4242);
        bool any_endpoint = false;
        for (double v : result) {
            if (v == 0.0 || v == 1.0) any_endpoint = true;
        }
        check(!any_endpoint, "generate_lhs_1d(10000, ...) never yields exactly 0.0 or 1.0");
    }

    // ---------------------------------------------------------------
    // generate_lhs_1d: same-seed identity — same (size, seed) reproduces an
    // identical vector, bit-for-bit.
    // ---------------------------------------------------------------
    {
        const std::vector<double> a = ecai_cpp2::lhs::generate_lhs_1d(50, 909090);
        const std::vector<double> b = ecai_cpp2::lhs::generate_lhs_1d(50, 909090);
        check(a == b, "generate_lhs_1d(50, 909090) called twice yields bit-for-bit identical vectors");
    }

    // ---------------------------------------------------------------
    // generate_lhs_1d: different-seed vector and rank difference, and no
    // duplicate sequence across named test streams (seed 2026 vs 2027).
    // ---------------------------------------------------------------
    {
        const std::vector<double> a = ecai_cpp2::lhs::generate_lhs_1d(5, 2026);
        const std::vector<double> b = ecai_cpp2::lhs::generate_lhs_1d(5, 2027);
        const std::vector<double> expected_a = {
            0.8511219367491099641, 0.49691936895701949162, 0.75183961652627206984,
            0.063499227159712359758, 0.330871453824236883,
        };
        const std::vector<double> expected_b = {
            0.4414859034314400632, 0.13699864410944159743, 0.83880108793292962321,
            0.73137464969212562238, 0.262894568130376316,
        };
        bool a_matches = a.size() == expected_a.size();
        for (std::size_t i = 0; a_matches && i < expected_a.size(); ++i) {
            a_matches = approx(a[i], expected_a[i], 1e-12);
        }
        bool b_matches = b.size() == expected_b.size();
        for (std::size_t i = 0; b_matches && i < expected_b.size(); ++i) {
            b_matches = approx(b[i], expected_b[i], 1e-12);
        }
        check(a_matches, "generate_lhs_1d(5, 2026) matches its known-answer vector");
        check(b_matches, "generate_lhs_1d(5, 2027) matches its known-answer vector");
        check(a != b, "generate_lhs_1d(5, 2026) and (5, 2027) produce different value vectors");
        check(argsort(a) != argsort(b),
              "generate_lhs_1d(5, 2026) and (5, 2027) produce different rank permutations");
    }

    // ---------------------------------------------------------------
    // generate_lhs_1d: zero rejection also confirmed for a range of
    // representative sizes/seeds is unnecessary (constructor already
    // rejects size 0 unconditionally regardless of seed) — confirm that
    // claim directly across a couple of distinct seeds.
    // ---------------------------------------------------------------
    {
        for (std::uint64_t seed : {0ULL, 1ULL, 999999ULL}) {
            bool threw = false;
            try {
                ecai_cpp2::lhs::generate_lhs_1d(0, seed);
            } catch (const std::invalid_argument&) {
                threw = true;
            }
            check(threw, "generate_lhs_1d(0, seed=" + std::to_string(seed) +
                              ") throws std::invalid_argument regardless of seed");
        }
    }

    // ---------------------------------------------------------------
    // Source-level compliance: no std::uniform_real_distribution or
    // std::shuffle anywhere in lhs.cpp (frozen-algorithm requirement,
    // section 11.2).
    // ---------------------------------------------------------------
    {
        std::ifstream source_file(ECAI_CPP2_LHS_SOURCE_PATH);
        check(source_file.good(), "lhs.cpp source file is readable for compliance scan");
        std::stringstream buffer;
        buffer << source_file.rdbuf();
        const std::string source = buffer.str();
        check(source.find("uniform_real_distribution") == std::string::npos,
              "lhs.cpp does not use std::uniform_real_distribution");
        check(source.find("std::shuffle") == std::string::npos,
              "lhs.cpp does not use std::shuffle");
    }

    if (g_failures > 0) {
        std::cerr << g_failures << " test(s) failed.\n";
        return 1;
    }
    std::cout << "All tests passed.\n";
    return 0;
}
