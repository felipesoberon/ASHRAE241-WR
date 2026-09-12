// Stage 5 tests for distributions: inverse CDFs transforming a single
// open-unit-interval LHS draw into a Jones et al. parameter distribution,
// per docs/MODEL_SPECIFICATION.md section 5 (Table 1) and section 11.5
// (validation tolerances), and IMPLEMENTATION_PLAN.md Section 7, Stage 5.
//
// Fixture provenance: every hard-coded reference value below (normal,
// location-scale normal, lognormal-from-GM/GSD, log10normal, uniform, beta,
// and binomial) was independently computed with SciPy 1.17.1 / Python 3.11.15
// (scipy.stats.norm/lognorm/beta/binom .ppf/.cdf) in a standalone script
// outside this repository (not by calling any production helper in this
// codebase), then hard-coded here as literals. The script and its output
// were deleted after transcription; it is not part of this commit.

#include "distributions.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

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

bool close(double actual, double expected, double abs_tol, double rel_tol) {
    if (!std::isfinite(actual) || !std::isfinite(expected)) return false;
    double diff = std::fabs(actual - expected);
    if (diff <= abs_tol) return true;
    double rel = diff / std::max(std::fabs(expected), 1e-300);
    return rel <= rel_tol;
}

template <typename Fn>
bool throws_invalid_argument(Fn&& fn) {
    try {
        fn();
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

}  // namespace

int main() {
    using namespace ecai_cpp2::distributions;

    const double kNormalAbs = 1e-9;
    const double kNormalRel = 1e-9;
    const double kUniformAbs = 1e-14;
    const double kBetaAbs = 1e-9;
    const double kBetaRel = 1e-8;

    // ================================================================
    // Standard normal ppf: dense grid, SciPy scipy.stats.norm.ppf.
    // ================================================================
    {
        struct Case { double u; double expected; };
        const Case cases[] = {
            {1e-08, -5.612001244174789},
            {1e-06, -4.753424308822899},
            {0.0001, -3.7190164854556804},
            {0.001, -3.090232306167813},
            {0.01, -2.3263478740408408},
            {0.025, -1.9599639845400545},
            {0.05, -1.6448536269514729},
            {0.1, -1.2815515655446004},
            {0.2, -0.8416212335729142},
            {0.3, -0.5244005127080409},
            {0.4, -0.2533471031357997},
            {0.5, 0.0},
            {0.6, 0.2533471031357997},
            {0.7, 0.5244005127080407},
            {0.8, 0.8416212335729143},
            {0.9, 1.2815515655446004},
            {0.95, 1.6448536269514722},
            {0.975, 1.959963984540054},
            {0.99, 2.3263478740408408},
            {0.999, 3.090232306167813},
            {0.9999, 3.719016485455709},
            {0.999999, 4.753424308817087},
            {0.99999999, 5.612001243305505},
        };
        for (const auto& c : cases) {
            double actual = normal_ppf(c.u);
            check(close(actual, c.expected, kNormalAbs, kNormalRel),
                  "normal_ppf(" + std::to_string(c.u) + ") matches SciPy norm.ppf");
        }
        // Symmetry: normal_ppf(u) == -normal_ppf(1-u) for a dense grid.
        const double sym_us[] = {1e-8, 1e-6, 1e-4, 0.001, 0.01, 0.1, 0.3, 0.49};
        for (double u : sym_us) {
            check(close(normal_ppf(u), -normal_ppf(1.0 - u), 1e-8, 1e-8),
                  "normal_ppf(" + std::to_string(u) + ") == -normal_ppf(1-u) (symmetry)");
        }
        // Monotonicity, dense near both tails.
        double prev = normal_ppf(1e-12);
        for (int i = 1; i <= 2000; ++i) {
            double u = 1e-12 + (0.02 - 1e-12) * (static_cast<double>(i) / 2000.0);
            double cur = normal_ppf(u);
            check(cur > prev, "normal_ppf strictly increasing near lower tail at i=" +
                                   std::to_string(i));
            prev = cur;
        }
        prev = normal_ppf(0.98);
        for (int i = 1; i <= 2000; ++i) {
            double u = 0.98 + (1.0 - 1e-12 - 0.98) * (static_cast<double>(i) / 2000.0);
            double cur = normal_ppf(u);
            check(cur > prev, "normal_ppf strictly increasing near upper tail at i=" +
                                   std::to_string(i));
            prev = cur;
        }
    }

    // ================================================================
    // Standard normal ppf: invalid inputs.
    // ================================================================
    {
        check(throws_invalid_argument([] { normal_ppf(0.0); }),
              "normal_ppf(0.0) throws std::invalid_argument");
        check(throws_invalid_argument([] { normal_ppf(1.0); }),
              "normal_ppf(1.0) throws std::invalid_argument");
        check(throws_invalid_argument([] { normal_ppf(-0.1); }),
              "normal_ppf(-0.1) throws std::invalid_argument");
        check(throws_invalid_argument([] { normal_ppf(1.1); }),
              "normal_ppf(1.1) throws std::invalid_argument");
        check(throws_invalid_argument([] { normal_ppf(std::numeric_limits<double>::quiet_NaN()); }),
              "normal_ppf(NaN) throws std::invalid_argument");
        check(throws_invalid_argument(
                  [] { normal_ppf(std::numeric_limits<double>::infinity()); }),
              "normal_ppf(+inf) throws std::invalid_argument");
        check(throws_invalid_argument(
                  [] { normal_ppf(-std::numeric_limits<double>::infinity()); }),
              "normal_ppf(-inf) throws std::invalid_argument");
    }

    // ================================================================
    // Location-scale normal ppf: mu=10, sigma=2.5, SciPy norm.ppf(u, loc, scale).
    // ================================================================
    {
        struct Case { double u; double expected; };
        const Case cases[] = {
            {1e-06, -1.8835607720572476},
            {0.001, 2.274419234580467},
            {0.01, 4.184130314897898},
            {0.1, 6.796121086138499},
            {0.3, 8.688998718229898},
            {0.5, 10.0},
            {0.7, 11.311001281770102},
            {0.9, 13.203878913861502},
            {0.99, 15.815869685102102},
            {0.999, 17.725580765419533},
            {0.999999, 21.88356077204272},
        };
        for (const auto& c : cases) {
            double actual = normal_ppf(c.u, 10.0, 2.5);
            check(close(actual, c.expected, kNormalAbs, kNormalRel),
                  "normal_ppf(" + std::to_string(c.u) + ", 10, 2.5) matches SciPy");
        }
        check(close(normal_ppf(0.5, 10.0, 2.5), 10.0, kNormalAbs, kNormalRel),
              "normal_ppf(0.5, mu, sigma) == mu");

        check(throws_invalid_argument([] { normal_ppf(0.5, 0.0, 0.0); }),
              "normal_ppf with sigma=0 throws std::invalid_argument");
        check(throws_invalid_argument([] { normal_ppf(0.5, 0.0, -1.0); }),
              "normal_ppf with negative sigma throws std::invalid_argument");
        check(throws_invalid_argument(
                  [] { normal_ppf(0.5, std::numeric_limits<double>::quiet_NaN(), 1.0); }),
              "normal_ppf with nonfinite mu throws std::invalid_argument");
        check(throws_invalid_argument([] { normal_ppf(0.0, 10.0, 2.5); }),
              "normal_ppf(0.0, mu, sigma) throws std::invalid_argument (open-u policy inherited)");
    }

    // ================================================================
    // Natural-lognormal from GM/GSD: lambda (Jones biological decay rate),
    // GM=0.52, GSD=1.9, SciPy lognorm.ppf(u, s=ln(GSD), scale=GM).
    // ================================================================
    {
        struct Case { double u; double expected; };
        const Case cases[] = {
            {1e-06, 0.02460193125926831},
            {0.0001, 0.04778747747412529},
            {0.001, 0.07154675819186283},
            {0.01, 0.11682223121570454},
            {0.05, 0.18092328812826466},
            {0.1, 0.22843688526552214},
            {0.3, 0.3713853878012194},
            {0.5, 0.52},
            {0.7, 0.7280846497512958},
            {0.9, 1.1836967558268985},
            {0.95, 1.4945560784209337},
            {0.99, 2.31462793670431},
            {0.999, 3.7793466375497244},
            {0.9999, 5.6583861357070795},
            {0.999999, 10.991007053445182},
        };
        for (const auto& c : cases) {
            double actual = lognormal_gm_gsd_ppf(c.u, 0.52, 1.9);
            check(close(actual, c.expected, kNormalAbs, kNormalRel),
                  "lognormal_gm_gsd_ppf(" + std::to_string(c.u) +
                      ", GM=0.52, GSD=1.9) matches SciPy");
        }
        check(close(lognormal_gm_gsd_ppf(0.5, 0.52, 1.9), 0.52, kNormalAbs, kNormalRel),
              "lognormal_gm_gsd_ppf(0.5, gm, gsd) == gm (median)");

        // Monotonicity, dense near both tails.
        double prev = lognormal_gm_gsd_ppf(1e-9, 0.52, 1.9);
        for (int i = 1; i <= 1000; ++i) {
            double u = 1e-9 + (0.01 - 1e-9) * (static_cast<double>(i) / 1000.0);
            double cur = lognormal_gm_gsd_ppf(u, 0.52, 1.9);
            check(cur > prev, "lognormal_gm_gsd_ppf strictly increasing near lower tail, i=" +
                                   std::to_string(i));
            prev = cur;
        }

        check(throws_invalid_argument([] { lognormal_gm_gsd_ppf(0.5, 0.0, 1.9); }),
              "lognormal_gm_gsd_ppf with gm=0 throws std::invalid_argument");
        check(throws_invalid_argument([] { lognormal_gm_gsd_ppf(0.5, -0.52, 1.9); }),
              "lognormal_gm_gsd_ppf with negative gm throws std::invalid_argument");
        check(throws_invalid_argument([] { lognormal_gm_gsd_ppf(0.5, 0.52, 1.0); }),
              "lognormal_gm_gsd_ppf with gsd=1.0 throws std::invalid_argument (degenerate)");
        check(throws_invalid_argument([] { lognormal_gm_gsd_ppf(0.5, 0.52, 0.9); }),
              "lognormal_gm_gsd_ppf with gsd<1 throws std::invalid_argument");
        check(throws_invalid_argument([] { lognormal_gm_gsd_ppf(0.0, 0.52, 1.9); }),
              "lognormal_gm_gsd_ppf(0.0, ...) throws std::invalid_argument");
        check(throws_invalid_argument([] { lognormal_gm_gsd_ppf(1.0, 0.52, 1.9); }),
              "lognormal_gm_gsd_ppf(1.0, ...) throws std::invalid_argument");
    }

    // ================================================================
    // Base-10 lognormal direct (GVL: Jones genomic viral load), mu=7.0,
    // sigma=1.4 (already log10-scale parameters, no ln() transform).
    // ================================================================
    {
        struct Case { double u; double expected; };
        const Case cases[] = {
            {1e-06, 2.2141445364803727},
            {0.0001, 62.14081149033561},
            {0.001, 471.7096613494902},
            {0.01, 5534.940751334816},
            {0.05, 49797.199811763676},
            {0.1, 160630.42492953505},
            {0.3, 1844332.7673042489},
            {0.5, 10000000.0},
            {0.7, 54220150.383254334},
            {0.9, 622547067.5550272},
            {0.95, 2008145043.8579967},
            {0.99, 18067040731.35448},
            {0.999, 211994809930.14917},
            {0.9999, 1609248376416.2239},
            {0.999999, 45164169886166.39},
        };
        for (const auto& c : cases) {
            double actual = log10normal_ppf(c.u, 7.0, 1.4);
            check(close(actual, c.expected, kNormalAbs, kNormalRel),
                  "log10normal_ppf(" + std::to_string(c.u) + ", mu=7, sigma=1.4) matches SciPy");
        }
        check(close(log10normal_ppf(0.5, 7.0, 1.4), 1e7, kNormalAbs, kNormalRel),
              "log10normal_ppf(0.5, mu, sigma) == 10^mu (median)");

        // A base-mistake guard: log10normal_ppf must NOT equal
        // lognormal_gm_gsd_ppf(u, exp(mu), exp(sigma)) -- these are two
        // different parameterizations (log10 direct vs. GM/GSD natural log)
        // and must not be silently interchanged.
        check(!close(log10normal_ppf(0.9, 7.0, 1.4),
                     lognormal_gm_gsd_ppf(0.9, std::exp(7.0), std::exp(1.4)), 1.0, 1e-3),
              "log10normal_ppf is not accidentally equal to a mis-based lognormal_gm_gsd_ppf "
              "call (parameter-order/base guard)");

        check(throws_invalid_argument([] { log10normal_ppf(0.5, 7.0, 0.0); }),
              "log10normal_ppf with sigma=0 throws std::invalid_argument");
        check(throws_invalid_argument([] { log10normal_ppf(0.5, 7.0, -1.0); }),
              "log10normal_ppf with negative sigma throws std::invalid_argument");
        check(throws_invalid_argument([] { log10normal_ppf(0.0, 7.0, 1.4); }),
              "log10normal_ppf(0.0, ...) throws std::invalid_argument");
        check(throws_invalid_argument([] { log10normal_ppf(1.0, 7.0, 1.4); }),
              "log10normal_ppf(1.0, ...) throws std::invalid_argument");
    }

    // ================================================================
    // Affine uniform: RTD (min=0.43, max=0.65) and DK (min=5, max=15).
    // ================================================================
    {
        struct Case { double u; double expected; };
        const Case rtd_cases[] = {
            {0.0, 0.43},        {1e-06, 0.43000022}, {0.001, 0.43022},
            {0.25, 0.485},      {0.5, 0.54},          {0.75, 0.595},
            {0.999, 0.64978},   {1.0, 0.65},
        };
        for (const auto& c : rtd_cases) {
            double actual = uniform_ppf(c.u, 0.43, 0.65);
            check(close(actual, c.expected, kUniformAbs, 0.0),
                  "uniform_ppf(" + std::to_string(c.u) + ", 0.43, 0.65) matches hand calc");
        }
        const Case dk_cases[] = {
            {0.0, 5.0}, {0.1, 6.0}, {0.5, 10.0}, {0.9, 14.0}, {1.0, 15.0},
        };
        for (const auto& c : dk_cases) {
            double actual = uniform_ppf(c.u, 5.0, 15.0);
            check(close(actual, c.expected, kUniformAbs, 0.0),
                  "uniform_ppf(" + std::to_string(c.u) + ", 5, 15) matches hand calc");
        }
        // Endpoints exact.
        check(uniform_ppf(0.0, 0.43, 0.65) == 0.43, "uniform_ppf(0.0, min, max) == min exactly");
        check(uniform_ppf(1.0, 0.43, 0.65) == 0.65, "uniform_ppf(1.0, min, max) == max exactly");

        // Monotonicity, dense grid.
        double prev = uniform_ppf(0.0, 0.42, 0.61);
        for (int i = 1; i <= 1000; ++i) {
            double u = static_cast<double>(i) / 1000.0;
            double cur = uniform_ppf(u, 0.42, 0.61);
            check(cur > prev, "uniform_ppf strictly increasing, i=" + std::to_string(i));
            prev = cur;
        }

        check(throws_invalid_argument([] { uniform_ppf(-1e-9, 0.43, 0.65); }),
              "uniform_ppf(u<0, ...) throws std::invalid_argument");
        check(throws_invalid_argument([] { uniform_ppf(1.0 + 1e-9, 0.43, 0.65); }),
              "uniform_ppf(u>1, ...) throws std::invalid_argument");
        check(throws_invalid_argument([] { uniform_ppf(0.5, 0.65, 0.43); }),
              "uniform_ppf with min > max throws std::invalid_argument");
        check(throws_invalid_argument([] { uniform_ppf(0.5, 0.43, 0.43); }),
              "uniform_ppf with min == max throws std::invalid_argument");
        check(throws_invalid_argument(
                  [] { uniform_ppf(std::numeric_limits<double>::quiet_NaN(), 0.43, 0.65); }),
              "uniform_ppf(NaN, ...) throws std::invalid_argument");
    }

    // ================================================================
    // Regularized incomplete beta: non-convergence diagnostic.
    // ================================================================
    {
        check(throws_invalid_argument([] { regularized_incomplete_beta(0.5, 0.0, 2.0); }),
              "regularized_incomplete_beta with a<=0 throws std::invalid_argument");
        check(throws_invalid_argument([] { regularized_incomplete_beta(0.5, 2.0, 0.0); }),
              "regularized_incomplete_beta with b<=0 throws std::invalid_argument");
        check(throws_invalid_argument([] { regularized_incomplete_beta(-0.1, 2.0, 2.0); }),
              "regularized_incomplete_beta with x<0 throws std::invalid_argument");
        check(throws_invalid_argument([] { regularized_incomplete_beta(1.1, 2.0, 2.0); }),
              "regularized_incomplete_beta with x>1 throws std::invalid_argument");
        check(regularized_incomplete_beta(0.0, 2.0, 2.0) == 0.0,
              "regularized_incomplete_beta(0, a, b) == 0 exactly");
        check(regularized_incomplete_beta(1.0, 2.0, 2.0) == 1.0,
              "regularized_incomplete_beta(1, a, b) == 1 exactly");

        bool threw_runtime = false;
        try {
            regularized_incomplete_beta(0.5, 1e300, 1e300);
        } catch (const std::runtime_error&) {
            threw_runtime = true;
        } catch (...) {
        }
        check(threw_runtime,
              "regularized_incomplete_beta with pathological extreme (a, b) throws "
              "std::runtime_error (non-convergence diagnostic) rather than returning a "
              "silently wrong value");
    }

    // ================================================================
    // Scaled beta: E (a=5, b=2, min=2, max=5), SciPy beta.ppf(u, a, b, loc, scale).
    // ================================================================
    {
        struct Case { double u; double expected; };
        const Case e_cases[] = {
            {0.0, 2.0},
            {1e-06, 2.133280855909455},
            {0.001, 2.544158399242598},
            {0.01, 2.8829410150408776},
            {0.1, 3.4690510803455252},
            {0.3, 3.918926928859939},
            {0.5, 4.20665005011302},
            {0.7, 4.454589586043152},
            {0.9, 4.722214223260614},
            {0.99, 4.919710426571735},
            {0.999, 4.97523352163676},
            {0.999999, 4.999225136511926},
            {1.0, 5.0},
        };
        for (const auto& c : e_cases) {
            double actual = beta_ppf(c.u, 5.0, 2.0, 2.0, 5.0);
            check(close(actual, c.expected, kBetaAbs, kBetaRel),
                  "beta_ppf(" + std::to_string(c.u) + ", a=5, b=2, 2, 5) matches SciPy");
        }

        struct VfCase { double u; double expected; };
        const VfCase vf_cases[] = {
            {0.0, 0.0001},
            {1e-06, 0.0005398268245012011},
            {0.001, 0.0018957227175005738},
            {0.01, 0.003013705349634896},
            {0.1, 0.0049478685651402335},
            {0.5, 0.007381945165372967},
            {0.9, 0.009083306936760027},
            {0.99, 0.009735044407686725},
            {0.999, 0.009918270621401307},
            {1.0, 0.01},
        };
        for (const auto& c : vf_cases) {
            double actual = beta_ppf(c.u, 5.0, 2.0, 1e-4, 1e-2);
            check(close(actual, c.expected, kBetaAbs, kBetaRel),
                  "beta_ppf(" + std::to_string(c.u) + ", a=5, b=2, 1e-4, 1e-2) matches SciPy");
        }

        // Dense tail-stress case: a=0.5, b=0.5 (U-shaped, harder for the
        // continued fraction near both endpoints), min=0, max=1.
        struct SymCase { double u; double expected; };
        const SymCase sym_cases[] = {
            {1e-06, 2.46740110027031e-12},
            {0.001, 2.4673990709169446e-06},
            {0.01, 0.00024671981713422146},
            {0.5, 0.4999999999999999},
            {0.99, 0.9997532801828658},
            {0.999, 0.999997532600929},
            {0.999999, 0.9999999999975326},
        };
        for (const auto& c : sym_cases) {
            double actual = beta_ppf(c.u, 0.5, 0.5, 0.0, 1.0);
            check(close(actual, c.expected, kBetaAbs, kBetaRel),
                  "beta_ppf(" + std::to_string(c.u) + ", a=0.5, b=0.5, 0, 1) matches SciPy");
        }
        // a=b symmetry about the midpoint for the symmetric-shape case.
        check(close(beta_ppf(0.3, 0.5, 0.5, 0.0, 1.0), 1.0 - beta_ppf(0.7, 0.5, 0.5, 0.0, 1.0),
                    kBetaAbs, kBetaRel),
              "beta_ppf(u, a, a, 0, 1) == 1 - beta_ppf(1-u, a, a, 0, 1) (a==b symmetry)");

        // Endpoints exact.
        check(beta_ppf(0.0, 5.0, 2.0, 2.0, 5.0) == 2.0, "beta_ppf(0.0, ...) == min exactly");
        check(beta_ppf(1.0, 5.0, 2.0, 2.0, 5.0) == 5.0, "beta_ppf(1.0, ...) == max exactly");

        // Monotonicity, dense grid.
        double prev = beta_ppf(0.0, 5.0, 2.0, 2.0, 5.0);
        for (int i = 1; i <= 1000; ++i) {
            double u = static_cast<double>(i) / 1000.0;
            double cur = beta_ppf(u, 5.0, 2.0, 2.0, 5.0);
            check(cur >= prev, "beta_ppf strictly increasing (E, a=5,b=2), i=" + std::to_string(i));
            prev = cur;
        }

        check(throws_invalid_argument([] { beta_ppf(0.5, 0.0, 2.0, 2.0, 5.0); }),
              "beta_ppf with a<=0 throws std::invalid_argument");
        check(throws_invalid_argument([] { beta_ppf(0.5, 5.0, -1.0, 2.0, 5.0); }),
              "beta_ppf with negative b throws std::invalid_argument");
        check(throws_invalid_argument([] { beta_ppf(0.5, 5.0, 2.0, 5.0, 2.0); }),
              "beta_ppf with min > max throws std::invalid_argument");
        check(throws_invalid_argument([] { beta_ppf(-0.1, 5.0, 2.0, 2.0, 5.0); }),
              "beta_ppf(u<0, ...) throws std::invalid_argument");
        check(throws_invalid_argument([] { beta_ppf(1.1, 5.0, 2.0, 2.0, 5.0); }),
              "beta_ppf(u>1, ...) throws std::invalid_argument");
    }

    // ================================================================
    // Binomial ppf: I0=400, p=0.01 and p=0.03 (Jones community rates),
    // SciPy binom.ppf(u, n, p).
    // ================================================================
    {
        struct Case { double u; long long expected; };
        const Case p01_cases[] = {
            {1e-06, 0},   {0.001, 0}, {0.01, 0},  {0.05, 1},    {0.1, 2},
            {0.25, 3},    {0.5, 4},   {0.75, 5},  {0.9, 7},     {0.95, 8},
            {0.99, 9},    {0.999, 11}, {0.9999, 13}, {0.99999999, 19},
        };
        for (const auto& c : p01_cases) {
            long long actual = binomial_ppf(c.u, 400, 0.01);
            check(actual == c.expected, "binomial_ppf(" + std::to_string(c.u) +
                                             ", n=400, p=0.01) == " + std::to_string(c.expected) +
                                             " (SciPy exact match)");
        }

        const Case p03_cases[] = {
            {1e-06, 0},   {0.001, 3}, {0.01, 5},  {0.05, 7},    {0.1, 8},
            {0.25, 10},   {0.5, 12},  {0.75, 14}, {0.9, 16},    {0.95, 18},
            {0.99, 21},   {0.999, 24}, {0.9999, 27}, {0.99999999, 35},
        };
        for (const auto& c : p03_cases) {
            long long actual = binomial_ppf(c.u, 400, 0.03);
            check(actual == c.expected, "binomial_ppf(" + std::to_string(c.u) +
                                             ", n=400, p=0.03) == " + std::to_string(c.expected) +
                                             " (SciPy exact match)");
        }
    }

    // ================================================================
    // Binomial ppf: boundary/tie stress at exact CDF steps, n=400, p=0.01.
    // SciPy binom.cdf(k, 400, 0.01) at k in {0,1,2,3,4,5,10,20}, then probed
    // just below / exactly at / just above each step.
    // ================================================================
    {
        struct StepCase { int k; double cdf; int expected_at_literal; };
        const StepCase steps[] = {
            {0, 0.017950553275045196, 0},  {1, 0.09047804125502519, 1},
            {2, 0.2366319185480154, 2},    {3, 0.43248795613256125, 3},
            {4, 0.6288385796706073, 4},    {5, 0.7859190785010428, 6},
            {10, 0.9973169239092451, 10},  {20, 0.9999999986484686, 20},
        };
        for (const auto& s : steps) {
            const double below = s.cdf - 1e-12;
            const double exact = s.cdf;
            const double above = s.cdf + 1e-12;
            check(binomial_ppf(below, 400, 0.01) == s.k,
                  "binomial_ppf(just below cdf(k=" + std::to_string(s.k) + ")) == k");
            // Pin the deterministic result for each literal. Only k=5 differs from the
            // mathematical boundary because this SciPy CDF literal lies about 4.4e-16 above
            // the mode-centered recurrence CDF. MODEL_SPECIFICATION.md section 11.5 permits
            // that documented one-step boundary discrepancy; all other rows remain exact.
            const long long exact_result = binomial_ppf(exact, 400, 0.01);
            check(exact_result == s.expected_at_literal,
                  "binomial_ppf(exactly cdf(k=" + std::to_string(s.k) +
                      ")) matches its pinned recurrence/SciPy-literal boundary result");
            check(binomial_ppf(above, 400, 0.01) == s.k + 1,
                  "binomial_ppf(just above cdf(k=" + std::to_string(s.k) + ")) == k+1");
        }
    }

    // Four independently observed SciPy/C++ one-step ties from the Stage 5
    // dense audit. Each literal is within 6.6e-16 of its CDF boundary, well
    // inside MODEL_SPECIFICATION.md section 11.5's documented 1e-12 band.
    {
        struct ToleratedTie { int n; double p; double u; long long cpp_result; };
        const ToleratedTie ties[] = {
            {1, 0.9, 0.10000000000000001, 1},
            {3, 0.9, 0.001, 1},
            {161, 0.5, 0.5, 80},
            {325, 0.5, 0.5, 162},
        };
        for (const auto& tie : ties) {
            check(binomial_ppf(tie.u, tie.n, tie.p) == tie.cpp_result,
                  "binomial_ppf documented tolerated CDF-boundary tie is pinned");
        }
    }

    // ================================================================
    // Binomial ppf: dense CDF fixture, n=10, p=0.3, exercising every k.
    // ================================================================
    {
        struct DenseCase { int k; double cdf; };
        const DenseCase cases[] = {
            {0, 0.02824752489999999}, {1, 0.1493083458999999}, {2, 0.38278278639999974},
            {3, 0.6496107184000001},  {4, 0.8497316674000001}, {5, 0.9526510126},
            {6, 0.9894079216},        {7, 0.9984096136},       {8, 0.9998563141},
            {9, 0.9999940951},        {10, 1.0},
        };
        int prev_k = -1;
        for (const auto& c : cases) {
            long long actual = binomial_ppf(c.cdf - 1e-13, 10, 0.3);
            check(actual == c.k,
                  "binomial_ppf(just below dense cdf(k=" + std::to_string(c.k) + ")) == k");
            check(actual > prev_k || c.k == 0,
                  "binomial_ppf dense n=10 p=0.3 result is non-decreasing across k");
            prev_k = static_cast<int>(actual);
        }
    }

    // ================================================================
    // Binomial ppf: boundaries n=0, p=0, p=1.
    // ================================================================
    {
        check(binomial_ppf(0.5, 0, 0.3) == 0, "binomial_ppf(u, n=0, p) == 0 for any valid u");
        check(binomial_ppf(1e-9, 0, 0.3) == 0, "binomial_ppf(tiny u, n=0, p) == 0");
        check(binomial_ppf(1.0, 0, 0.3) == 0, "binomial_ppf(1.0, n=0, p) == 0");

        check(binomial_ppf(0.999, 10, 0.0) == 0, "binomial_ppf(u, n=10, p=0.0) == 0 for any valid u");
        check(binomial_ppf(1.0, 10, 0.0) == 10,
              "binomial_ppf(1.0, n=10, p=0.0) uses the upper-support endpoint convention");

        check(binomial_ppf(0.001, 10, 1.0) == 10,
              "binomial_ppf(u, n=10, p=1.0) == n for any valid u");
        check(binomial_ppf(1.0, 10, 1.0) == 10, "binomial_ppf(1.0, n=10, p=1.0) == n");
        check(binomial_ppf(1.0, 10, 0.01) == 10,
              "binomial_ppf(1.0, n=10, p=0.01) == n through the general path");
        check(binomial_ppf(1.0, 100, 0.5) == 100,
              "binomial_ppf(1.0, n=100, p=0.5) == n through the general path");
        check(binomial_ppf(1.0, 400, 0.9) == 400,
              "binomial_ppf(1.0, n=400, p=0.9) == n through the general path");
    }

    // ================================================================
    // Binomial ppf: invalid inputs.
    // ================================================================
    {
        check(throws_invalid_argument([] { binomial_ppf(0.0, 400, 0.01); }),
              "binomial_ppf(0.0, ...) throws std::invalid_argument");
        check(throws_invalid_argument([] { binomial_ppf(-0.1, 400, 0.01); }),
              "binomial_ppf(u<0, ...) throws std::invalid_argument");
        check(throws_invalid_argument([] { binomial_ppf(1.1, 400, 0.01); }),
              "binomial_ppf(u>1, ...) throws std::invalid_argument");
        check(throws_invalid_argument(
                  [] { binomial_ppf(std::numeric_limits<double>::quiet_NaN(), 400, 0.01); }),
              "binomial_ppf(NaN, ...) throws std::invalid_argument");
        check(throws_invalid_argument([] { binomial_ppf(0.5, -1, 0.01); }),
              "binomial_ppf(u, n=-1, p) throws std::invalid_argument");
        check(throws_invalid_argument([] { binomial_ppf(0.5, 400, -0.01); }),
              "binomial_ppf(u, n, p<0) throws std::invalid_argument");
        check(throws_invalid_argument([] { binomial_ppf(0.5, 400, 1.01); }),
              "binomial_ppf(u, n, p>1) throws std::invalid_argument");
        check(throws_invalid_argument(
                  [] { binomial_ppf(0.5, 400, std::numeric_limits<double>::quiet_NaN()); }),
              "binomial_ppf(u, n, p=NaN) throws std::invalid_argument");
    }

    // ================================================================
    // C2 (Opus finding): normal/lognormal/log10normal at u = 1 - 1e-13 and
    // u = 1e-13, dense-tail SciPy literals (scipy.stats.norm/lognorm .ppf).
    // The Halley refinement previously suffered catastrophic cancellation
    // subtracting u (~1) from Phi(x) (~1) in the upper tail, missing the
    // abs-or-rel 1e-9 tolerance (MODEL_SPECIFICATION.md section 11.5).
    // ================================================================
    {
        const double u_hi = 1.0 - 1e-13;
        const double u_lo = 1e-13;
        check(close(normal_ppf(u_hi), 7.3487545403000425, kNormalAbs, kNormalRel),
              "normal_ppf(1-1e-13) matches SciPy (upper-tail cancellation guard)");
        check(close(normal_ppf(u_lo), -7.348796102800677, kNormalAbs, kNormalRel),
              "normal_ppf(1e-13) matches SciPy");
        check(close(lognormal_gm_gsd_ppf(u_hi, 0.52, 1.9), 58.14269180949784, kNormalAbs,
                    kNormalRel),
              "lognormal_gm_gsd_ppf(1-1e-13, GM=0.52, GSD=1.9) matches SciPy");
        check(close(lognormal_gm_gsd_ppf(u_lo, 0.52, 1.9), 0.004650503411626305, kNormalAbs,
                    kNormalRel),
              "lognormal_gm_gsd_ppf(1e-13, GM=0.52, GSD=1.9) matches SciPy");
        check(close(log10normal_ppf(u_hi, 7.0, 1.4), 1.9420318867064157e+17, kNormalAbs,
                    kNormalRel),
              "log10normal_ppf(1-1e-13, mu=7, sigma=1.4) matches SciPy");
        check(close(log10normal_ppf(u_lo, 7.0, 1.4), 0.0005148556180499769, kNormalAbs,
                    kNormalRel),
              "log10normal_ppf(1e-13, mu=7, sigma=1.4) matches SciPy");
    }

    // ================================================================
    // I1 (Opus finding): normal_ppf must return a finite value (or reject
    // consistently) at the smallest positive double (denorm_min), where the
    // Halley refinement's exp(x*x/2) intermediate overflows to +inf, and
    // remain accurate at the production LHS engine's actual minimum draw,
    // 2^-54 (docs/MODEL_SPECIFICATION.md / lhs.h never emit 0.0 or 1.0, and
    // 2^-54 is the smallest nonzero draw the exact-size 1D LHS engine can
    // produce).
    // ================================================================
    {
        const double denorm_min = std::numeric_limits<double>::denorm_min();
        double v = std::numeric_limits<double>::quiet_NaN();
        bool threw = false;
        try {
            v = normal_ppf(denorm_min);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(!threw && std::isfinite(v),
              "normal_ppf(denorm_min) returns a finite rational-approximation fallback");
        check(close(v, -38.467405617144344, kNormalAbs, kNormalRel),
              "normal_ppf(denorm_min) meets the frozen SciPy tolerance");
        check(close(normal_ppf(denorm_min, 10.0, 2.5), -86.16851404286086,
                    kNormalAbs, kNormalRel),
              "location-scale normal_ppf(denorm_min) meets the frozen SciPy tolerance");

        const double two_minus_54 = std::ldexp(1.0, -54);
        check(close(normal_ppf(two_minus_54), -8.292361075813597, kNormalAbs, kNormalRel),
              "normal_ppf(2^-54) matches SciPy (production LHS minimum draw)");
    }

    // ================================================================
    // C1 (Opus finding): binomial_ppf at n=400, p=0.9 and p=0.99 -- both
    // combinations where mean = n*p sits far from 0, so pmf(0) = (1-p)^n
    // underflows to exactly 0.0 (0.1^400 and 0.01^400 both underflow double
    // precision), and the old all-n forward-from-zero recurrence multiplies
    // zero by finite ratios forever, always falling through to the "return
    // n" fallback regardless of u. Dense grid including u=.001, .5, .999 and
    // tight neighborhoods, SciPy scipy.stats.binom.ppf(u, n, p) literals.
    // ================================================================
    {
        struct Case { double u; long long expected; };
        const Case p90_cases[] = {
            {1e-06, 329},     {1e-09, 320},     {0.00099, 340},   {0.001, 340},
            {0.00101, 340},   {0.0005, 339},    {0.0015, 341},    {0.01, 345},
            {0.05, 350},      {0.1, 352},        {0.25, 356},      {0.4999, 360},
            {0.5, 360},       {0.5001, 360},     {0.499, 360},     {0.501, 360},
            {0.75, 364},      {0.9, 368},        {0.95, 370},      {0.99, 373},
            {0.99899, 377},   {0.999, 377},      {0.99901, 377},   {0.9985, 377},
            {0.9995, 378},    {0.9999, 380},     {0.99999999, 389},
        };
        for (const auto& c : p90_cases) {
            long long actual = binomial_ppf(c.u, 400, 0.9);
            check(actual == c.expected,
                  "binomial_ppf(" + std::to_string(c.u) + ", n=400, p=0.9) == " +
                      std::to_string(c.expected) + " (SciPy exact match, mode-centered recurrence)");
        }

        const Case p99_cases[] = {
            {1e-06, 384},     {1e-09, 379},     {0.00099, 389},   {0.001, 389},
            {0.00101, 389},   {0.0005, 388},    {0.0015, 389},    {0.01, 391},
            {0.05, 392},      {0.1, 393},        {0.25, 395},      {0.4999, 396},
            {0.5, 396},       {0.5001, 396},     {0.499, 396},     {0.501, 396},
            {0.75, 397},      {0.9, 398},        {0.95, 399},      {0.99, 400},
            {0.99899, 400},   {0.999, 400},      {0.99901, 400},   {0.9985, 400},
            {0.9995, 400},    {0.9999, 400},     {0.99999999, 400},
        };
        for (const auto& c : p99_cases) {
            long long actual = binomial_ppf(c.u, 400, 0.99);
            check(actual == c.expected,
                  "binomial_ppf(" + std::to_string(c.u) + ", n=400, p=0.99) == " +
                      std::to_string(c.expected) + " (SciPy exact match, mode-centered recurrence)");
        }
    }

    // ================================================================
    // BinomialPpfTable: for a fixed (n, p), ppf(u) must equal the
    // single-shot binomial_ppf(u, n, p) exactly for every valid u --
    // dense grid, extreme tails, the pinned CDF-step / tie literals, and
    // the degenerate n=0 / p=0 / p=1 / u=1 conventions -- and its
    // validation and max-size guards must mirror binomial_ppf's.
    // ================================================================
    {
        struct NP { long long n; double p; };
        const NP params[] = {
            {400, 0.01}, {400, 0.03}, {400, 0.9},  {400, 0.99}, {10, 0.3},
            {1, 0.5},    {1, 0.9},    {3, 0.9},    {161, 0.5},  {325, 0.5},
            {2000, 0.5}, {50, 1e-6},  {50, 1.0 - 1e-6},
        };
        const double tails[] = {
            std::numeric_limits<double>::denorm_min(), 1e-300, 1e-15, 1e-12, 1e-9,
            std::ldexp(1.0, -54), 1.0 - 1e-9, 1.0 - 1e-12, 1.0 - std::ldexp(1.0, -53), 1.0,
        };
        const int grid = 20000;
        for (const NP& np : params) {
            const BinomialPpfTable table(np.n, np.p);
            bool all_equal = table.n() == np.n && table.p() == np.p;
            bool monotone = true;
            long long previous = -1;
            for (int i = 1; i <= grid && all_equal; ++i) {
                const double u = static_cast<double>(i) / static_cast<double>(grid);
                const long long actual = table.ppf(u);
                all_equal = actual == binomial_ppf(u, np.n, np.p);
                monotone = monotone && actual >= previous;
                previous = actual;
            }
            for (double u : tails) {
                all_equal = all_equal && table.ppf(u) == binomial_ppf(u, np.n, np.p);
            }
            check(all_equal, "BinomialPpfTable(n=" + std::to_string(np.n) + ", p=" +
                                 std::to_string(np.p) + ").ppf equals binomial_ppf over " +
                                 std::to_string(grid) + " dense u values and the tails");
            check(monotone, "BinomialPpfTable(n=" + std::to_string(np.n) + ", p=" +
                                std::to_string(np.p) + ").ppf is non-decreasing in u");
        }

        // Pinned boundary literals (exact CDF steps, tolerated ties, dense
        // n=10 p=0.3 steps): the table must resolve each one identically.
        {
            const BinomialPpfTable t400(400, 0.01);
            const double steps[] = {0.017950553275045196, 0.09047804125502519,
                                    0.2366319185480154,   0.43248795613256125,
                                    0.6288385796706073,   0.7859190785010428,
                                    0.9973169239092451,   0.9999999986484686};
            bool boundary_equal = true;
            for (double cdf : steps) {
                for (double u : {cdf - 1e-12, cdf, cdf + 1e-12}) {
                    boundary_equal = boundary_equal && t400.ppf(u) == binomial_ppf(u, 400, 0.01);
                }
            }
            check(boundary_equal,
                  "BinomialPpfTable(400, 0.01) matches binomial_ppf just below / exactly at / "
                  "just above every pinned CDF step");

            struct Tie { long long n; double p; double u; long long expected; };
            const Tie ties[] = {{1, 0.9, 0.10000000000000001, 1},
                                {3, 0.9, 0.001, 1},
                                {161, 0.5, 0.5, 80},
                                {325, 0.5, 0.5, 162}};
            for (const Tie& tie : ties) {
                check(BinomialPpfTable(tie.n, tie.p).ppf(tie.u) == tie.expected,
                      "BinomialPpfTable reproduces the pinned tolerated CDF-boundary tie");
            }

            const BinomialPpfTable t10(10, 0.3);
            const double dense[] = {0.02824752489999999, 0.1493083458999999, 0.38278278639999974,
                                    0.6496107184000001,  0.8497316674000001, 0.9526510126,
                                    0.9894079216,        0.9984096136,       0.9998563141,
                                    0.9999940951,        1.0};
            bool dense_equal = true;
            for (double cdf : dense) {
                const double u = cdf - 1e-13;
                dense_equal = dense_equal && t10.ppf(u) == binomial_ppf(u, 10, 0.3);
            }
            check(dense_equal,
                  "BinomialPpfTable(10, 0.3) matches binomial_ppf just below every dense CDF step");
        }

        // Degenerate conventions, answered without a table buffer.
        {
            const BinomialPpfTable n0(0, 0.3);
            check(n0.ppf(0.5) == 0 && n0.ppf(1e-9) == 0 && n0.ppf(1.0) == 0,
                  "BinomialPpfTable(n=0).ppf == 0 for every valid u");
            const BinomialPpfTable p0(10, 0.0);
            check(p0.ppf(0.999) == 0 && p0.ppf(1.0) == 10,
                  "BinomialPpfTable(p=0).ppf == 0 below u=1 and n at u=1 (SciPy endpoint "
                  "convention)");
            const BinomialPpfTable p1(10, 1.0);
            check(p1.ppf(0.001) == 10 && p1.ppf(1.0) == 10,
                  "BinomialPpfTable(p=1).ppf == n for every valid u");
            check(BinomialPpfTable(10, 0.01).ppf(1.0) == 10 &&
                      BinomialPpfTable(100, 0.5).ppf(1.0) == 100 &&
                      BinomialPpfTable(400, 0.9).ppf(1.0) == 400,
                  "BinomialPpfTable.ppf(1.0) == n through the general path");
        }

        // Validation mirrors binomial_ppf's contract.
        {
            check(throws_invalid_argument([] { BinomialPpfTable(-1, 0.01); }),
                  "BinomialPpfTable(n=-1, p) throws std::invalid_argument");
            check(throws_invalid_argument([] { BinomialPpfTable(400, -0.01); }),
                  "BinomialPpfTable(n, p<0) throws std::invalid_argument");
            check(throws_invalid_argument([] { BinomialPpfTable(400, 1.01); }),
                  "BinomialPpfTable(n, p>1) throws std::invalid_argument");
            check(throws_invalid_argument(
                      [] { BinomialPpfTable(400, std::numeric_limits<double>::quiet_NaN()); }),
                  "BinomialPpfTable(n, p=NaN) throws std::invalid_argument");
            const BinomialPpfTable table(400, 0.01);
            check(throws_invalid_argument([&] { table.ppf(0.0); }),
                  "BinomialPpfTable::ppf(0.0) throws std::invalid_argument");
            check(throws_invalid_argument([&] { table.ppf(-0.1); }),
                  "BinomialPpfTable::ppf(u<0) throws std::invalid_argument");
            check(throws_invalid_argument([&] { table.ppf(1.1); }),
                  "BinomialPpfTable::ppf(u>1) throws std::invalid_argument");
            check(throws_invalid_argument(
                      [&] { table.ppf(std::numeric_limits<double>::quiet_NaN()); }),
                  "BinomialPpfTable::ppf(NaN) throws std::invalid_argument");
        }

        // Max-size guard: an n whose (n+1)-double buffer cannot exist is
        // rejected with std::runtime_error before any allocation, by both
        // the single-shot function and the table constructor.
        {
            const long long huge_n = std::numeric_limits<long long>::max();
            bool single_threw = false;
            try {
                binomial_ppf(0.5, huge_n, 0.5);
            } catch (const std::runtime_error&) {
                single_threw = true;
            }
            check(single_threw, "binomial_ppf(u, n=LLONG_MAX, p) throws std::runtime_error");
            bool table_threw = false;
            try {
                BinomialPpfTable(huge_n, 0.5);
            } catch (const std::runtime_error&) {
                table_threw = true;
            }
            check(table_threw, "BinomialPpfTable(n=LLONG_MAX, p) throws std::runtime_error");
        }
    }

    // ================================================================
    // Finite-output requirement across the Jones-relevant open-tail range,
    // for every transform, including very close to (0,1).
    // ================================================================
    {
        const double tails[] = {1e-15, 1e-12, 1e-9, 1.0 - 1e-9, 1.0 - 1e-12, 1.0 - 1e-15};
        for (double u : tails) {
            check(std::isfinite(normal_ppf(u)), "normal_ppf finite at tail u=" + std::to_string(u));
            check(std::isfinite(lognormal_gm_gsd_ppf(u, 0.52, 1.9)),
                  "lognormal_gm_gsd_ppf finite at tail u=" + std::to_string(u));
            check(std::isfinite(log10normal_ppf(u, 7.0, 1.4)),
                  "log10normal_ppf finite at tail u=" + std::to_string(u));
            check(std::isfinite(uniform_ppf(u, 0.42, 0.61)),
                  "uniform_ppf finite at tail u=" + std::to_string(u));
            check(std::isfinite(beta_ppf(u, 5.0, 2.0, 2.0, 5.0)),
                  "beta_ppf finite at tail u=" + std::to_string(u));
        }
    }

    if (g_failures > 0) {
        std::cerr << g_failures << " test(s) failed.\n";
        return 1;
    }
    std::cout << "All tests passed.\n";
    return 0;
}
