// analysis/ecai_uncertainty_study.cpp
//
// Narrowly-scoped statistical study (NOT a unit test, NOT wired into CTest):
// For a given occupancy category (default "Patient"), estimate the 96th
// percentile of ECAi and its uncertainty using one of two selectable
// methods (--method replicate|bootstrap, default replicate):
//
// replicate mode (original behavior): run R independent replicates (default
// 25). Each replicate constructs its own RandomNumberManager(N) and performs
// exactly N simulations (default 10,000) using the same model call sequence
// as ecai.cpp's main loop (sample_parameters followed by compute_ECAi,
// target infection probability 0.001, require_infectors=false,
// category-default community rate, default Jones QER distribution),
// collects the N exact (non-rounded) ECAi values, sorts them, and computes
// the 96th percentile via the same linear interpolation ecai.cpp uses.
// Across the R replicate ECAi values this program reports sample mean,
// sample standard deviation, standard error, and a two-sided 95% Student-t
// confidence interval (df = R - 1).
//
// bootstrap mode: generate exactly N exact ECAi values once (same
// sample_parameters -> compute_ECAi path as replicate mode), compute their
// interpolated 96th percentile (the "original" estimate), then draw B
// bootstrap resamples (--resamples B), each of N values sampled with
// replacement from the original N values, computing the same interpolated
// 96th percentile for each resample. The percentile bootstrap 95% CI is the
// 2.5th and 97.5th percentiles (same interpolation) of the B bootstrap
// estimates. Bootstrap resampling uses its own std::mt19937_64 RNG seeded
// from std::random_device, independent of the RandomNumberManager used for
// the original N simulations.
//
// order-statistic mode: generate exactly N exact ECAi values once (same
// sample_parameters -> compute_ECAi path as the other modes), sort them,
// and report the interpolated 96th percentile as the point estimate (same
// interpolation as the other modes). A distribution-free two-sided 95%
// confidence interval for the population 96th percentile is then built
// directly from the sorted sample using the classical order-statistic
// (Hahn/Meeker "qbinom") construction: let K ~ Binomial(N, 0.96) be the
// number of the N sample values that fall at or below the true population
// 96th percentile. Let r = smallest k in [0,N] with P(K<=k) >= 0.025, and
// s = (smallest k in [0,N] with P(K<=k) >= 0.975) + 1; r and s are the
// 1-based ranks of the sorted sample (X_(1) <= ... <= X_(N), 1-based, i.e.
// X_(1) is the smallest of the N values) such that P(X_(r) <= xi_0.96 <=
// X_(s)) is approximately 0.95. If r < 1 or s > N, no valid rank bounds
// exist for the requested N (N is too small for a distribution-free 95% CI
// at the 96th percentile) and the program reports an error instead of a
// CI. Binomial tail probabilities are computed exactly via the regularized
// incomplete beta function (P(K<=k) = I_{1-p}(N-k, k+1)), not simulated.
//
// N is configurable via --simulations N (exact simulation count in all
// three modes); R (replicate mode only) via --repeats R; B (bootstrap mode
// only) via --resamples B. --repeats and --resamples are not valid outside
// their respective modes. With no arguments the defaults (method=replicate,
// N=10000, R=25) reproduce the original fixed-constant Patient-only
// behavior exactly. The occupancy category is configurable via --category
// NAME, accepting any category present in model.h's
// occupancy_params/category_order; the default remains "Patient".
//
// This file only reads from the existing ashrae241 library (model.h /
// random_manager.h); it does not modify simulation logic, defaults, the
// RandomNumberManager, or any of the four production executables.

#include "model.h"
#include "random_manager.h"

#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <limits>
#include <filesystem>
#include <iomanip>
#include <cstring>
#include <random>

namespace {

constexpr const char* kDefaultCategory = "Patient";
constexpr const char* kDefaultMethod   = "replicate";
constexpr int    kDefaultN          = 10000;
constexpr int    kDefaultReplicates = 25;
constexpr int    kDefaultResamples  = 2000;
constexpr double kTargetProb   = 0.001;
constexpr double kPercentile   = 0.96;
constexpr double kConfidenceLevel = 0.95;
constexpr double kAlpha = 1.0 - kConfidenceLevel;
constexpr bool   kRequireInfectors = false; // allows zero-infector draws (existing default)

// ---------------------------------------------------------------------------
// Student-t two-sided 97.5th-percentile critical value for arbitrary degrees
// of freedom, via the regularized incomplete beta function (continued
// fraction, Numerical Recipes formulation) and bisection on the t CDF.
// ---------------------------------------------------------------------------

double betacf(double a, double b, double x) {
    const int kMaxIter = 200;
    const double kEps = 3.0e-16;
    const double kTiny = 1.0e-300;

    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::fabs(d) < kTiny) d = kTiny;
    d = 1.0 / d;
    double h = d;

    for (int m = 1; m <= kMaxIter; m++) {
        int m2 = 2 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < kTiny) d = kTiny;
        c = 1.0 + aa / c;
        if (std::fabs(c) < kTiny) c = kTiny;
        d = 1.0 / d;
        h *= d * c;

        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < kTiny) d = kTiny;
        c = 1.0 + aa / c;
        if (std::fabs(c) < kTiny) c = kTiny;
        d = 1.0 / d;
        double del = d * c;
        h *= del;

        if (std::fabs(del - 1.0) < kEps) break;
    }
    return h;
}

double regularized_incomplete_beta(double a, double b, double x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    double ln_bt = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
                   a * std::log(x) + b * std::log(1.0 - x);
    double bt = std::exp(ln_bt);
    if (x < (a + 1.0) / (a + b + 2.0)) {
        return bt * betacf(a, b, x) / a;
    }
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

// ---------------------------------------------------------------------------
// Binomial CDF/quantile, used to build the order-statistic mode's rank
// bounds. Computed exactly via the regularized incomplete beta function
// (P(K<=k) = I_{1-p}(n-k, k+1) for a Binomial(n,p) random variable K),
// avoiding any risk of overflow/underflow from direct factorial terms even
// for n in the tens of thousands.
// ---------------------------------------------------------------------------

// P(K <= k) for K ~ Binomial(n, p), 0 <= k <= n.
double binomial_cdf(int k, int n, double p) {
    if (k < 0) return 0.0;
    if (k >= n) return 1.0;
    double a = static_cast<double>(n - k);
    double b = static_cast<double>(k + 1);
    return regularized_incomplete_beta(a, b, 1.0 - p);
}

// Smallest integer k in [0, n] such that P(K <= k) >= q, for K ~
// Binomial(n, p). Standard ("qbinom") quantile-function convention.
int binomial_quantile(double q, int n, double p) {
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (binomial_cdf(mid, n, p) >= q) hi = mid; else lo = mid + 1;
    }
    return lo;
}

// CDF of the Student-t distribution with df degrees of freedom.
double student_t_cdf(double t, double df) {
    double x = df / (df + t * t);
    double ib = regularized_incomplete_beta(df / 2.0, 0.5, x);
    return (t >= 0.0) ? (1.0 - 0.5 * ib) : (0.5 * ib);
}

// Two-sided 95% critical value (0.975 quantile) for the Student-t
// distribution with the given degrees of freedom, found by bisection.
double student_t_critical_975(int df) {
    double lo = 0.0, hi = 1.0e6;
    for (int i = 0; i < 200; i++) {
        double mid = 0.5 * (lo + hi);
        double p = student_t_cdf(mid, static_cast<double>(df));
        if (p < 0.975) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

// Interpolated percentile (same linear-interpolation convention ecai.cpp
// uses) of an already-sorted vector, at fraction p in [0, 1].
double interpolated_percentile(const std::vector<double>& sorted_values, double p) {
    const int n = static_cast<int>(sorted_values.size());
    double rank = p * (n - 1);
    int lo = static_cast<int>(std::floor(rank));
    int hi = static_cast<int>(std::ceil(rank));
    double frac = rank - lo;
    return sorted_values[lo] * (1 - frac) + sorted_values[hi] * frac;
}

// Generates n exact ECAi values for the given occupancy category via
// sample_parameters + compute_ECAi (same call sequence and argument
// semantics as ecai.cpp's main loop), using a fresh RandomNumberManager(n).
std::vector<double> generate_ecai_values(int n, const std::string& category) {
    RandomNumberManager rng(static_cast<std::size_t>(n));
    std::vector<double> ecai_list;
    ecai_list.reserve(n);
    for (int i = 0; i < n; i++) {
        auto par = sample_parameters(rng, category);
        auto [ecai_val, _infected_flag] = compute_ECAi(
            par, kTargetProb, rng, category,
            kRequireInfectors, /*override_community_rate=*/-1,
            QERDistribution{});
        ecai_list.push_back(ecai_val);
    }
    return ecai_list;
}

// Runs one replicate: n simulations of the given occupancy category.
// Returns the interpolated 96th percentile of the n exact ECAi values.
double run_replicate(int n, const std::string& category) {
    std::vector<double> ecai_list = generate_ecai_values(n, category);
    std::sort(ecai_list.begin(), ecai_list.end());
    return interpolated_percentile(ecai_list, kPercentile);
}

void print_usage(const char* prog) {
    std::fprintf(stderr,
        "Usage: %s [--method replicate|bootstrap|order-statistic] [--simulations N]\n"
        "           [--repeats R] [--resamples B] [--category NAME]\n"
        "  --method replicate|bootstrap|order-statistic\n"
        "                    Estimation method (default %s)\n"
        "  --simulations N   Simulations per replicate; original sample size in\n"
        "                    bootstrap mode; exact simulation count in\n"
        "                    order-statistic mode (positive integer, default %d)\n"
        "  --repeats R       Number of replicates; replicate mode only (positive\n"
        "                    integer >= 2, default %d)\n"
        "  --resamples B     Number of bootstrap resamples; bootstrap mode only\n"
        "                    (positive integer >= 2, default %d)\n"
        "  --category NAME   Occupancy category (default %s)\n",
        prog, kDefaultMethod, kDefaultN, kDefaultReplicates, kDefaultResamples, kDefaultCategory);
    std::fprintf(stderr, "Available categories:");
    for (const auto& k : category_order)
        std::fprintf(stderr, " %s", k.c_str());
    std::fprintf(stderr, "\n");
}

// Parses a positive integer from a string. Returns true on success.
bool parse_positive_int(const char* s, int* out) {
    if (s == nullptr || s[0] == '\0') return false;
    char* end = nullptr;
    long v = std::strtol(s, &end, 10);
    if (end == s || *end != '\0') return false;
    if (v <= 0 || v > std::numeric_limits<int>::max()) return false;
    *out = static_cast<int>(v);
    return true;
}

} // namespace

int main(int argc, char** argv) {
    int sim_count = kDefaultN;
    int replicate_count = kDefaultReplicates;
    int resample_count = kDefaultResamples;
    std::string category = kDefaultCategory;
    std::string method = kDefaultMethod;
    bool repeats_explicit = false;
    bool resamples_explicit = false;

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--category") == 0) {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "Error: --category requires an argument.\n");
                print_usage(argv[0]);
                return 1;
            }
            category = argv[++i];
        } else if (std::strcmp(argv[i], "--method") == 0) {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "Error: --method requires an argument.\n");
                print_usage(argv[0]);
                return 1;
            }
            method = argv[++i];
        } else if (std::strcmp(argv[i], "--simulations") == 0) {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "Error: --simulations requires an argument.\n");
                print_usage(argv[0]);
                return 1;
            }
            if (!parse_positive_int(argv[++i], &sim_count)) {
                std::fprintf(stderr, "Error: --simulations must be a positive integer, got '%s'.\n", argv[i]);
                print_usage(argv[0]);
                return 1;
            }
        } else if (std::strcmp(argv[i], "--repeats") == 0) {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "Error: --repeats requires an argument.\n");
                print_usage(argv[0]);
                return 1;
            }
            if (!parse_positive_int(argv[++i], &replicate_count)) {
                std::fprintf(stderr, "Error: --repeats must be a positive integer, got '%s'.\n", argv[i]);
                print_usage(argv[0]);
                return 1;
            }
            repeats_explicit = true;
        } else if (std::strcmp(argv[i], "--resamples") == 0) {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "Error: --resamples requires an argument.\n");
                print_usage(argv[0]);
                return 1;
            }
            if (!parse_positive_int(argv[++i], &resample_count)) {
                std::fprintf(stderr, "Error: --resamples must be a positive integer, got '%s'.\n", argv[i]);
                print_usage(argv[0]);
                return 1;
            }
            resamples_explicit = true;
        } else if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            std::fprintf(stderr, "Error: unrecognized argument '%s'.\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
    }

    if (method != "replicate" && method != "bootstrap" && method != "order-statistic") {
        std::fprintf(stderr, "Error: --method must be 'replicate', 'bootstrap', or 'order-statistic', got '%s'.\n", method.c_str());
        print_usage(argv[0]);
        return 1;
    }

    if (method != "replicate" && repeats_explicit) {
        std::fprintf(stderr, "Error: --repeats applies to replicate mode only; it is not valid with --method %s.\n", method.c_str());
        print_usage(argv[0]);
        return 1;
    }
    if (method != "bootstrap" && resamples_explicit) {
        std::fprintf(stderr, "Error: --resamples applies to bootstrap mode only; it is not valid with --method %s.\n", method.c_str());
        print_usage(argv[0]);
        return 1;
    }

    if (method == "replicate" && replicate_count < 2) {
        std::fprintf(stderr, "Error: --repeats must be at least 2 (need df >= 1 for the Student-t CI), got %d.\n",
                     replicate_count);
        return 1;
    }
    if (method == "bootstrap" && resample_count < 2) {
        std::fprintf(stderr, "Error: --resamples must be at least 2 (need >= 2 for a percentile CI), got %d.\n",
                     resample_count);
        return 1;
    }

    if (category.empty()) {
        std::fprintf(stderr, "Error: --category requires a non-empty value.\n");
        print_usage(argv[0]);
        return 1;
    }
    if (occupancy_params.find(category) == occupancy_params.end()) {
        std::fprintf(stderr, "Error: '%s' is not a valid category.\n", category.c_str());
        std::fprintf(stderr, "Available categories:");
        for (const auto& k : category_order)
            std::fprintf(stderr, " %s", k.c_str());
        std::fprintf(stderr, "\n");
        return 1;
    }

    const int kN = sim_count;
    const bool is_default_category = (category == kDefaultCategory);

    // Output directory and filename stems. Patient (the original default)
    // keeps the original replicate-mode filenames for backward
    // compatibility; any other category, or bootstrap mode, gets a
    // category/method-derived filename so results are never mistaken for
    // the original Patient replicate study.
    const std::string out_dir = "ecai_uncertainty_study_results";
    std::filesystem::create_directories(out_dir);

    if (method == "replicate") {
        const int kReplicates = replicate_count;
        const int df = kReplicates - 1;
        const double t_crit = student_t_critical_975(df);

        printf("\n%s ECAi Monte Carlo study (separate from production executables/tests)\n", category.c_str());
        printf("  Method:               replicate\n");
        printf("  Category:            %s\n", category.c_str());
        printf("  N per replicate:      %d\n", kN);
        printf("  Replicates:           %d\n", kReplicates);
        printf("  Target probability:   %.4f (%.2f%%)\n", kTargetProb, kTargetProb * 100);
        printf("  Percentile:           %.2f\n", kPercentile);
        printf("  Infector mode:        require_infectors=%s (zero-infector draws allowed)\n",
               kRequireInfectors ? "true" : "false");
        printf("  Community rate:       category default (%s)\n\n", category.c_str());

        std::vector<double> ecai_values;
        ecai_values.reserve(kReplicates);
        for (int r = 0; r < kReplicates; r++) {
            double ecai = run_replicate(kN, category);
            ecai_values.push_back(ecai);
            printf("  Replicate %2d/%d: ECAi_96th = %.6f L/s/person\n", r + 1, kReplicates, ecai);
        }

        // Sanity check: exactly kReplicates rows, all finite.
        bool all_finite = true;
        for (double v : ecai_values) {
            if (!std::isfinite(v)) { all_finite = false; break; }
        }
        bool row_count_ok = (static_cast<int>(ecai_values.size()) == kReplicates);

        double mean = std::numeric_limits<double>::quiet_NaN();
        double sd = std::numeric_limits<double>::quiet_NaN();
        double se = std::numeric_limits<double>::quiet_NaN();
        double ci_lo = std::numeric_limits<double>::quiet_NaN();
        double ci_hi = std::numeric_limits<double>::quiet_NaN();

        if (all_finite && row_count_ok) {
            double sum = 0;
            for (double v : ecai_values) sum += v;
            mean = sum / kReplicates;

            double ss = 0;
            for (double v : ecai_values) ss += (v - mean) * (v - mean);
            sd = std::sqrt(ss / (kReplicates - 1)); // sample std dev, df = n-1
            se = sd / std::sqrt(static_cast<double>(kReplicates));

            double margin = t_crit * se;
            ci_lo = mean - margin;
            ci_hi = mean + margin;
        }

        printf("\nSummary across %d replicates:\n", kReplicates);
        printf("  Sample mean:                    %.6f L/s/person\n", mean);
        printf("  Sample standard deviation:      %.6f L/s/person\n", sd);
        printf("  Standard error:                 %.6f L/s/person\n", se);
        printf("  95%% CI (two-sided t, df=%d):    [%.6f, %.6f]\n", df, ci_lo, ci_hi);
        printf("\nSanity check: %d/%d rows present, all finite = %s\n",
               static_cast<int>(ecai_values.size()), kReplicates,
               all_finite ? "true" : "false");

        const std::string replicates_stem = is_default_category
            ? "patient_ecai_96th_percentile_replicates"
            : category + "_ecai_96th_percentile_replicates";
        const std::string summary_stem = is_default_category
            ? "patient_ecai_96th_percentile_summary"
            : category + "_ecai_96th_percentile_summary";

        const std::string csv_path = out_dir + "/" + replicates_stem + ".csv";
        std::ofstream csv(csv_path);
        csv << std::fixed << std::setprecision(9);
        csv << "replicate,ECAi_96th_Lps\n";
        for (int r = 0; r < kReplicates; r++) {
            csv << (r + 1) << "," << ecai_values[r] << "\n";
        }
        csv.close();

        const std::string summary_path = out_dir + "/" + summary_stem + ".csv";
        std::ofstream summary(summary_path);
        summary << std::fixed << std::setprecision(9);
        summary << "parameter,value\n";
        summary << "method,replicate\n";
        summary << "category," << category << "\n";
        summary << "N," << kN << "\n";
        summary << "replicates," << kReplicates << "\n";
        summary << "target_probability," << kTargetProb << "\n";
        summary << "percentile," << kPercentile << "\n";
        summary << "require_infectors," << (kRequireInfectors ? "true" : "false") << "\n";
        summary << "sample_mean_Lps," << mean << "\n";
        summary << "sample_stddev_Lps," << sd << "\n";
        summary << "standard_error_Lps," << se << "\n";
        summary << "t_critical_df" << df << "_0975," << t_crit << "\n";
        summary << "ci95_lower_Lps," << ci_lo << "\n";
        summary << "ci95_upper_Lps," << ci_hi << "\n";
        summary << "row_count_ok," << (row_count_ok ? "true" : "false") << "\n";
        summary << "all_finite," << (all_finite ? "true" : "false") << "\n";
        summary.close();

        printf("\nResults written to:\n  %s\n  %s\n", csv_path.c_str(), summary_path.c_str());

        return 0;
    }

    if (method == "order-statistic") {
        // Rank bounds depend only on N and the target percentile, not on the
        // simulated data, so validate them before running any (potentially
        // expensive) simulations. 1-based ranks r and s into the sorted
        // sample X_(1) <= ... <= X_(N): X_(r) is the CI lower bound, X_(s)
        // is the CI upper bound. r/s are indices into the N sorted values,
        // one-based (X_(1) is the smallest of the N values, X_(N) the
        // largest); the code below converts to 0-based std::vector indices
        // (sorted[r - 1], sorted[s - 1]) when reading the bound values.
        const int kOrderN = sim_count;
        const int r = binomial_quantile(kAlpha / 2.0, kOrderN, kPercentile);
        const int r_hi = binomial_quantile(1.0 - kAlpha / 2.0, kOrderN, kPercentile);
        const int s = r_hi + 1;

        if (r < 1 || s > kOrderN) {
            std::fprintf(stderr,
                "Error: no valid order-statistic rank bounds for N=%d at the %.2f percentile with %.0f%% confidence.\n"
                "  Computed lower rank r=%d, upper rank s=%d (valid range is 1..%d).\n"
                "  N is too small for a distribution-free order-statistic CI at this percentile/confidence;\n"
                "  increase --simulations.\n",
                kOrderN, kPercentile, kConfidenceLevel * 100, r, s, kOrderN);
            return 1;
        }

        printf("\n%s ECAi Monte Carlo study (separate from production executables/tests)\n", category.c_str());
        printf("  Method:               order-statistic\n");
        printf("  Category:            %s\n", category.c_str());
        printf("  N (simulations):      %d\n", kOrderN);
        printf("  Target probability:   %.4f (%.2f%%)\n", kTargetProb, kTargetProb * 100);
        printf("  Percentile:           %.2f\n", kPercentile);
        printf("  Confidence level:     %.2f\n", kConfidenceLevel);
        printf("  Infector mode:        require_infectors=%s (zero-infector draws allowed)\n",
               kRequireInfectors ? "true" : "false");
        printf("  Community rate:       category default (%s)\n", category.c_str());
        printf("  Order-statistic ranks (1-based, into the N sorted values): r=%d (lower), s=%d (upper)\n\n",
               r, s);

        std::vector<double> values = generate_ecai_values(kOrderN, category);
        std::vector<double> sorted_values = values;
        std::sort(sorted_values.begin(), sorted_values.end());

        const double point_estimate = interpolated_percentile(sorted_values, kPercentile);
        const double ci_lo = sorted_values[static_cast<std::size_t>(r) - 1];
        const double ci_hi = sorted_values[static_cast<std::size_t>(s) - 1];

        bool all_finite = std::isfinite(point_estimate) && std::isfinite(ci_lo) && std::isfinite(ci_hi);
        bool row_count_ok = (static_cast<int>(sorted_values.size()) == kOrderN);

        printf("Point estimate (interpolated %.2f percentile): %.6f L/s/person\n", kPercentile, point_estimate);
        printf("95%% order-statistic CI: [%.6f, %.6f] L/s/person (X_(%d), X_(%d))\n",
               ci_lo, ci_hi, r, s);
        printf("\nSanity check: %d/%d values present, all finite = %s\n",
               static_cast<int>(sorted_values.size()), kOrderN, all_finite ? "true" : "false");

        const std::string values_stem = is_default_category
            ? "patient_ecai_96th_percentile_order_statistic_values"
            : category + "_ecai_96th_percentile_order_statistic_values";
        const std::string summary_stem = is_default_category
            ? "patient_ecai_96th_percentile_order_statistic_summary"
            : category + "_ecai_96th_percentile_order_statistic_summary";

        const std::string csv_path = out_dir + "/" + values_stem + ".csv";
        std::ofstream csv(csv_path);
        csv << std::fixed << std::setprecision(9);
        csv << "rank_1based,ECAi_Lps\n";
        for (int i = 0; i < kOrderN; i++) {
            csv << (i + 1) << "," << sorted_values[i] << "\n";
        }
        csv.close();

        const std::string summary_path = out_dir + "/" + summary_stem + ".csv";
        std::ofstream summary(summary_path);
        summary << std::fixed << std::setprecision(9);
        summary << "parameter,value\n";
        summary << "method,order-statistic\n";
        summary << "category," << category << "\n";
        summary << "N," << kOrderN << "\n";
        summary << "target_probability," << kTargetProb << "\n";
        summary << "percentile," << kPercentile << "\n";
        summary << "confidence_level," << kConfidenceLevel << "\n";
        summary << "require_infectors," << (kRequireInfectors ? "true" : "false") << "\n";
        summary << "point_estimate_Lps," << point_estimate << "\n";
        summary << "rank_lower_1based," << r << "\n";
        summary << "rank_upper_1based," << s << "\n";
        summary << "indexing_convention,1-based ranks into N ascending-sorted ECAi values (X_(1)..X_(N))\n";
        summary << "ci95_lower_Lps," << ci_lo << "\n";
        summary << "ci95_upper_Lps," << ci_hi << "\n";
        summary << "row_count_ok," << (row_count_ok ? "true" : "false") << "\n";
        summary << "all_finite," << (all_finite ? "true" : "false") << "\n";
        summary.close();

        printf("\nResults written to:\n  %s\n  %s\n", csv_path.c_str(), summary_path.c_str());

        return 0;
    }

    // --- Bootstrap mode -----------------------------------------------
    const int kResamples = resample_count;

    printf("\n%s ECAi Monte Carlo study (separate from production executables/tests)\n", category.c_str());
    printf("  Method:               bootstrap\n");
    printf("  Category:            %s\n", category.c_str());
    printf("  N (original sample): %d\n", kN);
    printf("  Bootstrap resamples:  %d\n", kResamples);
    printf("  Target probability:   %.4f (%.2f%%)\n", kTargetProb, kTargetProb * 100);
    printf("  Percentile:           %.2f\n", kPercentile);
    printf("  Infector mode:        require_infectors=%s (zero-infector draws allowed)\n",
           kRequireInfectors ? "true" : "false");
    printf("  Community rate:       category default (%s)\n\n", category.c_str());

    // Original N exact ECAi values, via the same sample_parameters ->
    // compute_ECAi path used by replicate mode.
    std::vector<double> original_values = generate_ecai_values(kN, category);
    std::vector<double> original_sorted = original_values;
    std::sort(original_sorted.begin(), original_sorted.end());
    double original_percentile = interpolated_percentile(original_sorted, kPercentile);

    printf("  Original ECAi_96th (N=%d): %.6f L/s/person\n\n", kN, original_percentile);

    // Independent bootstrap RNG, seeded consistently with the existing
    // random-device behavior used by RandomNumberManager, but kept
    // separate from it.
    std::mt19937_64 bootstrap_rng(std::random_device{}());
    std::uniform_int_distribution<std::size_t> index_dist(0, static_cast<std::size_t>(kN - 1));

    std::vector<double> bootstrap_estimates;
    bootstrap_estimates.reserve(kResamples);
    std::vector<double> resample_buf(kN);
    for (int b = 0; b < kResamples; b++) {
        for (int i = 0; i < kN; i++) {
            resample_buf[i] = original_values[index_dist(bootstrap_rng)];
        }
        std::sort(resample_buf.begin(), resample_buf.end());
        double est = interpolated_percentile(resample_buf, kPercentile);
        bootstrap_estimates.push_back(est);
        if ((b + 1) % std::max(1, kResamples / 10) == 0 || b + 1 == kResamples) {
            printf("  Bootstrap resample %4d/%d: ECAi_96th = %.6f L/s/person\n", b + 1, kResamples, est);
        }
    }

    // Sanity check: exactly kResamples estimates, all finite.
    bool all_finite = true;
    for (double v : bootstrap_estimates) {
        if (!std::isfinite(v)) { all_finite = false; break; }
    }
    bool row_count_ok = (static_cast<int>(bootstrap_estimates.size()) == kResamples);

    double ci_lo = std::numeric_limits<double>::quiet_NaN();
    double ci_hi = std::numeric_limits<double>::quiet_NaN();
    double boot_mean = std::numeric_limits<double>::quiet_NaN();
    double boot_sd = std::numeric_limits<double>::quiet_NaN();

    if (all_finite && row_count_ok) {
        std::vector<double> sorted_estimates = bootstrap_estimates;
        std::sort(sorted_estimates.begin(), sorted_estimates.end());
        ci_lo = interpolated_percentile(sorted_estimates, 0.025);
        ci_hi = interpolated_percentile(sorted_estimates, 0.975);

        double sum = 0;
        for (double v : bootstrap_estimates) sum += v;
        boot_mean = sum / kResamples;
        double ss = 0;
        for (double v : bootstrap_estimates) ss += (v - boot_mean) * (v - boot_mean);
        boot_sd = std::sqrt(ss / (kResamples - 1));
    }

    printf("\nSummary across %d bootstrap resamples:\n", kResamples);
    printf("  Original ECAi_96th:              %.6f L/s/person\n", original_percentile);
    printf("  Bootstrap mean of estimates:     %.6f L/s/person\n", boot_mean);
    printf("  Bootstrap std dev of estimates:  %.6f L/s/person\n", boot_sd);
    printf("  95%% percentile bootstrap CI:      [%.6f, %.6f]\n", ci_lo, ci_hi);
    printf("\nSanity check: %d/%d resamples present, all finite = %s\n",
           static_cast<int>(bootstrap_estimates.size()), kResamples,
           all_finite ? "true" : "false");

    const std::string estimates_stem = is_default_category
        ? "patient_ecai_96th_percentile_bootstrap_estimates"
        : category + "_ecai_96th_percentile_bootstrap_estimates";
    const std::string summary_stem = is_default_category
        ? "patient_ecai_96th_percentile_bootstrap_summary"
        : category + "_ecai_96th_percentile_bootstrap_summary";

    const std::string csv_path = out_dir + "/" + estimates_stem + ".csv";
    std::ofstream csv(csv_path);
    csv << std::fixed << std::setprecision(9);
    csv << "resample,ECAi_96th_Lps\n";
    for (int b = 0; b < kResamples; b++) {
        csv << (b + 1) << "," << bootstrap_estimates[b] << "\n";
    }
    csv.close();

    const std::string summary_path = out_dir + "/" + summary_stem + ".csv";
    std::ofstream summary(summary_path);
    summary << std::fixed << std::setprecision(9);
    summary << "parameter,value\n";
    summary << "method,bootstrap\n";
    summary << "category," << category << "\n";
    summary << "N," << kN << "\n";
    summary << "resamples," << kResamples << "\n";
    summary << "target_probability," << kTargetProb << "\n";
    summary << "percentile," << kPercentile << "\n";
    summary << "require_infectors," << (kRequireInfectors ? "true" : "false") << "\n";
    summary << "original_ecai_96th_Lps," << original_percentile << "\n";
    summary << "bootstrap_mean_Lps," << boot_mean << "\n";
    summary << "bootstrap_stddev_Lps," << boot_sd << "\n";
    summary << "ci95_lower_Lps," << ci_lo << "\n";
    summary << "ci95_upper_Lps," << ci_hi << "\n";
    summary << "row_count_ok," << (row_count_ok ? "true" : "false") << "\n";
    summary << "all_finite," << (all_finite ? "true" : "false") << "\n";
    summary.close();

    printf("\nResults written to:\n  %s\n  %s\n", csv_path.c_str(), summary_path.c_str());

    return 0;
}
