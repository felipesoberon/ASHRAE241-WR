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
// reproduction mode: a dedicated path for the reproduction study. It runs R independent
// replicates (default 500) of N simulations (default 10,000) per occupancy
// category, and for each replicate records the interpolated ECAi at every
// percentile in --percentiles (default 0.96,0.963). Unlike the other three
// modes it takes explicit community infection rates for the general and
// healthcare groups (defaults 0.01 and 0.03125), applies a category-specific
// mask condition (--mask-preset / --mask-efficiency), and seeds each
// replicate deterministically (--seed base; replicate r uses base + r).
// Across the R replicates it reports the exact unrounded mean, sample
// standard deviation, standard error and two-sided 95% Student-t CI for each
// requested percentile -- no 5 L/s/person ceiling is applied anywhere. It
// writes to its own results directory with its own filenames, so the
// replicate/bootstrap/order-statistic outputs are untouched.
//
// IMPORTANT -- LHS scheme. Reproduction mode uses the existing
// per-distribution-family LHS; there is no alternative sampling scheme.
// N is configurable via --simulations N (exact simulation count in all
// modes); R via --repeats R (replicate and reproduction modes only); B
// (bootstrap mode only) via --resamples B. --repeats and --resamples are not valid outside
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
#include <cerrno>
#include <cstdint>

namespace {

constexpr const char* kDefaultCategory = "Patient";
constexpr const char* kDefaultMethod   = "replicate";
constexpr int    kDefaultN          = 10000;
constexpr int    kDefaultReplicates = 25;
constexpr int    kDefaultResamples  = 2000;
constexpr int    kReproDefaultReplicates = 500;
constexpr double kReproDefaultRateGeneral    = 0.01;
constexpr double kReproDefaultRateHealthcare = 0.03125;
constexpr const char* kReproDefaultPercentiles = "0.96,0.963";
constexpr const char* kReproOutDir = "ecai_reproduction_study_results";
// One-line description of the sampling design actually used, recorded in
// every reproduction output so the deviation from the notebook is explicit.
constexpr const char* kFamilyLhsSchemeNote =
    "repository per-distribution-family LHS (one stratified buffer per distribution family)";
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

// ---------------------------------------------------------------------------
// Reproduction mode helpers
// ---------------------------------------------------------------------------

// Settings for one reproduction run. Every field is CLI-configurable; the
// defaults are the reproduction settings.
struct ReproductionConfig {
    int    n = kDefaultN;
    double community_rate_general     = kReproDefaultRateGeneral;
    double community_rate_healthcare  = kReproDefaultRateHealthcare;
    // mask_override=true forces mask_efficiency on every category (this is
    // how --mask-efficiency 0 is expressed). Otherwise greg_email_masks
    // selects between Greg's email condition and the model.h defaults.
    bool   mask_override    = false;
    double mask_efficiency  = 0.0;
    bool   greg_email_masks = true;
    std::uint64_t base_seed = 0;
};

// The reproduction email condition: no mask for Patient and Resident rooms; Exam,
// Group treatment and Waiting rooms keep their model.h mask efficiency.
// General (non-healthcare) categories have a 0 default already.
bool is_greg_unmasked_healthcare(const std::string& category) {
    return category == "Patient" || category == "Resident";
}

double resolve_mask_efficiency(const ReproductionConfig& cfg,
                               const std::string& category) {
    if (cfg.mask_override) return cfg.mask_efficiency;
    if (cfg.greg_email_masks && is_greg_unmasked_healthcare(category)) return 0.0;
    return get_occupancy_parameters(category).mask_eff;
}

double resolve_reproduction_rate(const ReproductionConfig& cfg,
                                 const std::string& category) {
    return is_healthcare_category(category) ? cfg.community_rate_healthcare
                                            : cfg.community_rate_general;
}

// One reproduction replicate: cfg.n simulations of `category` from a
// RandomNumberManager seeded with `seed`, with the community rate and mask
// efficiency overridden per the configuration. Returns the exact
// (non-rounded) ECAi values, ascending-sorted. The model call sequence is
// the same sample_parameters -> compute_ECAi path the other modes use; only
// par.mask_efficiency (a plain struct field) and the community-rate
// override argument differ, so no model API is changed.
std::vector<double> generate_reproduction_values(const ReproductionConfig& cfg,
                                                 const std::string& category,
                                                 std::uint64_t seed) {
    RandomNumberManager rng(static_cast<std::size_t>(cfg.n), seed);
    const double comm_rate = resolve_reproduction_rate(cfg, category);
    const double mask_eff  = resolve_mask_efficiency(cfg, category);
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(cfg.n));
    for (int i = 0; i < cfg.n; i++) {
        SimParameters par = sample_parameters(rng, category);
        par.mask_efficiency = mask_eff;
        auto [ecai_val, _infected_flag] = compute_ECAi(
            par, kTargetProb, rng, category,
            kRequireInfectors, comm_rate, QERDistribution{});
        values.push_back(ecai_val);
    }
    std::sort(values.begin(), values.end());
    return values;
}

// Compact tag for a percentile fraction, used in CSV column names:
// 0.96 -> "0.96", 0.963 -> "0.963".
std::string percentile_tag(double p) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%g", p);
    return std::string(buf);
}

void print_usage(const char* prog) {
    std::fprintf(stderr,
        "Usage: %s [--method replicate|bootstrap|order-statistic|reproduction]\n"
        "           [--simulations N] [--repeats R] [--resamples B] [--category NAME]\n"
        "           [--percentiles LIST] [--seed S] [--mask-preset PRESET]\n"
        "           [--mask-efficiency V] [--community-rate-general V]\n"
        "           [--community-rate-healthcare V]\n"
        "  --method replicate|bootstrap|order-statistic|reproduction\n"
        "                    Estimation method (default %s)\n"
        "  --simulations N   Simulations per replicate; original sample size in\n"
        "                    bootstrap mode; exact simulation count in\n"
        "                    order-statistic mode (positive integer, default %d)\n"
        "  --repeats R       Number of replicates; replicate and reproduction modes\n"
        "                    only (positive integer >= 2, default %d in replicate\n"
        "                    mode, %d in reproduction mode)\n"
        "  --resamples B     Number of bootstrap resamples; bootstrap mode only\n"
        "                    (positive integer >= 2, default %d)\n"
        "  --category NAME   Occupancy category (default %s); reproduction mode also\n"
        "                    accepts 'all' to sweep every category in order\n",
        prog, kDefaultMethod, kDefaultN, kDefaultReplicates,
        kReproDefaultReplicates, kDefaultResamples, kDefaultCategory);
    std::fprintf(stderr,
        "\nReproduction mode only (rejected with any other --method):\n"
        "  --percentiles LIST          Comma-separated percentile fractions in (0,1),\n"
        "                              e.g. 0.96,0.963 (default %s)\n"
        "  --seed S                    Base seed, non-negative integer (default 0);\n"
        "                              replicate r uses seed S + r\n"
        "  --mask-preset PRESET        greg-email (default): no mask for Patient and\n"
        "                              Resident, model.h default elsewhere;\n"
        "                              category-default: model.h values unchanged\n"
        "  --mask-efficiency V         Force mask efficiency V in [0,1] for every\n"
        "                              category; mutually exclusive with --mask-preset\n"
        "  --community-rate-general V     Rate in [0,1] for general spaces (default %g)\n"
        "  --community-rate-healthcare V  Rate in [0,1] for healthcare spaces\n"
        "                                 (Exam, Group, Patient, Resident, Waiting;\n"
        "                                 default %g)\n",
        kReproDefaultPercentiles, kReproDefaultRateGeneral,
        kReproDefaultRateHealthcare);
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

// Strict double parser: the whole string must be consumed and the value
// must be finite. Returns true on success.
bool parse_double_strict(const char* s, double* out) {
    if (s == nullptr || s[0] == '\0') return false;
    char* end = nullptr;
    errno = 0;
    double v = std::strtod(s, &end);
    if (end == s || *end != '\0') return false;
    if (errno == ERANGE || !std::isfinite(v)) return false;
    *out = v;
    return true;
}

// Strict unsigned 64-bit parser (for --seed). Rejects signs and trailing text.
bool parse_uint64_strict(const char* s, std::uint64_t* out) {
    if (s == nullptr || s[0] == '\0') return false;
    for (const char* c = s; *c != '\0'; c++) {
        if (*c < '0' || *c > '9') return false;
    }
    char* end = nullptr;
    errno = 0;
    unsigned long long v = std::strtoull(s, &end, 10);
    if (end == s || *end != '\0' || errno == ERANGE) return false;
    *out = static_cast<std::uint64_t>(v);
    return true;
}

// Parses a comma-separated list of percentile fractions, each strictly
// inside (0, 1). Rejects empty entries, duplicates, and out-of-range values.
// `error` receives a human-readable reason on failure.
bool parse_percentile_list(const std::string& spec, std::vector<double>* out,
                           std::string* error) {
    out->clear();
    std::size_t pos = 0;
    while (true) {
        std::size_t comma = spec.find(',', pos);
        std::string token = spec.substr(pos, comma == std::string::npos
                                             ? std::string::npos : comma - pos);
        if (token.empty()) {
            *error = "empty percentile entry";
            return false;
        }
        double v = 0;
        if (!parse_double_strict(token.c_str(), &v)) {
            *error = "'" + token + "' is not a number";
            return false;
        }
        if (!(v > 0.0 && v < 1.0)) {
            *error = "'" + token + "' is out of range; percentiles are "
                     "fractions strictly between 0 and 1 (use 0.963, not 96.3)";
            return false;
        }
        for (double existing : *out) {
            if (existing == v) {
                *error = "duplicate percentile '" + token + "'";
                return false;
            }
        }
        out->push_back(v);
        if (comma == std::string::npos) break;
        pos = comma + 1;
    }
    if (out->empty()) {
        *error = "no percentiles given";
        return false;
    }
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

    // Reproduction-mode-only settings (rejected with any other --method).
    std::string percentile_spec = kReproDefaultPercentiles;
    bool percentiles_explicit = false;
    std::uint64_t base_seed = 0;
    bool seed_explicit = false;
    std::string mask_preset = "greg-email";
    bool mask_preset_explicit = false;
    double mask_efficiency = 0.0;
    bool mask_efficiency_explicit = false;
    double rate_general = kReproDefaultRateGeneral;
    bool rate_general_explicit = false;
    double rate_healthcare = kReproDefaultRateHealthcare;
    bool rate_healthcare_explicit = false;

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
        } else if (std::strcmp(argv[i], "--lhs-method") == 0) {
            std::fprintf(stderr, "Error: --lhs-method is no longer supported; reproduction mode always uses the per-distribution-family LHS.\n");
            return 1;
        } else if (std::strcmp(argv[i], "--percentiles") == 0) {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "Error: --percentiles requires an argument.\n");
                print_usage(argv[0]);
                return 1;
            }
            percentile_spec = argv[++i];
            percentiles_explicit = true;
        } else if (std::strcmp(argv[i], "--seed") == 0) {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "Error: --seed requires an argument.\n");
                print_usage(argv[0]);
                return 1;
            }
            if (!parse_uint64_strict(argv[++i], &base_seed)) {
                std::fprintf(stderr, "Error: --seed must be a non-negative integer, got '%s'.\n", argv[i]);
                print_usage(argv[0]);
                return 1;
            }
            seed_explicit = true;
        } else if (std::strcmp(argv[i], "--mask-preset") == 0) {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "Error: --mask-preset requires an argument.\n");
                print_usage(argv[0]);
                return 1;
            }
            mask_preset = argv[++i];
            mask_preset_explicit = true;
        } else if (std::strcmp(argv[i], "--mask-efficiency") == 0) {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "Error: --mask-efficiency requires an argument.\n");
                print_usage(argv[0]);
                return 1;
            }
            if (!parse_double_strict(argv[++i], &mask_efficiency)) {
                std::fprintf(stderr, "Error: --mask-efficiency must be a number, got '%s'.\n", argv[i]);
                print_usage(argv[0]);
                return 1;
            }
            mask_efficiency_explicit = true;
        } else if (std::strcmp(argv[i], "--community-rate-general") == 0) {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "Error: --community-rate-general requires an argument.\n");
                print_usage(argv[0]);
                return 1;
            }
            if (!parse_double_strict(argv[++i], &rate_general)) {
                std::fprintf(stderr, "Error: --community-rate-general must be a number, got '%s'.\n", argv[i]);
                print_usage(argv[0]);
                return 1;
            }
            rate_general_explicit = true;
        } else if (std::strcmp(argv[i], "--community-rate-healthcare") == 0) {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "Error: --community-rate-healthcare requires an argument.\n");
                print_usage(argv[0]);
                return 1;
            }
            if (!parse_double_strict(argv[++i], &rate_healthcare)) {
                std::fprintf(stderr, "Error: --community-rate-healthcare must be a number, got '%s'.\n", argv[i]);
                print_usage(argv[0]);
                return 1;
            }
            rate_healthcare_explicit = true;
        } else if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            std::fprintf(stderr, "Error: unrecognized argument '%s'.\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
    }

    const bool is_reproduction = (method == "reproduction");

    if (method != "replicate" && method != "bootstrap" &&
        method != "order-statistic" && !is_reproduction) {
        std::fprintf(stderr, "Error: --method must be 'replicate', 'bootstrap', 'order-statistic', or 'reproduction', got '%s'.\n", method.c_str());
        print_usage(argv[0]);
        return 1;
    }

    if (method != "replicate" && !is_reproduction && repeats_explicit) {
        std::fprintf(stderr, "Error: --repeats applies to replicate and reproduction modes only; it is not valid with --method %s.\n", method.c_str());
        print_usage(argv[0]);
        return 1;
    }
    if (method != "bootstrap" && resamples_explicit) {
        std::fprintf(stderr, "Error: --resamples applies to bootstrap mode only; it is not valid with --method %s.\n", method.c_str());
        print_usage(argv[0]);
        return 1;
    }

    // Reproduction-only options must not silently do nothing in another mode.
    if (!is_reproduction) {
        struct { bool given; const char* name; } repro_only[] = {
            {percentiles_explicit,     "--percentiles"},
            {seed_explicit,            "--seed"},
            {mask_preset_explicit,     "--mask-preset"},
            {mask_efficiency_explicit, "--mask-efficiency"},
            {rate_general_explicit,    "--community-rate-general"},
            {rate_healthcare_explicit, "--community-rate-healthcare"},
        };
        for (const auto& opt : repro_only) {
            if (opt.given) {
                std::fprintf(stderr, "Error: %s applies to reproduction mode only; it is not valid with --method %s.\n",
                             opt.name, method.c_str());
                print_usage(argv[0]);
                return 1;
            }
        }
    }

    // Reproduction mode uses its own replicate default (Greg's 500).
    if (is_reproduction && !repeats_explicit) {
        replicate_count = kReproDefaultReplicates;
    }

    if ((method == "replicate" || is_reproduction) && replicate_count < 2) {
        std::fprintf(stderr, "Error: --repeats must be at least 2 (need df >= 1 for the Student-t CI), got %d.\n",
                     replicate_count);
        return 1;
    }

    // Reproduction-mode value validation.
    std::vector<double> repro_percentiles;
    if (is_reproduction) {
        if (mask_preset_explicit && mask_efficiency_explicit) {
            std::fprintf(stderr, "Error: --mask-preset and --mask-efficiency are mutually exclusive.\n");
            print_usage(argv[0]);
            return 1;
        }
        if (mask_preset != "greg-email" && mask_preset != "category-default") {
            std::fprintf(stderr, "Error: --mask-preset must be 'greg-email' or 'category-default', got '%s'.\n",
                         mask_preset.c_str());
            print_usage(argv[0]);
            return 1;
        }
        if (mask_efficiency_explicit && !(mask_efficiency >= 0.0 && mask_efficiency <= 1.0)) {
            std::fprintf(stderr, "Error: --mask-efficiency must be in [0, 1], got %g.\n", mask_efficiency);
            return 1;
        }
        if (!(rate_general >= 0.0 && rate_general <= 1.0)) {
            std::fprintf(stderr, "Error: --community-rate-general must be in [0, 1], got %g.\n", rate_general);
            return 1;
        }
        if (!(rate_healthcare >= 0.0 && rate_healthcare <= 1.0)) {
            std::fprintf(stderr, "Error: --community-rate-healthcare must be in [0, 1], got %g.\n", rate_healthcare);
            return 1;
        }
        std::string pct_error;
        if (!parse_percentile_list(percentile_spec, &repro_percentiles, &pct_error)) {
            std::fprintf(stderr, "Error: --percentiles: %s.\n", pct_error.c_str());
            print_usage(argv[0]);
            return 1;
        }
        // Guard against integer overflow of base_seed + (replicate_count - 1).
        if (base_seed > std::numeric_limits<std::uint64_t>::max() -
                        static_cast<std::uint64_t>(replicate_count - 1)) {
            std::fprintf(stderr, "Error: --seed %llu is too large for %d replicates (seed + repeats - 1 overflows).\n",
                         static_cast<unsigned long long>(base_seed), replicate_count);
            return 1;
        }
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
    // Reproduction mode additionally accepts the literal "all", which sweeps
    // every category in category_order.
    const bool repro_all_categories = (is_reproduction && category == "all");
    if (!repro_all_categories &&
        occupancy_params.find(category) == occupancy_params.end()) {
        std::fprintf(stderr, "Error: '%s' is not a valid category.\n", category.c_str());
        std::fprintf(stderr, "Available categories:");
        for (const auto& k : category_order)
            std::fprintf(stderr, " %s", k.c_str());
        if (is_reproduction) std::fprintf(stderr, " all");
        std::fprintf(stderr, "\n");
        return 1;
    }

    const int kN = sim_count;
    const bool is_default_category = (category == kDefaultCategory);

    // --- Reproduction mode -------------------------------------------
    // Handled before the legacy out_dir is created so that the existing
    // ecai_uncertainty_study_results/ directory and filenames are untouched.
    if (is_reproduction) {
        ReproductionConfig cfg;
        cfg.n = kN;
        cfg.community_rate_general    = rate_general;
        cfg.community_rate_healthcare = rate_healthcare;
        cfg.mask_override    = mask_efficiency_explicit;
        cfg.mask_efficiency  = mask_efficiency;
        cfg.greg_email_masks = (mask_preset == "greg-email");
        cfg.base_seed        = base_seed;

        const int kReplicates = replicate_count;
        const int df = kReplicates - 1;
        const double t_crit = student_t_critical_975(df);
        const std::size_t kNumPct = repro_percentiles.size();

        std::vector<std::string> categories;
        if (repro_all_categories) categories = category_order;
        else categories.push_back(category);

        const std::string mask_setting = mask_efficiency_explicit
            ? ("forced " + percentile_tag(mask_efficiency) + " for every category")
            : ("preset " + mask_preset);

        printf("\nECAi reproduction study\n");
        printf("  Method:                    reproduction\n");
        printf("  Categories:                %s\n",
               repro_all_categories ? "all (25, in model order)" : category.c_str());
        printf("  Simulations per replicate: %d\n", kN);
        printf("  Replicates:                %d\n", kReplicates);
        printf("  Seeds:                     %llu .. %llu (replicate r uses seed + r)\n",
               static_cast<unsigned long long>(base_seed),
               static_cast<unsigned long long>(base_seed + kReplicates - 1));
        printf("  Percentiles:              ");
        for (double q : repro_percentiles) printf(" %s", percentile_tag(q).c_str());
        printf("\n");
        printf("  Target probability:        %.4f (%.2f%%)\n", kTargetProb, kTargetProb * 100);
        printf("  Community rate general:    %g\n", rate_general);
        printf("  Community rate healthcare: %g\n", rate_healthcare);
        printf("  Mask efficiency:           %s\n", mask_setting.c_str());
        printf("  QER distribution:          jones (full 8-parameter calculation)\n");
        printf("  Infector mode:             require_infectors=%s (zero-infector draws allowed)\n",
               kRequireInfectors ? "true" : "false");
        printf("  Rounding:                  none (exact unrounded L/s/person; the\n");
        printf("                             5 L/s/person ceiling rule is NOT applied)\n");
        const char* lhs_note = kFamilyLhsSchemeNote;
        printf("  LHS scheme:                %s\n\n", lhs_note);

        const std::string repro_out_dir = kReproOutDir;
        std::filesystem::create_directories(repro_out_dir);

        // Combined per-category x per-percentile summary (also written for a
        // single category, since it is the ASHRAE 241 comparison table).
        const std::string combined_path =
            std::string(repro_out_dir) + "/reproduction_summary.csv";
        std::ofstream combined(combined_path);
        combined << std::fixed << std::setprecision(9);
        combined << "category,percentile,mean_Lps,stddev_Lps,standard_error_Lps,"
                    "ci95_lower_Lps,ci95_upper_Lps,min_Lps,max_Lps,"
                    "ashrae241_published_Lps,diff_vs_published_Lps,"
                    "abs_diff_le_5,mask_efficiency_applied,community_rate_applied,"
                    "simulations,replicates,seed_first,seed_last,all_finite\n";

        for (const std::string& cat : categories) {
            const double mask_applied = resolve_mask_efficiency(cfg, cat);
            const double rate_applied = resolve_reproduction_rate(cfg, cat);
            const double published = get_occupancy_parameters(cat).ECAi;

            // per_pct[j][r] = percentile j of replicate r
            std::vector<std::vector<double>> per_pct(kNumPct);
            for (auto& v : per_pct) v.reserve(static_cast<std::size_t>(kReplicates));

            printf("Category %s (mask %.4g, community rate %g, ASHRAE 241 published %g L/s/p)\n",
                   cat.c_str(), mask_applied, rate_applied, published);

            const int progress_every = std::max(1, kReplicates / 20);
            for (int r = 0; r < kReplicates; r++) {
                const std::uint64_t seed = base_seed + static_cast<std::uint64_t>(r);
                std::vector<double> sorted_values =
                    generate_reproduction_values(cfg, cat, seed);
                std::sort(sorted_values.begin(), sorted_values.end());
                for (std::size_t j = 0; j < kNumPct; j++) {
                    per_pct[j].push_back(
                        interpolated_percentile(sorted_values, repro_percentiles[j]));
                }
                const bool show = (kReplicates <= 25) || ((r + 1) % progress_every == 0)
                                  || (r + 1 == kReplicates);
                if (show) {
                    printf("  Replicate %4d/%d (seed %llu):", r + 1, kReplicates,
                           static_cast<unsigned long long>(seed));
                    for (std::size_t j = 0; j < kNumPct; j++) {
                        printf("  p%s = %.6f", percentile_tag(repro_percentiles[j]).c_str(),
                               per_pct[j].back());
                    }
                    printf(" L/s/person\n");
                }
            }

            bool all_finite = true;
            for (const auto& col : per_pct)
                for (double v : col)
                    if (!std::isfinite(v)) { all_finite = false; break; }
            const bool row_count_ok = std::all_of(
                per_pct.begin(), per_pct.end(),
                [&](const std::vector<double>& col) {
                    return static_cast<int>(col.size()) == kReplicates;
                });

            // Per-replicate CSV: one row per replicate, one column per percentile.
            const std::string replicates_path = repro_out_dir +
                "/reproduction_" + cat + "_replicates.csv";
            std::ofstream rep_csv(replicates_path);
            rep_csv << std::fixed << std::setprecision(9);
            rep_csv << "replicate,seed";
            for (double q : repro_percentiles)
                rep_csv << ",ECAi_p" << percentile_tag(q) << "_Lps";
            rep_csv << "\n";
            for (int r = 0; r < kReplicates; r++) {
                rep_csv << (r + 1) << ","
                        << (base_seed + static_cast<std::uint64_t>(r));
                for (std::size_t j = 0; j < kNumPct; j++)
                    rep_csv << "," << per_pct[j][r];
                rep_csv << "\n";
            }
            rep_csv.close();

            // Per-category summary CSV: run metadata then per-percentile stats.
            const std::string summary_path = repro_out_dir +
                "/reproduction_" + cat + "_summary.csv";
            std::ofstream summary(summary_path);
            summary << "parameter,value\n";
            summary << "method,reproduction\n";
            summary << "category," << cat << "\n";
            summary << "simulations_per_replicate," << kN << "\n";
            summary << "replicates," << kReplicates << "\n";
            summary << "seed_first," << base_seed << "\n";
            summary << "seed_last," << (base_seed + kReplicates - 1) << "\n";
            summary << "percentiles,";
            for (std::size_t j = 0; j < kNumPct; j++)
                summary << (j ? " " : "") << percentile_tag(repro_percentiles[j]);
            summary << "\n";
            summary << "target_probability," << kTargetProb << "\n";
            summary << "require_infectors," << (kRequireInfectors ? "true" : "false") << "\n";
            summary << "qer_distribution,jones\n";
            summary << "mask_setting," << mask_setting << "\n";
            summary << "mask_efficiency_applied," << mask_applied << "\n";
            summary << "community_rate_general," << rate_general << "\n";
            summary << "community_rate_healthcare," << rate_healthcare << "\n";
            summary << "community_rate_applied," << rate_applied << "\n";
            summary << "is_healthcare_category," << (is_healthcare_category(cat) ? "true" : "false") << "\n";
            summary << "ashrae241_published_Lps," << published << "\n";
            summary << "rounding,none (exact unrounded L/s/person)\n";
            summary << "lhs_scheme," << lhs_note << "\n";
            summary << "row_count_ok," << (row_count_ok ? "true" : "false") << "\n";
            summary << "all_finite," << (all_finite ? "true" : "false") << "\n";
            summary << "\n";
            summary << "percentile,mean_Lps,stddev_Lps,standard_error_Lps,"
                       "ci95_lower_Lps,ci95_upper_Lps,min_Lps,max_Lps,"
                       "ashrae241_published_Lps,diff_vs_published_Lps\n";
            summary << std::fixed << std::setprecision(9);

            printf("\n  Summary across %d replicates (exact, unrounded):\n", kReplicates);
            for (std::size_t j = 0; j < kNumPct; j++) {
                const auto& col = per_pct[j];
                double mean = std::numeric_limits<double>::quiet_NaN();
                double sd = mean, se = mean, ci_lo = mean, ci_hi = mean;
                double vmin = mean, vmax = mean;
                if (all_finite && row_count_ok) {
                    double sum = 0;
                    for (double v : col) sum += v;
                    mean = sum / kReplicates;
                    double ss = 0;
                    for (double v : col) ss += (v - mean) * (v - mean);
                    sd = std::sqrt(ss / (kReplicates - 1));
                    se = sd / std::sqrt(static_cast<double>(kReplicates));
                    ci_lo = mean - t_crit * se;
                    ci_hi = mean + t_crit * se;
                    vmin = *std::min_element(col.begin(), col.end());
                    vmax = *std::max_element(col.begin(), col.end());
                }
                const double diff = mean - published;
                const std::string tag = percentile_tag(repro_percentiles[j]);

                printf("    p%-6s mean = %.6f  sd = %.6f  se = %.6f  95%% CI (t, df=%d) = [%.6f, %.6f]\n",
                       tag.c_str(), mean, sd, se, df, ci_lo, ci_hi);
                printf("            min = %.6f  max = %.6f  published = %g  diff = %+.6f L/s/person\n",
                       vmin, vmax, published, diff);

                summary << tag << "," << mean << "," << sd << "," << se << ","
                        << ci_lo << "," << ci_hi << "," << vmin << "," << vmax << ","
                        << published << "," << diff << "\n";

                combined << cat << "," << tag << "," << mean << "," << sd << ","
                         << se << "," << ci_lo << "," << ci_hi << "," << vmin << ","
                         << vmax << "," << published << "," << diff << ","
                         << ((std::isfinite(diff) && std::fabs(diff) <= 5.0) ? "true" : "false")
                         << "," << mask_applied << "," << rate_applied << ","
                         << kN << "," << kReplicates << "," << base_seed << ","
                         << (base_seed + kReplicates - 1) << ","
                         << ((all_finite && row_count_ok) ? "true" : "false") << "\n";
            }
            summary.close();

            printf("  Sanity check: %d replicates per percentile, all finite = %s\n",
                   kReplicates, (all_finite && row_count_ok) ? "true" : "false");
            printf("  Written: %s\n           %s\n\n",
                   replicates_path.c_str(), summary_path.c_str());
        }

        combined.close();
        printf("Combined comparison table written to:\n  %s\n", combined_path.c_str());
        printf("\nNote: %s.\n"
               "This mode uses the repository family-LHS construction; values are\n"
               "not expected to match an external notebook exactly.\n",
               lhs_note);
        return 0;
    }

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
