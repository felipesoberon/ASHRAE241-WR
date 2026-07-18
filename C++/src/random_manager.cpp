// random_manager.cpp — LHS random number engine + inverse CDF implementations

#include "random_manager.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace {
constexpr double PI = 3.14159265358979323846;

// ----------------------------------------------------------------
// Inverse normal CDF — Acklam's algorithm + one Halley refinement
// Max relative error ~1.15e-9 (matches scipy.stats.norm.ppf closely)
// ----------------------------------------------------------------
double acklam_normal_ppf(double u) {
    static const double a[] = {
        -3.969683028665376e+01,  2.209460984245205e+02,
        -2.759285188461852e+02,  1.383577518672690e+02,
        -3.066479806614716e+01,  2.506628277453239e+00
    };
    static const double b[] = {
        -5.447609879822406e+01,  1.615858368580409e+02,
        -1.556989798598866e+02,  6.680131188771972e+01,
        -1.328068155288572e+01
    };
    static const double c[] = {
        -7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141118049e+00,  2.938163982698783e+00
    };
    static const double d[] = {
         7.784695709041462e-03,  3.224671290700398e-01,
         2.445134137142996e+00,  3.754408661907416e+00
    };
    const double plow = 0.02425, phigh = 1.0 - plow;
    double q, r, x;
    if (u < plow) {
        q = std::sqrt(-2.0 * std::log(u));
        x = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
            ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1);
    } else if (u <= phigh) {
        q = u - 0.5;
        r = q * q;
        x = (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q /
            (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1);
    } else {
        q = std::sqrt(-2.0 * std::log(1.0 - u));
        x = -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
            ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1);
    }
    // One Halley refinement step
    double e = 0.5 * std::erfc(-x / std::sqrt(2.0)) - u;
    double u1 = e * std::sqrt(2.0 * PI) * std::exp(x * x / 2.0);
    x = x - u1 / (1.0 + x * u1 / 2.0);
    return x;
}
} // anonymous namespace

double normal_ppf(double u) {
    if (u <= 0.0) return -std::numeric_limits<double>::infinity();
    if (u >= 1.0) return  std::numeric_limits<double>::infinity();
    return acklam_normal_ppf(u);
}

double lognormal_ppf(double u, double mean, double sigma) {
    // scipy lognorm.ppf(u, s=sigma, scale=exp(mean)) = exp(normal_ppf(u)*sigma + mean)
    return std::exp(normal_ppf(u) * sigma + mean);
}

double uniform_ppf(double u, double low, double high) {
    return low + (high - low) * u;
}

double log10normal_ppf(double u, double mu, double sigma) {
    return std::pow(10.0, normal_ppf(u) * sigma + mu);
}

// ----------------------------------------------------------------
// Incomplete beta function (continued fraction — Lentz's algorithm)
// Based on Numerical Recipes betacf
// ----------------------------------------------------------------
namespace {
double betacf(double a, double b, double x) {
    const int MAXIT = 200;
    const double EPS = 3e-16;
    const double FPMIN = 1e-300;
    double qab = a + b, qap = a + 1.0, qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::fabs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    double h = d;
    for (int m = 1, m2 = 2; m <= MAXIT; m++, m2 += 2) {
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (std::fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (std::fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (std::fabs(del - 1.0) < EPS) break;
    }
    return h;
}

double beta_cdf(double x, double a, double b) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    double bt = std::exp(std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b)
                         + a * std::log(x) + b * std::log(1.0 - x));
    if (x < (a + 1.0) / (a + b + 2.0))
        return bt * betacf(a, b, x) / a;
    else
        return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}
} // anonymous namespace

double beta_ppf(double u, double a, double b) {
    if (u <= 0.0) return 0.0;
    if (u >= 1.0) return 1.0;
    // Bisection on beta_cdf — robust, ~60 iterations for double precision
    double lo = 0.0, hi = 1.0, mid;
    for (int i = 0; i < 60; i++) {
        mid = 0.5 * (lo + hi);
        double cdf = beta_cdf(mid, a, b);
        if (cdf < u) lo = mid;
        else hi = mid;
    }
    return 0.5 * (lo + hi);
}

// ----------------------------------------------------------------
// Binomial inverse CDF — sequential CDF search (matches Python logic)
// Uses log-space PMF to avoid overflow for large n
// ----------------------------------------------------------------
int binomial_ppf(double u, int n, double p) {
    if (u <= 0.0) return 0;
    if (u >= 1.0) return n;
    double cdf = 0.0;
    double log_p = std::log(p);
    double log_1mp = std::log(1.0 - p);
    for (int k = 0; k <= n; k++) {
        double log_pmf = std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1)
                         + k * log_p + (n - k) * log_1mp;
        cdf += std::exp(log_pmf);
        if (u <= cdf)
            return k;
    }
    return n;
}

// ----------------------------------------------------------------
// LHS Random Number Manager
// ----------------------------------------------------------------
RandomNumberManager::RandomNumberManager(std::size_t batch_size)
    : batch_size_(batch_size), rng_(std::random_device{}())
{
    for (const auto& dist : {"beta", "normal", "lognormal", "uniform", "binomial", "log10normal", "qer_fitted"}) {
        refill(dist);
    }
}

void RandomNumberManager::refill(const std::string& dist) {
    std::vector<double> buffer(batch_size_);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    for (std::size_t i = 0; i < batch_size_; i++) {
        double lo = static_cast<double>(i) / batch_size_;
        double hi = static_cast<double>(i + 1) / batch_size_;
        buffer[i] = lo + uni(rng_) * (hi - lo);
    }
    std::shuffle(buffer.begin(), buffer.end(), rng_);
    buffers_[dist] = std::move(buffer);
    indices_[dist] = 0;
}

double RandomNumberManager::get(const std::string& dist) {
    auto idx_it = indices_.find(dist);
    if (idx_it == indices_.end() || idx_it->second >= batch_size_)
        refill(dist);
    return buffers_[dist][indices_[dist]++];
}

// ----------------------------------------------------------------
// Distribution sampling functions
// ----------------------------------------------------------------
double random_beta_lhs(RandomNumberManager& rng, double a, double b,
                       double loc, double scale) {
    double u = rng.get("beta");
    return beta_ppf(u, a, b) * scale + loc;
}

double random_normal_lhs(RandomNumberManager& rng, double mu, double sigma) {
    double u = rng.get("normal");
    return normal_ppf(u) * sigma + mu;
}

double random_lognormal_lhs(RandomNumberManager& rng, double mean, double sigma) {
    double u = rng.get("lognormal");
    return std::exp(normal_ppf(u) * sigma + mean);
}

double random_uniform_lhs(RandomNumberManager& rng, double low, double high) {
    double u = rng.get("uniform");
    return low + (high - low) * u;
}

int random_binomial_lhs(RandomNumberManager& rng, int n, double p) {
    double u = rng.get("binomial");
    return binomial_ppf(u, n, p);
}

double random_log10normal_lhs(RandomNumberManager& rng, double mu, double sigma) {
    double u = rng.get("log10normal");
    return std::pow(10.0, normal_ppf(u) * sigma + mu);
}