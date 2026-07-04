#pragma once
#include <string>
#include <vector>
#include <map>
#include <random>
#include <utility>
#include <cstdint>

// ------------------------------------------------------------------
// Inverse CDFs (percent point functions) — match scipy.stats.*.ppf
// ------------------------------------------------------------------
double normal_ppf(double u);
double lognormal_ppf(double u, double mean, double sigma);
double uniform_ppf(double u, double low, double high);
double beta_ppf(double u, double a, double b);
double log10normal_ppf(double u, double mu, double sigma);
int    binomial_ppf(double u, int n, double p);

// ------------------------------------------------------------------
// LHS Random Number Manager
// ------------------------------------------------------------------
class RandomNumberManager {
public:
    explicit RandomNumberManager(std::size_t batch_size = 1000000);
    double get(const std::string& dist);

private:
    std::size_t batch_size_;
    std::map<std::string, std::vector<double>> buffers_;
    std::map<std::string, std::size_t> indices_;
    std::mt19937_64 rng_;
    void refill(const std::string& dist);
};

// Distribution sampling functions (match Python API)
double random_beta_lhs(RandomNumberManager& rng, double a, double b,
                       double loc = 0, double scale = 1);
double random_normal_lhs(RandomNumberManager& rng, double mu, double sigma);
double random_lognormal_lhs(RandomNumberManager& rng, double mean, double sigma);
double random_uniform_lhs(RandomNumberManager& rng, double low, double high);
int    random_binomial_lhs(RandomNumberManager& rng, int n, double p);
double random_log10normal_lhs(RandomNumberManager& rng, double mu, double sigma);