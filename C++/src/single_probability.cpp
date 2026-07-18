// single_probability.cpp — port of singleProbability.py:
// Run infection probability simulations for a single occupancy category

#include "model.h"
#include "random_manager.h"

#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <map>
#include <cmath>

int main(int argc, char* argv[]) {
    int N = 10000;
    std::string category = "Classroom";
    double community_rate = -1;  // -1 means use default
    bool allow_zero_infectors = true;
    QERDistribution qer_distribution{};
    double ecai_override = -1;   // -1 means use default from params

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--N" && i + 1 < argc) {
            N = std::atoi(argv[++i]);
        } else if (arg == "--category" && i + 1 < argc) {
            category = argv[++i];
        } else if (arg == "--community_rate" && i + 1 < argc) {
            community_rate = std::atof(argv[++i]);
        } else if (arg == "--ecai" && i + 1 < argc) {
            ecai_override = std::atof(argv[++i]);
        } else if (arg == "--no_zero_infectors") {
            allow_zero_infectors = false;
        } else if (arg == "--use-fitted-qer") {
            qer_distribution = {QERDistributionKind::JonesFitted, {}};
        } else if (arg == "--qer-distribution" && i + 1 < argc) {
            if (!parse_qer_distribution(argv[++i], qer_distribution)) {
                fprintf(stderr, "Unknown QER distribution. See README for valid profile names.\n");
                return 2;
            }
        } else if (arg == "--show_plots") {
            // Ignored — no plotting in C++
        } else if (arg == "--help" || arg == "-h") {
            printf("Usage: single_probability [options]\n");
            printf("  --N <int>            Number of simulations (default 10000)\n");
            printf("  --category <name>    Occupancy category (default Classroom)\n");
            printf("  --community_rate <f> Community infection rate 0-1 (default: from params)\n");
            printf("  --ecai <float>       ECAi in L/s/person (default: from params)\n");
            printf("  --no_zero_infectors  Disallow zero infector simulations\n");
            printf("  --qer-distribution   jones, jones-fitted, or one of the 30\n");
            printf("                        mikszewski-<organism>-<activity> profiles\n");
            printf("                        (default jones)\n");
            printf("  --use-fitted-qer     Alias for --qer-distribution fitted\n");
            return 0;
        }
    }

    if (N < 1) { printf("Invalid N, using default 10000.\n"); N = 10000; }

    auto it = occupancy_params.find(category);
    if (it == occupancy_params.end()) {
        printf("Error: '%s' is not a valid category.\n", category.c_str());
        printf("Available categories:");
        for (const auto& k : category_order)
            printf(" %s", k.c_str());
        printf("\n");
        return 1;
    }

    double ECAi = (ecai_override >= 0) ? ecai_override : it->second.ECAi;
    double comm_rate = (community_rate >= 0) ? community_rate : 0;  // 0 means use default
    bool use_default_comm = (community_rate < 0);

    RandomNumberManager rng;
    std::vector<double> probabilities;
    probabilities.reserve(N);
    int zero_infectors_count = 0;

    for (int i = 0; i < N; i++) {
        auto par = sample_parameters(rng, category);
        double cr = use_default_comm ? 0 : comm_rate;
        auto [prob, infected_flag] = infection_probability(
            ECAi, par, rng, category, !allow_zero_infectors, cr,
 qer_distribution);
        probabilities.push_back(prob);
        if (infected_flag == 0) zero_infectors_count++;
    }

    std::sort(probabilities.begin(), probabilities.end());
    double rank = 0.96 * (N - 1);
    int lo = static_cast<int>(std::floor(rank));
    int hi = static_cast<int>(std::ceil(rank));
    double frac = rank - lo;
    double p96 = (probabilities[lo] * (1 - frac) + probabilities[hi] * frac) * 100;
    double zero_percent = 100.0 * zero_infectors_count / N;

    printf("\nSimulation results for category '%s':\n", category.c_str());
    printf("  Simulations: %d\n", N);
    printf("  ECAi: %g L/s/person\n", ECAi);
    if (use_default_comm) {
        printf("  Community infection rate: Default (%.4f)\n", it->second.community_rate);
    } else {
        printf("  Community infection rate: %.4f\n", comm_rate);
    }
    printf("  Allow zero infectors: %s\n", allow_zero_infectors ? "true" : "false");
    printf("  96th percentile of infection probability: %.3f%%\n", p96);
    printf("  Percentage of simulations with zero infectors: %.1f%%\n\n", zero_percent);

    return 0;
}