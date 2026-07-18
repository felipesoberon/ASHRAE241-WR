// probability_scan.cpp — port of probabilityScan.py:
// Scan ECAi values 5..100 (step 5), find minimum where 96th percentile < 0.1%

#include "model.h"
#include "random_manager.h"

#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <map>

int main(int argc, char* argv[]) {
    int N = 10000;
    QERDistribution qer_distribution{};
    double general_rate = -1;        // -1 = use each category's default (1%)
    double healthcare_rate = -1;     // -1 = use each category's default (3%)

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--use-fitted-qer") {
            qer_distribution = {QERDistributionKind::JonesFitted, {}};
        } else if (arg == "--qer-distribution" && i + 1 < argc) {
            if (!parse_qer_distribution(argv[++i], qer_distribution)) {
                fprintf(stderr, "Unknown QER distribution. See README for valid profile names.\n");
                return 2;
            }
        } else if (arg == "--community-rate-general" && i + 1 < argc) {
            general_rate = std::atof(argv[++i]);
        } else if (arg == "--community-rate-healthcare" && i + 1 < argc) {
            healthcare_rate = std::atof(argv[++i]);
        } else if (arg[0] != '-') {
            N = std::atoi(arg.c_str());
        }
    }
    if (N < 1) {
        printf("Invalid value for N, using default of 10000.\n");
        N = 10000;
    }

    double target_prob = 0.001;

    printf("\nSearching for ECAi (L/s/person) where 96th percentile of probability < %.2f%% for each category (N=%d):\n\n",
           target_prob * 100, N);
    printf("| %-20s | %28s |\n", "Category", "Min ECAi (L/s/person)");
    printf("|%s|%s|\n", std::string(22, '-').c_str(), std::string(30, '-').c_str());

    RandomNumberManager rng;
    struct Row {
        std::string category;
        std::string ecai_str;
    };
    std::vector<Row> results;

    for (const auto& category : category_order) {
        // Resolve this category's community rate from the group overrides;
        // -1 (unset) falls through to the category default (no impact).
        double cr = is_healthcare_category(category) ? healthcare_rate : general_rate;
        bool found = false;
        for (int ECAi = 5; ECAi <= 100; ECAi += 5) {
            std::vector<double> probabilities;
            probabilities.reserve(N);
            for (int i = 0; i < N; i++) {
                auto par = sample_parameters(rng, category);
                auto [prob, _f] = infection_probability(
                    static_cast<double>(ECAi), par, rng, category, false, cr,
                    qer_distribution);
                probabilities.push_back(prob);
            }
            std::sort(probabilities.begin(), probabilities.end());
            double rank = 0.96 * (N - 1);
            int lo = static_cast<int>(std::floor(rank));
            int hi = static_cast<int>(std::ceil(rank));
            double frac = rank - lo;
            double p96 = probabilities[lo] * (1 - frac) + probabilities[hi] * frac;

            if (p96 < target_prob) {
                printf("| %-20s | %28.2f |\n", category.c_str(), static_cast<double>(ECAi));
                results.push_back({category, std::to_string(ECAi)});
                found = true;
                break;
            }
        }
        if (!found) {
            printf("| %-20s | %28s |\n", category.c_str(), "Not found");
            results.push_back({category, "Not found"});
        }
    }

    std::string csv_filename = "ecai_min_for_p_lt_0.1pct.csv";
    std::ofstream csv(csv_filename);
    csv << "Category,ECAi_Lps_for_P_lt_0.1pct\n";
    for (const auto& r : results)
        csv << r.category << "," << r.ecai_str << "\n";
    csv.close();
    printf("\nSummary table saved to %s\n", csv_filename.c_str());

    return 0;
}