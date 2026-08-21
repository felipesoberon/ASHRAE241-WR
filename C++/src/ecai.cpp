// ecai.cpp — port of ecai.py: compute required ECAi at fixed infection probability target

#include "model.h"
#include "random_manager.h"

#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <map>

static const double LPS_TO_CFM = 2.11888;

int main(int argc, char* argv[]) {
    int N = 10000;
    bool require_infectors = false;  // default: match original behavior
    QERDistribution qer_distribution{};
    double general_rate = -1;        // -1 = use each category's default (1%)
    double healthcare_rate = -1;     // -1 = use each category's default (3%)

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--require-infectors") {
            require_infectors = true;
        } else if (arg == "--use-fitted-qer") {
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
    if (N <= 0) {
        fprintf(stderr, "Error: N must be a positive integer (got %d).\n", N);
        return 1;
    }

    double target_P = 0.001;
    printf("\nCalculating ECAi from fixed infection probability using %d samples for all categories...\n", N);
    printf("\nTarget Probability: %.3f%%\n\n", target_P * 100);
    printf("Round (L/s/p) is the ceiling to nearest multiple of 5.\n");
    printf("Round (CFM/p) is the ceiling to nearest multiple of 10.\n\n");

    std::string line = "+-----------------+---------------+---------------+-----------------+---------------+";
    printf("%s\n", line.c_str());
    printf("| %-15s | %13s | %13s | %15s | %13s |\n",
           "Category", "ECAi (L/s/p)", "Round (L/s/p)", "Round (CFM/p)", "Zero Inf. (%)");
    printf("%s\n", line.c_str());

    RandomNumberManager rng(static_cast<std::size_t>(N));
    struct Row {
        std::string category;
        std::string p96_str;
        int rounded_lps;
        int rounded_cfm;
        double zero_infected_percent;
    };
    std::vector<Row> results;
    int grand_zero_infected_count = 0;

    for (const auto& category : category_order) {
        std::vector<double> ECAi_list;
        ECAi_list.reserve(N);
        int zero_infected_count = 0;

        // Resolve this category's community rate from the group overrides;
        // -1 (unset) falls through to the category default (no impact).
        double cr = is_healthcare_category(category) ? healthcare_rate : general_rate;

        for (int i = 0; i < N; i++) {
            auto par = sample_parameters(rng, category);
            auto [ecai_val, infected_flag] = compute_ECAi(par, target_P, rng, category,
                require_infectors, cr, qer_distribution);
            ECAi_list.push_back(ecai_val);
            if (infected_flag == 0) zero_infected_count++;
        }

        std::sort(ECAi_list.begin(), ECAi_list.end());
        // 96th percentile via linear interpolation (matches numpy.percentile default)
        double percentile_96;
        {
            double rank = 0.96 * (N - 1);
            int lo = static_cast<int>(std::floor(rank));
            int hi = static_cast<int>(std::ceil(rank));
            double frac = rank - lo;
            percentile_96 = ECAi_list[lo] * (1 - frac) + ECAi_list[hi] * frac;
        }

        int rounded_lps = static_cast<int>(std::ceil(percentile_96 / 5.0)) * 5;
        double percentile_96_CFM = percentile_96 * LPS_TO_CFM;
        int rounded_cfm = static_cast<int>(std::ceil(percentile_96_CFM / 10.0)) * 10;
        double zero_infected_percent = 100.0 * zero_infected_count / N;

        printf("| %-15s | %13.2f | %13d | %15d | %12.1f  |\n",
               category.c_str(), percentile_96, rounded_lps, rounded_cfm, zero_infected_percent);

        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2) << percentile_96;
        results.push_back({category, oss.str(), rounded_lps, rounded_cfm, zero_infected_percent});
        grand_zero_infected_count += zero_infected_count;
    }

    printf("%s\n", line.c_str());

    // Grand total — matches Python ecai.py format
    int grand_total_simulations = static_cast<int>(category_order.size()) * N;
    double grand_percent = 100.0 * grand_zero_infected_count / grand_total_simulations;
    printf("\nGrand Total Simulations with Zero Infected: %d of %d (%.1f %%)\n",
           grand_zero_infected_count, grand_total_simulations, grand_percent);

    // Write CSV
    std::string csv_filename = "ecai_results.csv";
    std::ofstream csv(csv_filename);
    csv << "Category,Percentile_96_Lps,Rounded_Lps,Rounded_CFM,ZeroInfectedPercent\n";
    for (const auto& r : results) {
        csv << r.category << "," << r.p96_str << "," << r.rounded_lps
            << "," << r.rounded_cfm << "," << std::fixed << std::setprecision(1)
            << r.zero_infected_percent << "\n";
    }
    csv.close();
    printf("\nResults saved to %s\n", csv_filename.c_str());

    return 0;
}