// probability_ecai.cpp — port of probabilityECAi.py:
// Evaluate 96th percentile infection probability at ASHRAE 241 ECAi values

#include "model.h"
#include "random_manager.h"

#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <map>
#include <cstdint>

// Simple binary format for --save-all:
//   uint32_t magic = 0x41534852 ("ASHR")
//   uint32_t num_categories
//   For each category:
//     uint32_t name_len, char[name_len] name, uint32_t count, double[count] values

static void save_binary(const std::string& path,
                        const std::vector<std::pair<std::string, std::vector<float>>>& data) {
    std::ofstream out(path, std::ios::binary);
    uint32_t magic = 0x41534852;
    out.write(reinterpret_cast<char*>(&magic), sizeof(magic));
    uint32_t nc = static_cast<uint32_t>(data.size());
    out.write(reinterpret_cast<char*>(&nc), sizeof(nc));
    for (const auto& [cat, arr] : data) {
        uint32_t nl = static_cast<uint32_t>(cat.size());
        out.write(reinterpret_cast<char*>(&nl), sizeof(nl));
        out.write(cat.data(), nl);
        uint32_t cnt = static_cast<uint32_t>(arr.size());
        out.write(reinterpret_cast<char*>(&cnt), sizeof(cnt));
        out.write(reinterpret_cast<const char*>(arr.data()), cnt * sizeof(float));
    }
    out.close();
}

int main(int argc, char* argv[]) {
    int N = 10000;
    bool save_all = false;
    bool require_infectors = true;   // default: match original behavior
    double general_rate = -1;        // -1 = use each category's default (1%)
    double healthcare_rate = -1;     // -1 = use each category's default (3%)
    std::string outfile = "probabilityECAi_raw.bin";

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--save-all") {
            save_all = true;
        } else if (arg == "--outfile" && i + 1 < argc) {
            outfile = argv[++i];
        } else if (arg == "--no-require-infectors") {
            require_infectors = false;
        } else if (arg == "--community-rate-general" && i + 1 < argc) {
            general_rate = std::atof(argv[++i]);
        } else if (arg == "--community-rate-healthcare" && i + 1 < argc) {
            healthcare_rate = std::atof(argv[++i]);
        } else if (arg[0] != '-' && N == 10000) {
            // First positional = N
            N = std::atoi(arg.c_str());
        }
    }
    if (N < 1) { printf("Invalid value for N, using default of 10000.\n"); N = 10000; }

    printf("\nEvaluating the 96th percentile probability at ECAi from occupancy_params for each category (N=%d):\n\n", N);
    printf("| %-20s | %15s | %20s |\n", "Category", "ECAi (L/s/p)", "96th per P (%)");
    printf("|%s|%s|%s|\n", std::string(22, '-').c_str(), std::string(17, '-').c_str(), std::string(22, '-').c_str());

    RandomNumberManager rng;
    struct Row {
        std::string category;
        std::string ecai_str;
        std::string p96_str;
    };
    std::vector<Row> results;
    std::vector<std::pair<std::string, std::vector<float>>> raw_data;

    for (const auto& category : category_order) {
        double ECAi = occupancy_params.at(category).ECAi;
        std::vector<double> probabilities;
        probabilities.reserve(N);

        // Resolve this category's community rate from the group overrides;
        // -1 (unset) falls through to the category default (no impact).
        double cr = is_healthcare_category(category) ? healthcare_rate : general_rate;

        for (int i = 0; i < N; i++) {
            auto par = sample_parameters(rng, category);
            auto [prob, _] = infection_probability(ECAi, par, rng, category,
                require_infectors, cr);
            probabilities.push_back(prob);
        }

        std::sort(probabilities.begin(), probabilities.end());
        double rank = 0.96 * (N - 1);
        int lo = static_cast<int>(std::floor(rank));
        int hi = static_cast<int>(std::ceil(rank));
        double frac = rank - lo;
        double percentile_96 = (probabilities[lo] * (1 - frac) + probabilities[hi] * frac) * 100;

        if (save_all) {
            std::vector<float> arr(N);
            for (int i = 0; i < N; i++) arr[i] = static_cast<float>(probabilities[i]);
            raw_data.emplace_back(category, std::move(arr));
        }

        printf("| %-20s | %15.2f | %20.3f |\n", category.c_str(), ECAi, percentile_96);

        std::ostringstream oss_e, oss_p;
        oss_e << std::fixed << std::setprecision(2) << ECAi;
        // Match Python csv.DictWriter: write full float precision for P_96th_percentile
        oss_p << std::setprecision(15) << percentile_96;
        results.push_back({category, oss_e.str(), oss_p.str()});
    }

    // Write CSV
    std::string csv_filename = "ecai_ashrae241_96th_percentile.csv";
    std::ofstream csv(csv_filename);
    csv << "Category,ECAi_Lps,P_96th_percentile\n";
    for (const auto& r : results) {
        csv << r.category << "," << r.ecai_str << "," << r.p96_str << "\n";
    }
    csv.close();
    printf("\nSummary table saved to %s\n", csv_filename.c_str());

    if (save_all) {
        save_binary(outfile, raw_data);
        size_t total = 0;
        for (const auto& [_, arr] : raw_data) total += arr.size();
        printf("Raw data (%zu values, raw probability 0-1) saved to %s\n", total, outfile.c_str());
    }

    return 0;
}