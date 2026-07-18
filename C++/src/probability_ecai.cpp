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
//
// --save-inputs format (same structure, but values are SimInputs structs
// written as raw doubles, plus the sample_parameters fields):
//   uint32_t magic = 0x494E5054 ("INPT")
//   uint32_t num_categories
//   For each category:
//     uint32_t name_len, char[name_len] name, uint32_t count
//     double[count] PBR_sample, lambda_bio, gamma,  // from SimParameters
//     double[count] n_infected, phi, QER_sum, mask_factor, Q, P,
//     double[count] qer_first.PBR_qer, qer_first.C_drop, qer_first.d,
//     double[count] qer_first.E, qer_first.Vdrop, qer_first.GVL_ml,
//     double[count] qer_first.GVL_m3, qer_first.VF, qer_first.RTD,
//     double[count] qer_first.DK, qer_first.VER, qer_first.QER_val
// (17 doubles per simulation: 3 from SimParameters + 6 from SimInputs
//  + 12 from QERInputs = 21 doubles per sim, stored as 21 sequential
//  double arrays for cache-friendly writing)

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

// --- save_inputs_binary ---
// Writes per-simulation input parameters to a binary file.
// 21 fields per simulation, stored as 21 sequential double arrays
// per category (column-major: all N values of field 0, then all N of
// field 1, etc.) for efficient bulk writes.
//
// Field order (21 doubles per sim):
//   0:  PBR_sample    1:  lambda_bio   2:  gamma
//   3:  n_infected    4:  phi           5:  QER_sum
//   6:  mask_factor   7:  Q             8:  P
//   9:  qer.PBR_qer  10:  qer.C_drop   11:  qer.d
//  12:  qer.E        13:  qer.Vdrop    14:  qer.GVL_ml
//  15:  qer.GVL_m3   16:  qer.VF       17:  qer.RTD
//  18:  qer.DK       19:  qer.VER     20:  qer.QER_val

static const int N_INPUT_FIELDS = 21;

static void save_inputs_binary(const std::string& path,
    const std::vector<std::string>& cats,
    const std::vector<std::vector<double>>& inputs)
{
    std::ofstream out(path, std::ios::binary);
    uint32_t magic = 0x494E5054;  // "INPT"
    out.write(reinterpret_cast<char*>(&magic), sizeof(magic));
    uint32_t nc = static_cast<uint32_t>(cats.size());
    out.write(reinterpret_cast<char*>(&nc), sizeof(nc));
    for (size_t c = 0; c < cats.size(); c++) {
        const auto& cat = cats[c];
        uint32_t nl = static_cast<uint32_t>(cat.size());
        out.write(reinterpret_cast<char*>(&nl), sizeof(nl));
        out.write(cat.data(), nl);
        uint32_t count = static_cast<uint32_t>(inputs[c].size() / N_INPUT_FIELDS);
        out.write(reinterpret_cast<char*>(&count), sizeof(count));
        // Write all doubles at once (already laid out as 21 fields x count)
        out.write(reinterpret_cast<const char*>(inputs[c].data()),
                  inputs[c].size() * sizeof(double));
    }
    out.close();
}

int main(int argc, char* argv[]) {
    int N = 10000;
    bool save_all = false;
    bool save_inputs = false;
    bool require_infectors = true;   // default: match original behavior
    QERDistribution qer_distribution = QERDistribution::Jones;
    double general_rate = -1;        // -1 = use each category's default (1%)
    double healthcare_rate = -1;     // -1 = use each category's default (3%)
    std::string outfile = "probabilityECAi_raw.bin";
    std::string inputs_file = "probabilityECAi_inputs.bin";

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--save-all") {
            save_all = true;
        } else if (arg == "--save-inputs") {
            save_inputs = true;
        } else if (arg == "--outfile" && i + 1 < argc) {
            outfile = argv[++i];
        } else if (arg == "--inputs-file" && i + 1 < argc) {
            inputs_file = argv[++i];
        } else if (arg == "--no-require-infectors") {
            require_infectors = false;
        } else if (arg == "--use-fitted-qer") {
            qer_distribution = QERDistribution::JonesFitted;
        } else if (arg == "--qer-distribution" && i + 1 < argc) {
            if (!parse_qer_distribution(argv[++i], qer_distribution)) {
                fprintf(stderr, "Unknown QER distribution. Use jones, fitted, or mikszewski.\n");
                return 2;
            }
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
    // For --save-inputs: store 21 doubles per simulation per category
    std::vector<std::string> inputs_cats;
    std::vector<std::vector<double>> inputs_data;

    for (const auto& category : category_order) {
        double ECAi = occupancy_params.at(category).ECAi;
        std::vector<double> probabilities;
        probabilities.reserve(N);

        // For --save-inputs
        std::vector<double> cat_inputs;
        if (save_inputs)
            cat_inputs.reserve(static_cast<size_t>(N) * N_INPUT_FIELDS);

        // Resolve this category's community rate from the group overrides;
        // -1 (unset) falls through to the category default (no impact).
        double cr = is_healthcare_category(category) ? healthcare_rate : general_rate;

        for (int i = 0; i < N; i++) {
            auto par = sample_parameters(rng, category);
            if (save_inputs) {
                auto result = infection_probability_with_inputs(
                    ECAi, par, rng, category, require_infectors, cr,
                    qer_distribution);
                probabilities.push_back(result.P);

                // Store 21 fields: 3 from SimParameters + 6 from SimInputs
                // + 12 from QERInputs
                const auto& si = result.inputs;
                const auto& qi = si.qer_first;
                double fields[N_INPUT_FIELDS] = {
                    par.PBR, par.lambda_bio, par.gamma,
                    static_cast<double>(si.n_infected), si.phi, si.QER_sum,
                    si.mask_factor, si.Q, si.P,
                    qi.PBR_qer, qi.C_drop, qi.d,
                    qi.E, qi.Vdrop, qi.GVL_ml,
                    qi.GVL_m3, qi.VF, qi.RTD,
                    qi.DK, qi.VER, qi.QER_val
                };
                cat_inputs.insert(cat_inputs.end(), fields, fields + N_INPUT_FIELDS);
            } else {
                auto [prob, _] = infection_probability(ECAi, par, rng, category,
                    require_infectors, cr, qer_distribution);
                probabilities.push_back(prob);
            }
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

        if (save_inputs) {
            inputs_cats.push_back(category);
            inputs_data.push_back(std::move(cat_inputs));
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

    if (save_inputs) {
        save_inputs_binary(inputs_file, inputs_cats, inputs_data);
        size_t total_sims = 0;
        for (const auto& v : inputs_data) total_sims += v.size() / N_INPUT_FIELDS;
        printf("Input parameters (%zu simulations, %d fields each) saved to %s\n",
               total_sims, N_INPUT_FIELDS, inputs_file.c_str());
    }

    return 0;
}