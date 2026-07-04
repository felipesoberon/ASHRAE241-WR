// analysis/percentiles.cpp — port of analysis/percentiles.py:
// Percentile statistics and threshold search from binary raw data

#include "model.h"
#include "random_manager.h"

#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <utility>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <sstream>

// Ordered list of (category, values) — preserves file insertion order.
using DataSet = std::vector<std::pair<std::string, std::vector<double>>>;

static DataSet load_binary(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        fprintf(stderr, "Cannot open %s\n", path.c_str());
        exit(1);
    }
    uint32_t magic;
    in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    if (magic != 0x41534852) {
        fprintf(stderr, "Bad magic in %s\n", path.c_str());
        exit(1);
    }
    uint32_t nc;
    in.read(reinterpret_cast<char*>(&nc), sizeof(nc));
    DataSet data;
    for (uint32_t i = 0; i < nc; i++) {
        uint32_t nl;
        in.read(reinterpret_cast<char*>(&nl), sizeof(nl));
        std::string name(nl, '\0');
        in.read(&name[0], nl);
        uint32_t cnt;
        in.read(reinterpret_cast<char*>(&cnt), sizeof(cnt));
        std::vector<float> arr(cnt);
        in.read(reinterpret_cast<char*>(arr.data()), cnt * sizeof(float));
        data.emplace_back(name, std::vector<double>(arr.begin(), arr.end()));
    }
    return data;
}

static double percentile(const std::vector<double>& sorted, double p) {
    double rank = (p / 100.0) * (sorted.size() - 1);
    int lo = static_cast<int>(std::floor(rank));
    int hi = static_cast<int>(std::ceil(rank));
    double frac = rank - lo;
    return sorted[lo] * (1 - frac) + sorted[hi] * frac;
}

static void report_table(const DataSet& data,
                         const std::vector<int>& percentiles) {
    printf("%-15s", "Category");
    for (int p : percentiles) printf("%12s", (std::to_string(p) + "th (%)").c_str());
    printf("\n");
    printf("%s\n", std::string(15 + 12 * percentiles.size(), '-').c_str());

    for (const auto& [cat, arr] : data) {
        std::vector<double> sorted = arr;
        std::sort(sorted.begin(), sorted.end());
        printf("%-15s", cat.c_str());
        for (int p : percentiles) {
            double val = percentile(sorted, p) * 100;
            printf("%12.3f", val);
        }
        printf("\n");
    }
}

static double worst_at(const DataSet& data,
                       double p, std::string& worst_cat) {
    double worst = -1e9;
    for (const auto& [cat, arr] : data) {
        std::vector<double> sorted = arr;
        std::sort(sorted.begin(), sorted.end());
        double val = percentile(sorted, p) * 100;
        if (val > worst) { worst = val; worst_cat = cat; }
    }
    return worst;
}

static void report_threshold(const DataSet& data,
                              double target_pct) {
    printf("Coarse-to-fine search for highest percentile with all categories < %.2f%%\n\n", target_pct);
    int lo = 25, hi = 95;
    int best = -1;

    for (int step : {20, 10, 5, 1}) {
        std::vector<int> ps;
        for (int p = lo; p <= hi; p += step) ps.push_back(p);
        if (ps.back() != hi) ps.push_back(hi);

        printf("Scan step %d (%d..%d):\n", step, lo, hi);
        printf("%10s %12s %18s %14s\n", "Percentile", "Max P (%)", "Worst category", "All < target?");
        printf("%s\n", std::string(58, '-').c_str());

        std::vector<std::pair<int, bool>> results;
        for (int p : ps) {
            std::string wc;
            double worst = worst_at(data, p, wc);
            bool ok = worst < target_pct;
            results.push_back({p, ok});
            printf("%10d %12.4f %18s %14s\n", p, worst, wc.c_str(), ok ? "True" : "False");
        }
        printf("\n");

        // Find bracketing interval
        int interval_lo = -1, interval_hi = -1;
        for (size_t i = 0; i + 1 < results.size(); i++) {
            if (results[i].second && !results[i+1].second) {
                interval_lo = results[i].first;
                interval_hi = results[i+1].first;
                break;
            }
        }

        if (interval_lo == -1) {
            bool all_ok = true;
            for (auto& [_, ok] : results) if (!ok) { all_ok = false; break; }
            if (all_ok) {
                best = ps.back();
                printf("All categories stay below %.2f%% through the %dth percentile (top of %d..%d range); crossing not reached.\n",
                       target_pct, best, lo, hi);
            } else {
                best = -1;
                printf("Even the %dth percentile exceeds %.2f%%; crossing is below %d.\n", lo, target_pct, lo);
            }
            return;
        }

        lo = interval_lo;
        hi = interval_hi;
        best = lo;
    }
    printf("Highest percentile where ALL categories < %.2f%%: %d\n", target_pct, best);
}

int main(int argc, char* argv[]) {
    std::string infile = "probabilityECAi_raw.bin";
    std::string report = "both";
    std::vector<int> percentiles = {50, 75, 96};
    double target = 0.1;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--report" && i + 1 < argc) {
            report = argv[++i];
        } else if (arg == "--percentiles" && i + 1 < argc) {
            percentiles.clear();
            while (i + 1 < argc && argv[i+1][0] != '-') {
                percentiles.push_back(std::atoi(argv[++i]));
            }
        } else if (arg == "--target" && i + 1 < argc) {
            target = std::atof(argv[++i]);
        } else if (arg[0] != '-') {
            infile = arg;
        }
    }

    auto data = load_binary(infile);

    if (report == "table" || report == "both") {
        report_table(data, percentiles);
    }
    if (report == "both") {
        printf("\n");
    }
    if (report == "threshold" || report == "both") {
        report_threshold(data, target);
    }

    return 0;
}