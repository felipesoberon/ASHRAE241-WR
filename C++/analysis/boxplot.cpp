// analysis/boxplot.cpp — port of analysis/boxplot.py:
// Text-based horizontal box-and-whisker summary of raw probabilities

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

// Ordered list of (category, values) — preserves file insertion order.
using DataSet = std::vector<std::pair<std::string, std::vector<double>>>;

static DataSet load_binary(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) { fprintf(stderr, "Cannot open %s\n", path.c_str()); exit(1); }
    uint32_t magic;
    in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    if (magic != 0x41534852) { fprintf(stderr, "Bad magic\n"); exit(1); }
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

static std::string bar(double frac, int width) {
    int n = std::max(0, std::min(width, static_cast<int>(frac * width)));
    return std::string(n, '#') + std::string(width - n, '-');
}

int main(int argc, char* argv[]) {
    std::string infile = "probabilityECAi_raw.bin";
    double target = 0.1;
    double xmin = 0.001;
    double xmax = 100;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--target" && i + 1 < argc) target = std::atof(argv[++i]);
        else if (arg == "--xmin" && i + 1 < argc) xmin = std::atof(argv[++i]);
        else if (arg == "--xmax" && i + 1 < argc) xmax = std::atof(argv[++i]);
        else if (arg == "--save" && i + 1 < argc) {
            // Ignored — text output only
            i++;
        } else if (arg[0] != '-') {
            infile = arg;
        }
    }

    auto data = load_binary(infile);

    printf("\nInfection probability distribution by occupancy (ASHRAE 241 ECAi)\n");
    printf("Log scale: %.4f%% .. %.0f%%   target: %.3f%%\n\n", xmin, xmax, target);

    // Log scale range for the text bars
    double log_xmin = std::log10(xmin);
    double log_xmax = std::log10(xmax);
    int bar_width = 50;

    // Header
    printf("%-15s  %8s %8s %8s %8s %8s   %s\n",
           "Category", "Min", "Q1", "Med", "Q3", "P96", "  [log scale bar]");
    printf("%s\n", std::string(95, '-').c_str());

    for (const auto& [cat, arr] : data) {
        std::vector<double> sorted = arr;
        std::sort(sorted.begin(), sorted.end());

        // Convert to percent, drop zero for log display
        std::vector<double> pct;
        for (double v : sorted) {
            double p = v * 100;
            if (p > 0) pct.push_back(p);
        }
        if (pct.empty()) {
            printf("%-15s  (no positive values)\n", cat.c_str());
            continue;
        }

        std::sort(pct.begin(), pct.end());
        double q1 = percentile(pct, 25);
        double med = percentile(pct, 50);
        double q3 = percentile(pct, 75);
        double p96 = percentile(pct, 96);
        double mn = pct.front();

        // Bar: position of p96 on log scale
        double log_p96 = std::log10(std::max(p96, xmin));
        double frac96 = (log_p96 - log_xmin) / (log_xmax - log_xmin);

        printf("%-15s  %7.4f %8.3f %8.3f %8.3f %8.3f   [%s]\n",
               cat.c_str(), mn, q1, med, q3, p96, bar(frac96, bar_width).c_str());
    }

    // Target line
    double log_target = std::log10(target);
    double frac_t = (log_target - log_xmin) / (log_xmax - log_xmin);
    printf("\n  Target: %.3f%%  [%s]\n", target, bar(frac_t, bar_width).c_str());

    return 0;
}