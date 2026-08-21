// analysis/patient_ecai_monte_carlo_study.cpp
//
// Narrowly-scoped statistical study (NOT a unit test, NOT wired into CTest):
// For the "Patient" occupancy category, run 25 independent replicates. Each
// replicate constructs its own RandomNumberManager(N=10000) and performs
// exactly 10,000 simulations using the same model call sequence as
// ecai.cpp's main loop (sample_parameters followed by compute_ECAi, target
// infection probability 0.001, require_infectors=false, category-default
// community rate, default Jones QER distribution), collects the 10,000
// exact (non-rounded) ECAi values, sorts them, and computes the 96th
// percentile via the same linear interpolation ecai.cpp uses. The replicate
// result is that interpolated continuous ECAi value.
//
// Across the 25 replicate ECAi values this program reports sample mean,
// sample standard deviation, standard error, and a two-sided 95% Student-t
// confidence interval (df = 24).
//
// This file only reads from the existing ashrae241 library (model.h /
// random_manager.h); it does not modify simulation logic, defaults, the
// RandomNumberManager, or any of the four production executables.

#include "model.h"
#include "random_manager.h"

#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>
#include <limits>
#include <filesystem>
#include <iomanip>

namespace {

constexpr const char* kCategory = "Patient";
constexpr int    kN            = 10000;
constexpr int    kReplicates   = 25;
constexpr double kTargetProb   = 0.001;
constexpr double kPercentile   = 0.96;
constexpr bool   kRequireInfectors = false; // allows zero-infector draws (existing default)

// Two-sided 95% Student-t critical value for df = 24 (standard table value).
constexpr double kT_0975_df24 = 2.063898961579403;

// Runs one replicate: N simulations of the Patient category, each producing
// one exact ECAi value via sample_parameters + compute_ECAi (same call
// sequence and argument semantics as ecai.cpp's main loop). Returns the
// interpolated 96th percentile of the N exact ECAi values.
double run_replicate() {
    RandomNumberManager rng(static_cast<std::size_t>(kN));
    std::vector<double> ecai_list;
    ecai_list.reserve(kN);
    for (int i = 0; i < kN; i++) {
        auto par = sample_parameters(rng, kCategory);
        auto [ecai_val, _infected_flag] = compute_ECAi(
            par, kTargetProb, rng, kCategory,
            kRequireInfectors, /*override_community_rate=*/-1,
            QERDistribution{});
        ecai_list.push_back(ecai_val);
    }

    std::sort(ecai_list.begin(), ecai_list.end());
    double rank = kPercentile * (kN - 1);
    int lo = static_cast<int>(std::floor(rank));
    int hi = static_cast<int>(std::ceil(rank));
    double frac = rank - lo;
    return ecai_list[lo] * (1 - frac) + ecai_list[hi] * frac;
}

} // namespace

int main() {
    printf("\nPatient ECAi Monte Carlo study (separate from production executables/tests)\n");
    printf("  Category:            %s\n", kCategory);
    printf("  N per replicate:      %d\n", kN);
    printf("  Replicates:           %d\n", kReplicates);
    printf("  Target probability:   %.4f (%.2f%%)\n", kTargetProb, kTargetProb * 100);
    printf("  Percentile:           %.2f\n", kPercentile);
    printf("  Infector mode:        require_infectors=%s (zero-infector draws allowed)\n",
           kRequireInfectors ? "true" : "false");
    printf("  Community rate:       category default (Patient)\n\n");

    std::vector<double> ecai_values;
    ecai_values.reserve(kReplicates);
    for (int r = 0; r < kReplicates; r++) {
        double ecai = run_replicate();
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

        double margin = kT_0975_df24 * se;
        ci_lo = mean - margin;
        ci_hi = mean + margin;
    }

    printf("\nSummary across %d replicates:\n", kReplicates);
    printf("  Sample mean:                    %.6f L/s/person\n", mean);
    printf("  Sample standard deviation:      %.6f L/s/person\n", sd);
    printf("  Standard error:                 %.6f L/s/person\n", se);
    printf("  95%% CI (two-sided t, df=24):    [%.6f, %.6f]\n", ci_lo, ci_hi);
    printf("\nSanity check: %d/%d rows present, all finite = %s\n",
           static_cast<int>(ecai_values.size()), kReplicates,
           all_finite ? "true" : "false");

    // Write results artifact.
    const std::string out_dir = "patient_ecai_monte_carlo_study_results";
    std::filesystem::create_directories(out_dir);

    const std::string csv_path = out_dir + "/patient_ecai_96th_percentile_replicates.csv";
    std::ofstream csv(csv_path);
    csv << std::fixed << std::setprecision(9);
    csv << "replicate,ECAi_96th_Lps\n";
    for (int r = 0; r < kReplicates; r++) {
        csv << (r + 1) << "," << ecai_values[r] << "\n";
    }
    csv.close();

    const std::string summary_path = out_dir + "/patient_ecai_96th_percentile_summary.csv";
    std::ofstream summary(summary_path);
    summary << std::fixed << std::setprecision(9);
    summary << "parameter,value\n";
    summary << "category," << kCategory << "\n";
    summary << "N," << kN << "\n";
    summary << "replicates," << kReplicates << "\n";
    summary << "target_probability," << kTargetProb << "\n";
    summary << "percentile," << kPercentile << "\n";
    summary << "require_infectors," << (kRequireInfectors ? "true" : "false") << "\n";
    summary << "sample_mean_Lps," << mean << "\n";
    summary << "sample_stddev_Lps," << sd << "\n";
    summary << "standard_error_Lps," << se << "\n";
    summary << "t_critical_df24_0975," << kT_0975_df24 << "\n";
    summary << "ci95_lower_Lps," << ci_lo << "\n";
    summary << "ci95_upper_Lps," << ci_hi << "\n";
    summary << "row_count_ok," << (row_count_ok ? "true" : "false") << "\n";
    summary << "all_finite," << (all_finite ? "true" : "false") << "\n";
    summary.close();

    printf("\nResults written to:\n  %s\n  %s\n", csv_path.c_str(), summary_path.c_str());

    return 0;
}
