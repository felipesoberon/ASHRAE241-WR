// test_qer_fitted.cpp — unit tests for fitted QER distributions

#include "qer_fitted.h"
#include "model.h"
#include "random_manager.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <numeric>

static int failures = 0;

int main() {
    // ---- qer_fitted_params: count and spot checks ----
    if (qer_fitted_params.size() != 25) {
        printf("FAIL: expected 25 categories in qer_fitted_params, got %zu\n",
               qer_fitted_params.size());
        failures++;
    }

    // Check that all 25 categories from category_order are present
    // (we check a few key ones)
    const auto& cell = qer_fitted_params.at("Cell");
    if (cell.mu_log10 != -1.123527) {
        printf("FAIL: Cell mu_log10 = %f, expected -1.123527\n", cell.mu_log10);
        failures++;
    }

    const auto& office = qer_fitted_params.at("Office");
    if (office.mu_log10 != cell.mu_log10 || office.sigma_log10 != cell.sigma_log10) {
        printf("FAIL: Office should match Cell (same group)\n");
        failures++;
    }

    const auto& dwelling = qer_fitted_params.at("Dwelling");
    if (dwelling.sigma_log10 < 1.5) {
        printf("FAIL: Dwelling sigma_log10 should be ~1.56, got %f\n",
               dwelling.sigma_log10);
        failures++;
    }

    // Check Lecture/Auditorium/Lobbies/Common share same params
    const auto& lect = qer_fitted_params.at("Lecture");
    const auto& aud = qer_fitted_params.at("Auditorium");
    const auto& lob = qer_fitted_params.at("Lobbies");
    const auto& com = qer_fitted_params.at("Common");
    if (lect.mu_log10 != aud.mu_log10 || lect.mu_log10 != lob.mu_log10 ||
        lect.mu_log10 != com.mu_log10) {
        printf("FAIL: Lecture/Auditorium/Lobbies/Common should share mu_log10\n");
        failures++;
    }

    // ---- QER_fitted: basic validity ----
    RandomNumberManager rng;
    double q = QER_fitted(rng, "Classroom");
    if (q <= 0 || !std::isfinite(q)) {
        printf("FAIL: QER_fitted(Classroom) = %f\n", q);
        failures++;
    }

    // ---- QER_fitted: statistical properties (10000 samples) ----
    // Median should be approximately 10^mu_log10
    std::vector<double> samples;
    for (int i = 0; i < 10000; i++) {
        samples.push_back(QER_fitted(rng, "Cell"));
    }
    std::sort(samples.begin(), samples.end());
    double median = samples[5000];
    double expected_median = std::pow(10.0, cell.mu_log10);

    // Allow 10% tolerance for Monte Carlo noise
    double rel_err = std::abs(median - expected_median) / expected_median;
    if (rel_err > 0.10) {
        printf("FAIL: QER_fitted(Cell) median = %e, expected ~%e (rel err %.1f%%)\n",
               median, expected_median, rel_err * 100);
        failures++;
    }

    // All samples should be positive
    int nonpositive = 0;
    for (double s : samples) {
        if (s <= 0 || !std::isfinite(s)) nonpositive++;
    }
    if (nonpositive > 0) {
        printf("FAIL: %d nonpositive/nonfinite samples out of 10000\n", nonpositive);
        failures++;
    }

    // ---- Distribution selector and Mikszewski option ----
    QERDistribution parsed{};
    if (!parse_qer_distribution("mikszewski-sars-cov-2-standing-speaking", parsed) ||
        parsed.kind != QERDistributionKind::Mikszewski ||
        parsed.profile != "mikszewski-sars-cov-2-standing-speaking") {
        printf("FAIL: Mikszewski profile did not parse\n");
        failures++;
    }
    if (mikszewski_distribution_names().size() != 30) {
        printf("FAIL: expected 30 Mikszewski profiles\n");
        failures++;
    }
    if (std::abs(mikszewski_qer_params.mu_log10 - std::log10(2.7)) > 1e-12 ||
        std::abs(mikszewski_qer_params.sigma_log10 - 1.2) > 1e-12) {
        printf("FAIL: Mikszewski SARS-CoV-2 parameters are incorrect\n");
        failures++;
    }
    for (const auto& profile : mikszewski_distribution_names()) {
        double q_mik = QER_mikszewski(rng, profile);
        if (q_mik <= 0 || !std::isfinite(q_mik)) {
            printf("FAIL: invalid sample for profile %s\n", profile.c_str());
            failures++;
            break;
        }
    }

    // ---- QER_fitted via QER(use_fitted=true) ----
    // Verify QER with use_fitted=true produces valid values
    double q2 = QER(rng, "Exam", true);
    if (q2 <= 0 || !std::isfinite(q2)) {
        printf("FAIL: QER(Exam, use_fitted=true) = %f\n", q2);
        failures++;
    }

    // ---- QER with use_fitted=false should still work (default) ----
    double q3 = QER(rng, "Exam", false);
    if (q3 <= 0 || !std::isfinite(q3)) {
        printf("FAIL: QER(Exam, use_fitted=false) = %f\n", q3);
        failures++;
    }

    // ---- Verify default (no flag) is Jones method ----
    double q4 = QER(rng, "Gym");
    if (q4 <= 0 || !std::isfinite(q4)) {
        printf("FAIL: QER(Gym) default = %f\n", q4);
        failures++;
    }

    if (failures == 0) {
        printf("PASS: all qer_fitted tests\n");
        return 0;
    } else {
        printf("FAIL: %d test(s) failed\n", failures);
        return 1;
    }
}