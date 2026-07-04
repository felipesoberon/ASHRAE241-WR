// test_model.cpp — unit tests for occupancy parameters and model functions

#include "model.h"
#include "random_manager.h"
#include <cstdio>
#include <cmath>

static int failures = 0;

int main() {
    // ---- occupancy_params: count and spot checks ----
    if (occupancy_params.size() != 25) {
        printf("FAIL: expected 25 categories, got %zu\n", occupancy_params.size());
        failures++;
    }

    const auto& c = occupancy_params.at("Classroom");
    if (c.PBR_GM != 0.55 || c.I0 != 30 || c.volume_m3 != 320 || c.ECAi != 20) {
        printf("FAIL: Classroom params mismatch\n"); failures++;
    }

    const auto& g = occupancy_params.at("Gym");
    if (g.PBR_GM != 0.62 || g.I0 != 180 || g.volume_m3 != 1900 || g.ECAi != 40) {
        printf("FAIL: Gym params mismatch\n"); failures++;
    }

    const auto& dw = occupancy_params.at("Dwelling");
    if (dw.PBR_GSD != 4.3 || dw.Cdrop_GM != 40e4 || dw.I0 != 6 || dw.ECAi != 15) {
        printf("FAIL: Dwelling params mismatch\n"); failures++;
    }

    const auto& ex = occupancy_params.at("Exam");
    if (ex.mask_eff != 0.3 || ex.community_rate != 0.03 || ex.I0 != 3 || ex.volume_m3 != 41) {
        printf("FAIL: Exam params mismatch\n"); failures++;
    }

    const auto& mu = occupancy_params.at("Museum");
    if (mu.I0 != 400 || mu.volume_m3 != 10000 || mu.ECAi != 30) {
        printf("FAIL: Museum params mismatch\n"); failures++;
    }

    // ---- sample_parameters ----
    RandomNumberManager rng;
    auto par = sample_parameters(rng, "Classroom");
    if (par.D != 1.0) { printf("FAIL: D != 1.0\n"); failures++; }
    if (par.VOL != 320) { printf("FAIL: VOL != 320\n"); failures++; }
    if (par.I0 != 30) { printf("FAIL: I0 != 30\n"); failures++; }
    if (par.mask_efficiency != 0.0) { printf("FAIL: mask_eff != 0\n"); failures++; }
    if (par.community_rate != 0.01) { printf("FAIL: community_rate != 0.01\n"); failures++; }
    if (par.lambda_bio <= 0) { printf("FAIL: lambda_bio <= 0\n"); failures++; }
    if (par.gamma < 0.42 || par.gamma > 0.61) { printf("FAIL: gamma range\n"); failures++; }
    if (par.PBR <= 0) { printf("FAIL: PBR <= 0\n"); failures++; }

    // ---- QER ----
    double q = QER(rng, "Classroom");
    if (q <= 0 || !std::isfinite(q)) { printf("FAIL: QER = %f\n", q); failures++; }

    // ---- infection_probability ----
    auto [prob, flag] = infection_probability(20, par, rng, "Classroom", false, 0);
    if (prob < 0 || prob >= 1 || !std::isfinite(prob)) {
        printf("FAIL: P = %f\n", prob); failures++;
    }

    // ---- compute_ECAi ----
    auto [ecai, flag2] = compute_ECAi(par, 0.001, rng, "Classroom", false, 0);
    if (ecai < 0 || !std::isfinite(ecai)) {
        printf("FAIL: ECAi = %f\n", ecai); failures++;
    }

    // ---- infection_probability with require_infectors ----
    // With require_infectors=true, flag should always be 1
    for (int i = 0; i < 100; i++) {
        auto [p, f] = infection_probability(20, par, rng, "Classroom", true, 0);
        if (f != 1) { printf("FAIL: require_infectors gave flag=0\n"); failures++; break; }
    }

    // ---- compute_ECAi with require_infectors ----
    for (int i = 0; i < 100; i++) {
        auto [e, f] = compute_ECAi(par, 0.001, rng, "Classroom", true, 0);
        if (f != 1) { printf("FAIL: require_infectors (ECAi) gave flag=0\n"); failures++; break; }
    }

    // ---- Check mask_eff for healthcare categories ----
    auto par_exam = sample_parameters(rng, "Exam");
    if (par_exam.mask_efficiency != 0.3) {
        printf("FAIL: Exam mask_eff != 0.3\n"); failures++;
    }

    if (failures == 0) {
        printf("PASS: all model tests\n");
        return 0;
    } else {
        printf("FAIL: %d test(s) failed\n", failures);
        return 1;
    }
}