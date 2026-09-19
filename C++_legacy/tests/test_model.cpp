// test_model.cpp — unit tests for occupancy parameters and model functions

#include "model.h"
#include "random_manager.h"
#include <cstdio>
#include <cmath>
#include <stdexcept>

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
    // override_community_rate = -1 means "unset" (use category default).
    auto [prob, flag] = infection_probability(20, par, rng, "Classroom", false, -1);
    if (prob < 0 || prob >= 1 || !std::isfinite(prob)) {
        printf("FAIL: P = %f\n", prob); failures++;
    }

    // ---- compute_ECAi ----
    auto [ecai, flag2] = compute_ECAi(par, 0.001, rng, "Classroom", false, -1);
    if (ecai < 0 || !std::isfinite(ecai)) {
        printf("FAIL: ECAi = %f\n", ecai); failures++;
    }

    // ---- infection_probability with require_infectors ----
    // With require_infectors=true, flag should always be 1
    for (int i = 0; i < 100; i++) {
        auto [p, f] = infection_probability(20, par, rng, "Classroom", true, -1);
        if (f != 1) { printf("FAIL: require_infectors gave flag=0\n"); failures++; break; }
    }

    // ---- compute_ECAi with require_infectors ----
    for (int i = 0; i < 100; i++) {
        auto [e, f] = compute_ECAi(par, 0.001, rng, "Classroom", true, -1);
        if (f != 1) { printf("FAIL: require_infectors (ECAi) gave flag=0\n"); failures++; break; }
    }

    // ---- override_community_rate sentinel semantics ----
    // An explicit 0 override must be honored (no infectors ever drawn),
    // an omitted/default (-1) override must fall back to the category's
    // default community_rate, and a positive override must be honored too.
    {
        // Explicit zero: comm_rate = 0 means n_infected is always 0, so
        // infection_probability must return flag=0 and P=0 every time.
        bool saw_infected = false;
        for (int i = 0; i < 200; i++) {
            auto [p, f] = infection_probability(20, par, rng, "Classroom", false, 0.0);
            if (f != 0 || p != 0.0) { saw_infected = true; break; }
        }
        if (saw_infected) {
            printf("FAIL: explicit override_community_rate=0 was not honored\n");
            failures++;
        }
    }
    {
        // Omitted override (uses the default parameter value, -1) must fall
        // back to par.community_rate. Checked deterministically by forcing
        // that field to each extreme, instead of comparing two independent
        // random_device-seeded RNG streams.
        auto par_zero = par;
        par_zero.community_rate = 0.0;
        for (int i = 0; i < 50; i++) {
            auto [p, f] = infection_probability(20, par_zero, rng, "Classroom");
            if (f != 0 || p != 0.0) {
                printf("FAIL: omitted override did not use par.community_rate=0\n");
                failures++; break;
            }
        }
        auto par_one = par;
        par_one.community_rate = 1.0;
        for (int i = 0; i < 50; i++) {
            auto [p, f] = infection_probability(20, par_one, rng, "Classroom");
            if (f != 1 || !(p > 0.0) || !std::isfinite(p)) {
                printf("FAIL: omitted override did not use par.community_rate=1\n");
                failures++; break;
            }
        }
    }
    {
        // Positive override: with require_infectors=true and comm_rate=1.0,
        // every draw must yield n_infected == par.I0 infectors, so flag=1
        // is guaranteed and P should be > 0 (Classroom I0=30 infectors).
        auto [p, f] = infection_probability(20, par, rng, "Classroom", true, 1.0);
        if (f != 1 || !(p > 0.0) || !std::isfinite(p)) {
            printf("FAIL: positive override_community_rate=1.0 not honored (P=%f, flag=%d)\n", p, f);
            failures++;
        }
    }

    // ---- require_infectors with a zero community rate is invalid ----
    // Zero infectors are impossible, so all three model functions must throw
    // std::invalid_argument rather than retry forever.
    {
        bool threw = false;
        try {
            infection_probability(20, par, rng, "Classroom", true, 0.0);
        } catch (const std::invalid_argument&) { threw = true; }
        if (!threw) {
            printf("FAIL: infection_probability accepted require_infectors with rate 0\n");
            failures++;
        }
    }
    {
        bool threw = false;
        try {
            infection_probability_with_inputs(20, par, rng, "Classroom", true, 0.0);
        } catch (const std::invalid_argument&) { threw = true; }
        if (!threw) {
            printf("FAIL: infection_probability_with_inputs accepted require_infectors with rate 0\n");
            failures++;
        }
    }
    {
        bool threw = false;
        try {
            compute_ECAi(par, 0.001, rng, "Classroom", true, 0.0);
        } catch (const std::invalid_argument&) { threw = true; }
        if (!threw) {
            printf("FAIL: compute_ECAi accepted require_infectors with rate 0\n");
            failures++;
        }
    }
    {
        // Same rejection when the zero rate comes from par.community_rate
        // (omitted override) rather than an explicit override.
        auto par_zero_req = par;
        par_zero_req.community_rate = 0.0;
        bool threw = false;
        try {
            compute_ECAi(par_zero_req, 0.001, rng, "Classroom", true);
        } catch (const std::invalid_argument&) { threw = true; }
        if (!threw) {
            printf("FAIL: default zero community rate accepted with require_infectors\n");
            failures++;
        }
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