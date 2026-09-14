// Stage 10 tests for simulation.h/simulation.cpp: one replicate end to end
// (Stage 6 plan -> Stage 7 banks -> stable Stage 7 assignments -> Stage 8
// QER -> Stage 9 ECAi -> exactly N index-ordered outputs -> the frozen
// percentile of docs/MODEL_SPECIFICATION.md D8 / section 12), per
// IMPLEMENTATION_PLAN.md Section 7, Stage 10.
//
// Every expected value in the injected fixture below is derived in this file
// from the frozen equations (independently re-written here, never by calling
// the production QER/ECAi helpers), so the fixture reconciles the
// orchestration against the specification rather than against itself.

#include "simulation.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <optional>
#include <utility>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef ECAI_CPP2_SIMULATION_HEADER_PATH
#error "ECAI_CPP2_SIMULATION_HEADER_PATH must be defined by the build system"
#endif
#ifndef ECAI_CPP2_SIMULATION_SOURCE_PATH
#error "ECAI_CPP2_SIMULATION_SOURCE_PATH must be defined by the build system"
#endif

namespace {

int g_checks = 0;
int g_failures = 0;

void check(bool condition, const std::string& description) {
    ++g_checks;
    if (!condition) {
        std::cerr << "FAIL: " << description << "\n";
        ++g_failures;
    } else {
        std::cout << "PASS: " << description << "\n";
    }
}

bool approx(double actual, double expected, double abs_tol = 1e-12, double rel_tol = 1e-10) {
    if (!std::isfinite(actual) || !std::isfinite(expected)) return false;
    const double diff = std::fabs(actual - expected);
    if (diff <= abs_tol) return true;
    return diff <= rel_tol * std::fabs(expected);
}

// True only for a positive zero (sign bit clear); a zero-infector or clamped
// row must never carry -0.0 into the percentile stage.
bool is_positive_zero(double v) { return v == 0.0 && !std::signbit(v); }

// Bit-for-bit comparison (not ==, so a -0.0/+0.0 difference is caught).
bool bitwise_equal(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

bool bitwise_equal(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!bitwise_equal(a[i], b[i])) return false;
    }
    return true;
}

template <typename Callable>
bool throws_invalid_argument(Callable&& fn) {
    try {
        fn();
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

template <typename Callable>
bool throws_runtime_error(Callable&& fn) {
    try {
        fn();
    } catch (const std::runtime_error&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

using ecai_cpp2::OccupancyCategory;
using ecai_cpp2::lookup_occupancy_category;
using ecai_cpp2::sampling_plan::BankId;
using ecai_cpp2::sampling_plan::kAllBankIds;
using ecai_cpp2::sampling_plan::ParameterBank;
using ecai_cpp2::sampling_plan::ParameterBanks;
using ecai_cpp2::sampling_plan::SamplingPlan;
using ecai_cpp2::seed_manager::EffectiveConfig;
using ecai_cpp2::simulation::EffectiveOverrides;
using ecai_cpp2::simulation::EffectiveScenario;
using ecai_cpp2::simulation::EvaluationOrder;
using ecai_cpp2::simulation::InfectorPolicy;
using ecai_cpp2::simulation::OutputMode;
using ecai_cpp2::simulation::ReplicateOptions;
using ecai_cpp2::simulation::ReplicateResult;

constexpr double kPi = 3.14159265358979323846;

// ------------------------------------------------------------------
// Independent (test-owned) restatements of the frozen equations.
// ------------------------------------------------------------------

// Eq. 5 + Eq. 4 + Eq. 3 written out directly from
// docs/MODEL_SPECIFICATION.md section 3, with Table 3's printed x1e4 C_drop
// and x1e-6 d_bar scalings and Table 1's x1e6 GVL footnote applied here, in
// this file, rather than by calling model::compute_qer.
double expected_qer(double pbr, double c_drop, double d_bar, double e, double gvl_ml, double vf,
                    double rtd, double dk) {
    const double d_m = d_bar * 1e-6;
    const double hydrated = d_m * e;
    const double vdrop = (kPi / 6.0) * hydrated * hydrated * hydrated * (c_drop * 1e4);
    const double ver = pbr * vdrop * (gvl_ml * 1e6) * vf * 1.0;  // R = 1.0
    return rtd * ver / dk;
}

// The frozen ECAi inversion (section 3), clamped at exactly zero.
double expected_ecai(double pbr_s, double gamma, double lambda, double volume, long long i0,
                     double mu, double duration_h, double target_p, double qer_sum) {
    const double q_target = -std::log(1.0 - target_p);
    const double mask_factor = (1.0 - mu) * (1.0 - mu);
    const double term1 = pbr_s * duration_h * mask_factor * qer_sum / q_target;
    const double term2 = volume * (gamma + lambda);
    const double raw = (term1 - term2) / (3.6 * static_cast<double>(i0));
    return raw > 0.0 ? raw : 0.0;
}

// The forward probability evaluation (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md
// section 3.2), written out directly rather than via model::compute_probability.
double expected_probability(double pbr_s, double gamma, double lambda, double volume, long long i0,
                            double mu, double duration_h, double input_ecai, double qer_sum) {
    const double eca_total = 3.6 * static_cast<double>(i0) * input_ecai;
    const double removal = eca_total + volume * (gamma + lambda);
    const double mask_factor = (1.0 - mu) * (1.0 - mu);
    const double dose_numerator = pbr_s * duration_h * mask_factor * qer_sum;
    const double q = dose_numerator / removal;
    const double p_raw = -std::expm1(-q);
    return p_raw > 0.0 ? p_raw : 0.0;
}

// numpy.percentile "linear" (Type 7), written out independently here.
double expected_percentile(std::vector<double> values, double p) {
    std::sort(values.begin(), values.end());
    const double rank = p * static_cast<double>(values.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(rank));
    const std::size_t hi = static_cast<std::size_t>(std::ceil(rank));
    const double frac = rank - static_cast<double>(lo);
    return values[lo] * (1.0 - frac) + values[hi] * frac;
}

// ------------------------------------------------------------------
// The injected deterministic fixture: an explicit ninf vector and
// hand-chosen bank contents, so every expected row is algebraic.
// ------------------------------------------------------------------

// counts = {0, 1, 2, 0, 1}: N = 5, N_zero = 2, N_active = 3, K_total = 4.
SamplingPlan fixture_plan() {
    SamplingPlan plan;
    plan.counts = {0, 1, 2, 0, 1};
    plan.n = 5;
    plan.n_zero = 2;
    plan.n_active = 3;
    plan.k_total = 4;
    plan.susceptible_draws = 3;
    plan.infector_draws = 4;
    return plan;
}

// Susceptible bank contents, index-aligned (bank index 0 -> simulation 1,
// 1 -> simulation 2, 2 -> simulation 4, i.e. active simulations in order).
const std::vector<double> kSusceptiblePbr = {0.5, 0.6, 0.7};
const std::vector<double> kGamma = {0.42, 0.50, 0.61};
const std::vector<double> kLambda = {0.5, 0.4, 0.3};

// Infector bank contents (index 0 -> simulation 1's only infector,
// 1 and 2 -> simulation 2's two infectors, 3 -> simulation 4's infector).
const std::vector<double> kInfectorPbr = {1.0, 2.0, 3.0, 4.0};
const std::vector<double> kCDrop = {10.0, 20.0, 30.0, 40.0};
const std::vector<double> kDBar = {1.0, 2.0, 1.0, 2.0};
const std::vector<double> kE = {2.0, 3.0, 4.0, 5.0};
// Infector 0's GVL is deliberately three orders of magnitude below the rest
// so that simulation 1's raw ECAi is negative (exercising the clamp) while
// simulations 2 and 4 stay strictly positive.
const std::vector<double> kGvl = {1.0e6, 2.0e9, 1.0e10, 5.0e9};
const std::vector<double> kVf = {1.0e-4, 1.0e-3, 1.0e-2, 5.0e-3};
const std::vector<double> kRtd = {0.43, 0.50, 0.60, 0.65};
const std::vector<double> kDk = {5.0, 10.0, 15.0, 8.0};

const std::vector<double>& fixture_bank_values(BankId id) {
    switch (id) {
        case BankId::susceptible_pbr: return kSusceptiblePbr;
        case BankId::gamma: return kGamma;
        case BankId::lambda: return kLambda;
        case BankId::infector_pbr: return kInfectorPbr;
        case BankId::c_drop: return kCDrop;
        case BankId::d_bar: return kDBar;
        case BankId::E: return kE;
        case BankId::GVL: return kGvl;
        case BankId::VF: return kVf;
        case BankId::RTD: return kRtd;
        case BankId::DK: return kDk;
    }
    throw std::invalid_argument("fixture_bank_values: unknown bank");
}

ParameterBanks fixture_banks() {
    ParameterBanks banks;
    banks.n = 5;
    banks.n_active = 3;
    banks.k_total = 4;
    banks.budget = ecai_cpp2::sampling_plan::compute_bank_budget(5, 3, 4);
    for (BankId id : kAllBankIds) {
        ParameterBank bank;
        bank.id = id;
        bank.stream = ecai_cpp2::sampling_plan::bank_stream_id(id);
        bank.seed = 0;
        bank.generated = true;
        bank.values = fixture_bank_values(id);
        banks.banks[static_cast<std::size_t>(id)] = bank;
    }
    return banks;
}

// A deterministic all-active, one-infector-per-row plan/bank pair of size n,
// used only to pin the interleaved evaluation schedule and to compare it
// against the forward schedule's output. Values are hand-chosen distinct
// positive numbers in the same order of magnitude as the injected fixture
// above so every row's ECAi stays finite.
SamplingPlan all_active_plan(std::size_t n) {
    SamplingPlan plan;
    plan.counts.assign(n, 1);
    plan.n = n;
    plan.n_zero = 0;
    plan.n_active = n;
    plan.k_total = n;
    plan.susceptible_draws = n;
    plan.infector_draws = n;
    return plan;
}

double all_active_bank_base(BankId id, std::size_t i) {
    const double fi = static_cast<double>(i);
    switch (id) {
        case BankId::susceptible_pbr: return 0.5 + 0.1 * fi;
        case BankId::gamma: return 0.42 + 0.04 * fi;
        case BankId::lambda: return 0.3 + 0.05 * fi;
        case BankId::infector_pbr: return 1.0 + 0.5 * fi;
        case BankId::c_drop: return 10.0 + 5.0 * fi;
        case BankId::d_bar: return 1.0 + 0.2 * fi;
        case BankId::E: return 2.0 + 0.3 * fi;
        case BankId::GVL: return 1.0e9 + 1.0e8 * fi;
        case BankId::VF: return 1.0e-3 + 1.0e-4 * fi;
        case BankId::RTD: return 0.5 + 0.02 * fi;
        case BankId::DK: return 5.0 + 1.0 * fi;
    }
    throw std::invalid_argument("all_active_bank_base: unknown bank");
}

ParameterBanks all_active_banks(std::size_t n) {
    ParameterBanks banks;
    banks.n = n;
    banks.n_active = n;
    banks.k_total = n;
    banks.budget = ecai_cpp2::sampling_plan::compute_bank_budget(n, n, n);
    for (BankId id : kAllBankIds) {
        ParameterBank bank;
        bank.id = id;
        bank.stream = ecai_cpp2::sampling_plan::bank_stream_id(id);
        bank.seed = 0;
        bank.generated = true;
        bank.values.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            bank.values[i] = all_active_bank_base(id, i);
        }
        banks.banks[static_cast<std::size_t>(id)] = bank;
    }
    return banks;
}

EffectiveScenario fixture_scenario() {
    EffectiveScenario scenario;
    scenario.category_name = "fixture";
    scenario.occupants = 10;
    scenario.volume_m3 = 100.0;
    scenario.mask_efficiency = 0.0;
    scenario.community_infection_rate = 0.01;
    scenario.target_probability = 0.001;
    scenario.duration_hours = 1.0;
    scenario.percentile = 0.96;
    scenario.simulations_n = 5;
    return scenario;
}

// The three active rows' expected ECAi, computed from the fixture contents
// through the test-owned equations above.
std::vector<double> fixture_expected_outputs() {
    const double qer0 = expected_qer(kInfectorPbr[0], kCDrop[0], kDBar[0], kE[0], kGvl[0], kVf[0],
                                     kRtd[0], kDk[0]);
    const double qer1 = expected_qer(kInfectorPbr[1], kCDrop[1], kDBar[1], kE[1], kGvl[1], kVf[1],
                                     kRtd[1], kDk[1]);
    const double qer2 = expected_qer(kInfectorPbr[2], kCDrop[2], kDBar[2], kE[2], kGvl[2], kVf[2],
                                     kRtd[2], kDk[2]);
    const double qer3 = expected_qer(kInfectorPbr[3], kCDrop[3], kDBar[3], kE[3], kGvl[3], kVf[3],
                                     kRtd[3], kDk[3]);
    const EffectiveScenario s = fixture_scenario();
    std::vector<double> expected(5, 0.0);
    expected[1] = expected_ecai(kSusceptiblePbr[0], kGamma[0], kLambda[0], s.volume_m3,
                                s.occupants, s.mask_efficiency, s.duration_hours,
                                s.target_probability, qer0);
    expected[2] = expected_ecai(kSusceptiblePbr[1], kGamma[1], kLambda[1], s.volume_m3,
                                s.occupants, s.mask_efficiency, s.duration_hours,
                                s.target_probability, qer1 + qer2);
    expected[4] = expected_ecai(kSusceptiblePbr[2], kGamma[2], kLambda[2], s.volume_m3,
                                s.occupants, s.mask_efficiency, s.duration_hours,
                                s.target_probability, qer3);
    return expected;
}

// Same fixture, resolved into probability mode at a moderate fixed input
// ECAi chosen so that none of the three active rows saturates.
EffectiveScenario fixture_probability_scenario() {
    EffectiveScenario scenario = fixture_scenario();
    scenario.mode = OutputMode::probability;
    scenario.input_ecai_l_s_person = 5.0;
    return scenario;
}

// The three active rows' expected forward probability at fixture_probability_
// scenario()'s input ECAi, computed from the fixture contents through the
// test-owned equations above (same QER values as fixture_expected_outputs).
std::vector<double> fixture_expected_probabilities(double input_ecai) {
    const double qer0 = expected_qer(kInfectorPbr[0], kCDrop[0], kDBar[0], kE[0], kGvl[0], kVf[0],
                                     kRtd[0], kDk[0]);
    const double qer1 = expected_qer(kInfectorPbr[1], kCDrop[1], kDBar[1], kE[1], kGvl[1], kVf[1],
                                     kRtd[1], kDk[1]);
    const double qer2 = expected_qer(kInfectorPbr[2], kCDrop[2], kDBar[2], kE[2], kGvl[2], kVf[2],
                                     kRtd[2], kDk[2]);
    const double qer3 = expected_qer(kInfectorPbr[3], kCDrop[3], kDBar[3], kE[3], kGvl[3], kVf[3],
                                     kRtd[3], kDk[3]);
    const EffectiveScenario s = fixture_scenario();
    std::vector<double> expected(5, 0.0);
    expected[1] = expected_probability(kSusceptiblePbr[0], kGamma[0], kLambda[0], s.volume_m3,
                                       s.occupants, s.mask_efficiency, s.duration_hours,
                                       input_ecai, qer0);
    expected[2] = expected_probability(kSusceptiblePbr[1], kGamma[1], kLambda[1], s.volume_m3,
                                       s.occupants, s.mask_efficiency, s.duration_hours,
                                       input_ecai, qer1 + qer2);
    expected[4] = expected_probability(kSusceptiblePbr[2], kGamma[2], kLambda[2], s.volume_m3,
                                       s.occupants, s.mask_efficiency, s.duration_hours,
                                       input_ecai, qer3);
    return expected;
}

// ------------------------------------------------------------------
// Section 1: the injected deterministic replicate.
// ------------------------------------------------------------------

void test_injected_fixture() {
    ReplicateOptions options;
    options.collect_audit = true;
    const ReplicateResult result = ecai_cpp2::simulation::run_replicate_from_banks(
        fixture_plan(), fixture_banks(), fixture_scenario(), options);

    check(result.ecai_l_s_person.size() == 5, "injected fixture emits exactly N = 5 ECAi values");
    check(result.infector_counts == fixture_plan().counts,
          "injected fixture preserves the authoritative Stage 6 count vector exactly");

    const std::vector<double> expected = fixture_expected_outputs();
    for (std::size_t i = 0; i < expected.size(); ++i) {
        check(approx(result.ecai_l_s_person[i], expected[i]),
              "injected fixture row " + std::to_string(i) +
                  " equals its independently computed ECAi");
    }

    check(is_positive_zero(result.ecai_l_s_person[0]) &&
              is_positive_zero(result.ecai_l_s_person[3]),
          "zero-infector rows 0 and 3 are exactly +0.0");
    check(is_positive_zero(result.ecai_l_s_person[1]),
          "row 1's negative raw ECAi is clamped to exactly +0.0");
    check(result.ecai_l_s_person[2] > 0.0 && result.ecai_l_s_person[4] > 0.0,
          "rows 2 and 4 are strictly positive");

    check(result.counters.n == 5, "counter n == 5");
    check(result.counters.zero_infector_count == 2, "zero_infector_count == 2");
    check(result.counters.negative_clamp_count == 1, "negative_clamp_count == 1");
    check(result.counters.output_zero_count == 3,
          "output_zero_count == zero_infector_count + negative_clamp_count");
    check(result.counters.active_count == 3, "active_count == 3");
    check(result.counters.k_total == 4, "k_total == 4");
    check(result.counters.qer_evaluations == 4, "one QER evaluation per realized infector");
    check(result.counters.ecai_evaluations == 3, "one ECAi evaluation per active row");
    check(result.counters.susceptible_values_consumed == 9, "3 susceptible banks x 3 active rows");
    check(result.counters.infector_values_consumed == 32, "8 infector banks x 4 infectors");

    check(approx(result.percentile_l_s_person, expected_percentile(expected, 0.96)),
          "injected fixture percentile matches the independently interpolated value");

    // Cursor exhaustion and planned-length reconciliation.
    bool all_exhausted = result.cursors_exhausted;
    for (BankId id : kAllBankIds) {
        const std::size_t ordinal = static_cast<std::size_t>(id);
        const std::size_t planned =
            ecai_cpp2::sampling_plan::is_infector_bank(id) ? std::size_t{4} : std::size_t{3};
        if (result.planned_bank_lengths[ordinal] != planned) all_exhausted = false;
        if (result.final_cursors[ordinal] != planned) all_exhausted = false;
    }
    check(all_exhausted, "every bank cursor is exactly exhausted at its planned length");

    // Exact consumption reconciliation: every planned value consumed once.
    check(result.audit_collected, "audit was collected when requested");
    const std::set<std::uint64_t> consumed(result.audit.consumption_ids.begin(),
                                           result.audit.consumption_ids.end());
    check(result.audit.consumption_ids.size() == 41,
          "audit records 3*3 + 8*4 = 41 consumed values");
    check(consumed.size() == result.audit.consumption_ids.size(),
          "no consumption id appears twice");
    std::set<std::uint64_t> planned_ids;
    for (BankId id : kAllBankIds) {
        const std::size_t planned =
            ecai_cpp2::sampling_plan::is_infector_bank(id) ? std::size_t{4} : std::size_t{3};
        for (std::size_t index = 0; index < planned; ++index) {
            planned_ids.insert(ecai_cpp2::sampling_plan::consumption_id(id, index));
        }
    }
    check(consumed == planned_ids, "the consumed id set equals the planned id set exactly");
}

// ------------------------------------------------------------------
// Section 2: the seeded end-to-end replicate.
// ------------------------------------------------------------------

// Builds the effective configuration the resolver must agree with: category
// defaults, with any supplied override applied by this test independently of
// the production resolver.
EffectiveConfig make_config(const OccupancyCategory& category, std::uint64_t n,
                            const EffectiveOverrides& overrides = EffectiveOverrides{}) {
    EffectiveConfig config;
    config.category_name = category.canonical_name;
    config.i0 = overrides.occupants ? *overrides.occupants : category.i0;
    config.volume_m3 = overrides.volume_m3 ? *overrides.volume_m3 : category.volume_m3;
    config.mask_efficiency =
        overrides.mask_efficiency ? *overrides.mask_efficiency : category.mask_efficiency;
    config.community_infection_rate = overrides.community_infection_rate
                                          ? *overrides.community_infection_rate
                                          : category.community_infection_rate;
    config.target_probability = 0.001;
    config.duration_hours = 1.0;
    config.percentile = 0.96;
    config.confidence_level = 0.95;
    config.simulations_n = n;
    config.replicates_r = 1;
    return config;
}

void test_scenario_resolution() {
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");

    const EffectiveScenario defaults =
        ecai_cpp2::simulation::resolve_scenario(patient, make_config(patient, 16));
    check(defaults.category_name == patient.canonical_name,
          "resolved scenario keeps the canonical category name");
    check(defaults.occupants == patient.i0 && defaults.volume_m3 == patient.volume_m3,
          "category defaults supply I0 and VOL when nothing is overridden");
    check(defaults.mask_efficiency == patient.mask_efficiency &&
              defaults.community_infection_rate == patient.community_infection_rate,
          "category defaults supply mask efficiency and community rate");
    check(!defaults.custom_scenario, "an un-overridden run is not a custom scenario");

    EffectiveOverrides overrides;
    overrides.occupants = 7;
    overrides.volume_m3 = 250.0;
    overrides.mask_efficiency = 0.0;
    overrides.community_infection_rate = 0.05;
    const EffectiveScenario overridden = ecai_cpp2::simulation::resolve_scenario(
        patient, make_config(patient, 16, overrides), overrides);
    check(overridden.occupants == 7 && overridden.volume_m3 == 250.0 &&
              overridden.mask_efficiency == 0.0 && overridden.community_infection_rate == 0.05,
          "explicit overrides replace every category default they name");
    check(overridden.custom_scenario, "an overridden run is flagged as a custom scenario");
    check(overridden.occupants_overridden && overridden.volume_overridden &&
              overridden.mask_efficiency_overridden &&
              overridden.community_infection_rate_overridden,
          "each override is individually flagged");

    // The effective configuration (which fingerprints the seeds) must agree
    // with the resolved category+override values.
    EffectiveConfig mismatched = make_config(patient, 16);
    mismatched.volume_m3 = patient.volume_m3 + 1.0;
    check(throws_invalid_argument(
              [&] { ecai_cpp2::simulation::resolve_scenario(patient, mismatched); }),
          "a configuration whose VOL disagrees with the category is rejected");
    EffectiveConfig wrong_name = make_config(patient, 16);
    wrong_name.category_name = "Classroom";
    check(throws_invalid_argument(
              [&] { ecai_cpp2::simulation::resolve_scenario(patient, wrong_name); }),
          "a configuration naming a different category is rejected");
    EffectiveConfig missing_override = make_config(patient, 16);
    check(throws_invalid_argument([&] {
              ecai_cpp2::simulation::resolve_scenario(patient, missing_override, overrides);
          }),
          "overrides not reflected in the effective configuration are rejected");
}

void test_seeded_replicate_determinism() {
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");
    const EffectiveConfig config = make_config(patient, 64);

    const ReplicateResult a = ecai_cpp2::simulation::run_replicate(485, config, 0, patient);
    const ReplicateResult b = ecai_cpp2::simulation::run_replicate(485, config, 0, patient);
    const ReplicateResult indexed = ecai_cpp2::simulation::run_replicate(485, config, 7, patient);
    check(a.replicate_index == 0 && indexed.replicate_index == 7,
          "top-level replicate results preserve the exact requested replicate identity");
    check(a.infector_counts.size() == 64,
          "top-level replicate preserves exactly one authoritative infector count per row");
    const ecai_cpp2::seed_manager::SeedManager top_level_seeds(485);
    const SamplingPlan top_level_plan = ecai_cpp2::sampling_plan::build_sampling_plan(
        top_level_seeds, config, 0, 64, patient.i0, patient.community_infection_rate,
        ecai_cpp2::sampling_plan::ResourcePolicy{});
    check(a.infector_counts == top_level_plan.counts,
          "top-level replicate count vector exactly equals the authoritative Stage 6 plan");
    bool counts_ground_outputs = true;
    for (std::size_t row = 0; row < a.infector_counts.size(); ++row) {
        if (a.infector_counts[row] < 0 || a.infector_counts[row] > patient.i0 ||
            (a.infector_counts[row] == 0 && !is_positive_zero(a.ecai_l_s_person[row]))) {
            counts_ground_outputs = false;
        }
    }
    check(counts_ground_outputs,
          "top-level Stage 6 counts are bounded and every zero-infector row emits exact +0.0");
    check(a.ecai_l_s_person.size() == 64, "seeded replicate emits exactly N = 64 values");
    check(bitwise_equal(a.ecai_l_s_person, b.ecai_l_s_person),
          "the same base seed, configuration and replicate index reproduce bit-identical output");
    check(bitwise_equal(a.percentile_l_s_person, b.percentile_l_s_person),
          "the reproduced percentile is bit-identical");

    const ReplicateResult other_seed = ecai_cpp2::simulation::run_replicate(486, config, 0, patient);
    check(!bitwise_equal(a.ecai_l_s_person, other_seed.ecai_l_s_person),
          "a changed base seed changes the replicate output");

    const ReplicateResult other_replicate =
        ecai_cpp2::simulation::run_replicate(485, config, 1, patient);
    check(!bitwise_equal(a.ecai_l_s_person, other_replicate.ecai_l_s_person),
          "a changed replicate index changes the replicate output");

    // No state leakage: interleaving unrelated runs cannot perturb a repeat.
    const ReplicateResult repeat = ecai_cpp2::simulation::run_replicate(485, config, 0, patient);
    check(bitwise_equal(a.ecai_l_s_person, repeat.ecai_l_s_person),
          "no state leaks between replicates: a repeat after other runs is unchanged");

    // Category defaults vs. an explicit override of the same field.
    EffectiveOverrides overrides;
    overrides.community_infection_rate = 0.5;
    const ReplicateResult overridden = ecai_cpp2::simulation::run_replicate(
        485, make_config(patient, 64, overrides), 0, patient, overrides);
    check(overridden.counters.zero_infector_count < a.counters.zero_infector_count,
          "raising the community rate above the category default lowers the zero-infector count");
    check(overridden.counters.n == 64 && overridden.ecai_l_s_person.size() == 64,
          "the overridden run still emits exactly N values");

    check(a.counters.zero_infector_count + a.counters.active_count == 64,
          "zero-infector and active rows partition the replicate");
    check(a.counters.qer_evaluations == a.counters.k_total,
          "exactly one QER evaluation per realized infector");
    check(a.counters.ecai_evaluations == a.counters.active_count,
          "exactly one ECAi evaluation per active row");
    check(a.cursors_exhausted, "every cursor is exhausted in the seeded replicate");
}

void test_pinned_seeded_replicate() {
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");
    const EffectiveConfig config = make_config(patient, 8);
    const ReplicateResult result = ecai_cpp2::simulation::run_replicate(485, config, 0, patient);

    // Pinned once from this frozen WSL toolchain (g++ 13, -fno-fast-math
    // -ffp-contract=off) and thereafter treated as a regression fixture.
    const std::vector<double> pinned = {0.0, 0.0, 0.0, 0.0,
                                       0.0, 0.0, 0.0, 17.67698541527113};
    const double pinned_percentile = 12.727429498995209;
    check(result.ecai_l_s_person.size() == pinned.size(), "pinned vector length matches");
    bool matches = true;
    for (std::size_t i = 0; i < pinned.size() && i < result.ecai_l_s_person.size(); ++i) {
        if (!approx(result.ecai_l_s_person[i], pinned[i], 1e-12, 1e-12)) matches = false;
    }
    check(matches, "the pinned N = 8 seeded output vector is reproduced");
    check(approx(result.percentile_l_s_person, pinned_percentile, 1e-12, 1e-12),
          "the pinned seeded percentile is reproduced");
}

void test_all_zero_replicate() {
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");
    EffectiveOverrides overrides;
    overrides.community_infection_rate = 0.0;
    ReplicateOptions options;
    options.collect_audit = true;
    const ReplicateResult result = ecai_cpp2::simulation::run_replicate(
        485, make_config(patient, 32, overrides), 0, patient, overrides, options);

    check(result.counters.zero_infector_count == 32,
          "community rate 0 makes every simulation a zero-infector row");
    check(result.infector_counts == std::vector<long long>(32, 0),
          "an all-zero replicate retains 32 authoritative zero counts in stable order");
    check(result.counters.active_count == 0 && result.counters.k_total == 0,
          "an all-zero replicate has no active row and no realized infector");
    bool all_positive_zero = result.ecai_l_s_person.size() == 32;
    for (const double v : result.ecai_l_s_person) {
        if (!is_positive_zero(v)) all_positive_zero = false;
    }
    check(all_positive_zero, "every output of an all-zero replicate is exactly +0.0");
    check(result.counters.qer_evaluations == 0 && result.counters.ecai_evaluations == 0,
          "an all-zero replicate makes no QER or ECAi model request");
    bool none_generated = true;
    std::size_t total_planned = 0;
    for (BankId id : kAllBankIds) {
        const std::size_t ordinal = static_cast<std::size_t>(id);
        if (result.bank_generated[ordinal]) none_generated = false;
        total_planned += result.planned_bank_lengths[ordinal];
    }
    check(none_generated, "no non-binomial parameter bank is generated for an all-zero replicate");
    check(total_planned == 0, "every parameter bank is planned at length zero");
    check(result.audit.consumption_ids.empty(), "no parameter value is consumed");
    check(is_positive_zero(result.percentile_l_s_person),
          "the percentile of an all-zero replicate is exactly +0.0");
    check(result.cursors_exhausted, "empty cursors count as exhausted");
}

void test_person_identity_and_traversal_order() {
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");
    const EffectiveConfig config = make_config(patient, 48);
    ReplicateOptions forward;
    forward.collect_audit = true;
    forward.evaluation_order = EvaluationOrder::forward;
    ReplicateOptions reverse = forward;
    reverse.evaluation_order = EvaluationOrder::reverse;
    ReplicateOptions interleaved = forward;
    interleaved.evaluation_order = EvaluationOrder::interleaved;

    const ReplicateResult f = ecai_cpp2::simulation::run_replicate(485, config, 2, patient,
                                                                    EffectiveOverrides{}, forward);
    const ReplicateResult r = ecai_cpp2::simulation::run_replicate(485, config, 2, patient,
                                                                    EffectiveOverrides{}, reverse);
    const ReplicateResult i = ecai_cpp2::simulation::run_replicate(
        485, config, 2, patient, EffectiveOverrides{}, interleaved);

    check(f.counters.active_count > 1 && f.counters.k_total > 1,
          "the traversal-order fixture actually contains several active rows");
    check(bitwise_equal(f.ecai_l_s_person, r.ecai_l_s_person),
          "reverse evaluation order produces a bit-identical index-ordered output vector");
    check(bitwise_equal(f.ecai_l_s_person, i.ecai_l_s_person),
          "interleaved evaluation order produces a bit-identical output vector");
    check(bitwise_equal(f.percentile_l_s_person, r.percentile_l_s_person) &&
              bitwise_equal(f.percentile_l_s_person, i.percentile_l_s_person),
          "the percentile is unchanged by evaluation order");

    check(f.audit.rows.size() == 48 && r.audit.rows.size() == 48,
          "the audit records one row per simulation, in simulation index order");
    bool mapping_stable = true;
    bool sums_stable = true;
    bool indices_ordered = true;
    for (std::size_t row = 0; row < f.audit.rows.size(); ++row) {
        const auto& fr = f.audit.rows[row];
        const auto& rr = r.audit.rows[row];
        if (fr.simulation_index != row || rr.simulation_index != row) indices_ordered = false;
        if (fr.susceptible_ids != rr.susceptible_ids) mapping_stable = false;
        if (fr.infector_ids != rr.infector_ids) mapping_stable = false;
        if (!bitwise_equal(fr.qer_sum_per_h, rr.qer_sum_per_h)) sums_stable = false;
    }
    check(indices_ordered, "audit rows stay in stable simulation index order");
    check(mapping_stable, "person-to-value mapping is identical under a different traversal order");
    check(sums_stable, "per-row QER sums are bit-identical under a different traversal order");
    check(f.audit.evaluation_order != r.audit.evaluation_order,
          "the two runs really did evaluate their rows in different orders");

    // Distinct ids per infector, and no PBR shared with the susceptible.
    std::set<std::uint64_t> all_ids;
    bool distinct_within_infector = true;
    bool distinct_across_all = true;
    bool pbr_never_shared = true;
    std::size_t audited_infectors = 0;
    for (const auto& row : f.audit.rows) {
        if (!row.active) {
            check(row.infector_ids.empty(), "a zero-infector audit row consumes no infector value");
            continue;
        }
        for (const std::uint64_t id : row.susceptible_ids) {
            if (!all_ids.insert(id).second) distinct_across_all = false;
        }
        for (const auto& ids : row.infector_ids) {
            ++audited_infectors;
            const std::set<std::uint64_t> unique(ids.begin(), ids.end());
            if (unique.size() != 8) distinct_within_infector = false;
            for (const std::uint64_t id : ids) {
                if (!all_ids.insert(id).second) distinct_across_all = false;
            }
            if (ecai_cpp2::sampling_plan::consumption_id_bank(ids[0]) != BankId::infector_pbr ||
                ids[0] == row.susceptible_ids[0]) {
                pbr_never_shared = false;
            }
        }
    }
    check(audited_infectors == f.counters.k_total, "every realized infector is audited");
    check(distinct_within_infector, "every infector holds eight distinct consumption ids");
    check(distinct_across_all,
          "no consumption id is shared between two people or two parameters");
    check(pbr_never_shared,
          "an infector's PBR id comes from the infector PBR bank and never matches the "
          "susceptible PBR id");
    check(all_ids.size() == 3 * f.counters.active_count + 8 * f.counters.k_total,
          "the audited id count equals 3*N_active + 8*K_total");
}

void test_large_replicate_smoke() {
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");
    const EffectiveConfig config = make_config(patient, 10000);
    const ReplicateResult result = ecai_cpp2::simulation::run_replicate(485, config, 0, patient);

    check(result.ecai_l_s_person.size() == 10000, "N = 10000 smoke emits exactly 10000 values");
    bool all_finite = true;
    bool all_nonnegative = true;
    std::size_t zeros = 0;
    for (const double v : result.ecai_l_s_person) {
        if (!std::isfinite(v)) all_finite = false;
        if (v < 0.0) all_nonnegative = false;
        if (v == 0.0) {
            ++zeros;
            if (std::signbit(v)) all_nonnegative = false;
        }
    }
    check(all_finite, "every N = 10000 output is finite");
    check(all_nonnegative, "no N = 10000 output is negative or a negative zero");
    check(zeros == result.counters.output_zero_count, "output_zero_count matches the vector");
    check(result.counters.output_zero_count ==
              result.counters.zero_infector_count + result.counters.negative_clamp_count,
          "output zeros are exactly the zero-infector rows plus the clamped rows");
    check(result.counters.zero_infector_count > 0 && result.counters.negative_clamp_count > 0,
          "the smoke replicate exercises both zero sources");
    check(std::isfinite(result.percentile_l_s_person) && result.percentile_l_s_person > 0.0,
          "the N = 10000 percentile is finite and positive");
    check(result.cursors_exhausted, "the N = 10000 replicate exhausts every cursor");
}

// ------------------------------------------------------------------
// Section 3: the frozen percentile convention (D8 / section 12).
// ------------------------------------------------------------------

void test_percentile_convention() {
    using ecai_cpp2::simulation::percentile_linear;

    // N = 1: the single value is returned for every accepted p.
    const std::vector<double> single = {4.25};
    check(percentile_linear(single, 0.0) == 4.25 && percentile_linear(single, 0.5) == 4.25 &&
              percentile_linear(single, 0.96) == 4.25 && percentile_linear(single, 1.0) == 4.25,
          "N = 1 returns the single value for p = 0, 0.5, 0.96 and 1");

    // Odd N, fractional rank: rank = 0.96 * 4 = 3.84.
    const std::vector<double> odd = {1.0, 2.0, 3.0, 4.0, 5.0};
    check(approx(percentile_linear(odd, 0.96), 4.84, 1e-9),
          "odd N fractional rank interpolates 3.84 between 4 and 5");
    // Odd N, integer rank: rank = 0.5 * 4 = 2, i.e. the exact median element.
    check(percentile_linear(odd, 0.5) == 3.0, "odd N integer rank returns the exact element");

    // Even N: rank = 0.5 * 3 = 1.5, midway between 2 and 3.
    const std::vector<double> even = {1.0, 2.0, 3.0, 4.0};
    check(percentile_linear(even, 0.5) == 2.5, "even N median interpolates the two middle values");
    check(approx(percentile_linear(even, 0.96), 3.88, 1e-9),
          "even N fractional rank matches the independent Type-7 formula");
    check(approx(percentile_linear(even, 0.96), expected_percentile(even, 0.96), 1e-9),
          "the production percentile agrees with the test-owned Type-7 restatement");

    // p = 0 and p = 1 are the accepted endpoints of the frozen domain.
    check(percentile_linear(odd, 0.0) == 1.0, "p = 0 returns the minimum");
    check(percentile_linear(odd, 1.0) == 5.0, "p = 1 returns the maximum");

    // Duplicates and zeros, including an all-zero list.
    const std::vector<double> duplicates = {0.0, 0.0, 0.0, 5.0, 5.0};
    check(percentile_linear(duplicates, 0.5) == 0.0,
          "a duplicated value at the rank is returned unchanged");
    check(approx(percentile_linear(duplicates, 0.96), 5.0, 1e-9),
          "a rank inside a run of duplicates interpolates to that duplicate");
    const std::vector<double> zeros(7, 0.0);
    check(is_positive_zero(percentile_linear(zeros, 0.96)),
          "an all-zero list yields exactly +0.0");

    // Unsorted input is sorted internally, into a copy: the caller's
    // index-ordered vector is never mutated.
    std::vector<double> unsorted = {5.0, 1.0, 4.0, 2.0, 3.0};
    const std::vector<double> before = unsorted;
    const double from_unsorted = percentile_linear(unsorted, 0.96);
    check(bitwise_equal(unsorted, before),
          "percentile_linear does not mutate or reorder the caller's vector");
    check(from_unsorted == percentile_linear(odd, 0.96),
          "an unsorted list gives the same percentile as its sorted permutation");

    // Rejected inputs.
    check(throws_invalid_argument([&] { percentile_linear(std::vector<double>{}, 0.96); }),
          "an empty value list is rejected");
    const std::vector<double> with_nan = {1.0, std::numeric_limits<double>::quiet_NaN(), 3.0};
    check(throws_invalid_argument([&] { percentile_linear(with_nan, 0.96); }),
          "a NaN value is rejected");
    const std::vector<double> with_inf = {1.0, std::numeric_limits<double>::infinity(), 3.0};
    check(throws_invalid_argument([&] { percentile_linear(with_inf, 0.96); }),
          "an infinite value is rejected");
    check(throws_invalid_argument(
              [&] { percentile_linear(odd, std::numeric_limits<double>::quiet_NaN()); }),
          "a NaN percentile is rejected");
    check(throws_invalid_argument([&] { percentile_linear(odd, -0.01); }),
          "a percentile below 0 is rejected");
    check(throws_invalid_argument([&] { percentile_linear(odd, 1.01); }),
          "a percentile above 1 is rejected");
}

// ------------------------------------------------------------------
// Section 4: validation, resource limits, and bank non-exposure.
// ------------------------------------------------------------------

void test_rejects_malformed_inputs() {
    using ecai_cpp2::simulation::run_replicate_from_banks;

    {
        SamplingPlan plan = fixture_plan();
        plan.n = 6;  // disagrees with counts.size()
        check(throws_invalid_argument(
                  [&] { run_replicate_from_banks(plan, fixture_banks(), fixture_scenario()); }),
              "a plan whose N disagrees with its counts vector is rejected before any work");
    }
    {
        SamplingPlan plan = fixture_plan();
        plan.n_active = 2;  // disagrees with a fresh summary of counts
        check(throws_invalid_argument(
                  [&] { run_replicate_from_banks(plan, fixture_banks(), fixture_scenario()); }),
              "a plan whose accounting disagrees with its own counts is rejected");
    }
    {
        SamplingPlan plan = fixture_plan();
        plan.counts[0] = -1;
        check(throws_invalid_argument(
                  [&] { run_replicate_from_banks(plan, fixture_banks(), fixture_scenario()); }),
              "a negative infector count is rejected");
    }
    {
        ParameterBanks banks = fixture_banks();
        banks.banks[static_cast<std::size_t>(BankId::GVL)].values.pop_back();
        check(throws_invalid_argument(
                  [&] { run_replicate_from_banks(fixture_plan(), banks, fixture_scenario()); }),
              "a bank shorter than its planned length is rejected");
    }
    {
        ParameterBanks banks = fixture_banks();
        banks.k_total = 5;
        check(throws_invalid_argument(
                  [&] { run_replicate_from_banks(fixture_plan(), banks, fixture_scenario()); }),
              "banks sized for a different plan are rejected");
    }
    {
        EffectiveScenario scenario = fixture_scenario();
        scenario.simulations_n = 4;
        check(throws_invalid_argument(
                  [&] { run_replicate_from_banks(fixture_plan(), fixture_banks(), scenario); }),
              "a scenario whose simulation count disagrees with the plan is rejected");
    }
    const std::vector<std::pair<std::string, EffectiveScenario>> bad_scenarios = [] {
        std::vector<std::pair<std::string, EffectiveScenario>> cases;
        EffectiveScenario s = fixture_scenario();
        s.occupants = 0;
        cases.emplace_back("I0 = 0", s);
        s = fixture_scenario();
        s.volume_m3 = 0.0;
        cases.emplace_back("VOL = 0", s);
        s = fixture_scenario();
        s.mask_efficiency = 1.5;
        cases.emplace_back("mask efficiency above 1", s);
        s = fixture_scenario();
        s.target_probability = 0.0;
        cases.emplace_back("target probability 0", s);
        s = fixture_scenario();
        s.duration_hours = -1.0;
        cases.emplace_back("negative duration", s);
        s = fixture_scenario();
        s.percentile = 1.5;
        cases.emplace_back("percentile above 1", s);
        s = fixture_scenario();
        s.volume_m3 = std::numeric_limits<double>::quiet_NaN();
        cases.emplace_back("nonfinite VOL", s);
        return cases;
    }();
    for (const auto& entry : bad_scenarios) {
        check(throws_invalid_argument([&] {
                  run_replicate_from_banks(fixture_plan(), fixture_banks(), entry.second);
              }),
              "an out-of-domain scenario (" + entry.first + ") is rejected");
    }
}

void test_resource_limits_and_overflow() {
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");
    const EffectiveConfig config = make_config(patient, 64);

    ReplicateOptions tight_simulations;
    tight_simulations.policy.max_simulations = 8;
    check(throws_runtime_error([&] {
              ecai_cpp2::simulation::run_replicate(485, config, 0, patient, EffectiveOverrides{},
                                                    tight_simulations);
          }),
          "a run exceeding the simulation-count ceiling is rejected before allocation");

    ReplicateOptions tight_infectors;
    tight_infectors.policy.max_total_infectors = 1;
    check(throws_runtime_error([&] {
              ecai_cpp2::simulation::run_replicate(485, config, 0, patient, EffectiveOverrides{},
                                                    tight_infectors);
          }),
          "a run exceeding the total-infector ceiling is rejected before allocation");

    ReplicateOptions tight_bytes;
    tight_bytes.policy.max_bytes = 64;
    check(throws_runtime_error([&] {
              ecai_cpp2::simulation::run_replicate(485, config, 0, patient, EffectiveOverrides{},
                                                    tight_bytes);
          }),
          "a run exceeding the memory budget is rejected before allocation");

    EffectiveConfig zero_n = make_config(patient, 0);
    check(throws_invalid_argument(
              [&] { ecai_cpp2::simulation::run_replicate(485, zero_n, 0, patient); }),
          "a configuration with N = 0 is rejected");

    // A huge I0 override makes the checked N * I0 bound overflow the frozen
    // default ceilings, and must be rejected rather than attempted.
    EffectiveOverrides huge;
    huge.occupants = 1000000000LL;
    check(throws_runtime_error([&] {
              ecai_cpp2::simulation::run_replicate(485, make_config(patient, 64, huge), 0, patient,
                                                    huge);
          }),
          "an I0 override that overflows the K_total ceiling is rejected");
}

std::string read_file(const char* path) {
    std::ifstream file(path);
    if (!file) {
        return std::string();
    }
    std::ostringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

void test_banks_are_not_exposed_after_traversal() {
    const std::string header = read_file(ECAI_CPP2_SIMULATION_HEADER_PATH);
    const std::string source = read_file(ECAI_CPP2_SIMULATION_SOURCE_PATH);
    check(!header.empty() && !source.empty(), "the simulation header and source are readable");

    const std::size_t start = header.find("struct ReplicateResult {");
    const std::size_t end = header.find("};", start);
    check(start != std::string::npos && end != std::string::npos,
          "ReplicateResult is declared in the header");
    const std::string result_body = header.substr(start, end - start);
    check(result_body.find("ParameterBank") == std::string::npos,
          "ReplicateResult exposes no parameter bank");
    check(result_body.find("ParameterCursor") == std::string::npos,
          "ReplicateResult exposes no cursor");
    check(header.find("ParameterCursor") == std::string::npos,
          "the simulation header exposes no cursor type at all");

    for (const std::string& token : {std::string("refill"), std::string("rewind"),
                                     std::string("reset")}) {
        check(source.find(token) == std::string::npos,
              "simulation.cpp contains no \"" + token + "\" operation");
    }
    check(source.find("const sampling_plan::ParameterBanks banks =") != std::string::npos,
          "the banks are a const local of the run, destroyed when it returns");
}

// ------------------------------------------------------------------
// Section 5: preflight rejections that must happen before allocation.
// ------------------------------------------------------------------

void test_invalid_evaluation_order_rejected() {
    using ecai_cpp2::simulation::run_replicate_from_banks;
    ReplicateOptions options;
    options.evaluation_order = static_cast<EvaluationOrder>(999);
    check(throws_invalid_argument([&] {
              run_replicate_from_banks(fixture_plan(), fixture_banks(), fixture_scenario(),
                                       options);
          }),
          "an out-of-range EvaluationOrder is rejected before any output allocation");
}

void test_run_replicate_from_plan_context_mismatch() {
    using ecai_cpp2::simulation::run_replicate_from_plan;
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");
    const EffectiveConfig config = make_config(patient, 16);
    const ecai_cpp2::seed_manager::SeedManager seeds(485);
    const EffectiveScenario scenario = ecai_cpp2::simulation::resolve_scenario(patient, config);
    const SamplingPlan plan = ecai_cpp2::sampling_plan::build_sampling_plan(
        seeds, config, 0, scenario.simulations_n, scenario.occupants,
        scenario.community_infection_rate, ecai_cpp2::sampling_plan::ResourcePolicy{});

    check(!throws_invalid_argument([&] {
              run_replicate_from_plan(seeds, config, 0, plan, patient, scenario);
          }),
          "a config/category/scenario/plan that all agree is accepted");

    {
        EffectiveScenario bad = scenario;
        bad.category_name = "Classroom";
        check(throws_invalid_argument(
                  [&] { run_replicate_from_plan(seeds, config, 0, plan, patient, bad); }),
              "a scenario category name disagreeing with the category record is rejected before "
              "build_parameter_banks");
    }
    {
        EffectiveConfig bad = config;
        bad.category_name = "Classroom";
        check(throws_invalid_argument(
                  [&] { run_replicate_from_plan(seeds, bad, 0, plan, patient, scenario); }),
              "a config category name disagreeing with the category record is rejected");
    }
    {
        EffectiveConfig bad = config;
        bad.i0 = config.i0 + 1;
        check(throws_invalid_argument(
                  [&] { run_replicate_from_plan(seeds, bad, 0, plan, patient, scenario); }),
              "a config I0 disagreeing with the scenario is rejected");
    }
    {
        EffectiveConfig bad = config;
        bad.volume_m3 = config.volume_m3 + 1.0;
        check(throws_invalid_argument(
                  [&] { run_replicate_from_plan(seeds, bad, 0, plan, patient, scenario); }),
              "a config VOL disagreeing with the scenario is rejected");
    }
    {
        EffectiveConfig bad = config;
        bad.mask_efficiency = config.mask_efficiency < 0.5 ? config.mask_efficiency + 0.1
                                                            : config.mask_efficiency - 0.1;
        check(throws_invalid_argument(
                  [&] { run_replicate_from_plan(seeds, bad, 0, plan, patient, scenario); }),
              "a config mask efficiency disagreeing with the scenario is rejected");
    }
    {
        EffectiveConfig bad = config;
        bad.community_infection_rate = config.community_infection_rate < 0.5
                                            ? config.community_infection_rate + 0.1
                                            : config.community_infection_rate - 0.1;
        check(throws_invalid_argument(
                  [&] { run_replicate_from_plan(seeds, bad, 0, plan, patient, scenario); }),
              "a config community infection rate disagreeing with the scenario is rejected");
    }
    {
        EffectiveConfig bad = config;
        bad.target_probability = config.target_probability + 0.0001;
        check(throws_invalid_argument(
                  [&] { run_replicate_from_plan(seeds, bad, 0, plan, patient, scenario); }),
              "a config target_probability disagreeing with the scenario is rejected");
    }
    {
        EffectiveConfig bad = config;
        bad.duration_hours = config.duration_hours + 1.0;
        check(throws_invalid_argument(
                  [&] { run_replicate_from_plan(seeds, bad, 0, plan, patient, scenario); }),
              "a config duration_hours disagreeing with the scenario is rejected");
    }
    {
        EffectiveConfig bad = config;
        bad.percentile = config.percentile < 0.5 ? config.percentile + 0.01 : config.percentile - 0.01;
        check(throws_invalid_argument(
                  [&] { run_replicate_from_plan(seeds, bad, 0, plan, patient, scenario); }),
              "a config percentile disagreeing with the scenario is rejected");
    }
    {
        EffectiveConfig bad = config;
        bad.simulations_n = config.simulations_n + 1;
        check(throws_invalid_argument(
                  [&] { run_replicate_from_plan(seeds, bad, 0, plan, patient, scenario); }),
              "a config simulations_n disagreeing with the plan's N is rejected");
    }
}

void test_interleaved_schedule_pinned() {
    using ecai_cpp2::simulation::run_replicate_from_banks;
    for (const std::size_t n : {std::size_t{4}, std::size_t{5}}) {
        const SamplingPlan plan = all_active_plan(n);
        const ParameterBanks banks = all_active_banks(n);
        EffectiveScenario scenario = fixture_scenario();
        scenario.simulations_n = n;

        ReplicateOptions forward;
        forward.collect_audit = true;
        forward.evaluation_order = EvaluationOrder::forward;
        ReplicateOptions interleaved = forward;
        interleaved.evaluation_order = EvaluationOrder::interleaved;

        const ReplicateResult f = run_replicate_from_banks(plan, banks, scenario, forward);
        const ReplicateResult i = run_replicate_from_banks(plan, banks, scenario, interleaved);

        const std::vector<std::size_t> expected_order =
            (n == 4) ? std::vector<std::size_t>{0, 3, 1, 2}
                     : std::vector<std::size_t>{0, 4, 1, 3, 2};
        check(i.audit.evaluation_order == expected_order,
              "interleaved schedule for N = " + std::to_string(n) +
                  " active rows matches the pinned permutation exactly");

        const std::set<std::size_t> unique_indices(i.audit.evaluation_order.begin(),
                                                    i.audit.evaluation_order.end());
        check(unique_indices.size() == n,
              "interleaved schedule for N = " + std::to_string(n) + " has no repeated index");
        bool all_present = true;
        for (std::size_t idx = 0; idx < n; ++idx) {
            if (unique_indices.find(idx) == unique_indices.end()) all_present = false;
        }
        check(all_present, "interleaved schedule for N = " + std::to_string(n) +
                                " is a complete permutation of every active row");

        check(bitwise_equal(f.ecai_l_s_person, i.ecai_l_s_person),
              "interleaved schedule for N = " + std::to_string(n) +
                  " active rows reproduces the forward output exactly");
    }
}

// ------------------------------------------------------------------
// Section 6: probability mode (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md
// section 4.2/4.3).
// ------------------------------------------------------------------

void test_probability_injected_fixture() {
    ReplicateOptions options;
    options.collect_audit = true;
    const EffectiveScenario scenario = fixture_probability_scenario();
    const ReplicateResult result = ecai_cpp2::simulation::run_replicate_from_banks(
        fixture_plan(), fixture_banks(), scenario, options);

    check(result.mode == OutputMode::probability,
          "probability-mode replicate reports mode == probability");
    check(result.probability.size() == 5,
          "probability-mode fixture emits exactly N = 5 probabilities");
    check(result.ecai_l_s_person.empty(),
          "probability-mode replicate leaves ecai_l_s_person empty");
    check(result.percentile_l_s_person == 0.0,
          "probability-mode replicate leaves percentile_l_s_person at 0");

    const std::vector<double> expected =
        fixture_expected_probabilities(scenario.input_ecai_l_s_person);
    for (std::size_t i = 0; i < expected.size(); ++i) {
        check(approx(result.probability[i], expected[i]),
              "probability fixture row " + std::to_string(i) +
                  " equals its independently computed probability");
    }

    check(is_positive_zero(result.probability[0]) && is_positive_zero(result.probability[3]),
          "zero-infector rows 0 and 3 are exactly +0.0 in probability mode");
    check(result.probability[1] > 0.0 && result.probability[2] > 0.0 &&
              result.probability[4] > 0.0,
          "active rows produce strictly positive probabilities (no negative clamping in "
          "probability mode)");
    bool all_in_range = true;
    for (const double v : result.probability) {
        if (!std::isfinite(v) || !(v >= 0.0) || !(v <= 1.0)) all_in_range = false;
    }
    check(all_in_range, "every probability lies in the finite closed range [0, 1]");

    check(result.counters.negative_clamp_count == 0, "probability mode never negative-clamps");
    check(result.counters.ecai_evaluations == 0, "probability mode makes no ECAi evaluation");
    check(result.counters.probability_evaluations == 3,
          "one probability evaluation per active row");
    check(result.counters.saturated_count == 0,
          "the fixture's moderate input ECAi does not saturate any row");
    check(result.counters.output_zero_count == 2,
          "probability mode has no clamp source, so output zeros equal zero-infector rows only");
    check(result.counters.qer_evaluations == 4, "one QER evaluation per realized infector");

    check(approx(result.percentile_probability, expected_percentile(expected, 0.96)),
          "probability fixture percentile matches the independently interpolated value");
    check(result.percentile == 0.96, "percentile echoes the scenario's requested percentile");

    check(result.audit.rows[1].probability > 0.0 && !result.audit.rows[1].saturated,
          "audit rows record the per-row probability and saturation flag");
    check(is_positive_zero(result.audit.rows[0].probability) && !result.audit.rows[0].saturated,
          "a zero-infector row's audit entry carries +0.0 and is not saturated");
}

void test_cross_mode_identity() {
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");
    const EffectiveConfig config = make_config(patient, 64);

    ReplicateOptions options;
    options.collect_audit = true;

    EffectiveOverrides ecai_overrides;  // default mode == ecai
    EffectiveOverrides probability_overrides;
    probability_overrides.mode = OutputMode::probability;

    const ReplicateResult ecai_result =
        ecai_cpp2::simulation::run_replicate(485, config, 3, patient, ecai_overrides, options);
    const ReplicateResult probability_result = ecai_cpp2::simulation::run_replicate(
        485, config, 3, patient, probability_overrides, options);

    check(ecai_result.mode == OutputMode::ecai, "the ECAi-mode run reports mode == ecai");
    check(probability_result.mode == OutputMode::probability,
          "the probability-mode run reports mode == probability");

    // Neither mode nor input ECAi is a fingerprint field, so identical seeds
    // must yield identical plans, banks, assignments and cursors
    // (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md section 6).
    check(ecai_result.infector_counts == probability_result.infector_counts,
          "both modes realize the identical Stage 6 infector-count vector");
    check(ecai_result.planned_bank_lengths == probability_result.planned_bank_lengths &&
              ecai_result.final_cursors == probability_result.final_cursors &&
              ecai_result.bank_generated == probability_result.bank_generated,
          "both modes plan, fill and exhaust identical bank lengths and cursors");
    check(ecai_result.cursors_exhausted && probability_result.cursors_exhausted,
          "both modes exhaust every cursor exactly");
    check(ecai_result.counters.active_count == probability_result.counters.active_count &&
              ecai_result.counters.k_total == probability_result.counters.k_total &&
              ecai_result.counters.zero_infector_count ==
                  probability_result.counters.zero_infector_count,
          "both modes agree on active/k_total/zero-infector counts");
    check(ecai_result.counters.qer_evaluations == probability_result.counters.qer_evaluations,
          "both modes evaluate the identical number of QERs");

    check(ecai_result.audit.consumption_ids == probability_result.audit.consumption_ids,
          "both modes consume the identical parameter values in the identical order");
    check(ecai_result.audit.evaluation_order == probability_result.audit.evaluation_order,
          "both modes evaluate the same active rows in the same schedule");

    bool rows_match = ecai_result.audit.rows.size() == probability_result.audit.rows.size();
    bool qer_sums_match = true;
    for (std::size_t row = 0; rows_match && row < ecai_result.audit.rows.size(); ++row) {
        const auto& e = ecai_result.audit.rows[row];
        const auto& p = probability_result.audit.rows[row];
        if (e.simulation_index != p.simulation_index || e.infector_count != p.infector_count ||
            e.active != p.active || e.susceptible_ids != p.susceptible_ids ||
            e.infector_ids != p.infector_ids) {
            rows_match = false;
        }
        if (!bitwise_equal(e.qer_sum_per_h, p.qer_sum_per_h)) qer_sums_match = false;
    }
    check(rows_match, "both modes assign identical people (ids) to identical simulation rows");
    check(qer_sums_match, "both modes compute bit-identical per-row QER sums");

    // Mode isolation: each mode populates only its own pair of fields.
    check(ecai_result.probability.empty() && ecai_result.percentile_probability == 0.0 &&
              ecai_result.counters.probability_evaluations == 0 &&
              ecai_result.counters.saturated_count == 0,
          "the ECAi-mode result carries no probability-mode field");
    check(probability_result.ecai_l_s_person.empty() &&
              probability_result.percentile_l_s_person == 0.0 &&
              probability_result.counters.negative_clamp_count == 0 &&
              probability_result.counters.ecai_evaluations == 0,
          "the probability-mode result carries no ECAi-mode field");

    check(!probability_result.probability.empty() &&
              probability_result.probability.size() == ecai_result.ecai_l_s_person.size(),
          "the probability-mode result emits as many probabilities as the ECAi-mode result "
          "emits ECAi values");

    bool zero_infectors_are_positive_zero = true;
    bool all_in_range = true;
    for (std::size_t row = 0; row < probability_result.infector_counts.size(); ++row) {
        const double v = probability_result.probability[row];
        if (probability_result.infector_counts[row] == 0 && !is_positive_zero(v)) {
            zero_infectors_are_positive_zero = false;
        }
        if (!std::isfinite(v) || !(v >= 0.0) || !(v <= 1.0)) all_in_range = false;
    }
    check(zero_infectors_are_positive_zero,
          "every zero-infector row is exactly +0.0 in probability mode too");
    check(all_in_range, "every probability-mode output lies in the finite closed range [0, 1]");

    check(approx(probability_result.percentile_probability,
                 expected_percentile(probability_result.probability, 0.96)),
          "the probability-mode percentile matches the independently interpolated Type-7 value");
}

// ------------------------------------------------------------------
// Section 7: risk-averse infector policy simulation integration
// (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md Phase 5). Phase 3/4 already cover
// CLI parsing/resolution and the Stage 6 plan-level conditional-PPF
// semantics in their own test files; these tests exercise the policy through
// the simulation-layer entry points (resolve_scenario, run_replicate,
// run_replicate_from_plan, run_replicate_from_banks, run_replicates).
// ------------------------------------------------------------------

EffectiveOverrides risk_averse_overrides() {
    EffectiveOverrides overrides;
    overrides.infector_policy = InfectorPolicy::risk_averse;
    overrides.infector_policy_explicit = true;
    return overrides;
}

// run_replicate must forward the resolved policy into build_sampling_plan
// (previously it silently always built a community plan, ignoring
// scenario.infector_policy entirely -- the defect this phase fixes). Proven
// by reconstructing both policies' plans independently, with the same seed
// and config, and checking run_replicate's realized counts against the
// risk-averse one specifically, not the community one.
void test_risk_averse_policy_forwarded_into_plan() {
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");
    const EffectiveConfig config = make_config(patient, 200);
    const EffectiveOverrides ra = risk_averse_overrides();

    const ReplicateResult result = ecai_cpp2::simulation::run_replicate(485, config, 0, patient, ra);
    check(result.infector_policy == InfectorPolicy::risk_averse,
          "a risk-averse run reports result.infector_policy == risk_averse");

    const ecai_cpp2::seed_manager::SeedManager seeds(485);
    const SamplingPlan community_plan = ecai_cpp2::sampling_plan::build_sampling_plan(
        seeds, config, 0, 200, patient.i0, patient.community_infection_rate,
        ecai_cpp2::sampling_plan::ResourcePolicy{}, ecai_cpp2::sampling_plan::InfectorPolicy::community);
    const SamplingPlan risk_averse_plan = ecai_cpp2::sampling_plan::build_sampling_plan(
        seeds, config, 0, 200, patient.i0, patient.community_infection_rate,
        ecai_cpp2::sampling_plan::ResourcePolicy{},
        ecai_cpp2::sampling_plan::InfectorPolicy::risk_averse);

    check(community_plan.n_zero > 0,
          "sanity check: Patient room's small I0/low community rate really does produce "
          "zero-infector rows under community policy for this fixture");
    check(result.infector_counts == risk_averse_plan.counts,
          "run_replicate's realized counts match the independently rebuilt risk-averse plan, "
          "not the community plan");
    check(result.infector_counts != community_plan.counts,
          "run_replicate's risk-averse counts differ from what a community plan would have "
          "produced for the identical seed/config");

    check(risk_averse_plan.lhs_generated == 200 && risk_averse_plan.lhs_consumed == 200,
          "the risk-averse plan generates and consumes exactly N count-stream uniforms");
    check(community_plan.lhs_generated == 200 && community_plan.lhs_consumed == 200,
          "the community plan generates and consumes exactly N count-stream uniforms");

    // Identical count-stream uniforms (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md
    // Phase 4/7): the same seed and LHS list feed both policies, so a
    // risk-averse count is never lower than the community count the same u
    // produced, and risk-averse never contains a zero.
    bool pointwise_at_least = true;
    bool risk_averse_never_zero = true;
    for (std::size_t i = 0; i < risk_averse_plan.counts.size(); ++i) {
        if (risk_averse_plan.counts[i] < community_plan.counts[i]) pointwise_at_least = false;
        if (risk_averse_plan.counts[i] < 1) risk_averse_never_zero = false;
    }
    check(pointwise_at_least,
          "every risk-averse count is at least its matched community count (identical "
          "count-stream uniforms)");
    check(risk_averse_never_zero, "no risk-averse count is zero");
}

// Risk-averse ECAi and probability modes must both produce no zero-infector
// rows, with every simulation active and every count in [1, I0].
void test_risk_averse_no_zero_rows_ecai_and_probability() {
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");
    const EffectiveConfig config = make_config(patient, 300);

    EffectiveOverrides ecai_ra = risk_averse_overrides();
    const ReplicateResult ecai_result =
        ecai_cpp2::simulation::run_replicate(485, config, 0, patient, ecai_ra);
    check(ecai_result.mode == OutputMode::ecai, "the risk-averse ECAi run resolves to ECAi mode");
    check(ecai_result.counters.zero_infector_count == 0,
          "risk-averse ECAi contains no zero-infector rows");
    check(ecai_result.counters.active_count == 300 && ecai_result.counters.k_total >= 300,
          "risk-averse ECAi marks every row active with K_total >= N");
    bool ecai_counts_bounded = true;
    for (const long long c : ecai_result.infector_counts) {
        if (c < 1 || c > patient.i0) ecai_counts_bounded = false;
    }
    check(ecai_counts_bounded, "every risk-averse ECAi count lies in [1, I0]");

    EffectiveOverrides prob_ra = risk_averse_overrides();
    prob_ra.mode = OutputMode::probability;
    const ReplicateResult prob_result =
        ecai_cpp2::simulation::run_replicate(485, config, 0, patient, prob_ra);
    check(prob_result.mode == OutputMode::probability,
          "the risk-averse probability run resolves to probability mode");
    check(prob_result.counters.zero_infector_count == 0,
          "risk-averse probability contains no zero-infector rows");
    check(prob_result.counters.active_count == 300 && prob_result.counters.k_total >= 300,
          "risk-averse probability marks every row active with K_total >= N");
    bool prob_counts_bounded = true;
    for (const long long c : prob_result.infector_counts) {
        if (c < 1 || c > patient.i0) prob_counts_bounded = false;
    }
    check(prob_counts_bounded, "every risk-averse probability count lies in [1, I0]");
    bool all_in_range = true;
    for (const double v : prob_result.probability) {
        if (!std::isfinite(v) || !(v >= 0.0) || !(v <= 1.0)) all_in_range = false;
    }
    check(all_in_range, "every risk-averse probability output is finite and in [0, 1]");
}

// Assignments, bank cursors and evaluation counters must reconcile exactly
// for a risk-averse replicate, exactly as they do for community.
void test_risk_averse_assignments_bank_cursors_and_qer_reconcile() {
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");
    const EffectiveConfig config = make_config(patient, 150);
    ReplicateOptions options;
    options.collect_audit = true;
    const EffectiveOverrides ra = risk_averse_overrides();

    const ReplicateResult result =
        ecai_cpp2::simulation::run_replicate(485, config, 0, patient, ra, options);

    check(result.audit.rows.size() == 150, "the audit records one row per simulation");
    bool every_row_active = true;
    bool infector_counts_match_ids = true;
    for (std::size_t row = 0; row < result.audit.rows.size(); ++row) {
        const auto& audit_row = result.audit.rows[row];
        if (!audit_row.active) every_row_active = false;
        if (audit_row.infector_ids.size() !=
            static_cast<std::size_t>(result.infector_counts[row])) {
            infector_counts_match_ids = false;
        }
    }
    check(every_row_active, "every simulation is active under risk-averse policy");
    check(infector_counts_match_ids,
          "every row's audited infector-id count equals its authoritative infector count");

    check(result.cursors_exhausted, "every bank cursor is exactly exhausted");
    bool bank_lengths_reconcile = true;
    for (BankId id : kAllBankIds) {
        const std::size_t ordinal = static_cast<std::size_t>(id);
        const std::size_t expected =
            ecai_cpp2::sampling_plan::is_infector_bank(id) ? result.counters.k_total
                                                            : result.counters.active_count;
        if (result.planned_bank_lengths[ordinal] != expected ||
            result.final_cursors[ordinal] != expected || !result.bank_generated[ordinal]) {
            bank_lengths_reconcile = false;
        }
    }
    check(bank_lengths_reconcile,
          "every bank is planned, generated and exhausted at N_active/K_total exactly");

    check(result.counters.qer_evaluations == result.counters.k_total,
          "one QER evaluation per realized infector");
    check(result.counters.ecai_evaluations == result.counters.active_count &&
              result.counters.active_count == 150,
          "one ECAi evaluation per active row, and every row is active");
    check(result.counters.susceptible_values_consumed == 3 * result.counters.active_count,
          "3 susceptible banks x N_active");
    check(result.counters.infector_values_consumed == 8 * result.counters.k_total,
          "8 infector banks x K_total");
}

// Risk-averse probability outputs stay bounded and can genuinely saturate;
// no zero-infector short-circuit ever produces a spurious +0.0/1.0.
void test_risk_averse_probability_saturation_and_bounds() {
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");
    EffectiveConfig config = make_config(patient, 100);
    config.duration_hours = 1.0e6;  // forces saturation for at least some rows
    EffectiveOverrides ra = risk_averse_overrides();
    ra.mode = OutputMode::probability;

    const ReplicateResult result = ecai_cpp2::simulation::run_replicate(485, config, 0, patient, ra);
    check(result.counters.zero_infector_count == 0,
          "the saturating risk-averse probability run still has no zero-infector rows");
    bool all_finite_in_range = true;
    bool no_negative_zero = true;
    for (const double v : result.probability) {
        if (!std::isfinite(v) || !(v >= 0.0) || !(v <= 1.0)) all_finite_in_range = false;
        if (v == 0.0 && std::signbit(v)) no_negative_zero = false;
    }
    check(all_finite_in_range, "every saturating risk-averse probability stays finite and in [0, 1]");
    check(no_negative_zero, "no risk-averse probability output is a negative zero");
    check(result.counters.saturated_count > 0,
          "the extreme-duration fixture actually exercises numerical saturation");
    check(result.counters.negative_clamp_count == 0,
          "probability mode never negative-clamps, risk-averse included");
}

// The forward probability calculation must stay independent of
// target_probability (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md Phase 5): it
// is an ECAi-mode comparator input only, and must not leak into probability
// mode under either infector policy.
void test_probability_independent_of_target_probability() {
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");

    EffectiveConfig config_a = make_config(patient, 80);
    config_a.target_probability = 0.001;
    EffectiveConfig config_b = config_a;
    config_b.target_probability = 0.5;

    EffectiveOverrides prob_ra = risk_averse_overrides();
    prob_ra.mode = OutputMode::probability;
    const ReplicateResult ra_a =
        ecai_cpp2::simulation::run_replicate(485, config_a, 0, patient, prob_ra);
    const ReplicateResult ra_b =
        ecai_cpp2::simulation::run_replicate(485, config_b, 0, patient, prob_ra);
    check(ra_a.mode == OutputMode::probability && ra_b.mode == OutputMode::probability,
          "target_probability changes do not change the resolved probability mode");
    check(ra_a.infector_policy == InfectorPolicy::risk_averse &&
              ra_b.infector_policy == InfectorPolicy::risk_averse,
          "target_probability changes do not change the resolved risk-averse policy");
    check(ra_a.infector_counts.size() == 80 && ra_b.infector_counts.size() == 80,
          "target_probability changes preserve exactly N probability rows");
    check(ra_a.probability.size() == 80 && ra_b.probability.size() == 80,
          "target_probability changes preserve probability vector size");
    check(ra_a.counters.zero_infector_count == 0 && ra_b.counters.zero_infector_count == 0,
          "target_probability changes preserve the risk-averse no-zero invariant");
    check(ra_a.counters.qer_evaluations == ra_a.counters.k_total &&
              ra_b.counters.qer_evaluations == ra_b.counters.k_total,
          "target_probability changes preserve one QER evaluation per realized infector");

    EffectiveOverrides prob_community;
    prob_community.mode = OutputMode::probability;
    const ReplicateResult community_a =
        ecai_cpp2::simulation::run_replicate(485, config_a, 0, patient, prob_community);
    const ReplicateResult community_b =
        ecai_cpp2::simulation::run_replicate(485, config_b, 0, patient, prob_community);
    check(community_a.mode == OutputMode::probability &&
              community_b.mode == OutputMode::probability,
          "target_probability changes preserve community probability mode");
    check(community_a.infector_counts.size() == 80 && community_b.infector_counts.size() == 80,
          "target_probability changes preserve community probability row count");
    check(community_a.counters.qer_evaluations == community_a.counters.k_total &&
              community_b.counters.qer_evaluations == community_b.counters.k_total,
          "target_probability changes preserve community QER accounting");
    check(ra_a.probability != ra_b.probability || community_a.probability != community_b.probability,
          "target_probability is allowed to change streams through the current configuration fingerprint");
    // The forward API has no target_probability field; its direct separation is
    // covered by model tests. These run-level checks deliberately test policy,
    // mode, row-count, and accounting invariants rather than asserting stream
    // identity across changed fingerprints.
}

// Policy disagreement between a resolved scenario and an independently built
// plan (a tampered or stale pairing, or a mixed-policy misuse) must be
// rejected before build_parameter_banks allocates anything, and again if a
// caller bypasses run_replicate_from_plan and hands a mismatched plan/bank
// pair directly to run_replicate_from_banks.
void test_policy_scenario_plan_mismatch_rejected() {
    using ecai_cpp2::simulation::run_replicate_from_banks;
    using ecai_cpp2::simulation::run_replicate_from_plan;
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");
    const ecai_cpp2::seed_manager::SeedManager seeds(485);

    // Community-resolved scenario paired with a risk-averse-built plan.
    {
        const EffectiveConfig config = make_config(patient, 32);
        const EffectiveScenario scenario = ecai_cpp2::simulation::resolve_scenario(patient, config);
        const SamplingPlan mismatched_plan = ecai_cpp2::sampling_plan::build_sampling_plan(
            seeds, config, 0, 32, patient.i0, patient.community_infection_rate,
            ecai_cpp2::sampling_plan::ResourcePolicy{},
            ecai_cpp2::sampling_plan::InfectorPolicy::risk_averse);
        check(throws_invalid_argument([&] {
                  run_replicate_from_plan(seeds, config, 0, mismatched_plan, patient, scenario);
              }),
              "a community scenario paired with a risk-averse plan is rejected by "
              "run_replicate_from_plan before build_parameter_banks");

        const ecai_cpp2::sampling_plan::ParameterBanks banks =
            ecai_cpp2::sampling_plan::build_parameter_banks(seeds, config, 0, mismatched_plan,
                                                              patient);
        check(throws_invalid_argument(
                  [&] { run_replicate_from_banks(mismatched_plan, banks, scenario); }),
              "the same mismatch is independently rejected by run_replicate_from_banks");
    }

    // Risk-averse-resolved scenario paired with a community-built plan.
    {
        const EffectiveConfig config = make_config(patient, 32);
        const EffectiveScenario scenario =
            ecai_cpp2::simulation::resolve_scenario(patient, config, risk_averse_overrides());
        const SamplingPlan mismatched_plan = ecai_cpp2::sampling_plan::build_sampling_plan(
            seeds, config, 0, 32, patient.i0, patient.community_infection_rate,
            ecai_cpp2::sampling_plan::ResourcePolicy{},
            ecai_cpp2::sampling_plan::InfectorPolicy::community);
        check(throws_invalid_argument([&] {
                  run_replicate_from_plan(seeds, config, 0, mismatched_plan, patient, scenario);
              }),
              "a risk-averse scenario paired with a community plan is rejected by "
              "run_replicate_from_plan before build_parameter_banks");

        const ecai_cpp2::sampling_plan::ParameterBanks banks =
            ecai_cpp2::sampling_plan::build_parameter_banks(seeds, config, 0, mismatched_plan,
                                                              patient);
        check(throws_invalid_argument(
                  [&] { run_replicate_from_banks(mismatched_plan, banks, scenario); }),
              "the same mismatch is independently rejected by run_replicate_from_banks");
    }

    // An unknown policy value on the scenario is rejected by validate_scenario
    // before any plan/policy cross-check even runs.
    {
        EffectiveScenario bad = fixture_scenario();
        bad.infector_policy = static_cast<InfectorPolicy>(99);
        check(throws_invalid_argument(
                  [&] { run_replicate_from_banks(fixture_plan(), fixture_banks(), bad); }),
              "an unknown InfectorPolicy on the scenario is rejected");
    }
}

// Risk-averse's own domain (I0 >= 1, community rate > 0) is enforced
// authoritatively at the simulation layer, not only by the CLI or the Stage 6
// plan builder, and rejects before any allocation.
void test_risk_averse_invalid_domain_rejected_at_simulation_layer() {
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");

    {
        EffectiveOverrides overrides = risk_averse_overrides();
        overrides.community_infection_rate = 0.0;
        const EffectiveConfig config = make_config(patient, 16, overrides);
        check(throws_invalid_argument([&] {
                  ecai_cpp2::simulation::resolve_scenario(patient, config, overrides);
              }),
              "risk-averse with community_infection_rate == 0 (p = 0) is rejected at "
              "resolve_scenario");
        check(throws_invalid_argument([&] {
                  ecai_cpp2::simulation::run_replicate(485, config, 0, patient, overrides);
              }),
              "risk-averse with p = 0 is rejected by run_replicate before any allocation");
    }
    {
        EffectiveOverrides overrides = risk_averse_overrides();
        overrides.occupants = 0;
        const EffectiveConfig config = make_config(patient, 16, overrides);
        check(throws_invalid_argument([&] {
                  ecai_cpp2::simulation::resolve_scenario(patient, config, overrides);
              }),
              "risk-averse with I0 == 0 is rejected at resolve_scenario");
        check(throws_invalid_argument([&] {
                  ecai_cpp2::simulation::run_replicate(485, config, 0, patient, overrides);
              }),
              "risk-averse with I0 = 0 is rejected by run_replicate before any allocation");
    }
}

// run_replicates' batch-level policy consistency: every replicate in a
// risk-averse batch resolves to the same policy and carries no zero-infector
// row, and an all-community batch remains an unaffected regression baseline.
void test_run_replicates_batch_policy_consistency() {
    const OccupancyCategory patient = lookup_occupancy_category("Patient room");

    ecai_cpp2::simulation::EffectiveConfiguration ra_config;
    ra_config.category_name = patient.canonical_name;
    ra_config.i0 = patient.i0;
    ra_config.volume_m3 = patient.volume_m3;
    ra_config.mask_efficiency = patient.mask_efficiency;
    ra_config.community_infection_rate = patient.community_infection_rate;
    ra_config.target_probability = 0.001;
    ra_config.duration_hours = 1.0;
    ra_config.percentile = 0.96;
    ra_config.simulations_n = 50;
    ra_config.replicates_r = 5;
    const ecai_cpp2::simulation::ReplicateBatch ra_batch = ecai_cpp2::simulation::run_replicates(
        485, ra_config, patient, risk_averse_overrides());
    check(ra_batch.replicates.size() == 5, "a risk-averse batch preserves exactly R = 5 results");
    bool every_replicate_risk_averse = true;
    bool every_replicate_zero_free = true;
    for (const ReplicateResult& replicate : ra_batch.replicates) {
        if (replicate.infector_policy != InfectorPolicy::risk_averse) {
            every_replicate_risk_averse = false;
        }
        if (replicate.counters.zero_infector_count != 0) every_replicate_zero_free = false;
    }
    check(every_replicate_risk_averse,
          "every replicate in the risk-averse batch resolves to InfectorPolicy::risk_averse");
    check(every_replicate_zero_free,
          "every replicate in the risk-averse batch has zero zero-infector rows");
    check(std::isfinite(ra_batch.summary.mean_l_s_person),
          "the risk-averse batch summary mean is finite");

    ecai_cpp2::simulation::EffectiveConfiguration community_config;
    community_config.category_name = patient.canonical_name;
    community_config.i0 = patient.i0;
    community_config.volume_m3 = patient.volume_m3;
    community_config.mask_efficiency = patient.mask_efficiency;
    community_config.community_infection_rate = patient.community_infection_rate;
    community_config.target_probability = 0.001;
    community_config.duration_hours = 1.0;
    community_config.percentile = 0.96;
    community_config.simulations_n = 50;
    community_config.replicates_r = 3;
    const ecai_cpp2::simulation::ReplicateBatch community_batch =
        ecai_cpp2::simulation::run_replicates(485, community_config, patient);
    check(community_batch.replicates.size() == 3,
          "an unaffected community batch preserves exactly R = 3 results (regression baseline)");
    bool every_replicate_community = true;
    for (const ReplicateResult& replicate : community_batch.replicates) {
        if (replicate.infector_policy != InfectorPolicy::community) {
            every_replicate_community = false;
        }
    }
    check(every_replicate_community,
          "every replicate in an un-overridden batch resolves to InfectorPolicy::community");
}

}  // namespace

int main() {
    test_injected_fixture();
    test_scenario_resolution();
    test_seeded_replicate_determinism();
    test_pinned_seeded_replicate();
    test_all_zero_replicate();
    test_person_identity_and_traversal_order();
    test_large_replicate_smoke();
    test_percentile_convention();
    test_rejects_malformed_inputs();
    test_resource_limits_and_overflow();
    test_banks_are_not_exposed_after_traversal();
    test_invalid_evaluation_order_rejected();
    test_run_replicate_from_plan_context_mismatch();
    test_interleaved_schedule_pinned();
    test_probability_injected_fixture();
    test_cross_mode_identity();
    test_risk_averse_policy_forwarded_into_plan();
    test_risk_averse_no_zero_rows_ecai_and_probability();
    test_risk_averse_assignments_bank_cursors_and_qer_reconcile();
    test_risk_averse_probability_saturation_and_bounds();
    test_probability_independent_of_target_probability();
    test_policy_scenario_plan_mismatch_rejected();
    test_risk_averse_invalid_domain_rejected_at_simulation_layer();
    test_run_replicates_batch_policy_consistency();

    std::cout << g_checks << " checks executed.\n";
    if (g_checks == 0) {
        std::cerr << "FAIL: no checks executed.\n";
        return 1;
    }
    if (g_failures > 0) {
        std::cerr << g_failures << " test(s) failed.\n";
        return 1;
    }
    std::cout << "All tests passed.\n";
    return 0;
}
