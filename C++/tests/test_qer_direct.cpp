// Work package 1 tests for qer_direct: typed QER family/profile selection,
// the immutable Jones-fitted (25 records) and Mikszewski (30 profiles)
// parameter tables, their explicit legacy-short-to-canonical mappings, and
// the pure direct-QER log10-normal transform.
//
// Jones-fitted mu_log10/sigma_log10 values and Mikszewski medians/sigmas are
// typed directly from C++_legacy/src/qer_fitted.cpp (the legacy source of
// truth), not copied from qer_direct.cpp, so this test independently
// verifies the implementation against that source. Canonical category names
// are typed from C++/include/category_data.h's publication-order table.

#include "qer_direct.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "distributions.h"

namespace {

int g_failures = 0;

void check(bool condition, const std::string& description) {
    if (!condition) {
        std::cerr << "FAIL: " << description << "\n";
        ++g_failures;
    } else {
        std::cout << "PASS: " << description << "\n";
    }
}

bool approx(double a, double b, double tol = 1e-9) {
    return std::fabs(a - b) <= tol;
}

template <typename Fn>
bool throws_invalid_argument(Fn&& fn) {
    try {
        fn();
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

}  // namespace

int main() {
    using namespace ecai_cpp2::qer_direct;

    // ================================================================
    // QerFamily: canonical rendering and strict parsing.
    // ================================================================
    {
        check(to_string(QerFamily::Jones) == "jones", "QerFamily::Jones renders \"jones\"");
        check(to_string(QerFamily::JonesFitted) == "jones-fitted",
              "QerFamily::JonesFitted renders \"jones-fitted\"");
        check(to_string(QerFamily::Mikszewski) == "mikszewski",
              "QerFamily::Mikszewski renders \"mikszewski\"");

        check(parse_qer_family("jones") == QerFamily::Jones, "parse \"jones\"");
        check(parse_qer_family("jones-fitted") == QerFamily::JonesFitted,
              "parse \"jones-fitted\"");
        check(parse_qer_family("mikszewski") == QerFamily::Mikszewski, "parse \"mikszewski\"");

        check(throws_invalid_argument([] { parse_qer_family("Jones"); }),
              "parse_qer_family rejects wrong case");
        check(throws_invalid_argument([] { parse_qer_family("jones_fitted"); }),
              "parse_qer_family rejects underscore variant");
        check(throws_invalid_argument([] { parse_qer_family(""); }),
              "parse_qer_family rejects empty string");
        check(throws_invalid_argument([] { parse_qer_family("unknown"); }),
              "parse_qer_family rejects unknown family");
    }

    // ================================================================
    // Jones-fitted table: exactly 25 records, literal mu/sigma from
    // C++/src/qer_fitted.cpp, in category_data.h publication order.
    // ================================================================
    {
        struct Expected {
            std::string canonical_category;
            double mu_log10;
            double sigma_log10;
        };
        const Expected kExpected[25] = {
            {"Cell", -1.123527, 1.428904},
            {"Dayroom", -0.992113, 1.429109},
            {"Food and beverage facilities", -1.093306, 1.427494},
            {"Gym", -1.298071, 1.432531},
            {"Office", -1.123527, 1.428904},
            {"Retail", -1.299721, 1.429381},
            {"Transportation waiting", -1.026098, 1.429956},
            {"Classroom", -1.039436, 1.427156},
            {"Lecture hall", -1.349078, 1.428965},
            {"Manufacturing", -0.970826, 1.430100},
            {"Sorting, packing, light assembly", -0.844935, 1.431221},
            {"Warehouse", -0.947311, 1.432394},
            {"Exam room", -0.645766, 1.430231},
            {"Group treatment area", -1.260304, 1.428842},
            {"Patient room", -0.673900, 1.430657},
            {"Resident room", -0.827809, 1.430225},
            {"Waiting room", -1.299721, 1.429381},
            {"Auditorium", -1.349078, 1.428965},
            {"Place of religious worship", -1.389396, 1.427938},
            {"Museum", -1.399522, 1.429133},
            {"Convention", -1.399522, 1.429133},
            {"Spectator area", -1.093306, 1.427494},
            {"Lobbies", -1.349078, 1.428965},
            {"Common space", -1.349078, 1.428965},
            {"Dwelling unit", -0.509175, 1.561887},
        };

        const auto& table = jones_fitted_table();
        check(table.size() == 25, "jones_fitted_table has exactly 25 records");

        for (std::size_t i = 0; i < table.size() && i < 25; ++i) {
            const std::string tag = kExpected[i].canonical_category;
            check(table[i].canonical_category == kExpected[i].canonical_category,
                  "publication order index " + std::to_string(i) + " is " + tag);
            check(approx(table[i].params.mu_log10, kExpected[i].mu_log10),
                  tag + ": mu_log10");
            check(approx(table[i].params.sigma_log10, kExpected[i].sigma_log10),
                  tag + ": sigma_log10");
        }

        for (const auto& expected : kExpected) {
            const auto& params =
                jones_fitted_params_for_canonical_category(expected.canonical_category);
            check(approx(params.mu_log10, expected.mu_log10),
                  expected.canonical_category + ": lookup mu_log10");
            check(approx(params.sigma_log10, expected.sigma_log10),
                  expected.canonical_category + ": lookup sigma_log10");
        }

        check(throws_invalid_argument(
                  [] { jones_fitted_params_for_canonical_category("Nonexistent Category"); }),
              "jones_fitted_params_for_canonical_category rejects unknown category");
        check(throws_invalid_argument(
                  [] { jones_fitted_params_for_canonical_category(""); }),
              "jones_fitted_params_for_canonical_category rejects empty string");
    }

    // ================================================================
    // Jones-fitted legacy-short-to-canonical mapping: complete (25 entries,
    // one per canonical category, no unknown-category fallback).
    // ================================================================
    {
        struct Alias {
            std::string legacy_short_name;
            std::string canonical_category;
        };
        const Alias kAliases[25] = {
            {"Cell", "Cell"},
            {"Dayroom", "Dayroom"},
            {"Food", "Food and beverage facilities"},
            {"Gym", "Gym"},
            {"Office", "Office"},
            {"Retail", "Retail"},
            {"Transportation", "Transportation waiting"},
            {"Classroom", "Classroom"},
            {"Lecture", "Lecture hall"},
            {"Manufacturing", "Manufacturing"},
            {"Sorting", "Sorting, packing, light assembly"},
            {"Warehouse", "Warehouse"},
            {"Exam", "Exam room"},
            {"Group", "Group treatment area"},
            {"Patient", "Patient room"},
            {"Resident", "Resident room"},
            {"Waiting", "Waiting room"},
            {"Auditorium", "Auditorium"},
            {"Place", "Place of religious worship"},
            {"Museum", "Museum"},
            {"Convention", "Convention"},
            {"Spectator", "Spectator area"},
            {"Lobbies", "Lobbies"},
            {"Common", "Common space"},
            {"Dwelling", "Dwelling unit"},
        };

        const auto& aliases = jones_fitted_legacy_aliases();
        check(aliases.size() == 25, "jones_fitted_legacy_aliases has exactly 25 entries");

        for (const auto& expected : kAliases) {
            bool found = false;
            for (const auto& alias : aliases) {
                if (alias.legacy_short_name == expected.legacy_short_name) {
                    found = true;
                    check(alias.canonical_category == expected.canonical_category,
                          "legacy alias \"" + expected.legacy_short_name + "\" maps to \"" +
                              expected.canonical_category + "\"");
                    break;
                }
            }
            check(found, "legacy alias \"" + expected.legacy_short_name + "\" is present");
        }

        // Every canonical category referenced by the table has exactly one
        // legacy short-name entry mapping to it (no incomplete mapping).
        const auto& table = jones_fitted_table();
        for (const auto& record : table) {
            int match_count = 0;
            for (const auto& alias : aliases) {
                if (alias.canonical_category == record.canonical_category) ++match_count;
            }
            check(match_count == 1, record.canonical_category +
                                         ": exactly one legacy alias maps to this category");
        }

        for (const auto& expected : kAliases) {
            const auto& params = jones_fitted_params_for_legacy_short_name(expected.legacy_short_name);
            const auto& canonical_params =
                jones_fitted_params_for_canonical_category(expected.canonical_category);
            check(approx(params.mu_log10, canonical_params.mu_log10) &&
                      approx(params.sigma_log10, canonical_params.sigma_log10),
                  "legacy short name \"" + expected.legacy_short_name +
                      "\" resolves to the same params as \"" + expected.canonical_category + "\"");
        }

        check(throws_invalid_argument(
                  [] { jones_fitted_params_for_legacy_short_name("Nonexistent"); }),
              "jones_fitted_params_for_legacy_short_name rejects unknown short name");
        check(throws_invalid_argument(
                  [] { jones_fitted_params_for_legacy_short_name("Food and beverage facilities"); }),
              "jones_fitted_params_for_legacy_short_name does not fall back to canonical lookup");
    }

    // ================================================================
    // Mikszewski profile table: exactly 30 profiles (10 organisms x 3
    // activities), literal medians/sigmas from C++/src/qer_fitted.cpp,
    // mu_log10 = log10(median) computed without rounding.
    // ================================================================
    {
        struct Expected {
            MikszewskiOrganism organism;
            MikszewskiActivity activity;
            std::string canonical_name;
            double median;
            double sigma_log10;
        };
        const Expected kExpected[30] = {
            {MikszewskiOrganism::SarsCov1, MikszewskiActivity::Resting,
             "mikszewski-sars-cov-1-resting", 0.0084, 1.30},
            {MikszewskiOrganism::SarsCov1, MikszewskiActivity::StandingSpeaking,
             "mikszewski-sars-cov-1-standing-speaking", 0.042, 1.30},
            {MikszewskiOrganism::SarsCov1, MikszewskiActivity::LightSpeakingLoudly,
             "mikszewski-sars-cov-1-light-speaking-loudly", 0.71, 1.30},
            {MikszewskiOrganism::Mers, MikszewskiActivity::Resting,
             "mikszewski-mers-resting", 0.011, 1.60},
            {MikszewskiOrganism::Mers, MikszewskiActivity::StandingSpeaking,
             "mikszewski-mers-standing-speaking", 0.056, 1.60},
            {MikszewskiOrganism::Mers, MikszewskiActivity::LightSpeakingLoudly,
             "mikszewski-mers-light-speaking-loudly", 0.96, 1.60},
            {MikszewskiOrganism::TbOnTreatment, MikszewskiActivity::Resting,
             "mikszewski-tb-on-treatment-resting", 0.020, 1.40},
            {MikszewskiOrganism::TbOnTreatment, MikszewskiActivity::StandingSpeaking,
             "mikszewski-tb-on-treatment-standing-speaking", 0.098, 1.40},
            {MikszewskiOrganism::TbOnTreatment, MikszewskiActivity::LightSpeakingLoudly,
             "mikszewski-tb-on-treatment-light-speaking-loudly", 1.7, 1.40},
            {MikszewskiOrganism::Influenza, MikszewskiActivity::Resting,
             "mikszewski-influenza-resting", 0.035, 0.84},
            {MikszewskiOrganism::Influenza, MikszewskiActivity::StandingSpeaking,
             "mikszewski-influenza-standing-speaking", 0.17, 0.84},
            {MikszewskiOrganism::Influenza, MikszewskiActivity::LightSpeakingLoudly,
             "mikszewski-influenza-light-speaking-loudly", 3.0, 0.84},
            {MikszewskiOrganism::Coxsackievirus, MikszewskiActivity::Resting,
             "mikszewski-coxsackievirus-resting", 0.062, 1.10},
            {MikszewskiOrganism::Coxsackievirus, MikszewskiActivity::StandingSpeaking,
             "mikszewski-coxsackievirus-standing-speaking", 0.31, 1.10},
            {MikszewskiOrganism::Coxsackievirus, MikszewskiActivity::LightSpeakingLoudly,
             "mikszewski-coxsackievirus-light-speaking-loudly", 5.2, 1.10},
            {MikszewskiOrganism::Rhinovirus, MikszewskiActivity::Resting,
             "mikszewski-rhinovirus-resting", 0.21, 0.83},
            {MikszewskiOrganism::Rhinovirus, MikszewskiActivity::StandingSpeaking,
             "mikszewski-rhinovirus-standing-speaking", 1.0, 0.83},
            {MikszewskiOrganism::Rhinovirus, MikszewskiActivity::LightSpeakingLoudly,
             "mikszewski-rhinovirus-light-speaking-loudly", 18.0, 0.83},
            {MikszewskiOrganism::SarsCov2, MikszewskiActivity::Resting,
             "mikszewski-sars-cov-2-resting", 0.55, 1.20},
            {MikszewskiOrganism::SarsCov2, MikszewskiActivity::StandingSpeaking,
             "mikszewski-sars-cov-2-standing-speaking", 2.7, 1.20},
            {MikszewskiOrganism::SarsCov2, MikszewskiActivity::LightSpeakingLoudly,
             "mikszewski-sars-cov-2-light-speaking-loudly", 46.0, 1.20},
            {MikszewskiOrganism::TbUntreated, MikszewskiActivity::Resting,
             "mikszewski-tb-untreated-resting", 0.62, 1.30},
            {MikszewskiOrganism::TbUntreated, MikszewskiActivity::StandingSpeaking,
             "mikszewski-tb-untreated-standing-speaking", 3.1, 1.30},
            {MikszewskiOrganism::TbUntreated, MikszewskiActivity::LightSpeakingLoudly,
             "mikszewski-tb-untreated-light-speaking-loudly", 52.0, 1.30},
            {MikszewskiOrganism::Adenovirus, MikszewskiActivity::Resting,
             "mikszewski-adenovirus-resting", 0.78, 0.95},
            {MikszewskiOrganism::Adenovirus, MikszewskiActivity::StandingSpeaking,
             "mikszewski-adenovirus-standing-speaking", 3.9, 0.95},
            {MikszewskiOrganism::Adenovirus, MikszewskiActivity::LightSpeakingLoudly,
             "mikszewski-adenovirus-light-speaking-loudly", 66.0, 0.95},
            {MikszewskiOrganism::Measles, MikszewskiActivity::Resting,
             "mikszewski-measles-resting", 3.1, 1.60},
            {MikszewskiOrganism::Measles, MikszewskiActivity::StandingSpeaking,
             "mikszewski-measles-standing-speaking", 15.0, 1.60},
            {MikszewskiOrganism::Measles, MikszewskiActivity::LightSpeakingLoudly,
             "mikszewski-measles-light-speaking-loudly", 260.0, 1.60},
        };

        const auto& table = mikszewski_profile_table();
        check(table.size() == 30, "mikszewski_profile_table has exactly 30 records");

        for (std::size_t i = 0; i < table.size() && i < 30; ++i) {
            const std::string tag = kExpected[i].canonical_name;
            check(table[i].profile.organism == kExpected[i].organism &&
                      table[i].profile.activity == kExpected[i].activity,
                  tag + ": organism/activity");
            check(approx(table[i].median_quanta_per_hour, kExpected[i].median),
                  tag + ": literal median");
            check(canonical_mikszewski_profile_name(table[i].profile) == kExpected[i].canonical_name,
                  tag + ": canonical profile name");
            const double expected_mu = std::log10(kExpected[i].median);
            check(approx(table[i].params.mu_log10, expected_mu, 1e-12),
                  tag + ": mu_log10 == log10(median), unrounded");
            check(approx(table[i].params.sigma_log10, kExpected[i].sigma_log10),
                  tag + ": sigma_log10");
        }

        // Parsing and rendering round-trip for every profile.
        for (const auto& expected : kExpected) {
            const MikszewskiProfile parsed = parse_mikszewski_profile(expected.canonical_name);
            check(parsed.organism == expected.organism && parsed.activity == expected.activity,
                  "parse_mikszewski_profile(\"" + expected.canonical_name + "\") round-trips");

            const auto& params = mikszewski_params_for_profile(parsed);
            check(approx(params.mu_log10, std::log10(expected.median), 1e-12) &&
                      approx(params.sigma_log10, expected.sigma_log10),
                  expected.canonical_name + ": mikszewski_params_for_profile matches table");

            const auto& params_by_name = mikszewski_params_for_canonical_name(expected.canonical_name);
            check(approx(params_by_name.mu_log10, params.mu_log10) &&
                      approx(params_by_name.sigma_log10, params.sigma_log10),
                  expected.canonical_name + ": mikszewski_params_for_canonical_name matches");
        }

        check(throws_invalid_argument([] { parse_mikszewski_profile("mikszewski"); }),
              "parse_mikszewski_profile rejects bare \"mikszewski\" (no default resolution)");
        check(throws_invalid_argument([] { parse_mikszewski_profile("sars-cov-2-standing-speaking"); }),
              "parse_mikszewski_profile rejects profile name missing the mikszewski- prefix");
        check(throws_invalid_argument(
                  [] { parse_mikszewski_profile("mikszewski-sars-cov-2-standing"); }),
              "parse_mikszewski_profile rejects an incomplete activity suffix");
        check(throws_invalid_argument([] { parse_mikszewski_profile(""); }),
              "parse_mikszewski_profile rejects empty string");

        // All 30 canonical names and legacy dictionary keys are distinct.
        std::vector<std::string> names;
        for (const auto& record : table) {
            names.push_back(canonical_mikszewski_profile_name(record.profile));
        }
        bool all_unique = true;
        for (std::size_t i = 0; i < names.size() && all_unique; ++i) {
            for (std::size_t j = i + 1; j < names.size(); ++j) {
                if (names[i] == names[j]) all_unique = false;
            }
        }
        check(all_unique, "all 30 canonical Mikszewski profile names are unique");
    }

    // ================================================================
    // Positivity, finiteness, and parameter-table immutability sanity.
    // ================================================================
    {
        for (const auto& record : jones_fitted_table()) {
            check(std::isfinite(record.params.mu_log10) && std::isfinite(record.params.sigma_log10) &&
                      record.params.sigma_log10 > 0.0,
                  record.canonical_category + ": finite mu_log10 and positive finite sigma_log10");
        }
        for (const auto& record : mikszewski_profile_table()) {
            check(record.median_quanta_per_hour > 0.0 && std::isfinite(record.median_quanta_per_hour),
                  canonical_mikszewski_profile_name(record.profile) + ": positive finite median");
            check(std::isfinite(record.params.mu_log10) && std::isfinite(record.params.sigma_log10) &&
                      record.params.sigma_log10 > 0.0,
                  canonical_mikszewski_profile_name(record.profile) +
                      ": finite mu_log10 and positive finite sigma_log10");
        }
    }

    // ================================================================
    // Pure direct-QER log10-normal transform: identical to
    // distributions::log10normal_ppf(u, mu_log10, sigma_log10), reusing the
    // existing C++ inverse-normal machinery, and rejecting the same
    // out-of-domain u and parameter values.
    // ================================================================
    {
        const DirectQerParams params{-1.0, 1.4};
        const double u_values[] = {1e-6, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 1.0 - 1e-6};
        for (double u : u_values) {
            const double actual = sample_direct_qer_log10normal(u, params);
            const double expected =
                ecai_cpp2::distributions::log10normal_ppf(u, params.mu_log10, params.sigma_log10);
            check(approx(actual, expected, 1e-9),
                  "sample_direct_qer_log10normal(" + std::to_string(u) +
                      ") matches distributions::log10normal_ppf");
            check(actual > 0.0 && std::isfinite(actual),
                  "sample_direct_qer_log10normal(" + std::to_string(u) + ") is positive and finite");
        }

        // Every table record's params round-trip through the transform.
        for (const auto& record : jones_fitted_table()) {
            const double value = sample_direct_qer_log10normal(0.5, record.params);
            check(value > 0.0 && std::isfinite(value),
                  record.canonical_category + ": transform at u=0.5 is positive and finite");
        }
        for (const auto& record : mikszewski_profile_table()) {
            const double value = sample_direct_qer_log10normal(0.5, record.params);
            check(value > 0.0 && std::isfinite(value),
                  canonical_mikszewski_profile_name(record.profile) +
                      ": transform at u=0.5 is positive and finite");
        }

        check(throws_invalid_argument([&] { sample_direct_qer_log10normal(0.0, params); }),
              "sample_direct_qer_log10normal rejects u == 0");
        check(throws_invalid_argument([&] { sample_direct_qer_log10normal(1.0, params); }),
              "sample_direct_qer_log10normal rejects u == 1");
        check(throws_invalid_argument([&] {
                  sample_direct_qer_log10normal(0.5, DirectQerParams{-1.0, -1.0});
              }),
              "sample_direct_qer_log10normal rejects non-positive sigma_log10");
        check(throws_invalid_argument([&] {
                  sample_direct_qer_log10normal(0.5, DirectQerParams{
                                                          std::numeric_limits<double>::quiet_NaN(), 1.4});
              }),
              "sample_direct_qer_log10normal rejects nonfinite mu_log10");
    }

    // ================================================================
    // Table-version provenance: non-empty, distinct identifiers.
    // ================================================================
    {
        check(!std::string(kJonesFittedTableVersion).empty(), "Jones-fitted table version is non-empty");
        check(!std::string(kMikszewskiTableVersion).empty(), "Mikszewski table version is non-empty");
        check(std::string(kJonesFittedTableVersion) != std::string(kMikszewskiTableVersion),
              "Jones-fitted and Mikszewski table versions are distinct");
    }

    if (g_failures > 0) {
        std::cerr << g_failures << " test(s) failed.\n";
        return 1;
    }
    std::cout << "All tests passed.\n";
    return 0;
}
