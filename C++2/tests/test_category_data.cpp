// Stage 2 tests for category_data: strict, literal-value verification of all
// 25 Jones et al. occupancy categories against docs/MODEL_SPECIFICATION.md
// (Tables 1-4, D2, D7, and the frozen-conventions summary), plus canonical
// lookup, normalization, alias, and rejection behavior per
// IMPLEMENTATION_PLAN.md Section 7, Stage 2.
//
// Expected values below are typed directly from MODEL_SPECIFICATION.md, not
// copied from category_data.cpp, so this test independently verifies the
// implementation against the frozen specification.

#include "category_data.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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

struct Expected {
    std::string canonical_name;
    std::string group;
    bool is_health_care_group;
    int i0;
    double volume_m3;
    double probability_susceptible_present;
    double people_outdoor_airflow_Lsp;
    double area_outdoor_airflow_Lsm2;
    double baseline_outdoor_airflow_Lsp;
    double community_infection_rate;
    double mask_efficiency;
    double pbr_gm_m3h;
    double pbr_gsd;
    double activity_breathing_pct;
    double activity_talking_pct;
    double activity_vocalizing_pct;
    double cdrop_gm_e4_per_m3;
    double cdrop_gsd;
    double d_gm_e6_m;
    double d_gsd;
    double published_ecai_Lsp;
    double published_rr;
};

// Table 2 (geometry/occupancy, p.20), Table 3 (breathing/droplet, p.21),
// Table 4 (published adjusted ECAi/RR, p.22), Table 1 (community rate/mask,
// p.15-16), and D2/D7 adjudications, in exact publication order.
const Expected kExpected[25] = {
    {"Cell", "Correctional Facilities", false, 30, 320, 0.25, 2.5, 0.6, 4.9,
     0.01, 0.0, 0.55, 1.2, 80, 20, 0, 14, 1.1, 1.9, 1.1, 15, 0.44},
    {"Dayroom", "Correctional Facilities", false, 30, 600, 0.25, 2.5, 0.3, 3.5,
     0.01, 0.0, 0.55, 1.2, 75, 25, 0, 19, 1.1, 1.9, 1.1, 20, 0.35},
    {"Food and beverage facilities", "Commercial/Retail", false, 50, 360, 0.38,
     3.8, 0.9, 5.6, 0.01, 0.0, 0.55, 1.2, 75, 25, 0, 15, 1.1, 1.9, 1.1, 30,
     0.24},
    {"Gym", "Commercial/Retail", false, 180, 1900, 0.83, 10, 0.3, 11, 0.01,
     0.0, 0.62, 1.3, 100, 0, 0, 9.8, 1.1, 1.8, 1.1, 40, 0.32},
    {"Office", "Commercial/Retail", false, 50, 2700, 0.38, 2.5, 0.3, 8.5, 0.01,
     0.0, 0.55, 1.2, 80, 20, 0, 14, 1.1, 1.9, 1.1, 15, 0.78},
    {"Retail", "Commercial/Retail", false, 150, 6000, 0.77, 3.8, 0.6, 7.8,
     0.01, 0.0, 0.55, 1.2, 95, 5, 0, 11, 1.1, 1.8, 1.1, 20, 0.61},
    {"Transportation waiting", "Commercial/Retail", false, 100, 10000, 0.62,
     3.8, 0.3, 6.8, 0.01, 0.0, 0.64, 1.2, 75, 25, 0, 15, 1.1, 1.9, 1.1, 30,
     0.60},
    {"Classroom", "Educational Facilities", false, 30, 320, 0.25, 5.0, 0.6,
     7.4, 0.01, 0.0, 0.55, 1.2, 70, 30, 0, 17, 1.1, 1.9, 1.1, 20, 0.45},
    {"Lecture hall", "Educational Facilities", false, 150, 1800, 0.77, 2.5,
     0.3, 2.9, 0.01, 0.0, 0.55, 1.2, 100, 0, 0, 9.8, 1.1, 1.8, 1.1, 25, 0.22},
    {"Manufacturing", "Industrial", false, 70, 12000, 0.50, 5.0, 0.3, 9.3,
     0.01, 0.0, 0.73, 1.3, 75, 25, 0, 15, 1.1, 1.9, 1.1, 25, 0.78},
    {"Sorting, packing, light assembly", "Industrial", false, 20, 4800, 0.17,
     5.0, 0.3, 11, 0.01, 0.0, 0.77, 1.3, 75, 25, 0, 19, 1.1, 1.9, 1.1, 10,
     1.0},
    {"Warehouse", "Industrial", false, 20, 1200, 0.17, 5.0, 0.3, 6.5, 0.01,
     0.0, 0.77, 1.3, 75, 25, 0, 15, 1.1, 1.9, 1.1, 10, 0.87},
    {"Exam room", "Health Care", true, 3, 41, 0.057, 7.5, 0.0, 7.5, 0.03, 0.30,
     0.62, 1.3, 100, 0, 0, 32, 1.1, 2.0, 1.1, 20, 0.47},
    {"Group treatment area", "Health Care", true, 20, 270, 0.43, 7.5, 0.0, 7.5,
     0.03, 0.30, 0.55, 1.2, 90, 10, 0, 12, 1.1, 1.8, 1.1, 35, 0.29},
    {"Patient room", "Health Care", true, 3, 81, 0.057, 15, 0.0, 15, 0.03,
     0.30, 0.62, 1.3, 100, 0, 0, 30, 1.1, 2.0, 1.1, 35, 0.53},
    {"Resident room", "Health Care", true, 3, 81, 0.057, 15, 0.0, 15, 0.03,
     0.30, 0.62, 1.3, 100, 0, 0, 21, 1.1, 2.0, 1.1, 25, 0.69},
    {"Waiting room", "Health Care", true, 30, 270, 0.57, 5.0, 0.0, 5.0, 0.03,
     0.30, 0.55, 1.2, 90, 10, 0, 11, 1.1, 1.8, 1.1, 45, 0.16},
    {"Auditorium", "Public Assembly/Sports/Entertainment", false, 150, 4600,
     0.77, 2.5, 0.3, 3.5, 0.01, 0.0, 0.55, 1.2, 90, 10, 0, 9.8, 1.1, 1.8, 1.1,
     25, 0.36},
    {"Place of religious worship", "Public Assembly/Sports/Entertainment",
     false, 180, 2400, 0.83, 2.5, 0.3, 2.8, 0.01, 0.0, 0.50, 1.2, 100, 0, 0,
     9.8, 1.1, 1.8, 1.1, 25, 0.23},
    {"Museum", "Public Assembly/Sports/Entertainment", false, 400, 10000, 0.97,
     3.8, 0.3, 4.6, 0.01, 0.0, 0.49, 1.2, 100, 0, 0, 9.8, 1.1, 1.8, 1.1, 30,
     0.31},
    {"Convention", "Public Assembly/Sports/Entertainment", false, 400, 10000,
     0.97, 3.8, 0.3, 4.6, 0.01, 0.0, 0.49, 1.2, 100, 0, 0, 9.8, 1.1, 1.8, 1.1,
     30, 0.31},
    {"Spectator area", "Public Assembly/Sports/Entertainment", false, 100,
     6300, 0.62, 3.8, 0.3, 4.7, 0.01, 0.0, 0.55, 1.2, 70, 30, 0, 15, 1.1, 1.9,
     1.1, 25, 0.52},
    {"Lobbies", "Public Assembly/Sports/Entertainment", false, 150, 4600, 0.77,
     2.5, 0.3, 3.5, 0.01, 0.0, 0.55, 1.2, 90, 10, 0, 9.8, 1.1, 1.8, 1.1, 25,
     0.34},
    {"Common space", "Residential", false, 150, 4600, 0.77, 2.5, 0.3, 1.0,
     0.01, 0.0, 0.55, 1.2, 90, 10, 0, 9.8, 1.1, 1.8, 1.1, 25, 0.26},
    {"Dwelling unit", "Residential", false, 6, 540, 0.049, 6.8, 0.0, 6.8, 0.01,
     0.0, 0.51, 4.3, 15, 70, 15, 40, 1.1, 2.2, 1.1, 15, 0.80},
};

const std::string kHealthCareNames[5] = {
    "Exam room", "Group treatment area", "Patient room", "Resident room",
    "Waiting room",
};

const std::pair<std::string, std::string> kAliases[10] = {
    {"Food", "Food and beverage facilities"},
    {"Transportation", "Transportation waiting"},
    {"Group", "Group treatment area"},
    {"Patient", "Patient room"},
    {"Resident", "Resident room"},
    {"Place", "Place of religious worship"},
    {"Spectator", "Spectator area"},
    {"Common", "Common space"},
    {"Dwelling", "Dwelling unit"},
    {"Sorting", "Sorting, packing, light assembly"},
};

void check_category_fields(const ecai_cpp2::OccupancyCategory& actual,
                            const Expected& expected) {
    const std::string tag = expected.canonical_name;
    check(actual.canonical_name == expected.canonical_name, tag + ": canonical_name");
    check(actual.group == expected.group, tag + ": group");
    check(actual.is_health_care_group == expected.is_health_care_group,
          tag + ": is_health_care_group");
    check(actual.i0 == expected.i0, tag + ": i0");
    check(approx(actual.volume_m3, expected.volume_m3), tag + ": volume_m3");
    check(approx(actual.probability_susceptible_present,
                 expected.probability_susceptible_present),
          tag + ": probability_susceptible_present");
    check(approx(actual.people_outdoor_airflow_Lsp,
                 expected.people_outdoor_airflow_Lsp),
          tag + ": people_outdoor_airflow_Lsp");
    check(approx(actual.area_outdoor_airflow_Lsm2,
                 expected.area_outdoor_airflow_Lsm2),
          tag + ": area_outdoor_airflow_Lsm2");
    check(approx(actual.baseline_outdoor_airflow_Lsp,
                 expected.baseline_outdoor_airflow_Lsp),
          tag + ": baseline_outdoor_airflow_Lsp");
    check(approx(actual.community_infection_rate,
                 expected.community_infection_rate),
          tag + ": community_infection_rate");
    check(approx(actual.mask_efficiency, expected.mask_efficiency),
          tag + ": mask_efficiency");
    check(approx(actual.pbr_gm_m3h, expected.pbr_gm_m3h), tag + ": pbr_gm_m3h");
    check(approx(actual.pbr_gsd, expected.pbr_gsd), tag + ": pbr_gsd");
    check(approx(actual.activity_breathing_pct, expected.activity_breathing_pct),
          tag + ": activity_breathing_pct");
    check(approx(actual.activity_talking_pct, expected.activity_talking_pct),
          tag + ": activity_talking_pct");
    check(approx(actual.activity_vocalizing_pct,
                 expected.activity_vocalizing_pct),
          tag + ": activity_vocalizing_pct");
    check(approx(actual.cdrop_gm_e4_per_m3, expected.cdrop_gm_e4_per_m3),
          tag + ": cdrop_gm_e4_per_m3");
    check(approx(actual.cdrop_gsd, expected.cdrop_gsd), tag + ": cdrop_gsd");
    check(approx(actual.d_gm_e6_m, expected.d_gm_e6_m), tag + ": d_gm_e6_m");
    check(approx(actual.d_gsd, expected.d_gsd), tag + ": d_gsd");
    check(approx(actual.published_ecai_Lsp, expected.published_ecai_Lsp),
          tag + ": published_ecai_Lsp");
    check(approx(actual.published_rr, expected.published_rr),
          tag + ": published_rr");
}

}  // namespace

int main() {
    const auto& categories = ecai_cpp2::all_occupancy_categories();

    check(categories.size() == 25, "exactly 25 categories");

    for (std::size_t i = 0; i < categories.size() && i < 25; ++i) {
        check(categories[i].canonical_name == kExpected[i].canonical_name,
              "publication order index " + std::to_string(i) + " is " +
                  kExpected[i].canonical_name);
        check_category_fields(categories[i], kExpected[i]);
    }

    // All five Table 2 Health Care categories carry 3% community rate and
    // 0.30 mask efficiency; every other category carries 1% and 0.0.
    for (const auto& category : categories) {
        bool is_expected_health_care = false;
        for (const auto& name : kHealthCareNames) {
            if (category.canonical_name == name) {
                is_expected_health_care = true;
                break;
            }
        }
        check(category.is_health_care_group == is_expected_health_care,
              category.canonical_name + ": is_health_care_group matches Table 2 Health Care group");
        if (is_expected_health_care) {
            check(approx(category.community_infection_rate, 0.03),
                  category.canonical_name + ": health care community rate is 3%");
            check(approx(category.mask_efficiency, 0.30),
                  category.canonical_name + ": health care mask efficiency is 0.30");
        } else {
            check(approx(category.community_infection_rate, 0.01),
                  category.canonical_name + ": non-health-care community rate is 1%");
            check(approx(category.mask_efficiency, 0.0),
                  category.canonical_name + ": non-health-care mask efficiency is 0.0");
        }
    }

    // Museum and Convention share the verified maximum occupancy of 400.
    {
        const std::string museum_name = "Museum";
        const std::string convention_name = "Convention";
        const auto& museum = ecai_cpp2::lookup_occupancy_category(museum_name);
        const auto& convention = ecai_cpp2::lookup_occupancy_category(convention_name);
        check(museum.i0 == 400, "Museum I0 == 400");
        check(convention.i0 == 400, "Convention I0 == 400");
    }

    // Dwelling unit PBR GM uses the publication's Table 3 value, 0.51.
    {
        const std::string dwelling_name = "Dwelling unit";
        const auto& dwelling = ecai_cpp2::lookup_occupancy_category(dwelling_name);
        check(approx(dwelling.pbr_gm_m3h, 0.51),
              "Dwelling unit PBR GM is the publication value 0.51");
    }

    // Canonical normalization: case-insensitive; spaces, commas, and hyphens
    // are equivalent, interchangeable delimiters; no other punctuation is
    // normalized.
    {
        const std::string forms[] = {
            "Food and beverage facilities",
            "food and beverage facilities",
            "FOOD AND BEVERAGE FACILITIES",
            "Food-and-beverage-facilities",
            "Food, and beverage facilities",
            "food   and-beverage,facilities",
            "  Food and beverage facilities  ",
        };
        for (const auto& form : forms) {
            try {
                const auto& category = ecai_cpp2::lookup_occupancy_category(form);
                check(category.canonical_name == "Food and beverage facilities",
                      "normalized form resolves: \"" + form + "\"");
            } catch (const std::exception&) {
                check(false, "normalized form resolves: \"" + form + "\"");
            }
        }

        try {
            const std::string hyphenated_form = "Place-of-religious-worship";
            const auto& category = ecai_cpp2::lookup_occupancy_category(hyphenated_form);
            check(category.canonical_name == "Place of religious worship",
                  "hyphenated multi-word normalized form resolves");
        } catch (const std::exception&) {
            check(false, "hyphenated multi-word normalized form resolves");
        }
    }

    // Every one of the 10 frozen legacy aliases resolves to its correct
    // canonical full record (per D7).
    for (const auto& alias_pair : kAliases) {
        const std::string& alias = alias_pair.first;
        const std::string& expected_canonical = alias_pair.second;
        try {
            const auto& category = ecai_cpp2::lookup_occupancy_category(alias);
            check(category.canonical_name == expected_canonical,
                  "alias \"" + alias + "\" resolves to \"" + expected_canonical + "\"");
        } catch (const std::exception&) {
            check(false, "alias \"" + alias + "\" resolves to \"" + expected_canonical + "\"");
        }
    }

    // Unknown/near-miss names are rejected without fallback -- no fuzzy,
    // prefix, or edit-distance matching.
    {
        const std::string unknown_names[] = {
            "Nonexistent Category",
            "Patient",  // handled above as a valid alias; excluded below
            "Pat",
            "Museum Hall",
            "Musem",
            "Food and beverage facilitie",
            "Food and beverage facilities extra",
            "",
            "Patient room2",
        };
        for (const auto& name : unknown_names) {
            if (name == "Patient") continue;  // valid alias, not a rejection case
            bool threw = false;
            try {
                ecai_cpp2::lookup_occupancy_category(name);
            } catch (const std::invalid_argument&) {
                threw = true;
            } catch (const std::exception&) {
                threw = true;
            }
            check(threw, "unknown/near-miss name rejected: \"" + name + "\"");
        }
    }

    // Unique canonical keys: every canonical name and alias normalizes to a
    // distinct key with no collisions.
    {
        std::vector<std::string> keys;
        for (const auto& category : categories) {
            keys.push_back(ecai_cpp2::normalize_category_key(category.canonical_name));
        }
        for (const auto& alias_pair : kAliases) {
            keys.push_back(ecai_cpp2::normalize_category_key(alias_pair.first));
        }
        bool all_unique = true;
        for (std::size_t i = 0; i < keys.size() && all_unique; ++i) {
            for (std::size_t j = i + 1; j < keys.size(); ++j) {
                if (keys[i] == keys[j]) {
                    all_unique = false;
                    std::cerr << "Collision: \"" << keys[i] << "\" (indices " << i
                              << ", " << j << ")\n";
                    break;
                }
            }
        }
        check(all_unique, "all 25 canonical keys + 10 alias keys are unique");
        check(keys.size() == 35, "35 total keys checked (25 canonical + 10 alias)");
    }

    if (g_failures > 0) {
        std::cerr << g_failures << " test(s) failed.\n";
        return 1;
    }
    std::cout << "All tests passed.\n";
    return 0;
}
