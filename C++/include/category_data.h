// Stage 2: OccupancyCategory data and strict category lookup.
//
// Source of truth: docs/MODEL_SPECIFICATION.md Tables 1-4 (Jones et al. 2025,
// printed pp. 15, 20, 21, 22), Discrepancy D2 (Dwelling unit PBR GM), and
// Discrepancy D7 (canonical naming / normalization / legacy aliases).
// IMPLEMENTATION_PLAN.md Section 7, Stage 2.

#ifndef ECAI_CPP2_CATEGORY_DATA_H
#define ECAI_CPP2_CATEGORY_DATA_H

#include <array>
#include <string>

namespace ecai_cpp2 {

// One row of Jones et al. Tables 2-4 plus the Table 1 community-rate/mask
// assignment, fully typed for later stages (Stage 5+ distribution sampling,
// Stage 8/9 QER/ECAi calculation). Numeric fields are stored exactly as
// printed in the publication; unit/scale notes are given per field below
// where the printed value is not yet in SI base units.
struct OccupancyCategory {
    std::string canonical_name;  // Publication full name (Tables 2/3/4).
    std::string group;           // Table 2 group heading, e.g. "Health Care".
    bool is_health_care_group;   // Explicit per-category flag (Table 2's
                                  // Health Care group of exactly 5 rows), not
                                  // inferred from `group` at lookup time.

    int i0;             // Max occupants, Table 2.
    double volume_m3;   // Space volume VOL, Table 2, m^3.

    // Table 2 columns not consumed by the QER/ECAi calculation itself;
    // retained here for provenance/completeness per data/reference/README.md.
    double probability_susceptible_present;   // Table 2, dimensionless.
    double people_outdoor_airflow_Lsp;        // Table 2, L/s per person.
    double area_outdoor_airflow_Lsm2;         // Table 2, L/s per m^2.
    double baseline_outdoor_airflow_Lsp;      // Table 4 baseline column,
                                               // L/s per person (RR comparator).

    double community_infection_rate;  // Table 1: 0.03 Health Care, else 0.01.
    double mask_efficiency;           // Sec. 3.3, p.16: 0.30 Health Care, else 0.0.

    double pbr_gm_m3h;  // Table 3, breathing rate geometric mean, m^3/h.
    double pbr_gsd;     // Table 3, breathing rate geometric standard deviation.
    double activity_breathing_pct;   // Table 3 activity split (already folded
    double activity_talking_pct;     // into PBR/C_drop/d GM/GSD per p.19; not
    double activity_vocalizing_pct;  // separately sampled), percent, sums to 100.

    double cdrop_gm_e4_per_m3;  // Table 3, C_drop GM as printed (x10^4 m^-3).
    double cdrop_gsd;           // Table 3, C_drop geometric standard deviation.

    double d_gm_e6_m;  // Table 3, droplet diameter GM as printed (x10^-6 m).
    double d_gsd;       // Table 3, droplet diameter geometric standard deviation.

    double published_ecai_Lsp;  // Table 4, adjusted publication comparator, L/s/person.
    double published_rr;        // Table 4, relative risk.
};

// All 25 Jones et al. occupancy categories, in exact publication order
// (Table 2, p.20).
const std::array<OccupancyCategory, 25>& all_occupancy_categories();

// Strict, case-insensitive lookup by canonical full publication name or by
// one of the 10 frozen legacy short aliases (D7). Runs of spaces, commas, and
// hyphens in `name` are treated as equivalent, interchangeable delimiters; no
// other punctuation is normalized, and no fuzzy/prefix/edit-distance matching
// is performed. Always returns the full canonical record. Throws
// std::invalid_argument if `name` does not exactly match (after
// normalization) any canonical name or enumerated alias.
const OccupancyCategory& lookup_occupancy_category(const std::string& name);

// Exposed for testing: applies the D7 normalization rule (lowercase; runs of
// ' ', ',', '-' collapsed to a single space; leading/trailing delimiters
// trimmed) used internally by lookup_occupancy_category.
std::string normalize_category_key(const std::string& name);

}  // namespace ecai_cpp2

#endif  // ECAI_CPP2_CATEGORY_DATA_H
