// Work package 1: typed direct-QER family/profile selection and the
// immutable Jones-fitted/Mikszewski parameter tables. See
// include/qer_direct.h for the full contract and provenance.

#include "qer_direct.h"

#include <cmath>
#include <stdexcept>

#include "distributions.h"

namespace ecai_cpp2 {
namespace qer_direct {

namespace {

[[noreturn]] void throw_invalid(const std::string& message) {
    throw std::invalid_argument(message);
}

}  // namespace

std::string to_string(QerFamily family) {
    switch (family) {
        case QerFamily::Jones:
            return "jones";
        case QerFamily::JonesFitted:
            return "jones-fitted";
        case QerFamily::Mikszewski:
            return "mikszewski";
    }
    throw_invalid("to_string(QerFamily): unknown family value");
}

QerFamily parse_qer_family(const std::string& name) {
    if (name == "jones") return QerFamily::Jones;
    if (name == "jones-fitted") return QerFamily::JonesFitted;
    if (name == "mikszewski") return QerFamily::Mikszewski;
    throw_invalid("parse_qer_family: unknown QER family: \"" + name + "\"");
}

// ------------------------------------------------------------------
// Jones-fitted table.
// ------------------------------------------------------------------

const char* const kJonesFittedTableVersion =
    "jones-fitted-v1 (C++/src/qer_fitted.cpp qer_fitted_params, "
    "fit_qer_lognormal.py / qer_lognormal_fit_results.csv)";

namespace {

// The 25 records, in category_data.h publication order. mu_log10/sigma_log10
// are copied exactly from C++/src/qer_fitted.cpp's `qer_fitted_params` map;
// canonical_category is the corresponding category_data.h canonical_name.
const std::array<JonesFittedRecord, kJonesFittedRecordCount>& build_jones_fitted_table() {
    static const std::array<JonesFittedRecord, kJonesFittedRecordCount> table = {{
        {"Cell", {-1.123527, 1.428904}},
        {"Dayroom", {-0.992113, 1.429109}},
        {"Food and beverage facilities", {-1.093306, 1.427494}},
        {"Gym", {-1.298071, 1.432531}},
        {"Office", {-1.123527, 1.428904}},
        {"Retail", {-1.299721, 1.429381}},
        {"Transportation waiting", {-1.026098, 1.429956}},
        {"Classroom", {-1.039436, 1.427156}},
        {"Lecture hall", {-1.349078, 1.428965}},
        {"Manufacturing", {-0.970826, 1.430100}},
        {"Sorting, packing, light assembly", {-0.844935, 1.431221}},
        {"Warehouse", {-0.947311, 1.432394}},
        {"Exam room", {-0.645766, 1.430231}},
        {"Group treatment area", {-1.260304, 1.428842}},
        {"Patient room", {-0.673900, 1.430657}},
        {"Resident room", {-0.827809, 1.430225}},
        {"Waiting room", {-1.299721, 1.429381}},
        {"Auditorium", {-1.349078, 1.428965}},
        {"Place of religious worship", {-1.389396, 1.427938}},
        {"Museum", {-1.399522, 1.429133}},
        {"Convention", {-1.399522, 1.429133}},
        {"Spectator area", {-1.093306, 1.427494}},
        {"Lobbies", {-1.349078, 1.428965}},
        {"Common space", {-1.349078, 1.428965}},
        {"Dwelling unit", {-0.509175, 1.561887}},
    }};
    return table;
}

// The complete legacy-short-to-canonical mapping, exactly 25 entries: one per
// canonical category in build_jones_fitted_table(), in the same order as
// C++/src/qer_fitted.cpp's `qer_fitted_params` map.
const std::array<LegacyJonesFittedAlias, kJonesFittedRecordCount>&
build_jones_fitted_legacy_aliases() {
    static const std::array<LegacyJonesFittedAlias, kJonesFittedRecordCount> aliases = {{
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
    }};
    return aliases;
}

}  // namespace

const std::array<JonesFittedRecord, kJonesFittedRecordCount>& jones_fitted_table() {
    return build_jones_fitted_table();
}

const DirectQerParams& jones_fitted_params_for_canonical_category(
    const std::string& canonical_category_name) {
    for (const auto& record : build_jones_fitted_table()) {
        if (record.canonical_category == canonical_category_name) return record.params;
    }
    throw_invalid(
        "jones_fitted_params_for_canonical_category: unknown canonical category: \"" +
        canonical_category_name + "\"");
}

const std::array<LegacyJonesFittedAlias, kJonesFittedRecordCount>& jones_fitted_legacy_aliases() {
    return build_jones_fitted_legacy_aliases();
}

const DirectQerParams& jones_fitted_params_for_legacy_short_name(
    const std::string& legacy_short_name) {
    for (const auto& alias : build_jones_fitted_legacy_aliases()) {
        if (alias.legacy_short_name == legacy_short_name) {
            return jones_fitted_params_for_canonical_category(alias.canonical_category);
        }
    }
    throw_invalid(
        "jones_fitted_params_for_legacy_short_name: unknown legacy short name: \"" +
        legacy_short_name + "\"");
}

// ------------------------------------------------------------------
// Mikszewski table.
// ------------------------------------------------------------------

const char* const kMikszewskiTableVersion =
    "mikszewski-2022-v1 (Mikszewski et al. 2022 Tables 1-2, via "
    "C++/src/qer_fitted.cpp mikszewski_params)";

namespace {

const char* organism_key(MikszewskiOrganism organism) {
    switch (organism) {
        case MikszewskiOrganism::SarsCov1:
            return "sars-cov-1";
        case MikszewskiOrganism::Mers:
            return "mers";
        case MikszewskiOrganism::TbOnTreatment:
            return "tb-on-treatment";
        case MikszewskiOrganism::Influenza:
            return "influenza";
        case MikszewskiOrganism::Coxsackievirus:
            return "coxsackievirus";
        case MikszewskiOrganism::Rhinovirus:
            return "rhinovirus";
        case MikszewskiOrganism::SarsCov2:
            return "sars-cov-2";
        case MikszewskiOrganism::TbUntreated:
            return "tb-untreated";
        case MikszewskiOrganism::Adenovirus:
            return "adenovirus";
        case MikszewskiOrganism::Measles:
            return "measles";
    }
    throw_invalid("organism_key: unknown MikszewskiOrganism value");
}

const char* activity_key(MikszewskiActivity activity) {
    switch (activity) {
        case MikszewskiActivity::Resting:
            return "resting";
        case MikszewskiActivity::StandingSpeaking:
            return "standing-speaking";
        case MikszewskiActivity::LightSpeakingLoudly:
            return "light-speaking-loudly";
    }
    throw_invalid("activity_key: unknown MikszewskiActivity value");
}

DirectQerParams make_mikszewski_params(double median, double sigma_log10) {
    return DirectQerParams{std::log10(median), sigma_log10};
}

// The 30 records, in C++/src/qer_fitted.cpp `mikszewski_params` map order:
// organism-major, resting/standing-speaking/light-speaking-loudly within
// each organism. Literal medians (quanta/h) and organism-specific
// sigma_log10 are copied exactly from that map; mu_log10 is computed here as
// log10(median), without rounding.
const std::array<MikszewskiProfileRecord, kMikszewskiProfileCount>&
build_mikszewski_profile_table() {
    static const std::array<MikszewskiProfileRecord, kMikszewskiProfileCount> table = {{
        {{MikszewskiOrganism::SarsCov1, MikszewskiActivity::Resting}, 0.0084,
         make_mikszewski_params(0.0084, 1.30)},
        {{MikszewskiOrganism::SarsCov1, MikszewskiActivity::StandingSpeaking}, 0.042,
         make_mikszewski_params(0.042, 1.30)},
        {{MikszewskiOrganism::SarsCov1, MikszewskiActivity::LightSpeakingLoudly}, 0.71,
         make_mikszewski_params(0.71, 1.30)},

        {{MikszewskiOrganism::Mers, MikszewskiActivity::Resting}, 0.011,
         make_mikszewski_params(0.011, 1.60)},
        {{MikszewskiOrganism::Mers, MikszewskiActivity::StandingSpeaking}, 0.056,
         make_mikszewski_params(0.056, 1.60)},
        {{MikszewskiOrganism::Mers, MikszewskiActivity::LightSpeakingLoudly}, 0.96,
         make_mikszewski_params(0.96, 1.60)},

        {{MikszewskiOrganism::TbOnTreatment, MikszewskiActivity::Resting}, 0.020,
         make_mikszewski_params(0.020, 1.40)},
        {{MikszewskiOrganism::TbOnTreatment, MikszewskiActivity::StandingSpeaking}, 0.098,
         make_mikszewski_params(0.098, 1.40)},
        {{MikszewskiOrganism::TbOnTreatment, MikszewskiActivity::LightSpeakingLoudly}, 1.7,
         make_mikszewski_params(1.7, 1.40)},

        {{MikszewskiOrganism::Influenza, MikszewskiActivity::Resting}, 0.035,
         make_mikszewski_params(0.035, 0.84)},
        {{MikszewskiOrganism::Influenza, MikszewskiActivity::StandingSpeaking}, 0.17,
         make_mikszewski_params(0.17, 0.84)},
        {{MikszewskiOrganism::Influenza, MikszewskiActivity::LightSpeakingLoudly}, 3.0,
         make_mikszewski_params(3.0, 0.84)},

        {{MikszewskiOrganism::Coxsackievirus, MikszewskiActivity::Resting}, 0.062,
         make_mikszewski_params(0.062, 1.10)},
        {{MikszewskiOrganism::Coxsackievirus, MikszewskiActivity::StandingSpeaking}, 0.31,
         make_mikszewski_params(0.31, 1.10)},
        {{MikszewskiOrganism::Coxsackievirus, MikszewskiActivity::LightSpeakingLoudly}, 5.2,
         make_mikszewski_params(5.2, 1.10)},

        {{MikszewskiOrganism::Rhinovirus, MikszewskiActivity::Resting}, 0.21,
         make_mikszewski_params(0.21, 0.83)},
        {{MikszewskiOrganism::Rhinovirus, MikszewskiActivity::StandingSpeaking}, 1.0,
         make_mikszewski_params(1.0, 0.83)},
        {{MikszewskiOrganism::Rhinovirus, MikszewskiActivity::LightSpeakingLoudly}, 18.0,
         make_mikszewski_params(18.0, 0.83)},

        {{MikszewskiOrganism::SarsCov2, MikszewskiActivity::Resting}, 0.55,
         make_mikszewski_params(0.55, 1.20)},
        {{MikszewskiOrganism::SarsCov2, MikszewskiActivity::StandingSpeaking}, 2.7,
         make_mikszewski_params(2.7, 1.20)},
        {{MikszewskiOrganism::SarsCov2, MikszewskiActivity::LightSpeakingLoudly}, 46.0,
         make_mikszewski_params(46.0, 1.20)},

        {{MikszewskiOrganism::TbUntreated, MikszewskiActivity::Resting}, 0.62,
         make_mikszewski_params(0.62, 1.30)},
        {{MikszewskiOrganism::TbUntreated, MikszewskiActivity::StandingSpeaking}, 3.1,
         make_mikszewski_params(3.1, 1.30)},
        {{MikszewskiOrganism::TbUntreated, MikszewskiActivity::LightSpeakingLoudly}, 52.0,
         make_mikszewski_params(52.0, 1.30)},

        {{MikszewskiOrganism::Adenovirus, MikszewskiActivity::Resting}, 0.78,
         make_mikszewski_params(0.78, 0.95)},
        {{MikszewskiOrganism::Adenovirus, MikszewskiActivity::StandingSpeaking}, 3.9,
         make_mikszewski_params(3.9, 0.95)},
        {{MikszewskiOrganism::Adenovirus, MikszewskiActivity::LightSpeakingLoudly}, 66.0,
         make_mikszewski_params(66.0, 0.95)},

        {{MikszewskiOrganism::Measles, MikszewskiActivity::Resting}, 3.1,
         make_mikszewski_params(3.1, 1.60)},
        {{MikszewskiOrganism::Measles, MikszewskiActivity::StandingSpeaking}, 15.0,
         make_mikszewski_params(15.0, 1.60)},
        {{MikszewskiOrganism::Measles, MikszewskiActivity::LightSpeakingLoudly}, 260.0,
         make_mikszewski_params(260.0, 1.60)},
    }};
    return table;
}

}  // namespace

const std::array<MikszewskiProfileRecord, kMikszewskiProfileCount>& mikszewski_profile_table() {
    return build_mikszewski_profile_table();
}

std::string canonical_mikszewski_profile_name(MikszewskiProfile profile) {
    return std::string("mikszewski-") + organism_key(profile.organism) + "-" +
           activity_key(profile.activity);
}

MikszewskiProfile parse_mikszewski_profile(const std::string& canonical_name) {
    for (const auto& record : build_mikszewski_profile_table()) {
        if (canonical_mikszewski_profile_name(record.profile) == canonical_name) {
            return record.profile;
        }
    }
    throw_invalid("parse_mikszewski_profile: unknown Mikszewski profile: \"" + canonical_name + "\"");
}

const DirectQerParams& mikszewski_params_for_profile(MikszewskiProfile profile) {
    for (const auto& record : build_mikszewski_profile_table()) {
        if (record.profile.organism == profile.organism &&
            record.profile.activity == profile.activity) {
            return record.params;
        }
    }
    throw_invalid("mikszewski_params_for_profile: profile not present in mikszewski_profile_table()");
}

const DirectQerParams& mikszewski_params_for_canonical_name(const std::string& canonical_name) {
    return mikszewski_params_for_profile(parse_mikszewski_profile(canonical_name));
}

// ------------------------------------------------------------------
// Pure direct-QER log10-normal transform.
// ------------------------------------------------------------------

double sample_direct_qer_log10normal(double u, const DirectQerParams& params) {
    return distributions::log10normal_ppf(u, params.mu_log10, params.sigma_log10);
}

}  // namespace qer_direct
}  // namespace ecai_cpp2
