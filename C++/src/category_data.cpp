#include "category_data.h"

#include <cctype>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ecai_cpp2 {

namespace {

constexpr double kHealthCareCommunityRate = 0.03;
constexpr double kHealthCareMaskEfficiency = 0.30;
constexpr double kGeneralCommunityRate = 0.01;
constexpr double kGeneralMaskEfficiency = 0.0;

const std::array<OccupancyCategory, 25>& build_categories() {
    static const std::array<OccupancyCategory, 25> categories = {{
        {"Cell", "Correctional Facilities", false,
         30, 320,
         0.25, 2.5, 0.6, 4.9,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.55, 1.2, 80, 20, 0,
         14, 1.1,
         1.9, 1.1,
         15, 0.44},
        {"Dayroom", "Correctional Facilities", false,
         30, 600,
         0.25, 2.5, 0.3, 3.5,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.55, 1.2, 75, 25, 0,
         19, 1.1,
         1.9, 1.1,
         20, 0.35},
        {"Food and beverage facilities", "Commercial/Retail", false,
         50, 360,
         0.38, 3.8, 0.9, 5.6,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.55, 1.2, 75, 25, 0,
         15, 1.1,
         1.9, 1.1,
         30, 0.24},
        {"Gym", "Commercial/Retail", false,
         180, 1900,
         0.83, 10, 0.3, 11,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.62, 1.3, 100, 0, 0,
         9.8, 1.1,
         1.8, 1.1,
         40, 0.32},
        {"Office", "Commercial/Retail", false,
         50, 2700,
         0.38, 2.5, 0.3, 8.5,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.55, 1.2, 80, 20, 0,
         14, 1.1,
         1.9, 1.1,
         15, 0.78},
        {"Retail", "Commercial/Retail", false,
         150, 6000,
         0.77, 3.8, 0.6, 7.8,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.55, 1.2, 95, 5, 0,
         11, 1.1,
         1.8, 1.1,
         20, 0.61},
        {"Transportation waiting", "Commercial/Retail", false,
         100, 10000,
         0.62, 3.8, 0.3, 6.8,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.64, 1.2, 75, 25, 0,
         15, 1.1,
         1.9, 1.1,
         30, 0.60},
        {"Classroom", "Educational Facilities", false,
         30, 320,
         0.25, 5.0, 0.6, 7.4,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.55, 1.2, 70, 30, 0,
         17, 1.1,
         1.9, 1.1,
         20, 0.45},
        {"Lecture hall", "Educational Facilities", false,
         150, 1800,
         0.77, 2.5, 0.3, 2.9,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.55, 1.2, 100, 0, 0,
         9.8, 1.1,
         1.8, 1.1,
         25, 0.22},
        {"Manufacturing", "Industrial", false,
         70, 12000,
         0.50, 5.0, 0.3, 9.3,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.73, 1.3, 75, 25, 0,
         15, 1.1,
         1.9, 1.1,
         25, 0.78},
        {"Sorting, packing, light assembly", "Industrial", false,
         20, 4800,
         0.17, 5.0, 0.3, 11,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.77, 1.3, 75, 25, 0,
         19, 1.1,
         1.9, 1.1,
         10, 1.0},
        {"Warehouse", "Industrial", false,
         20, 1200,
         0.17, 5.0, 0.3, 6.5,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.77, 1.3, 75, 25, 0,
         15, 1.1,
         1.9, 1.1,
         10, 0.87},
        {"Exam room", "Health Care", true,
         3, 41,
         0.057, 7.5, 0.0, 7.5,
         kHealthCareCommunityRate, kHealthCareMaskEfficiency,
         0.62, 1.3, 100, 0, 0,
         32, 1.1,
         2.0, 1.1,
         20, 0.47},
        {"Group treatment area", "Health Care", true,
         20, 270,
         0.43, 7.5, 0.0, 7.5,
         kHealthCareCommunityRate, kHealthCareMaskEfficiency,
         0.55, 1.2, 90, 10, 0,
         12, 1.1,
         1.8, 1.1,
         35, 0.29},
        {"Patient room", "Health Care", true,
         3, 81,
         0.057, 15, 0.0, 15,
         kHealthCareCommunityRate, kHealthCareMaskEfficiency,
         0.62, 1.3, 100, 0, 0,
         30, 1.1,
         2.0, 1.1,
         35, 0.53},
        {"Resident room", "Health Care", true,
         3, 81,
         0.057, 15, 0.0, 15,
         kHealthCareCommunityRate, kHealthCareMaskEfficiency,
         0.62, 1.3, 100, 0, 0,
         21, 1.1,
         2.0, 1.1,
         25, 0.69},
        {"Waiting room", "Health Care", true,
         30, 270,
         0.57, 5.0, 0.0, 5.0,
         kHealthCareCommunityRate, kHealthCareMaskEfficiency,
         0.55, 1.2, 90, 10, 0,
         11, 1.1,
         1.8, 1.1,
         45, 0.16},
        {"Auditorium", "Public Assembly/Sports/Entertainment", false,
         150, 4600,
         0.77, 2.5, 0.3, 3.5,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.55, 1.2, 90, 10, 0,
         9.8, 1.1,
         1.8, 1.1,
         25, 0.36},
        {"Place of religious worship", "Public Assembly/Sports/Entertainment", false,
         180, 2400,
         0.83, 2.5, 0.3, 2.8,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.50, 1.2, 100, 0, 0,
         9.8, 1.1,
         1.8, 1.1,
         25, 0.23},
        {"Museum", "Public Assembly/Sports/Entertainment", false,
         400, 10000,
         0.97, 3.8, 0.3, 4.6,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.49, 1.2, 100, 0, 0,
         9.8, 1.1,
         1.8, 1.1,
         30, 0.31},
        {"Convention", "Public Assembly/Sports/Entertainment", false,
         400, 10000,
         0.97, 3.8, 0.3, 4.6,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.49, 1.2, 100, 0, 0,
         9.8, 1.1,
         1.8, 1.1,
         30, 0.31},
        {"Spectator area", "Public Assembly/Sports/Entertainment", false,
         100, 6300,
         0.62, 3.8, 0.3, 4.7,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.55, 1.2, 70, 30, 0,
         15, 1.1,
         1.9, 1.1,
         25, 0.52},
        {"Lobbies", "Public Assembly/Sports/Entertainment", false,
         150, 4600,
         0.77, 2.5, 0.3, 3.5,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.55, 1.2, 90, 10, 0,
         9.8, 1.1,
         1.8, 1.1,
         25, 0.34},
        {"Common space", "Residential", false,
         150, 4600,
         0.77, 2.5, 0.3, 1.0,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.55, 1.2, 90, 10, 0,
         9.8, 1.1,
         1.8, 1.1,
         25, 0.26},
        {"Dwelling unit", "Residential", false,
         6, 540,
         0.049, 6.8, 0.0, 6.8,
         kGeneralCommunityRate, kGeneralMaskEfficiency,
         0.51, 4.3, 15, 70, 15,
         40, 1.1,
         2.2, 1.1,
         15, 0.80},
    }};
    return categories;
}

struct AliasEntry {
    const char* alias;
    const char* canonical_name;
};

// The 10 frozen legacy short aliases (D7), matching the existing
// C++/Python `occupancy_params` map's shortened keys.
constexpr AliasEntry kLegacyAliases[10] = {
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

bool is_delimiter(char c) {
    return c == ' ' || c == ',' || c == '-';
}

const std::unordered_map<std::string, std::size_t>& build_lookup_index() {
    static const std::unordered_map<std::string, std::size_t> index = [] {
        std::unordered_map<std::string, std::size_t> map;
        const auto& categories = build_categories();
        for (std::size_t i = 0; i < categories.size(); ++i) {
            map.emplace(normalize_category_key(categories[i].canonical_name), i);
        }
        for (const auto& alias : kLegacyAliases) {
            for (std::size_t i = 0; i < categories.size(); ++i) {
                if (categories[i].canonical_name == alias.canonical_name) {
                    map.emplace(normalize_category_key(alias.alias), i);
                    break;
                }
            }
        }
        return map;
    }();
    return index;
}

}  // namespace

std::string normalize_category_key(const std::string& name) {
    std::string out;
    out.reserve(name.size());
    bool pending_delimiter = false;
    bool started = false;
    for (unsigned char c : name) {
        if (is_delimiter(static_cast<char>(c))) {
            if (started) {
                pending_delimiter = true;
            }
            continue;
        }
        if (pending_delimiter) {
            out.push_back(' ');
            pending_delimiter = false;
        }
        out.push_back(static_cast<char>(std::tolower(c)));
        started = true;
    }
    return out;
}

const std::array<OccupancyCategory, 25>& all_occupancy_categories() {
    return build_categories();
}

const OccupancyCategory& lookup_occupancy_category(const std::string& name) {
    const auto& index = build_lookup_index();
    const auto it = index.find(normalize_category_key(name));
    if (it == index.end()) {
        throw std::invalid_argument(
            "lookup_occupancy_category: unknown occupancy category name: \"" +
            name + "\"");
    }
    return build_categories()[it->second];
}

}  // namespace ecai_cpp2
