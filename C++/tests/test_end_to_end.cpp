// Stage 14 bounded deterministic regression, wiring, and accounting checks.
#include "model.h"

#include <array>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/wait.h>
#include <vector>

#ifndef ECAI_CPP2_SOURCE_DIR
#error "ECAI_CPP2_SOURCE_DIR must be defined by CMake"
#endif
#ifndef ECAI_CPP2_EXECUTABLE
#error "ECAI_CPP2_EXECUTABLE must be defined by CMake"
#endif

namespace fs = std::filesystem;
namespace {
int checks = 0;
int failures = 0;

void check(bool condition, const std::string& message) {
    ++checks;
    if (condition) std::cout << "PASS: " << message << '\n';
    else { ++failures; std::cerr << "FAIL: " << message << '\n'; }
}

bool approx(double actual, double expected, double abs_tol = 1e-12,
            double rel_tol = 1e-10) {
    if (!std::isfinite(actual) || !std::isfinite(expected)) return false;
    const double difference = std::fabs(actual - expected);
    return difference <= abs_tol || difference <= rel_tol * std::fabs(expected);
}

std::string read_text(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot read " + path.string());
    std::ostringstream text;
    text << input.rdbuf();
    return text.str();
}

std::string collapse_whitespace(const std::string& text) {
    std::string result;
    bool pending_space = false;
    for (const unsigned char ch : text) {
        if (std::isspace(ch)) {
            pending_space = !result.empty();
        } else {
            if (pending_space) result.push_back(' ');
            result.push_back(static_cast<char>(ch));
            pending_space = false;
        }
    }
    return result;
}

std::vector<std::vector<std::string>> parse_csv(const std::string& text) {
    std::vector<std::vector<std::string>> rows(1);
    std::string field;
    bool quoted = false;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char ch = text[i];
        if (quoted) {
            if (ch == '"') {
                if (i + 1 < text.size() && text[i + 1] == '"') {
                    field.push_back('"');
                    ++i;
                } else {
                    quoted = false;
                }
            } else {
                field.push_back(ch);
            }
        } else if (ch == '"' && field.empty()) {
            quoted = true;
        } else if (ch == ',') {
            rows.back().push_back(field);
            field.clear();
        } else if (ch == '\n') {
            if (!field.empty() && field.back() == '\r') field.pop_back();
            rows.back().push_back(field);
            field.clear();
            rows.emplace_back();
        } else {
            field.push_back(ch);
        }
    }
    if (quoted) throw std::runtime_error("unterminated quoted CSV field");
    if (!field.empty() || !rows.back().empty()) rows.back().push_back(field);
    if (!rows.empty() && rows.back().empty()) rows.pop_back();
    return rows;
}

std::map<std::string, std::string> named_row(
    const std::vector<std::string>& header, const std::vector<std::string>& row) {
    if (header.size() != row.size()) throw std::runtime_error("CSV row width mismatch");
    std::map<std::string, std::string> names;
    for (std::size_t i = 0; i < header.size(); ++i) names.emplace(header[i], row[i]);
    return names;
}

double finite_number(const std::string& text) {
    std::size_t used = 0;
    const double value = std::stod(text, &used);
    if (used != text.size() || !std::isfinite(value))
        throw std::runtime_error("invalid finite number: " + text);
    return value;
}

bool bounded_reference_matches(const std::string& field, const std::string& actual,
                               const std::string& expected) {
    if (field != "raw_mean_l_s_person" && field != "raw_ci_lower_l_s_person" &&
        field != "raw_ci_upper_l_s_person") {
        return actual == expected;
    }
    try {
        return std::fabs(finite_number(actual) - finite_number(expected)) <= 1e-9;
    } catch (const std::exception&) {
        return false;
    }
}

std::string shell_quote(const std::string& text) {
    std::string result = "'";
    for (const char ch : text) result += ch == '\'' ? "'\\''" : std::string(1, ch);
    result += "'";
    return result;
}

struct CommandResult {
    int exit_code = -1;
    std::string output;
};

CommandResult run_command(const std::string& arguments) {
    const std::string command = shell_quote(ECAI_CPP2_EXECUTABLE) + " " + arguments + " 2>&1";
    FILE* pipe = popen(command.c_str(), "r");
    if (pipe == nullptr) throw std::runtime_error("failed to launch executable");
    std::array<char, 4096> buffer{};
    std::string output;
    std::size_t count = 0;
    while ((count = std::fread(buffer.data(), 1, buffer.size(), pipe)) != 0)
        output.append(buffer.data(), count);
    const int status = pclose(pipe);
    CommandResult result;
    result.output = output;
    if (WIFEXITED(status)) result.exit_code = WEXITSTATUS(status);
    return result;
}

CommandResult run_to_directory(const std::string& category, std::size_t n, std::size_t r,
                               const fs::path& directory, bool write_replicates = false,
                               const std::string& extra = "") {
    std::ostringstream args;
    args << "--category " << shell_quote(category)
         << " --simulations " << n << " --replicates " << r
         << " --output-directory " << shell_quote(directory.string());
    if (write_replicates) args << " --write-replicates";
    if (!extra.empty()) args << ' ' << extra;
    return run_command(args.str());
}

struct TempTree {
    fs::path path;
    ~TempTree() {
        std::error_code ignored;
        fs::remove_all(path, ignored);
    }
};

void check_documentation(const fs::path& source) {
    const std::string validation = read_text(source / "docs/VALIDATION.md");
    const std::string semantics = collapse_whitespace(validation);
    check(validation.find("doi:10.1016/j.buildenv.2025.113318") != std::string::npos,
          "validation report cites the publication DOI");
    check(validation.find("printed pp. 10 and 19") != std::string::npos &&
              validation.find("printed p. 12") != std::string::npos &&
              validation.find("p. 22") != std::string::npos,
          "validation report cites printed source pages for PBR, LHS, and Table 4");
    check(validation.find("Supplementary Materials") != std::string::npos &&
              validation.find("unavailable") != std::string::npos &&
              validation.find("not used by the Project Committee") != std::string::npos,
          "validation report qualifies the unavailable supplementary CDF");
    check(validation.find("implementation/invariant/regression gates: **PASS**") !=
              std::string::npos &&
              validation.find("Patient publication comparator: `FAIL_GT_ONE_INCREMENT`") !=
              std::string::npos &&
              validation.find("Museum publication comparator: `NEAR_WITHIN_ONE_INCREMENT`") !=
              std::string::npos &&
              validation.find("overall publication reproduction: **NOT PASS**") !=
              std::string::npos,
          "validation status separates deterministic gates from failed publication reproduction");
    check(validation.find("required investigations are complete") != std::string::npos &&
              validation.find("no production bug") != std::string::npos &&
              validation.find("frozen alternative design") != std::string::npos &&
              validation.find("implementation may advance") != std::string::npos &&
              validation.find("comparator remains explicitly failed") != std::string::npos,
          "validation report records the plan 864-867 investigation disposition");
    check(semantics.find("Only the Table 4 values and source parameters") != std::string::npos &&
              semantics.find("deterministic C++ implementation baselines") != std::string::npos &&
              semantics.find("not independent scientific-correctness evidence") !=
              std::string::npos,
          "validation report distinguishes publication references from regression baselines");
    check(validation.find("person-specific separate PBR") != std::string::npos &&
              validation.find("explicit identical PBR") != std::string::npos &&
              validation.find("matrix design is underspecified") != std::string::npos &&
              validation.find("not demonstrated equivalent") != std::string::npos,
          "validation report limits the architecture conclusion to source-supported distinctions");
    check(semantics.find("approximate randomized-LHS Monte Carlo Student-t CI for the mean") !=
              std::string::npos &&
              semantics.find("not a CI for a population percentile, individual variation, or total model uncertainty") !=
              std::string::npos &&
              semantics.find("publication supplies no corresponding CI") != std::string::npos,
          "validation report defines and limits the reported Monte Carlo CI");
    check(validation.find("all_replicates_zero_infector_count") != std::string::npos &&
              validation.find("all_replicates_negative_clamp_count") != std::string::npos &&
              validation.find("all_replicates_output_zero_count") != std::string::npos &&
              validation.find("aggregate metadata from the same bounded runs") != std::string::npos,
          "validation report makes zero-mass claims traceable to retained aggregate counters");
    check(validation.find("Empirical percentile of published adjusted comparator") ==
              std::string::npos &&
              validation.find("Conditional-active p96") == std::string::npos &&
              validation.find("Shared-PBR diagnostic") == std::string::npos &&
              validation.find("Legacy family-buffer diagnostic") == std::string::npos &&
              validation.find("positive/active subset") == std::string::npos,
          "validation report excludes diagnostics not reproducible from retained artifacts");
    check(validation.find("10,000 simulations × 500 replicates was **not run**") !=
              std::string::npos && validation.find("pending explicit user authorization") !=
              std::string::npos,
          "validation report leaves the full default run unauthorized and pending");
}

}  // namespace

int main() {
    using ecai_cpp2::model::ECAiInputs;
    using ecai_cpp2::sampling_plan::InfectorParameters;
    const fs::path source(ECAI_CPP2_SOURCE_DIR);

    check(bounded_reference_matches("raw_mean_l_s_person", "1.0000000005", "1"),
          "bounded floating-point regression accepts a difference below 1e-9");
    check(!bounded_reference_matches("raw_mean_l_s_person", "1.000000002", "1"),
          "bounded floating-point regression rejects a difference above 1e-9");
    check(!bounded_reference_matches("raw_mean_l_s_person", "nan", "1"),
          "bounded floating-point regression rejects non-finite values");
    check(!bounded_reference_matches("all_replicates_n", "100001", "100000"),
          "bounded regression still compares integer counters exactly");

    // Hard-coded deterministic algebra regression fixtures through production pure functions.
    InfectorParameters infector;
    infector.pbr_m3h = 2.0;
    infector.c_drop = 10.0;
    infector.d_bar = 2.0;
    infector.evaporation_factor = 3.0;
    infector.gvl_copies_per_ml = 1.0;
    infector.viable_fraction = 0.01;
    infector.rtd = 0.5;
    infector.dose_constant = 10.0;
    const auto qer = ecai_cpp2::model::compute_qer(infector);
    check(approx(qer.vdrop, 1.1309733552923255e-11), "tiny fixture Vdrop");
    check(approx(qer.gvl_copies_per_m3, 1000000.0), "tiny fixture GVL");
    check(approx(qer.ver_virions_per_h, 2.2619467105846509e-7), "tiny fixture VER");
    check(approx(qer.qer_per_h, 1.1309733552923255e-8), "tiny fixture QER");

    ECAiInputs inputs;
    inputs.susceptible_pbr_m3h = 0.6;
    inputs.gamma_per_h = 0.5;
    inputs.lambda_per_h = 0.5;
    inputs.volume_m3 = 50.0;
    inputs.occupants = 10;
    inputs.mask_efficiency = 0.0;
    inputs.duration_h = 1.0;
    inputs.target_probability = 0.5;
    const auto ecai = ecai_cpp2::model::compute_ecai(inputs, std::vector<double>{100.0});
    check(approx(ecai.ecai_l_s_person, 1.01560284592605), "tiny fixture ECAi");

    const std::vector<std::string> reference_header = {
        "source_citation", "category", "source_i0", "source_volume_m3",
        "source_community_rate", "source_mask_efficiency", "effective_i0",
        "effective_volume_m3", "effective_community_rate", "effective_mask_efficiency",
        "simulations_n", "replicates_r", "base_seed", "percentile", "confidence_level",
        "target_probability", "duration_hours", "all_replicates_n",
        "all_replicates_zero_infector_count", "all_replicates_negative_clamp_count",
        "all_replicates_output_zero_count", "raw_mean_l_s_person",
        "raw_ci_lower_l_s_person", "raw_ci_upper_l_s_person", "ceiling_to_5_l_s_person",
        "published_adjusted_l_s_person", "signed_difference_l_s_person",
        "absolute_difference_l_s_person", "comparator_status", "sampling_design_version",
        "qualification"};
    const auto reference = parse_csv(read_text(source / "data/reference/publication_cases.csv"));
    check(reference.size() == 3, "reference CSV has exactly one header and two data rows");
    check(!reference.empty() && reference[0] == reference_header,
          "reference CSV header and order are exact");
    if (reference.size() != 3 || reference[0] != reference_header) {
        std::cerr << failures << " of " << checks << " checks failed\n";
        return 1;
    }

    std::vector<std::map<std::string, std::string>> cases;
    for (std::size_t row_index = 1; row_index < reference.size(); ++row_index) {
        check(reference[row_index].size() == reference_header.size(),
              "reference row has exact width " + std::to_string(row_index));
        if (reference[row_index].size() != reference_header.size()) continue;
        cases.push_back(named_row(reference_header, reference[row_index]));
    }
    check(cases.size() == 2 && cases[0].at("category") == "Patient room" &&
              cases[1].at("category") == "Museum",
          "reference rows use the two canonical categories in frozen order");

    const std::vector<std::string> numeric_reference_fields = {
        "source_i0", "source_volume_m3", "source_community_rate", "source_mask_efficiency",
        "effective_i0", "effective_volume_m3", "effective_community_rate",
        "effective_mask_efficiency", "simulations_n", "replicates_r", "base_seed",
        "percentile", "confidence_level", "target_probability", "duration_hours",
        "all_replicates_n", "all_replicates_zero_infector_count",
        "all_replicates_negative_clamp_count", "all_replicates_output_zero_count",
        "raw_mean_l_s_person", "raw_ci_lower_l_s_person", "raw_ci_upper_l_s_person",
        "ceiling_to_5_l_s_person", "published_adjusted_l_s_person",
        "signed_difference_l_s_person", "absolute_difference_l_s_person",
        "sampling_design_version"};
    for (const auto& item : cases) {
        bool finite = true;
        try {
            for (const auto& field : numeric_reference_fields) (void)finite_number(item.at(field));
        } catch (const std::exception&) {
            finite = false;
        }
        check(finite, item.at("category") + " reference numeric fields are finite");
        check(item.at("source_i0") == item.at("effective_i0") &&
                  item.at("source_volume_m3") == item.at("effective_volume_m3") &&
                  item.at("source_community_rate") == item.at("effective_community_rate") &&
                  item.at("source_mask_efficiency") == item.at("effective_mask_efficiency"),
              item.at("category") + " retains unmodified source conditions");
        const double mean = finite_number(item.at("raw_mean_l_s_person"));
        const double ceiling = finite_number(item.at("ceiling_to_5_l_s_person"));
        const double published = finite_number(item.at("published_adjusted_l_s_person"));
        const double signed_difference = finite_number(item.at("signed_difference_l_s_person"));
        const double absolute_difference = finite_number(item.at("absolute_difference_l_s_person"));
        check(ceiling == 5.0 * std::ceil(mean / 5.0),
              item.at("category") + " ceiling-to-5 formula is exact");
        check(signed_difference == ceiling - published &&
                  absolute_difference == std::fabs(signed_difference),
              item.at("category") + " signed and absolute differences reconcile");
        const std::string expected_status = absolute_difference == 0.0 ? "EXACT_MATCH" :
            (absolute_difference == 5.0 ? "NEAR_WITHIN_ONE_INCREMENT" : "FAIL_GT_ONE_INCREMENT");
        check(item.at("comparator_status") == expected_status,
              item.at("category") + " comparison policy classification is exact");
        check(item.at("qualification") ==
                  "BOUNDED_CPP2_DIAGNOSTIC_NOT_PUBLICATION_EQUIVALENT",
              item.at("category") + " is explicitly qualified as non-publication-equivalent");
        const double all_n = finite_number(item.at("all_replicates_n"));
        const double zero_infector = finite_number(item.at("all_replicates_zero_infector_count"));
        const double negative_clamp = finite_number(item.at("all_replicates_negative_clamp_count"));
        const double output_zero = finite_number(item.at("all_replicates_output_zero_count"));
        check(all_n == finite_number(item.at("simulations_n")) *
                           finite_number(item.at("replicates_r")),
              item.at("category") + " retained aggregate row count equals N times R");
        check(zero_infector >= 0.0 && negative_clamp >= 0.0 && output_zero >= 0.0 &&
                  zero_infector <= all_n && negative_clamp <= all_n && output_zero <= all_n,
              item.at("category") + " retained aggregate counters are bounded by all rows");
        check(output_zero == zero_infector + negative_clamp,
              item.at("category") + " retained zero outputs reconcile to zero-infector plus clamp counts");
    }
    check(cases[0].at("source_i0") == "3" && cases[0].at("source_volume_m3") == "81" &&
              cases[0].at("source_community_rate") == "0.029999999999999999" &&
              cases[0].at("source_mask_efficiency") == "0.29999999999999999" &&
              cases[1].at("source_i0") == "400" &&
              cases[1].at("source_volume_m3") == "10000" &&
              cases[1].at("source_community_rate") == "0.01" &&
              cases[1].at("source_mask_efficiency") == "0",
          "Patient and Museum source conditions are frozen explicitly");

    check_documentation(source);

    TempTree temporary{fs::temp_directory_path() /
        ("ecai_cpp2_stage14_" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()))};
    std::error_code error;
    fs::remove_all(temporary.path, error);

    // Small N/R path: executable, simulation, summary writer and replicate writer.
    const fs::path small_a = temporary.path / "small-a";
    const fs::path small_b = temporary.path / "small-b";
    const CommandResult small_result = run_to_directory("Exam room", 17, 3, small_a, true);
    check(small_result.exit_code == 0, "small N/R full executable path succeeds");
    const auto small_summary = parse_csv(read_text(small_a / "summary.csv"));
    const auto small_replicates = parse_csv(read_text(small_a / "replicates.csv"));
    check(small_summary.size() == 2 && small_summary[0].size() == 89 &&
              small_summary[1].size() == 89,
          "small case summary has one row and exact width 89");
    bool replicate_widths = small_replicates.size() == 4;
    for (const auto& row : small_replicates) replicate_widths = replicate_widths && row.size() == 19;
    check(replicate_widths, "small case replicate CSV has exact width 19 and R rows");
    if (small_summary.size() == 2 && small_summary[0].size() == 89 &&
        small_summary[1].size() == 89) {
        const auto names = named_row(small_summary[0], small_summary[1]);
        check(names.at("category") == "Exam room" && names.at("simulations_n") == "17" &&
                  names.at("replicates_r") == "3" && names.at("base_seed") == "485",
              "small case metadata identifies canonical category and applied N/R/seed");
        check(names.at("all_replicates_n") == "51" && names.at("counters_reconcile") == "true" &&
                  names.at("all_cursors_exhausted") == "true" &&
                  names.at("all_replicates_finite") == "true" &&
                  names.at("banks_reconcile") == "true",
              "small case accounting, cursor, finiteness, and bank invariants pass");
    }
    const CommandResult repeat_result = run_to_directory("Exam room", 17, 3, small_b, true);
    check(repeat_result.exit_code == 0 &&
              read_text(small_a / "summary.csv") == read_text(small_b / "summary.csv") &&
              read_text(small_a / "replicates.csv") == read_text(small_b / "replicates.csv"),
          "same toolchain/configuration/seed produces byte-identical CSV files");

    const fs::path custom = temporary.path / "custom";
    const CommandResult custom_result = run_to_directory(
        "Exam room", 17, 3, custom, false, "--volume-m3 50");
    const auto custom_csv = parse_csv(read_text(custom / "summary.csv"));
    bool custom_suppressed = custom_result.exit_code == 0 && custom_csv.size() == 2;
    if (custom_suppressed) {
        const auto names = named_row(custom_csv[0], custom_csv[1]);
        custom_suppressed = names.at("custom_scenario") == "true" &&
            names.at("comparator_eligible") == "false" &&
            names.at("published_comparator_l_s_person").empty() &&
            names.at("comparator_signed_difference_l_s_person").empty() &&
            names.at("comparator_abs_difference_l_s_person").empty() &&
            names.at("comparator_status") == "NOT_APPLICABLE_CUSTOM";
    }
    check(custom_suppressed, "custom condition override suppresses publication comparator");

    // Publication-condition bounded regressions: only R is bounded; model defaults are unchanged.
    for (std::size_t index = 0; index < cases.size(); ++index) {
        const auto& item = cases[index];
        const std::size_t n = static_cast<std::size_t>(finite_number(item.at("simulations_n")));
        const std::size_t r = static_cast<std::size_t>(finite_number(item.at("replicates_r")));
        const fs::path output_directory = temporary.path / (index == 0 ? "patient" : "museum");
        const CommandResult result = run_to_directory(item.at("category"), n, r, output_directory);
        check(result.exit_code == 0, item.at("category") + " bounded executable run succeeds");
        const auto generated = parse_csv(read_text(output_directory / "summary.csv"));
        bool reference_match = generated.size() == 2 && generated[0].size() == 89 &&
                           generated[1].size() == 89;
        if (reference_match) {
            const auto names = named_row(generated[0], generated[1]);
            const std::vector<std::pair<std::string, std::string>> mappings = {
                {"category", "category"}, {"source_i0", "source_i0"},
                {"source_volume_m3", "source_volume_m3"},
                {"source_community_rate", "source_community_rate"},
                {"source_mask_efficiency", "source_mask_efficiency"},
                {"effective_i0", "effective_i0"}, {"effective_volume_m3", "effective_volume_m3"},
                {"effective_community_rate", "effective_community_rate"},
                {"effective_mask_efficiency", "effective_mask_efficiency"},
                {"simulations_n", "simulations_n"}, {"replicates_r", "replicates_r"},
                {"base_seed", "base_seed"}, {"percentile", "percentile"},
                {"confidence_level", "confidence_level"},
                {"target_probability", "target_probability"}, {"duration_hours", "duration_hours"},
                {"all_replicates_n", "all_replicates_n"},
                {"all_replicates_zero_infector_count", "all_replicates_zero_infector_count"},
                {"all_replicates_negative_clamp_count", "all_replicates_negative_clamp_count"},
                {"all_replicates_output_zero_count", "all_replicates_output_zero_count"},
                {"raw_mean_l_s_person", "raw_mean_l_s_person"},
                {"raw_ci_lower_l_s_person", "raw_ci_lower_l_s_person"},
                {"raw_ci_upper_l_s_person", "raw_ci_upper_l_s_person"},
                {"ceiling_to_5_l_s_person", "ceiling_to_5_l_s_person"},
                {"published_comparator_l_s_person", "published_adjusted_l_s_person"},
                {"comparator_signed_difference_l_s_person", "signed_difference_l_s_person"},
                {"comparator_abs_difference_l_s_person", "absolute_difference_l_s_person"},
                {"comparator_status", "comparator_status"},
                {"sampling_design_version", "sampling_design_version"}};
            for (const auto& mapping : mappings)
                reference_match = reference_match && bounded_reference_matches(
                    mapping.first, names.at(mapping.first), item.at(mapping.second));
            reference_match = reference_match && names.at("custom_scenario") == "false" &&
                names.at("comparator_eligible") == "true" &&
                names.at("counters_reconcile") == "true" &&
                names.at("all_cursors_exhausted") == "true" &&
                names.at("all_replicates_finite") == "true" &&
                names.at("banks_reconcile") == "true";
        }
        check(reference_match, item.at("category") +
              " generated output matches the bounded regression and exact accounting reference");
    }

    fs::remove_all(temporary.path, error);
    check(!error && !fs::exists(temporary.path), "all Stage 14 generated output is removed");

    if (failures != 0) {
        std::cerr << failures << " of " << checks << " checks failed\n";
        return 1;
    }
    std::cout << "All " << checks
              << " Stage 14 deterministic regression, wiring, and accounting checks passed.\n";
    return 0;
}
