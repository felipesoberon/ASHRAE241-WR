#include "output.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <random>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <utility>

#include "model.h"
#include "seed_manager.h"

#ifndef ECAI_CPP2_PROTECTED_LEGACY_TREE
#define ECAI_CPP2_PROTECTED_LEGACY_TREE ""
#endif

namespace fs = std::filesystem;

namespace ecai_cpp2::output {
namespace {

detail::WriteTestHook g_write_test_hook;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::invalid_argument("output: " + message);
}

void require_finite(double value, const std::string& name) {
    require(std::isfinite(value), "non-finite " + name);
}

bool same_double(double left, double right) {
    return std::memcmp(&left, &right, sizeof(double)) == 0;
}

bool same_optional(const std::optional<double>& left, const std::optional<double>& right) {
    if (left.has_value() != right.has_value()) return false;
    return !left || same_double(*left, *right);
}

std::string format_double(double value) {
    require_finite(value, "numeric value");
    std::ostringstream text;
    text.imbue(std::locale::classic());
    text << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return text.str();
}

std::string format_optional(const std::optional<double>& value) {
    return value ? format_double(*value) : "";
}

std::string boolean(bool value) { return value ? "true" : "false"; }

std::string quote_csv_field(const std::string& field) {
    if (field.find_first_of(",\"\r\n") == std::string::npos) return field;
    std::string quoted = "\"";
    for (const char ch : field) quoted += ch == '"' ? "\"\"" : std::string(1, ch);
    quoted += '"';
    return quoted;
}

void validate_document(const OutputDocument& document) {
    require(document.summary_rows.empty() || !document.summary_columns.empty(),
            "summary has rows but no schema");
    require(document.replicate_rows.empty() || !document.replicate_columns.empty(),
            "replicates have rows but no schema");
    for (const auto& row : document.summary_rows)
        require(row.size() == document.summary_columns.size(), "summary row width mismatch");
    for (const auto& row : document.replicate_rows)
        require(row.size() == document.replicate_columns.size(), "replicate row width mismatch");
}

std::size_t parse_size(const std::string& text, const char* name) {
    std::size_t value = 0;
    const char* const begin = text.data();
    const char* const end = begin + text.size();
    const auto parsed = std::from_chars(begin, end, value);
    require(begin != end && parsed.ec == std::errc{} && parsed.ptr == end,
            std::string("invalid ") + name + " metadata");
    return value;
}

std::string serialize(const std::vector<std::string>& header,
                      const std::vector<std::vector<std::string>>& rows) {
    require(!header.empty(), "empty CSV schema");
    std::ostringstream output;
    auto write_row = [&](const std::vector<std::string>& row) {
        require(row.size() == header.size(), "CSV row width mismatch");
        for (std::size_t i = 0; i < row.size(); ++i) {
            if (i != 0) output << ',';
            output << quote_csv_field(row[i]);
        }
        output << '\n';
    };
    write_row(header);
    for (const auto& row : rows) write_row(row);
    return output.str();
}

std::size_t checked_add(std::size_t left, std::size_t right, const char* name) {
    require(right <= std::numeric_limits<std::size_t>::max() - left,
            std::string("counter overflow: ") + name);
    return left + right;
}

std::size_t checked_multiply(std::size_t left, std::size_t right, const char* name) {
    require(left == 0 || right <= std::numeric_limits<std::size_t>::max() / left,
            std::string("counter overflow: ") + name);
    return left * right;
}

std::string absolute_path(const std::string& raw) {
    require(!raw.empty(), "output path is empty");
    for (const unsigned char ch : raw)
        require(ch >= 32 && ch != 127, "output path contains control characters");
    std::error_code error;
    const fs::path path = fs::absolute(fs::path(raw), error);
    require(!error, "cannot make output path absolute");
    return path.lexically_normal().string();
}

bool same_or_child(const fs::path& child, const fs::path& root) {
    std::error_code error;
    const fs::path canonical_child = fs::weakly_canonical(child, error);
    require(!error, "cannot resolve output destination");
    const fs::path canonical_root = fs::weakly_canonical(root, error);
    require(!error, "cannot resolve protected tree");
    auto child_part = canonical_child.begin();
    for (auto root_part = canonical_root.begin(); root_part != canonical_root.end();
         ++root_part, ++child_part) {
        if (child_part == canonical_child.end() || *child_part != *root_part) return false;
    }
    return true;
}

void validate_destination(const fs::path& destination, const std::string& protected_tree) {
    std::error_code error;
    const bool exists = fs::exists(destination, error);
    require(!error, "cannot inspect destination " + destination.string());
    require(!exists, "refusing to overwrite existing destination " + destination.string());
    if (!protected_tree.empty()) {
        require(!same_or_child(destination, fs::path(absolute_path(protected_tree))),
                "destination is equal to or inside the protected legacy C++ tree");
    }
}

std::vector<std::string> summary_header() {
    return {
        "schema_version", "model_version", "sampling_design_version", "category", "group",
        "is_health_care_group", "source_i0", "source_volume_m3",
        "source_probability_susceptible_present", "source_people_outdoor_airflow_l_s_person",
        "source_area_outdoor_airflow_l_s_m2", "source_baseline_outdoor_airflow_l_s_person",
        "source_community_rate", "source_mask_efficiency", "source_pbr_gm_m3_h",
        "source_pbr_gsd", "source_activity_breathing_pct", "source_activity_talking_pct",
        "source_activity_vocalizing_pct", "source_cdrop_gm_e4_per_m3", "source_cdrop_gsd",
        "source_d_gm_e6_m", "source_d_gsd", "source_published_ecai_l_s_person",
        "source_published_rr", "effective_i0", "effective_volume_m3",
        "effective_community_rate", "effective_mask_efficiency", "target_probability",
        "duration_hours", "target_probability_over_duration", "simulations_n", "replicates_r",
        "percentile", "confidence_level", "base_seed", "seed_policy_version",
        "zero_infector_policy", "negative_clamping_policy", "custom_scenario",
        "comparator_eligible", "published_comparator_l_s_person",
        "comparator_signed_difference_l_s_person", "comparator_abs_difference_l_s_person",
        "comparator_status", "dispersion_status", "raw_mean_l_s_person",
        "raw_degrees_of_freedom", "raw_sd_l_s_person", "raw_se_l_s_person", "raw_t_critical",
        "raw_ci_lower_l_s_person", "raw_ci_upper_l_s_person", "raw_min_l_s_person",
        "raw_max_l_s_person", "ceiling_to_5_l_s_person", "summary_row_count",
        "replicate_csv_row_count", "all_replicates_n", "all_replicates_zero_infector_count",
        "all_replicates_active_count", "all_replicates_k_total",
        "all_replicates_negative_clamp_count", "all_replicates_output_zero_count",
        "all_replicates_qer_evaluations", "all_replicates_ecai_evaluations",
        "all_replicates_susceptible_values_consumed",
        "all_replicates_infector_values_consumed", "counters_reconcile",
        "all_cursors_exhausted", "all_replicates_finite", "banks_reconcile",
        "explicit_category", "explicit_simulations", "explicit_replicates", "explicit_percentile",
        "explicit_confidence", "explicit_target_probability", "explicit_duration",
        "explicit_base_seed", "explicit_community_rate", "explicit_occupants", "explicit_volume_m3",
        "explicit_mask_efficiency", "explicit_output_directory", "explicit_summary_file",
        "explicit_replicates_file", "explicit_write_replicates"};
}

std::vector<std::string> replicate_header() {
    return {"replicate_index", "requested_percentile", "raw_percentile_l_s_person", "n",
            "zero_infector_count", "active_count", "k_total", "negative_clamp_count",
            "output_zero_count", "qer_evaluations", "ecai_evaluations",
            "susceptible_values_consumed", "infector_values_consumed", "planned_bank_lengths",
            "final_cursors", "counters_reconcile", "cursors_exhausted", "finite_values",
            "replicate_csv_row_count"};
}

// Probability-mode summary schema (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md
// section 5.2), 96 ordered columns. Columns 8-26 and 79-94 are byte-identical
// to ECAi summary_header() columns 7-25 and 74-89 respectively.
std::vector<std::string> probability_summary_header() {
    return {
        "schema_version", "model_version", "sampling_design_version", "mode", "category",
        "group", "is_health_care_group", "source_i0", "source_volume_m3",
        "source_probability_susceptible_present", "source_people_outdoor_airflow_l_s_person",
        "source_area_outdoor_airflow_l_s_m2", "source_baseline_outdoor_airflow_l_s_person",
        "source_community_rate", "source_mask_efficiency", "source_pbr_gm_m3_h",
        "source_pbr_gsd", "source_activity_breathing_pct", "source_activity_talking_pct",
        "source_activity_vocalizing_pct", "source_cdrop_gm_e4_per_m3", "source_cdrop_gsd",
        "source_d_gm_e6_m", "source_d_gsd", "source_published_ecai_l_s_person",
        "source_published_rr", "effective_i0", "effective_volume_m3",
        "effective_community_rate", "effective_mask_efficiency", "input_ecai_l_s_person",
        "input_ecai_source", "effective_eca_total_m3_h", "target_probability",
        "duration_hours", "target_probability_over_duration", "simulations_n", "replicates_r",
        "percentile", "confidence_level", "base_seed", "seed_policy_version",
        "zero_infector_policy", "probability_policy", "saturation_policy", "custom_scenario",
        "comparator_eligible", "comparator_target_probability", "comparator_signed_difference",
        "comparator_ratio", "comparator_status", "dispersion_status", "raw_mean_probability",
        "raw_degrees_of_freedom", "raw_sd_probability", "raw_se_probability", "raw_t_critical",
        "raw_ci_lower_probability", "raw_ci_upper_probability", "raw_min_probability",
        "raw_max_probability", "summary_row_count", "replicate_csv_row_count",
        "all_replicates_n", "all_replicates_zero_infector_count", "all_replicates_active_count",
        "all_replicates_k_total", "all_replicates_output_zero_count",
        "all_replicates_saturated_count", "all_replicates_qer_evaluations",
        "all_replicates_probability_evaluations", "all_replicates_susceptible_values_consumed",
        "all_replicates_infector_values_consumed", "counters_reconcile", "all_cursors_exhausted",
        "all_replicates_finite", "banks_reconcile",
        // docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md section 5.2 mirrors the ECAi tail
        // (columns 73-89) byte-for-byte, so "banks_reconcile" deliberately appears twice:
        // once as this schema's own reconciliation flag, once copied from ECAi's column 73.
        "banks_reconcile", "explicit_category", "explicit_simulations", "explicit_replicates",
        "explicit_percentile", "explicit_confidence", "explicit_target_probability",
        "explicit_duration", "explicit_base_seed", "explicit_community_rate",
        "explicit_occupants", "explicit_volume_m3", "explicit_mask_efficiency",
        "explicit_output_directory", "explicit_summary_file", "explicit_replicates_file",
        "explicit_write_replicates", "explicit_mode", "explicit_input_ecai"};
}

// Probability-mode replicate schema, 20 ordered columns.
std::vector<std::string> probability_replicate_header() {
    return {"replicate_index", "mode", "requested_percentile", "raw_percentile_probability", "n",
            "zero_infector_count", "active_count", "k_total", "saturated_count",
            "output_zero_count", "qer_evaluations", "probability_evaluations",
            "susceptible_values_consumed", "infector_values_consumed", "planned_bank_lengths",
            "final_cursors", "counters_reconcile", "cursors_exhausted", "finite_values",
            "replicate_csv_row_count"};
}

void validate_probability_stage13_document(const OutputDocument& document) {
    const std::vector<std::string> expected_summary = probability_summary_header();
    const std::vector<std::string> expected_replicates = probability_replicate_header();
    require(document.summary_columns == expected_summary,
            "probability summary schema is not the authoritative 96-column Stage 13 header");
    require(document.replicate_columns == expected_replicates,
            "probability replicate schema is not the authoritative 20-column Stage 13 header");
    validate_document(document);

    const auto& summary = document.summary_rows[0];
    require(summary[61] == "1", "summary_row_count metadata must equal one");
    const std::size_t requested_replicates = parse_size(summary[37], "replicates_r");
    const std::size_t declared_replicates = parse_size(summary[62], "replicate_csv_row_count");
    require(document.replicate_rows.size() == requested_replicates &&
                document.replicate_rows.size() == declared_replicates,
            "probability replicate rows disagree with summary replicate count metadata");
    for (const auto& row : document.replicate_rows) {
        require(parse_size(row[19], "replicate row count") == document.replicate_rows.size(),
                "probability replicate row-count metadata mismatch");
    }
}

void validate_stage13_document(const OutputDocument& document) {
    require(document.summary_rows.size() == 1, "summary must contain exactly one row");
    require(!document.summary_columns.empty() && document.summary_columns[0] == "schema_version" &&
                !document.summary_rows[0].empty(),
            "summary schema missing a leading schema_version column/value");
    if (document.summary_rows[0][0] == kProbabilitySchemaVersion) {
        validate_probability_stage13_document(document);
        return;
    }

    const std::vector<std::string> expected_summary = summary_header();
    const std::vector<std::string> expected_replicates = replicate_header();
    require(document.summary_columns == expected_summary,
            "summary schema is not the authoritative 89-column Stage 13 header");
    require(document.replicate_columns == expected_replicates,
            "replicate schema is not the authoritative 19-column Stage 13 header");
    validate_document(document);

    const auto& summary = document.summary_rows[0];
    require(summary[57] == "1", "summary_row_count metadata must equal one");
    const std::size_t requested_replicates = parse_size(summary[33], "replicates_r");
    const std::size_t declared_replicates = parse_size(summary[58], "replicate_csv_row_count");
    require(document.replicate_rows.size() == requested_replicates &&
                document.replicate_rows.size() == declared_replicates,
            "replicate rows disagree with summary replicate count metadata");
    for (const auto& row : document.replicate_rows) {
        require(parse_size(row[18], "replicate row count") == document.replicate_rows.size(),
                "replicate row-count metadata mismatch");
    }
}

std::string bank_sizes(const std::array<std::size_t, 11>& values) {
    std::ostringstream text;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i != 0) text << ';';
        text << values[i];
    }
    return text.str();
}

void require_same_summary(const statistics::ReplicateSummary& supplied,
                          const statistics::ReplicateSummary& authoritative) {
    require(supplied.replicate_count == authoritative.replicate_count,
            "summary replicate count mismatch");
    require(supplied.replicate_estimates_l_s_person.size() ==
                authoritative.replicate_estimates_l_s_person.size(),
            "summary raw vector length mismatch");
    for (std::size_t i = 0; i < authoritative.replicate_estimates_l_s_person.size(); ++i)
        require(same_double(supplied.replicate_estimates_l_s_person[i],
                            authoritative.replicate_estimates_l_s_person[i]),
                "summary raw vector/order mismatch");
    require(same_double(supplied.confidence_level, authoritative.confidence_level),
            "summary confidence mismatch");
    require(supplied.dispersion_status == authoritative.dispersion_status,
            "summary dispersion status mismatch");
    require(same_double(supplied.mean_l_s_person, authoritative.mean_l_s_person),
            "summary mean mismatch");
    require(same_double(supplied.min_l_s_person, authoritative.min_l_s_person),
            "summary minimum mismatch");
    require(same_double(supplied.max_l_s_person, authoritative.max_l_s_person),
            "summary maximum mismatch");
    require(same_optional(supplied.degrees_of_freedom, authoritative.degrees_of_freedom),
            "summary degrees of freedom mismatch");
    require(same_optional(supplied.sample_sd_l_s_person, authoritative.sample_sd_l_s_person),
            "summary sample SD mismatch");
    require(same_optional(supplied.standard_error_l_s_person,
                          authoritative.standard_error_l_s_person),
            "summary standard error mismatch");
    require(same_optional(supplied.t_critical, authoritative.t_critical),
            "summary t critical mismatch");
    require(same_optional(supplied.ci_lower_l_s_person, authoritative.ci_lower_l_s_person),
            "summary CI lower mismatch");
    require(same_optional(supplied.ci_upper_l_s_person, authoritative.ci_upper_l_s_person),
            "summary CI upper mismatch");
    require(same_double(supplied.adjusted_ceiling_to_5_l_s_person,
                        authoritative.adjusted_ceiling_to_5_l_s_person),
            "summary ceiling mismatch");
}

void require_same_run(const cli::ResolvedRun& supplied, const cli::ResolvedRun& authoritative) {
    require(supplied.category == authoritative.category,
            "resolved category pointer/name disagrees with CLI options");
    const auto& left = supplied.config;
    const auto& right = authoritative.config;
    require(left.category_name == right.category_name && left.i0 == right.i0 &&
                same_double(left.volume_m3, right.volume_m3) &&
                same_double(left.mask_efficiency, right.mask_efficiency) &&
                same_double(left.community_infection_rate, right.community_infection_rate) &&
                same_double(left.target_probability, right.target_probability) &&
                same_double(left.duration_hours, right.duration_hours) &&
                same_double(left.percentile, right.percentile) &&
                same_double(left.confidence_level, right.confidence_level) &&
                left.simulations_n == right.simulations_n && left.replicates_r == right.replicates_r,
            "resolved effective configuration disagrees with CLI options");
    require(supplied.overrides.occupants == authoritative.overrides.occupants &&
                same_optional(supplied.overrides.volume_m3, authoritative.overrides.volume_m3) &&
                same_optional(supplied.overrides.mask_efficiency,
                              authoritative.overrides.mask_efficiency) &&
                same_optional(supplied.overrides.community_infection_rate,
                              authoritative.overrides.community_infection_rate) &&
                supplied.overrides.mode == authoritative.overrides.mode &&
                same_optional(supplied.overrides.input_ecai_l_s_person,
                              authoritative.overrides.input_ecai_l_s_person),
            "resolved overrides disagree with CLI options");
    const auto& ls = supplied.scenario;
    const auto& rs = authoritative.scenario;
    require(ls.category_name == rs.category_name && ls.occupants == rs.occupants &&
                same_double(ls.volume_m3, rs.volume_m3) &&
                same_double(ls.mask_efficiency, rs.mask_efficiency) &&
                same_double(ls.community_infection_rate, rs.community_infection_rate) &&
                same_double(ls.target_probability, rs.target_probability) &&
                same_double(ls.duration_hours, rs.duration_hours) &&
                same_double(ls.percentile, rs.percentile) && ls.simulations_n == rs.simulations_n &&
                ls.occupants_overridden == rs.occupants_overridden &&
                ls.volume_overridden == rs.volume_overridden &&
                ls.mask_efficiency_overridden == rs.mask_efficiency_overridden &&
                ls.community_infection_rate_overridden == rs.community_infection_rate_overridden &&
                ls.custom_scenario == rs.custom_scenario && ls.mode == rs.mode &&
                same_double(ls.input_ecai_l_s_person, rs.input_ecai_l_s_person) &&
                ls.input_ecai_overridden == rs.input_ecai_overridden,
            "resolved effective scenario disagrees with CLI options");
    require(supplied.base_seed == authoritative.base_seed &&
                supplied.custom_scenario == authoritative.custom_scenario &&
                supplied.target_is_hourly == authoritative.target_is_hourly,
            "resolved seed/custom/hourly provenance disagrees with CLI options");
}

void record_missing_parents(const fs::path& parent, std::vector<fs::path>& owned) {
    std::vector<fs::path> missing;
    fs::path current = parent;
    std::error_code error;
    while (!current.empty() && !fs::exists(current, error)) {
        require(!error, "cannot inspect output parent directory");
        missing.push_back(current);
        current = current.parent_path();
    }
    require(!error, "cannot inspect output parent directory");
    std::reverse(missing.begin(), missing.end());
    for (const auto& path : missing) {
        if (std::find(owned.begin(), owned.end(), path) == owned.end()) owned.push_back(path);
    }
}

bool cleanup_temporary_payloads(const std::vector<fs::path>& directories,
                                const std::vector<fs::path>& payload_paths) noexcept {
    bool complete = true;
    for (std::size_t i = 0; i < directories.size(); ++i) {
        std::error_code error;
        if (i < payload_paths.size()) {
            fs::remove(payload_paths[i], error);
            if (error) complete = false;
        }
        error.clear();
        fs::remove(directories[i], error);  // removes only an empty owned directory
        if (error) {
            complete = false;
        } else {
            const bool remains = fs::exists(directories[i], error);
            if (error || remains) complete = false;
        }
    }
    return complete;
}

fs::path reserve_temp_directory(const fs::path& destination, std::random_device& random) {
    for (unsigned attempt = 0; attempt < 64; ++attempt) {
        std::ostringstream name;
        name << destination.filename().string() << ".tmp-" << std::hex << random() << random()
             << attempt;
        const fs::path candidate = destination.parent_path() / name.str();
        std::error_code error;
        if (fs::create_directory(candidate, error)) return candidate;
        if (error && error != std::errc::file_exists)
            require(false, "cannot reserve temporary directory: " + error.message());
    }
    require(false, "cannot reserve unique temporary directory");
    return {};
}

}  // namespace

namespace detail {
void set_write_test_hook(WriteTestHook hook) { g_write_test_hook = std::move(hook); }
}  // namespace detail

std::string serialize_summary_csv(const OutputDocument& document) {
    validate_document(document);
    return serialize(document.summary_columns, document.summary_rows);
}

std::string serialize_replicates_csv(const OutputDocument& document) {
    validate_document(document);
    return serialize(document.replicate_columns, document.replicate_rows);
}

OutputDocument build_probability_document(const cli::ResolvedRun& run,
                                          const simulation::ReplicateBatch& batch,
                                          const cli::CliOptions& options);

OutputDocument build_document(const cli::ResolvedRun& run,
                              const simulation::ReplicateBatch& batch,
                              const cli::CliOptions& options,
                              const std::string&) {
    if (run.scenario.mode == simulation::OutputMode::probability) {
        return build_probability_document(run, batch, options);
    }
    require(run.overrides.mode == simulation::OutputMode::ecai,
            "ECAi document requested for a non-ECAi batch");
    for (const simulation::ReplicateResult& replicate : batch.replicates) {
        require(replicate.mode == simulation::OutputMode::ecai &&
                    replicate.probability.empty() && replicate.percentile_probability == 0.0 &&
                    replicate.counters.probability_evaluations == 0 &&
                    replicate.counters.saturated_count == 0,
                "ECAi document requested for a non-ECAi batch");
    }

    const cli::ResolvedRun authoritative_run = cli::resolve_run(options);
    require_same_run(run, authoritative_run);
    require(run.category != nullptr, "null category");
    const OccupancyCategory& category = *run.category;
    const simulation::EffectiveScenario& scenario = run.scenario;

    require(run.category == &lookup_occupancy_category(category.canonical_name),
            "category pointer does not identify the canonical source record");
    require(category.canonical_name == run.config.category_name &&
                run.config.category_name == scenario.category_name,
            "category/config/scenario disagreement");
    const bool category_custom = options.occupants.has_value() || options.volume_m3.has_value() ||
                                 options.mask_efficiency.has_value() ||
                                 options.community_infection_rate.has_value();
    require(run.custom_scenario == options.custom_scenario &&
                scenario.custom_scenario == category_custom,
            "scenario provenance disagreement");
    require(run.base_seed == options.base_seed, "base-seed provenance disagreement");
    require(run.overrides.occupants == options.occupants &&
                run.overrides.volume_m3 == options.volume_m3 &&
                run.overrides.mask_efficiency == options.mask_efficiency &&
                run.overrides.community_infection_rate == options.community_infection_rate,
            "override provenance disagreement");
    require(scenario.occupants_overridden == options.occupants.has_value() &&
                scenario.volume_overridden == options.volume_m3.has_value() &&
                scenario.mask_efficiency_overridden == options.mask_efficiency.has_value() &&
                scenario.community_infection_rate_overridden ==
                    options.community_infection_rate.has_value(),
            "per-field scenario provenance disagreement");
    require(run.target_is_hourly == options.target_is_hourly &&
                run.target_is_hourly == (scenario.duration_hours == 1.0) &&
                run.config.i0 == scenario.occupants && run.config.volume_m3 == scenario.volume_m3 &&
                run.config.mask_efficiency == scenario.mask_efficiency &&
                run.config.community_infection_rate == scenario.community_infection_rate &&
                run.config.simulations_n == scenario.simulations_n &&
                run.config.simulations_n == options.simulations_n &&
                run.config.replicates_r == options.replicates_r &&
                run.config.target_probability == scenario.target_probability &&
                run.config.target_probability == options.target_probability &&
                run.config.duration_hours == scenario.duration_hours &&
                run.config.duration_hours == options.duration_hours &&
                run.config.percentile == scenario.percentile &&
                run.config.percentile == options.percentile &&
                run.config.confidence_level == options.confidence_level,
            "configuration/statistics provenance disagreement");
    require(batch.replicates.size() == options.replicates_r,
            "batch replicate count disagrees with requested R");

    std::array<std::size_t, 10> totals{};
    std::vector<double> raw_percentiles;
    raw_percentiles.reserve(batch.replicates.size());
    for (std::size_t index = 0; index < batch.replicates.size(); ++index) {
        const simulation::ReplicateResult& replicate = batch.replicates[index];
        const simulation::ReplicateCounters& counters = replicate.counters;
        require(replicate.replicate_index == index, "replicate identity/order mismatch");
        require(same_double(replicate.percentile, scenario.percentile),
                "replicate requested percentile mismatch");
        require(counters.n == scenario.simulations_n,
                "replicate counter N disagrees with effective N");
        require(replicate.ecai_l_s_person.size() == scenario.simulations_n,
                "replicate output count disagrees with effective N");
        require(replicate.infector_counts.size() == scenario.simulations_n,
                "replicate infector-count vector disagrees with effective N");
        std::size_t grounded_zero = 0;
        std::size_t grounded_active = 0;
        std::size_t grounded_k_total = 0;
        std::size_t zero_count = 0;
        for (std::size_t row = 0; row < replicate.ecai_l_s_person.size(); ++row) {
            const long long count = replicate.infector_counts[row];
            const double value = replicate.ecai_l_s_person[row];
            require(count >= 0 && count <= scenario.occupants,
                    "replicate infector count lies outside [0, effective I0]");
            if (count == 0) {
                grounded_zero = checked_add(grounded_zero, 1, "grounded zero infector");
                require(value == 0.0 && !std::signbit(value),
                        "zero-infector row output is not exact positive +0.0");
            } else {
                grounded_active = checked_add(grounded_active, 1, "grounded active");
                grounded_k_total = checked_add(grounded_k_total, static_cast<std::size_t>(count),
                                               "grounded K_total");
            }
            require(std::isfinite(value), "replicate contains non-finite ECAi output");
            require(value >= 0.0, "replicate contains negative ECAi output");
            if (value == 0.0) zero_count = checked_add(zero_count, 1, "output zero");
        }
        require(counters.zero_infector_count == grounded_zero &&
                    counters.active_count == grounded_active && counters.k_total == grounded_k_total,
                "replicate counters disagree with authoritative Stage 6 counts");
        require(counters.output_zero_count == zero_count, "output zero count mismatch");
        require(counters.active_count <= counters.n &&
                    counters.zero_infector_count == counters.n - counters.active_count,
                "active/zero-infector partition mismatch");
        require(counters.negative_clamp_count <= counters.active_count,
                "negative clamp count exceeds active count");
        require(counters.output_zero_count >= counters.zero_infector_count &&
                    counters.negative_clamp_count <=
                        counters.output_zero_count - counters.zero_infector_count,
                "output zero sources do not reconcile");
        require(counters.qer_evaluations == counters.k_total,
                "QER count disagrees with K_total");
        require(counters.ecai_evaluations == counters.active_count,
                "ECAi count disagrees with active count");
        require(counters.susceptible_values_consumed ==
                    checked_multiply(3, counters.active_count, "susceptible"),
                "susceptible consumption mismatch");
        require(counters.infector_values_consumed ==
                    checked_multiply(8, counters.k_total, "infector"),
                "infector consumption mismatch");
        for (std::size_t bank = 0; bank < 11; ++bank) {
            const std::size_t expected = bank < 3 ? counters.active_count : counters.k_total;
            require(replicate.planned_bank_lengths[bank] == expected,
                    "planned bank length mismatch");
            require(replicate.final_cursors[bank] == expected, "final bank cursor mismatch");
            require(replicate.bank_generated[bank] == (expected > 0),
                    "generated-bank flag mismatch");
        }
        require(replicate.cursors_exhausted, "replicate cursors are not exhausted");
        const double percentile =
            simulation::percentile_linear(replicate.ecai_l_s_person, scenario.percentile);
        require(same_double(replicate.percentile_l_s_person, percentile),
                "stored replicate percentile mismatch");
        raw_percentiles.push_back(percentile);

        totals[0] = checked_add(totals[0], counters.n, "n");
        totals[1] = checked_add(totals[1], counters.zero_infector_count, "zero infector");
        totals[2] = checked_add(totals[2], counters.active_count, "active");
        totals[3] = checked_add(totals[3], counters.k_total, "K_total");
        totals[4] = checked_add(totals[4], counters.negative_clamp_count, "negative clamp");
        totals[5] = checked_add(totals[5], counters.output_zero_count, "output zero");
        totals[6] = checked_add(totals[6], counters.qer_evaluations, "QER");
        totals[7] = checked_add(totals[7], counters.ecai_evaluations, "ECAi");
        totals[8] = checked_add(totals[8], counters.susceptible_values_consumed, "susceptible");
        totals[9] = checked_add(totals[9], counters.infector_values_consumed, "infector");
    }

    const statistics::ReplicateSummary authoritative =
        statistics::summarize_replicates(raw_percentiles, options.confidence_level);
    require_same_summary(batch.summary, authoritative);

    const bool comparator_eligible = !run.custom_scenario;
    std::string published;
    std::string signed_difference;
    std::string absolute_difference;
    std::string comparator_status = "NOT_APPLICABLE_CUSTOM";
    if (comparator_eligible) {
        require(std::fmod(authoritative.adjusted_ceiling_to_5_l_s_person, 5.0) == 0.0 &&
                    std::fmod(category.published_ecai_Lsp, 5.0) == 0.0,
                "publication comparison values must lie on the exact 5 L/s/person grid");
        const double difference =
            authoritative.adjusted_ceiling_to_5_l_s_person - category.published_ecai_Lsp;
        published = format_double(category.published_ecai_Lsp);
        signed_difference = format_double(difference);
        absolute_difference = format_double(std::abs(difference));
        comparator_status = difference == 0.0
                                ? "EXACT_MATCH"
                                : (std::abs(difference) == 5.0
                                       ? "NEAR_WITHIN_ONE_INCREMENT"
                                       : "FAIL_GT_ONE_INCREMENT");
    }

    const auto& q = authoritative;
    const std::string dispersion_status =
        q.dispersion_status == statistics::DispersionStatus::available
            ? "AVAILABLE"
            : "UNAVAILABLE_SINGLE_REPLICATE";
    std::vector<std::string> row = {
        kOutputSchemaVersion, kModelVersion, std::to_string(seed_manager::kDesignVersion),
        category.canonical_name, category.group, boolean(category.is_health_care_group),
        std::to_string(category.i0), format_double(category.volume_m3),
        format_double(category.probability_susceptible_present),
        format_double(category.people_outdoor_airflow_Lsp),
        format_double(category.area_outdoor_airflow_Lsm2),
        format_double(category.baseline_outdoor_airflow_Lsp),
        format_double(category.community_infection_rate), format_double(category.mask_efficiency),
        format_double(category.pbr_gm_m3h), format_double(category.pbr_gsd),
        format_double(category.activity_breathing_pct), format_double(category.activity_talking_pct),
        format_double(category.activity_vocalizing_pct), format_double(category.cdrop_gm_e4_per_m3),
        format_double(category.cdrop_gsd), format_double(category.d_gm_e6_m),
        format_double(category.d_gsd), format_double(category.published_ecai_Lsp),
        format_double(category.published_rr), std::to_string(scenario.occupants),
        format_double(scenario.volume_m3), format_double(scenario.community_infection_rate),
        format_double(scenario.mask_efficiency), format_double(scenario.target_probability),
        format_double(scenario.duration_hours),
        scenario.duration_hours == 1.0 ? "hourly_target" : "target_over_selected_duration",
        std::to_string(scenario.simulations_n), std::to_string(options.replicates_r),
        format_double(scenario.percentile), format_double(options.confidence_level),
        std::to_string(run.base_seed), seed_manager::kSeedPolicyVersion,
        "zero_infector_short_circuit", "negative_ecai_clamped_to_zero",
        boolean(run.custom_scenario), boolean(comparator_eligible), published, signed_difference,
        absolute_difference, comparator_status, dispersion_status, format_double(q.mean_l_s_person),
        format_optional(q.degrees_of_freedom), format_optional(q.sample_sd_l_s_person),
        format_optional(q.standard_error_l_s_person), format_optional(q.t_critical),
        format_optional(q.ci_lower_l_s_person), format_optional(q.ci_upper_l_s_person),
        format_double(q.min_l_s_person), format_double(q.max_l_s_person),
        format_double(q.adjusted_ceiling_to_5_l_s_person), "1",
        std::to_string(batch.replicates.size()), std::to_string(totals[0]),
        std::to_string(totals[1]), std::to_string(totals[2]), std::to_string(totals[3]),
        std::to_string(totals[4]), std::to_string(totals[5]), std::to_string(totals[6]),
        std::to_string(totals[7]), std::to_string(totals[8]), std::to_string(totals[9]),
        "true", "true", "true", "true", boolean(options.category_supplied),
        boolean(options.simulations_explicit), boolean(options.replicates_explicit),
        boolean(options.percentile_explicit), boolean(options.confidence_explicit),
        boolean(options.target_probability_explicit), boolean(options.duration_explicit),
        boolean(options.base_seed_explicit), boolean(options.community_infection_rate.has_value()),
        boolean(options.occupants.has_value()), boolean(options.volume_m3.has_value()),
        boolean(options.mask_efficiency.has_value()), boolean(options.output_directory.has_value()),
        boolean(options.summary_file.has_value()), boolean(options.replicates_file.has_value()),
        boolean(options.write_replicates)};

    OutputDocument document;
    document.summary_columns = summary_header();
    require(document.summary_columns.size() == 89, "summary schema width changed");
    require(row.size() == document.summary_columns.size(), "summary schema construction mismatch");
    document.summary_rows.push_back(std::move(row));
    document.replicate_columns = replicate_header();
    require(document.replicate_columns.size() == 19, "replicate schema width changed");
    for (const auto& replicate : batch.replicates) {
        const auto& counters = replicate.counters;
        document.replicate_rows.push_back({
            std::to_string(replicate.replicate_index), format_double(replicate.percentile),
            format_double(replicate.percentile_l_s_person), std::to_string(counters.n),
            std::to_string(counters.zero_infector_count), std::to_string(counters.active_count),
            std::to_string(counters.k_total), std::to_string(counters.negative_clamp_count),
            std::to_string(counters.output_zero_count), std::to_string(counters.qer_evaluations),
            std::to_string(counters.ecai_evaluations),
            std::to_string(counters.susceptible_values_consumed),
            std::to_string(counters.infector_values_consumed),
            bank_sizes(replicate.planned_bank_lengths), bank_sizes(replicate.final_cursors),
            "true", "true", "true", std::to_string(batch.replicates.size())});
    }
    return document;
}

OutputDocument build_probability_document(const cli::ResolvedRun& run,
                                          const simulation::ReplicateBatch& batch,
                                          const cli::CliOptions& options) {
    const cli::ResolvedRun authoritative_run = cli::resolve_run(options);
    require_same_run(run, authoritative_run);
    require(run.category != nullptr, "null category");
    const OccupancyCategory& category = *run.category;
    const simulation::EffectiveScenario& scenario = run.scenario;

    require(run.category == &lookup_occupancy_category(category.canonical_name),
            "category pointer does not identify the canonical source record");
    require(category.canonical_name == run.config.category_name &&
                run.config.category_name == scenario.category_name,
            "category/config/scenario disagreement");
    require(run.overrides.mode == simulation::OutputMode::probability &&
                scenario.mode == simulation::OutputMode::probability &&
                options.mode == simulation::OutputMode::probability,
            "build_probability_document requires a probability-mode resolved run");
    const bool category_custom = options.occupants.has_value() || options.volume_m3.has_value() ||
                                 options.mask_efficiency.has_value() ||
                                 options.community_infection_rate.has_value();
    require(run.custom_scenario == options.custom_scenario &&
                scenario.custom_scenario == category_custom,
            "scenario provenance disagreement");
    require(run.base_seed == options.base_seed, "base-seed provenance disagreement");
    require(run.overrides.occupants == options.occupants &&
                run.overrides.volume_m3 == options.volume_m3 &&
                run.overrides.mask_efficiency == options.mask_efficiency &&
                run.overrides.community_infection_rate == options.community_infection_rate,
            "override provenance disagreement");
    require(scenario.occupants_overridden == options.occupants.has_value() &&
                scenario.volume_overridden == options.volume_m3.has_value() &&
                scenario.mask_efficiency_overridden == options.mask_efficiency.has_value() &&
                scenario.community_infection_rate_overridden ==
                    options.community_infection_rate.has_value(),
            "per-field scenario provenance disagreement");
    require(scenario.input_ecai_overridden == options.input_ecai_l_s_person.has_value(),
            "input ECAi provenance disagreement");
    require(same_optional(run.overrides.input_ecai_l_s_person, options.input_ecai_l_s_person),
            "input ECAi override disagreement");
    require(run.target_is_hourly == options.target_is_hourly &&
                run.target_is_hourly == (scenario.duration_hours == 1.0) &&
                run.config.i0 == scenario.occupants && run.config.volume_m3 == scenario.volume_m3 &&
                run.config.mask_efficiency == scenario.mask_efficiency &&
                run.config.community_infection_rate == scenario.community_infection_rate &&
                run.config.simulations_n == scenario.simulations_n &&
                run.config.simulations_n == options.simulations_n &&
                run.config.replicates_r == options.replicates_r &&
                run.config.target_probability == scenario.target_probability &&
                run.config.target_probability == options.target_probability &&
                run.config.duration_hours == scenario.duration_hours &&
                run.config.duration_hours == options.duration_hours &&
                run.config.percentile == scenario.percentile &&
                run.config.percentile == options.percentile &&
                run.config.confidence_level == options.confidence_level,
            "configuration/statistics provenance disagreement");
    require(batch.replicates.size() == options.replicates_r,
            "batch replicate count disagrees with requested R");

    std::array<std::size_t, 10> totals{};
    std::vector<double> raw_percentiles;
    raw_percentiles.reserve(batch.replicates.size());
    for (std::size_t index = 0; index < batch.replicates.size(); ++index) {
        const simulation::ReplicateResult& replicate = batch.replicates[index];
        const simulation::ReplicateCounters& counters = replicate.counters;
        require(replicate.replicate_index == index, "replicate identity/order mismatch");
        require(replicate.mode == simulation::OutputMode::probability,
                "replicate did not resolve to probability mode");
        require(same_double(replicate.percentile, scenario.percentile),
                "replicate requested percentile mismatch");
        require(counters.n == scenario.simulations_n,
                "replicate counter N disagrees with effective N");
        require(replicate.ecai_l_s_person.empty(),
                "probability-mode replicate must not carry ECAi-mode output");
        require(replicate.percentile_l_s_person == 0.0,
                "probability-mode replicate must not carry an ECAi-mode percentile");
        require(counters.negative_clamp_count == 0,
                "probability-mode replicate must not carry ECAi-mode negative clamps");
        require(counters.ecai_evaluations == 0,
                "probability-mode replicate must not carry ECAi-mode evaluation counts");
        require(replicate.probability.size() == scenario.simulations_n,
                "replicate output count disagrees with effective N");
        require(replicate.infector_counts.size() == scenario.simulations_n,
                "replicate infector-count vector disagrees with effective N");
        std::size_t grounded_zero = 0;
        std::size_t grounded_active = 0;
        std::size_t grounded_k_total = 0;
        std::size_t zero_output_count = 0;
        std::size_t saturated_count = 0;
        for (std::size_t row = 0; row < replicate.probability.size(); ++row) {
            const long long count = replicate.infector_counts[row];
            const double value = replicate.probability[row];
            require(count >= 0 && count <= scenario.occupants,
                    "replicate infector count lies outside [0, effective I0]");
            if (count == 0) {
                grounded_zero = checked_add(grounded_zero, 1, "grounded zero infector");
                require(value == 0.0 && !std::signbit(value),
                        "zero-infector row output is not exact positive +0.0");
            } else {
                grounded_active = checked_add(grounded_active, 1, "grounded active");
                grounded_k_total = checked_add(grounded_k_total, static_cast<std::size_t>(count),
                                               "grounded K_total");
            }
            require(std::isfinite(value), "replicate contains a non-finite probability");
            require(value >= 0.0 && value <= 1.0,
                    "replicate contains a probability outside the closed unit interval");
            if (value == 0.0) zero_output_count = checked_add(zero_output_count, 1, "output zero");
            if (value == 1.0) saturated_count = checked_add(saturated_count, 1, "saturated");
        }
        require(counters.zero_infector_count == grounded_zero &&
                    counters.active_count == grounded_active && counters.k_total == grounded_k_total,
                "replicate counters disagree with authoritative Stage 6 counts");
        require(counters.output_zero_count == zero_output_count, "output zero count mismatch");
        require(counters.saturated_count == saturated_count, "saturated count mismatch");
        require(counters.active_count <= counters.n &&
                    counters.zero_infector_count == counters.n - counters.active_count,
                "active/zero-infector partition mismatch");
        require(counters.saturated_count <= counters.active_count,
                "saturated count exceeds active count");
        require(counters.output_zero_count >= counters.zero_infector_count,
                "output zero sources do not reconcile");
        require(counters.qer_evaluations == counters.k_total,
                "QER count disagrees with K_total");
        require(counters.probability_evaluations == counters.active_count,
                "probability evaluation count disagrees with active count");
        require(counters.susceptible_values_consumed ==
                    checked_multiply(3, counters.active_count, "susceptible"),
                "susceptible consumption mismatch");
        require(counters.infector_values_consumed ==
                    checked_multiply(8, counters.k_total, "infector"),
                "infector consumption mismatch");
        for (std::size_t bank = 0; bank < 11; ++bank) {
            const std::size_t expected = bank < 3 ? counters.active_count : counters.k_total;
            require(replicate.planned_bank_lengths[bank] == expected,
                    "planned bank length mismatch");
            require(replicate.final_cursors[bank] == expected, "final bank cursor mismatch");
            require(replicate.bank_generated[bank] == (expected > 0),
                    "generated-bank flag mismatch");
        }
        require(replicate.cursors_exhausted, "replicate cursors are not exhausted");
        const double percentile =
            simulation::percentile_linear(replicate.probability, scenario.percentile);
        require(same_double(replicate.percentile_probability, percentile),
                "stored replicate percentile mismatch");
        raw_percentiles.push_back(percentile);

        totals[0] = checked_add(totals[0], counters.n, "n");
        totals[1] = checked_add(totals[1], counters.zero_infector_count, "zero infector");
        totals[2] = checked_add(totals[2], counters.active_count, "active");
        totals[3] = checked_add(totals[3], counters.k_total, "K_total");
        totals[4] = checked_add(totals[4], counters.output_zero_count, "output zero");
        totals[5] = checked_add(totals[5], counters.saturated_count, "saturated");
        totals[6] = checked_add(totals[6], counters.qer_evaluations, "QER");
        totals[7] = checked_add(totals[7], counters.probability_evaluations, "probability");
        totals[8] = checked_add(totals[8], counters.susceptible_values_consumed, "susceptible");
        totals[9] = checked_add(totals[9], counters.infector_values_consumed, "infector");
    }

    const statistics::ReplicateSummary authoritative =
        statistics::summarize_replicates(raw_percentiles, options.confidence_level);
    require_same_summary(batch.summary, authoritative);

    const bool comparator_eligible = !run.custom_scenario;
    std::string comparator_target;
    std::string comparator_signed_difference;
    std::string comparator_ratio_text;
    std::string comparator_status = "NOT_APPLICABLE_CUSTOM";
    if (comparator_eligible) {
        const double target = scenario.target_probability;
        const double mean = authoritative.mean_l_s_person;
        const double difference = mean - target;
        const double ratio = mean / target;
        comparator_target = format_double(target);
        comparator_signed_difference = format_double(difference);
        comparator_ratio_text = format_double(ratio);
        // The [0.9, 1.1] band is closed on both ends (PROBABILITY_MODE_IMPLEMENTATION_PLAN.md
        // §9.3 T9), but `ratio = mean / target` is a single IEEE-754 division: when `mean` is
        // mathematically exactly `0.9 * target` (or `1.1 * target`), rounding the division can
        // land one ULP to the wrong side of the literal 0.9/1.1 comparison (e.g. mean ==
        // 0.0009, target == 0.001 rounds to ratio == 0.89999999999999991...). Tolerate a few ULPs
        // of rounding noise around each boundary so the documented closed interval is honored;
        // this tolerance is orders of magnitude smaller than any ratio that is genuinely outside
        // the band.
        constexpr double kBoundaryUlpTolerance = 4.0 * std::numeric_limits<double>::epsilon();
        const auto near_boundary = [](double value, double boundary) {
            return std::abs(value - boundary) <= kBoundaryUlpTolerance * boundary;
        };
        const bool at_or_above_lower = ratio >= 0.9 || near_boundary(ratio, 0.9);
        const bool at_or_below_upper = ratio <= 1.1 || near_boundary(ratio, 1.1);
        comparator_status = (at_or_above_lower && at_or_below_upper)
                                ? "AT_TARGET_WITHIN_10PCT"
                                : (!at_or_above_lower ? "BELOW_TARGET" : "ABOVE_TARGET");
    }

    const auto& q = authoritative;
    const std::string dispersion_status =
        q.dispersion_status == statistics::DispersionStatus::available
            ? "AVAILABLE"
            : "UNAVAILABLE_SINGLE_REPLICATE";
    const std::string input_ecai_source =
        scenario.input_ecai_overridden ? "explicit" : "published_table4";
    const double eca_total = model::kM3PerHourPerLitrePerSecond *
                             static_cast<double>(scenario.occupants) *
                             scenario.input_ecai_l_s_person;
    const std::string banks_reconcile_value = "true";

    std::vector<std::string> row = {
        kProbabilitySchemaVersion, kModelVersion, std::to_string(seed_manager::kDesignVersion),
        "probability", category.canonical_name, category.group,
        boolean(category.is_health_care_group), std::to_string(category.i0),
        format_double(category.volume_m3), format_double(category.probability_susceptible_present),
        format_double(category.people_outdoor_airflow_Lsp),
        format_double(category.area_outdoor_airflow_Lsm2),
        format_double(category.baseline_outdoor_airflow_Lsp),
        format_double(category.community_infection_rate), format_double(category.mask_efficiency),
        format_double(category.pbr_gm_m3h), format_double(category.pbr_gsd),
        format_double(category.activity_breathing_pct), format_double(category.activity_talking_pct),
        format_double(category.activity_vocalizing_pct), format_double(category.cdrop_gm_e4_per_m3),
        format_double(category.cdrop_gsd), format_double(category.d_gm_e6_m),
        format_double(category.d_gsd), format_double(category.published_ecai_Lsp),
        format_double(category.published_rr), std::to_string(scenario.occupants),
        format_double(scenario.volume_m3), format_double(scenario.community_infection_rate),
        format_double(scenario.mask_efficiency), format_double(scenario.input_ecai_l_s_person),
        input_ecai_source, format_double(eca_total), format_double(scenario.target_probability),
        format_double(scenario.duration_hours),
        scenario.duration_hours == 1.0 ? "hourly_target" : "target_over_selected_duration",
        std::to_string(scenario.simulations_n), std::to_string(options.replicates_r),
        format_double(scenario.percentile), format_double(options.confidence_level),
        std::to_string(run.base_seed), seed_manager::kSeedPolicyVersion,
        "zero_infector_short_circuit", "forward_wells_riley_expm1",
        "probability_one_is_numerical_saturation", boolean(run.custom_scenario),
        boolean(comparator_eligible), comparator_target, comparator_signed_difference,
        comparator_ratio_text, comparator_status, dispersion_status,
        format_double(q.mean_l_s_person), format_optional(q.degrees_of_freedom),
        format_optional(q.sample_sd_l_s_person), format_optional(q.standard_error_l_s_person),
        format_optional(q.t_critical), format_optional(q.ci_lower_l_s_person),
        format_optional(q.ci_upper_l_s_person), format_double(q.min_l_s_person),
        format_double(q.max_l_s_person), "1", std::to_string(batch.replicates.size()),
        std::to_string(totals[0]), std::to_string(totals[1]), std::to_string(totals[2]),
        std::to_string(totals[3]), std::to_string(totals[4]), std::to_string(totals[5]),
        std::to_string(totals[6]), std::to_string(totals[7]), std::to_string(totals[8]),
        std::to_string(totals[9]), "true", "true", "true", banks_reconcile_value,
        banks_reconcile_value, boolean(options.category_supplied),
        boolean(options.simulations_explicit), boolean(options.replicates_explicit),
        boolean(options.percentile_explicit), boolean(options.confidence_explicit),
        boolean(options.target_probability_explicit), boolean(options.duration_explicit),
        boolean(options.base_seed_explicit), boolean(options.community_infection_rate.has_value()),
        boolean(options.occupants.has_value()), boolean(options.volume_m3.has_value()),
        boolean(options.mask_efficiency.has_value()), boolean(options.output_directory.has_value()),
        boolean(options.summary_file.has_value()), boolean(options.replicates_file.has_value()),
        boolean(options.write_replicates), boolean(options.mode_explicit),
        boolean(options.input_ecai_l_s_person.has_value())};

    OutputDocument document;
    document.summary_columns = probability_summary_header();
    require(document.summary_columns.size() == 96, "probability summary schema width changed");
    require(row.size() == document.summary_columns.size(),
            "probability summary schema construction mismatch");
    document.summary_rows.push_back(std::move(row));
    document.replicate_columns = probability_replicate_header();
    require(document.replicate_columns.size() == 20, "probability replicate schema width changed");
    for (const auto& replicate : batch.replicates) {
        const auto& counters = replicate.counters;
        document.replicate_rows.push_back({
            std::to_string(replicate.replicate_index), "probability",
            format_double(replicate.percentile), format_double(replicate.percentile_probability),
            std::to_string(counters.n), std::to_string(counters.zero_infector_count),
            std::to_string(counters.active_count), std::to_string(counters.k_total),
            std::to_string(counters.saturated_count), std::to_string(counters.output_zero_count),
            std::to_string(counters.qer_evaluations),
            std::to_string(counters.probability_evaluations),
            std::to_string(counters.susceptible_values_consumed),
            std::to_string(counters.infector_values_consumed),
            bank_sizes(replicate.planned_bank_lengths), bank_sizes(replicate.final_cursors),
            "true", "true", "true", std::to_string(batch.replicates.size())});
    }
    return document;
}

OutputPaths resolve_paths(const cli::CliOptions& options, const std::string& protected_tree) {
    OutputPaths paths;
    paths.requested = options.output_directory.has_value() || options.summary_file.has_value() ||
                      options.replicates_file.has_value() || options.write_replicates;
    if (!paths.requested) return paths;
    require(!(options.replicates_file && !options.write_replicates),
            "--replicates-file requires --write-replicates");
    require(!options.write_replicates || options.replicates_file || options.output_directory,
            "--write-replicates requires --replicates-file or --output-directory");
    if (options.summary_file) {
        paths.summary_requested = true;
        paths.summary_path = absolute_path(*options.summary_file);
    } else if (options.output_directory) {
        paths.summary_requested = true;
        paths.summary_path = absolute_path((fs::path(*options.output_directory) / "summary.csv").string());
    }
    if (options.write_replicates) {
        paths.replicates_requested = true;
        paths.replicates_path = absolute_path(
            options.replicates_file ? fs::path(*options.replicates_file).string()
                                    : (fs::path(*options.output_directory) / "replicates.csv").string());
    }
    if (paths.summary_requested && paths.replicates_requested)
        require(!(same_or_child(paths.summary_path, paths.replicates_path) &&
                  same_or_child(paths.replicates_path, paths.summary_path)),
                "summary and replicate destinations are identical or aliased");
    const std::string tree = protected_tree.empty() ? ECAI_CPP2_PROTECTED_LEGACY_TREE : protected_tree;
    if (paths.summary_requested) validate_destination(paths.summary_path, tree);
    if (paths.replicates_requested) validate_destination(paths.replicates_path, tree);
    return paths;
}

WriteResult write_outputs(const cli::CliOptions& options,
                          const OutputDocument& document,
                          const std::string& protected_tree) {
    validate_stage13_document(document);
    const OutputPaths paths = resolve_paths(options, protected_tree);
    const std::string summary = paths.summary_requested ? serialize_summary_csv(document) : "";
    const std::string replicates =
        paths.replicates_requested ? serialize_replicates_csv(document) : "";
    WriteResult result;
    result.paths = paths;
    result.summary_rows = document.summary_rows.size();
    result.replicate_rows = document.replicate_rows.size();
    if (!paths.requested) return result;

    std::vector<fs::path> destinations;
    std::vector<const std::string*> payloads;
    destinations.reserve(2);
    payloads.reserve(2);
    if (paths.summary_requested) {
        destinations.emplace_back(paths.summary_path);
        payloads.push_back(&summary);
    }
    if (paths.replicates_requested) {
        destinations.emplace_back(paths.replicates_path);
        payloads.push_back(&replicates);
    }

    // Complete all path discovery and reserve every ownership-tracking slot
    // before the first directory or file can be created.
    std::vector<fs::path> planned_parents;
    for (const auto& destination : destinations)
        record_missing_parents(destination.parent_path(), planned_parents);
    std::vector<fs::path> owned_parents;
    std::vector<fs::path> owned_temp_directories;
    std::vector<fs::path> payload_paths;
    std::vector<std::size_t> committed_indices;
    owned_parents.reserve(planned_parents.size());
    owned_temp_directories.reserve(destinations.size());
    payload_paths.reserve(destinations.size());
    committed_indices.reserve(destinations.size());

    try {
        for (const auto& parent : planned_parents) {
            std::error_code error;
            const bool created = fs::create_directory(parent, error);
            require(!error, "cannot create output parent directory " + parent.string());
            if (created) {
                owned_parents.push_back(parent);
            } else {
                require(fs::is_directory(parent, error) && !error,
                        "output parent path is not a directory " + parent.string());
            }
        }

        std::random_device random;
        for (std::size_t i = 0; i < destinations.size(); ++i) {
            fs::path temp_directory = reserve_temp_directory(destinations[i], random);
            owned_temp_directories.push_back(std::move(temp_directory));
            payload_paths.push_back(owned_temp_directories.back() / "payload");
            if (g_write_test_hook)
                g_write_test_hook("before_payload_write", i, owned_temp_directories.back());
            std::ofstream output(payload_paths.back(), std::ios::binary | std::ios::out);
            require(static_cast<bool>(output), "cannot create temporary payload");
            output.write(payloads[i]->data(), static_cast<std::streamsize>(payloads[i]->size()));
            require(static_cast<bool>(output), "failed writing temporary payload");
            output.flush();
            require(static_cast<bool>(output), "failed flushing temporary payload");
            if (g_write_test_hook)
                g_write_test_hook("before_payload_close", i, owned_temp_directories.back());
            output.close();
            require(!output.fail(), "failed closing temporary payload");
        }

        const std::string tree =
            protected_tree.empty() ? ECAI_CPP2_PROTECTED_LEGACY_TREE : protected_tree;
        for (std::size_t i = 0; i < destinations.size(); ++i) {
            if (g_write_test_hook) g_write_test_hook("before_commit", i, owned_temp_directories[i]);
            validate_destination(destinations[i], tree);
            std::error_code error;
            fs::create_hard_link(payload_paths[i], destinations[i], error);
            require(!error, "atomic no-replace hard-link commit failed: " + error.message());
            committed_indices.push_back(i);
        }
    } catch (...) {
        for (auto it = committed_indices.rbegin(); it != committed_indices.rend(); ++it) {
            const std::size_t i = *it;
            std::error_code error;
            const bool still_owned = fs::equivalent(destinations[i], payload_paths[i], error);
            if (!error && still_owned) {
                fs::remove(destinations[i], error);
            }
        }
        (void)cleanup_temporary_payloads(owned_temp_directories, payload_paths);
        for (auto it = owned_parents.rbegin(); it != owned_parents.rend(); ++it) {
            std::error_code ignored;
            fs::remove(*it, ignored);  // succeeds only when still empty
        }
        throw;
    }

    // Committed destinations are success even if best-effort removal of the
    // retained payload links cannot complete. Surface that state without
    // reporting a write failure that would contradict the filesystem.
    result.cleanup_complete = cleanup_temporary_payloads(owned_temp_directories, payload_paths);
    return result;
}

}  // namespace ecai_cpp2::output
