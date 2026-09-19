#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "output.h"
#include "seed_manager.h"

#ifndef ECAI_CPP2_EXECUTABLE
#define ECAI_CPP2_EXECUTABLE "ecai_cpp2"
#endif

namespace fs = std::filesystem;
using ecai_cpp2::simulation::ReplicateBatch;
static int checks = 0;
static int failures = 0;

void check(bool ok, const std::string& message) {
    ++checks;
    if (ok) std::cout << "PASS: " << message << '\n';
    else { ++failures; std::cerr << "FAIL: " << message << '\n'; }
}

template <class F> bool rejects(F&& f) {
    try { f(); } catch (const std::exception&) { return true; }
    return false;
}

std::string read_text(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    std::ostringstream text;
    text << input.rdbuf();
    return text.str();
}

std::string round_trip_text(double value) {
    std::ostringstream text;
    text << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return text.str();
}

std::vector<std::vector<std::string>> parse_csv(const std::string& text) {
    std::vector<std::vector<std::string>> rows(1);
    std::string field;
    bool quoted = false;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char ch = text[i];
        if (quoted) {
            if (ch == '"') {
                if (i + 1 < text.size() && text[i + 1] == '"') { field.push_back('"'); ++i; }
                else quoted = false;
            } else field.push_back(ch);
        } else if (ch == '"' && field.empty()) quoted = true;
        else if (ch == ',') { rows.back().push_back(field); field.clear(); }
        else if (ch == '\n') { rows.back().push_back(field); field.clear(); rows.emplace_back(); }
        else field.push_back(ch);
    }
    if (!field.empty() || !rows.back().empty()) rows.back().push_back(field);
    if (!rows.empty() && rows.back().empty()) rows.pop_back();
    if (quoted) throw std::runtime_error("unterminated CSV field");
    return rows;
}

std::map<std::string, std::string> named_row(const ecai_cpp2::output::OutputDocument& document) {
    std::map<std::string, std::string> result;
    for (std::size_t i = 0; i < document.summary_columns.size(); ++i)
        result.emplace(document.summary_columns[i], document.summary_rows.at(0).at(i));
    return result;
}

ReplicateBatch constant_batch(ReplicateBatch batch, double value, double confidence) {
    for (std::size_t i = 0; i < batch.replicates.size(); ++i) {
        auto& r = batch.replicates[i];
        r.replicate_index = i;
        r.infector_counts.assign(r.ecai_l_s_person.size(), 1);
        std::fill(r.ecai_l_s_person.begin(), r.ecai_l_s_person.end(), value);
        r.percentile_l_s_person = value;
        auto& c = r.counters;
        c.zero_infector_count = 0;
        c.active_count = c.n;
        c.k_total = c.n;
        c.negative_clamp_count = 0;
        c.output_zero_count = value == 0.0 ? c.n : 0;
        c.qer_evaluations = c.k_total;
        c.ecai_evaluations = c.active_count;
        c.susceptible_values_consumed = 3 * c.active_count;
        c.infector_values_consumed = 8 * c.k_total;
        for (std::size_t bank = 0; bank < 11; ++bank) {
            const std::size_t length = bank < 3 ? c.active_count : c.k_total;
            r.planned_bank_lengths[bank] = length;
            r.final_cursors[bank] = length;
            r.bank_generated[bank] = length > 0;
        }
        r.cursors_exhausted = true;
    }
    std::vector<double> raw(batch.replicates.size(), value);
    batch.summary = ecai_cpp2::statistics::summarize_replicates(raw, confidence);
    return batch;
}

// Probability-mode analogue of constant_batch(): every row in every replicate
// becomes an active row (exactly one infector) whose probability is `value`,
// so the batch remains internally consistent under build_document's Stage 6
// counter/bank reconciliation (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md
// section 5.1). Used only to exercise the comparator/saturation boundaries at
// a chosen mean, not to model a realistic sample.
ReplicateBatch constant_probability_batch(ReplicateBatch batch, double value, double confidence) {
    for (std::size_t i = 0; i < batch.replicates.size(); ++i) {
        auto& r = batch.replicates[i];
        r.replicate_index = i;
        r.infector_counts.assign(r.probability.size(), 1);
        std::fill(r.probability.begin(), r.probability.end(), value);
        r.percentile_probability = value;
        auto& c = r.counters;
        c.zero_infector_count = 0;
        c.active_count = c.n;
        c.k_total = c.n;
        c.saturated_count = value == 1.0 ? c.n : 0;
        c.output_zero_count = value == 0.0 ? c.n : 0;
        c.qer_evaluations = c.k_total;
        c.probability_evaluations = c.active_count;
        c.susceptible_values_consumed = 3 * c.active_count;
        c.infector_values_consumed = 8 * c.k_total;
        for (std::size_t bank = 0; bank < 11; ++bank) {
            const std::size_t length = bank < 3 ? c.active_count : c.k_total;
            r.planned_bank_lengths[bank] = length;
            r.final_cursors[bank] = length;
            r.bank_generated[bank] = length > 0;
        }
        r.cursors_exhausted = true;
    }
    std::vector<double> raw(batch.replicates.size(), value);
    batch.summary = ecai_cpp2::statistics::summarize_replicates(raw, confidence);
    return batch;
}

int main() {
    using namespace ecai_cpp2;
    using output::OutputDocument;

    // Real CSV parsing, including every character that triggers field quoting.
    OutputDocument quote_doc;
    quote_doc.summary_columns = {"plain", "quoted", "empty"};
    quote_doc.summary_rows = {{"a", "comma, quote \" CR\r LF\nend", ""}};
    const auto quote_rows = parse_csv(output::serialize_summary_csv(quote_doc));
    check(quote_rows.size() == 2 && quote_rows[1] == quote_doc.summary_rows[0],
          "CSV fields round-trip comma, quote, CR and LF with LF record endings");

    const cli::CliOptions options = cli::parse_cli(
        {"--category", "Exam room", "--simulations", "4", "--replicates", "2"});
    const cli::ResolvedRun run = cli::resolve_run(options);
    const ReplicateBatch batch = simulation::run_replicates(
        run.base_seed, run.config, *run.category, run.overrides);
    const OutputDocument document = output::build_document(run, batch, options);
    const std::vector<std::string> expected_summary_columns = {
        "schema_version","model_version","sampling_design_version","category","group","is_health_care_group",
        "source_i0","source_volume_m3","source_probability_susceptible_present",
        "source_people_outdoor_airflow_l_s_person","source_area_outdoor_airflow_l_s_m2",
        "source_baseline_outdoor_airflow_l_s_person","source_community_rate","source_mask_efficiency",
        "source_pbr_gm_m3_h","source_pbr_gsd","source_activity_breathing_pct","source_activity_talking_pct",
        "source_activity_vocalizing_pct","source_cdrop_gm_e4_per_m3","source_cdrop_gsd","source_d_gm_e6_m",
        "source_d_gsd","source_published_ecai_l_s_person","source_published_rr","effective_i0",
        "effective_volume_m3","effective_community_rate","effective_mask_efficiency","target_probability",
        "duration_hours","target_probability_over_duration","simulations_n","replicates_r","percentile",
        "confidence_level","base_seed","seed_policy_version","zero_infector_policy","negative_clamping_policy",
        "custom_scenario","comparator_eligible","published_comparator_l_s_person",
        "comparator_signed_difference_l_s_person","comparator_abs_difference_l_s_person","comparator_status",
        "dispersion_status","raw_mean_l_s_person","raw_degrees_of_freedom","raw_sd_l_s_person",
        "raw_se_l_s_person","raw_t_critical","raw_ci_lower_l_s_person","raw_ci_upper_l_s_person",
        "raw_min_l_s_person","raw_max_l_s_person","ceiling_to_5_l_s_person","summary_row_count",
        "replicate_csv_row_count","all_replicates_n","all_replicates_zero_infector_count",
        "all_replicates_active_count","all_replicates_k_total","all_replicates_negative_clamp_count",
        "all_replicates_output_zero_count","all_replicates_qer_evaluations","all_replicates_ecai_evaluations",
        "all_replicates_susceptible_values_consumed","all_replicates_infector_values_consumed",
        "counters_reconcile","all_cursors_exhausted","all_replicates_finite","banks_reconcile",
        "explicit_category","explicit_simulations","explicit_replicates","explicit_percentile",
        "explicit_confidence","explicit_target_probability","explicit_duration","explicit_base_seed",
        "explicit_community_rate","explicit_occupants","explicit_volume_m3","explicit_mask_efficiency",
        "explicit_output_directory","explicit_summary_file","explicit_replicates_file","explicit_write_replicates"};
    const std::vector<std::string> expected_replicate_columns = {
        "replicate_index","requested_percentile","raw_percentile_l_s_person","n","zero_infector_count",
        "active_count","k_total","negative_clamp_count","output_zero_count","qer_evaluations",
        "ecai_evaluations","susceptible_values_consumed","infector_values_consumed","planned_bank_lengths",
        "final_cursors","counters_reconcile","cursors_exhausted","finite_values","replicate_csv_row_count"};
    check(document.summary_columns == expected_summary_columns && document.summary_rows.size() == 1 &&
              document.summary_rows[0].size() == 89,
          "summary schema is frozen at exactly 89 ordered columns and one row");
    check(document.replicate_columns == expected_replicate_columns && document.replicate_rows.size() == 2,
          "replicate schema is frozen at exactly 19 ordered columns and R rows");
    check(document.summary_columns[19] == "source_cdrop_gm_e4_per_m3",
          "source c-drop units are unambiguous");
    check(std::find(document.summary_columns.begin(), document.summary_columns.end(),
                    "planned_bank_lengths") == document.summary_columns.end(),
          "summary does not present first-replicate bank data as aggregate metadata");
    const auto names = named_row(document);
    check(names.at("sampling_design_version") == std::to_string(seed_manager::kDesignVersion) &&
              names.at("seed_policy_version") == seed_manager::kSeedPolicyVersion,
          "output versions are emitted from their authoritative seed-manager owners");
    check(names.at("summary_row_count") == "1" && names.at("replicate_csv_row_count") == "2",
          "summary and replicate CSV row counts are explicit and unambiguous");

    const auto parsed_summary = parse_csv(output::serialize_summary_csv(document));
    const auto parsed_replicates = parse_csv(output::serialize_replicates_csv(document));
    check(parsed_summary.size() == 2 && parsed_summary[0].size() == 89 &&
              parsed_summary[1].size() == 89,
          "serialized summary retains exact header and row widths");
    check(parsed_replicates.size() == 3 && parsed_replicates[0].size() == 19 &&
              parsed_replicates[1].size() == 19 && parsed_replicates[2].size() == 19,
          "serialized replicate CSV contains exactly one header and R fixed-width rows");
    for (std::size_t i = 0; i < batch.replicates.size(); ++i) {
        const std::size_t column = std::distance(document.replicate_columns.begin(),
            std::find(document.replicate_columns.begin(), document.replicate_columns.end(),
                      "raw_percentile_l_s_person"));
        const double round_trip = std::stod(document.replicate_rows[i][column]);
        check(std::memcmp(&round_trip, &batch.replicates[i].percentile_l_s_person,
                          sizeof(double)) == 0,
              "max_digits10 replicate values round-trip bitwise");
    }

    auto expect_batch_reject = [&](const std::string& label,
                                   const std::function<void(ReplicateBatch&)>& mutate) {
        ReplicateBatch bad = batch;
        mutate(bad);
        check(rejects([&] { (void)output::build_document(run, bad, options); }), label);
    };

    expect_batch_reject("reordered stable replicate identity is rejected", [](auto& b) {
        std::swap(b.replicates[0], b.replicates[1]);
    });
    expect_batch_reject("tampered replicate identity is rejected", [](auto& b) {
        ++b.replicates[0].replicate_index;
    });
    expect_batch_reject("tampered raw ECAi value is rejected", [](auto& b) {
        b.replicates[0].ecai_l_s_person[0] = std::nextafter(
            b.replicates[0].ecai_l_s_person[0], std::numeric_limits<double>::infinity());
    });
    expect_batch_reject("negative raw ECAi value is rejected", [](auto& b) {
        b.replicates[0].ecai_l_s_person[0] = -1.0;
    });
    expect_batch_reject("nonfinite raw ECAi value is rejected", [](auto& b) {
        b.replicates[0].ecai_l_s_person[0] = std::numeric_limits<double>::infinity();
    });
    expect_batch_reject("wrong raw ECAi vector length is rejected", [](auto& b) {
        b.replicates[0].ecai_l_s_person.pop_back();
    });
    expect_batch_reject("stored replicate percentile mismatch is rejected", [](auto& b) {
        b.replicates[0].percentile_l_s_person = std::nextafter(
            b.replicates[0].percentile_l_s_person, std::numeric_limits<double>::infinity());
    });
    expect_batch_reject("stored requested percentile mismatch is rejected", [](auto& b) {
        b.replicates[0].percentile = 0.5;
    });

    const std::vector<std::function<void(ReplicateBatch&)>> counter_mutations = {
        [](auto& b){++b.replicates[0].counters.n;},
        [](auto& b){++b.replicates[0].counters.zero_infector_count;},
        [](auto& b){++b.replicates[0].counters.active_count;},
        [](auto& b){++b.replicates[0].counters.k_total;},
        [](auto& b){++b.replicates[0].counters.negative_clamp_count;},
        [](auto& b){++b.replicates[0].counters.output_zero_count;},
        [](auto& b){++b.replicates[0].counters.qer_evaluations;},
        [](auto& b){++b.replicates[0].counters.ecai_evaluations;},
        [](auto& b){++b.replicates[0].counters.susceptible_values_consumed;},
        [](auto& b){++b.replicates[0].counters.infector_values_consumed;}
    };
    for (std::size_t i = 0; i < counter_mutations.size(); ++i)
        expect_batch_reject("tampered replicate counter " + std::to_string(i) + " is rejected",
                            counter_mutations[i]);

    expect_batch_reject("tampered susceptible bank plan is rejected", [](auto& b) {
        ++b.replicates[0].planned_bank_lengths[0];
    });
    expect_batch_reject("tampered infector bank plan is rejected", [](auto& b) {
        ++b.replicates[0].planned_bank_lengths[3];
    });
    expect_batch_reject("tampered final bank cursor is rejected", [](auto& b) {
        ++b.replicates[0].final_cursors[0];
    });
    expect_batch_reject("tampered generated-bank flag is rejected", [](auto& b) {
        b.replicates[0].bank_generated[0] = !b.replicates[0].bank_generated[0];
    });
    expect_batch_reject("false cursor-exhaustion claim is rejected", [](auto& b) {
        b.replicates[0].cursors_exhausted = false;
    });
    expect_batch_reject("wrong authoritative infector-count vector length is rejected", [](auto& b) {
        b.replicates[0].infector_counts.pop_back();
    });
    expect_batch_reject("infector count above effective I0 is rejected", [&](auto& b) {
        b.replicates[0].infector_counts[0] = run.scenario.occupants + 1;
    });
    expect_batch_reject("zero-infector output negative zero is rejected", [](auto& b) {
        const auto zero = std::find(b.replicates[0].infector_counts.begin(),
                                    b.replicates[0].infector_counts.end(), 0);
        if (zero == b.replicates[0].infector_counts.end())
            throw std::runtime_error("fixture unexpectedly has no zero-infector row");
        const std::size_t row = static_cast<std::size_t>(
            std::distance(b.replicates[0].infector_counts.begin(), zero));
        b.replicates[0].ecai_l_s_person[row] = -0.0;
    });
    {
        ReplicateBatch forged = constant_batch(batch, 10.0, options.confidence_level);
        auto& r = forged.replicates[0];
        r.ecai_l_s_person[0] = 0.0;
        ++r.counters.zero_infector_count;
        --r.counters.active_count;
        ++r.counters.output_zero_count;
        --r.counters.ecai_evaluations;
        r.counters.susceptible_values_consumed -= 3;
        for (std::size_t bank = 0; bank < 3; ++bank) {
            --r.planned_bank_lengths[bank];
            --r.final_cursors[bank];
        }
        check(rejects([&] { (void)output::build_document(run, forged, options); }),
              "coherently forged active/zero/ECAi/susceptible-bank accounting is rejected by Stage 6 counts");
    }
    {
        ReplicateBatch forged = constant_batch(batch, 10.0, options.confidence_level);
        auto& r = forged.replicates[0];
        ++r.counters.k_total;
        ++r.counters.qer_evaluations;
        r.counters.infector_values_consumed += 8;
        for (std::size_t bank = 3; bank < 11; ++bank) {
            ++r.planned_bank_lengths[bank];
            ++r.final_cursors[bank];
        }
        check(rejects([&] { (void)output::build_document(run, forged, options); }),
              "coherently forged K/QER/infector-bank accounting is rejected by Stage 6 counts");
    }

    {
        cli::CliOptions unprovenanced = options;
        unprovenanced.simulations_explicit = false;
        check(rejects([&] { (void)output::build_document(run, batch, unprovenanced); }),
              "build_document revalidates unprovenanced scalar mutations through resolve_run");
    }
    {
        cli::CliOptions wrong_category = options;
        wrong_category.category_name = "Classroom";
        check(rejects([&] { (void)output::build_document(run, batch, wrong_category); }),
              "build_document rejects a category name whose authoritative resolved run differs");
    }
    {
        cli::CliOptions mismatched_custom = options;
        mismatched_custom.volume_m3 = run.category->volume_m3;
        check(rejects([&] { (void)output::build_document(run, batch, mismatched_custom); }),
              "build_document rejects a category override/custom-provenance mismatch");
    }

    const std::vector<std::function<void(ReplicateBatch&)>> summary_mutations = {
        [](auto& b){++b.summary.replicate_count;},
        [](auto& b){b.summary.replicate_estimates_l_s_person[0] =
            std::nextafter(b.summary.replicate_estimates_l_s_person[0], INFINITY);},
        [](auto& b){b.summary.confidence_level = std::nextafter(b.summary.confidence_level, 1.0);},
        [](auto& b){b.summary.mean_l_s_person = std::nextafter(b.summary.mean_l_s_person, INFINITY);},
        [](auto& b){b.summary.min_l_s_person = std::nextafter(b.summary.min_l_s_person, INFINITY);},
        [](auto& b){b.summary.max_l_s_person = std::nextafter(b.summary.max_l_s_person, INFINITY);},
        [](auto& b){b.summary.dispersion_status = statistics::DispersionStatus::unavailable_single_replicate;},
        [](auto& b){b.summary.degrees_of_freedom = std::nextafter(*b.summary.degrees_of_freedom, INFINITY);},
        [](auto& b){b.summary.sample_sd_l_s_person = std::nextafter(*b.summary.sample_sd_l_s_person, INFINITY);},
        [](auto& b){b.summary.standard_error_l_s_person = std::nextafter(*b.summary.standard_error_l_s_person, INFINITY);},
        [](auto& b){b.summary.t_critical = std::nextafter(*b.summary.t_critical, INFINITY);},
        [](auto& b){b.summary.ci_lower_l_s_person = std::nextafter(*b.summary.ci_lower_l_s_person, INFINITY);},
        [](auto& b){b.summary.ci_upper_l_s_person = std::nextafter(*b.summary.ci_upper_l_s_person, INFINITY);},
        [](auto& b){b.summary.adjusted_ceiling_to_5_l_s_person =
            std::nextafter(b.summary.adjusted_ceiling_to_5_l_s_person, INFINITY);}
    };
    for (std::size_t i = 0; i < summary_mutations.size(); ++i)
        expect_batch_reject("tampered summary field " + std::to_string(i) + " is rejected",
                            summary_mutations[i]);
    {
        ReplicateBatch ordered = constant_batch(batch, 10.0, options.confidence_level);
        std::fill(ordered.replicates[1].ecai_l_s_person.begin(),
                  ordered.replicates[1].ecai_l_s_person.end(), 15.0);
        ordered.replicates[1].percentile_l_s_person = 15.0;
        ordered.summary = statistics::summarize_replicates({10.0, 15.0}, options.confidence_level);
        (void)output::build_document(run, ordered, options);
        std::swap(ordered.summary.replicate_estimates_l_s_person[0],
                  ordered.summary.replicate_estimates_l_s_person[1]);
        check(rejects([&]{ (void)output::build_document(run, ordered, options); }),
              "summary raw replicate vector order is authoritative");
    }

    const cli::CliOptions r1_options = cli::parse_cli(
        {"--category", "Exam room", "--simulations", "1", "--replicates", "1"});
    const auto r1_run = cli::resolve_run(r1_options);
    const auto r1_batch = simulation::run_replicates(
        r1_run.base_seed, r1_run.config, *r1_run.category, r1_run.overrides);
    const auto r1_names = named_row(output::build_document(r1_run, r1_batch, r1_options));
    for (const std::string field : {"raw_degrees_of_freedom", "raw_sd_l_s_person",
                                    "raw_se_l_s_person", "raw_t_critical",
                                    "raw_ci_lower_l_s_person", "raw_ci_upper_l_s_person"})
        check(r1_names.at(field).empty(), "R=1 leaves named dispersion field blank: " + field);

    const double published = run.category->published_ecai_Lsp;
    for (const auto& expected : std::vector<std::pair<double, std::string>>{
             {published, "EXACT_MATCH"}, {published + 5.0, "NEAR_WITHIN_ONE_INCREMENT"},
             {published + 10.0, "FAIL_GT_ONE_INCREMENT"}}) {
        const auto comparator = named_row(output::build_document(
            run, constant_batch(batch, expected.first, options.confidence_level), options));
        const double difference = expected.first - published;
        check(comparator.at("comparator_eligible") == "true" &&
                  comparator.at("comparator_status") == expected.second &&
                  comparator.at("comparator_signed_difference_l_s_person") == round_trip_text(difference) &&
                  comparator.at("comparator_abs_difference_l_s_person") == round_trip_text(std::abs(difference)),
              "eligible comparator emits status " + expected.second);
    }
    auto custom_options = cli::parse_cli({"--category", "Exam room", "--simulations", "4",
                                          "--replicates", "2", "--volume-m3", "50"});
    auto custom_run = cli::resolve_run(custom_options);
    auto custom_batch = simulation::run_replicates(custom_run.base_seed, custom_run.config,
                                                    *custom_run.category, custom_run.overrides);
    const auto custom_names = named_row(output::build_document(custom_run, custom_batch, custom_options));
    check(custom_names.at("comparator_eligible") == "false" &&
              custom_names.at("comparator_status") == "NOT_APPLICABLE_CUSTOM" &&
              custom_names.at("published_comparator_l_s_person").empty(),
          "custom scenario suppresses comparator values and status");
    for (const double raw : {0.0, 10.0, std::nextafter(10.0, 0.0),
                             std::nextafter(10.0, INFINITY)}) {
        const ReplicateBatch boundary_batch = constant_batch(batch, raw, options.confidence_level);
        const auto boundary = named_row(output::build_document(run, boundary_batch, options));
        check(boundary.at("ceiling_to_5_l_s_person") ==
                  round_trip_text(boundary_batch.summary.adjusted_ceiling_to_5_l_s_person),
              "ceiling boundary is serialized from authoritative statistics output");
    }

    const fs::path root = fs::temp_directory_path() / "ecai_cpp2_stage13_output_test";
    std::error_code ec;
    fs::remove_all(root, ec);

    // resolve_paths is a public, pure path-resolution API. Exercise it
    // directly rather than obtaining its result only through write_outputs.
    const output::OutputPaths no_paths = output::resolve_paths(cli::CliOptions{});
    check(!no_paths.requested && !no_paths.summary_requested &&
              !no_paths.replicates_requested && no_paths.summary_path.empty() &&
              no_paths.replicates_path.empty(),
          "resolve_paths: default options request no output paths");

    const fs::path resolve_dir = root / "resolve-only";
    cli::CliOptions directory_paths;
    directory_paths.output_directory = resolve_dir.string();
    directory_paths.write_replicates = true;
    const output::OutputPaths resolved_directory = output::resolve_paths(directory_paths);
    check(resolved_directory.requested && resolved_directory.summary_requested &&
              resolved_directory.replicates_requested &&
              resolved_directory.summary_path ==
                  fs::absolute(resolve_dir / "summary.csv").lexically_normal().string() &&
              resolved_directory.replicates_path ==
                  fs::absolute(resolve_dir / "replicates.csv").lexically_normal().string() &&
              !fs::exists(root),
          "resolve_paths: output directory resolves exact default filenames without filesystem mutation");

    cli::CliOptions mixed_paths = directory_paths;
    mixed_paths.summary_file = (root / "explicit" / "named-summary.csv").string();
    const output::OutputPaths resolved_mixed = output::resolve_paths(mixed_paths);
    check(resolved_mixed.summary_path ==
              fs::absolute(*mixed_paths.summary_file).lexically_normal().string() &&
              resolved_mixed.replicates_path ==
                  fs::absolute(resolve_dir / "replicates.csv").lexically_normal().string() &&
              !fs::exists(root),
          "resolve_paths: explicit summary overrides the directory default while replicates use it");

    cli::CliOptions invalid_paths;
    invalid_paths.replicates_file = (root / "invalid" / "replicates.csv").string();
    check(rejects([&] { (void)output::resolve_paths(invalid_paths); }) && !fs::exists(root),
          "resolve_paths: replicates file without write flag is rejected without filesystem mutation");

    fs::create_directories(root, ec);
    const fs::path existing = root / "foreign.csv";
    { std::ofstream f(existing); f << "FOREIGN\n"; }
    auto existing_options = cli::parse_cli({"--category", "Exam room", "--summary-file", existing.string()});
    check(rejects([&]{ (void)output::write_outputs(existing_options, document); }) &&
              read_text(existing) == "FOREIGN\n",
          "pre-existing foreign destination is never overwritten or removed");

    auto expect_invalid_document_no_mutation = [&](const std::string& label,
                                                    const std::function<void(OutputDocument&)>& mutate) {
        const fs::path parent = root / ("invalid-" + std::to_string(checks));
        auto invalid_options = cli::parse_cli(
            {"--category", "Exam room", "--summary-file", (parent / "summary.csv").string(),
             "--replicates-file", (parent / "replicates.csv").string(), "--write-replicates"});
        OutputDocument bad = document;
        mutate(bad);
        check(rejects([&] { (void)output::write_outputs(invalid_options, bad); }) &&
                  !fs::exists(parent), label);
    };
    expect_invalid_document_no_mutation("writer rejects a wrong summary header before parent creation",
                                        [](auto& d) { d.summary_columns[0] = "wrong"; });
    expect_invalid_document_no_mutation("writer rejects zero summary rows before parent creation",
                                        [](auto& d) { d.summary_rows.clear(); });
    expect_invalid_document_no_mutation("writer rejects two summary rows before parent creation",
                                        [](auto& d) { d.summary_rows.push_back(d.summary_rows[0]); });
    expect_invalid_document_no_mutation("writer rejects a wrong replicate header before parent creation",
                                        [](auto& d) { d.replicate_columns[0] = "wrong"; });
    expect_invalid_document_no_mutation("writer rejects replicate row-count metadata before parent creation",
                                        [](auto& d) { d.summary_rows[0][58] = "3"; });
    expect_invalid_document_no_mutation("writer rejects the wrong number of replicate rows before parent creation",
                                        [](auto& d) { d.replicate_rows.pop_back(); });

    const fs::path race_dir = root / "race" / "nested";
    const fs::path first = race_dir / "summary.csv";
    const fs::path second = race_dir / "replicates.csv";
    auto race_options = cli::parse_cli({"--category", "Exam room", "--summary-file", first.string(),
                                        "--replicates-file", second.string(), "--write-replicates"});
    output::detail::set_write_test_hook([&](const std::string& stage, std::size_t index,
                                             const fs::path& temp_directory) {
        (void)temp_directory;
        if (stage == "before_commit" && index == 1) {
            std::error_code ignored;
            fs::remove(first, ignored);
            std::ofstream replacement(first);
            replacement << "FOREIGN-FIRST\n";
            std::ofstream foreign(second);
            foreign << "RACE-WINNER\n";
        }
    });
    check(rejects([&]{ (void)output::write_outputs(race_options, document); }),
          "deterministic second-destination race fails the two-file commit");
    output::detail::set_write_test_hook({});
    check(read_text(first) == "FOREIGN-FIRST\n" && read_text(second) == "RACE-WINNER\n",
          "rollback preserves foreign replacements at both destination paths");
    bool temp_left = false;
    for (const auto& entry : fs::directory_iterator(race_dir))
        if (entry.is_directory() && entry.path().filename().string().find(".tmp-") != std::string::npos)
            temp_left = true;
    check(!temp_left, "all owned temporary directories are cleaned after commit failure");

    const fs::path write_fail_parent = root / "write-fail" / "owned";
    const fs::path write_fail_dest = write_fail_parent / "summary.csv";
    auto write_fail_options = cli::parse_cli(
        {"--category", "Exam room", "--summary-file", write_fail_dest.string()});
    output::detail::set_write_test_hook([](const std::string& stage, std::size_t,
                                            const fs::path& temp_directory) {
        if (stage == "before_payload_write") {
            std::error_code ignored;
            fs::remove_all(temp_directory, ignored);
        }
    });
    check(rejects([&]{ (void)output::write_outputs(write_fail_options, document); }),
          "deterministic payload open/write failure is reported");
    output::detail::set_write_test_hook({});
    check(!fs::exists(write_fail_dest) && !fs::exists(root / "write-fail"),
          "write failure leaves no destination, temp directory, or owned empty parent");

    const fs::path pre_close_dest = root / "pre-close-exception" / "summary.csv";
    auto pre_close_options = cli::parse_cli(
        {"--category", "Exam room", "--summary-file", pre_close_dest.string()});
    output::detail::set_write_test_hook([](const std::string& stage, std::size_t,
                                            const fs::path&) {
        if (stage == "before_payload_close")
            throw std::runtime_error("injected pre-close exception");
    });
    check(rejects([&]{ (void)output::write_outputs(pre_close_options, document); }),
          "deterministic injected pre-close exception is reported");
    output::detail::set_write_test_hook({});
    check(!fs::exists(pre_close_dest) && !fs::exists(root / "pre-close-exception"),
          "injected pre-close exception leaves no destination, temp directory, or owned empty parent");

    const fs::path protected_tree = root / "protected";
    fs::create_directories(protected_tree / "child", ec);
    const fs::path alias = root / "protected-link";
    fs::create_directory_symlink(protected_tree, alias, ec);
    for (const fs::path& candidate : {protected_tree / "exact.csv", protected_tree / "child" / "x.csv",
                                     protected_tree / "child" / ".." / "dotdot.csv", alias / "link.csv"}) {
        auto blocked = cli::parse_cli({"--category", "Exam room", "--summary-file", candidate.string()});
        check(rejects([&]{ (void)output::write_outputs(blocked, document, protected_tree.string()); }),
              "protected exact/descendant/dotdot/symlink destination is rejected");
    }
    auto aliases = cli::parse_cli({"--category", "Exam room", "--summary-file", (root/"same.csv").string(),
                                   "--replicates-file", (root/"."/"same.csv").string(), "--write-replicates"});
    check(rejects([&]{ (void)output::write_outputs(aliases, document); }),
          "identical lexical aliases are rejected before filesystem mutation");
    const fs::path real_alias_parent = root / "real-alias-parent";
    const fs::path link_alias_parent = root / "link-alias-parent";
    fs::create_directories(real_alias_parent, ec);
    fs::create_directory_symlink(real_alias_parent, link_alias_parent, ec);
    auto symlink_aliases = cli::parse_cli(
        {"--category", "Exam room", "--summary-file", (real_alias_parent/"same.csv").string(),
         "--replicates-file", (link_alias_parent/"same.csv").string(), "--write-replicates"});
    check(rejects([&]{ (void)output::write_outputs(symlink_aliases, document); }),
          "identical symlink aliases are rejected before filesystem mutation");

    const fs::path success_dir = root / "success";
    auto success_options = cli::parse_cli(
        {"--category", "Exam room", "--summary-file", (success_dir/"summary.csv").string(),
         "--replicates-file", (success_dir/"replicates.csv").string(), "--write-replicates"});
    const output::WriteResult success = output::write_outputs(success_options, document);
    const auto disk_summary = parse_csv(read_text(success_dir / "summary.csv"));
    const auto disk_replicates = parse_csv(read_text(success_dir / "replicates.csv"));
    check(success.cleanup_complete && disk_summary.size() == 2 && disk_summary[0] == expected_summary_columns &&
              disk_replicates.size() == 3 && disk_replicates[0] == expected_replicate_columns,
          "successful writer returns cleanup complete and exact 89/19 headers with one/R rows on disk");
    const auto disk_names = [&] {
        std::map<std::string, std::string> values;
        for (std::size_t i = 0; i < disk_summary[0].size(); ++i)
            values.emplace(disk_summary[0][i], disk_summary[1][i]);
        return values;
    }();
    check(disk_names.at("source_i0") == std::to_string(run.category->i0) &&
              disk_names.at("effective_i0") == std::to_string(run.scenario.occupants) &&
              disk_names.at("explicit_simulations") == "true" &&
              disk_names.at("raw_mean_l_s_person") == names.at("raw_mean_l_s_person") &&
              disk_names.at("counters_reconcile") == "true",
          "on-disk summary retains representative source/effective/provenance/stat/counter fields");
    bool success_temp_left = false;
    for (const auto& entry : fs::directory_iterator(success_dir))
        if (entry.path().filename().string().find(".tmp-") != std::string::npos) success_temp_left = true;
    check(!success_temp_left, "successful writer leaves no temporary directory");

    // --- Probability mode (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md section 5) ---
    const cli::CliOptions prob_options = cli::parse_cli(
        {"--category", "Exam room", "--simulations", "6", "--replicates", "2", "--mode",
         "probability"});
    const cli::ResolvedRun prob_run = cli::resolve_run(prob_options);
    const ReplicateBatch prob_batch = simulation::run_replicates(
        prob_run.base_seed, prob_run.config, *prob_run.category, prob_run.overrides);
    const OutputDocument prob_document = output::build_document(prob_run, prob_batch, prob_options);
    check(prob_document.summary_columns.size() == 96 && prob_document.summary_rows.size() == 1 &&
              prob_document.summary_rows[0].size() == 96,
          "probability summary schema is frozen at exactly 96 ordered columns and one row");
    check(prob_document.replicate_columns.size() == 20 && prob_document.replicate_rows.size() == 2 &&
              prob_document.replicate_rows[0].size() == 20,
          "probability replicate schema is frozen at exactly 20 ordered columns and R rows");
    const auto prob_names = named_row(prob_document);
    check(prob_names.at("schema_version") == "stage13-prob-v1" &&
              prob_names.at("mode") == "probability",
          "probability document carries the probability schema version and mode marker");
    check(prob_names.at("comparator_eligible") == "true" &&
              (prob_names.at("comparator_status") == "AT_TARGET_WITHIN_10PCT" ||
               prob_names.at("comparator_status") == "BELOW_TARGET" ||
               prob_names.at("comparator_status") == "ABOVE_TARGET"),
          "an uncustomized probability run reports an eligible comparator with a valid status");
    check(!prob_names.at("raw_mean_probability").empty() &&
              std::stod(prob_names.at("raw_min_probability")) >= 0.0 &&
              std::stod(prob_names.at("raw_max_probability")) <= 1.0,
          "probability raw statistics are bounded to the closed unit interval");
    check(prob_names.at("input_ecai_source") == "published_table4" &&
              prob_names.at("input_ecai_l_s_person") ==
                  round_trip_text(prob_run.category->published_ecai_Lsp),
          "an unoverridden probability run resolves the published Table 4 ECAi as its input");

    const auto prob_parsed_summary = parse_csv(output::serialize_summary_csv(prob_document));
    const auto prob_parsed_replicates = parse_csv(output::serialize_replicates_csv(prob_document));
    check(prob_parsed_summary.size() == 2 && prob_parsed_summary[0].size() == 96 &&
              prob_parsed_summary[1].size() == 96,
          "serialized probability summary retains exact header and row widths");
    check(prob_parsed_replicates.size() == 3 && prob_parsed_replicates[0].size() == 20 &&
              prob_parsed_replicates[1].size() == 20 && prob_parsed_replicates[2].size() == 20,
          "serialized probability replicate CSV contains exactly one header and R fixed-width rows");

    {
        const cli::CliOptions custom_prob_options = cli::parse_cli(
            {"--category", "Exam room", "--simulations", "4", "--replicates", "2", "--mode",
             "probability", "--input-ecai", "40"});
        const cli::ResolvedRun custom_prob_run = cli::resolve_run(custom_prob_options);
        const ReplicateBatch custom_prob_batch = simulation::run_replicates(
            custom_prob_run.base_seed, custom_prob_run.config, *custom_prob_run.category,
            custom_prob_run.overrides);
        const auto custom_prob_names = named_row(
            output::build_document(custom_prob_run, custom_prob_batch, custom_prob_options));
        check(custom_prob_names.at("input_ecai_source") == "explicit" &&
                  custom_prob_names.at("input_ecai_l_s_person") == round_trip_text(40.0) &&
                  custom_prob_names.at("comparator_eligible") == "false" &&
                  custom_prob_names.at("comparator_status") == "NOT_APPLICABLE_CUSTOM" &&
                  custom_prob_names.at("comparator_target_probability").empty(),
              "explicit --input-ecai marks the run custom and suppresses the comparator");
    }
    {
        ReplicateBatch tampered = prob_batch;
        tampered.replicates[0].mode = simulation::OutputMode::ecai;
        check(rejects([&] { (void)output::build_document(prob_run, tampered, prob_options); }),
              "build_probability_document rejects a replicate that did not resolve to probability mode");
    }
    {
        const cli::CliOptions ecai_options = cli::parse_cli(
            {"--category", "Exam room", "--simulations", "4", "--replicates", "2"});
        const cli::ResolvedRun ecai_run = cli::resolve_run(ecai_options);
        ReplicateBatch tampered = simulation::run_replicates(
            ecai_run.base_seed, ecai_run.config, *ecai_run.category, ecai_run.overrides);
        tampered.replicates[0].mode = simulation::OutputMode::probability;
        check(rejects([&] { (void)output::build_document(ecai_run, tampered, ecai_options); }),
              "ECAi document construction rejects a replicate that resolved to probability mode");
    }
    {
        ReplicateBatch bad = prob_batch;
        auto& counts = bad.replicates[0].infector_counts;
        const auto zero_it = std::find(counts.begin(), counts.end(), 0LL);
        check(zero_it != counts.end(), "probability fixture has a zero-infector row");
        const std::size_t zero_row = static_cast<std::size_t>(std::distance(counts.begin(), zero_it));
        bad.replicates[0].probability[zero_row] = -0.0;
        check(rejects([&] { (void)output::build_document(prob_run, bad, prob_options); }),
              "probability document rejects a non-canonical zero-infector output");
    }
    {
        ReplicateBatch bad = prob_batch;
        auto& counts = bad.replicates[0].infector_counts;
        const auto active_it = std::find_if(counts.begin(), counts.end(),
                                            [](long long c) { return c > 0; });
        check(active_it != counts.end(), "probability fixture has an active row");
        const std::size_t active_row =
            static_cast<std::size_t>(std::distance(counts.begin(), active_it));
        bad.replicates[0].probability[active_row] = 2.0;
        check(rejects([&] { (void)output::build_document(prob_run, bad, prob_options); }),
              "probability document rejects a probability above 1");
    }
    {
        ReplicateBatch bad = prob_batch;
        ++bad.replicates[0].counters.saturated_count;
        check(rejects([&] { (void)output::build_document(prob_run, bad, prob_options); }),
              "probability document rejects a tampered saturated_count counter");
    }

    const fs::path prob_dir = root / "probability-success";
    const auto prob_write_options = cli::parse_cli(
        {"--category", "Exam room", "--simulations", "6", "--replicates", "2", "--mode",
         "probability", "--summary-file", (prob_dir / "summary.csv").string(), "--replicates-file",
         (prob_dir / "replicates.csv").string(), "--write-replicates"});
    const output::WriteResult prob_written = output::write_outputs(prob_write_options, prob_document);
    const auto prob_disk_summary = parse_csv(read_text(prob_dir / "summary.csv"));
    const auto prob_disk_replicates = parse_csv(read_text(prob_dir / "replicates.csv"));
    check(prob_written.cleanup_complete && prob_disk_summary.size() == 2 &&
              prob_disk_summary[0].size() == 96 && prob_disk_summary[1].size() == 96 &&
              prob_disk_replicates.size() == 3 && prob_disk_replicates[0].size() == 20,
          "probability writer returns cleanup complete and exact 96/20 headers with one/R rows on disk");

    // --- Probability mode: aggregate counter invariants on a real batch -----
    check(prob_names.at("all_replicates_active_count") ==
                  std::to_string(std::stoull(prob_names.at("all_replicates_n")) -
                                  std::stoull(prob_names.at("all_replicates_zero_infector_count"))),
          "probability aggregate active + zero-infector counts partition all_replicates_n");
    check(std::stoull(prob_names.at("all_replicates_output_zero_count")) >=
              std::stoull(prob_names.at("all_replicates_zero_infector_count")),
          "probability aggregate output-zero count is at least the zero-infector count");
    check(std::stoull(prob_names.at("all_replicates_saturated_count")) <=
              std::stoull(prob_names.at("all_replicates_active_count")),
          "probability aggregate saturated count never exceeds the active count");
    check(std::stoull(prob_names.at("all_replicates_qer_evaluations")) ==
              std::stoull(prob_names.at("all_replicates_k_total")),
          "probability aggregate QER evaluations equal K_total");
    check(std::stoull(prob_names.at("all_replicates_probability_evaluations")) ==
              std::stoull(prob_names.at("all_replicates_active_count")),
          "probability aggregate evaluation count equals the active count");
    check(prob_names.at("counters_reconcile") == "true" &&
              prob_names.at("all_cursors_exhausted") == "true" &&
              prob_names.at("all_replicates_finite") == "true" &&
              prob_names.at("banks_reconcile") == "true",
          "probability document reports every counter/cursor/finiteness/bank invariant as true");
    check(prob_names.at("probability_policy") == "forward_wells_riley_expm1" &&
              prob_names.at("saturation_policy") == "probability_one_is_numerical_saturation" &&
              prob_names.at("zero_infector_policy") == "zero_infector_short_circuit",
          "probability document reports the frozen policy literals");

    // --- Probability mode: arbitrary --percentile ----------------------------
    {
        const cli::CliOptions pct_options = cli::parse_cli(
            {"--category", "Exam room", "--simulations", "6", "--replicates", "2", "--mode",
             "probability", "--percentile", "0.5"});
        const cli::ResolvedRun pct_run = cli::resolve_run(pct_options);
        const ReplicateBatch pct_batch = simulation::run_replicates(
            pct_run.base_seed, pct_run.config, *pct_run.category, pct_run.overrides);
        const OutputDocument pct_document = output::build_document(pct_run, pct_batch, pct_options);
        const auto pct_names = named_row(pct_document);
        check(pct_names.at("percentile") == round_trip_text(0.5),
              "an arbitrary --percentile is reported exactly in the probability summary");
        for (const auto& row : pct_document.replicate_rows) {
            const std::size_t column = std::distance(
                pct_document.replicate_columns.begin(),
                std::find(pct_document.replicate_columns.begin(),
                          pct_document.replicate_columns.end(), "requested_percentile"));
            check(row[column] == round_trip_text(0.5),
                  "every probability replicate row echoes the arbitrary requested percentile");
        }
    }

    // --- Probability mode: comparator boundary classification ---------------
    // Means chosen exactly per docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md
    // section 9.3 T9 against the default target_probability = 0.001: 0.9 and
    // 1.1 are the inclusive AT_TARGET_WITHIN_10PCT ratio boundaries.
    for (const auto& expected : std::vector<std::pair<double, std::string>>{
             {0.0009, "AT_TARGET_WITHIN_10PCT"}, {0.0011, "AT_TARGET_WITHIN_10PCT"},
             {0.00089, "BELOW_TARGET"}, {0.00111, "ABOVE_TARGET"}, {0.0, "BELOW_TARGET"}}) {
        const ReplicateBatch boundary_batch =
            constant_probability_batch(prob_batch, expected.first, prob_options.confidence_level);
        const auto boundary = named_row(output::build_document(prob_run, boundary_batch, prob_options));
        const double target = std::stod(boundary.at("target_probability"));
        const double ratio = expected.first / target;
        check(boundary.at("comparator_eligible") == "true" &&
                  boundary.at("comparator_status") == expected.second &&
                  boundary.at("comparator_ratio") == round_trip_text(ratio),
              "probability comparator classifies mean " + round_trip_text(expected.first) +
                  " as " + expected.second);
    }

    // --- Probability mode: saturation counters flow into aggregate totals ---
    {
        const ReplicateBatch saturated_batch =
            constant_probability_batch(prob_batch, 1.0, prob_options.confidence_level);
        const auto saturated_names =
            named_row(output::build_document(prob_run, saturated_batch, prob_options));
        check(saturated_names.at("all_replicates_saturated_count") ==
                  saturated_names.at("all_replicates_active_count") &&
                  saturated_names.at("all_replicates_active_count") ==
                      saturated_names.at("all_replicates_n") &&
                  std::stod(saturated_names.at("raw_max_probability")) == 1.0,
              "a fully saturated batch reports saturated_count equal to every active row");
        for (const auto& row : output::build_document(prob_run, saturated_batch, prob_options)
                                    .replicate_rows) {
            const auto& columns =
                output::build_document(prob_run, saturated_batch, prob_options).replicate_columns;
            const std::size_t saturated_column = std::distance(
                columns.begin(), std::find(columns.begin(), columns.end(), "saturated_count"));
            const std::size_t active_column = std::distance(
                columns.begin(), std::find(columns.begin(), columns.end(), "active_count"));
            check(row[saturated_column] == row[active_column],
                  "each fully saturated replicate row reports saturated_count == active_count");
        }
    }

    // --- Package 4a: direct-QER provenance/schema variants -------------------
    {
        const cli::CliOptions direct_options = cli::parse_cli(
            {"--category", "Exam room", "--simulations", "4", "--replicates", "2",
             "--qer-distribution", "jones-fitted"});
        const cli::ResolvedRun direct_run = cli::resolve_run(direct_options);
        const ReplicateBatch direct_batch = simulation::run_replicates(
            direct_run.base_seed, direct_run.config, *direct_run.category, direct_run.overrides);
        const OutputDocument direct_document =
            output::build_document(direct_run, direct_batch, direct_options);
        const auto direct_names = named_row(direct_document);
        check(direct_document.summary_columns.size() == 98 &&
                  direct_document.replicate_columns.size() == 28,
              "Jones-fitted community ECAi uses the versioned 98/28-column direct schema");
        check(direct_names.at("schema_version") == output::kDirectOutputSchemaVersion &&
                  direct_names.at("qer_family") == "jones-fitted" &&
                  direct_names.at("qer_profile").empty() &&
                  direct_names.at("qer_units") == "quanta/h per infector" &&
                  direct_names.at("qer_transform") ==
                      "10**(mu_log10+sigma_log10*normal_ppf(u))" &&
                  direct_names.at("qer_parameter_table_version") ==
                      qer_direct::kJonesFittedTableVersion &&
                  direct_names.at("qer_direct_bank") ==
                      "direct-qer-bank-v1;stream=direct_qer;layout=K_total" &&
                  direct_names.at("qer_seed_policy") == seed_manager::kSeedPolicyVersion,
              "Jones-fitted metadata identifies family, transform, units, table, bank and seed policy");
        check(direct_names.at("comparator_eligible") == "false" &&
                  direct_names.at("comparator_status") == "NOT_APPLICABLE_QER_FAMILY",
              "direct-family ECAi suppresses the full-Jones comparator");
        for (const auto& row : direct_document.replicate_rows) {
            check(row.back() == seed_manager::kSeedPolicyVersion,
                  "direct replicate rows carry the same run-level QER provenance");
        }
    }
    {
        const cli::CliOptions direct_options = cli::parse_cli(
            {"--category", "Exam room", "--simulations", "4", "--replicates", "2", "--mode",
             "probability", "--qer-distribution", "mikszewski-sars-cov-2-standing-speaking",
             "--infector-policy", "risk-averse"});
        const cli::ResolvedRun direct_run = cli::resolve_run(direct_options);
        const ReplicateBatch direct_batch = simulation::run_replicates(
            direct_run.base_seed, direct_run.config, *direct_run.category, direct_run.overrides);
        const OutputDocument direct_document =
            output::build_document(direct_run, direct_batch, direct_options);
        const auto direct_names = named_row(direct_document);
        check(direct_document.summary_columns.size() == 112 &&
                  direct_document.replicate_columns.size() == 31,
              "Mikszewski risk-averse probability uses the versioned direct schema");
        check(direct_names.at("schema_version") == output::kDirectRiskAverseProbabilitySchemaVersion &&
                  direct_names.at("qer_family") == "mikszewski" &&
                  direct_names.at("qer_profile") ==
                      "mikszewski-sars-cov-2-standing-speaking" &&
                  direct_names.at("qer_parameter_table_version") ==
                      qer_direct::kMikszewskiTableVersion &&
                  direct_names.at("comparator_eligible") == "false" &&
                  direct_names.at("comparator_status") == "NOT_APPLICABLE_QER_FAMILY",
              "Mikszewski risk-averse metadata carries canonical profile and suppresses comparator");
    }
    {
        const cli::CliOptions direct_options = cli::parse_cli(
            {"--category", "Exam room", "--simulations", "2", "--replicates", "1",
             "--qer-distribution", "jones-fitted"});
        const cli::ResolvedRun direct_run = cli::resolve_run(direct_options);
        const ReplicateBatch direct_batch = simulation::run_replicates(
            direct_run.base_seed, direct_run.config, *direct_run.category, direct_run.overrides);
        cli::CliOptions mismatched = direct_options;
        mismatched.qer_family = qer_direct::QerFamily::Mikszewski;
        mismatched.qer_mikszewski_profile =
            qer_direct::parse_mikszewski_profile("mikszewski-sars-cov-2-resting");
        check(rejects([&] { (void)output::build_document(direct_run, direct_batch, mismatched); }),
              "scenario/options QER family-profile mismatch is rejected before output composition");
    }

    fs::remove_all(root, ec);
    if (failures) {
        std::cerr << failures << " of " << checks << " checks failed.\n";
        return 1;
    }
    std::cout << "All " << checks << " checks passed.\n";
    return 0;
}
