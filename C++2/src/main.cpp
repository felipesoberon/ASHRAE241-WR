// CLI entry point for ecai_cpp2 (Stage 13).
//
// Thin by construction: it parses the argument vector, resolves the run, runs
// the replicate batch exactly once, and prints one console line. Every rule
// lives in include/cli.h; every statistic lives in Stage 11. Stage 13 output
// composition and writing are delegated to output.cpp after the single run.

#include <algorithm>
#include <exception>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

#include "cli.h"
#include "simulation.h"
#include "output.h"

#ifndef ECAI_CPP2_PROTECTED_LEGACY_TREE
#define ECAI_CPP2_PROTECTED_LEGACY_TREE ""
#endif

namespace {

std::string program_name(int argc, char** argv) {
    if (argc > 0 && argv != nullptr && argv[0] != nullptr && argv[0][0] != '\0') {
        return argv[0];
    }
    return "ecai_cpp2";
}

// Console-only convenience marker (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md
// Phase 6): lets a human running the executable see at a glance which
// infector policy a run resolved to, without opening the CSV. Not part of
// any Stage 13 document/schema.
const char* infector_policy_name(ecai_cpp2::simulation::InfectorPolicy policy) {
    return policy == ecai_cpp2::simulation::InfectorPolicy::risk_averse ? "risk_averse"
                                                                        : "community";
}

}  // namespace

int main(int argc, char** argv) {
    const std::string program = program_name(argc, argv);
    try {
        std::vector<std::string> args;
        for (int index = 1; index < argc; ++index) {
            args.emplace_back(argv[index] != nullptr ? argv[index] : "");
        }

        const ecai_cpp2::cli::CliOptions options = ecai_cpp2::cli::parse_cli(args);
        if (options.help_requested) {
            // Sole-argument --help: no category lookup, no simulation.
            std::cout << ecai_cpp2::cli::help_text();
            return 0;
        }

        const ecai_cpp2::cli::ResolvedRun run = ecai_cpp2::cli::resolve_run(options);
        const ecai_cpp2::simulation::ReplicateBatch batch = ecai_cpp2::simulation::run_replicates(
            run.base_seed, run.config, *run.category, run.overrides);

        const ecai_cpp2::output::OutputDocument document =
            ecai_cpp2::output::build_document(run, batch, options,
                                              ECAI_CPP2_PROTECTED_LEGACY_TREE);
        const ecai_cpp2::output::WriteResult written =
            ecai_cpp2::output::write_outputs(options, document,
                                             ECAI_CPP2_PROTECTED_LEGACY_TREE);
        std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
        if (run.scenario.mode == ecai_cpp2::simulation::OutputMode::probability) {
            // The probability-mode fields are pulled from the already-validated
            // document rather than recomputed here, so main.cpp stays a thin
            // presenter of Stage 13 output (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md
            // section 4.5).
            const auto field = [&](const char* name) -> const std::string& {
                const auto it = std::find(document.summary_columns.begin(),
                                          document.summary_columns.end(), std::string(name));
                const auto column = static_cast<std::size_t>(
                    std::distance(document.summary_columns.begin(), it));
                return document.summary_rows.at(0).at(column);
            };
            std::cout << "console_summary (not a CSV file) mode=probability category=\""
                      << run.config.category_name << "\" input_ecai_l_s_person="
                      << field("input_ecai_l_s_person") << " percentile=" << field("percentile")
                      << " raw_mean_probability=" << field("raw_mean_probability")
                      << " target_probability=" << field("target_probability")
                      << " comparator_status=" << field("comparator_status");
        } else {
            std::cout << "console_summary (not a CSV file) category=\""
                      << run.config.category_name
                      << "\" raw_mean_l_s_person=" << batch.summary.mean_l_s_person
                      << " ceiling_to_5_l_s_person=" << batch.summary.adjusted_ceiling_to_5_l_s_person;
        }
        std::cout << " infector_policy=" << infector_policy_name(run.scenario.infector_policy);
        if (written.paths.summary_requested) std::cout << " summary_file=\"" << written.paths.summary_path << "\"";
        if (written.paths.replicates_requested) std::cout << " replicates_file=\"" << written.paths.replicates_path << "\"";
        std::cout << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ecai_cpp2: error: " << error.what() << "\n"
                  << "Try '" << program << " --help' for usage.\n";
        return 1;
    }
}
