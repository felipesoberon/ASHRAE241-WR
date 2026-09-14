// Stage 12 tests for cli.h/cli.cpp and the CLI entry point in src/main.cpp,
// per IMPLEMENTATION_PLAN.md Section 7, Stage 12 and sections 5.1/5.2, and
// docs/MODEL_SPECIFICATION.md D7 (category naming/normalization) and section
// 12 (frozen conventions).
//
// Every expected value below is taken from the frozen defaults in
// IMPLEMENTATION_PLAN.md section 5.1, from the published category table, or
// from an existing frozen constant (statistics::kMaxSupportedReplicateCount,
// sampling_plan::ResourcePolicy::max_simulations); never read back from this
// implementation.

#include "cli.h"

#include <sys/stat.h>
#include <sys/wait.h>

#include <array>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "category_data.h"
#include "sampling_plan.h"
#include "simulation.h"
#include "statistics.h"

#ifndef ECAI_CPP2_EXECUTABLE
#error "ECAI_CPP2_EXECUTABLE must be defined by the build system"
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

using Args = std::vector<std::string>;

// A small, bounded scenario used everywhere a parse must also be resolvable:
// "Exam room" has I0 = 3 and VOL = 41 m^3 (Jones et al. Table 2).
Args bounded_args() {
    return Args{"--category", "Exam room", "--simulations", "8", "--replicates", "2"};
}

// Every rejection must be a std::invalid_argument carrying one concise line.
void expect_reject(const Args& args, const std::string& description) {
    try {
        ecai_cpp2::cli::parse_cli(args);
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        check(!message.empty() && message.find('\n') == std::string::npos,
              "single-line diagnostic: " + description);
        return;
    } catch (...) {
        check(false, "rejects with std::invalid_argument: " + description);
        return;
    }
    check(false, "rejects " + description);
}

void expect_accept(const Args& args, const std::string& description) {
    try {
        ecai_cpp2::cli::parse_cli(args);
    } catch (const std::exception& error) {
        check(false, "accepts " + description + " (got: " + error.what() + ")");
        return;
    }
    check(true, "accepts " + description);
}

// bounded_args() plus one flag/value pair.
Args with_flag(const std::string& flag, const std::string& value) {
    Args args = bounded_args();
    args.push_back(flag);
    args.push_back(value);
    return args;
}

// A category plus exactly one flag/value pair, for flags bounded_args()
// already supplies (--simulations, --replicates), which may not be repeated.
Args only_flag(const std::string& flag, const std::string& value) {
    return Args{"--category", "Exam room", flag, value};
}

// --- Executable command matrix helpers ------------------------------------

struct CommandResult {
    int exit_code = -1;
    std::string output;
};

CommandResult run_command(const std::string& args) {
    const std::string command = std::string("\"") + ECAI_CPP2_EXECUTABLE + "\" " + args + " 2>&1";
    FILE* pipe = popen(command.c_str(), "r");
    if (pipe == nullptr) {
        std::cerr << "Failed to launch: " << command << "\n";
        std::exit(1);
    }
    std::string output;
    std::array<char, 4096> buffer{};
    std::size_t bytes_read = 0;
    while ((bytes_read = std::fread(buffer.data(), 1, buffer.size(), pipe)) > 0) {
        output.append(buffer.data(), bytes_read);
    }
    const int status = pclose(pipe);
    CommandResult result;
    result.output = output;
    if (WIFEXITED(status)) {
        result.exit_code = WEXITSTATUS(status);
    }
    return result;
}

std::size_t count_lines(const std::string& text) {
    std::size_t lines = 0;
    for (const char c : text) {
        if (c == '\n') ++lines;
    }
    return lines;
}

std::size_t count_occurrences(const std::string& text, const std::string& needle) {
    std::size_t count = 0;
    for (std::size_t at = text.find(needle); at != std::string::npos;
         at = text.find(needle, at + needle.size())) {
        ++count;
    }
    return count;
}

// Reads back the value printed after `label`, so the printed precision can be
// compared bit-for-bit with an in-process result.
double field_value(const std::string& text, const std::string& label) {
    const std::size_t at = text.find(label);
    if (at == std::string::npos) return std::numeric_limits<double>::quiet_NaN();
    const std::size_t start = at + label.size();
    std::size_t end = start;
    while (end < text.size() && text[end] != ' ' && text[end] != '\n') ++end;
    return std::strtod(text.substr(start, end - start).c_str(), nullptr);
}

bool path_exists(const std::string& path) {
    struct stat info {};
    return stat(path.c_str(), &info) == 0;
}

}  // namespace

int main() {
    using namespace ecai_cpp2;

    // --- Frozen defaults (IMPLEMENTATION_PLAN.md section 5.1) ---------------
    {
        const cli::CliOptions options = cli::parse_cli(Args{"--category", "Exam room"});
        check(!options.help_requested, "a category invocation does not request help");
        check(options.category_supplied, "--category is recorded as supplied");
        check(options.category_name == "Exam room", "--category preserves its string literally");
        check(options.simulations_n == 10000, "default simulations N is 10000");
        check(options.replicates_r == 500, "default replicates R is 500");
        check(options.percentile == 0.96, "default percentile is 0.96");
        check(options.confidence_level == 0.95, "default confidence level is 0.95");
        check(options.target_probability == 0.001, "default target probability is 0.001");
        check(options.duration_hours == 1.0, "default duration is 1 hour");
        check(options.base_seed == 485, "default base seed is 485");
        check(!options.occupants.has_value(), "occupants defaults to the category value");
        check(!options.volume_m3.has_value(), "volume defaults to the category value");
        check(!options.mask_efficiency.has_value(),
              "mask efficiency defaults to the category value");
        check(!options.community_infection_rate.has_value(),
              "community rate defaults to the category value");
        check(!options.output_directory.has_value(), "no output directory by default");
        check(!options.summary_file.has_value(), "no summary file by default");
        check(!options.replicates_file.has_value(), "no replicates file by default");
        check(!options.write_replicates, "replicate writing is off by default");
        check(!options.custom_scenario, "a bare category run is not a custom scenario");
        check(options.target_is_hourly, "the default 1 h duration makes the target hourly");
        check(!options.simulations_explicit && !options.replicates_explicit &&
                  !options.percentile_explicit && !options.confidence_explicit &&
                  !options.target_probability_explicit && !options.duration_explicit &&
                  !options.base_seed_explicit,
              "no per-field explicit flag is set for a bare category run");
    }
    {
        const cli::CliOptions options = cli::parse_cli(Args{"--category", "Patient"});
        check(options.category_name == "Patient",
              "a legacy alias is preserved literally by the parser");
    }

    // --- --help -------------------------------------------------------------
    {
        const cli::CliOptions options = cli::parse_cli(Args{"--help"});
        check(options.help_requested, "--help sets help_requested");
        check(!options.category_supplied,
              "--help parses successfully without a category, which is otherwise required");
    }
    {
        const std::string help = cli::help_text();
        check(!help.empty(), "help text is not empty");
        check(help.find("ecai_cpp2") != std::string::npos,
              "help text identifies the ecai_cpp2 executable");
        check(help.find("C++2") != std::string::npos,
              "help text identifies the C++2 implementation");

        // Bidirectional: every accepted flag is documented, and every
        // long-option token in the help text is an accepted flag.
        const std::vector<std::string> planned = {
            "--category",           "--simulations",     "--replicates",
            "--percentile",         "--confidence",      "--target-probability",
            "--duration-hours",     "--base-seed",       "--community-rate",
            "--occupants",          "--volume-m3",       "--mask-efficiency",
            "--mode",               "--input-ecai",      "--infector-policy",
            "--output-directory",   "--summary-file",    "--replicates-file",
            "--write-replicates",   "--help"};
        check(planned.size() == cli::accepted_flags().size(),
              "the parser accepts exactly the 20 planned flags");
        for (const std::string& flag : planned) {
            bool known = false;
            for (const std::string& parser_flag : cli::accepted_flags()) {
                if (parser_flag == flag) known = true;
            }
            check(known, "the parser accepts the planned flag " + flag);
            check(help.find(flag) != std::string::npos, "help text documents " + flag);
        }
        for (std::size_t i = 0; i + 1 < help.size(); ++i) {
            if (help[i] != '-' || help[i + 1] != '-') continue;
            if (i > 0 &&
                (help[i - 1] == '-' || std::isalnum(static_cast<unsigned char>(help[i - 1])) != 0)) {
                continue;
            }
            std::size_t end = i + 2;
            while (end < help.size() &&
                   (std::isalnum(static_cast<unsigned char>(help[end])) != 0 || help[end] == '-')) {
                ++end;
            }
            const std::string token = help.substr(i, end - i);
            bool known = false;
            for (const std::string& parser_flag : cli::accepted_flags()) {
                if (parser_flag == token) known = true;
            }
            check(known, "help text mentions only accepted flags (found " + token + ")");
            i = end - 1;
        }

        // Constraints stated in the help text must be the frozen ones.
        check(help.find(std::to_string(sampling_plan::ResourcePolicy{}.max_simulations)) !=
                  std::string::npos,
              "help text states the frozen maximum simulation count");
        check(help.find(std::to_string(statistics::kMaxSupportedReplicateCount)) !=
                  std::string::npos,
              "help text states the frozen maximum replicate count");
        check(help.find("0.96") != std::string::npos, "help text states the default percentile");
        check(help.find("0.95") != std::string::npos, "help text states the default confidence");
        check(help.find("0.001") != std::string::npos,
              "help text states the default target probability");
        check(help.find("10000") != std::string::npos, "help text states the default N");
        check(help.find("500") != std::string::npos, "help text states the default R");
        check(help.find("485") != std::string::npos, "help text states the default base seed");
    }

    // --- Every flag accepts a value and records its provenance --------------
    {
        const cli::CliOptions options = cli::parse_cli(Args{
            "--category", "exam-room", "--simulations", "8", "--replicates", "2",
            "--percentile", "0.9", "--confidence", "0.9", "--target-probability", "0.002",
            "--duration-hours", "2.5", "--base-seed", "17", "--community-rate", "0.05",
            "--occupants", "4", "--volume-m3", "50", "--mask-efficiency", "0.25",
            "--output-directory", "out dir", "--summary-file", "summary one.csv",
            "--replicates-file", "reps one.csv", "--write-replicates"});
        check(options.category_name == "exam-room", "delimiter-variant category text is preserved");
        check(options.simulations_n == 8, "--simulations is parsed");
        check(options.replicates_r == 2, "--replicates is parsed");
        check(options.percentile == 0.9, "--percentile is parsed");
        check(options.confidence_level == 0.9, "--confidence is parsed");
        check(options.target_probability == 0.002, "--target-probability is parsed");
        check(options.duration_hours == 2.5, "--duration-hours is parsed");
        check(options.base_seed == 17, "--base-seed is parsed");
        check(options.community_infection_rate.has_value() &&
                  *options.community_infection_rate == 0.05,
              "--community-rate is parsed as an override");
        check(options.occupants.has_value() && *options.occupants == 4,
              "--occupants is parsed as an override");
        check(options.volume_m3.has_value() && *options.volume_m3 == 50.0,
              "--volume-m3 is parsed as an override");
        check(options.mask_efficiency.has_value() && *options.mask_efficiency == 0.25,
              "--mask-efficiency is parsed as an override");
        check(options.output_directory.has_value() && *options.output_directory == "out dir",
              "--output-directory preserves a path containing a space");
        check(options.summary_file.has_value() && *options.summary_file == "summary one.csv",
              "--summary-file preserves a path containing a space");
        check(options.replicates_file.has_value() && *options.replicates_file == "reps one.csv",
              "--replicates-file preserves a path containing a space");
        check(options.write_replicates, "--write-replicates is a standalone boolean flag");
        check(options.simulations_explicit && options.replicates_explicit &&
                  options.percentile_explicit && options.confidence_explicit &&
                  options.target_probability_explicit && options.duration_explicit &&
                  options.base_seed_explicit,
              "every supplied scalar flag records an explicit-value provenance flag");
        check(options.custom_scenario, "explicit scientific overrides mark the run custom");
        check(!options.target_is_hourly, "a 2.5 h duration is never labelled an hourly target");
    }

    // --- custom_scenario provenance rules -----------------------------------
    {
        Args args = bounded_args();
        args.push_back("--base-seed");
        args.push_back("99");
        args.push_back("--output-directory");
        args.push_back("out");
        args.push_back("--write-replicates");
        args.push_back("--confidence");
        args.push_back("0.99");
        const cli::CliOptions options = cli::parse_cli(args);
        check(!options.custom_scenario,
              "N, R, base seed, confidence and output controls alone are not a custom scenario");
        check(options.target_is_hourly, "the untouched 1 h default remains an hourly target");
    }
    {
        Args args = bounded_args();
        args.push_back("--duration-hours");
        args.push_back("1");
        const cli::CliOptions options = cli::parse_cli(args);
        check(options.custom_scenario,
              "explicitly supplying the default duration still counts as an override");
        check(options.duration_explicit, "the explicit duration flag is recorded");
        check(options.target_is_hourly,
              "an explicit duration of exactly 1 h may still be called hourly");
    }
    {
        const std::vector<std::string> scientific = {"--target-probability", "--percentile",
                                                     "--community-rate", "--mask-efficiency"};
        for (const std::string& flag : scientific) {
            Args args = bounded_args();
            args.push_back(flag);
            args.push_back("0.5");
            check(cli::parse_cli(args).custom_scenario,
                  "explicit " + flag + " marks the run as a custom scenario");
        }
        {
            Args args = bounded_args();
            args.push_back("--occupants");
            args.push_back("2");
            check(cli::parse_cli(args).custom_scenario,
                  "explicit --occupants marks the run as a custom scenario");
        }
        {
            Args args = bounded_args();
            args.push_back("--volume-m3");
            args.push_back("30");
            check(cli::parse_cli(args).custom_scenario,
                  "explicit --volume-m3 marks the run as a custom scenario");
        }
    }

    // --- Output mode: --mode / --input-ecai parsing (additive; see
    // docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md section 9.3 T8) -----------
    {
        const cli::CliOptions options = cli::parse_cli(bounded_args());
        check(options.mode == simulation::OutputMode::ecai,
              "the default output mode is ecai");
        check(!options.mode_explicit, "the default output mode is not marked explicit");
        check(!options.input_ecai_l_s_person.has_value(),
              "there is no input ECAi override by default");
    }
    {
        const cli::CliOptions options = cli::parse_cli(with_flag("--mode", "probability"));
        check(options.mode == simulation::OutputMode::probability,
              "--mode probability parses");
        check(options.mode_explicit, "an explicit --mode is recorded as explicit");
        check(!options.custom_scenario, "--mode alone never marks the run a custom scenario");
    }
    {
        const cli::CliOptions options = cli::parse_cli(with_flag("--mode", "ecai"));
        check(options.mode == simulation::OutputMode::ecai,
              "an explicit --mode ecai parses");
        check(options.mode_explicit, "an explicit --mode ecai is still recorded as explicit");
        check(!options.custom_scenario, "an explicit --mode ecai does not mark the run custom");
    }
    expect_reject(with_flag("--mode", "Probability"),
                  "a wrong-case --mode value (case sensitive, no aliases)");
    expect_reject(with_flag("--mode", "ECAI"), "a wrong-case --mode value");
    expect_reject(with_flag("--mode", "foo"), "an unknown --mode value");
    expect_reject(with_flag("--mode", ""), "an empty --mode value");
    {
        Args args = bounded_args();
        args.push_back("--mode");
        args.push_back("probability");
        args.push_back("--mode");
        args.push_back("ecai");
        expect_reject(args, "a duplicate --mode flag");
    }
    expect_reject(with_flag("--input-ecai", "30"),
                  "--input-ecai without --mode probability");
    {
        Args args = with_flag("--mode", "probability");
        args.push_back("--input-ecai");
        args.push_back("30");
        const cli::CliOptions options = cli::parse_cli(args);
        check(options.input_ecai_l_s_person.has_value() && *options.input_ecai_l_s_person == 30.0,
              "--input-ecai parses alongside --mode probability");
        check(options.custom_scenario, "an explicit --input-ecai marks the run a custom scenario");
    }
    {
        Args args = with_flag("--mode", "probability");
        args.push_back("--input-ecai");
        args.push_back("-1");
        expect_reject(args, "a negative --input-ecai value");
    }
    {
        Args args = with_flag("--mode", "probability");
        args.push_back("--input-ecai");
        args.push_back("nan");
        expect_reject(args, "a non-finite --input-ecai value");
    }
    {
        Args args = with_flag("--mode", "probability");
        args.push_back("--input-ecai");
        args.push_back("0");
        const cli::CliOptions options = cli::parse_cli(args);
        check(options.input_ecai_l_s_person.has_value() && *options.input_ecai_l_s_person == 0.0,
              "--input-ecai accepts exactly 0");
    }

    // --- Infector policy: --infector-policy parsing (additive; see
    // docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md Phase 3) -----------------------
    {
        const cli::CliOptions options = cli::parse_cli(bounded_args());
        check(options.infector_policy == simulation::InfectorPolicy::community,
              "the default infector policy is community");
        check(!options.infector_policy_explicit,
              "the default infector policy is not marked explicit");
    }
    {
        const cli::CliOptions options =
            cli::parse_cli(with_flag("--infector-policy", "risk-averse"));
        check(options.infector_policy == simulation::InfectorPolicy::risk_averse,
              "--infector-policy risk-averse parses");
        check(options.infector_policy_explicit,
              "an explicit --infector-policy is recorded as explicit");
        check(!options.custom_scenario,
              "--infector-policy alone never marks the run a custom scenario");
    }
    {
        const cli::CliOptions options =
            cli::parse_cli(with_flag("--infector-policy", "community"));
        check(options.infector_policy == simulation::InfectorPolicy::community,
              "an explicit --infector-policy community parses");
        check(options.infector_policy_explicit,
              "an explicit --infector-policy community is still recorded as explicit");
        check(!options.custom_scenario,
              "an explicit --infector-policy community does not mark the run custom");
    }
    expect_reject(with_flag("--infector-policy", "Risk-Averse"),
                  "a wrong-case --infector-policy value (case sensitive, no aliases)");
    expect_reject(with_flag("--infector-policy", "COMMUNITY"),
                  "a wrong-case --infector-policy value");
    expect_reject(with_flag("--infector-policy", "risk_averse"),
                  "an underscore-spelled --infector-policy value (no aliases)");
    expect_reject(with_flag("--infector-policy", "foo"),
                  "an unknown --infector-policy value");
    expect_reject(with_flag("--infector-policy", ""), "an empty --infector-policy value");
    {
        Args args = bounded_args();
        args.push_back("--infector-policy");
        args.push_back("risk-averse");
        args.push_back("--infector-policy");
        args.push_back("community");
        expect_reject(args, "a duplicate --infector-policy flag");
    }
    // --- risk-averse requires a strictly positive community rate, checked
    // early for an explicit conflict --------------------------------------
    {
        Args args = only_flag("--infector-policy", "risk-averse");
        args.push_back("--community-rate");
        args.push_back("0");
        expect_reject(args,
                      "--infector-policy risk-averse combined with an explicit "
                      "--community-rate of exactly 0");
    }
    {
        Args args = only_flag("--infector-policy", "risk-averse");
        args.push_back("--community-rate");
        args.push_back("-0.01");
        expect_reject(args,
                      "--infector-policy risk-averse combined with a negative "
                      "--community-rate (rejected by the generic bound first)");
    }
    {
        Args args = only_flag("--infector-policy", "risk-averse");
        args.push_back("--community-rate");
        args.push_back("0.01");
        expect_accept(args,
                      "--infector-policy risk-averse with an explicit positive "
                      "--community-rate");
    }
    {
        // Exam room's published community rate (0.03) is already strictly
        // positive, so risk-averse resolves without any explicit
        // --community-rate override.
        const cli::ResolvedRun run =
            cli::resolve_run(cli::parse_cli(only_flag("--infector-policy", "risk-averse")));
        check(run.scenario.infector_policy == simulation::InfectorPolicy::risk_averse,
              "risk-averse resolves using the category's positive default community rate");
    }

    // --- Infector policy: resolved scenario storage and provenance ----------
    {
        const cli::ResolvedRun run = cli::resolve_run(cli::parse_cli(bounded_args()));
        check(run.scenario.infector_policy == simulation::InfectorPolicy::community,
              "a bare run resolves to the community infector policy");
        check(!run.scenario.infector_policy_explicit,
              "a bare run reports no explicit infector-policy override");
        check(run.overrides.infector_policy == simulation::InfectorPolicy::community,
              "the effective overrides record the community policy by default");
        check(!run.overrides.infector_policy_explicit,
              "the effective overrides record no explicit infector-policy provenance by default");
    }
    {
        const cli::ResolvedRun run =
            cli::resolve_run(cli::parse_cli(with_flag("--infector-policy", "risk-averse")));
        check(run.scenario.infector_policy == simulation::InfectorPolicy::risk_averse,
              "--infector-policy risk-averse resolves to the risk-averse scenario policy");
        check(run.scenario.infector_policy_explicit,
              "an explicit --infector-policy reaches the resolved scenario as explicit");
        check(run.overrides.infector_policy == simulation::InfectorPolicy::risk_averse &&
                  run.overrides.infector_policy_explicit,
              "the effective overrides carry the explicit risk-averse policy");
        check(!run.custom_scenario && !run.scenario.custom_scenario,
              "an infector-policy override alone is not a custom scenario (it is an estimand "
              "change, not a four-field scientific override)");
    }

    // --- Infector policy combines with both output modes --------------------
    {
        Args args = with_flag("--mode", "ecai");
        args.push_back("--infector-policy");
        args.push_back("risk-averse");
        const cli::ResolvedRun run = cli::resolve_run(cli::parse_cli(args));
        check(run.scenario.mode == simulation::OutputMode::ecai &&
                  run.scenario.infector_policy == simulation::InfectorPolicy::risk_averse,
              "risk-averse combines with --mode ecai");
    }
    {
        Args args = with_flag("--mode", "probability");
        args.push_back("--infector-policy");
        args.push_back("risk-averse");
        const cli::ResolvedRun run = cli::resolve_run(cli::parse_cli(args));
        check(run.scenario.mode == simulation::OutputMode::probability &&
                  run.scenario.infector_policy == simulation::InfectorPolicy::risk_averse,
              "risk-averse combines with --mode probability");
    }

    // --- Infector-policy tamper check, mirroring the mode tamper block above ---
    {
        cli::CliOptions options = cli::parse_cli(bounded_args());
        options.infector_policy = simulation::InfectorPolicy::risk_averse;
        bool rejected = false;
        try {
            cli::resolve_run(options);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check(rejected,
              "resolve_run rejects an unprovenanced change to the infector policy");
    }

    // --- Direct scenario construction rejects invalid risk-averse inputs,
    // even when the CLI is bypassed entirely --------------------------------
    {
        const cli::ResolvedRun base = cli::resolve_run(cli::parse_cli(bounded_args()));
        simulation::EffectiveOverrides overrides;
        overrides.infector_policy = simulation::InfectorPolicy::risk_averse;
        overrides.infector_policy_explicit = true;
        overrides.community_infection_rate = 0.0;
        simulation::EffectiveConfiguration config = base.config;
        config.community_infection_rate = 0.0;
        bool rejected = false;
        try {
            simulation::resolve_scenario(*base.category, config, overrides);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check(rejected,
              "direct resolve_scenario rejects risk-averse with community_infection_rate <= 0, "
              "bypassing the CLI");
    }
    {
        const cli::ResolvedRun base = cli::resolve_run(cli::parse_cli(bounded_args()));
        simulation::EffectiveOverrides overrides;
        overrides.infector_policy = simulation::InfectorPolicy::risk_averse;
        overrides.infector_policy_explicit = true;
        overrides.occupants = 0;
        simulation::EffectiveConfiguration config = base.config;
        config.i0 = 0;
        bool rejected = false;
        try {
            simulation::resolve_scenario(*base.category, config, overrides);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check(rejected,
              "direct resolve_scenario rejects risk-averse with occupants (I0) == 0, bypassing "
              "the CLI");
    }
    {
        // A valid direct construction (no CLI involved at all) still resolves.
        const cli::ResolvedRun base = cli::resolve_run(cli::parse_cli(bounded_args()));
        simulation::EffectiveOverrides overrides;
        overrides.infector_policy = simulation::InfectorPolicy::risk_averse;
        overrides.infector_policy_explicit = true;
        const simulation::EffectiveScenario scenario =
            simulation::resolve_scenario(*base.category, base.config, overrides);
        check(scenario.infector_policy == simulation::InfectorPolicy::risk_averse,
              "a valid direct resolve_scenario call resolves the risk-averse policy");
    }

    // --- Output mode: scenario resolution (additive) ------------------------
    {
        const cli::ResolvedRun run = cli::resolve_run(cli::parse_cli(bounded_args()));
        check(run.scenario.mode == simulation::OutputMode::ecai,
              "a bare run resolves to ecai mode");
        check(run.scenario.input_ecai_l_s_person == 0.0,
              "ecai mode carries no probability-mode input ECAi");
        check(!run.scenario.input_ecai_overridden,
              "ecai mode never reports an input-ECAi override");
        check(run.overrides.mode == simulation::OutputMode::ecai,
              "the effective overrides record ecai mode by default");
        check(!run.overrides.input_ecai_l_s_person.has_value(),
              "the effective overrides carry no input ECAi by default");
    }
    {
        const cli::ResolvedRun run =
            cli::resolve_run(cli::parse_cli(with_flag("--mode", "probability")));
        check(run.scenario.mode == simulation::OutputMode::probability,
              "--mode probability resolves to probability mode");
        check(run.scenario.input_ecai_l_s_person == run.category->published_ecai_Lsp,
              "with no --input-ecai override, probability mode defaults to the category's "
              "published Table 4 ECAi");
        check(!run.scenario.input_ecai_overridden,
              "no explicit --input-ecai leaves input_ecai_overridden false");
        check(!run.custom_scenario,
              "--mode probability alone (no --input-ecai) is not a custom scenario");
    }
    {
        Args args = with_flag("--mode", "probability");
        args.push_back("--input-ecai");
        args.push_back("42.5");
        const cli::ResolvedRun run = cli::resolve_run(cli::parse_cli(args));
        check(run.scenario.mode == simulation::OutputMode::probability,
              "--mode probability with an explicit --input-ecai resolves to probability mode");
        check(run.scenario.input_ecai_l_s_person == 42.5,
              "an explicit --input-ecai reaches the resolved scenario");
        check(run.scenario.input_ecai_overridden,
              "an explicit --input-ecai sets input_ecai_overridden");
        check(run.custom_scenario, "an explicit --input-ecai marks the resolved run custom");
        check(run.overrides.mode == simulation::OutputMode::probability &&
                  run.overrides.input_ecai_l_s_person.has_value() &&
                  *run.overrides.input_ecai_l_s_person == 42.5,
              "the effective overrides carry the explicit mode and input ECAi");
    }
    {
        // Tamper check, mirroring the existing scalar-provenance tamper block
        // above: an out-of-band mode change without mode_explicit must be
        // rejected, not silently accepted.
        cli::CliOptions options = cli::parse_cli(bounded_args());
        options.mode = simulation::OutputMode::probability;
        bool rejected = false;
        try {
            cli::resolve_run(options);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check(rejected, "resolve_run rejects an unprovenanced change to the output mode");
    }

    // --- Structural rejections ----------------------------------------------
    {
        expect_reject(Args{}, "an invocation with no arguments at all");
        expect_reject(Args{"--simulations", "8"}, "an invocation with no --category");
        expect_reject(Args{"--category", ""}, "an empty category string");
        expect_reject(Args{"--category", "Exam room", "--this-flag-does-not-exist"},
                      "an unknown long flag");
        expect_reject(Args{"--category", "Exam room", "-x"}, "an unknown short flag");
        expect_reject(Args{"--category", "Exam room", "extra"}, "a positional argument");
        expect_reject(Args{"Exam room"}, "a leading positional argument");
        expect_reject(Args{"--category"}, "a flag whose value is missing at the end");
        expect_reject(Args{"--category", "--simulations", "8"},
                      "a value that is itself a flag");
        expect_reject(Args{"--category", "Exam room", "--simulations"},
                      "a trailing flag with no value");
        expect_reject(Args{"--category", "Exam room", "--category", "Patient room"},
                      "a duplicate --category");
        expect_reject(Args{"--category", "Exam room", "--simulations", "8", "--simulations", "9"},
                      "a duplicate --simulations");
        expect_reject(Args{"--category", "Exam room", "--write-replicates", "--write-replicates"},
                      "a repeated boolean flag");
        expect_reject(Args{"--help", "--category", "Exam room"},
                      "--help combined with another flag");
        expect_reject(Args{"--category", "Exam room", "--help"},
                      "--help supplied after another flag");
        expect_reject(Args{"--help", "--help"}, "a duplicated --help");
        expect_reject(Args{"--category", "Exam room", "--output-directory", ""},
                      "an empty --output-directory path");
        expect_reject(Args{"--category", "Exam room", "--summary-file", ""},
                      "an empty --summary-file path");
        expect_reject(Args{"--category", "Exam room", "--replicates-file", ""},
                      "an empty --replicates-file path");
    }

    // --- Strict numeric parsing ---------------------------------------------
    {
        const std::vector<std::string> malformed = {"",      " 0.5",  "0.5 ", "0.5x", "x0.5",
                                                    "abc",   "0.5.5", "1,5",  "--",   "0x1"};
        for (const std::string& text : malformed) {
            expect_reject(with_flag("--percentile", text),
                          "a malformed --percentile value \"" + text + "\"");
        }
        const std::vector<std::string> nonfinite = {"nan",  "NaN",       "NAN", "inf",
                                                    "INF",  "Infinity",  "-inf"};
        for (const std::string& text : nonfinite) {
            expect_reject(with_flag("--volume-m3", text),
                          "a non-finite --volume-m3 value \"" + text + "\"");
        }
        expect_reject(with_flag("--volume-m3", "1e400"), "a --volume-m3 value that overflows double");
        expect_reject(with_flag("--volume-m3", "1e-400"),
                      "a --volume-m3 value that underflows double");
        expect_reject(with_flag("--percentile", "+0.5"),
                      "a leading plus sign, which from_chars does not accept");
        expect_reject(only_flag("--simulations", "+8"),
                      "a leading plus sign on an integer field");
        expect_reject(only_flag("--simulations", "8.5"), "a non-integer --simulations value");
        expect_reject(only_flag("--simulations", "8 "), "a trailing space on --simulations");
        expect_reject(only_flag("--replicates", "2.0"), "a non-integer --replicates value");
        expect_reject(with_flag("--occupants", "3.0"), "a non-integer --occupants value");
        expect_reject(with_flag("--occupants", "1e2"), "an exponent form for --occupants");
        expect_reject(with_flag("--base-seed", "-1"), "a negative --base-seed");
        expect_reject(with_flag("--base-seed", "18446744073709551616"),
                      "a --base-seed above the unsigned 64 bit range");
        expect_accept(with_flag("--base-seed", "18446744073709551615"),
                      "the largest unsigned 64 bit --base-seed");
        expect_accept(with_flag("--base-seed", "0"), "a zero --base-seed");
        expect_reject(with_flag("--occupants", "9223372036854775808"),
                      "an --occupants value above the signed 64 bit range");
    }

    // --- Domain bounds -------------------------------------------------------
    {
        expect_reject(only_flag("--simulations", "0"), "N = 0");
        expect_accept(only_flag("--simulations", "1"), "N = 1");
        const std::size_t max_n = sampling_plan::ResourcePolicy{}.max_simulations;
        expect_accept(only_flag("--simulations", std::to_string(max_n)),
                      "N at the frozen resource-policy maximum");
        expect_reject(only_flag("--simulations", std::to_string(max_n + 1)),
                      "N above the frozen resource-policy maximum");

        expect_reject(only_flag("--replicates", "0"), "R = 0");
        expect_accept(only_flag("--replicates", "1"), "R = 1");
        expect_accept(only_flag("--replicates", std::to_string(
                                    statistics::kMaxSupportedReplicateCount)),
                      "R at the frozen Student-t maximum");
        expect_reject(only_flag("--replicates", std::to_string(
                                    statistics::kMaxSupportedReplicateCount + 1)),
                      "R above the frozen Student-t maximum");

        expect_reject(with_flag("--occupants", "0"), "I0 = 0");
        expect_reject(with_flag("--occupants", "-3"), "a negative I0");
        expect_accept(with_flag("--occupants", "1"), "I0 = 1");

        expect_reject(with_flag("--community-rate", "-0.001"), "a community rate below 0");
        expect_reject(with_flag("--community-rate", "1.001"), "a community rate above 1");
        expect_accept(with_flag("--community-rate", "0"), "a community rate of exactly 0");
        expect_accept(with_flag("--community-rate", "1"), "a community rate of exactly 1");

        expect_reject(with_flag("--mask-efficiency", "-0.001"), "a mask efficiency below 0");
        expect_reject(with_flag("--mask-efficiency", "1.001"), "a mask efficiency above 1");
        expect_accept(with_flag("--mask-efficiency", "0"), "a mask efficiency of exactly 0");
        expect_accept(with_flag("--mask-efficiency", "1"), "a mask efficiency of exactly 1");

        expect_reject(with_flag("--target-probability", "0"), "a target probability of exactly 0");
        expect_reject(with_flag("--target-probability", "1"), "a target probability of exactly 1");
        expect_reject(with_flag("--target-probability", "-0.1"),
                      "a negative target probability");
        expect_accept(with_flag("--target-probability", "0.5"),
                      "a target probability inside (0, 1)");

        expect_reject(with_flag("--confidence", "0"), "a confidence level of exactly 0");
        expect_reject(with_flag("--confidence", "1"), "a confidence level of exactly 1");
        expect_accept(with_flag("--confidence", "0.99"), "a confidence level inside (0, 1)");

        expect_reject(with_flag("--percentile", "-0.001"), "a percentile below 0");
        expect_reject(with_flag("--percentile", "1.001"), "a percentile above 1");
        expect_accept(with_flag("--percentile", "0"), "a percentile of exactly 0");
        expect_accept(with_flag("--percentile", "1"), "a percentile of exactly 1");

        expect_reject(with_flag("--volume-m3", "0"), "a volume of exactly 0");
        expect_reject(with_flag("--volume-m3", "-1"), "a negative volume");
        expect_accept(with_flag("--volume-m3", "0.5"), "a small positive volume");

        expect_reject(with_flag("--duration-hours", "0"), "a duration of exactly 0");
        expect_reject(with_flag("--duration-hours", "-1"), "a negative duration");
        expect_accept(with_flag("--duration-hours", "0.25"), "a fractional positive duration");
    }

    // --- Category resolution (docs/MODEL_SPECIFICATION.md D7) ---------------
    {
        const cli::ResolvedRun run = cli::resolve_run(cli::parse_cli(bounded_args()));
        check(run.category != nullptr, "resolve_run retains the occupancy category record");
        check(run.category->canonical_name == "Exam room",
              "the retained record carries the canonical publication name");
        check(run.config.category_name == "Exam room",
              "the effective configuration names the canonical category");
        check(run.scenario.category_name == "Exam room",
              "the resolved scenario names the canonical category");
        // Jones et al. Table 2 (I0, VOL) and Table 1 / section 3.3 (rate, mask).
        check(run.category->i0 == 3, "Exam room source default I0 is 3");
        check(run.category->volume_m3 == 41.0, "Exam room source default VOL is 41 m^3");
        check(run.category->community_infection_rate == 0.03,
              "Exam room source default community rate is the health care 0.03");
        check(run.category->mask_efficiency == 0.30,
              "Exam room source default mask efficiency is the health care 0.30");
        check(run.config.i0 == run.category->i0 &&
                  run.config.volume_m3 == run.category->volume_m3 &&
                  run.config.community_infection_rate == run.category->community_infection_rate &&
                  run.config.mask_efficiency == run.category->mask_efficiency,
              "with no override, every effective category input equals its source default");
        check(!run.overrides.occupants && !run.overrides.volume_m3 &&
                  !run.overrides.mask_efficiency && !run.overrides.community_infection_rate,
              "no EffectiveOverrides entry is set without an explicit flag");
        check(run.config.simulations_n == 8 && run.config.replicates_r == 2,
              "N and R reach the effective configuration");
        check(run.config.percentile == 0.96 && run.config.confidence_level == 0.95 &&
                  run.config.target_probability == 0.001 && run.config.duration_hours == 1.0,
              "the frozen scientific defaults reach the effective configuration");
        check(run.base_seed == 485, "the default base seed reaches the resolved run");
        check(!run.custom_scenario, "an unmodified category run is not custom");
        check(!run.scenario.custom_scenario,
              "the scenario reports no category-input override either");
        check(run.target_is_hourly, "the 1 h default duration is an hourly target");
    }
    {
        const std::vector<std::string> spellings = {"Exam room", "exam room", "EXAM ROOM",
                                                    "exam-room", "Exam,room", "  exam   room  "};
        for (const std::string& spelling : spellings) {
            const cli::ResolvedRun run =
                cli::resolve_run(cli::parse_cli(Args{"--category", spelling}));
            check(run.category != nullptr && run.category->canonical_name == "Exam room",
                  "category spelling [" + spelling + "] resolves to Exam room");
        }
        const std::vector<std::pair<std::string, std::string>> aliases = {
            {"Patient", "Patient room"},
            {"Dwelling", "Dwelling unit"},
            {"Sorting", "Sorting, packing, light assembly"},
            {"place-of-religious-worship", "Place of religious worship"}};
        for (const auto& entry : aliases) {
            const cli::ResolvedRun run =
                cli::resolve_run(cli::parse_cli(Args{"--category", entry.first}));
            check(run.category != nullptr && run.category->canonical_name == entry.second,
                  "[" + entry.first + "] resolves to the canonical [" + entry.second + "]");
        }
    }
    {
        const std::vector<std::string> unknown = {"Exam", "Exam roomm", "exam_room", "Pat",
                                                  "Patient roo", "Not a category"};
        for (const std::string& name : unknown) {
            bool rejected = false;
            try {
                cli::resolve_run(cli::parse_cli(Args{"--category", name}));
            } catch (const std::invalid_argument&) {
                rejected = true;
            }
            check(rejected, "no fuzzy or prefix match for the category [" + name + "]");
        }
    }

    // --- Override precedence and preserved source defaults -------------------
    {
        const cli::CliOptions options = cli::parse_cli(Args{
            "--category", "Exam room", "--simulations", "8", "--replicates", "2", "--occupants",
            "5", "--volume-m3", "60", "--mask-efficiency", "0.1", "--community-rate", "0.02"});
        const cli::ResolvedRun run = cli::resolve_run(options);
        check(run.config.i0 == 5 && run.config.volume_m3 == 60.0 &&
                  run.config.mask_efficiency == 0.1 &&
                  run.config.community_infection_rate == 0.02,
              "every explicit category override wins over the category default");
        check(run.scenario.occupants == 5 && run.scenario.volume_m3 == 60.0 &&
                  run.scenario.mask_efficiency == 0.1 &&
                  run.scenario.community_infection_rate == 0.02,
              "the resolved scenario carries the overridden values");
        check(run.overrides.occupants && *run.overrides.occupants == 5 &&
                  run.overrides.volume_m3 && *run.overrides.volume_m3 == 60.0 &&
                  run.overrides.mask_efficiency && *run.overrides.mask_efficiency == 0.1 &&
                  run.overrides.community_infection_rate &&
                  *run.overrides.community_infection_rate == 0.02,
              "each override is carried in EffectiveOverrides");
        check(run.category->i0 == 3 && run.category->volume_m3 == 41.0 &&
                  run.category->mask_efficiency == 0.30 &&
                  run.category->community_infection_rate == 0.03,
              "the source category defaults are preserved, not overwritten");
        check(run.scenario.occupants_overridden && run.scenario.volume_overridden &&
                  run.scenario.mask_efficiency_overridden &&
                  run.scenario.community_infection_rate_overridden,
              "per-field override provenance reaches the scenario");
        check(run.custom_scenario && run.scenario.custom_scenario,
              "category-input overrides make the run a custom scenario");
    }
    {
        const cli::ResolvedRun run =
            cli::resolve_run(cli::parse_cli(with_flag("--target-probability", "0.002")));
        check(run.config.target_probability == 0.002,
              "an explicit target probability reaches the effective configuration");
        check(run.custom_scenario,
              "a non-category scientific override still makes the run custom");
        check(!run.scenario.custom_scenario,
              "the scenario custom flag covers only the four category inputs");
    }
    {
        const cli::ResolvedRun run =
            cli::resolve_run(cli::parse_cli(with_flag("--duration-hours", "4")));
        check(run.config.duration_hours == 4.0 && run.scenario.duration_hours == 4.0,
              "an explicit duration reaches both the configuration and the scenario");
        check(!run.target_is_hourly, "a 4 h duration is never labelled an hourly target");
        check(run.custom_scenario, "a duration override makes the run custom");
    }
    {
        const cli::ResolvedRun run =
            cli::resolve_run(cli::parse_cli(with_flag("--duration-hours", "1")));
        check(run.target_is_hourly, "an explicit duration of exactly 1 h stays an hourly target");
        check(run.custom_scenario, "and is still recorded as custom");
    }

    // --- Provenance agreement is required, not assumed -----------------------
    {
        cli::CliOptions options = cli::parse_cli(with_flag("--volume-m3", "60"));
        options.custom_scenario = false;  // tampered: an override is present
        bool rejected = false;
        try {
            cli::resolve_run(options);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check(rejected,
              "resolve_run rejects a custom_scenario flag that contradicts the overrides");
    }
    {
        cli::CliOptions options = cli::parse_cli(with_flag("--duration-hours", "3"));
        options.target_is_hourly = true;  // tampered: the duration is not 1 h
        bool rejected = false;
        try {
            cli::resolve_run(options);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check(rejected, "resolve_run rejects an hourly label that contradicts the duration");
    }
    {
        struct ScalarTamper {
            const char* name;
            void (*mutate)(cli::CliOptions&);
        };
        const ScalarTamper tampers[] = {
            {"simulations", [](cli::CliOptions& o) { o.simulations_n = 9; }},
            {"replicates", [](cli::CliOptions& o) { o.replicates_r = 3; }},
            {"confidence", [](cli::CliOptions& o) { o.confidence_level = 0.9; }},
            {"percentile", [](cli::CliOptions& o) { o.percentile = 0.9; }},
            {"target probability", [](cli::CliOptions& o) { o.target_probability = 0.002; }},
            {"duration", [](cli::CliOptions& o) { o.duration_hours = 2.0; }},
            {"base seed", [](cli::CliOptions& o) { o.base_seed = 486; }},
        };
        for (const ScalarTamper& tamper : tampers) {
            cli::CliOptions options = cli::parse_cli(Args{"--category", "Exam room"});
            tamper.mutate(options);
            bool rejected = false;
            try {
                cli::resolve_run(options);
            } catch (const std::invalid_argument&) {
                rejected = true;
            }
            check(rejected, std::string("resolve_run rejects an unprovenanced change to ") +
                                tamper.name);
        }
    }
    {
        cli::CliOptions options = cli::parse_cli(bounded_args());
        options.category_supplied = false;
        bool rejected = false;
        try {
            cli::resolve_run(options);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check(rejected, "resolve_run refuses to resolve without a category");
    }
    {
        expect_reject(Args{"--category", "Exam room", "--percentile", "0.5\nx"},
                      "a numeric argument containing a line break");
        bool single_line = false;
        try {
            cli::resolve_run(cli::parse_cli(Args{"--category", "Bad\nName"}));
        } catch (const std::invalid_argument& error) {
            single_line = std::string(error.what()).find('\n') == std::string::npos &&
                          std::string(error.what()).find('\r') == std::string::npos;
        }
        check(single_line, "an unknown category containing a line break has a single-line diagnostic");
    }
    {
        bool rejected = false;
        try {
            cli::resolve_run(cli::parse_cli(Args{"--help"}));
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check(rejected, "resolve_run refuses to resolve a --help invocation");
    }

    // --- Executable command matrix ------------------------------------------
    // Every run below uses an explicitly small N and R: the CLI tests never
    // execute the 10000 x 500 defaults.
    {
        const CommandResult result = run_command("--help");
        check(result.exit_code == 0, "`--help` exits 0");
        check(result.output == cli::help_text(),
              "`--help` prints exactly the current complete help text");
        check(result.output.find("console_summary") == std::string::npos,
              "`--help` runs no simulation and prints no summary");
    }
    {
        const CommandResult result = run_command("");
        check(result.exit_code != 0, "an invocation with no arguments exits nonzero");
        check(result.output.find("--category") != std::string::npos,
              "the no-argument error names the required --category flag");
        check(result.output.find("Try ") != std::string::npos,
              "the no-argument error points at --help");
        check(count_lines(result.output) == 2,
              "the no-argument error is one diagnostic line plus one Try line");
        check(result.output.find("console_summary") == std::string::npos,
              "a failed invocation prints no scientific output");
    }
    {
        const std::vector<std::pair<std::string, std::string>> failures = {
            {"--this-flag-does-not-exist", "an unknown flag"},
            {"--category", "a flag with a missing value"},
            {"--category \"Not a category\"", "an unknown category"},
            {"--category \"Exam room\" --simulations 0", "N = 0"},
            {"--category \"Exam room\" --percentile 1.5", "a percentile above 1"},
            {"--category \"Exam room\" --duration-hours nan", "a NaN duration"},
            {"--category \"Exam room\" extra", "a positional argument"},
            {"--help --category \"Exam room\"", "--help combined with another flag"}};
        for (const auto& entry : failures) {
            const CommandResult result = run_command(entry.first);
            check(result.exit_code != 0, "the executable exits nonzero for " + entry.second);
            check(count_lines(result.output) == 2,
                  "the executable prints exactly two error lines for " + entry.second);
            check(result.output.find("console_summary") == std::string::npos,
                  "no partial scientific output for " + entry.second);
        }
    }
    {
        const std::string command =
            "--category \"Exam room\" --simulations 4 --replicates 2";
        const CommandResult result = run_command(command);
        check(result.exit_code == 0, "a valid bounded invocation exits 0");
        check(count_lines(result.output) == 1, "a valid run prints exactly one summary line");
        check(result.output.find("console_summary") != std::string::npos,
              "the summary line is marked as a console summary");
        check(result.output.find("not a CSV") != std::string::npos,
              "the summary line states that it is not a CSV file");
        check(result.output.find("Exam room") != std::string::npos,
              "the summary line names the canonical category");
        check(result.output.find("raw_mean_l_s_person") != std::string::npos,
              "the summary line reports the raw mean in L/s per person");
        check(result.output.find("ceiling_to_5_l_s_person") != std::string::npos,
              "the summary line reports the ceiling-to-5 value in L/s per person");

        // The printed values must round-trip exactly to the values an
        // in-process run of the same command produces.
        const cli::ResolvedRun run = cli::resolve_run(cli::parse_cli(
            Args{"--category", "Exam room", "--simulations", "4", "--replicates", "2"}));
        const simulation::ReplicateBatch batch = simulation::run_replicates(
            run.base_seed, run.config, *run.category, run.overrides);
        check(field_value(result.output, "raw_mean_l_s_person=") ==
                  batch.summary.mean_l_s_person,
              "the printed raw mean round-trips bit-exactly");
        check(field_value(result.output, "ceiling_to_5_l_s_person=") ==
                  batch.summary.adjusted_ceiling_to_5_l_s_person,
              "the printed ceiling-to-5 value round-trips bit-exactly");
        check(count_occurrences(result.output, "raw_mean_l_s_person") == 1 &&
                  count_occurrences(result.output, "ceiling_to_5_l_s_person") == 1,
              "no statistic is reported twice");

        const CommandResult repeat = run_command(command);
        check(repeat.output == result.output, "the same command produces identical output");
    }
    {
        const CommandResult result =
            run_command("--category Patient --simulations 4 --replicates 1");
        check(result.exit_code == 0, "a legacy alias runs successfully");
        check(result.output.find("Patient room") != std::string::npos,
              "the summary emits the canonical name, never the alias");
    }
    {
        std::error_code cleanup_error;
        std::filesystem::remove_all("stage12 out dir", cleanup_error);
        const CommandResult result = run_command(
            "--category exam-room --simulations 4 --replicates 2 --percentile 0.9 "
            "--confidence 0.9 --target-probability 0.002 --duration-hours 2 "
            "--base-seed 17 --community-rate 0.02 --occupants 2 --volume-m3 50 "
            "--mask-efficiency 0.2 --output-directory \"stage12 out dir\" "
            "--summary-file \"stage12 out dir/summary one.csv\" "
            "--replicates-file \"stage12 out dir/reps one.csv\" --write-replicates");
        check(result.exit_code == 0, "every flag together runs successfully");
        check(result.output.find("Exam room") != std::string::npos,
              "a delimiter-variant category still reports the canonical name");
        check(path_exists("stage12 out dir"),
              "Stage 13 creates the explicitly requested output directory");
        check(path_exists("stage12 out dir/summary one.csv"),
              "Stage 13 writes the explicitly requested summary file");
        check(path_exists("stage12 out dir/reps one.csv"),
              "Stage 13 writes the explicitly requested per-replicate file");
        std::filesystem::remove_all("stage12 out dir", cleanup_error);
    }

    if (g_failures > 0) {
        std::cerr << g_failures << " of " << g_checks << " check(s) failed.\n";
        return 1;
    }
    std::cout << "All " << g_checks << " checks passed.\n";
    return 0;
}
