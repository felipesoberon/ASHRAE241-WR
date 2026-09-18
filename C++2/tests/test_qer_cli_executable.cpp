// Package 3b regression test: proves the --qer-distribution CLI selector
// actually propagates into the executed run, not merely into the CLI-layer
// ResolvedQerSelection descriptor (see tests/test_qer_cli_scenario.cpp,
// which is deliberately scoped to CLI/scenario resolution only and does not
// invoke the built executable or exercise sampling_plan/simulation
// integration).
//
// This file invokes the actual built ecai_cpp2 executable (mirroring the
// tests/test_cli.cpp "Executable command matrix" pattern) and reads the
// console_summary line's "qer_family=" marker
// (docs/QER_DISTRIBUTION_IMPLEMENTATION_PLAN.md Phase 5), which main.cpp
// prints from the EffectiveScenario the run actually executed under -- the
// same field simulation::resolve_scenario/run_replicate branch on to select
// the full-Jones versus direct-QER-family code path. This is deliberately
// not a CSV/output-schema assertion: no document field is added or read
// here, only the existing non-schema console line, exactly like the
// existing "infector_policy=" marker it sits next to.
//
// Before this test, --qer-distribution parsed and resolved correctly in
// isolation (test_qer_cli_scenario.cpp), but the resolved selector was never
// copied into simulation::EffectiveOverrides in cli::resolve_run, so every
// run -- regardless of the flag -- silently executed as full-Jones. This
// test is written to fail red against that defect and pass once
// cli::resolve_run actually propagates the selector through
// EffectiveOverrides -> EffectiveScenario -> the executed run.

#include <array>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

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
#ifdef WIFEXITED
    if (WIFEXITED(status)) {
        result.exit_code = WEXITSTATUS(status);
    }
#else
    result.exit_code = status;
#endif
    return result;
}

// Extracts the exact token following `label` up to the next whitespace, so a
// substring match cannot confuse "jones" with "jones-fitted" or a Mikszewski
// profile name that merely starts with "mikszewski".
std::string field_token(const std::string& text, const std::string& label) {
    const std::size_t at = text.find(label);
    if (at == std::string::npos) return "";
    const std::size_t start = at + label.size();
    std::size_t end = text.find(' ', start);
    if (end == std::string::npos) end = text.find('\n', start);
    if (end == std::string::npos) end = text.size();
    return text.substr(start, end - start);
}

std::string bounded_command(const std::string& extra_flags = "") {
    return "--category \"Exam room\" --simulations 4 --replicates 1 --base-seed 7" +
           (extra_flags.empty() ? "" : " " + extra_flags);
}

}  // namespace

int main() {
    // Omitted --qer-distribution: the executed run resolves to full Jones.
    {
        const CommandResult result = run_command(bounded_command());
        check(result.exit_code == 0, "omitted --qer-distribution: executable exits 0");
        check(field_token(result.output, "qer_family=") == "jones",
              "omitted --qer-distribution: executed run's qer_family marker reads exactly "
              "\"jones\"");
    }

    // --qer-distribution jones-fitted: the executed run reaches JonesFitted,
    // not the silently-defaulted full Jones path.
    {
        const CommandResult result =
            run_command(bounded_command("--qer-distribution jones-fitted"));
        check(result.exit_code == 0, "--qer-distribution jones-fitted: executable exits 0");
        check(field_token(result.output, "qer_family=") == "jones-fitted",
              "--qer-distribution jones-fitted: executed run's qer_family marker reads exactly "
              "\"jones-fitted\"");
    }

    // --qer-distribution mikszewski-influenza-resting: the executed run
    // reaches Mikszewski, not the silently-defaulted full Jones path.
    {
        const CommandResult result =
            run_command(bounded_command("--qer-distribution mikszewski-influenza-resting"));
        check(result.exit_code == 0,
              "--qer-distribution mikszewski-influenza-resting: executable exits 0");
        check(field_token(result.output, "qer_family=") == "mikszewski",
              "--qer-distribution mikszewski-influenza-resting: executed run's qer_family "
              "marker reads exactly \"mikszewski\"");
    }

    if (g_failures > 0) {
        std::cerr << g_failures << " of " << g_checks << " check(s) failed.\n";
        return 1;
    }
    std::cout << "All " << g_checks << " checks passed.\n";
    return 0;
}
