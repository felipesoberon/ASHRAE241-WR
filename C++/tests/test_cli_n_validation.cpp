// test_cli_n_validation.cpp — CLI entry points must reject N <= 0 with a
// nonzero exit status and an explanatory stderr message, and must accept a
// small positive N. Runs the built executables as subprocesses from a
// scratch subdirectory of the build dir so side-effect output files
// (CSV/bin) don't clobber other tests' fixtures.

#include <cstdio>
#include <cstdlib>
#include <string>
#include <array>
#include <filesystem>
#include <sys/wait.h>

namespace fs = std::filesystem;

static int failures = 0;
static int checks = 0;

struct RunResult {
    bool ran = false;
    int exit_code = -1;
    std::string stderr_output;
};

// Runs `../<exe> <args>` with cwd = scratch dir, capturing only stderr.
static RunResult run(const std::string& scratch, const std::string& exe,
                      const std::string& args) {
    RunResult result;
    std::string cmd = "cd '" + scratch + "' && ../" + exe + " " + args +
                       " 2>&1 1>/dev/null";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return result;
    std::array<char, 512> buf{};
    while (fgets(buf.data(), buf.size(), pipe)) result.stderr_output += buf.data();
    int status = pclose(pipe);
    result.ran = true;
    if (WIFEXITED(status)) result.exit_code = WEXITSTATUS(status);
    return result;
}

static void expect_rejected(const std::string& scratch, const std::string& exe,
                             const std::string& args, const std::string& label) {
    checks++;
    RunResult r = run(scratch, exe, args);
    if (!r.ran) {
        printf("FAIL: %s — could not launch subprocess\n", label.c_str());
        failures++;
        return;
    }
    if (r.exit_code == 0) {
        printf("FAIL: %s — expected nonzero exit for N<=0, got 0\n", label.c_str());
        failures++;
        return;
    }
    if (r.stderr_output.find("Error: N must be a positive integer") == std::string::npos) {
        printf("FAIL: %s — stderr missing expected message. Got: %s\n",
               label.c_str(), r.stderr_output.c_str());
        failures++;
    }
}

static void expect_accepted(const std::string& scratch, const std::string& exe,
                             const std::string& args, const std::string& label) {
    checks++;
    RunResult r = run(scratch, exe, args);
    if (!r.ran) {
        printf("FAIL: %s — could not launch subprocess\n", label.c_str());
        failures++;
        return;
    }
    if (r.exit_code != 0) {
        printf("FAIL: %s — expected exit 0 for small positive N, got %d (stderr: %s)\n",
               label.c_str(), r.exit_code, r.stderr_output.c_str());
        failures++;
    }
}

int main() {
    fs::path scratch = fs::path("cli_n_validation_tmp");
    std::error_code ec;
    fs::create_directories(scratch, ec);

    // ---- ecai: positional N ----
    expect_rejected(scratch.string(), "ecai", "0", "ecai N=0");
    expect_accepted(scratch.string(), "ecai", "5", "ecai N=5");

    // ---- probability_scan: positional N ----
    expect_rejected(scratch.string(), "probability_scan", "0", "probability_scan N=0");
    expect_accepted(scratch.string(), "probability_scan", "5", "probability_scan N=5");

    // ---- probability_ecai: positional N ----
    expect_rejected(scratch.string(), "probability_ecai", "0", "probability_ecai N=0");
    expect_accepted(scratch.string(), "probability_ecai", "5", "probability_ecai N=5");

    // ---- single_probability: --N flag (also accepts negative values) ----
    expect_rejected(scratch.string(), "single_probability", "--N 0", "single_probability --N 0");
    expect_rejected(scratch.string(), "single_probability", "--N -5", "single_probability --N -5");
    expect_accepted(scratch.string(), "single_probability", "--N 5", "single_probability --N 5");

    fs::remove_all(scratch, ec);

    if (failures == 0) {
        printf("PASS: all CLI N-validation tests (%d checks)\n", checks);
        return 0;
    } else {
        printf("FAIL: %d/%d CLI N-validation check(s) failed\n", failures, checks);
        return 1;
    }
}
