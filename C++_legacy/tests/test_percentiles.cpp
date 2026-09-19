// test_percentiles.cpp — exercises the `percentiles` analysis executable's
// binary-report path (analysis/percentiles.cpp) for:
//   1. valid interpolation on a normal series
//   2. out-of-range percentiles (negative / >100) clamped instead of
//      indexing out of bounds
//   3. an empty series handled without a crash
//   4. the threshold search reporting no-data (nonzero exit) rather than
//      "all categories pass" when every category is empty
//
// Writes small synthetic .bin fixtures in the ASHR binary format that
// load_binary() in analysis/percentiles.cpp expects, then runs the built
// `percentiles` executable as a subprocess and checks its exit code and
// output, mirroring the subprocess pattern used by test_cli_n_validation.cpp.

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <filesystem>
#include <fstream>

#ifdef _WIN32
#define WR_POPEN _popen
#define WR_PCLOSE _pclose
#else
#include <sys/wait.h>
#define WR_POPEN popen
#define WR_PCLOSE pclose
#endif

namespace fs = std::filesystem;

static int failures = 0;
static int checks = 0;

static void write_category(std::ofstream& out, const std::string& name,
                            const std::vector<float>& values) {
    uint32_t nl = static_cast<uint32_t>(name.size());
    out.write(reinterpret_cast<const char*>(&nl), sizeof(nl));
    out.write(name.data(), nl);
    uint32_t cnt = static_cast<uint32_t>(values.size());
    out.write(reinterpret_cast<const char*>(&cnt), sizeof(cnt));
    if (cnt) out.write(reinterpret_cast<const char*>(values.data()), cnt * sizeof(float));
}

static void write_header(std::ofstream& out, uint32_t nc) {
    uint32_t magic = 0x41534852;
    out.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
    out.write(reinterpret_cast<const char*>(&nc), sizeof(nc));
}

// Fixture: "Known" has 5 values (valid interpolation target) and "Empty"
// has zero values (must not crash the percentile calculation).
static void write_fixture(const fs::path& path) {
    std::ofstream out(path, std::ios::binary);
    write_header(out, 2);
    write_category(out, "Known", {0.5f, 0.1f, 0.3f, 0.4f, 0.2f});
    write_category(out, "Empty", {});
}

// Fixture where every category is empty — the threshold search has nothing
// to evaluate and must not claim all categories are below the target.
static void write_empty_fixture(const fs::path& path) {
    std::ofstream out(path, std::ios::binary);
    write_header(out, 2);
    write_category(out, "EmptyA", {});
    write_category(out, "EmptyB", {});
}

struct RunResult {
    bool ran = false;
    int exit_code = -1;
    std::string output;
};

// Runs `../percentiles <args>` with cwd = scratch dir, capturing stdout+stderr.
static RunResult run(const std::string& scratch, const std::string& args) {
    RunResult result;
#ifdef _WIN32
    std::string cmd = "cd /d \"" + scratch + "\" && ..\\percentiles.exe " + args + " 2>&1";
#else
    std::string cmd = "cd '" + scratch + "' && ../percentiles " + args + " 2>&1";
#endif
    FILE* pipe = WR_POPEN(cmd.c_str(), "r");
    if (!pipe) return result;
    std::array<char, 512> buf{};
    while (fgets(buf.data(), buf.size(), pipe)) result.output += buf.data();
    int status = WR_PCLOSE(pipe);
    result.ran = true;
#ifdef _WIN32
    result.exit_code = status;
#else
    if (WIFEXITED(status)) result.exit_code = WEXITSTATUS(status);
#endif
    return result;
}

static void check(bool cond, const std::string& label, const std::string& detail) {
    checks++;
    if (!cond) {
        printf("FAIL: %s — %s\n", label.c_str(), detail.c_str());
        failures++;
    }
}

int main() {
    fs::path scratch = fs::path("percentiles_tmp");
    std::error_code ec;
    fs::create_directories(scratch, ec);
    write_fixture(scratch / "fixture.bin");
    write_empty_fixture(scratch / "all_empty.bin");

    // Valid interpolation: sorted Known = [0.1,0.2,0.3,0.4,0.5]; p=50 -> rank
    // 2.0 -> sorted[2]=0.3 -> reported as 30.000.
    // Out-of-range percentiles clamp: p=-10 -> 0th -> 10.000; p=150 -> 100th
    // -> 50.000. Empty series must not crash the process.
    RunResult r = run(scratch.string(), "fixture.bin --report table --percentiles 50 -10 150");

    check(r.ran, "percentiles subprocess", "could not launch");
    if (r.ran) {
        check(r.exit_code == 0, "exit code",
              "expected 0, got " + std::to_string(r.exit_code) + " (output: " + r.output + ")");
        check(r.output.find("30.000") != std::string::npos, "valid interpolation (p=50)",
              "expected 30.000 in output: " + r.output);
        check(r.output.find("10.000") != std::string::npos, "clamped low (p=-10 -> 0th)",
              "expected 10.000 in output: " + r.output);
        check(r.output.find("50.000") != std::string::npos, "clamped high (p=150 -> 100th)",
              "expected 50.000 in output: " + r.output);
        check(r.output.find("Empty") != std::string::npos, "empty series row present",
              "expected Empty category row in output: " + r.output);
    }

    // All categories empty: worst_at() has no finite value to report, so the
    // threshold search must fail with a no-data message rather than declaring
    // that every category stays below the target.
    RunResult e = run(scratch.string(), "all_empty.bin --report threshold");

    check(e.ran, "all-empty subprocess", "could not launch");
    if (e.ran) {
        check(e.exit_code != 0, "all-empty exit code",
              "expected nonzero for no-data threshold search, got " +
                  std::to_string(e.exit_code) + " (output: " + e.output + ")");
        check(e.output.find("no data") != std::string::npos, "all-empty no-data message",
              "expected a no-data message in output: " + e.output);
        check(e.output.find("stay below") == std::string::npos, "all-empty must not report a pass",
              "unexpected pass claim in output: " + e.output);
    }

    fs::remove_all(scratch, ec);

    if (failures == 0) {
        printf("PASS: all percentiles tests (%d checks)\n", checks);
        return 0;
    } else {
        printf("FAIL: %d/%d percentiles check(s) failed\n", failures, checks);
        return 1;
    }
}
