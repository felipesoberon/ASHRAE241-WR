// Stage 1 test harness for ecai_cpp2. Exercises the built executable as a
// subprocess so behavior is verified end-to-end, not via disabled assert().

#include <array>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sys/wait.h>

#ifndef ECAI_CPP2_EXECUTABLE
#error "ECAI_CPP2_EXECUTABLE must be defined by the build system"
#endif

namespace {

struct CommandResult {
    int exit_code;
    std::string output;
};

CommandResult run_command(const std::string& args) {
    const std::string command =
        std::string("\"") + ECAI_CPP2_EXECUTABLE + "\" " + args + " 2>&1";

    FILE* pipe = popen(command.c_str(), "r");
    if (pipe == nullptr) {
        std::cerr << "Failed to launch: " << command << "\n";
        std::exit(1);
    }

    std::string output;
    std::array<char, 4096> buffer{};
    size_t bytes_read = 0;
    while ((bytes_read = std::fread(buffer.data(), 1, buffer.size(), pipe)) > 0) {
        output.append(buffer.data(), bytes_read);
    }

    const int status = pclose(pipe);
    int exit_code = -1;
    if (WIFEXITED(status)) {
        exit_code = WEXITSTATUS(status);
    }
    return CommandResult{exit_code, output};
}

int g_failures = 0;

void check(bool condition, const std::string& description) {
    if (!condition) {
        std::cerr << "FAIL: " << description << "\n";
        ++g_failures;
    } else {
        std::cout << "PASS: " << description << "\n";
    }
}

}  // namespace

int main() {
    {
        const CommandResult result = run_command("--help");
        check(result.exit_code == 0, "`--help` exits with status 0");
        check(result.output.find("ecai_cpp2") != std::string::npos,
              "`--help` output identifies the ecai_cpp2 executable");
        check(result.output.find("C++2") != std::string::npos,
              "`--help` output identifies the C++2 implementation");
    }
    {
        const CommandResult result = run_command("--this-flag-does-not-exist");
        check(result.exit_code != 0,
              "unsupported invocation exits with a nonzero status");
        check(result.output.find("ECAi") == std::string::npos,
              "unsupported invocation does not report any ECAi/model output");
    }
    {
        const CommandResult result = run_command("");
        check(result.exit_code != 0,
              "invocation with no arguments exits with a nonzero status");
    }

    if (g_failures > 0) {
        std::cerr << g_failures << " test(s) failed.\n";
        return 1;
    }
    std::cout << "All tests passed.\n";
    return 0;
}
