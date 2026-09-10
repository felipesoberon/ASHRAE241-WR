#ifndef ECAI_CPP2_OUTPUT_H
#define ECAI_CPP2_OUTPUT_H

#include <cstddef>
#include <filesystem>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "cli.h"
#include "simulation.h"

namespace ecai_cpp2::output {

inline constexpr const char* kModelVersion = "ASHRAE241-Jones-2025-Cpp2-v1";
inline constexpr const char* kOutputSchemaVersion = "stage13-v1";

struct OutputDocument {
    std::vector<std::string> summary_columns;
    std::vector<std::vector<std::string>> summary_rows;
    std::vector<std::string> replicate_columns;
    std::vector<std::vector<std::string>> replicate_rows;
};

struct OutputPaths {
    bool requested = false;
    bool summary_requested = false;
    bool replicates_requested = false;
    std::string summary_path;
    std::string replicates_path;
};

struct WriteResult {
    OutputPaths paths;
    std::size_t summary_rows = 0;
    std::size_t replicate_rows = 0;
    // False only when destinations committed but best-effort temporary
    // payload cleanup could not be completed.
    bool cleanup_complete = true;
};

// Pure CSV serializers using RFC4180-style field quoting. Record endings are
// deliberately LF (not RFC4180's CRLF requirement).
std::string serialize_summary_csv(const OutputDocument& document);
std::string serialize_replicates_csv(const OutputDocument& document);

// Builds and validates the typed Stage 13 document without touching the file
// system. protected_tree is explicit so tests can exercise the legacy-tree
// refusal without relying on a build definition.
OutputDocument build_document(const cli::ResolvedRun& run,
                              const simulation::ReplicateBatch& batch,
                              const cli::CliOptions& options,
                              const std::string& protected_tree = {});

// Resolves and validates all output controls and destination paths. It does
// not create directories or files.
OutputPaths resolve_paths(const cli::CliOptions& options,
                          const std::string& protected_tree = {});

// Requires the exact Stage 13 89/19-column schemas and row-count metadata,
// then writes requested files atomically while refusing overwrite. Parent
// directories are created only after all validation passes.
WriteResult write_outputs(const cli::CliOptions& options,
                          const OutputDocument& document,
                          const std::string& protected_tree = {});

namespace detail {
// Deterministic failure/race seam used only by writer tests. Production leaves
// it empty. Stages are "before_payload_write", "before_payload_close", and
// "before_commit".
using WriteTestHook = std::function<void(const std::string&, std::size_t,
                                         const std::filesystem::path&)>;
void set_write_test_hook(WriteTestHook hook);
}  // namespace detail

}  // namespace ecai_cpp2::output

#endif
