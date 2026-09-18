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
// Probability-mode document schema
// (docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md section 5). The ECAi schema
// above is unaffected: build_document dispatches to build_probability_document
// only when the resolved run is in probability mode.
inline constexpr const char* kProbabilitySchemaVersion = "stage13-prob-v1";
// Risk-averse infector-policy document schemas
// (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md Phase 6). Distinct from, and never
// substituted for, the two community schemas above: build_document dispatches
// to one of these only when the resolved run's infector policy is
// risk_averse, and the community/probability schemas and their serialized
// bytes are otherwise emitted completely unchanged.
inline constexpr const char* kRiskAverseSchemaVersion = "stage13-ra-v1";
inline constexpr const char* kRiskAverseProbabilitySchemaVersion = "stage13-prob-ra-v1";
// Direct-QER schemas are deliberately distinct from every full-Jones schema.
// The suffix columns are run/replicate metadata, not per-simulation fields.
inline constexpr const char* kDirectOutputSchemaVersion = "stage13-qer-direct-v1";
inline constexpr const char* kDirectProbabilitySchemaVersion = "stage13-prob-qer-direct-v1";
inline constexpr const char* kDirectRiskAverseSchemaVersion = "stage13-ra-qer-direct-v1";
inline constexpr const char* kDirectRiskAverseProbabilitySchemaVersion =
    "stage13-prob-ra-qer-direct-v1";

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
