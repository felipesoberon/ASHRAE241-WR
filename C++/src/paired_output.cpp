#include "paired_output.h"

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef ECAI_SOURCE_REVISION
#define ECAI_SOURCE_REVISION "unknown"
#endif
#ifndef ECAI_SOURCE_TRACKED_DIRTY
#define ECAI_SOURCE_TRACKED_DIRTY "unknown"
#endif
#ifndef ECAI_SOURCE_CONTENT_HASH
#define ECAI_SOURCE_CONTENT_HASH "unknown"
#endif
#ifndef ECAI_COMPILER_ID
#define ECAI_COMPILER_ID "unknown"
#endif

namespace ecai_cpp2::paired_output {
namespace {
namespace fs = std::filesystem;

void write_u64(std::ostream& output, std::uint64_t value) {
    char bytes[8];
    for (unsigned i = 0; i < 8; ++i) bytes[i] = static_cast<char>((value >> (i * 8)) & 0xff);
    output.write(bytes, 8);
    if (!output) throw std::runtime_error("failed writing paired results");
}

void write_double(std::ostream& output, double value) {
    static_assert(sizeof(double) == sizeof(std::uint64_t), "binary64 required");
    static_assert(std::numeric_limits<double>::is_iec559, "IEEE-754 required");
    std::uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    write_u64(output, bits);
}

bool contains_path(const fs::path& root, const fs::path& child) {
    auto a = root.begin();
    auto b = child.begin();
    for (; a != root.end(); ++a, ++b) {
        if (b == child.end() || *a != *b) return false;
    }
    return true;
}
}  // namespace

void run_and_write(const cli::ResolvedRun& run, const cli::CliOptions& options,
                   const char* protected_legacy_tree) {
    const auto& targets = options.target_probabilities;
    if (targets.empty() || !options.raw_output_directory) {
        throw std::invalid_argument("paired output requires targets and a destination");
    }
    // Bound the additional value matrices to 1 GiB per replicate, in addition
    // to the existing sampler's independent memory budget.
    constexpr std::uint64_t max_paired_values = (1ULL << 30) / sizeof(double);
    if (run.config.simulations_n > max_paired_values / targets.size()) {
        throw std::invalid_argument("paired result vectors exceed the 1 GiB per-replicate limit");
    }
    const fs::path destination = fs::absolute(*options.raw_output_directory).lexically_normal();
    if (destination.filename().empty() || !fs::is_directory(destination.parent_path())) {
        throw std::invalid_argument("paired output requires an existing parent and a new directory");
    }
    const fs::path resolved = fs::weakly_canonical(destination);
    if (protected_legacy_tree && *protected_legacy_tree &&
        contains_path(fs::weakly_canonical(protected_legacy_tree), resolved)) {
        throw std::invalid_argument("paired output cannot be written inside the legacy tree");
    }
    if (!fs::create_directory(destination)) {
        throw std::runtime_error("paired output directory already exists or cannot be created");
    }
    // A failed run is deliberately left as a partial directory. It has no
    // manifest.txt and cannot be mistaken for a completed dataset.
    const fs::path data_partial = destination / "ecai.bin.partial";
    std::ofstream data(data_partial, std::ios::binary | std::ios::out);
    if (!data) throw std::runtime_error("cannot create paired data file");
    data.write("ECAIP01\n", 8);
    write_u64(data, static_cast<std::uint64_t>(targets.size()));
    write_u64(data, run.config.simulations_n);
    write_u64(data, run.config.replicates_r);
    for (double target : targets) write_double(data, target);

    simulation::ReplicateOptions replicate_options;
    replicate_options.paired_target_probabilities.assign(targets.begin() + 1, targets.end());
    for (std::uint64_t replicate = 0; replicate < run.config.replicates_r; ++replicate) {
        const auto result = simulation::run_replicate(run.base_seed, run.config, replicate,
                                                      *run.category, run.overrides,
                                                      replicate_options);
        for (double value : result.ecai_l_s_person) write_double(data, value);
        for (const auto& values : result.paired_ecai_l_s_person) {
            for (double value : values) write_double(data, value);
        }
    }
    data.flush();
    if (!data) throw std::runtime_error("failed flushing paired data");
    data.close();
    if (!data) throw std::runtime_error("failed closing paired data");
    fs::rename(data_partial, destination / "ecai.bin");

    std::ofstream metadata(destination / "manifest.txt.partial", std::ios::out);
    if (!metadata) throw std::runtime_error("cannot create paired manifest");
    metadata << std::setprecision(std::numeric_limits<double>::max_digits10)
             << "format=ECAIP01\ncategory=" << std::quoted(run.config.category_name)
             << "\nunits=L/s/person\ninfector_policy="
             << (run.scenario.infector_policy == simulation::InfectorPolicy::community ? "community" : "risk-averse")
             << "\nqer_distribution=" << run.qer_selection.canonical_name
             << "\nbase_seed=" << run.base_seed
             << "\nseed_policy=" << seed_manager::kSeedPolicyVersion
             << "\nsource_revision_at_configure=" << ECAI_SOURCE_REVISION
             << "\ntracked_source_dirty_at_configure=" << ECAI_SOURCE_TRACKED_DIRTY
             << "\nsource_content_sha256_at_configure=" << ECAI_SOURCE_CONTENT_HASH
             << "\ncompiler=" << ECAI_COMPILER_ID
             << "\nseed_target=first_listed_target"
             << "\noccupants=" << run.config.i0
             << "\nvolume_m3=" << run.config.volume_m3
             << "\nmask_efficiency=" << run.config.mask_efficiency
             << "\ncommunity_rate=" << run.config.community_infection_rate
             << "\nduration_hours=" << run.config.duration_hours
             << "\nseed_percentile=" << run.config.percentile
             << "\nconfidence=" << run.config.confidence_level
             << "\nsimulations=" << run.config.simulations_n
             << "\nreplicates=" << run.config.replicates_r << "\ntargets=";
    for (std::size_t i = 0; i < targets.size(); ++i) {
        if (i) metadata << ',';
        metadata << targets[i];
    }
    metadata << "\nlayout=replicate,target,simulation\n";
    metadata.flush();
    if (!metadata) throw std::runtime_error("failed writing paired manifest");
    metadata.close();
    if (!metadata) throw std::runtime_error("failed closing paired manifest");
    fs::rename(destination / "manifest.txt.partial", destination / "manifest.txt");
}
}  // namespace ecai_cpp2::paired_output
