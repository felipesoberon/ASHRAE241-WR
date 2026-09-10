// Stage 4 tests for seed_manager: deterministic per-replicate/per-stream
// seed derivation, per docs/MODEL_SPECIFICATION.md section 11.3 and
// docs/SEED_POLICY.md, and IMPLEMENTATION_PLAN.md Section 7, Stage 4.
//
// All known-answer values below (FNV-1a, SplitMix64, canonical serialization
// strings, and full-tuple derived seeds) were independently derived by
// hand-assembling the exact frozen algorithm from section 11.3 in a
// throwaway reference program (std::to_chars-based %.17g-equivalent
// formatting, FNV-1a-64, SplitMix64, and the normative fold) and running it
// separately from seed_manager.cpp, then hard-coded here — not copied from
// seed_manager.cpp after the fact. splitmix64(0) == 0xe220a8397b1dcdaf is
// also the well-known public first output of Vigna's reference splitmix64.c
// seeded at 0, independently cross-checking the mixing step.

#include "seed_manager.h"

#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

void check(bool condition, const std::string& description) {
    if (!condition) {
        std::cerr << "FAIL: " << description << "\n";
        ++g_failures;
    } else {
        std::cout << "PASS: " << description << "\n";
    }
}

using ecai_cpp2::seed_manager::EffectiveConfig;
using ecai_cpp2::seed_manager::StreamId;

EffectiveConfig make_config1() {
    EffectiveConfig config;
    config.category_name = "Patient";
    config.i0 = 1;
    config.volume_m3 = 25.0;
    config.mask_efficiency = 0.3;
    config.community_infection_rate = 0.03;
    config.target_probability = 0.001;
    config.duration_hours = 1.0;
    config.percentile = 0.96;
    config.confidence_level = 0.95;
    config.simulations_n = 10000;
    config.replicates_r = 500;
    return config;
}

EffectiveConfig make_config2() {
    EffectiveConfig config;
    config.category_name = "Museum, general";
    config.i0 = 400;
    config.volume_m3 = 3000.0;
    config.mask_efficiency = 0.0;
    config.community_infection_rate = 0.01;
    config.target_probability = 0.001;
    config.duration_hours = 1.0;
    config.percentile = 0.96;
    config.confidence_level = 0.95;
    config.simulations_n = 10000;
    config.replicates_r = 500;
    return config;
}

}  // namespace

int main() {
    using namespace ecai_cpp2::seed_manager;

    // ---------------------------------------------------------------
    // StreamId: unique, exhaustive, name-mapped.
    // ---------------------------------------------------------------
    {
        check(kAllStreamIds.size() == 12, "kAllStreamIds contains exactly 12 entries");

        std::set<std::uint64_t> raw_values;
        for (StreamId id : kAllStreamIds) {
            raw_values.insert(static_cast<std::uint64_t>(id));
        }
        check(raw_values.size() == kAllStreamIds.size(), "every StreamId raw value is unique");

        const std::uint64_t expected_values[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
        for (int i = 0; i < 12; ++i) {
            check(static_cast<std::uint64_t>(kAllStreamIds[static_cast<std::size_t>(i)]) ==
                      expected_values[i],
                  "kAllStreamIds[" + std::to_string(i) + "] has its frozen fixed value");
        }

        const std::string expected_names[12] = {
            "infector_count", "susceptible_pbr", "gamma", "lambda",
            "infector_pbr",   "c_drop",          "d_bar", "E",
            "GVL",            "VF",              "RTD",   "DK",
        };
        std::set<std::string> names;
        for (std::size_t i = 0; i < kAllStreamIds.size(); ++i) {
            const std::string name = stream_id_name(kAllStreamIds[i]);
            check(name == expected_names[i],
                  "stream_id_name(kAllStreamIds[" + std::to_string(i) + "]) == \"" +
                      expected_names[i] + "\"");
            names.insert(name);
        }
        check(names.size() == kAllStreamIds.size(), "every StreamId name is unique");
    }

    // ---------------------------------------------------------------
    // stream_id_name: rejects a StreamId value not in kAllStreamIds
    // (constructed by casting an out-of-range integer).
    // ---------------------------------------------------------------
    {
        bool threw = false;
        try {
            stream_id_name(static_cast<StreamId>(12));
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "stream_id_name rejects an unrecognized StreamId (raw value 12) "
                     "with std::invalid_argument");

        threw = false;
        try {
            stream_id_name(static_cast<StreamId>(0xFFFFFFFFFFFFFFFFULL));
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "stream_id_name rejects an unrecognized StreamId (raw value UINT64_MAX) "
                     "with std::invalid_argument");
    }

    // ---------------------------------------------------------------
    // fnv1a_64: known-answer values, independently fixed.
    // ---------------------------------------------------------------
    {
        check(fnv1a_64("") == 0xcbf29ce484222325ULL,
              "fnv1a_64(\"\") equals the FNV-1a-64 offset basis");
        check(fnv1a_64("a") == 0xaf63dc4c8601ec8cULL, "fnv1a_64(\"a\") matches the known-answer");
        check(fnv1a_64("CPP2SEED") == 0xecf972d99d663045ULL,
              "fnv1a_64(\"CPP2SEED\") matches the known-answer");
    }

    // ---------------------------------------------------------------
    // splitmix64: known-answer values, independently fixed. The seed==0
    // case additionally cross-checks against the public reference
    // splitmix64.c first output.
    // ---------------------------------------------------------------
    {
        check(splitmix64(0ULL) == 0xe220a8397b1dcdafULL,
              "splitmix64(0) matches the known-answer (also the public splitmix64.c "
              "reference first output)");
        check(splitmix64(1ULL) == 0x910a2dec89025cc1ULL, "splitmix64(1) matches the known-answer");
        check(splitmix64(42ULL) == 0xbdd732262feb6e95ULL,
              "splitmix64(42) matches the known-answer");
        check(splitmix64(0x9E3779B97F4A7C15ULL) == 0x6e789e6aa1b965f4ULL,
              "splitmix64(golden-gamma constant) matches the known-answer");
        check(splitmix64(0xFFFFFFFFFFFFFFFFULL) == 0xe4d971771b652c20ULL,
              "splitmix64(UINT64_MAX) matches the known-answer");

        check(splitmix64(0ULL) == splitmix64(0ULL), "splitmix64 is a pure function: same input, "
                                                      "same output across repeated calls");
    }

    // ---------------------------------------------------------------
    // canonical_config_string: exact serialization strings, fixed field
    // order, ASCII '|' delimiter, integers decimal, doubles %.17g-equivalent.
    // ---------------------------------------------------------------
    {
        const std::string expected1 =
            "Patient|1|25|0.29999999999999999|0.029999999999999999|0.001|1|"
            "0.95999999999999996|0.94999999999999996|10000|500";
        check(canonical_config_string(make_config1()) == expected1,
              "canonical_config_string(config1) matches the known-answer string exactly");

        const std::string expected2 =
            "Museum, general|400|3000|0|0.01|0.001|1|0.95999999999999996|"
            "0.94999999999999996|10000|500";
        check(canonical_config_string(make_config2()) == expected2,
              "canonical_config_string(config2) matches the known-answer string exactly");
    }

    // ---------------------------------------------------------------
    // canonical_config_string / config_fingerprint: negative-zero policy —
    // a field of -0.0 canonicalizes identically to +0.0, so the serialized
    // string and fingerprint are unaffected by the sign of a zero field.
    // ---------------------------------------------------------------
    {
        EffectiveConfig zero = make_config1();
        zero.mask_efficiency = 0.0;
        EffectiveConfig negzero = make_config1();
        negzero.mask_efficiency = -0.0;

        const std::string serialized_zero = canonical_config_string(zero);
        const std::string serialized_negzero = canonical_config_string(negzero);
        check(serialized_zero == serialized_negzero,
              "canonical_config_string treats a -0.0 field identically to +0.0");
        check(serialized_zero.find("|0|0.029999999999999999|") != std::string::npos,
              "canonical_config_string emits a bare '0' for a zero-valued double field, "
              "regardless of the sign of zero");

        check(config_fingerprint(zero) == config_fingerprint(negzero),
              "config_fingerprint is unaffected by the sign of a zero-valued field");
        check(config_fingerprint(zero) == 0xf9b0a06ca885a949ULL,
              "config_fingerprint(config1 with mask_efficiency = 0.0) matches the "
              "known-answer fingerprint");
    }

    // ---------------------------------------------------------------
    // canonical_config_string / config_fingerprint: deterministic rejection
    // of nonfinite double fields (NaN and +/-infinity), and of a category
    // name containing the '|' delimiter (invalid configuration).
    // ---------------------------------------------------------------
    {
        const double nonfinite_values[3] = {
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
        };
        const char* labels[3] = {"NaN", "+infinity", "-infinity"};
        for (int i = 0; i < 3; ++i) {
            EffectiveConfig config = make_config1();
            config.target_probability = nonfinite_values[i];
            bool threw = false;
            try {
                canonical_config_string(config);
            } catch (const std::invalid_argument&) {
                threw = true;
            }
            check(threw, std::string("canonical_config_string rejects a ") + labels[i] +
                             " double field with std::invalid_argument");

            bool fp_threw = false;
            try {
                config_fingerprint(config);
            } catch (const std::invalid_argument&) {
                fp_threw = true;
            }
            check(fp_threw, std::string("config_fingerprint rejects a ") + labels[i] +
                                 " double field with std::invalid_argument");
        }

        EffectiveConfig bad_name = make_config1();
        bad_name.category_name = "Patient|room";
        bool threw_delim = false;
        try {
            canonical_config_string(bad_name);
        } catch (const std::invalid_argument&) {
            threw_delim = true;
        }
        check(threw_delim, "canonical_config_string rejects a category_name containing "
                            "the '|' delimiter with std::invalid_argument");
    }

    // ---------------------------------------------------------------
    // fold_stream_seed: several full-tuple known-answer seeds, independent
    // of EffectiveConfig (raw fingerprint inputs).
    // ---------------------------------------------------------------
    {
        check(fold_stream_seed(1, 1, 100, 0, 0, 10000) == 0x9280e0cf84f0022bULL,
              "fold_stream_seed(base=1, dv=1, fp=100, rep=0, stream=0, len=10000) "
              "matches the known-answer seed");
        check(fold_stream_seed(485, 1, 100, 3, 11, 250) == 0xd7d91197f25503f5ULL,
              "fold_stream_seed(base=485, dv=1, fp=100, rep=3, stream=11, len=250) "
              "matches the known-answer seed");
    }

    // ---------------------------------------------------------------
    // fold_stream_seed: each tuple component (base_seed, design_version,
    // fingerprint, replicate_index, stream_id, list_length) changes the
    // derived seed, holding the other five fixed. Every "changed" value
    // below is itself an independently fixed known-answer value, not merely
    // compared for inequality against the baseline.
    // ---------------------------------------------------------------
    {
        const std::uint64_t baseline = fold_stream_seed(1, 1, 100, 0, 0, 10000);
        check(baseline == 0x9280e0cf84f0022bULL, "fold_stream_seed baseline matches its "
                                                  "known-answer seed");

        const std::uint64_t changed_base_seed = fold_stream_seed(2, 1, 100, 0, 0, 10000);
        check(changed_base_seed == 0x54b28fcabcac4423ULL,
              "fold_stream_seed with base_seed changed matches its known-answer seed");
        check(changed_base_seed != baseline, "changing base_seed changes the derived seed");

        const std::uint64_t changed_design_version = fold_stream_seed(1, 2, 100, 0, 0, 10000);
        check(changed_design_version == 0x40077849ddeba2b0ULL,
              "fold_stream_seed with design_version changed matches its known-answer seed");
        check(changed_design_version != baseline,
              "changing design_version changes the derived seed");

        const std::uint64_t changed_fingerprint = fold_stream_seed(1, 1, 101, 0, 0, 10000);
        check(changed_fingerprint == 0x8bbb2d71437044c6ULL,
              "fold_stream_seed with fingerprint changed matches its known-answer seed");
        check(changed_fingerprint != baseline,
              "changing the configuration fingerprint changes the derived seed");

        const std::uint64_t changed_replicate = fold_stream_seed(1, 1, 100, 1, 0, 10000);
        check(changed_replicate == 0x0fe801b99446b664ULL,
              "fold_stream_seed with replicate_index changed matches its known-answer seed");
        check(changed_replicate != baseline,
              "changing replicate_index changes the derived seed");

        const std::uint64_t changed_stream = fold_stream_seed(1, 1, 100, 0, 1, 10000);
        check(changed_stream == 0x127c7ed594e6204dULL,
              "fold_stream_seed with stream_id changed matches its known-answer seed");
        check(changed_stream != baseline, "changing stream_id changes the derived seed");

        const std::uint64_t changed_length = fold_stream_seed(1, 1, 100, 0, 0, 9999);
        check(changed_length == 0x0749c43e72d34d1eULL,
              "fold_stream_seed with list_length changed matches its known-answer seed");
        check(changed_length != baseline, "changing list_length changes the derived seed");
    }

    // ---------------------------------------------------------------
    // fold_stream_seed: same tuple is deterministic (called twice).
    // ---------------------------------------------------------------
    {
        const std::uint64_t a = fold_stream_seed(485, 1, 12345, 7, 3, 8000);
        const std::uint64_t b = fold_stream_seed(485, 1, 12345, 7, 3, 8000);
        check(a == b, "fold_stream_seed(same tuple) called twice yields an identical seed");
    }

    // ---------------------------------------------------------------
    // SeedManager::derive_seed: several full-tuple known-answer seeds built
    // from real EffectiveConfig values (base_seed=485, default design
    // version), exercising the full config -> fingerprint -> fold pipeline.
    // ---------------------------------------------------------------
    {
        SeedManager manager(485);
        const EffectiveConfig config1 = make_config1();
        const EffectiveConfig config2 = make_config2();

        const std::uint64_t seed_a =
            manager.derive_seed(config1, 0, StreamId::infector_count, 10000);
        check(seed_a == 0xa8dc5d1018e8dc93ULL,
              "SeedManager::derive_seed(config1, rep=0, infector_count, 10000) matches "
              "the known-answer seed");

        const std::uint64_t seed_b =
            manager.derive_seed(config1, 1, StreamId::infector_count, 10000);
        check(seed_b == 0x6d04c5fa311ce8bcULL,
              "SeedManager::derive_seed(config1, rep=1, infector_count, 10000) matches "
              "the known-answer seed");
        check(seed_b != seed_a, "changing only replicate_index changes SeedManager::derive_seed's "
                                 "output");

        const std::uint64_t seed_c =
            manager.derive_seed(config1, 0, StreamId::susceptible_pbr, 9800);
        check(seed_c == 0x2705070f7f765ae3ULL,
              "SeedManager::derive_seed(config1, rep=0, susceptible_pbr, 9800) matches "
              "the known-answer seed");
        check(seed_c != seed_a, "changing stream_id and list_length changes "
                                 "SeedManager::derive_seed's output");

        const std::uint64_t seed_d =
            manager.derive_seed(config2, 0, StreamId::infector_count, 10000);
        check(seed_d == 0x116e350d18a38d84ULL,
              "SeedManager::derive_seed(config2, rep=0, infector_count, 10000) matches "
              "the known-answer seed");
        check(seed_d != seed_a, "changing only the effective configuration changes "
                                 "SeedManager::derive_seed's output");

        // Same tuple is deterministic across repeated calls and independent
        // SeedManager instances with the same base seed.
        SeedManager manager_again(485);
        check(manager_again.derive_seed(config1, 0, StreamId::infector_count, 10000) == seed_a,
              "SeedManager::derive_seed(same tuple) is deterministic across independent "
              "SeedManager instances with the same base seed");
    }

    // ---------------------------------------------------------------
    // SeedManager: an explicit design_version constructor argument changes
    // derive_seed's output relative to the default design version, all
    // else held fixed.
    // ---------------------------------------------------------------
    {
        SeedManager default_version(1);
        SeedManager explicit_version_2(1, 2);
        const EffectiveConfig config1 = make_config1();

        const std::uint64_t seed_default =
            default_version.derive_seed(config1, 0, StreamId::infector_count, 10000);
        const std::uint64_t seed_v2 =
            explicit_version_2.derive_seed(config1, 0, StreamId::infector_count, 10000);
        check(seed_default != seed_v2,
              "SeedManager constructed with a different design_version yields a different "
              "derived seed, all other tuple components held fixed");

        SeedManager explicit_version_1(1, 1);
        check(explicit_version_1.derive_seed(config1, 0, StreamId::infector_count, 10000) ==
                  seed_default,
              "SeedManager(base, kDesignVersion) matches SeedManager(base) (the default "
              "design_version argument)");
    }

    // ---------------------------------------------------------------
    // SeedManager::derive_seed: rejects an unrecognized StreamId before
    // ever consulting the configuration.
    // ---------------------------------------------------------------
    {
        SeedManager manager(1);
        bool threw = false;
        try {
            manager.derive_seed(make_config1(), 0, static_cast<StreamId>(99), 10000);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "SeedManager::derive_seed rejects an unrecognized StreamId with "
                     "std::invalid_argument");
    }

    // ---------------------------------------------------------------
    // SeedManager::derive_seed: rejects an invalid configuration (nonfinite
    // double field), surfacing the same std::invalid_argument as
    // canonical_config_string.
    // ---------------------------------------------------------------
    {
        SeedManager manager(1);
        EffectiveConfig invalid_config = make_config1();
        invalid_config.percentile = std::numeric_limits<double>::quiet_NaN();
        bool threw = false;
        try {
            manager.derive_seed(invalid_config, 0, StreamId::infector_count, 10000);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "SeedManager::derive_seed rejects an invalid (nonfinite-field) "
                     "configuration with std::invalid_argument");
    }

    // ---------------------------------------------------------------
    // No seed collisions across 500 replicates x all 12 defined stream IDs,
    // using a realistic configuration and realistic per-stream list lengths
    // (N for the binomial stream, plausible N_active/K_total-scale lengths
    // for the rest). 500 * 12 = 6000 derived seeds, all pairwise distinct.
    // ---------------------------------------------------------------
    {
        SeedManager manager(485);
        const EffectiveConfig config = make_config1();

        auto list_length_for = [](StreamId id) -> std::uint64_t {
            switch (id) {
                case StreamId::infector_count:
                    return 10000;
                case StreamId::susceptible_pbr:
                case StreamId::gamma:
                case StreamId::lambda:
                    return 9800;  // representative N_active
                default:
                    return 40000;  // representative K_total (aggregate infector demand)
            }
        };

        std::set<std::uint64_t> seeds;
        std::size_t total = 0;
        for (std::uint64_t replicate_index = 0; replicate_index < 500; ++replicate_index) {
            for (StreamId id : kAllStreamIds) {
                const std::uint64_t seed =
                    manager.derive_seed(config, replicate_index, id, list_length_for(id));
                seeds.insert(seed);
                ++total;
            }
        }
        check(total == 500 * 12, "generated exactly 500 * 12 = 6000 seeds for the "
                                  "collision sweep");
        check(seeds.size() == total,
              "no seed collisions across 500 replicates x all 12 defined stream IDs, "
              "using realistic per-stream list lengths");
    }

    // ---------------------------------------------------------------
    // Separately: no two streams in that same 500 x 12 supported set produce
    // an identical first-8-output std::mt19937_64 sequence. Seed uniqueness
    // does not by itself guarantee sequence uniqueness (a finite 64-bit hash
    // is not a mathematical collision-impossibility claim), so this is
    // checked independently, per docs/SEED_POLICY.md.
    // ---------------------------------------------------------------
    {
        SeedManager manager(485);
        const EffectiveConfig config = make_config1();

        auto list_length_for = [](StreamId id) -> std::uint64_t {
            switch (id) {
                case StreamId::infector_count:
                    return 10000;
                case StreamId::susceptible_pbr:
                case StreamId::gamma:
                case StreamId::lambda:
                    return 9800;
                default:
                    return 40000;
            }
        };

        std::set<std::vector<std::uint64_t>> prefixes;
        std::size_t total = 0;
        for (std::uint64_t replicate_index = 0; replicate_index < 500; ++replicate_index) {
            for (StreamId id : kAllStreamIds) {
                const std::uint64_t seed =
                    manager.derive_seed(config, replicate_index, id, list_length_for(id));
                std::mt19937_64 engine(seed);
                std::vector<std::uint64_t> prefix;
                prefix.reserve(8);
                for (int i = 0; i < 8; ++i) {
                    prefix.push_back(engine());
                }
                prefixes.insert(std::move(prefix));
                ++total;
            }
        }
        check(total == 500 * 12, "generated exactly 500 * 12 = 6000 mt19937_64 streams for "
                                  "the sequence-prefix sweep");
        check(prefixes.size() == total,
              "no two streams in the supported 500-replicate x 12-stream set share an "
              "identical first-8-output std::mt19937_64 sequence");
    }

    if (g_failures > 0) {
        std::cerr << g_failures << " test(s) failed.\n";
        return 1;
    }
    std::cout << "All tests passed.\n";
    return 0;
}
