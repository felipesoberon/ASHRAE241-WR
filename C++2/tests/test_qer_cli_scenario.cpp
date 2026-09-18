// Work package 2 tests for the QER family/profile CLI selector and its
// immutable resolved-descriptor propagation, per
// docs/QER_DISTRIBUTION_IMPLEMENTATION_PLAN.md Phase 5 (scoped to CLI/
// scenario resolution only: no sampling_plan/simulation/output integration
// is exercised or expected here).
//
// This file is deliberately independent of tests/test_cli.cpp (which
// exercises the full flag surface and the executable command matrix): it
// covers only the new --qer-distribution selector, its strict parsing rules,
// and the ResolvedQerSelection descriptor resolve_run builds from the work
// package 1 qer_direct API.

#include "cli.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "qer_direct.h"

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

using Args = std::vector<std::string>;
using ecai_cpp2::cli::CliOptions;
using ecai_cpp2::cli::ResolvedRun;

// A small, bounded scenario used everywhere a parse must also be resolvable:
// "Exam room" has I0 = 3 and VOL = 41 m^3 (Jones et al. Table 2), matching
// tests/test_cli.cpp's bounded_args().
Args bounded_args() {
    return Args{"--category", "Exam room", "--simulations", "8", "--replicates", "2"};
}

Args with_flag(const std::string& flag, const std::string& value) {
    Args args = bounded_args();
    args.push_back(flag);
    args.push_back(value);
    return args;
}

void expect_parse_reject(const Args& args, const std::string& description) {
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

void expect_resolve_reject(const CliOptions& options, const std::string& description) {
    try {
        ecai_cpp2::cli::resolve_run(options);
    } catch (const std::invalid_argument&) {
        check(true, "rejects " + description);
        return;
    } catch (...) {
        check(false, "rejects with std::invalid_argument: " + description);
        return;
    }
    check(false, "rejects " + description);
}

}  // namespace

int main() {
    using ecai_cpp2::qer_direct::DirectQerParams;
    using ecai_cpp2::qer_direct::MikszewskiActivity;
    using ecai_cpp2::qer_direct::MikszewskiOrganism;
    using ecai_cpp2::qer_direct::MikszewskiProfile;
    using ecai_cpp2::qer_direct::QerFamily;

    // ================================================================
    // accepted_flags(): the selector is registered exactly once.
    // ================================================================
    {
        int occurrences = 0;
        for (const std::string& flag : ecai_cpp2::cli::accepted_flags()) {
            if (flag == "--qer-distribution") ++occurrences;
        }
        check(occurrences == 1, "accepted_flags() lists --qer-distribution exactly once");
    }

    // ================================================================
    // parse_cli: default (no flag) leaves the Jones default untouched and
    // does not mark the run a custom scenario by itself.
    // ================================================================
    {
        const CliOptions options = ecai_cpp2::cli::parse_cli(bounded_args());
        check(options.qer_family == QerFamily::Jones,
              "omitted --qer-distribution defaults to Jones");
        check(!options.qer_family_explicit, "omitted --qer-distribution is not explicit");
        check(!options.qer_mikszewski_profile.has_value(),
              "omitted --qer-distribution resolves no Mikszewski profile");
        check(!options.custom_scenario,
              "omitted --qer-distribution does not mark the run a custom scenario");
    }

    // ================================================================
    // parse_cli: explicit "jones" still marks the run a custom scenario,
    // because it explicitly names the QER mechanism rather than being
    // silently defaulted.
    // ================================================================
    {
        const CliOptions options = ecai_cpp2::cli::parse_cli(with_flag("--qer-distribution", "jones"));
        check(options.qer_family == QerFamily::Jones, "explicit \"jones\" resolves to QerFamily::Jones");
        check(options.qer_family_explicit, "explicit \"jones\" sets qer_family_explicit");
        check(!options.qer_mikszewski_profile.has_value(), "explicit \"jones\" carries no profile");
        check(options.custom_scenario,
              "explicit --qer-distribution jones marks the run a custom scenario");
    }

    // ================================================================
    // parse_cli: "jones-fitted".
    // ================================================================
    {
        const CliOptions options =
            ecai_cpp2::cli::parse_cli(with_flag("--qer-distribution", "jones-fitted"));
        check(options.qer_family == QerFamily::JonesFitted, "\"jones-fitted\" resolves correctly");
        check(options.qer_family_explicit, "\"jones-fitted\" sets qer_family_explicit");
        check(!options.qer_mikszewski_profile.has_value(), "\"jones-fitted\" carries no profile");
        check(options.custom_scenario, "\"jones-fitted\" marks the run a custom scenario");
    }

    // ================================================================
    // parse_cli: a full canonical Mikszewski profile name.
    // ================================================================
    {
        const CliOptions options = ecai_cpp2::cli::parse_cli(
            with_flag("--qer-distribution", "mikszewski-sars-cov-2-standing-speaking"));
        check(options.qer_family == QerFamily::Mikszewski,
              "a full Mikszewski profile name resolves to QerFamily::Mikszewski");
        check(options.qer_family_explicit, "a full Mikszewski profile name sets qer_family_explicit");
        check(options.qer_mikszewski_profile.has_value(),
              "a full Mikszewski profile name resolves a profile");
        if (options.qer_mikszewski_profile) {
            check(options.qer_mikszewski_profile->organism == MikszewskiOrganism::SarsCov2,
                  "the resolved profile has the correct organism");
            check(options.qer_mikszewski_profile->activity == MikszewskiActivity::StandingSpeaking,
                  "the resolved profile has the correct activity");
        }
        check(options.custom_scenario, "a full Mikszewski profile name marks the run a custom scenario");
    }

    // Every one of the 30 canonical Mikszewski profile names round-trips
    // through the CLI parser to the same profile qer_direct itself reports.
    {
        bool all_round_trip = true;
        for (const auto& record : ecai_cpp2::qer_direct::mikszewski_profile_table()) {
            const std::string canonical_name =
                ecai_cpp2::qer_direct::canonical_mikszewski_profile_name(record.profile);
            const CliOptions options =
                ecai_cpp2::cli::parse_cli(with_flag("--qer-distribution", canonical_name));
            if (options.qer_family != QerFamily::Mikszewski ||
                !options.qer_mikszewski_profile.has_value() ||
                options.qer_mikszewski_profile->organism != record.profile.organism ||
                options.qer_mikszewski_profile->activity != record.profile.activity) {
                all_round_trip = false;
                break;
            }
        }
        check(all_round_trip, "every one of the 30 canonical Mikszewski profile names round-trips "
                               "through the CLI parser");
    }

    // ================================================================
    // parse_cli: strict rejections.
    // ================================================================
    expect_parse_reject(with_flag("--qer-distribution", "mikszewski"),
                         "bare \"mikszewski\" with no explicit profile");
    expect_parse_reject(with_flag("--qer-distribution", "Jones"),
                         "case-mismatched \"Jones\" (no case folding)");
    expect_parse_reject(with_flag("--qer-distribution", "jones fitted"),
                         "a space instead of the canonical hyphen");
    expect_parse_reject(with_flag("--qer-distribution", "mikszewski-sars-cov-2"),
                         "a Mikszewski organism with no activity suffix");
    expect_parse_reject(with_flag("--qer-distribution", "mikszewski-not-a-real-organism-resting"),
                         "an unknown Mikszewski organism");
    expect_parse_reject(with_flag("--qer-distribution", ""), "an empty selector value");
    expect_parse_reject(with_flag("--qer-distribution", "not-a-real-family"),
                         "an unrecognized selector value");
    {
        // Duplicate --qer-distribution flag (even with the same value) is a
        // duplicate flag, exactly like every other repeated flag.
        Args args = with_flag("--qer-distribution", "jones-fitted");
        args.push_back("--qer-distribution");
        args.push_back("jones-fitted");
        expect_parse_reject(args, "a repeated --qer-distribution flag");
    }

    // ================================================================
    // resolve_run: the immutable descriptor for Jones (default family).
    // ================================================================
    {
        const CliOptions options = ecai_cpp2::cli::parse_cli(bounded_args());
        const ResolvedRun run = ecai_cpp2::cli::resolve_run(options);
        check(run.qer_selection.family == QerFamily::Jones,
              "resolve_run: default descriptor family is Jones");
        check(run.qer_selection.canonical_name == "jones",
              "resolve_run: default descriptor canonical_name is \"jones\"");
        check(!run.qer_selection.mikszewski_profile.has_value(),
              "resolve_run: default descriptor carries no Mikszewski profile");
    }

    // ================================================================
    // resolve_run: the immutable descriptor for Jones-fitted uses the
    // resolved category's canonical name (via the work package 1 API), and
    // matches jones_fitted_params_for_canonical_category exactly.
    // ================================================================
    {
        const CliOptions options =
            ecai_cpp2::cli::parse_cli(with_flag("--qer-distribution", "jones-fitted"));
        const ResolvedRun run = ecai_cpp2::cli::resolve_run(options);
        check(run.qer_selection.family == QerFamily::JonesFitted,
              "resolve_run: Jones-fitted descriptor family");
        check(run.qer_selection.canonical_name == "jones-fitted",
              "resolve_run: Jones-fitted descriptor canonical_name");
        check(run.category != nullptr, "resolve_run: category is resolved");
        if (run.category != nullptr) {
            const DirectQerParams& expected =
                ecai_cpp2::qer_direct::jones_fitted_params_for_canonical_category(
                    run.category->canonical_name);
            check(run.qer_selection.direct_params.mu_log10 == expected.mu_log10 &&
                      run.qer_selection.direct_params.sigma_log10 == expected.sigma_log10,
                  "resolve_run: Jones-fitted descriptor params match "
                  "jones_fitted_params_for_canonical_category for the resolved category");
        }
    }

    // Every one of the 25 categories resolves a Jones-fitted descriptor
    // without throwing (the category_data.h <-> qer_direct.h canonical-name
    // mapping is complete in both directions).
    {
        bool all_resolved = true;
        for (const auto& record : ecai_cpp2::qer_direct::jones_fitted_table()) {
            try {
                const CliOptions options = ecai_cpp2::cli::parse_cli(
                    Args{"--category", record.canonical_category, "--qer-distribution",
                         "jones-fitted"});
                const ResolvedRun run = ecai_cpp2::cli::resolve_run(options);
                if (run.qer_selection.direct_params.mu_log10 != record.params.mu_log10 ||
                    run.qer_selection.direct_params.sigma_log10 != record.params.sigma_log10) {
                    all_resolved = false;
                    break;
                }
            } catch (const std::exception&) {
                all_resolved = false;
                break;
            }
        }
        check(all_resolved,
              "every one of the 25 canonical categories resolves a matching Jones-fitted descriptor");
    }

    // ================================================================
    // resolve_run: the immutable descriptor for a Mikszewski profile.
    // ================================================================
    {
        const CliOptions options = ecai_cpp2::cli::parse_cli(
            with_flag("--qer-distribution", "mikszewski-influenza-resting"));
        const ResolvedRun run = ecai_cpp2::cli::resolve_run(options);
        check(run.qer_selection.family == QerFamily::Mikszewski,
              "resolve_run: Mikszewski descriptor family");
        check(run.qer_selection.canonical_name == "mikszewski-influenza-resting",
              "resolve_run: Mikszewski descriptor canonical_name round-trips");
        check(run.qer_selection.mikszewski_profile.has_value(),
              "resolve_run: Mikszewski descriptor carries the resolved profile");
        if (run.qer_selection.mikszewski_profile) {
            const DirectQerParams& expected = ecai_cpp2::qer_direct::mikszewski_params_for_profile(
                *run.qer_selection.mikszewski_profile);
            check(run.qer_selection.direct_params.mu_log10 == expected.mu_log10 &&
                      run.qer_selection.direct_params.sigma_log10 == expected.sigma_log10,
                  "resolve_run: Mikszewski descriptor params match mikszewski_params_for_profile");
        }
    }

    // ================================================================
    // resolve_run: determinism. Two independent resolutions of the same
    // options produce an identical descriptor.
    // ================================================================
    {
        const CliOptions options = ecai_cpp2::cli::parse_cli(
            with_flag("--qer-distribution", "mikszewski-measles-light-speaking-loudly"));
        const ResolvedRun first = ecai_cpp2::cli::resolve_run(options);
        const ResolvedRun second = ecai_cpp2::cli::resolve_run(options);
        check(first.qer_selection.family == second.qer_selection.family &&
                  first.qer_selection.canonical_name == second.qer_selection.canonical_name &&
                  first.qer_selection.direct_params.mu_log10 ==
                      second.qer_selection.direct_params.mu_log10 &&
                  first.qer_selection.direct_params.sigma_log10 ==
                      second.qer_selection.direct_params.sigma_log10,
              "resolve_run is deterministic: two resolutions of the same options agree exactly");
    }

    // ================================================================
    // resolve_run: provenance disagreement is refused rather than silently
    // resolved (mirrors the existing scalar-provenance checks).
    // ================================================================
    {
        CliOptions options = ecai_cpp2::cli::parse_cli(bounded_args());
        options.qer_family = QerFamily::JonesFitted;  // edited out of band, not through the parser
        expect_resolve_reject(options,
                               "a QER family edited out of band without qer_family_explicit");
    }
    {
        CliOptions options = ecai_cpp2::cli::parse_cli(bounded_args());
        options.qer_mikszewski_profile =
            MikszewskiProfile{MikszewskiOrganism::Measles, MikszewskiActivity::Resting};
        expect_resolve_reject(
            options, "a Mikszewski profile set without QerFamily::Mikszewski and without "
                     "qer_family_explicit");
    }
    {
        CliOptions options =
            ecai_cpp2::cli::parse_cli(with_flag("--qer-distribution", "jones-fitted"));
        options.custom_scenario = false;  // edited out of band
        expect_resolve_reject(
            options, "custom_scenario cleared out of band after an explicit --qer-distribution");
    }
    {
        // Mikszewski family without a resolved profile: parse_cli can never
        // produce this, but resolve_run's domain check must still catch it.
        CliOptions options = ecai_cpp2::cli::parse_cli(bounded_args());
        options.qer_family = QerFamily::Mikszewski;
        options.qer_family_explicit = true;
        options.custom_scenario = true;
        expect_resolve_reject(options, "QerFamily::Mikszewski without a resolved profile");
    }
    {
        // Out-of-band QerFamily value (not one of the three named
        // enumerators): unreachable via parse_cli, but resolve_run's
        // descriptor build must still fail closed instead of falling
        // through its switch with an empty canonical_name and
        // zero-initialized direct_params.
        CliOptions options = ecai_cpp2::cli::parse_cli(bounded_args());
        options.qer_family = static_cast<QerFamily>(99);
        options.qer_family_explicit = true;
        options.custom_scenario = true;
        expect_resolve_reject(options, "an out-of-band QerFamily enum value");
    }
    {
        // An out-of-band QerFamily must be rejected before the category
        // lookup -- not merely rejected eventually. Pairing it with an
        // equally invalid category name distinguishes the two: if
        // resolve_run validated the family only where the switch in the
        // descriptor-build step currently sits (after
        // lookup_occupancy_category), the category lookup would fail first
        // and this would surface as an "unknown category" message instead.
        CliOptions options = ecai_cpp2::cli::parse_cli(bounded_args());
        options.category_name = "not-a-real-category-xyz";  // also invalid, deliberately
        options.qer_family = static_cast<QerFamily>(99);
        options.qer_family_explicit = true;
        options.custom_scenario = true;
        try {
            ecai_cpp2::cli::resolve_run(options);
            check(false, "rejects an out-of-band QerFamily before category lookup");
        } catch (const std::invalid_argument& error) {
            const std::string message = error.what();
            check(message.find("category") == std::string::npos,
                  "out-of-band QerFamily is rejected before the category lookup runs "
                  "(error must not be the \"unknown category\" message)");
            check(message.find("QER family") != std::string::npos,
                  "out-of-band QerFamily rejection names the QER family as the cause");
        } catch (...) {
            check(false, "rejects an out-of-band QerFamily with std::invalid_argument, "
                          "before category lookup");
        }
    }

    std::cout << (g_failures == 0 ? "All qer_cli_scenario checks passed.\n"
                                   : std::to_string(g_failures) + " qer_cli_scenario check(s) failed.\n");
    return g_failures == 0 ? 0 : 1;
}
