// Stage 12: strict command-line parsing and run resolution.
// See include/cli.h for the contract this file implements.

#include "cli.h"

#include <charconv>
#include <cmath>
#include <cstddef>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include "qer_direct.h"
#include "sampling_plan.h"
#include "statistics.h"

namespace ecai_cpp2 {
namespace cli {

namespace {

// Every diagnostic is a single line with no trailing newline, so a caller can
// print it verbatim on one line next to a "Try --help" hint.
[[noreturn]] void reject(const std::string& message) {
    std::string single_line;
    single_line.reserve(message.size());
    for (const unsigned char character : message) {
        switch (character) {
            case '\n': single_line += "\\n"; break;
            case '\r': single_line += "\\r"; break;
            default:
                if (character < 0x20 || character == 0x7f) {
                    single_line += '?';
                } else {
                    single_line.push_back(static_cast<char>(character));
                }
                break;
        }
    }
    throw std::invalid_argument(single_line);
}

const std::string& max_simulations_text() {
    static const std::string text =
        std::to_string(sampling_plan::ResourcePolicy{}.max_simulations);
    return text;
}

const std::string& max_replicates_text() {
    static const std::string text = std::to_string(statistics::kMaxSupportedReplicateCount);
    return text;
}

// Strict, locale-independent, full-string numeric parsing. std::from_chars
// never consults the locale, never skips leading whitespace, and never
// accepts a leading '+'; anything it leaves unconsumed is a trailing
// character and is rejected here.
double parse_double_strict(const std::string& text, const std::string& flag) {
    if (text.empty()) {
        reject("flag " + flag + " requires a nonempty numeric value");
    }
    double value = 0.0;
    const char* const first = text.data();
    const char* const last = first + text.size();
    const std::from_chars_result result = std::from_chars(first, last, value);
    if (result.ec == std::errc::result_out_of_range) {
        reject("flag " + flag + " value \"" + text +
               "\" cannot be represented as a double (overflow or underflow)");
    }
    if (result.ec != std::errc{}) {
        reject("flag " + flag + " expects a number, but got \"" + text + "\"");
    }
    if (result.ptr != last) {
        reject("flag " + flag + " expects a complete number with no trailing characters, but got \"" +
               text + "\"");
    }
    if (!std::isfinite(value)) {
        reject("flag " + flag + " expects a finite number, but got \"" + text + "\"");
    }
    return value;
}

std::uint64_t parse_uint64_strict(const std::string& text, const std::string& flag) {
    if (text.empty()) {
        reject("flag " + flag + " requires a nonempty integer value");
    }
    if (text.front() == '-') {
        reject("flag " + flag + " expects an unsigned integer, but got the signed value \"" + text +
               "\"");
    }
    std::uint64_t value = 0;
    const char* const first = text.data();
    const char* const last = first + text.size();
    const std::from_chars_result result = std::from_chars(first, last, value, 10);
    if (result.ec == std::errc::result_out_of_range) {
        reject("flag " + flag + " value \"" + text + "\" does not fit in an unsigned 64 bit integer");
    }
    if (result.ec != std::errc{}) {
        reject("flag " + flag + " expects a decimal integer, but got \"" + text + "\"");
    }
    if (result.ptr != last) {
        reject("flag " + flag +
               " expects a whole integer with no trailing characters, but got \"" + text + "\"");
    }
    return value;
}

long long parse_longlong_strict(const std::string& text, const std::string& flag) {
    if (text.empty()) {
        reject("flag " + flag + " requires a nonempty integer value");
    }
    long long value = 0;
    const char* const first = text.data();
    const char* const last = first + text.size();
    const std::from_chars_result result = std::from_chars(first, last, value, 10);
    if (result.ec == std::errc::result_out_of_range) {
        reject("flag " + flag + " value \"" + text + "\" does not fit in a 64 bit signed integer");
    }
    if (result.ec != std::errc{}) {
        reject("flag " + flag + " expects a decimal integer, but got \"" + text + "\"");
    }
    if (result.ptr != last) {
        reject("flag " + flag +
               " expects a whole integer with no trailing characters, but got \"" + text + "\"");
    }
    return value;
}

std::string require_identifier(const std::string& text, const std::string& flag) {
    if (text.empty()) {
        reject("flag " + flag + " requires a nonempty value");
    }
    return text;
}

bool is_accepted_flag(const std::string& token) {
    for (const std::string& flag : accepted_flags()) {
        if (flag == token) return true;
    }
    return false;
}

// Consumes the value that follows args[index]. A value that begins with "--"
// is always a mistake: no accepted flag takes a value in that form.
const std::string& take_value(const std::vector<std::string>& args, std::size_t& index,
                               const std::string& flag) {
    if (index + 1 >= args.size()) {
        reject("flag " + flag + " requires a value");
    }
    const std::string& value = args[index + 1];
    if (value.rfind("--", 0) == 0) {
        reject("flag " + flag + " requires a value, but the next argument is \"" + value + "\"");
    }
    ++index;
    return value;
}

// Every domain check that can be made without the category record. Run once,
// after the whole argument vector has been consumed, so a rejection never
// depends on flag order.
void validate_domains(const CliOptions& options) {
    if (!options.target_probabilities.empty()) {
        if (options.target_probabilities.size() > 64 ||
            options.target_probabilities.front() != options.target_probability ||
            !options.target_probability_explicit ||
            options.mode != simulation::OutputMode::ecai || !options.raw_output_directory ||
            options.output_directory || options.summary_file || options.replicates_file ||
            options.write_replicates) {
            reject("paired targets require ECAi mode and --raw-output-directory, without other output options");
        }
        std::set<double> unique;
        for (double value : options.target_probabilities) {
            if (!std::isfinite(value) || !(value > 0.0 && value < 1.0) ||
                !unique.insert(value).second) {
                reject("paired targets must be unique finite probabilities strictly between 0 and 1");
            }
        }
    } else if (options.raw_output_directory) {
        reject("--raw-output-directory requires --target-probabilities");
    }
    if (options.simulations_n < 1) {
        reject("flag --simulations requires N >= 1");
    }
    if (options.simulations_n > sampling_plan::ResourcePolicy{}.max_simulations) {
        reject("flag --simulations requires N <= " + max_simulations_text() +
               " (the frozen resource policy maximum)");
    }
    if (options.replicates_r < 1) {
        reject("flag --replicates requires R >= 1");
    }
    if (options.replicates_r > statistics::kMaxSupportedReplicateCount) {
        reject("flag --replicates requires R <= " + max_replicates_text() +
               " (the supported Student-t degrees-of-freedom range)");
    }
    if (!(options.percentile >= 0.0) || !(options.percentile <= 1.0)) {
        reject("flag --percentile requires 0 <= P <= 1");
    }
    if (!(options.confidence_level > 0.0) || !(options.confidence_level < 1.0)) {
        reject("flag --confidence requires 0 < LEVEL < 1");
    }
    if (!(options.target_probability > 0.0) || !(options.target_probability < 1.0)) {
        reject("flag --target-probability requires 0 < P < 1");
    }
    if (!(options.duration_hours > 0.0)) {
        reject("flag --duration-hours requires a finite duration > 0");
    }
    if (options.occupants && *options.occupants < 1) {
        reject("flag --occupants requires I0 >= 1");
    }
    if (options.volume_m3 && !(*options.volume_m3 > 0.0)) {
        reject("flag --volume-m3 requires a finite volume > 0");
    }
    if (options.community_infection_rate &&
        (!(*options.community_infection_rate >= 0.0) ||
         !(*options.community_infection_rate <= 1.0))) {
        reject("flag --community-rate requires 0 <= RATE <= 1");
    }
    if (options.mask_efficiency &&
        (!(*options.mask_efficiency >= 0.0) || !(*options.mask_efficiency <= 1.0))) {
        reject("flag --mask-efficiency requires 0 <= FRACTION <= 1");
    }
    if (options.input_ecai_l_s_person) {
        if (!(*options.input_ecai_l_s_person >= 0.0)) {
            reject("flag --input-ecai requires a finite value >= 0");
        }
        if (options.mode != simulation::OutputMode::probability) {
            reject("flag --input-ecai requires --mode probability");
        }
    }
    // Risk-averse infector policy (docs/RISK_AVERSE_IMPLEMENTATION_PLAN.md
    // Phase 3): reject an explicit --community-rate <= 0 combined with
    // --infector-policy risk-averse as early as possible, without a category
    // lookup. An unsupplied --community-rate still resolves to the category
    // default, which is checked authoritatively later in
    // simulation::resolve_scenario (validate_scenario), once that default is
    // known; this early check only ever catches an explicit conflict.
    if (options.infector_policy == simulation::InfectorPolicy::risk_averse &&
        options.community_infection_rate && !(*options.community_infection_rate > 0.0)) {
        reject("flag --infector-policy risk-averse requires --community-rate > 0");
    }
    // QER family/profile selector (docs/QER_DISTRIBUTION_IMPLEMENTATION_PLAN.md
    // Phase 5). Validate the enum domain itself first, and as early as
    // possible: parse_cli can only ever produce one of the three named
    // enumerators, so this only ever catches options built or edited out of
    // band, but it must run before the category lookup, scenario
    // resolution, allocation or output creation in resolve_run, not merely
    // be caught eventually by the descriptor-build switch further down.
    switch (options.qer_family) {
        case qer_direct::QerFamily::Jones:
        case qer_direct::QerFamily::JonesFitted:
        case qer_direct::QerFamily::Mikszewski:
            break;
        default:
            reject("QER family selector is out of range; expected Jones, JonesFitted, or Mikszewski");
    }
    // Mikszewski always carries a resolved profile once parsing accepts it;
    // Jones and JonesFitted never do. A mismatch here means the options were
    // edited out of band rather than produced by parse_cli.
    const bool has_mikszewski_profile = options.qer_mikszewski_profile.has_value();
    if ((options.qer_family == qer_direct::QerFamily::Mikszewski) != has_mikszewski_profile) {
        reject("QER family/profile provenance is inconsistent: Mikszewski requires exactly one "
               "resolved profile, and no other family may carry one");
    }
}

// Require non-explicit scalar values to retain their parser defaults. This
// prevents later output from attaching false provenance to an options object
// that was edited after parsing.
void validate_scalar_provenance(const CliOptions& options) {
    if (!options.simulations_explicit && options.simulations_n != kDefaultSimulations) {
        reject("simulation count differs from its unoverridden default");
    }
    if (!options.replicates_explicit && options.replicates_r != kDefaultReplicates) {
        reject("replicate count differs from its unoverridden default");
    }
    if (!options.confidence_explicit && options.confidence_level != kDefaultConfidenceLevel) {
        reject("confidence level differs from its unoverridden default");
    }
    if (!options.percentile_explicit && options.percentile != kDefaultPercentile) {
        reject("percentile differs from its unoverridden default");
    }
    if (!options.target_probability_explicit &&
        options.target_probability != kDefaultTargetProbability) {
        reject("target probability differs from its unoverridden default");
    }
    if (!options.duration_explicit && options.duration_hours != kDefaultDurationHours) {
        reject("duration differs from its unoverridden default");
    }
    if (!options.base_seed_explicit && options.base_seed != kDefaultBaseSeed) {
        reject("base seed differs from its unoverridden default");
    }
    if (!options.mode_explicit && options.mode != kDefaultMode) {
        reject("output mode differs from its unoverridden default");
    }
    if (!options.infector_policy_explicit &&
        options.infector_policy != simulation::InfectorPolicy::community) {
        reject("infector policy differs from its unoverridden default");
    }
    if (!options.qer_family_explicit &&
        (options.qer_family != qer_direct::QerFamily::Jones ||
         options.qer_mikszewski_profile.has_value())) {
        reject("QER distribution selector differs from its unoverridden default");
    }
}

}  // namespace

const std::vector<std::string>& accepted_flags() {
    static const std::vector<std::string> flags = {
        "--category",       "--simulations",       "--replicates",
        "--percentile",     "--confidence",        "--target-probability",
        "--duration-hours", "--base-seed",         "--community-rate",
        "--occupants",      "--volume-m3",         "--mask-efficiency",
        "--mode",           "--input-ecai",        "--infector-policy",
        "--qer-distribution",
        "--target-probabilities", "--raw-output-directory",
        "--output-directory", "--summary-file",    "--replicates-file",
        "--write-replicates", "--help"};
    return flags;
}

std::string help_text() {
    return
        "ecai_cpp2 (C++) : CLI-only implementation of the ASHRAE 241 aggregate-demand\n"
        "1D-LHS ECAi model of Jones et al. (2025). Independent of, and not linked\n"
        "against, the legacy C++_legacy implementation.\n"
        "\n"
        "Usage:\n"
        "  ecai_cpp2 --category NAME [options]\n"
        "  ecai_cpp2 --help\n"
        "\n"
        "One occupancy category per invocation. --category is required; there is no\n"
        "default category.\n"
        "\n"
        "Category selection:\n"
        "  --category NAME             Occupancy category, given as a full publication\n"
        "                              name or as one of the 10 frozen legacy aliases.\n"
        "                              Matching is case insensitive and treats spaces,\n"
        "                              commas and hyphens as interchangeable delimiters.\n"
        "                              No fuzzy, prefix or typo tolerant matching.\n"
        "\n"
        "Run controls (these never make a run a custom scenario):\n"
        "  --simulations N             Simulations per replicate. Integer,\n"
        "                              1 <= N <= " + max_simulations_text() + ". Default 10000.\n"
        "  --replicates R              Replicates. Integer,\n"
        "                              1 <= R <= " + max_replicates_text() +
        " (the supported Student-t\n"
        "                              degrees-of-freedom range, R - 1). Default 500.\n"
        "  --base-seed SEED            Unsigned 64 bit base seed. Default 485.\n"
        "  --confidence LEVEL          Confidence level for the across-replicate\n"
        "                              Student-t interval, 0 < LEVEL < 1. Default 0.95.\n"
        "\n"
        "Scientific inputs (any explicit value here marks the run a custom scenario):\n"
        "  --percentile P              Percentile of the ECAi distribution taken within\n"
        "                              each replicate, 0 <= P <= 1. Default 0.96.\n"
        "  --target-probability P      Target infection probability over the effective\n"
        "                              duration, 0 < P < 1. Default 0.001.\n"
        "  --target-probabilities LIST Comma-separated distinct targets, 0 < P < 1.\n"
        "                              The first determines the existing seed fingerprint.\n"
        "                              Requires --raw-output-directory; ECAi mode only.\n"
        "  --duration-hours HOURS      Exposure duration in hours, finite and > 0.\n"
        "                              Default 1. The target probability always applies\n"
        "                              over the effective duration and is an hourly\n"
        "                              target only when that duration is exactly 1 hour.\n"
        "  --community-rate RATE       Community infection rate, 0 <= RATE <= 1.\n"
        "                              Default: the selected category value.\n"
        "  --occupants I0              Maximum occupants, integer >= 1.\n"
        "                              Default: the selected category value.\n"
        "  --volume-m3 VOLUME          Space volume in m^3, finite and > 0.\n"
        "                              Default: the selected category value.\n"
        "  --mask-efficiency FRACTION  Mask efficiency, 0 <= FRACTION <= 1.\n"
        "                              Default: the selected category value.\n"
        "\n"
        "Output mode (--mode never marks the run a custom scenario by itself):\n"
        "  --mode MODE                 \"ecai\" (default) or \"probability\", case sensitive,\n"
        "                              no aliases. ecai mode inverts for the ECAi that\n"
        "                              reproduces --target-probability, exactly as before.\n"
        "                              probability mode instead evaluates the forward\n"
        "                              infection probability at a fixed input ECAi.\n"
        "  --input-ecai ECAI           Fixed input ECAi for probability mode, in\n"
        "                              L/s/person, finite and >= 0. Requires\n"
        "                              --mode probability; an explicit value marks the\n"
        "                              run a custom scenario. Default: the selected\n"
        "                              category's published Table 4 ECAi.\n"
        "\n"
        "Infector policy (--infector-policy never marks the run a custom scenario by\n"
        "itself, and applies identically to both output modes above):\n"
        "  --infector-policy POLICY    \"community\" (default) or \"risk-averse\", case\n"
        "                              sensitive, no aliases. community retains the\n"
        "                              unconditional per-simulation binomial infector\n"
        "                              count, including zero-infector rows. risk-averse\n"
        "                              conditions every simulation on at least one\n"
        "                              infector and requires a strictly positive\n"
        "                              effective community rate (I0 >= 1 is already\n"
        "                              required in both policies).\n"
        "\n"
        "QER distribution selector (docs/QER_DISTRIBUTION_IMPLEMENTATION_PLAN.md; any explicit\n"
        "value here, including an explicit \"jones\", marks the run a custom scenario, because\n"
        "it changes the QER calculation mechanism):\n"
        "  --qer-distribution NAME     \"jones\" (default, full mechanistic Jones et al.\n"
        "                              calculation), \"jones-fitted\", or a full canonical\n"
        "                              Mikszewski profile name, e.g.\n"
        "                              \"mikszewski-sars-cov-2-standing-speaking\". Bare\n"
        "                              \"mikszewski\" is rejected: no default profile is\n"
        "                              resolved for it. Case sensitive, no aliases. Selecting\n"
        "                              a direct family selects the family-aware sampling and\n"
        "                              simulation path and emits a versioned provenance\n"
        "                              schema that is not a full-Jones comparator output.\n"
        "\n"
        "Output controls (Stage 12 parses and records these only; no directory or file\n"
        "is created, opened, written or checked until Stage 13):\n"
        "  --output-directory PATH     Directory for Stage 13 output.\n"
        "  --raw-output-directory PATH New directory for paired per-simulation data.\n"
        "  --summary-file PATH         Summary file path for Stage 13 output.\n"
        "  --replicates-file PATH      Per-replicate file path for Stage 13 output.\n"
        "  --write-replicates          Request per-replicate output in Stage 13.\n"
        "\n"
        "  --help                      Print this help and exit. Must be the only\n"
        "                              argument; it runs no lookup and no simulation.\n"
        "\n"
        "Parsing rules:\n"
        "  Numeric values are parsed strictly, locale independently, over the whole\n"
        "  string. Empty values, leading or trailing characters or whitespace, NaN,\n"
        "  infinity, values that overflow or underflow their field, negative or\n"
        "  non-integer values for integer fields, duplicate flags, unknown flags,\n"
        "  positional arguments and missing values are all rejected. Every rejection\n"
        "  exits nonzero with a single line diagnostic and produces no scientific\n"
        "  output.\n"
        "\n"
        "Explicit category overrides win directly over the selected category defaults.\n"
        "The category source defaults are retained separately from the effective\n"
        "values, and a run using any explicit scientific override is reported as a\n"
        "custom scenario rather than as a published category run.\n";
}

CliOptions parse_cli(const std::vector<std::string>& args) {
    for (const std::string& arg : args) {
        if (arg == "--help") {
            if (args.size() != 1) {
                reject("--help must be the only argument");
            }
            CliOptions help_options;
            help_options.help_requested = true;
            return help_options;
        }
    }

    CliOptions options;
    std::set<std::string> seen;

    for (std::size_t index = 0; index < args.size(); ++index) {
        const std::string& arg = args[index];
        if (arg.rfind("--", 0) != 0) {
            if (!arg.empty() && arg.front() == '-') {
                reject("unknown flag \"" + arg + "\"");
            }
            reject("unexpected positional argument \"" + arg + "\"");
        }
        if (!is_accepted_flag(arg)) {
            reject("unknown flag \"" + arg + "\"");
        }
        if (!seen.insert(arg).second) {
            reject("flag " + arg + " was supplied more than once");
        }

        if (arg == "--write-replicates") {
            options.write_replicates = true;
            continue;
        }

        const std::string value = take_value(args, index, arg);
        if (arg == "--category") {
            options.category_name = require_identifier(value, arg);
            options.category_supplied = true;
        } else if (arg == "--simulations") {
            options.simulations_n = parse_uint64_strict(value, arg);
            options.simulations_explicit = true;
        } else if (arg == "--replicates") {
            options.replicates_r = parse_uint64_strict(value, arg);
            options.replicates_explicit = true;
        } else if (arg == "--percentile") {
            options.percentile = parse_double_strict(value, arg);
            options.percentile_explicit = true;
        } else if (arg == "--confidence") {
            options.confidence_level = parse_double_strict(value, arg);
            options.confidence_explicit = true;
        } else if (arg == "--target-probability") {
            options.target_probability = parse_double_strict(value, arg);
            options.target_probability_explicit = true;
        } else if (arg == "--target-probabilities") {
            std::size_t start = 0;
            while (true) {
                const auto end = value.find(',', start);
                options.target_probabilities.push_back(parse_double_strict(
                    value.substr(start, end == std::string::npos ? end : end - start), arg));
                if (end == std::string::npos) break;
                start = end + 1;
            }
        } else if (arg == "--raw-output-directory") {
            options.raw_output_directory = require_identifier(value, arg);
        } else if (arg == "--duration-hours") {
            options.duration_hours = parse_double_strict(value, arg);
            options.duration_explicit = true;
        } else if (arg == "--base-seed") {
            options.base_seed = parse_uint64_strict(value, arg);
            options.base_seed_explicit = true;
        } else if (arg == "--community-rate") {
            options.community_infection_rate = parse_double_strict(value, arg);
        } else if (arg == "--occupants") {
            options.occupants = parse_longlong_strict(value, arg);
        } else if (arg == "--volume-m3") {
            options.volume_m3 = parse_double_strict(value, arg);
        } else if (arg == "--mask-efficiency") {
            options.mask_efficiency = parse_double_strict(value, arg);
        } else if (arg == "--mode") {
            if (value == "ecai") {
                options.mode = simulation::OutputMode::ecai;
            } else if (value == "probability") {
                options.mode = simulation::OutputMode::probability;
            } else {
                reject("flag --mode expects \"ecai\" or \"probability\", but got \"" + value +
                       "\"");
            }
            options.mode_explicit = true;
        } else if (arg == "--input-ecai") {
            options.input_ecai_l_s_person = parse_double_strict(value, arg);
        } else if (arg == "--infector-policy") {
            if (value == "community") {
                options.infector_policy = simulation::InfectorPolicy::community;
            } else if (value == "risk-averse") {
                options.infector_policy = simulation::InfectorPolicy::risk_averse;
            } else {
                reject("flag --infector-policy expects \"community\" or \"risk-averse\", but got \"" +
                       value + "\"");
            }
            options.infector_policy_explicit = true;
        } else if (arg == "--qer-distribution") {
            if (value == "jones") {
                options.qer_family = qer_direct::QerFamily::Jones;
                options.qer_mikszewski_profile.reset();
            } else if (value == "jones-fitted") {
                options.qer_family = qer_direct::QerFamily::JonesFitted;
                options.qer_mikszewski_profile.reset();
            } else if (value == "mikszewski") {
                reject(
                    "flag --qer-distribution requires an explicit Mikszewski profile, e.g. "
                    "\"mikszewski-sars-cov-2-standing-speaking\"; bare \"mikszewski\" is not "
                    "accepted");
            } else {
                try {
                    options.qer_mikszewski_profile = qer_direct::parse_mikszewski_profile(value);
                } catch (const std::invalid_argument&) {
                    reject("flag --qer-distribution expects \"jones\", \"jones-fitted\", or a full "
                           "Mikszewski profile name, but got \"" +
                           value + "\"");
                }
                options.qer_family = qer_direct::QerFamily::Mikszewski;
            }
            options.qer_family_explicit = true;
        } else if (arg == "--output-directory") {
            options.output_directory = require_identifier(value, arg);
        } else if (arg == "--summary-file") {
            options.summary_file = require_identifier(value, arg);
        } else if (arg == "--replicates-file") {
            options.replicates_file = require_identifier(value, arg);
        } else {
            reject("unknown flag \"" + arg + "\"");
        }
    }

    if (!options.category_supplied) {
        reject("--category NAME is required; there is no default category");
    }

    if (!options.target_probabilities.empty()) {
        if (seen.count("--target-probability")) {
            reject("--target-probability cannot be combined with --target-probabilities");
        }
        options.target_probability = options.target_probabilities.front();
        options.target_probability_explicit = true;
    }

    options.custom_scenario = options.target_probability_explicit || options.duration_explicit ||
                               options.percentile_explicit ||
                               options.community_infection_rate.has_value() ||
                               options.occupants.has_value() || options.volume_m3.has_value() ||
                               options.mask_efficiency.has_value() ||
                               options.input_ecai_l_s_person.has_value() ||
                               options.qer_family_explicit;
    options.target_is_hourly = (options.duration_hours == 1.0);
    validate_domains(options);
    return options;
}

ResolvedRun resolve_run(const CliOptions& options) {
    if (options.help_requested) {
        reject("--help performs no category lookup and no run");
    }
    if (!options.category_supplied || options.category_name.empty()) {
        reject("--category NAME is required; there is no default category");
    }

    validate_scalar_provenance(options);

    // The parser's provenance is re-derived here and must agree: a
    // disagreement means the options were built or edited out of band, and
    // the run is refused rather than mislabelled.
    const bool expected_custom =
        options.target_probability_explicit || options.duration_explicit ||
        options.percentile_explicit || options.community_infection_rate.has_value() ||
        options.occupants.has_value() || options.volume_m3.has_value() ||
        options.mask_efficiency.has_value() || options.input_ecai_l_s_person.has_value() ||
        options.qer_family_explicit;
    if (expected_custom != options.custom_scenario) {
        reject("custom-scenario provenance disagrees with the supplied scientific overrides");
    }
    if ((options.duration_hours == 1.0) != options.target_is_hourly) {
        reject("hourly-target provenance disagrees with the effective duration");
    }
    validate_domains(options);

    ResolvedRun run;
    try {
        run.category = &lookup_occupancy_category(options.category_name);
    } catch (const std::invalid_argument&) {
        reject("unknown category \"" + options.category_name +
               "\"; give a full publication name or a frozen legacy alias (no fuzzy matching)");
    }

    run.base_seed = options.base_seed;
    run.custom_scenario = options.custom_scenario;
    run.target_is_hourly = options.target_is_hourly;

    // Explicit CLI values win directly over the category defaults; the
    // category record itself is retained untouched behind run.category.
    run.overrides.occupants = options.occupants;
    run.overrides.volume_m3 = options.volume_m3;
    run.overrides.mask_efficiency = options.mask_efficiency;
    run.overrides.community_infection_rate = options.community_infection_rate;
    run.overrides.mode = options.mode;
    run.overrides.input_ecai_l_s_person = options.input_ecai_l_s_person;
    run.overrides.infector_policy = options.infector_policy;
    run.overrides.infector_policy_explicit = options.infector_policy_explicit;

    // Build the immutable QER selector descriptor once, entirely through the
    // work package 1 qer_direct API, and copy it into EffectiveOverrides
    // before resolve_scenario runs below: resolve_scenario copies these
    // fields through unchanged into EffectiveScenario (Phase 5), so the
    // descriptor must exist first. Jones-fitted params are looked up by the
    // resolved category's canonical name, not by CLI category text, so any
    // legacy alias still resolves to the correct fitted record.
    run.qer_selection.family = options.qer_family;
    run.qer_selection.mikszewski_profile = options.qer_mikszewski_profile;
    switch (options.qer_family) {
        case qer_direct::QerFamily::Jones:
            run.qer_selection.canonical_name = qer_direct::to_string(qer_direct::QerFamily::Jones);
            break;
        case qer_direct::QerFamily::JonesFitted:
            run.qer_selection.canonical_name =
                qer_direct::to_string(qer_direct::QerFamily::JonesFitted);
            run.qer_selection.direct_params = qer_direct::jones_fitted_params_for_canonical_category(
                run.category->canonical_name);
            break;
        case qer_direct::QerFamily::Mikszewski:
            if (!options.qer_mikszewski_profile) {
                reject("QER family Mikszewski resolved without a profile");
            }
            run.qer_selection.canonical_name =
                qer_direct::canonical_mikszewski_profile_name(*options.qer_mikszewski_profile);
            run.qer_selection.direct_params =
                qer_direct::mikszewski_params_for_profile(*options.qer_mikszewski_profile);
            break;
        default:
            reject("QER family selector is out of range; expected Jones, JonesFitted, or Mikszewski");
    }

    run.overrides.qer_family = run.qer_selection.family;
    run.overrides.qer_family_explicit = options.qer_family_explicit;
    run.overrides.qer_mikszewski_profile = run.qer_selection.mikszewski_profile;
    run.overrides.qer_direct_params = run.qer_selection.direct_params;

    run.config.category_name = run.category->canonical_name;
    run.config.i0 = options.occupants ? *options.occupants
                                      : static_cast<std::int64_t>(run.category->i0);
    run.config.volume_m3 = options.volume_m3 ? *options.volume_m3 : run.category->volume_m3;
    run.config.mask_efficiency =
        options.mask_efficiency ? *options.mask_efficiency : run.category->mask_efficiency;
    run.config.community_infection_rate = options.community_infection_rate
                                              ? *options.community_infection_rate
                                              : run.category->community_infection_rate;
    run.config.target_probability = options.target_probability;
    run.config.duration_hours = options.duration_hours;
    run.config.percentile = options.percentile;
    run.config.confidence_level = options.confidence_level;
    run.config.simulations_n = options.simulations_n;
    run.config.replicates_r = options.replicates_r;

    // Stage 10 re-checks that the fingerprinted configuration and the
    // category/override resolution agree exactly, and applies every scenario
    // domain check, before anything is allocated or sampled.
    run.scenario = simulation::resolve_scenario(*run.category, run.config, run.overrides);

    return run;
}

}  // namespace cli
}  // namespace ecai_cpp2
