# C++2 Seed Policy

## Scope

This policy defines deterministic seed derivation for the C++2 **aggregate-demand-conditioned marginal one-dimensional LHS** design. It is a C++2 engineering decision, not a Jones et al. scientific assumption. The bank construction, assignment order, and non-joint limitations are defined in [`SAMPLING_DESIGN.md`](SAMPLING_DESIGN.md).

## Fixed stream identifiers

The numeric assignments are permanent and must not be renumbered or reused:

| Value | Stream |
|---:|---|
| 0 | `infector_count` |
| 1 | `susceptible_pbr` |
| 2 | `gamma` |
| 3 | `lambda` |
| 4 | `infector_pbr` |
| 5 | `c_drop` |
| 6 | `d_bar` |
| 7 | `E` |
| 8 | `GVL` |
| 9 | `VF` |
| 10 | `RTD` |
| 11 | `DK` |

A removed stream leaves its identifier reserved. An unrecognized identifier is rejected.

## Effective-configuration serialization

The effective configuration is serialized as one ASCII string with `|` delimiters and this fixed field order:

1. canonical publication category name;
2. `I0`;
3. volume in m³;
4. mask efficiency;
5. community infection rate;
6. target infection probability;
7. duration in hours;
8. percentile as a fraction;
9. confidence level as a fraction;
10. simulations per replicate;
11. replicate count.

Integers use unsigned or signed decimal notation according to their field type. Doubles use the locale-independent `std::to_chars` general format with precision 17, equivalent to C-locale `%.17g` for the supported toolchain. Negative zero is normalized to positive zero and serialized as `0`. NaN and positive or negative infinity are rejected. A category name containing `|` is rejected. The canonical category name itself is already resolved by the strict category-data lookup before this serialization is used.

The configuration fingerprint is 64-bit FNV-1a over the serialized bytes:

- offset basis: `0xcbf29ce484222325`;
- prime: `0x100000001b3`;
- for each byte in order: XOR the byte into the state, then multiply by the prime with defined unsigned 64-bit wraparound.

`std::hash` is prohibited because its result is not specified to be stable across implementations or releases.

## SplitMix64 transform

`splitmix64(x)` is the following exact unsigned 64-bit transform, including the initial addition:

    x += 0x9E3779B97F4A7C15
    z = x
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
    z = (z ^ (z >> 27)) * 0x94D049BB133111EB
    z = z ^ (z >> 31)
    return z

All additions and multiplications use defined modulo-2⁶⁴ unsigned arithmetic.

## Stream-seed derivation

The domain tag is `0x4350503253454544` (`CPP2SEED`). For base seed `b`, initialize:

    state = splitmix64(b ^ 0x4350503253454544)

Then fold this ordered tuple:

    (design_version,
     effective_configuration_fingerprint,
     replicate_index,
     stream_id,
     list_length)

For each tuple element `e`, apply exactly:

    state = splitmix64(state ^ splitmix64(e))

The final state seeds that stream's `std::mt19937_64`. No component may derive seeds independently with `seed + index`, and no alternate additive or XOR fold is permitted.

The current sampling-design version is `1`. Any future change that alters stochastic stream construction, assignment, or consumption must increment the design version after review.

Every serialized field is seed-relevant. In particular, changing the confidence level changes the effective-configuration fingerprint and therefore all model-stream seeds, even though confidence level is otherwise a reporting setting. Comparisons intended to reuse identical stochastic samples must keep the complete effective configuration, including confidence level, unchanged.

`std::random_device` is prohibited for model streams. Reproducibility begins with the explicit user/default base seed.

## Verified uniqueness scope

Tests exhaustively cover the supported default study domain of 500 replicate indices and all 12 stream identifiers for a representative effective configuration and realistic stream lengths. They separately verify:

- all 6,000 derived seeds are distinct;
- all 6,000 first-eight-output `std::mt19937_64` prefixes are distinct;
- changing each tuple component changes the checked derived seed;
- repeated identical tuples produce identical seeds.

These are empirical supported-domain checks, not a mathematical proof that a finite 64-bit derivation can never collide. Seed uniqueness and generated-sequence uniqueness are separate tested properties.

## Probability mode

Probability mode (`--mode probability`, `include/simulation.h` `OutputMode::probability`) evaluates the forward Wells-Riley probability at a fixed input ECAi instead of inverting for ECAi. It reuses this seed policy unchanged:

- `mode` and the resolved `input_ecai_l_s_person` are **not** fields of `EffectiveConfig` and are not serialized into the effective-configuration string or folded into the fingerprint. Both modes therefore derive identical stream seeds, banks, replicate assignments, and cursor positions for the same `(category, overrides, N, R, percentile, confidence, target, duration, base seed)`.
- `target_probability` remains in the fingerprint even though probability mode uses it only as an output-stage comparator against `raw_mean_probability`, never as a model input. Consequently sweeping `--target-probability` in probability mode changes the derived seeds and samples exactly as it does in ECAi mode; it is not a same-sample sweep.
- `--input-ecai` is not a fingerprint field, so sweeping it reuses one sample set under one seed. `--mode` itself is not a fingerprint field either, so `--mode ecai` explicit derives byte-identical seeds to omitting `--mode`.
- The resolved `OutputMode` (`ReplicateResult::mode`, set from the resolved `EffectiveScenario::mode`) determines only which deterministic per-row calculation and output fields are populated (`probability`/`percentile_probability` vs. `ecai_l_s_person`/`percentile_l_s_person`); it does not alter stream construction, assignment order, or QER consumption, which are identical between modes.
- The sampling-design version is not incremented for probability mode: it changes the deterministic per-row map, not stochastic stream construction, assignment, or consumption.
- An explicit `--input-ecai` and any override of the four category scientific inputs mark the run `custom_scenario` and suppress the comparator, but do not change fingerprint behavior beyond what those overrides already did in ECAi mode.

See `docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md` §6 for the full rationale and cross-mode test coverage.

## Probability mode

`--mode probability` evaluates the forward Wells-Riley probability at a fixed input ECAi instead
of inverting for ECAi. Neither `mode` nor `input_ecai_l_s_person` is a field of `EffectiveConfig`
and neither is serialized into the fingerprint. Consequently, for an identical
`(category, overrides, N, R, percentile, confidence, target, duration, base seed)`, ECAi mode and
probability mode derive identical stream seeds and therefore identical `infector_counts`,
parameter banks, per-replicate assignment order, and cursor state. Every susceptible and infector
row consumes the same QER and other sampled values in the same order in both modes; only the
deterministic per-row formula applied to those already-sampled values differs. The sampling-design
version stays `1`: probability mode is a new consumer of the existing streams, not a change to
stream construction, assignment, or consumption order.

Dispatch onto the two modes is by the *resolved* `scenario.mode`, never by the caller-supplied
`overrides.mode` alone. `run_replicate_from_banks` sets `result.mode = scenario.mode`, and every
downstream consumer (`run_replicates`, `main.cpp`, and both `output.cpp` document builders)
re-checks that `overrides.mode`, `scenario.mode`, and every replicate's `mode` agree before
proceeding, rejecting the run otherwise.

`target_probability` remains a fingerprint field in both modes. In probability mode the target is
only a comparator applied to `raw_mean_probability` at the output stage; it is never read by the
forward evaluation itself (`ProbabilityInputs` carries the seven physical inputs and the fixed
input ECAi only, no target). Because the target is still fingerprinted, an explicit
`--target-probability` in probability mode changes the derived seeds and therefore the sample, so
sweeping the target across runs is not a same-sample sweep — comparing several thresholds against
one sample set requires holding the target at its default and comparing `raw_percentile_probability`
values offline. An explicit `--target-probability` also marks the run custom, like any other
scientific-input override, and suppresses the comparator, so the only comparator-eligible target is
the frozen default. `--input-ecai` and `--mode` are not fingerprint fields, so sweeping the input
ECAi, or toggling an explicit `--mode ecai` against omitting `--mode` entirely, reuses one sample
set.

## Portability boundary

The fixed integer seed fold, `std::mt19937_64` engine outputs, threshold-rejection bounded integers, and explicit Fisher–Yates index permutation use standard-defined integer behavior and explicit algorithms. Bit identity for integer-to-open-double conversion, `nextafter` corrections, binary64 stratum arithmetic, or shuffled double values additionally requires supported IEEE-754 binary64 semantics and the same relevant compiler, standard and mathematical libraries, architecture, and floating-point flags; C++17 conformance alone is not sufficient.

That claim does not extend automatically to inverse distribution transforms or ECAi calculations because they use platform mathematical-library functions. End-to-end bit reproducibility requires the frozen compiler, standard library, mathematical library, architecture, and floating-point flags. On other supported toolchains, transformed values and model outputs must agree within the tolerances in `MODEL_SPECIFICATION.md` §11.5. Fast-math and uncontrolled floating-point contraction are prohibited.

## Infector policy fingerprinting

`community` remains the default and its canonical seed fingerprint/design version is unchanged. `risk-averse` is resolved separately from output mode and is included in risk-averse provenance; it must not be inferred from raw CLI flags. Runs with the same resolved risk-averse scenario, banks, assignments, and seeds are reproducible subject to the portability boundary above. Changing other fingerprinted scientific inputs, including target probability, can legitimately change the stream.