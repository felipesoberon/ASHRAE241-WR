# ecai_cpp2 CLI reference

This document describes the command-line interface as it exists now. It is a
current-state reference for `include/cli.h`, `src/cli.cpp` and `src/main.cpp`,
and it matches the text printed by `--help`.

The CLI parses arguments, resolves one occupancy category and its effective
inputs, runs the replicate batch once, and optionally writes validated summary
and per-replicate CSV files. Without output controls it prints only one console
summary line.

## Invocation

```text
ecai_cpp2 --category NAME [options]
ecai_cpp2 --help
```

One occupancy category per invocation. `--category` is required for any normal
execution; there is no default category and none is inferred.

## Flags

The interface accepts exactly 19 flags:

- `--category NAME` — occupancy category. Required. No default.
- `--simulations N` — simulations per replicate. Default `10000`.
- `--replicates R` — replicate count. Default `500`.
- `--percentile P` — percentile taken within each replicate. Default `0.96`.
- `--confidence LEVEL` — confidence level of the across-replicate interval. Default `0.95`.
- `--target-probability P` — target infection probability over the effective duration. Default `0.001`.
- `--duration-hours HOURS` — exposure duration in hours. Default `1`.
- `--base-seed SEED` — unsigned 64-bit base seed. Default `485`.
- `--community-rate RATE` — community infection rate. Default: the selected category value.
- `--occupants I0` — maximum occupants. Default: the selected category value.
- `--volume-m3 VOLUME` — space volume in m³. Default: the selected category value.
- `--mask-efficiency FRACTION` — mask efficiency. Default: the selected category value.
- `--mode MODE` — output mode, `ecai` or `probability`, case sensitive, no aliases. Default `ecai`.
- `--input-ecai ECAI` — fixed input ECAi for probability mode, L/s/person, finite and `>= 0`.
  Requires `--mode probability`. Default: the selected category's published Table 4 ECAi.
- `--output-directory PATH` — output directory; requests `summary.csv`.
- `--summary-file PATH` — literal summary CSV destination override.
- `--replicates-file PATH` — literal per-replicate CSV destination override.
- `--write-replicates` — boolean; requests per-replicate output.
- `--help` — print the help text and exit. Must be the only argument.

## Output mode

`ecai_cpp2` supports two output modes, selected with `--mode`:

- `ecai` (default): inverts Eqs. 1, 2, 6 and 7 for the ECAi that reproduces
  `--target-probability` over the effective duration. Identical to every invocation before
  `--mode` existed; omitting `--mode` or passing `--mode ecai` explicitly is byte-identical.
- `probability`: evaluates the forward Wells-Riley infection probability (Eq. 1 forward) for
  every sampled row at a fixed **input** ECAi, and reports the configurable within-replicate
  percentile (default `0.96`) of those per-simulation probabilities alongside the existing
  across-replicate statistics. See `docs/MODEL_SPECIFICATION.md` §3.1 for the forward equation,
  row rules, saturation, and comparator semantics, and `docs/SEED_POLICY.md` for what is and is
  not fingerprinted across modes.

`--input-ecai` is only accepted together with `--mode probability`; omitting it resolves to the
selected category's published Table 4 ECAi (`input_ecai_source=published_table4`). An explicit
value (`input_ecai_source=explicit`) marks the run a custom scenario. `--mode` alone — including
an explicit `--mode ecai` — never marks a run custom.

## Category selection

Category matching follows the frozen naming rule in
`docs/MODEL_SPECIFICATION.md` D7.

- The canonical identifier is the full publication name, for example
  `Exam room`, `Food and beverage facilities`, `Sorting, packing, light assembly`.
- Matching is case-insensitive, and runs of spaces, commas and hyphens are
  interchangeable delimiters. `Exam room`, `exam room`, `EXAM ROOM`,
  `exam-room` and `Exam,room` all select the same category.
- Ten explicitly enumerated legacy aliases are accepted: `Food`,
  `Transportation`, `Group`, `Patient`, `Resident`, `Place`, `Spectator`,
  `Common`, `Dwelling`, `Sorting`.
- There is no fuzzy, prefix or edit-distance matching. An input that does not
  match a canonical name or an enumerated alias after normalization is
  rejected, never guessed.
- Any reported category identifier is the canonical full publication name,
  never the alias that was typed.

## Defaults and domains

Defaults applied when a flag is omitted:

| Input | Default | Accepted domain |
| --- | --- | --- |
| `--simulations N` | 10000 | integer, `1 <= N <= 536870912` |
| `--replicates R` | 500 | integer, `1 <= R <= 10000001` |
| `--percentile P` | 0.96 | `0 <= P <= 1` |
| `--confidence LEVEL` | 0.95 | `0 < LEVEL < 1` |
| `--target-probability P` | 0.001 | `0 < P < 1` |
| `--duration-hours HOURS` | 1 | finite, `> 0` |
| `--base-seed SEED` | 485 | unsigned 64-bit integer, `0` to `18446744073709551615` |
| `--community-rate RATE` | category value | `0 <= RATE <= 1` |
| `--occupants I0` | category value | integer, `>= 1`, fits a 64-bit signed integer |
| `--volume-m3 VOLUME` | category value | finite, `> 0` |
| `--mask-efficiency FRACTION` | category value | `0 <= FRACTION <= 1` |
| `--mode MODE` | `ecai` | `ecai` or `probability`, case sensitive, no aliases |
| `--input-ecai ECAI` | category's published Table 4 ECAi (probability mode only) | finite, `>= 0`; only accepted with `--mode probability` |

The 20 general categories use community rate `0.01` and mask efficiency `0`.
The five categories in the Health Care group use community rate `0.03` and mask
efficiency `0.30`. Occupancy and volume remain category-specific values from the
25-record source table.

The maximum `N` is the frozen `sampling_plan::ResourcePolicy` ceiling
(`max_simulations`, 2^29 = 536870912). The maximum `R` is
`statistics::kMaxSupportedReplicateCount` = 10000001, which is the largest
replicate count whose integer degrees of freedom `R - 1` stay inside the
verified Student-t domain.

## Parsing rules

Numeric values are parsed strictly, locale-independently, over the whole
string, using `std::from_chars`.

- An empty value is rejected.
- Leading or trailing characters, including whitespace, are rejected
  (`0.5 `, ` 0.5`, `0.5x`, `1,5`, `0x1`).
- A leading `+` is rejected, because `from_chars` does not accept it.
- `NaN` and infinity are rejected in any spelling or case.
- Values that overflow or underflow the target type are rejected
  (`1e400`, `1e-400`, a base seed above the unsigned 64-bit range).
- Integer fields reject non-integer text, including `8.5`, `2.0` and `1e2`.
- A negative value for an unsigned field, such as `--base-seed -1`, is rejected.
- Path and category strings are preserved literally, including spaces; an empty
  string is rejected.

Structural rules:

- Unknown flags are rejected.
- Positional arguments are rejected.
- A repeated flag is rejected, including a repeated boolean flag such as
  `--write-replicates --write-replicates`.
- A flag whose value is missing is rejected, and so is a value that begins with
  `--`, since no accepted flag takes a value in that form.
- `--help` must be the only argument; combining it with anything else is
  rejected.
- Domain checks run once, after the whole argument vector is consumed, so a
  rejection never depends on flag order.

## Source defaults, effective values and custom scenarios

The selected occupancy-category record is retained alongside the resolved run,
so the published source defaults remain available next to the effective values;
the record is never overwritten by an override. An explicit CLI value wins
directly over the corresponding category default.

A run is recorded as a custom scenario when at least one explicit scientific
override is supplied:

- `--target-probability`
- `--duration-hours`
- `--percentile`
- `--community-rate`
- `--occupants`
- `--volume-m3`
- `--mask-efficiency`
- `--input-ecai`

Explicitly supplying a value equal to the default still counts as an override
for provenance purposes; for example `--duration-hours 1` marks the run custom.

These flags never make a run custom on their own:

- `--simulations`
- `--replicates`
- `--base-seed`
- `--confidence`
- `--mode`
- `--output-directory`, `--summary-file`, `--replicates-file`, `--write-replicates`

Per-field provenance is recorded separately for the four category inputs
(occupants, volume, mask efficiency, community rate) and for the run-level
scientific inputs, so a report can distinguish a published-category run from a
modified one.

## Duration and target semantics

The target probability always applies over the effective duration. It is
labelled an hourly target only when the effective duration is exactly 1.0 hour.
A duration of any other value is never labelled hourly. An explicitly supplied
duration of exactly 1 hour remains an hourly target and is simultaneously
recorded as a custom scenario.

## Statistics reported

Statistics are computed from the raw per-replicate percentile estimates.

- The raw across-replicate mean of the per-replicate percentile estimates, in
  L/s per person, is reported unrounded and unadjusted.
- The ceiling-to-5 value, `5 * ceil(raw_mean / 5)`, is a separately labelled
  comparator. It is never an input to the mean, sample SD, standard error,
  confidence interval, minimum or maximum.
- The confidence interval produced at `--confidence LEVEL` is an approximate
  randomized-LHS Monte Carlo Student-t interval for the mean of the replicate
  percentile estimator. It is not an interval for the population percentile of
  the ECAi distribution, not an interval for individual ECAi values, and not a
  characterization of model-input uncertainty.

## Console output

A valid run prints exactly one line to standard output, at round-trip
precision. In `ecai` mode (default):

```text
console_summary (not a CSV file) category="Patient room" raw_mean_l_s_person=41.606977504215571 ceiling_to_5_l_s_person=45
```

In `probability` mode the line instead reports the resolved input ECAi, the
requested percentile, the raw mean probability, the target, and the
comparator status, for example:

```text
console_summary (not a CSV file) mode=probability category="Patient room" input_ecai_l_s_person=35 percentile=0.96 raw_mean_probability=0.66694130367595861 target_probability=0.001 comparator_status=ABOVE_TARGET
```

Probability-mode comparator status is one of `AT_TARGET_WITHIN_10PCT`,
`BELOW_TARGET`, `ABOVE_TARGET` (comparator-eligible runs), or
`NOT_APPLICABLE_CUSTOM` (any custom scenario, including an explicit
`--input-ecai`) — a distinct vocabulary from the ECAi-mode comparator status
below.

The percentile is always printed as a number (never a `p96`-style label), so
any `--percentile` value renders unambiguously. These are console summaries,
not CSV output. Each statistic appears once. The category is the canonical
publication name. CSV output is opt-in and is written only after complete
validation of the document and destinations.

## Output files

Output is opt-in. An output directory writes `summary.csv`; with
`--write-replicates` it also writes `replicates.csv`. Explicit summary and
replicate paths are literal overrides for their respective destination; a
literal summary path takes precedence over the directory's `summary.csv`, and a
literal replicate path takes precedence over the directory's `replicates.csv`.
Replicate output requires
`--write-replicates`, destinations refuse overwrite, and the legacy `C++_legacy`
source tree is protected. Summary has one data row; replicate output has exactly
R rows and does not contain the N raw ECAi values.

The supported combinations are:

- no output flags: no files;
- `--output-directory DIR`: `DIR/summary.csv`;
- `--output-directory DIR --write-replicates`: `DIR/summary.csv` and
  `DIR/replicates.csv`;
- `--summary-file FILE`: summary only at the literal path;
- `--replicates-file FILE --write-replicates`: replicates only at the literal path;
- both literal file flags plus `--write-replicates`: both literal paths.

`--replicates-file` without `--write-replicates`, or `--write-replicates`
without either a replicate file or output directory, is rejected.

All document/schema and path validation precedes filesystem mutation. Each
payload is written, flushed and closed inside an exclusively reserved
same-parent temporary directory. Commit uses a no-replace hard link and retains
the temporary payload link until every requested destination is committed. On a
two-file failure, rollback removes a committed destination only when
`filesystem::equivalent` still proves that it names the invocation's retained
payload inode; a foreign replacement is preserved. Only empty owned temporary
and parent directories are removed. A filesystem that cannot create the
required hard link fails without creating that destination.

`WriteResult::cleanup_complete` reports whether best-effort removal of all
retained temporary payloads and directories completed. Once all destinations
are committed, cleanup failure does not turn the successful write into a
reported exception; committed paths are returned with
`cleanup_complete=false`. Normal successful writes report `true` and leave no
temporary directory.

Protected-tree and alias checks use resolved paths before mutation and are
repeated for each destination immediately before its no-replace commit. These
checks and inode-verified rollback protect the tested replacement races. The
C++17 `std::filesystem` pathname API does not provide directory-handle-relative
operations on every supported platform, so this is not a claim that every
hostile concurrent parent-directory or symlink swap is eliminated portably.

CSV uses RFC4180-style field quoting for commas, quotes, CR and LF, with LF
record endings. Floating-point fields use `max_digits10` round-trip precision.

### CSV schemas

The summary schema has exactly **89 columns**, in the order required by
`tests/test_output.cpp`, and exactly one row (`summary_row_count=1`). It includes
source category data, effective inputs, all run/output provenance, authoritative
replicate statistics, comparator fields, and aggregate counters. Bank plans and
cursors are intentionally not copied from the first replicate into the summary;
they remain per-replicate fields. `replicate_csv_row_count=R` states the separate
replicate-file row count. The c-drop source header is
`source_cdrop_gm_e4_per_m3`.

The replicate schema has exactly **19 columns** in frozen order and exactly R
rows. `replicate_index` is stable identity, not a row number synthesized by the
writer. Each row retains its planned bank lengths, final cursors, counters and
reconciliation flags.

For R=1, `dispersion_status=UNAVAILABLE_SINGLE_REPLICATE` and the six named
fields `raw_degrees_of_freedom`, `raw_sd_l_s_person`, `raw_se_l_s_person`,
`raw_t_critical`, `raw_ci_lower_l_s_person`, and `raw_ci_upper_l_s_person` are
empty. They are not written as zero or NaN.

For an unmodified published-category run, `comparator_eligible=true` and the
status is `EXACT_MATCH`, `NEAR_WITHIN_ONE_INCREMENT` (absolute difference 5),
or `FAIL_GT_ONE_INCREMENT`. Signed and absolute differences are both emitted.
For a custom scenario, `comparator_eligible=false`, status is
`NOT_APPLICABLE_CUSTOM`, and the published comparator and difference fields are
empty.

### Probability-mode CSV schema

A `--mode probability` run writes its own document, never the ECAi schema above; the two are
never mixed in one file. `schema_version=stage13-prob-v1`.

The probability summary schema has exactly **96 columns** and one row. Most source-category
columns (7–25) and explicit-provenance columns (78–94) are byte-identical in name and order to
the corresponding ECAi summary columns; the remaining columns are probability-mode specific,
including `mode`, `input_ecai_l_s_person`, `input_ecai_source` (`published_table4` or
`explicit`), `effective_eca_total_m3_h`, `zero_infector_policy`, `probability_policy`
(`forward_wells_riley_expm1`), `saturation_policy` (`probability_one_is_numerical_saturation`),
`comparator_target_probability`, `comparator_signed_difference`, `comparator_ratio`,
`comparator_status`, the `raw_*_probability` statistics, `all_replicates_saturated_count`, and
`all_replicates_probability_evaluations`. Columns that would be misleading in this mode are
deliberately absent: `negative_clamping_policy`, `all_replicates_negative_clamp_count`,
`all_replicates_ecai_evaluations`, `published_comparator_l_s_person`,
`comparator_abs_difference_l_s_person`, `ceiling_to_5_l_s_person` — probability mode never
negative-clamps and has no ceiling-to-5 comparator.

The probability replicate schema has exactly **20 columns** and R rows: the ECAi replicate schema
plus a `mode` column and `raw_percentile_probability`/`saturated_count`/
`probability_evaluations` in place of `raw_percentile_l_s_person`/`negative_clamp_count`/
`ecai_evaluations`.

For R=1, the same dispersion-unavailable convention applies: `dispersion_status=
UNAVAILABLE_SINGLE_REPLICATE` and the six named dispersion fields are empty, not zero or NaN.

No field ever contains a signed zero in either schema; a zero-infector row's probability is
always exact positive `+0.0`.

## Errors

Any parse, lookup, configuration or simulation failure exits with a nonzero
status and prints exactly two lines to standard error: one concise diagnostic
line and one `Try '<program> --help' for usage.` line. No partial scientific
output is produced.

`--help` exits 0, performs no category lookup and runs no simulation.

## Examples

Small explicit `N` and `R` keep these runs fast; they are illustrative, not the
publication configuration.

```text
ecai_cpp2 --help
ecai_cpp2 --category "Exam room" --simulations 4 --replicates 2
ecai_cpp2 --category Patient --simulations 4 --replicates 1
ecai_cpp2 --category exam-room --simulations 8 --replicates 2 --duration-hours 2
ecai_cpp2 --category "Patient room" --simulations 8 --replicates 2 --occupants 2 --volume-m3 50
```

The default invocation, `--category NAME` with no other flags, runs
10000 simulations across 500 replicates.

## Risk-averse infector policy

`--infector-policy community` is the default and preserves zero-infector community-mixture behavior. `--infector-policy risk-averse` conditions both ECAi and probability calculations on at least one infector, using a conditional PPF in CDF space and exactly N count-stream LHS values. It rejects `I0==0` and nonpositive community rates before output creation.

Risk-averse output uses `stage13-ra-v1` for ECAi and `stage13-prob-ra-v1` for probability. These schemas record policy, zero-mass and conditional-PPF provenance, enforce zero-infector count 0 and active count N, and use `NOT_APPLICABLE_RISK_AVERSE` with population comparator fields blank. C++ keeps community as the default in both modes; this differs from the legacy probability executable's historical require-infectors default.

## QER distribution selector

`--qer-distribution NAME` selects the QER (quanta emission rate) calculation family
(`docs/QER_DISTRIBUTION_IMPLEMENTATION_PLAN.md`). Accepted canonical values are:

- `jones` — full mechanistic Jones et al. calculation. This is the default when
  `--qer-distribution` is omitted.
- `jones-fitted` — the direct log10-normal fit to Jones et al. output.
- a full canonical Mikszewski profile name, e.g.
  `mikszewski-sars-cov-2-standing-speaking`. Bare `mikszewski` is rejected: no
  default profile is resolved for it.

Matching is case sensitive with no aliases. Any other value, including a
partial or malformed Mikszewski profile name, is rejected strictly at parse
time. Any explicit `--qer-distribution` value, including an explicit `jones`,
marks the run a custom scenario, because it changes the QER calculation
mechanism itself.

`--qer-distribution` is resolved to an immutable descriptor before sampling
begins. The selector controls the selected QER family/profile and the
direct-QER bank where applicable, propagates through the simulation loop, and
is recorded in versioned output provenance. Full Jones remains the default
when `--qer-distribution` is omitted.
