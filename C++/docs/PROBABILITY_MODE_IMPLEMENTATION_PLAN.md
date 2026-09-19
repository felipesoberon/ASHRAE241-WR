# Infection-probability-percentile mode — implementation plan

Status: plan only. Nothing in this document has been implemented, built, or run. Base commit for all line references: `7b93348` (`perf: optimize C++2 sampling pipeline`) — the C++2 tree, now `C++/`. Revision 2 (2026-09-13) incorporates the Codex review; the disposition of each review point is in §12.

## 1. Objective

Add a second output mode to `ecai_cpp2`. In the new mode the executable takes the selected category's published ASHRAE 241 ECAi value (Jones et al. 2025, Table 4, `OccupancyCategory::published_ecai_Lsp`, `include/category_data.h:54`) as a fixed input, runs the existing sampling pipeline unchanged, evaluates the **forward** Wells-Riley probability of infection for every simulation row instead of inverting for ECAi, and reports the configurable within-replicate percentile of those probabilities (default `0.96`) with the existing across-replicate statistics.

Primary acceptance criterion: **the existing ECAi mode is bit-for-bit unchanged** — same console line, same CSV bytes, same stderr diagnostics and exit codes, same seeds, same bank lengths and cursors, same schemas, same tests, same performance envelope — for every invocation that does not contain `--mode`, `--input-ecai` or `--help`. The `--help` text is the one deliberate, additive exception (§9.1).

### 1.1 Terminology (fixed for this document and the implementation)

- **Percentile** in this codebase means the `p`-quantile with `p ∈ [0, 1]` under the frozen Type 7 linear convention (`src/simulation.cpp:176-208`, rank `= p (N-1)`). `--percentile 0.96` is the 96th percentile = the 0.96-quantile. The word "percentile" is kept in all new names for consistency with the ECAi schema; no new name uses "quantile".
- **Replicate percentile** (`raw_percentile_probability`): the percentile taken *within* one replicate over its `N` per-simulation probabilities (zero-infector rows included as `+0.0`).
- **`raw_mean_probability`, `raw_sd_probability`, `raw_ci_*_probability`, `raw_min/max_probability`**: statistics *across* the `R` replicate percentiles, exactly as `raw_mean_l_s_person` etc. are across-replicate statistics of the replicate ECAi percentiles. The CI is a Student-t interval for the mean of the replicate-percentile estimator, not for the population percentile.
- **Per-simulation probability** `P`: the forward Wells-Riley probability for one susceptible in one sampled simulation row over duration `D` at the fixed input ECAi (§3.2). It is a property of the sampled row, not of the category.

## 2. What is reused unchanged

| Component | File | Reuse |
| --- | --- | --- |
| Category lookup, Table 4 published ECAi | `include/category_data.h:54`, `src/category_data.cpp` | `published_ecai_Lsp` is the default probability-mode input. No change. |
| Seed derivation, fingerprint | `include/seed_manager.h:58-70,101-126`, `src/seed_manager.cpp:90-146` | Untouched. The mode and the input ECAi are deliberately **not** fingerprint fields (§6). |
| Infector-count plan, banks, cursor, assignment | `include/sampling_plan.h`, `src/sampling_plan.cpp` | Untouched. Probability mode consumes exactly the same values in the same order. |
| Per-infector QER (Eq. 3-5) | `src/model.cpp:114-152` `compute_qer` | Reused as-is. |
| Percentile convention (Type 7 linear) | `src/simulation.cpp:176-208` `percentile_linear` | Reused on the probability vector. |
| Across-replicate statistics | `src/statistics.cpp:228` `summarize_replicates` | Reused; it accepts any finite non-negative estimates, so `[0, 1]` probabilities are in-domain. Its `ReplicateSummary` field names keep the `_l_s_person` suffix (no rename); in probability mode they hold probabilities and only the new document names them correctly. `adjusted_ceiling_to_5_l_s_person` is still computed (it is `5` for any mean in `(0, 1]`, `0` for a zero mean) but is **never reported** in probability mode. |
| Atomic CSV writer, path policy | `src/output.cpp:618-755` | Reused unchanged; only the document builder/validator gains a branch. |
| Strict numeric parsing, flag structure | `src/cli.cpp:55-157` | Reused for the two new flags. |

## 3. Numerical model

### 3.1 Existing inversion (unchanged, `src/model.cpp:155-206`)

```text
Q_t      = -log1p(-P_t)                                   (Eq. 1 inverted)
m        = (1 - mu)^2                                     (Eq. 2 mask factor)
term1    = PBR_s * D * m * QER_sum / Q_t                  [m^3/h]
term2    = VOL * (gamma + lambda)                         [m^3/h]
T_ECAi   = term1 - term2                                  [m^3/h]
ECAi     = max(T_ECAi, 0) / (3.6 * I0)                    [L/s/person]
```

### 3.2 New forward evaluation (probability mode)

Given the input `ECAi_in` [L/s/person] and one active row's sampled values:

```text
ECA_total = 3.6 * I0 * ECAi_in                            [m^3/h]   (Eq. 7 forward)
Lambda    = ECA_total + VOL * (gamma + lambda)            [m^3/h]   (total removal)
dose_num  = PBR_s * D * (1 - mu)^2 * QER_sum              [m^3 quanta / h]
Q         = dose_num / Lambda                             [quanta, dimensionless]
P_raw     = -expm1(-Q)                                    (Eq. 1 forward)
P         = (P_raw > 0) ? P_raw : +0.0                    (sign normalization, mirrors model.cpp:204)
```

Row rules:

- **Zero-infector row:** `P = +0.0` exactly, assigned before any bank value is requested (same short-circuit as ECAi mode, `src/simulation.cpp:283-292`). No equation runs.
- **No division by zero, no clamp:** `Lambda > 0` always (`ECAi_in >= 0`, `VOL > 0`, `gamma + lambda > 0`). `negative_clamp_count` is structurally zero in this mode and is not reported.
- **Range:** `Q >= 0` and finite ⇒ `P ∈ [0, 1]` **closed**. `P == 1.0` exactly is reachable: IEEE-754 double `-expm1(-Q)` rounds to `1.0` once `exp(-Q) < 2^-53`, i.e. for `Q` greater than roughly `36.7` (`ln 2^53`); the exact boundary is a property of the platform `expm1` and is not asserted. This is **numerical saturation**, not an error. It is recorded per row as `saturated = (P == 1.0)`, counted per replicate (`saturated_count`), and accepted by every validator (`0 <= P <= 1`). A saturated row is not an overflow: `Q` itself must still be finite (below).
- **`P == 0` for an active row** occurs only when `dose_num == 0`, i.e. `mu == 1` (mask factor zero) or an underflowed `QER_sum`; `compute_qer` guarantees `QER_i > 0` for in-domain inputs, so this is the only active-row zero source. `output_zero_count` therefore counts `zero_infector_count` plus these rows.
- **Precision:** `expm1` is used (not `1 - exp(-Q)`) to retain relative precision for the small `Q` (~1e-3) the default target implies.
- **Finiteness (every intermediate, in order):** `ECA_total`, `Lambda`, `dose_num`, `Q`, `P_raw` are each checked finite with a `require_finite_probability_intermediate(value, name)` helper (appended, §4.1) that throws `std::runtime_error("probability <name> is not finite (overflow or invalid intermediate)")`. Additionally `Q >= 0` is asserted (a negative `Q` is impossible after validation and would indicate a defect, so it is a `runtime_error`, not silently clamped). `ECA_total` can overflow for a pathological `--input-ecai` (e.g. `1e308`); that is rejected by the `ECA_total` check with a diagnostic, never propagated. `simulation.cpp` re-checks `std::isfinite(P)` after the call exactly as it does for ECAi (`:388-391`).

Dimensional check: `(m^3/h)(h)(-)(quanta/h)/(m^3/h) = quanta`. Consistent with `docs/MODEL_SPECIFICATION.md` §10.

### 3.3 Exact algebraic relations used as tests

1. **Round trip.** For an active row with `T_ECAi > 0` and `P_t < 1` (always true, `P_t ∈ (0, 1)`), evaluating §3.2 with `ECAi_in = ECAi_raw` (unclamped) returns `P = P_t` to relative `1e-12`. Applies only to non-saturated rows; the default `P_t = 0.001` never saturates.
2. **Indicator identity.** For any `X >= 0`, row-wise `P_i(X) > P_t  ⇔  ECAi_raw_i > X`. Hence `#{P_i(X) > P_t} = #{ECAi_raw_i > X}` for the same sampled rows. Numerically the two sides may disagree by one row when `|ECAi_raw_i - X| <= 1e-9 * max(1, |X|)`; the test excludes rows inside that band and requires exact equality on the rest.
3. **Order-statistic identity.** Choose `N` with `p (N-1)` an exact integer in double arithmetic (verify at runtime that `percentile_linear`'s `frac == 0`; `N = 26`, `p = 0.96` gives rank `24.0` exactly) and a case whose ECAi-mode percentile `X = sorted_ECAi[24] > 0` is unique (no tie). Then in probability mode with `ECAi_in = X`, `sorted_P[24] = P_t` to relative `1e-10`. This follows from (1) and (2) — `24` rows have `P < P_t`, `1` row has `P > P_t`, the row at `X` has `P = P_t` — and **not** from rank preservation: because `PBR_s`, `gamma`, `lambda` differ per row, the sort order of `P_i(X)` is not in general the sort order of `ECAi_raw_i`. The test must not assert full order preservation.
4. **Monotonicity.** `P` is strictly decreasing in `ECAi_in` for any row with `dose_num > 0` and `P < 1` (a saturated row is flat at `1.0` until `ECAi_in` is large enough to leave saturation).

These hold only because both modes consume identical samples (§6).

## 4. Public API changes

All additions are new symbols or new fields with defaults that reproduce current behaviour. No existing signature changes. `OutputMode::ecai` is enumerator `0`, so a value-initialized `OutputMode{}` and every defaulted struct member are ECAi mode; a `static_assert(static_cast<int>(OutputMode::ecai) == 0)` pins this.

### 4.1 `include/model.h` (append after line 146) / `src/model.cpp` (append after line 206)

The forward evaluation is **decoupled from `target_probability`**: it does not take, validate or read the target. The target is only a comparator at the output stage (§5.3).

```cpp
// Seven physical inputs (ECAiInputs minus target_probability) plus the fixed ECAi.
struct ProbabilityInputs {
    double susceptible_pbr_m3h = 0.0;
    double gamma_per_h = 0.0;
    double lambda_per_h = 0.0;
    double volume_m3 = 0.0;
    long long occupants = 0;
    double mask_efficiency = 0.0;
    double duration_h = 0.0;
    double input_ecai_l_s_person = 0.0;   // finite, >= 0
};
struct ProbabilityResult {
    double eca_total_m3_per_h = 0.0, removal_m3_per_h = 0.0, mask_factor = 0.0,
           qer_sum_per_h = 0.0, dose_numerator = 0.0, quanta_inhaled = 0.0,
           probability = 0.0;             // in [0, 1] closed, never -0.0
    bool saturated = false;               // probability == 1.0
};
ProbabilityResult compute_probability_from_qer_sum(const ProbabilityInputs&, double qer_sum_per_h);
ProbabilityResult compute_probability(const ProbabilityInputs&, const std::vector<double>& qer_per_h);
```

Implementation constraints:

- Everything is appended after `src/model.cpp:206`, including a **second anonymous-namespace block** holding `validate_probability_inputs` (the seven checks of `validate_ecai_inputs`, `:95-104`, with `"ProbabilityInputs."` message prefix, plus `input_ecai_l_s_person` finite and `>= 0`) and `require_finite_probability_intermediate`. The existing anonymous namespace (`:10-112`) is not edited; `require_qer_value` (`:70-77`) is reused by name from the appended code (same translation unit).
- `compute_probability` sums in supplied order with the same overflow check as `compute_ecai` (`:215-222`).
- `compute_ecai*` bodies are not touched.

### 4.2 `include/simulation.h`

- New `enum class OutputMode { ecai = 0, probability = 1 };` (before line 54).
- `EffectiveScenario` (lines 54-75): add `OutputMode mode = OutputMode::ecai; double input_ecai_l_s_person = 0.0; bool input_ecai_overridden = false;`. **`custom_scenario` (line 74) keeps its current meaning** — "one of the four category scientific inputs was overridden" — and is *not* ORed with `input_ecai_overridden`; that keeps the ECAi `build_document` provenance check (`src/output.cpp:396-401`) valid without edits. An explicit `--input-ecai` marks the run custom at the CLI level (`CliOptions::custom_scenario`, like `--target-probability`), §4.4.
- `EffectiveOverrides` (lines 81-86): add `OutputMode mode = OutputMode::ecai; std::optional<double> input_ecai_l_s_person;`. Defaults keep every existing call site identical.
- `ReplicateCounters` (lines 143-154): add `std::size_t probability_evaluations = 0;` (one per active row in probability mode) and `std::size_t saturated_count = 0;` (active rows with `P == 1.0`). Both are zero in ECAi mode; `negative_clamp_count` and `ecai_evaluations` are zero in probability mode.
- `RowAudit` (lines 158-173): add `double probability = 0.0; bool saturated = false;`.
- `ReplicateResult` (lines 190-211): add `OutputMode mode = OutputMode::ecai; std::vector<double> probability; double percentile_probability = 0.0;`. In ECAi mode these stay `ecai`/empty/zero; in probability mode `ecai_l_s_person` stays empty and `percentile_l_s_person` stays `0.0`. Existing fields are never reinterpreted.
- `run_replicates` (lines 313-319): branch on the **resolved** mode carried back in `result.mode` (set from `scenario.mode` in `run_replicate_from_banks`), not on `overrides.mode`: `raw_percentiles.push_back(result.mode == OutputMode::probability ? result.percentile_probability : result.percentile_l_s_person)`. All replicates in a batch must report the same mode (checked; a disagreement is a `runtime_error`).

### 4.3 `src/simulation.cpp`

- `validate_scenario` (lines 97-111): append: `mode` must be one of the two enumerators (same defensive switch as `validate_evaluation_order`, `:51-59`); if `mode == probability`, `input_ecai_l_s_person` must be finite and `>= 0`; if `mode == ecai`, `input_ecai_overridden` must be `false` and `input_ecai_l_s_person == 0.0` (an ECAi scenario carrying an input ECAi is a defect).
- `resolve_scenario` (lines 210-262): after line 230 set `scenario.mode = overrides.mode`; if probability mode, `scenario.input_ecai_l_s_person = overrides.input_ecai_l_s_person ? *… : category.published_ecai_Lsp` and `input_ecai_overridden = has_value()`; in ECAi mode both stay at their defaults regardless of `overrides.input_ecai_l_s_person` (which the CLI has already rejected). Reject a category record whose `published_ecai_Lsp` is not finite or is negative (defensive; none expected). Fingerprint agreement checks (lines 242-250) are unchanged — the new fields are not in `EffectiveConfig`.
- `validate_plan_context` (lines 66-94): unchanged.
- `run_replicate_from_banks`: `result.mode = scenario.mode` after line 279. The ECAi path is preserved by **dispatch, not by editing**:
  - Step 1 (line 287): `std::vector<double>& outputs = (scenario.mode == OutputMode::probability) ? result.probability : result.ecai_l_s_person; outputs.assign(plan.n, 0.0);` — the only line of the ECAi flow that changes text.
  - Step 3 (lines 349-400): the existing loop body is wrapped verbatim in `if (scenario.mode == OutputMode::ecai) { … } else { evaluate_probability_schedule(…); }` where the probability loop lives in a new anonymous-namespace function above `run_replicate_from_banks`. That loop is identical to the ECAi loop up to the `compute_qer` accumulation (lines 372-377), then builds `ProbabilityInputs` (no target), calls `compute_probability`, increments `probability_evaluations` (and `saturated_count` when `saturated`), re-checks `std::isfinite`, writes `result.probability[row]`, and fills `RowAudit` (`qer_sum_per_h`, `probability`, `saturated`; `ecai_*` and `clamped` left at defaults).
  - Zero count (lines 402-406) and Step 4 (line 422) iterate/percentile over `outputs`; Step 4 stores into `percentile_probability` or `percentile_l_s_person` by mode.
  - `git diff -w 7b93348 -- src/simulation.cpp` must show only: the three `outputs` lines, the wrapping `if/else`, the `result.mode` assignment, and inserted code. The ECAi loop body is byte-identical modulo indentation.

### 4.4 `include/cli.h` / `src/cli.cpp`

- Defaults (`cli.h:42-48`): add `inline constexpr simulation::OutputMode kDefaultMode = simulation::OutputMode::ecai;` (no default constant for `--input-ecai`: it resolves from the category).
- `CliOptions` (`cli.h:61-120`): add `simulation::OutputMode mode = kDefaultMode; bool mode_explicit = false; std::optional<double> input_ecai_l_s_person;`.
- `accepted_flags()` (`cli.cpp:236-245`): append `"--mode"`, `"--input-ecai"` → **19 flags**.
- `help_text()` (`cli.cpp:247-319`): insert one new block "Output mode:" between the "Scientific inputs" and "Output controls" blocks, documenting `--mode MODE` (`ecai` | `probability`, default `ecai`, never custom) and `--input-ecai ECAI` (probability mode only; L/s/person, finite, `>= 0`; default: the selected category's published Table 4 ECAi; an explicit value marks the run custom), plus one caveat sentence (§10). **No existing line of the help text is edited, removed or reordered**: the pre-change help text, split on `\n`, must be an in-order subsequence of the post-change text (T8b, §9.3). The help text does not state a flag count, so none needs changing.
- `parse_cli` (`cli.cpp:356-398`): `--mode` accepts exactly `ecai` or `probability` (case-sensitive, no aliases), sets `mode_explicit`; `--input-ecai` uses `parse_double_strict`.
- `validate_domains` (`cli.cpp:162-204`): `--input-ecai` requires `>= 0`; `--input-ecai` without `--mode probability` is rejected (`"flag --input-ecai requires --mode probability"`); `--mode probability` with `--target-probability` is allowed (it is the comparator, §5.3).
- `custom_scenario` (`cli.cpp:404-408` and `:427-431`): also OR `input_ecai_l_s_person.has_value()`. `--mode` alone never marks custom.
- `validate_scalar_provenance` (`cli.cpp:209-232`): add `if (!mode_explicit && mode != kDefaultMode) reject(...)`.
- `resolve_run` (`cli.cpp:414-480`): copy `mode`/`input_ecai` into `run.overrides` after line 457. `run.config` assembly (lines 459-473) is **not** modified, so the fingerprint string and every derived seed are identical across modes.

### 4.5 `src/main.cpp` (lines 59-66)

Branch on `run.scenario.mode` (the resolved scenario; `run.overrides.mode` and every `batch.replicates[i].mode` must agree, which `build_document` checks first). ECAi branch: existing statements verbatim. Probability branch prints one line:

```text
console_summary (not a CSV file) mode=probability category="Patient room" input_ecai_l_s_person=35 percentile=0.96 raw_mean_probability=<double> target_probability=0.001 comparator_status=<status>
```

followed by the same `summary_file=`/`replicates_file=` suffixes. The percentile is printed as a number, never formatted into a label (`p96`), so any `--percentile` value renders unambiguously. Doubles use `max_digits10` as today.

## 5. Output / schema changes

Rule: the ECAi summary header (`src/output.cpp:155-188`, 89 columns) and replicate header (`:190-197`, 19 columns), `kOutputSchemaVersion = "stage13-v1"`, and the ECAi `build_document` body (`:381-616`) are **not edited**. Probability mode gets its own document.

### 5.1 Code changes

- `include/output.h:17`: add `kProbabilitySchemaVersion = "stage13-prob-v1"`.
- `src/output.cpp`: add `probability_summary_header()`, `probability_replicate_header()` (§5.2) and `build_probability_document(run, batch, options)`. `build_document` gains a **leading guard block** before its first statement (line 385): `if (run.scenario.mode == probability) return build_probability_document(...);` followed by `require(run.overrides.mode == ecai && every replicate.mode == ecai && replicate.probability.empty() && replicate.percentile_probability == 0.0 && counters.probability_evaluations == 0 && counters.saturated_count == 0, "ECAi document requested for a non-ECAi batch")`. Nothing after the guard changes.
- `build_probability_document` mirrors the ECAi builder's provenance checks and additionally requires: `run.overrides.mode == run.scenario.mode == probability`; every `replicate.mode == probability`; `ecai_l_s_person.empty()`, `percentile_l_s_person == 0.0`, `negative_clamp_count == 0`, `ecai_evaluations == 0`; `probability_evaluations == active_count`; `saturated_count <= active_count`; every probability finite, `>= 0`, `<= 1`, and `+0.0` (not `-0.0`) on zero-infector rows; `saturated_count == #{P == 1.0}`; `output_zero_count == #{P == 0}`; recomputed `percentile_linear` bit-equal to `percentile_probability`; `summarize_replicates` bit-equal to `batch.summary`.
- `validate_stage13_document` (`:199-220`): dispatch on `summary_columns[0] == "schema_version"` and `summary_rows[0][0]` (`stage13-v1` → existing body verbatim; `stage13-prob-v1` → new validator using the §5.2 indices).
- `require_same_run` (`:270-311`): add equality checks for `overrides.mode`, `overrides.input_ecai_l_s_person` (`same_optional`), `scenario.mode`, `scenario.input_ecai_l_s_person` (`same_double`), `scenario.input_ecai_overridden`. Pure equality checks that hold trivially at defaults.

### 5.2 Probability schemas (complete, frozen at implementation)

Summary (`stage13-prob-v1`), one row, **96 columns**, 1-based positions:

```text
 1 schema_version                      2 model_version
 3 sampling_design_version             4 mode                      (= "probability")
 5 category                            6 group
 7 is_health_care_group
 8-26  the 19 source_* columns, byte-identical names and order to ECAi columns 7-25
27 effective_i0                       28 effective_volume_m3
29 effective_community_rate           30 effective_mask_efficiency
31 input_ecai_l_s_person              32 input_ecai_source         (published_table4 | explicit)
33 effective_eca_total_m3_h           (= 3.6 * effective_i0 * input_ecai)
34 target_probability                 35 duration_hours
36 target_probability_over_duration   (hourly_target | target_over_selected_duration, as ECAi)
37 simulations_n                      38 replicates_r
39 percentile                         40 confidence_level
41 base_seed                          42 seed_policy_version
43 zero_infector_policy               (= zero_infector_short_circuit)
44 probability_policy                 (= forward_wells_riley_expm1)
45 saturation_policy                  (= probability_one_is_numerical_saturation)
46 custom_scenario                    47 comparator_eligible
48 comparator_target_probability      49 comparator_signed_difference
50 comparator_ratio                   51 comparator_status
52 dispersion_status
53 raw_mean_probability               54 raw_degrees_of_freedom
55 raw_sd_probability                 56 raw_se_probability
57 raw_t_critical                     58 raw_ci_lower_probability
59 raw_ci_upper_probability           60 raw_min_probability
61 raw_max_probability
62 summary_row_count                  63 replicate_csv_row_count
64 all_replicates_n                   65 all_replicates_zero_infector_count
66 all_replicates_active_count        67 all_replicates_k_total
68 all_replicates_output_zero_count   69 all_replicates_saturated_count
70 all_replicates_qer_evaluations     71 all_replicates_probability_evaluations
72 all_replicates_susceptible_values_consumed
73 all_replicates_infector_values_consumed
74 counters_reconcile                 75 all_cursors_exhausted
76 all_replicates_finite              77 banks_reconcile
78-94 the 17 explicit_* columns, byte-identical names and order to ECAi columns 73-89
95 explicit_mode                      96 explicit_input_ecai
```

Deliberately absent (they would be misleading in this mode): `negative_clamping_policy`, `all_replicates_negative_clamp_count`, `all_replicates_ecai_evaluations`, `published_comparator_l_s_person`, `comparator_abs_difference_l_s_person`, `ceiling_to_5_l_s_person`. Columns 48-51 are empty strings when `comparator_eligible` is `false` (status `NOT_APPLICABLE_CUSTOM`), as in the ECAi schema.

Replicate (`stage13-prob-v1`), `R` rows, **20 columns**:

```text
 1 replicate_index          2 mode                       3 requested_percentile
 4 raw_percentile_probability                            5 n
 6 zero_infector_count      7 active_count               8 k_total
 9 saturated_count         10 output_zero_count         11 qer_evaluations
12 probability_evaluations 13 susceptible_values_consumed
14 infector_values_consumed                             15 planned_bank_lengths
16 final_cursors           17 counters_reconcile        18 cursors_exhausted
19 finite_values           20 replicate_csv_row_count
```

Validator metadata indices (0-based): `replicates_r = 37`, `summary_row_count = 61`, `replicate_csv_row_count = 62`, replicate-row count `= 19`. The width constants `96`/`20` are pinned in `build_probability_document` and in T9 exactly as `89`/`19` are today (`:596,600`).

### 5.3 Target comparator (probability mode)

- **Quantity compared:** `raw_mean_probability` (column 53), the across-replicate mean of the replicate percentiles, against `target_probability` (column 34). No ceiling, no rounding, no 5 L/s/person grid.
- **Eligibility:** `comparator_eligible = !run.custom_scenario`. Any of `--target-probability`, `--duration-hours`, `--percentile`, `--community-rate`, `--occupants`, `--volume-m3`, `--mask-efficiency`, `--input-ecai` makes the run custom and suppresses the comparator. Consequently an eligible comparison always uses the frozen target `0.001`, `D = 1 h` (`hourly_target`), `p = 0.96`, the category's own `I0`/`VOL`/`mu`/rate and the published Table 4 ECAi.
- **Duration:** the target applies over the effective `D`; because an explicit duration is custom, every eligible comparison is hourly. A custom-duration run still reports column 36 so the reader can see which target the printed probabilities relate to.
- **Statistic and boundaries (closed, evaluated on doubles):** `comparator_signed_difference = mean - target`; `comparator_ratio = mean / target` (`target ∈ (0, 1)`, so never a division by zero; `mean == 0` gives ratio `0`). `comparator_status = AT_TARGET_WITHIN_10PCT` iff `0.9 <= ratio <= 1.1`; `BELOW_TARGET` iff `ratio < 0.9`; `ABOVE_TARGET` iff `ratio > 1.1`. The 10 % band is a self-consistency tolerance, not a standard's acceptance criterion.
- **Zero-infector behaviour:** zero-infector rows enter each replicate's percentile as `P = +0.0`. When the zero-infector fraction of a replicate is `>= p` the replicate percentile is exactly `0` (e.g. Patient ≈ 91 % zero-infector at `p = 0.96` sits just off that boundary and its percentile is governed by the sparse active tail). A mean of `0` yields `BELOW_TARGET`; that is the correct reading of the sample, not a defect. `all_replicates_zero_infector_count` and the per-replicate `raw_percentile_probability` values let the reader see it.
- **Interpretation:** the comparator asks "at the published ECAi, does the C++ sample's 96th-percentile per-simulation infection probability reproduce the 0.001 target?" It is a self-consistency diagnostic of the C++ sample; it is not the probability that the standard's target is met and does not validate either implementation (§10).

## 6. Seed and sampling semantics

- The seed tuple, `canonical_config_string` (`src/seed_manager.cpp:90-120`), `kDesignVersion = 1` and `kSeedPolicyVersion` are unchanged. Mode and input ECAi are **not** serialized. Consequence: for identical `(category, overrides, N, R, percentile, confidence, target, duration, base seed)`, both modes derive identical stream seeds, identical `infector_counts`, identical banks, identical consumption order. This is intentional and is the basis of §3.3 and the cross-mode tests.
- Design version is **not** incremented: stochastic stream construction, assignment and consumption are unchanged; only the deterministic per-row map differs. Record this decision in `docs/SEED_POLICY.md` (after the "Verified uniqueness scope" section, line 93) and `docs/SAMPLING_DESIGN.md` (line 115) during implementation.
- **`target_probability` stays in the fingerprint** (`EffectiveConfig`, `seed_manager.h:64`). Removing or conditionally omitting it would change every ECAi-mode seed and is therefore excluded. Implications to document in `docs/SEED_POLICY.md` and `docs/CLI_REFERENCE.md`:
  1. In probability mode the target is only a comparator, yet `--target-probability 0.002` derives different seeds than `0.001`: a target sweep in probability mode is **not** a same-sample sweep. Users wanting to test several thresholds against one sample set should run once at the default target and compare `raw_percentile_probability` offline.
  2. An explicit target also makes the run custom (comparator suppressed), so the only comparator-eligible target is the frozen default.
  3. A probability run and an ECAi run with the same target/duration/percentile/N/R/seed share samples exactly (§3.3); a probability run at a different target shares samples with the ECAi run at that same target, not with the default ECAi run.
  4. `--input-ecai` does not affect seeds, so sweeping input ECAi over the same seed reuses one sample set (used by §3.3-2, -3, -4).
  5. `--mode` does not affect seeds; `--mode ecai` explicit is byte-identical to omitting it.

## 7. Published-ECAi resolution

1. `--input-ecai X` explicit → `X` (finite, `>= 0`), `input_ecai_source = explicit`, run is custom.
2. Otherwise → `category.published_ecai_Lsp` (Table 4 adjusted value, multiples of 5), `input_ecai_source = published_table4`.
3. `resolve_scenario` rejects a category record whose `published_ecai_Lsp` is not finite or negative (defensive; none expected).
4. `ECA_total` uses the **effective** `I0` (`scenario.occupants`). An `--occupants` override therefore rescales total ECA; the run is custom and the comparator is suppressed.

## 8. Implementation steps and commit structure

Two commits on a feature branch (merging/pushing is a separate decision):

**Commit A — `test: freeze ECAi-mode baseline before probability mode`** (no production code):

0. Build at `7b93348` (Release, §9.4 configuration), run the §9.1 command set, record the artefacts into `docs/ECAI_MODE_BASELINE.sha256` (new file, same header style as `docs/ISOLATION_BASELINE.sha256`: purpose, base HEAD, UTC timestamp, exact generating command). Add the fingerprint known-answer test (§9.2, additive lines in `tests/test_seed_manager.cpp`) and the ECAi-mode invariant checks (§9.2, additive lines in `tests/test_simulation.cpp` / `tests/test_output.cpp`, written only against symbols that already exist). This commit must pass `ctest` against the unchanged code — that is what proves the pins are correct.

**Commit B — `feat: add probability output mode`** (everything else):

1. `model.h/.cpp`: §4.1. Unit tests T1-T4.
2. `simulation.h/.cpp`: §4.2-4.3. Tests T5-T7, T10.
3. `cli.h/.cpp`: §4.4. Update the 17-flag pin in `tests/test_cli.cpp:225-232` (planned-flag vector + message) to 19 — the **only** modification of an existing assertion in the whole change. The `--help` byte comparison at `test_cli.cpp:716` compares against `help_text()` and needs no edit. Tests T8, T8b.
4. `output.h/.cpp`: §5. Tests T9.
5. `main.cpp`: §4.5. Test T11.
6. `CMakeLists.txt`: register `tests/test_probability_mode.cpp` as its own executable/`add_test` in the existing pattern (lines 196-228), linking the same sources as `ecai_cpp2_output_tests`.
7. Docs: `README.md:42` (17 → 19, mode paragraph), `docs/CLI_REFERENCE.md:24` (19 flags, new schema section, §6 implications), `docs/MODEL_SPECIFICATION.md` (new §3 subsection "Forward probability evaluation used by probability mode", including saturation), `docs/SEED_POLICY.md`, `docs/SAMPLING_DESIGN.md`, `docs/VALIDATION.md` (probability-mode bounded results and the §9.4 benchmark record; clearly not a publication pass).
8. Run every §9 gate. Only then commit B.

Test policy: existing assertions may not be modified or deleted except step 3. **Additive** tests in existing test files are permitted (and required by commit A); each addition must be a new `check(...)` block that does not reorder or reuse existing fixtures in a way that changes their observed values.

## 9. Preservation gates for ECAi mode (all must pass)

### 9.1 Baseline artefacts

Run before (commit A) and after commit B from the `C++` directory as `./build/ecai_cpp2 …` so `argv[0]` (echoed in the stderr "Try '<program> --help'" line) is stable. `T` is a fresh empty directory per case, deleted between cases. For each case record separately: (i) exit code, (ii) SHA-256 of stdout after normalizing `summary_file="…"`/`replicates_file="…"` to `<PATH>`, (iii) SHA-256 of stderr, (iv) the sorted list of files present under `T` afterwards, (v) SHA-256 of each file. Byte identity is required for every item of cases B-M; case A is handled separately.

```text
A  --help                                    RECORDED, NOT hash-compared: help text changes by
                                             design (new flags). Gate = exit 0, stdout == help_text(),
                                             no "console_summary", and T8b subsequence check.
B  --category "Exam room" --simulations 17 --replicates 3 --output-directory T --write-replicates
C  --category "Patient room" --simulations 10000 --replicates 10 --output-directory T --write-replicates
D  --category Museum --simulations 10000 --replicates 5 --output-directory T --write-replicates
E  --category exam-room --simulations 8 --replicates 2 --duration-hours 2 --output-directory T --write-replicates
F  --category Patient --simulations 4 --replicates 1 --output-directory T          (R=1: empty dispersion fields)
G  --category "Exam room" --simulations 6 --replicates 2 --base-seed 0 --output-directory T --write-replicates
H  --category "Exam room" --bogus-flag                                                (stderr, exit 1, T absent)
I  --category "Exam room" --simulations 4 --replicates 1 --replicates-file T/r.csv   (stderr "requires --write-replicates", exit 1, T absent)
J  --category Nope --simulations 4 --replicates 1 --output-directory T               (unknown category, exit 1, T absent)
K  case F run twice into the same T                                                  (2nd: "refusing to overwrite", exit 1, 1st file bytes unchanged)
L  --category "Exam room" --simulations 4 --replicates 1 --volume-m3 50 --output-directory T   (custom: comparator NOT_APPLICABLE_CUSTOM)
M  --category "Exam room" --simulations 4 --replicates 1 --mode probability          (baseline: "unknown flag", exit 1 — RECORDED only;
                                                                                     after commit B this is a valid run and is not compared)
```

Additional assertions on every CSV produced by B-G, K, L (before and after): no field matches `(^|,)-0(\.0*)?(,|$)` (no signed zero ever reaches an output); header lines byte-identical to baseline. C and D must additionally still match `data/reference/publication_cases.csv` exactly (existing focused test). Failing cases (H, I, J, second K) must leave `T` non-existent or unchanged: the writer creates parents only after validation (`src/output.cpp:691-702`) and that must remain observable.

### 9.2 Existing test-suite gates

- `ctest` full suite passes with zero modified expectations except the flag-count pin (§8 step 3).
- `tests/test_output.cpp` expected 89/19 headers unchanged; `tests/test_end_to_end.cpp` regression fixtures (`ECAi 1.01560284592605` at `:265`, counters, byte-identical same-seed repeat at `:396-399`) unchanged.
- Commit A additive pins: (a) `tests/test_seed_manager.cpp`: known-answer test for `canonical_config_string` and `config_fingerprint` of the Exam-room default configuration (`10000 × 500`, target `0.001`, `D = 1`, `p = 0.96`, confidence `0.95`), values observed at `7b93348`; (b) `tests/test_simulation.cpp`: after a default `run_replicate` (Exam room, `N = 17`, seed `485`, index `0`), pin `infector_counts`, `planned_bank_lengths`, `final_cursors`, `bank_generated`, `counters`, and the bit pattern of `percentile_l_s_person` and of every `ecai_l_s_person[i]`; (c) `tests/test_output.cpp`: pin the SHA-256 (or the full text) of `serialize_summary_csv`/`serialize_replicates_csv` for the same batch. After commit B these pins must not move; they are the in-process complement of §9.1.

### 9.3 Focused new tests (`tests/test_probability_mode.cpp`)

- **T1** round trip (§3.3-1), rel `1e-12`, over all 25 categories using a fixed hand-built `InfectorParameters`; for each, also assert `!saturated`.
- **T2** zero dose: empty QER vector → `P == 0.0 && !signbit`; `mu = 1` with positive QER → `P == 0.0 && !signbit`.
- **T2b** saturation: inputs giving `Q = 1e3` → `P == 1.0 && saturated`; `Q = 50` → `P == 1.0`; `Q = 20` → `P < 1.0 && !saturated`. Inputs giving overflowing `ECA_total` (`input_ecai = 1e308`, `I0 = 1000`) or overflowing `dose_num` → `std::runtime_error` naming the intermediate; `input_ecai` negative / NaN / Inf → `std::invalid_argument`.
- **T3** monotonicity in `input_ecai` (§3.3-4) on non-saturated rows; **T4** domain rejections (each of the seven physical fields, with `ProbabilityInputs.` message prefix; unchanged `ECAiInputs` rejections still fire, and `ProbabilityInputs` has no target field to reject).
- **T5** same-seed cross-mode identity (mode isolation): for Exam room `N = 17, R = 3, seed 485` and Museum `N = 200, R = 2`, resolve both CLI invocations (with and without `--mode probability`) and assert (a) `run.config` equal field-by-field (bitwise for doubles) and `config_fingerprint` equal; (b) `build_sampling_plan` → equal `counts`; (c) `build_parameter_banks` → for all 11 banks equal `seed`, `generated`, and bitwise-equal `values`; (d) `assign_replicate` → equal consumption ids and `final_cursors`; (e) `run_replicate` with `collect_audit = true` in both modes → equal `infector_counts`, `planned_bank_lengths`, `final_cursors`, `bank_generated`, `cursors_exhausted`, `audit.consumption_ids`, `audit.evaluation_order`, per-row `susceptible_ids`/`infector_ids`, per-row `qer_sum_per_h` bitwise, and all counters except `negative_clamp_count`, `ecai_evaluations`, `probability_evaluations`, `saturated_count`. (Seeds cannot differ by construction because neither mode nor input ECAi is fingerprinted; the test is a regression guard against a future fingerprint edit.)
- **T6** indicator identity (§3.3-2) on audit rows for three `X` values with the tolerance band; **T7** order-statistic identity (§3.3-3) at `N = 26`, Museum, with the runtime `frac == 0` and no-tie preconditions asserted.
- **T8** CLI: defaults (`mode == ecai`, `!mode_explicit`, not custom); `--mode probability` parses; `--mode Probability`, `--mode foo`, `--mode` twice, `--input-ecai` without probability mode, `--input-ecai -1`, `--input-ecai nan` rejected with single-line diagnostics; `--mode ecai` explicit does not mark custom; `--mode probability --input-ecai 30` is custom; `accepted_flags().size() == 19`.
- **T8b** help text: the `7b93348` help text (embedded verbatim in the test) split on `\n` is an in-order subsequence of `help_text()`; every new `--` token is an accepted flag (existing bidirectional check covers it).
- **T9** output: probability document has width 96/20 and one row / `R` rows; column names exactly as §5.2; ECAi `build_document` on a probability batch is refused and `build_probability_document` on an ECAi batch is refused; probability-mode same-seed repeat is byte-identical; destination refusal, overwrite refusal and protected-tree checks behave identically in both modes; no `-0` field; `raw_min_probability >= 0`, `raw_max_probability <= 1`; comparator columns empty and status `NOT_APPLICABLE_CUSTOM` when custom; `AT_TARGET_WITHIN_10PCT`/`BELOW_TARGET`/`ABOVE_TARGET` boundaries exercised via a hand-built `ReplicateBatch` with means `0.0009`, `0.0011`, `0.00089`, `0.00111`, `0`.
- **T10** mode isolation of results: `run_replicate` in ECAi mode leaves `mode == ecai`, `probability` empty, `percentile_probability == 0`, `probability_evaluations == 0`, `saturated_count == 0`; probability mode leaves `ecai_l_s_person` empty, `percentile_l_s_person == 0`, `negative_clamp_count == 0`, `ecai_evaluations == 0`; `run_replicates` with mixed-mode results is impossible through the public API (documented, not tested).
- **T11** executable matrix (same `run_command` harness as `test_cli.cpp:105-130`): `--mode probability --category "Exam room" --simulations 17 --replicates 3 --output-directory T --write-replicates` → exit 0, console line starts with `console_summary (not a CSV file) mode=probability`, files exist with the 96/20 headers; the ECAi console line for the same args without `--mode` is byte-identical to baseline B.

### 9.4 Benchmark protocol (reproducible)

- Build: `cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release` with the default WSL toolchain (record `c++ --version` and the exact CMake cache flags in the result). The same build directory, compiler and flags for before/after; before = `7b93348`, after = commit B.
- Cases: (i) Exam room `N = 10000, R = 10`; (ii) Patient `N = 10000, R = 10`; (iii) Patient `N = 10000, R = 10 --write-replicates --output-directory T` (I/O included). Each run in ECAi mode; (i) and (ii) additionally in probability mode.
- Procedure per case: one discarded warm-up run, then 7 timed runs of `/usr/bin/time -f %e ./build-bench/ecai_cpp2 …` alternating before/after binaries (A B A B …) to share thermal/background conditions; no other user workload; report median and min/max. If the baseline's own max/min spread exceeds 5 %, repeat with 15 runs before judging.
- Gate: ECAi-mode median after > 1.05 × median before fails. Probability-mode median must be within ±10 % of the ECAi-mode median at the same settings (one `expm1` replaces one `log1p`; no allocation beyond the output vector).
- Record the table (commit hashes, toolchain, 7 raw timings per cell, medians) in `docs/VALIDATION.md` as part of commit B.

### 9.5 Schema compatibility

- ECAi `summary.csv` header line and `replicates.csv` header line byte-identical to baseline; `schema_version` still `stage13-v1`.
- No column added, renamed, or reordered in the ECAi schemas. Probability output is only ever written when `--mode probability` was parsed, resolved into `scenario.mode`, and carried by every replicate.

### 9.6 Mode isolation review checklist

- `git diff 7b93348 -- src/model.cpp` shows only appended code after line 206 (including the second anonymous-namespace block).
- `git diff 7b93348 -- src/seed_manager.cpp include/seed_manager.h src/sampling_plan.cpp include/sampling_plan.h src/lhs.cpp src/distributions.cpp src/statistics.cpp include/statistics.h src/category_data.cpp include/category_data.h` is empty.
- `git diff -w 7b93348 -- src/simulation.cpp src/output.cpp src/main.cpp` shows the ECAi branches as unchanged text (whitespace-insensitive), with only the dispatch/guard lines listed in §4.3, §5.1, §4.5 touching the ECAi flow.
- `run.config` construction in `cli.cpp:459-473` unchanged; `canonical_config_string` unchanged.
- `git diff --stat 7b93348` lists only: `include/{model,simulation,cli,output}.h`, `src/{model,simulation,cli,output,main}.cpp`, `CMakeLists.txt`, `tests/test_probability_mode.cpp` (new), `tests/test_cli.cpp` (flag pin only), the additive test files from commit A, `docs/*` and `README.md`.
- Every `switch (mode)` has no `default:` label so `-Wswitch` flags a future third enumerator.

### 9.7 Rollback criteria

Any of: a §9.1 mismatch on cases B-M (items i-v); any pre-existing ctest failure; any commit-A pin change; ECAi benchmark regression `> 5 %`; any ECAi schema byte change. Action: `git revert` commit B (single commit behind a default-off flag). Commit A stays as the standing baseline.

## 10. Scientific caveats (must appear in `docs/VALIDATION.md` and help text)

- Probability mode inherits every limitation of C++: aggregate-demand-conditioned marginal 1D-LHS, person-specific PBR, not demonstrated publication-equivalent. Results are not "the publication's probability".
- The reported quantity is the `p`-th percentile (0.96-quantile by default) of per-simulation infection probability at a fixed ECAi, not `1 - p` of anything, and not the probability that the standard's target is met.
- Comparing the mean replicate percentile with `target_probability` is a self-consistency diagnostic of the C++ sample at the published ECAi. Because Table 4 values are ceilinged to multiples of 5 and were produced by an unpublished sampling matrix, agreement or disagreement with `0.001` does not validate or invalidate either implementation.
- Zero-infector rows carry `P = 0` exactly; for categories with a large zero mass (Patient ≈ 91 % zero-infector) the `p = 0.96` percentile sits close to the point-mass boundary described in `docs/VALIDATION.md` and is sensitive to the sparse positive tail. Probability mode has fewer zero outputs than ECAi mode for the same sample (no negative clamping), so the two modes' zero fractions differ by design.
- `P == 1.0` is numerical saturation of `-expm1(-Q)` in double precision (`Q` beyond roughly 36.7 quanta), not a model statement that infection is certain; it can only occur at very small input ECAi or extreme sampled emission.
- Steady-state, well-mixed, constant-emission assumptions of Eqs. 1-7 apply. Duration `D` enters `Q` linearly; the target is hourly only at `D = 1`.
- The across-replicate CI is an approximate randomized-LHS Monte Carlo Student-t interval for the mean of the replicate-percentile estimator, not for the population percentile or for model-input uncertainty.
- `--occupants` overrides rescale total ECA (per-person input × effective `I0`) and suppress the comparator. `--target-probability` changes the sample (§6) as well as the comparator.

## 11. Open items to settle at implementation time (not blocking)

- Whether `input_ecai_source` should also distinguish "published but `--occupants` rescaled" (currently reported only through `custom_scenario` and `effective_eca_total_m3_h`).
- Whether to expose `saturated_count` in the console line. Current decision: CSV only.

## 12. Codex review disposition

| # | Point | Disposition |
| --- | --- | --- |
| 1 | Exclude `--help` from byte baseline | Incorporated: §9.1 case A recorded but not hash-compared; replaced by the T8b in-order subsequence check plus exit-code/no-summary gates. |
| 2 | Allow/define `P == 1` saturation | Incorporated: §3.2 closed range `[0, 1]`, `saturated` flag, `saturated_count`, `saturation_policy` column, T2b, caveat in §10. |
| 3 | Finite/overflow checks on all intermediates | Incorporated: §3.2 lists `ECA_total`, `Lambda`, `dose_num`, `Q`, `P_raw` each checked, `Q >= 0` asserted, overflow tests in T2b. |
| 4 | Decouple forward API from `target_probability` | Incorporated: §4.1 `ProbabilityInputs` has no target; validation is a separate appended helper; the target is used only by the comparator (§5.3). |
| 5 | Branch on resolved `scenario.mode` | Incorporated: §4.2 `run_replicates` uses `result.mode` (from `scenario.mode`); `main.cpp` and `build_document` branch on `run.scenario.mode` and verify `overrides.mode` and every replicate agree. |
| 6 | Complete CSV schema; no misleading counters | Incorporated: §5.2 enumerates all 96/20 columns with indices; clamp/ECAi-evaluation/ceiling columns removed from the probability schema; `output_zero_count` semantics defined in §3.2. |
| 7 | Mode-isolation checks for seeds, banks, assignments, cursors, QER | Incorporated: T5 (a)-(e) covers config/fingerprint, plan, bank seeds and values, assignment ids/cursors, and per-row QER sums. |
| 8 | Document target-probability fingerprint implications | Incorporated: §6 items 1-5; kept in fingerprint because removal would change ECAi seeds. |
| 9 | Clarify quantile terminology | Incorporated: §1.1; console prints `percentile=<number>` instead of a `p96` label. |
| 10 | Define comparator boundaries, duration, zero-infector behaviour, interpretation | Incorporated: §5.3. |
| 11 | stderr / exit-code / output-absence / signed-zero regression checks | Incorporated: §9.1 records five artefacts per case, adds H-M, `-0` regex, and `T`-absence checks. |
| 12 | Reproducible benchmark protocol | Incorporated: §9.4 (build type, toolchain record, warm-up, 7 alternating runs, spread rule, recorded table). |
| 13 | Resolve baseline-vs-feature commit structure | Incorporated: §8 two commits (A baseline+pins, B feature); §9.7 reverts B only. |
| 14 | Allow additive tests where required | Incorporated: §8 test policy permits additive `check` blocks in existing files; only the 17→19 pin is a modification. |

No review point was rejected.
