# Risk-Averse Infector Mode Implementation Plan

> For Hermes: implement this plan in small, dependency-ordered Sonnet phases after independent plan review.

Goal: Add an optional risk-averse infector mode to C++ that conditions every simulation on at least one infector while preserving the current community-infection behavior as the default.

Architecture: Add an explicit resolved scenario policy for infector handling, independent of the existing ECAi/probability output mode. In risk-averse mode, precompute the binomial zero mass before consuming LHS values, apply a conditional PPF directly in CDF space, and generate exactly one infector-count result per requested simulation without rejection loops or variable LHS consumption. Apply the same policy to both ECAi and infection-probability calculations.

Scope and constraints:

- Work only in the C++ tree: `C++/include`, `C++/src`, `C++/tests`, and relevant C++ documentation.
- Read the legacy C++ implementation for behavioral reference only; do not modify the legacy C++ tree.
- Preserve the current community mode and existing ECAi/probability output modes by default.
- Do not commit or push until implementation and verification are complete.
- Keep scientific validation and audit counters; do not optimize away checks that establish correct conditional sampling.

## Required semantics

1. Add an explicit infector policy with the current community behavior as the default and a risk-averse/require-infectors option.
2. In community mode, retain the current binomial distribution, including zero-infector rows and their existing `+0.0` handling.
3. In risk-averse mode, for category parameters `I0` and prevalence `p`, determine before sampling the zero mass using the same binomial PPF table that will perform inversion:

   `zero_mass = CDF_table[0] = P(X = 0)`

   Risk-averse mode requires `I0 >= 1` and rejects `p <= 0` before LHS generation or allocation. Accept `0 < p <= 1`. For `p=1`, return `I0` deterministically. If the zero mass rounds to `1.0` in binary64 for a representable positive `p`, use an explicit limiting branch that returns `1` and records the condition; do not allow this result to arise accidentally from a `0 >= 0` comparison. Community mode retains its existing domain and behavior.
4. Do not form an affine `u_conditional` in floating-point arithmetic. Define the conditional PPF mathematically as the smallest `k` in `[1, I0]` satisfying:

   `CDF_table[k] - zero_mass >= u * (1 - zero_mass)`

   Implement that definition using retained conditional-tail cumulative masses `T[k] = sum(pmf[1..k]) / sum(pmf[1..I0])`, built from the same PMF terms with checked/compensated accumulation, and return the first `k` with `T[k] >= u`. Do not subtract separately rounded CDF values near `zero_mass==1`. Start the search at `k=1`, retain a checked `k=I0` fallback, and make the helper a `BinomialPpfTable` member with a `zero_mass()` accessor. For `p=1`, short-circuit before CDF access, record `zero_mass=0.0`, set the deterministic-branch flag, and return `I0`. This structurally prevents zero-infector results and avoids lower-cutoff and upper-endpoint rounding hazards. The mathematically equivalent conditional-CDF interval is `(zero_mass, 1)`; its lower endpoint is exclusive under the existing `CDF >= u` convention.
5. Use exactly the requested number of LHS values for each category and simulation count. Do not implement a retry loop that consumes an unknown number of LHS values. The conditional PPF must handle `u=denorm_min`, the first/last generated strata, and `nextafter(1,0)` without constructing `u=1` or silently returning `I0` because `zero_mass` rounded to 1.0.
6. Apply the policy before the shared ECAi or forward-probability calculation. The forward probability API must remain independent of `target_probability`.
7. Risk-averse probability results must still be bounded in `[0,1]`, permit numerical saturation at exactly `1.0`, reject non-finite intermediates, and preserve signed-zero normalization for any zero output.
8. Preserve auditability: record the policy, cutoff/provenance, effective infector counts, bank cursors, QER evaluations, and invariant results in output metadata without silently reusing community-mode counters.

## Implementation phases

### Phase 1: Read-only source and behavior audit

Files to inspect:

- `C++/docs/PROBABILITY_MODE_IMPLEMENTATION_PLAN.md`
- `C++/include/cli.h`, `C++/src/cli.cpp`
- `C++/include/simulation.h`, `C++/src/simulation.cpp`
- `C++/include/sampling_plan.h`, `C++/src/sampling_plan.cpp`
- `C++/include/distributions.h`, `C++/src/distributions.cpp`
- `C++/include/output.h`, `C++/src/output.cpp`
- `C++/tests/test_sampling_plan.cpp`, `C++/tests/test_simulation.cpp`, `C++/tests/test_cli.cpp`, `C++/tests/test_output.cpp`
- Corresponding legacy C++ binomial/risk-averse implementation, read-only.

Deliverable: a short implementation map before coding, confirming that `build_sampling_plan` owns count/LHS construction, that `BinomialPpfTable` has no current conditional helper, and that `SamplingPlan`, `EffectiveScenario`, and `ReplicateResult` need policy/cutoff provenance. Confirm the legacy defaults: legacy ECAi defaults to allowing zero infectors, while legacy `probability_ecai` defaults to requiring infectors; C++ must retain its current community default for both modes. No source changes.

### Phase 2: Pure cutoff and conditional-PPF helper (TDD)

Likely files:

- Modify `C++/include/distributions.h` and `C++/src/distributions.cpp`; make the helper a `BinomialPpfTable` member with `zero_mass()` and conditional-tail accessors.
- Add focused tests in `C++/tests/test_distributions.cpp` (use the existing target; do not create a new unregistered target).
- Update `kBinomialTableBuffers`/peak-byte accounting in `C++/src/sampling_plan.cpp` if retained conditional-tail storage adds a simultaneously live buffer; if storage is risk-averse-only, account for that branch explicitly.

Tests first:

- `p=0` is rejected before simulation in risk-averse mode; `p=1` and zero-mass-underflow cases have explicit deterministic/limiting tests.
- The conditional transform is monotonic and remains in the valid open/closed endpoint policy.
- Values at the first and last LHS strata are finite and produce counts at least one.
- The inverse-CDF result agrees with direct conditional-binomial probabilities at small fixtures.
- Large valid `I0`/`p` cases do not underflow or overflow silently.

Implement checked arithmetic, stable zero-mass evaluation, and a single shared `conditional_ppf_at_least_one` helper rather than duplicate formulas in CLI, sampling, and simulation code. The helper must use the PPF table's own `cdf_[0]`; do not calculate a separate `(1-p)^I0` value for the inversion boundary. Avoid unconditional subtraction of separately rounded `CDF[k] - zero_mass` near `zero_mass==1`; retain conditional-tail cumulative masses or use an explicitly justified limiting branch, and compare against higher-precision conditional-PMF/oracle fixtures.

### Phase 3: CLI and resolved scenario policy

Likely files:

- `C++/include/cli.h`, `C++/src/cli.cpp`
- `C++/include/simulation.h` and scenario-resolution code
- `C++/tests/test_cli.cpp`

Add an explicit option such as `--infector-policy community|risk-averse`, with community as the default. Resolve the policy once in the scenario and branch on the resolved policy, not raw CLI flags. Validate `--infector-policy risk-averse --community-rate 0` before simulation and again during authoritative scenario resolution; rejected commands must create no output files. Do not add legacy aliases unless source review establishes an existing C++ compatibility requirement.

Add authoritative validation below the CLI: the policy-aware sampling-plan entry point must independently reject `I0==0`, `p<=0`, non-finite values, and impossible policy/scenario combinations before LHS generation or allocation. Additive tests must verify:

- default parsing resolves to community;
- explicit risk-averse parsing resolves correctly;
- invalid policy values are rejected;
- help and accepted-flag enumeration agree;
- explicit policy provenance is recorded;
- the option works with both `--mode ecai` and `--mode probability`;
- unrelated scientific overrides retain their existing custom-scenario behavior;
- direct plan construction rejects invalid risk-averse inputs even when the CLI is bypassed.

### Phase 4: Sampling-plan integration

Likely files:

- `C++/include/sampling_plan.h`, `C++/src/sampling_plan.cpp`
- `C++/tests/test_sampling_plan.cpp`

Construct the conditional LHS/binomial plan before simulation rows consume parameter banks. For each category and requested `N`, use the PPF table's conditional operation and transform exactly `N` base LHS values. The plan must expose or validate an explicit count-stream `lhs_generated == N` and `lhs_consumed == N` invariant; `counts.size()==N` alone is insufficient and must remain distinct from parameter-bank/QER consumption counters. Verify:

- exactly `N` infector counts are produced and exactly `N` count-stream uniforms are generated/consumed;
- every risk-averse count is in `[1,I0]`;
- community counts remain unchanged for identical seeds and inputs;
- risk-averse bank sizing uses the realized positive counts and checked sums;
- zero-infector rows are absent from risk-averse plans, while community plans preserve them;
- susceptible/infector bank cursors and QER demand match the planned counts;
- repeat construction with the same seed is identical;
- the infector-count stream seed is unchanged between policies, so each risk-averse count is no lower than its matched community count; bank streams may differ because their lengths differ;
- policy is stored in the plan and checked against the resolved scenario during plan validation and replicate execution.

Add policy, cutoff, and numerical-branch provenance to `SamplingPlan`, `EffectiveScenario`, and `ReplicateResult`. Define the provenance fields and meanings explicitly: policy enum; stored `zero_mass`; `zero_mass_rounded_to_zero`; `zero_mass_rounded_to_one`; conditional-PPF algorithm/version; and whether a deterministic/limiting branch was used. Use a zero-valued enum member for community so existing value initialization remains safe. Validate plan policy against scenario policy in `validate_plan`, `validate_plan_context`, `run_replicate_from_plan`, and the batch runner. Add a batch-consistency check rejecting mixed policies. Preserve the canonical seed fingerprint and design version unchanged as specified above.

### Phase 5: ECAi and probability simulation integration

Likely files:

- `C++/include/simulation.h`, `C++/src/simulation.cpp`
- `C++/tests/test_simulation.cpp`

Branch on the resolved infector policy and resolved output mode. Use the same downstream equations and parameter-bank order after the infector-count plan is established.

Additive tests must verify:

- community ECAi behavior remains numerically and structurally unchanged;
- risk-averse ECAi contains no zero-infector rows;
- risk-averse probability contains no zero-infector rows;
- both modes consume exactly the planned LHS values;
- identical seeds and category inputs produce reproducible results within each policy;
- all risk-averse probabilities are finite and in `[0,1]`;
- saturation and `+0.0` normalization remain correct;
- QER sums, assignments, cursors, and evaluation counters reconcile;
- target probability is not passed into the forward probability calculation;
- risk-averse output suppresses population-level published comparators and reports a distinct `NOT_APPLICABLE_RISK_AVERSE` status;
- exact conditional-boundary tests cover `u` equal to each conditional cumulative boundary and `nextafter` immediately below and above it;
- invalid cases include unknown policy, risk-averse with `p=0`, risk-averse with `I0=0`, and tampered scenario/plan/output policy disagreement.

The policy is an estimand change, not a scientific custom-input override: do not fold it into the existing four-field `custom_scenario` calculation.

### Phase 6: Output schema and end-to-end behavior

Likely files:

- `C++/include/output.h`, `C++/src/output.cpp`, `C++/src/main.cpp`
- `C++/tests/test_output.cpp`, `C++/tests/test_end_to_end.cpp`
- C++ CLI/model documentation.

Add policy fields and cutoff provenance to distinct versioned risk-averse schemas; do not change the legacy community ECAi schema or its serialized bytes. Freeze these schema versions before implementation: `stage13-ra-v1` for risk-averse ECAi summary/replicate output and `stage13-prob-ra-v1` for risk-averse probability summary/replicate output. Community schemas remain `stage13-v1` and `stage13-prob-v1`. For each of the four risk-averse documents, enumerate exact ordered headers, column widths, and blank/value rules before coding; added provenance fields include policy, `zero_mass`, `zero_mass_rounded_to_zero`, `zero_mass_rounded_to_one`, conditional-PPF algorithm version, and limiting/deterministic-branch flag. For `p=1` no CDF table exists, so store `zero_mass=0.0` with the deterministic-branch flag set and validate that convention. `validate_stage13_document` must dispatch by schema version and require policy agreement across `EffectiveScenario`, `SamplingPlan`, `ReplicateResult`, and output metadata. Risk-averse documents must assert `zero_infector_count == 0`, `active_count == n`, `k_total >= n`, and `min(infector_counts) >= 1`. Comparator precedence is explicit: risk-averse status `NOT_APPLICABLE_RISK_AVERSE` wins over custom status because the estimand mismatch is fundamental; all population-level comparator fields are blank. Custom provenance remains independently recorded and still controls comparator eligibility for community runs. Test exact headers, widths, row counts, finite values, policy markers, and absence of files after rejected invocations.

Run bounded real-executable cases for:

- community + ECAi;
- risk-averse + ECAi;
- community + probability;
- risk-averse + probability;
- explicit input ECAi and arbitrary percentile;
- invalid policy/scenario cases: unknown policy, risk-averse with `p=0`, risk-averse with `I0=0`, and tampered scenario/plan/output policy disagreement.

### Phase 7: Documentation and independent review

Update only relevant C++ documentation after behavior is verified:

- `C++/docs/CLI_REFERENCE.md`
- `C++/docs/MODEL_SPECIFICATION.md`
- `C++/docs/SAMPLING_DESIGN.md`
- `C++/docs/SEED_POLICY.md`: policy is not added to the canonical seed fingerprint or design version, preserving community seed/output regression behavior; the infector-count stream is therefore matched across policies, while bank lengths and downstream bank streams may differ.
- `C++/docs/VALIDATION.md`: document the distinct conditional estimand, endpoint/underflow safeguards, no-zero invariant, and legacy probability executable's different default (legacy `probability_ecai` requires infectors by default; C++ does not).

Document the cutoff equation, conditional LHS transformation, endpoint policy, exact-N guarantee, distinction from redraw-until-positive, default community behavior, risk-averse interpretation, seed fingerprint implications, limitations, and output metadata.

Dispatch an independent read-only review agent against the implementation, plan, legacy reference, tests, and generated outputs. Require it to report exact files read, unresolved scientific issues, and whether the community path is preserved.

## Verification gates

After each phase:

1. Inspect the agent report and `git diff`/`git status`.
2. Build the existing CMake build directory without deleting or recreating it.
3. Run the focused test target.
4. Run the full registered CTest suite before advancing when the phase changes shared code.
5. Verify no legacy C++ files changed.
6. Check `git diff --check`.

Final acceptance requires:

- all existing C++ tests pass;
- new focused cutoff, conditional-binomial, policy, simulation, output, and end-to-end tests pass;
- community ECAi output is regression-equivalent for representative fixed-seed cases, excluding intentionally changed help text/schema metadata where applicable;
- risk-averse ECAi and probability runs produce no zero-infector simulations;
- exact-N LHS consumption is demonstrated by cursor and plan invariants;
- cutoff and conditional-CDF comparison values are finite and reproducible;
- exact count-stream generated/consumed invariants equal `N`; parameter-bank consumption is checked separately;
- invalid invocations create no output files;
- generated output documents validate completely;
- benchmark results, if added, state build type, compiler/toolchain, warm-up, timing method, and machine-load assumptions;
- only intended C++ files are staged; unrelated artifacts such as `../masks/` remain excluded.

## Risks and open questions for review

- The exact CLI spelling should match any existing `require_infectors` or related option in the repository; do not add a duplicate alias without source evidence.
- The lower cutoff endpoint and LHS open-interval conversion must be defined precisely to avoid a zero-infector inverse-CDF result caused by rounding; the implementation must use conditional CDF-space inversion rather than affine endpoint construction.
- A transformed one-dimensional LHS is conditional stratification, not a full joint LHS over the binomial distribution; documentation must not overstate its statistical properties.
- The risk-averse mode changes the estimand from community prevalence to conditional risk given at least one infector; comparisons between policies must be labelled accordingly.
- The policy is deliberately excluded from the canonical seed fingerprint and design version to preserve community streams; this must be tested, documented, and kept consistent with the matched infector-count stream invariant.
