# C++2 QER Distribution Implementation Plan

Status: revised plan for review; no implementation changes made.

## Goal

Add an explicit C++2 option selecting one of three QER calculation families:

1. Full mechanistic Jones et al. calculation (current default).
2. Jones-fitted direct QER distributions using the fitted parameters in the
   legacy C++ implementation.
3. Mikszewski et al. direct QER distributions and profiles.

The existing full-Jones, community-prevalence, and risk-averse behavior must
remain unchanged when no new option is supplied.

## Required scientific distinction

Full Jones samples the existing eight infector-side QER inputs and computes
QER mechanistically.

Jones-fitted and Mikszewski modes use one direct log10-normal QER distribution
per realized infector:

```text
QER_i = 10 ** (mu_log10 + sigma_log10 * normal_ppf(u))
```

The direct values are in quanta/h per infector. They replace the infector-side
Jones emission calculation; susceptible-side inputs and the existing ECAi and
probability equations remain active.

Direct modes must use exactly one LHS value per realized infector. They must
not sample, construct, or validate the eight unused full-Jones QER inputs.

## Decisions required before implementation

- Choose one strict CLI form: one `--qer-distribution` selector, or separate
  family/profile options.
- Decide which legacy aliases are retained. Canonical names are preferred;
  aliases must be explicit and tested if retained.
- Decide whether bare `mikszewski` resolves to
  `mikszewski-sars-cov-2-standing-speaking`, matching legacy behavior, or
  requires an explicit profile.
- Decide whether QER family/profile is included in the seed fingerprint. The
  selected policy must preserve omitted-selector full-Jones reproducibility and
  prevent unintended direct-family stream collisions.
- Decide whether direct modes receive new schema variants or a versioned
  extension while preserving existing Jones documents unchanged.

## Phase 1: Baseline and source contract

Before coding:

- Freeze the current C++2 full-Jones default outputs, seed fingerprints, bank
  layouts, resource formulas, and representative serialized documents.
- Read and verify the legacy tables in:
  `C++/src/qer_fitted.cpp`.
- Read the C++2 sampling, seed, simulation, CLI, output, and documentation
  paths listed in the source audit.
- Do not modify the legacy C++ tree.

## Phase 2: Typed QER selection and parameter tables

Add a typed C++2 QER selection layer containing:

- `Jones` as the zero/default family value;
- `JonesFitted`;
- `Mikszewski` plus a canonical profile identifier;
- strict validation and canonical rendering;
- immutable parameter-table/version provenance.

Add the complete Jones-fitted table:

- exactly 25 records;
- `(mu_log10, sigma_log10)` values copied exactly from the legacy source;
- explicit mapping from legacy short names to C++2 canonical categories;
- no unknown-category fallback;
- rejection of incomplete mappings.

Required mapping includes, among others:

- `Food` -> `Food and beverage facilities`;
- `Transportation` -> `Transportation waiting`;
- `Lecture` -> `Lecture hall`;
- `Sorting` -> `Sorting, packing, light assembly`;
- `Exam` -> `Exam room`;
- `Group` -> `Group treatment area`;
- `Patient` -> `Patient room`;
- `Resident` -> `Resident room`;
- `Waiting` -> `Waiting room`;
- `Place` -> `Place of religious worship`;
- `Spectator` -> `Spectator area`;
- `Common` -> `Common space`;
- `Dwelling` -> `Dwelling unit`.

Add all 30 Mikszewski profiles: ten organisms crossed with resting,
standing-speaking, and light-speaking-loudly activity. Store the literal
median quanta/h and organism-specific `sigma_log10`; calculate
`mu_log10 = log10(median)` at initialization without rounding.

## Phase 3: Direct-QER sampling architecture

The sampling plan must be family-aware before allocation.

Full Jones:

- retain the current three susceptible banks and eight infector QER banks;
- retain the current bank order, seeds, cursors, assignments, and equations;
- retain the current hot path and default output behavior.

Jones-fitted and Mikszewski:

- allocate the three susceptible banks;
- allocate exactly one direct-QER LHS bank of length `K_total`;
- do not allocate or iterate over the eight full-Jones QER banks;
- transform each direct LHS value once while generating the direct bank;
- consume one direct value per realized infector in stable infector-list order.

Do not represent direct mode as eight empty Jones banks plus a side path. Use a
separate family-aware bank/assignment representation so unused structures are
not constructed or traversed.

Family-specific persistent-bank formulas are:

```text
Full Jones:        N + 3*N_active + 8*K_total
Direct QER modes:  N + 3*N_active + K_total
```

Resource checks must also include transient LHS vectors, transform scratch,
assignment storage, and any audit buffers. Validate the selected layout before
allocation, with checked arithmetic and no partial allocation.

Community zero-infector rows consume no susceptible or QER values. Risk-averse
runs retain `N_active == N`, no zero rows, and exactly one direct QER draw or
eight Jones draws per infector as appropriate.

## Phase 4: Seed and performance policy

Resolve family, profile, table version, and stream identity once before
sampling. Do not perform string lookup, map lookup, selector validation,
provenance construction, or schema checks inside the simulation hot loop.

The hot loop should perform only the selected indexed bank access, required
accumulation, and existing numerical checks. Direct log10-normal transforms
should occur during bank generation, not during ECAi or probability evaluation.

Add a dedicated direct-QER stream or equivalent direct-only seed domain. The
seed policy must explicitly define whether family/profile participates in the
fingerprint. The omitted-selector full-Jones path must retain its approved
seed and output compatibility; changing the seed design version or effective
configuration must not silently alter existing Jones runs.

Record provenance once per run/replicate metadata object. Keep per-value
consumption counters and final cursor checks, but make detailed consumption
IDs debug/test data unless the existing output contract requires them.

## Phase 5: CLI, scenario, and simulation integration

Add strict selector parsing and propagation through:

- CLI options and help text;
- `EffectiveOverrides` and resolved scenario;
- sampling-plan construction;
- replicate/batch consistency checks;
- simulation execution;
- console summaries and output metadata.

Reject unknown families, unknown profiles, family/profile combinations that are
inconsistent, duplicate conflicting selectors, and invalid categories before
allocation or output creation.

Branch once after plan construction:

- Full Jones dispatches through the existing `compute_qer()` path.
- Direct modes dispatch through the direct-QER bank path and never construct
  `InfectorParameters` or call full-Jones `compute_qer()`.

Apply the selected QER values identically in ECAi and Probability modes. Keep
infector-count sampling independent of QER-family selection.

## Phase 6: Output and comparator handling

Preserve existing full-Jones schemas and serialized output where compatibility
is required. Add deliberately versioned direct-family provenance containing:

- QER family;
- canonical Mikszewski profile, where applicable;
- parameter-table/version identifier;
- transform and units;
- direct-bank layout and stream/version;
- seed/fingerprint policy;
- direct values generated, consumed, and evaluated.

Do not duplicate expensive or unnecessary provenance fields per simulation row.

Existing Jones publication comparators are not valid for fitted or Mikszewski
alternatives. Mark alternative-family comparator results as not applicable
unless a separate comparator is explicitly approved.

## Phase 7: Tests before implementation acceptance

Add focused tests for every new function and data element:

- strict selector parsing, aliases, conflicts, and rejection;
- exactly 25 fitted records and complete canonical mapping;
- exactly 30 Mikszewski profiles and all organism/activity combinations;
- literal medians, sigmas, units, known `u` transforms, positivity, finiteness,
  and endpoint behavior;
- direct-bank length `K_total` and one LHS/transform per infector;
- no direct bank for `K_total == 0`;
- no unused Jones-bank generation in direct modes;
- full-Jones three-plus-eight bank layout unchanged;
- zero-infector short-circuit and risk-averse `N_active == N`;
- cursor exhaustion, overrun, under-consumption, and resource limits;
- same-seed reproducibility and approved seed migration behavior;
- evaluation-order-independent direct assignment;
- CLI/scenario propagation and mixed-family rejection;
- ECAi and Probability operation for all three families;
- output schemas, provenance, comparator suppression, and rejected-run cleanup.

Run focused tests before accepting implementation, then the complete registered
C++2 CTest suite and a bounded CLI smoke matrix covering default Jones,
explicit Jones, Jones-fitted, representative Mikszewski profiles, both output
modes, both infector policies, invalid selectors, and resource rejection.

## Final review gate

After implementation and tests, obtain an independent read-only review covering:

- all 25 fitted records and mappings;
- all 30 Mikszewski profiles;
- direct versus full-Jones bank allocation;
- hot-loop overhead;
- seed compatibility;
- community and risk-averse behavior;
- output schemas and comparator handling;
- absence of legacy C++ changes.

No merge, commit, push, or production-scale sweep is authorized until the
focused tests, full CTest suite, bounded smoke matrix, and independent review
pass.

## Source references

Primary legacy references:

- `C++/src/qer_distribution.cpp`
- `C++/src/qer_fitted.cpp`
- `C++/src/model.cpp`
- `C++/tests/test_qer_fitted.cpp`

Primary C++2 references:

- `C++2/src/sampling_plan.cpp`
- `C++2/src/simulation.cpp`
- `C++2/src/model.cpp`
- `C++2/src/cli.cpp`
- `C++2/src/output.cpp`
- `C++2/docs/SEED_POLICY.md`
- `C++2/docs/SAMPLING_DESIGN.md`

This plan incorporates the Opus 5 performance review correction: direct modes
must be a genuinely separate low-overhead sampling and evaluation path, not a
full-Jones plan with unused banks or objects retained.
