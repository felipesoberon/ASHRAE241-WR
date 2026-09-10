# C++2 Aggregate-Demand 1D-LHS ASHRAE 241 Model Implementation Plan

> **For Hermes:** GPT-5.6-sol is the coordinating brain. Implement this plan task-by-task with independent background GPT-5.6-luna or Claude Sonnet agents. Every implementation task requires a separate read-only verification agent before integration proceeds.

**Goal:** Build a new, CLI-only C++ implementation of the Jones et al. ASHRAE 241 Wells–Riley ECAi model that pre-samples infector counts and then generates exact-size, one-dimensional LHS lists for every required susceptible and infector parameter, so every generated value is consumed exactly once and no sampled value is shared between people.

**Architecture:** Each replicate first generates a one-dimensional LHS of size `N` and transforms it through the binomial inverse CDF to obtain one infector count per simulation. Zero-infector simulations are assigned ECAi = 0 immediately and consume no other sampled values. The program then computes aggregate demand across all active simulations, generates a separate one-dimensional LHS list of exactly that size for every named stochastic parameter, and assigns a distinct transformed value to each applicable person. Every infector receives a unique complete infector parameter set, while the susceptible receives a separate susceptible parameter set; pulmonary breathing rate and all other sampled values are never shared between susceptible and infector(s). This is aggregate-demand-conditioned marginal stratification, not a grouped or multidimensional joint LHS.

**Tech stack:** C++17 or later, CMake, CTest, C++ standard library only unless an external dependency is explicitly justified and approved. No Python runtime, GUI, or dependency on the existing `C++` binaries.

**Status:** Planning only. No model implementation or execution is authorized by this document.

---

## 1. Scope and isolation

### New implementation location

All new implementation work must remain under the `C++2` directory.

The existing implementation under the sibling `C++` directory

is read-only reference material. No agent may edit, format, rename, delete, rebuild, or generate output inside `C++` while implementing `C++2`.

### Included in the first version

- Pure C++ command-line program; no GUI and no Python production code.
- Jones et al. full mechanistic QER calculation.
- One occupancy category per invocation initially, with strict category lookup.
- Jones et al. defaults when flags are omitted.
- Default 10,000 simulations per replicate.
- Default 500 replicates.
- Default 96th percentile.
- Default hourly target infection probability of 0.001 (0.1%).
- Category-specific Jones community infection rate: 1% general and 3% healthcare.
- Jones category occupancy, volume, breathing-rate, droplet, mask, and ECAi inputs.
- Zero-infector simulations assigned ECAi = 0 without sampling any other parameter.
- Negative calculated ECAi clamped to zero.
- Mean, sample standard deviation, standard error, and confidence interval of replicate percentile estimates.
- Both the raw ECAi percentile in L/s/person and its ceiling to the next multiple of 5 L/s/person; statistical summaries remain based on unrounded raw values.
- Reproducible, collision-resistant seed derivation for every replicate and every LHS parameter list.
- CSV and concise console output with complete run metadata.
- Comprehensive unit, property, integration, CLI, statistical, and end-to-end tests.

### Explicitly excluded from the first version

- GUI code.
- Python production code.
- Changes to `C++`.
- Alternative QER distributions such as Jones-fitted or Mikszewski.
- Probability-at-fixed-ECAi calculations unless later requested.
- Bootstrap and order-statistic methods unless later requested.
- Multi-category batch execution inside the first executable unless later requested.
- Claims that conditional LHS is equivalent to the Jones publication’s original sampling architecture.
- Commit or push unless separately authorized.

---

## 2. Sources of truth and provenance

### Primary scientific source

The authoritative publication file is:

`2025 - Jones et al. - Risk modeling for ASHRAE Standard 241-2023 - control of infectious aerosols.pdf`

This publication is the source of truth for:

- Wells–Riley equations and units;
- ECAi inversion;
- target infection probability;
- one-hour exposure assumption;
- probabilistic input distributions;
- category values and category order;
- mask assumptions;
- community infection rates;
- percentile definition;
- published ECAi comparison values;
- any constant or transformation attributed to Jones et al.

Before implementing each scientific component, the implementation agent must inspect the relevant publication section/table/equation. It must cite the page, table, or equation in its report. A separate verifier must independently inspect the same source location.

### Secondary implementation references

Read-only reference files:

- `C++/src/model.cpp`
- `C++/src/model.h`
- `C++/src/random_manager.cpp`
- `C++/src/random_manager.h`
- `C++/src/ecai.cpp`
- `C++/analysis/ecai_uncertainty_study.cpp`
- `C++/tests/`
- `Python/model.py`
- `Python/randomManager.py`

These files may clarify the existing implementation but do not override the publication. Any discrepancy between the publication and existing code must be recorded and adjudicated by GPT-5.6-sol before implementation continues.

### Experimental design note

The preliminary conditional-LHS concept is documented in `conditional-lhs-by-infector-count.md`.

This note defines the proposed sampling architecture, not the Jones et al. scientific equations. It is not a source for physical constants or publication claims.

### Source hierarchy

1. Jones et al. publication for scientific model and defaults.
2. Published supplementary material, if referenced and available.
3. Existing Python reference implementation for implementation comparison.
4. Existing C++ implementation for implementation comparison.
5. Correspondence and temporary experiments for hypotheses only.

### Explicit design decision on person-specific values

Regardless of ambiguous wording in the publication or behavior in prior implementations, C++2 must not share a realized stochastic parameter between people. For every active simulation:

- the susceptible receives its own `PBR_s` draw;
- every infector receives its own complete infector parameter set, including its own `PBR_i` draw;
- no infector shares a realized parameter value with another infector by assignment;
- no infector parameter value is reused for the susceptible;
- using the same source distribution for two roles does not permit them to consume the same LHS value.

The frozen model specification must identify this as an explicit C++2 sampling-design decision wherever it differs from, or clarifies ambiguity in, the publication. Scientific equations and marginal distributions still follow Jones et al.

Do not silently resolve conflicting sources. Record the discrepancy, pause that task, and let GPT-5.6-sol adjudicate.

---

## 3. Required directory and file organization

The following directories have been created under `C++2`:

```text
C++2/
├── src/
├── include/
├── tests/
├── data/
│   └── reference/
├── analysis/
├── docs/
└── results/
```

Planned files and responsibilities:

```text
C++2/
├── CMakeLists.txt
├── README.md
├── IMPLEMENTATION_PLAN.md
├── .gitignore
├── include/
│   ├── category_data.h          # OccupancyCategory structure and strict lookup API
│   ├── cli.h                    # CLI options, parser, usage, validation
│   ├── distributions.h          # Inverse CDF interfaces and Jones transformations
│   ├── lhs.h                    # Exact-size 1D LHS generator and validation helpers
│   ├── model.h                  # QER and ECAi equations using supplied values
│   ├── sampling_plan.h          # Infector pre-sampling, exact draw counts, cursor accounting
│   ├── seed_manager.h           # Deterministic per-replicate/per-stream seed derivation
│   ├── simulation.h             # Replicate and study orchestration
│   ├── statistics.h             # Percentile, SD, SE, confidence interval
│   └── output.h                 # Console and CSV writers
├── src/
│   ├── main.cpp                 # Thin CLI entry point only
│   ├── cli.cpp
│   ├── category_data.cpp
│   ├── distributions.cpp
│   ├── lhs.cpp
│   ├── model.cpp
│   ├── sampling_plan.cpp
│   ├── seed_manager.cpp
│   ├── simulation.cpp
│   ├── statistics.cpp
│   └── output.cpp
├── tests/
│   ├── test_main.cpp            # Minimal test harness if no test framework is added
│   ├── test_category_data.cpp
│   ├── test_cli.cpp
│   ├── test_distributions.cpp
│   ├── test_lhs.cpp
│   ├── test_model.cpp
│   ├── test_sampling_plan.cpp
│   ├── test_seed_manager.cpp
│   ├── test_simulation.cpp
│   ├── test_statistics.cpp
│   ├── test_output.cpp
│   └── test_end_to_end.cpp
├── data/
│   └── reference/
│       ├── README.md             # Provenance, units, generation rules
│       └── publication_cases.csv # Only independently verified reference cases
├── analysis/
│   ├── README.md                 # Analysis executable purpose and reproducibility
│   └── sampling_audit.cpp        # Optional C++-only audit executable; no Python
├── docs/
│   ├── MODEL_SPECIFICATION.md    # Current equations, units, defaults, provenance
│   ├── SAMPLING_DESIGN.md        # Conditional-LHS design and limitations
│   ├── CLI_REFERENCE.md
│   ├── SEED_POLICY.md
│   └── VALIDATION.md
└── results/
    └── .gitkeep                  # Generated outputs ignored except placeholder, if desired
```

Rules:

- Header files contain public types and interfaces; implementations remain in `src`.
- Scientific constants and category data belong in `category_data`, not `main.cpp`.
- Distribution transforms belong in `distributions`, not the simulation loop.
- LHS generation must be independent of the physical distribution transforms.
- Seed derivation must be centralized; no component may improvise `seed + index` independently.
- `results/` contains generated data and must be ignored by Git by default.
- `data/reference/` contains only small, reviewed fixtures with documented provenance.
- Analysis code must be C++ and must read model outputs rather than rerun simulations invisibly.

---

## 4. Proposed scientific and sampling workflow

### 4.1 Replicate setup

For each replicate:

1. Derive a unique binomial-LHS seed from:
   - user-provided base seed;
   - replicate index;
   - fixed stream identifier for `ninf`.
2. Generate exactly `N` one-dimensional LHS uniform values strictly in `(0,1)`.
3. Transform each uniform value using the binomial inverse CDF:

   `ninf[i] = BinomialPPF(u[i], I0, C_i)`

4. Store one infector count for each simulation index.
5. Compute:
   - `N_zero = count(ninf == 0)`;
   - `N_active = count(ninf > 0)`;
   - `K_total = sum(ninf)`.

Required invariant:

`N_zero + N_active == N`

### 4.2 Exact parameter-list sizes

Zero-infector simulations require no further parameter values and receive ECAi = 0.

For active simulations, generate exact-size one-dimensional LHS lists based on aggregate demand across the replicate. These are not separate LHS designs for each realized `ninf = k` group.

Susceptible/simulation-level stochastic parameters use list size `N_active`:

- susceptible pulmonary breathing rate `PBR_s`;
- biological decay/removal parameter `lambda`;
- deposition/removal parameter `gamma`.

Infector-specific QER stochastic parameters use list size `K_total`:

- infector pulmonary breathing rate `PBR_i`;
- droplet concentration `C_drop`;
- droplet diameter `d`;
- multiplicative diameter factor `E`;
- genomic viral load `GVL`;
- viable fraction `VF`;
- respiratory tract deposition `RTD`;
- dose constant `DK`.

Publication constants such as exposure duration `D` and the dimensionless fractional viability after emission `R`, if fixed by Jones et al., consume no LHS values. `R` must not be described as respiratory deposition; respiratory-tract deposition is represented separately by `RTD`.

Each named physical parameter receives its own one-dimensional LHS list and its own derived seed. Parameters that use the same distribution family still receive separate lists. No generic distribution-family buffer is shared among named parameters.

Each list is marginally stratified. Independent stream permutations form simulation-level parameter combinations, within-infector parameter vectors, susceptible–infector pairings, and multi-infector QER sums; those combinations are not a multidimensional joint LHS and must not be represented as such.

### 4.3 Assignment and consumption

- Maintain one cursor per parameter list.
- Iterate simulations in stable original index order.
- If `ninf == 0`, write ECAi = 0 and advance no susceptible or infector cursors.
- If `ninf > 0`, consume exactly one value from each susceptible list.
- Consume exactly `ninf` values from each infector-specific list.
- Never reuse a value.
- Never refill a list during a replicate.
- Abort with a diagnostic if a cursor would exceed its planned list size.
- At replicate completion, require every cursor to equal its list size.
- Stable simulation-index order and contiguous cursor consumption are the normative assignment algorithm. If calculation is later parallelized or traversed in another order, parameter assignments must first be fixed by simulation/person identity so traversal order cannot change results.

Core acceptance invariants:

- Every generated parameter value is consumed exactly once.
- No parameter value is consumed twice.
- No non-binomial parameter value is generated for a zero-infector simulation.
- Susceptible parameters are independent draws from infector parameters.
- Every infector receives its own distinct value from every infector-specific bank; values are not shared between infectors.
- Different named parameters do not use the same LHS sequence or seed.
- Replicates do not use the same LHS sequence or seed.

### 4.4 QER and ECAi equations

The implementation must reproduce the publication equations and units. The existing C++ implementation currently represents the Jones calculation as:

- `Vdrop = (pi / 6) * (d * E)^3 * C_drop`
- `GVL_m3 = GVL_ml * 1e6`
- `VER = PBR_i * Vdrop * GVL_m3 * VF * R`
- `QER_i = RTD * VER / DK`
- `QER_sum = sum(QER_i)` across realized infectors
- `Q_target = -log(1 - target_probability)`
- `mask_factor = (1 - mask_efficiency)^2`
- `term1 = (PBR_s * D * mask_factor / Q_target) * QER_sum`
- `term2 = VOL * (gamma + lambda)`
- `ECAi = (term1 - term2) / (3.6 * I0)`
- `ECAi = max(ECAi, 0)`

These expressions are provisional until independently checked against the publication. The frozen specification must define `R` using the publication's terminology and document a dimensional derivation showing `Vdrop` in m³ respiratory fluid per m³ exhaled air, `VER` in viable virions/h/person, `QER` in quanta/h/person, both ECAi inversion terms in m³/h, and division by `3.6 * I0` yielding L/s/person. Every unit conversion, especially `3.6`, `1e6`, hours, metres, litres, and ECAi units, must be documented and tested.

### 4.5 Percentile and replicate summary

Within each replicate:

- retain all `N` ECAi values, including zeros;
- compute the requested percentile using a documented interpolation convention;
- default to `p = 0.96`;
- do not round before calculating or aggregating.
- preserve and report the raw percentile in L/s/person;
- after obtaining the across-replicate mean raw percentile, also report `5 * ceil(mean_raw_percentile / 5)` in L/s/person as the ceiling-to-5 result;
- do not average ceiling-adjusted replicate values and do not use adjusted values to calculate SD, SE, or confidence intervals.

Across `R` replicates:

- mean of replicate percentile estimates;
- sample standard deviation using denominator `R - 1`;
- standard error `SD / sqrt(R)`;
- default two-sided 95% confidence interval for the mean using Student-t with `R - 1` degrees of freedom;
- minimum and maximum replicate percentile;
- exact replicate count and finite-value checks.

The confidence interval must be labelled as an approximate randomized-LHS Monte Carlo confidence interval for the mean of the replicate percentile estimator. It is not a confidence interval for the population percentile, an interval for individual ECAi values, model-input uncertainty, or protection against percentile-estimator bias. The approximate-normality assumption for replicate estimates must be documented and assessed diagnostically.

---

## 5. Default inputs and CLI design

### 5.1 Defaults

Publication defaults must be verified before encoding. Planned defaults:

- category: `Patient` only if the publication/current workflow supports this as a sensible single-category default; otherwise require `--category` explicitly;
- simulations: `10000`;
- replicates: `500`;
- percentile: `0.96`;
- confidence level: `0.95`;
- target infection probability: `0.001`;
- duration: `1.0` hour in Jones mode; any override makes the run a custom-duration scenario whose target probability applies over the selected duration and must not be labelled an hourly target;
- community rate: category default, 0.01 general and 0.03 healthcare;
- base seed: publication-compatible documented value, provisionally `485`, subject to source verification;
- negative ECAi clamping: enabled and not optional in version 1;
- zero-infector short circuit: enabled and not optional in version 1;
- Jones full QER distributions only.

### 5.2 Planned flags

```text
--category NAME
--simulations N
--replicates R
--percentile P
--confidence LEVEL
--target-probability P
--duration-hours HOURS
--base-seed SEED
--community-rate RATE
--occupants I0
--volume-m3 VOLUME
--mask-efficiency FRACTION
--output-directory PATH
--summary-file PATH
--replicates-file PATH
--write-replicates
--help
```

Rules:

- Strict full-string numeric parsing; reject trailing text.
- Reject `N <= 0`, `R <= 0`, nonpositive or non-integer occupants, unknown categories, community rates outside `[0,1]`, target probabilities outside `(0,1)`, confidence levels outside `(0,1)`, percentiles outside the frozen supported interval, volume `<= 0`, nonpositive duration, mask efficiencies outside `[0,1]`, non-finite values, derived-size overflow, and seed overflow.
- `--community-rate` overrides the selected category directly.
- Category defaults provide the Jones general or healthcare rate; no overlapping general/healthcare CLI override flags are needed for a one-category invocation.
- Any override of category scientific inputs marks the run as a custom scenario. Output must retain the category-default values, report effective values separately, and qualify or suppress direct comparison with the published category ECAi.
- Output paths must be explicit; no hidden writes into the source tree.
- Console output remains concise and machine-readable enough for automation.

Distribution parameters should not all become CLI flags in version 1 unless requested. Category-level overrides listed above are sufficient for the first implementation; adding every distribution parameter prematurely would enlarge the validation surface and weaken provenance.

---

## 6. Seed architecture

### Requirements

- Same complete CLI configuration and base seed must reproduce bit-for-bit identical C++2 outputs on the same supported implementation/toolchain.
- Different replicate indices must produce different binomial and parameter LHS sequences.
- Different named parameter streams within a replicate must produce different sequences.
- Different effective configurations and list lengths must not silently reuse the same stream seed.
- Thread scheduling must not affect results if parallelism is added later.

### Proposed policy

Create a `SeedManager` that derives each stream seed from a stable tuple:

```text
(base_seed, design_version, effective_configuration_fingerprint,
 replicate_index, stream_id, list_length)
```

Use a documented 64-bit mixing function such as SplitMix64 for derivation. Give every named stochastic parameter a fixed enum `stream_id`. Define a canonical, versioned encoding of effective configuration inputs instead of `std::hash<string>`, whose stability is not guaranteed. Do not use the same RNG object implicitly across unrelated modules. Seed uniqueness is verified over the supported domain rather than claimed as a mathematical collision impossibility.

Planned stream identifiers include:

- binomial infector counts;
- susceptible PBR;
- gamma;
- lambda;
- infector PBR;
- C_drop;
- d;
- E;
- GVL;
- VF;
- RTD;
- DK.

Tests must enumerate all seeds for at least 500 replicates and assert uniqueness for every supported `(configuration, replicate, stream, list_length)` test tuple. Tests must separately compare generated vectors and stratum permutations because different seed integers do not mathematically guarantee different finite sequences. Before coding the LHS engine, Stage 0 must freeze either (a) a fully specified random engine, open-interval integer-to-double conversion, bounded-integer generation, and Fisher–Yates shuffle for cross-toolchain reproducibility, or (b) a narrower reproducibility promise tied to named compiler and standard-library versions.

---

## 7. Detailed staged implementation plan

Each stage follows this gate:

1. GPT-5.6-sol supplies the bounded task and verified source evidence.
2. One GPT-5.6-luna or Sonnet background agent writes failing tests first.
3. The same agent implements only that task and runs focused tests.
4. GPT-5.6-sol inspects the actual diff and test output.
5. A different GPT-5.6-luna or Sonnet agent performs a read-only specification and code-quality review.
6. If the verifier reports a defect, GPT-5.6-sol adjudicates against the publication and dispatches a focused repair agent.
7. A fresh verifier reruns the affected tests.
8. Work proceeds only after an explicit PASS.

No agent may commit or push without separate authorization.

### Stage 0: Source extraction and frozen specification

**Objective:** Produce a verified, citation-based model specification before coding equations.

**Files:**

- Create: `docs/MODEL_SPECIFICATION.md`
- Create: `data/reference/README.md`

**Implementation-agent tasks:**

1. Read the publication’s methods, equations, tables, appendices, and relevant supplementary references.
2. Extract every category input, stochastic distribution, constant, unit, mask rule, community rate, percentile, and target probability.
3. Compare the publication values with `C++/src/model.cpp` and `Python/model.py` read-only.
4. Record exact discrepancies without resolving them silently.
5. Produce an explicit parameter table and equation list with page/table/equation citations.
6. Record the C++2 decision that susceptible PBR and every infector's PBR are separate draws, and that every infector receives a separate complete parameter set.
7. Correctly define publication symbol `R` as fractional viability after emission and distinguish it from `RTD`.
8. Freeze a dimensional-analysis table for every QER and ECAi intermediate.
9. Freeze the percentile interpolation convention and the formula `ceiling_5 = 5 * ceil(raw / 5)`.
10. Freeze the portable RNG/open-interval/shuffle policy and its reproducibility scope.
11. Create a validation acceptance table before scientific coding begins. At minimum it must specify tolerances for deterministic algebra, every inverse CDF, binomial PPF, percentile/statistical fixtures, and publication comparisons.
12. Classify Jones Table 4 values as adjusted publication comparators rather than exact raw Monte Carlo fixtures, and record the publication's adjustment rule with a citation.

**Independent verification:**

- A separate agent checks every extracted value against the PDF.
- GPT-5.6-sol adjudicates discrepancies.
- Gate: no scientific code until all required defaults have verified provenance and the equation, units, random-number policy, percentile convention, rounding rule, and quantitative validation tolerances are frozen.

### Stage 1: CMake project and test harness

**Objective:** Create a buildable C++2 skeleton with executable and test targets, without model logic.

**Files:**

- Create: `CMakeLists.txt`
- Create: `.gitignore`
- Create: `src/main.cpp`
- Create: `tests/test_main.cpp`
- Create: `README.md`

**Tests first:**

- Test executable runs and reports zero failures.
- Main executable supports `--help` only.
- Build does not reference any `C++` source or library.
- Before creating code, record the repository Git status and a file inventory plus hashes for the existing `C++` tree in the verification report; repeat the comparison after every stage so C++2-attributable changes can be detected.

**Commands planned:**

Set `CPP2` to the local C++2 directory, then run these commands in WSL:

```bash
cmake -S "$CPP2" -B "$CPP2/build"
cmake --build "$CPP2/build" -j2
ctest --test-dir "$CPP2/build" --output-on-failure
```

**Independent verification:** Confirm source isolation, CMake target scope, and absence of GUI/Python production files.

### Stage 2: Category data and defaults

**Objective:** Implement strict category lookup using the frozen publication specification.

**Files:**

- Create: `include/category_data.h`
- Create: `src/category_data.cpp`
- Create: `tests/test_category_data.cpp`

**Tests first:**

- Exactly 25 expected categories in publication order.
- Every category’s I0, volume, community rate, mask, breathing-rate parameters, droplet parameters, and published ECAi match the frozen specification.
- Museum and Convention use the verified maximum occupancy of 400 if confirmed by the publication.
- Healthcare detection is explicit, not inferred from arbitrary thresholds.
- Unknown category is rejected; no silent fallback category.

**Independent verification:** Recheck a sample from every category family and all healthcare rows against the PDF.

### Stage 3: Exact-size one-dimensional LHS engine

**Objective:** Implement a reusable 1D LHS generator independent of distributions.

**Files:**

- Create: `include/lhs.h`
- Create: `src/lhs.cpp`
- Create: `tests/test_lhs.cpp`

**Algorithm:** For size `n`, generate one jittered point strictly inside every stratum `(i/n, (i+1)/n)`, then apply the frozen deterministic Fisher–Yates shuffle. Construct jitter from RNG integers using the Stage 0 open-interval mapping; do not rely on implementation-dependent `std::uniform_real_distribution` if cross-toolchain known-answer reproducibility is selected.

**Tests first:**

- Reject size zero.
- Output length exactly `n`.
- Every value finite and strictly in `(0,1)`.
- `floor(u*n)` covers every stratum exactly once.
- Same seed gives identical values.
- Different seeds give different values and different rank permutations.
- Small sizes `1`, `2`, `3`; normal size `10,000`; odd sizes.
- No duplicate sequence across named test streams.
- Endpoint handling never yields exactly 0 or 1.
- Known-answer tests cover the first and last strata and the smallest/largest RNG integer mappings.
- The frozen shuffle known-answer vector passes on every supported toolchain.

**Independent verification:** Compare structural properties with SciPy `qmc.LatinHypercube(d=1).random(n)` using a temporary external test only; do not require bit-for-bit SciPy equivalence unless later specified.

### Stage 4: Seed manager

**Objective:** Derive reproducible streams and verify distinct seeds and sequences across the supported configurations, replicates, parameters, and list lengths.

**Files:**

- Create: `include/seed_manager.h`
- Create: `src/seed_manager.cpp`
- Create: `tests/test_seed_manager.cpp`
- Create: `docs/SEED_POLICY.md`

**Tests first:**

- Stable known-answer vectors for several tuples.
- No collisions across 500 replicates and all defined stream IDs.
- Changing base seed, replicate, or stream changes the derived seed.
- Changing effective configuration, design version, or list length changes the derived seed.
- Same tuple is deterministic.
- Stream IDs are unique and exhaustive.

**Independent verification:** Audit integer overflow semantics, canonical configuration encoding, and confirm no platform-dependent `std::hash` use. Distinguish verified seed uniqueness from verified sequence uniqueness and avoid mathematical guarantees that a finite hash cannot provide.

### Stage 5: Distribution inverse CDFs

**Objective:** Transform LHS uniforms into Jones et al. parameter distributions.

**Files:**

- Create: `include/distributions.h`
- Create: `src/distributions.cpp`
- Create: `tests/test_distributions.cpp`

**Functions:**

- normal inverse CDF;
- lognormal inverse CDF;
- base-10 lognormal inverse CDF;
- uniform inverse CDF;
- beta inverse CDF;
- binomial inverse CDF.

**Tests first:**

- Reject `u <= 0` and `u >= 1` for unbounded continuous transforms; test values approaching both endpoints.
- Invalid parameter rejection.
- Known quantiles and reference vectors generated from authoritative SciPy functions.
- Monotonicity for every inverse CDF.
- Binomial smallest-`k` convention: return smallest `k` with `CDF(k) >= u`.
- Binomial boundaries `p=0`, `p=1`, `n=0`, `I0=400`, and rates 0.01/0.03.
- Numerical stability and finiteness in Jones tail ranges.
- Explicit iteration limits and non-convergence diagnostics for iterative beta, binomial, and Student-t support calculations.
- Dense monotonicity grids near both tails and independently generated tail-reference fixtures.
- Function-specific absolute/relative tolerances from the frozen Stage 0 validation table.

**Independent verification:** Compare a bounded grid against SciPy reference outputs and document tolerances.

### Stage 6: Binomial pre-sampling and draw-count plan

**Objective:** Generate all `ninf` values first and calculate exact parameter-list sizes.

**Files:**

- Create: `include/sampling_plan.h`
- Create: `src/sampling_plan.cpp`
- Create: `tests/test_sampling_plan.cpp`

**Tests first:**

- Exactly `N` binomial LHS values and `N` infector counts.
- Every count in `[0, I0]`.
- `N_zero + N_active == N`.
- `K_total == sum(ninf)`.
- All-zero synthetic case: susceptible and infector draw counts both zero.
- Mixed known case with expected counts.
- All-infected case.
- 10,000 simulations, I0=400, p=0.01.
- Draw plan sizes: each susceptible stream `N_active`; each infector stream `K_total`.
- Determinism and cross-replicate difference.
- Checked summation for `K_total`, checked `N * I0`, checked vector byte counts, and graceful rejection before allocation when a configured resource limit or `size_t` bound would be exceeded.
- Explicit practical CLI limits or a documented pre-allocation memory-budget check; no reliance on publication defaults to make arbitrary overrides safe.

**Independent verification:** Confirm binomial transform convention and manually reconcile several small fixtures.

### Stage 7: Exact-size parameter banks and cursor accounting

**Objective:** Generate and consume every required LHS value exactly once.

**Files:**

- Extend: `include/sampling_plan.h`
- Extend: `src/sampling_plan.cpp`
- Extend: `tests/test_sampling_plan.cpp`

**Tests first:**

- One exact-size vector per named physical parameter.
- Correct transformed distribution per vector.
- Distinct derived seed per vector.
- No exact vector/rank-sequence equality among same-size parameter vectors.
- Cursor underflow/overflow detected.
- Zero-infector row advances no cursor.
- Active row advances each susceptible cursor once.
- Active row with `k` infectors advances each infector cursor exactly `k` times.
- All cursors exhausted exactly at replicate completion.
- No refill path exists.
- Instrumented consumption IDs show no duplicate use and no unused IDs.

**Independent verification:** Audit the complete cursor map and run property tests over randomized synthetic `ninf` vectors.

### Stage 8: Jones QER calculation

**Objective:** Implement QER as a pure function of supplied sampled values.

**Files:**

- Create: `include/model.h`
- Create: `src/model.cpp`
- Create: `tests/test_model.cpp`

**Tests first:**

- Known hand-calculated QER fixture.
- Each factor’s unit/conversion applied once.
- `Vdrop`, `GVL_m3`, `VER`, and `QER_i` intermediate values.
- Positive finite output over distribution boundary fixtures.
- Scaling tests: doubling appropriate multiplicative inputs doubles QER.
- No RNG or LHS generation inside model equations.

**Independent verification:** Compare algebra term-by-term with the publication and existing Python implementation.

### Stage 9: ECAi calculation and zero shortcut

**Objective:** Calculate one active simulation’s ECAi from supplied parameters and infector QER values, while keeping zero-infector consumption control in orchestration.

**Files:**

- Extend: `include/model.h`
- Extend: `src/model.cpp`
- Extend: `tests/test_model.cpp`

**Tests first:**

- The pure algebra helper returns exactly zero for an explicitly supplied zero-QER case.
- One-infector hand calculation.
- Multiple-infector QER sum.
- Negative raw ECAi clamps to exactly zero.
- Positive ECAi remains unchanged.
- Mask factor applied according to publication.
- Target probability boundary validation.
- Duration and unit-conversion tests.
- Finite output requirement.

The model equations do not own parameter-bank cursors and therefore cannot prove that no values were requested. `simulation.cpp`/`sampling_plan.cpp` must short-circuit a zero-infector row before requesting susceptible or infector parameters. In an all-zero replicate, no non-binomial parameter bank is allocated. In a mixed replicate, aggregate banks exist only for active rows and no cursor advances for a zero row.

**Independent verification:** Publication equation audit and direct comparison with a deterministic fixture evaluated independently.

### Stage 10: Replicate simulation

**Objective:** Combine pre-sampling, parameter banks, model calculation, and percentile extraction.

**Files:**

- Create: `include/simulation.h`
- Create: `src/simulation.cpp`
- Create: `tests/test_simulation.cpp`

**Tests first:**

- Small deterministic replicate with explicit `ninf` fixture.
- Output contains exactly `N` ECAi values.
- Zero count in ECAi vector equals zero-infector count plus any clamped negative results, tracked separately.
- All planned parameter values consumed exactly once.
- Every infector receives distinct consumption IDs for every infector-specific parameter, and none matches the susceptible consumption ID.
- 96th-percentile interpolation known-answer fixtures.
- Same seed/config identical; changed seed different.
- No state leakage between replicates.
- Category defaults and explicit overrides.
- Stable index mapping is recorded before calculation; changing computational traversal order does not change person-to-value assignments or outputs.

**Independent verification:** Inspect run metadata, counters, cursor exhaustion, and one small real replicate.

### Stage 11: Across-replicate statistics

**Objective:** Run multiple replicates and summarize their percentile estimates.

**Files:**

- Create: `include/statistics.h`
- Create: `src/statistics.cpp`
- Extend: `include/simulation.h`
- Extend: `src/simulation.cpp`
- Create: `tests/test_statistics.cpp`

**Tests first:**

- Percentile convention on odd/even and duplicate-valued vectors.
- Mean and sample SD known-answer vectors.
- Standard error.
- Student-t confidence interval known-answer fixtures.
- `R=1` policy explicitly handled: SD/CI unavailable or rejected, never fabricated.
- Default `R=500`, confidence 95%.
- No rounding before statistics.
- The approximate-normality assumption is documented and replicate-estimate diagnostics are available without changing the reported estimator.
- Raw mean percentile, raw SD/SE/CI, and `5 * ceil(raw_mean / 5)` known-answer fixtures; no adjusted value enters SD, SE, or CI calculations.

**Independent verification:** Compare bounded numeric fixtures with an independent statistical implementation.

### Stage 12: Strict CLI

**Objective:** Expose Jones defaults and safe overrides.

**Files:**

- Create: `include/cli.h`
- Create: `src/cli.cpp`
- Extend: `src/main.cpp`
- Create: `tests/test_cli.cpp`
- Create: `docs/CLI_REFERENCE.md`

**Tests first:**

- Every accepted flag and default.
- Unknown flag rejection.
- Missing value rejection.
- Too few/trailing characters, NaN, infinity, overflow, negative, and zero boundary cases.
- Category case policy.
- Override precedence.
- Duration override marks the target as applying over the selected duration and prevents it from being labelled hourly unless duration is exactly one hour.
- Every scientific override marks the run as custom, preserves source-category defaults in metadata, and reports effective values separately.
- `--help` succeeds without running simulations.
- Error exits are nonzero and concise.

**Independent verification:** Run a matrix of valid/invalid commands and compare help text with parser behavior.

### Stage 13: Output and metadata

**Objective:** Produce reproducible summaries without filling the repository with raw data.

**Files:**

- Create: `include/output.h`
- Create: `src/output.cpp`
- Create: `tests/test_output.cpp`
- Update: `.gitignore`

**Default summary fields:**

- model/sampling-design version;
- category, source-category defaults, effective inputs, and Jones-default/custom-scenario status;
- simulations and replicates;
- percentile and confidence level;
- target probability, duration, and explicit `target_probability_over_duration` semantics;
- community rate and mask efficiency;
- base seed and seed-policy version;
- zero-infector policy and negative-clamping policy;
- raw mean percentile in L/s/person, raw SD, raw SE, raw CI bounds, raw min, and raw max;
- ceiling-to-5 result in L/s/person, calculated once as `5 * ceil(raw_mean_percentile / 5)`;
- published adjusted ECAi comparator and comparison status only for an unmodified Jones category run;
- finite-value and row-count checks;
- per-replicate zero-infector and total-infector diagnostics summarized appropriately.

**Tests first:**

- CSV quoting and column order.
- Exact row count.
- Round-trip parsing.
- Full-precision numeric output.
- Raw and ceiling-to-5 fields have unambiguous names and units.
- Ceiling boundary fixtures include zero, an exact multiple of 5, and values immediately above/below a multiple of 5.
- Confidence statistics are calculated from raw values only.
- No output created on invalid configuration.
- Raw replicate CSV only with explicit `--write-replicates`.
- Output directory cannot accidentally be `C++`.

**Independent verification:** Read generated files back and validate schema programmatically.

### Stage 14: End-to-end scientific validation

**Objective:** Verify the complete C++2 executable without claiming sampling equivalence to existing models.

**Files:**

- Create: `tests/test_end_to_end.cpp`
- Create: `data/reference/publication_cases.csv`
- Create: `docs/VALIDATION.md`

**Validation sequence:**

1. Tiny deterministic synthetic fixture.
2. One category, small N and R, to exercise every path.
3. Patient category at publication defaults with bounded replicates for smoke verification.
4. Museum or Convention to exercise high occupancy and larger infector totals.
5. Full default 10,000 x 500 only after all focused gates pass and the user authorizes execution.
6. Compare scientific equations and marginal distributions with the existing implementations.
7. Compare results while explicitly labelling the sampling architecture as different.
8. For every unmodified Jones category selected for publication validation, report the raw 96th-percentile mean, raw CI, ceiling-to-5 result, published adjusted value, absolute difference, and comparison classification.

**Quantitative acceptance policy:**

- Deterministic algebra fixtures must match independently calculated double-precision references within both the absolute and relative tolerances frozen in Stage 0; the initial target is absolute `1e-12` or relative `1e-10`, tightened or changed only with documented numerical evidence before implementation.
- Inverse CDF, binomial PPF, percentile, and Student-t fixtures use function-specific tolerances frozen in `docs/MODEL_SPECIFICATION.md`; tolerances may not be loosened after seeing C++2 results without documented root-cause evidence and independent approval.
- LHS structure, cursor accounting, seeds, output row counts, clamping, and zero-infector behavior are exact pass/fail invariants with no statistical tolerance.
- Jones Table 4 values are adjusted publication comparators, not exact raw Monte Carlo fixtures. The reproduction target is equality between C++2's ceiling-to-5 value and the published adjusted ECAi for an unmodified Jones case.
- A publication difference of one 5 L/s/person reporting increment is classified as `NEAR_WITHIN_ONE_INCREMENT` and requires investigation and documentation; it is not silently called an exact reproduction.
- A publication difference greater than 5 L/s/person is a failed scientific-validation gate. Implementation completion is blocked until the equation, inputs, masks, percentile convention, stochastic transforms, and sampling-design effect are independently investigated and GPT-5.6-sol adjudicates the cause.
- Algebraic correctness and distribution-transform correctness must pass even if a publication comparator differs because of the intentionally different sampling architecture.

**Acceptance checks:**

- No changes under `C++`.
- All CTest targets pass.
- Sanitizer build passes for small tests if supported.
- No cursor overrun, under-consumption, seed collision, NaN, or infinity.
- Zero-infector rows consume no parameter-bank values.
- All generated LHS values are consumed exactly once.
- Different replicates and parameter streams have different seeds and sequences.
- Same configuration is reproducible.
- CLI metadata exactly describes applied settings.
- Each active susceptible and each individual infector has its own consumption IDs and complete applicable parameter set.
- Raw and ceiling-to-5 ECAi fields satisfy the frozen formula exactly, while all variability statistics use raw values.
- Checked-size/resource failures occur before allocation and leave no partial output.
- The final Git status and `C++` inventory/hash baseline match the initial isolation baseline.

**Independent verification:** At least two distinct read-only agents:

- scientific-equation and source-provenance audit;
- sampling, seeds, consumption accounting, statistics, and CLI/output audit.

GPT-5.6-sol resolves discrepancies and requests focused repairs before sign-off.

### Stage 15: Documentation and cleanup

**Objective:** Leave C++2 understandable, reproducible, and free of temporary artifacts.

**Files:**

- Finalize: `README.md`
- Finalize: `docs/MODEL_SPECIFICATION.md`
- Finalize: `docs/SAMPLING_DESIGN.md`
- Finalize: `docs/CLI_REFERENCE.md`
- Finalize: `docs/SEED_POLICY.md`
- Finalize: `docs/VALIDATION.md`
- Finalize: `analysis/README.md`

**Checks:**

- Documentation describes current behavior, not implementation history.
- Aggregate-demand-conditioned marginal 1D LHS is labelled as a proposed alternative design and is never called a joint, grouped, or publication-equivalent LHS.
- Publication citations and source filenames are included.
- Generated logs, temporary CSVs, debug binaries, and scratch scripts are removed.
- Only intended source, tests, documentation, and approved small reference fixtures remain.
- No commit or push unless the user explicitly requests it.

---

## 8. Comprehensive testing matrix

Every testable function must have direct unit coverage. Integration tests do not replace unit tests.

### Deterministic unit tests

- LHS stratum construction and shuffling.
- Open-interval RNG integer-to-double conversion and deterministic bounded-integer generation.
- Seed derivation known-answer values.
- Every inverse CDF.
- Category lookup and overrides.
- Draw-count planning.
- Cursor accounting.
- QER intermediate equations.
- ECAi inversion and clamping.
- Percentile interpolation.
- Mean, sample SD, SE, and Student-t CI.
- CSV serialization.
- CLI parsing and validation.

### Property and invariant tests

- One value per LHS stratum.
- Every generated value consumed once.
- No generated value consumed twice.
- Every cursor equals planned length at completion.
- Sum of infector counts equals infector-bank lengths.
- Zero-infector simulations consume nothing beyond the binomial draw.
- All distribution outputs satisfy support constraints.
- All ECAi outputs finite and nonnegative.
- Same seed reproducibility.
- Cross-stream and cross-replicate seed uniqueness.
- Distinct person-level consumption IDs and traversal-order-invariant assignment.
- Checked arithmetic and graceful pre-allocation resource rejection.

### Statistical reference tests

Statistical tests must use deterministic tolerances and known reference vectors. Avoid flaky assertions such as requiring a random sample mean to fall within a narrow range.

- Compare inverse CDF points with SciPy-generated fixtures.
- Compare binomial count frequencies with exact probability intervals where the 1D LHS implies deterministic count bounds.
- Compare percentile and CI calculations with frozen known-answer vectors.
- Compare raw-to-ceiling-5 conversion with exact boundary fixtures.
- Use histogram/correlation diagnostics as audit outputs, not proof of independence.

### Integration tests

- Full single replicate with synthetic parameters.
- Multiple replicates with distinct seeds.
- CLI to summary CSV.
- Invalid CLI creates no output.
- High-occupancy case.
- Zero-community-rate case: all zeros without parameter banks.
- Community-rate override.
- Mask override.
- Publication-default category cases.
- Custom-scenario metadata for each scientific override, including duration semantics.
- Deliberately excessive inputs rejected before allocation without partial output.

### Build variants

Where supported:

- Debug build with assertions.
- Release build.
- Address/undefined behavior sanitizers for focused small tests.
- Warnings enabled (`-Wall -Wextra -Wpedantic`, with justified exceptions documented).

---

## 9. Risks and required design caveats

1. **Aggregate-demand-conditioned marginal stratification is not a standard global joint LHS.** It stratifies each named one-dimensional parameter bank over exactly the values that are used. It does not separately stratify `ninf=k` groups, simulation-level parameter tuples, within-infector vectors, susceptible–infector pairings, or QER sums. This changes dependence structure and must not be called publication-equivalent without evidence.
2. **Pre-sampling `ninf` changes call ordering.** Model equations may remain identical while outputs differ from the current family-buffer implementation.
3. **Separate 1D LHS lists do not prove statistical independence.** Distinct seeds and permutations prevent accidental identity; they do not mathematically prove independence.
4. **Discrete binomial output repeats counts.** Stratification applies to the input uniform values, not to unique integer infector counts.
5. **Tail percentiles are sensitive to dependence.** Validation must separate algebraic correctness from sampling-architecture differences.
6. **Zero-infector shortcut changes consumption but not ECAi.** Tests must prove it returns the same zero while consuming no unrelated values.
7. **Exact-size banks may be large.** Publication defaults may be manageable, but arbitrary CLI overrides require checked sizes, a memory-budget policy, and graceful rejection before allocation.
8. **Confidence intervals characterize the mean replicate percentile estimator.** They are approximate randomized-LHS Monte Carlo intervals and must not be presented as population-percentile intervals, individual-ECAi variation, or model-input uncertainty.
9. **Cross-platform random engines may differ.** Reproducibility claims must state supported compiler/library assumptions or use a fully specified RNG and shuffle algorithm.
10. **Publication versus correspondence settings may differ.** Version 1 defaults follow the publication; later reproduction presets must be explicit alternatives.
11. **Ceiling adjustment discards information.** Preserve the raw L/s/person result and calculate the ceiling-to-5 value only as a separately labelled reporting transformation; never use it in SD, SE, or CI calculations.
12. **Person-specific draws are intentional.** Every infector and susceptible receives separate applicable sampled values, including separate PBR draws, even where source wording could be read as sharing an occupant value.

---

## 10. Open questions to resolve before implementation

These are gates, not permission to guess:

- Confirm the exact Jones percentile interpolation convention.
- Confirm the publication’s mask values for all healthcare categories, especially Patient and Resident, rather than using later correspondence adjustments.
- Confirm whether the default base seed is specified by the publication or should be a documented implementation default.
- Confirm whether a category must be mandatory or whether `Patient` should be the CLI default.
- Verify the publication one-hour default; retain `--duration-hours` as an explicitly labelled custom-scenario override whose target probability applies over that duration.
- Use a 95% confidence level by default unless the user overrides `--confidence`.
- Confirm whether per-replicate percentile rows should be written by default or only on request; this plan proposes only on request.
- Confirm the exact output directory naming.

GPT-5.6-sol must resolve these through the source or user instruction before the relevant task begins.

---

## 11. Final acceptance criteria

The first C++2 version is complete only when all of the following are verified:

- Pure C++ CLI implementation exists only under `C++2`.
- Existing `C++` has no modifications attributable to C++2 work.
- Jones et al. equations, distributions, category inputs, masks, and defaults are cited and independently verified.
- Default run configuration is 10,000 simulations, 500 replicates, and the 96th percentile.
- One LHS-derived binomial count exists per simulation.
- Zero-infector simulations return zero and consume no other random values.
- Every named susceptible and infector parameter uses an exact-size 1D LHS list.
- Every generated list value is consumed exactly once.
- Every active susceptible receives its own susceptible parameter values, every individual infector receives its own complete infector parameter set, and no value—including PBR—is shared between susceptible and infector(s) or between infectors.
- No parameter stream or replicate shares an accidental seed or sequence.
- Negative ECAi is clamped to zero.
- Mean, sample SD, SE, and confidence interval are correct.
- The raw percentile mean and raw uncertainty statistics are reported in L/s/person, together with the separately labelled `5 * ceil(raw_mean / 5)` result.
- Adjusted values never enter the SD, SE, or confidence-interval calculations.
- Command-line defaults and overrides are strict, documented, and tested.
- Every testable function has unit coverage.
- Every stage has passed a separate read-only verification agent.
- The frozen quantitative validation table is satisfied, and any publication difference greater than 5 L/s/person has blocked completion pending adjudication.
- All focused tests, full CTest, smoke tests, and approved end-to-end tests pass.
- Generated results and temporary artifacts do not clutter the repository.
- No commit or push occurs without explicit user authorization.

---

## 12. Agent orchestration and adjudication protocol

### Coordinating brain

GPT-5.6-sol owns:

- task sequencing;
- source-of-truth interpretation;
- acceptance criteria;
- actual diff inspection;
- comparison of agent reports with tool output;
- adjudication of disagreements;
- authorization to proceed to the next stage;
- final completion assessment.

### Implementation agents

Use one fresh GPT-5.6-luna or Claude Sonnet background agent per digestible stage. Give each agent:

- exact file scope;
- exact source pages/tables/equations already verified;
- failing tests to write first;
- exact build/test commands;
- prohibition on modifying `C++`;
- prohibition on commits/pushes;
- required structured report.

### Verification agents

After every implementation stage, use a different GPT-5.6-luna or Sonnet agent in read-only mode. Require:

- source/specification compliance;
- code-quality review;
- meaningful test assessment;
- execution of the focused test target;
- exact file/line evidence for defects;
- explicit PASS or FAIL.

A verifier must not merely restate the implementer’s report.

### Discrepancy handling

If reports disagree:

1. GPT-5.6-sol reads the actual source, diff, tests, and authoritative Jones source.
2. Separate scientific-source disagreement from implementation disagreement.
3. State the adjudicated conclusion and evidence.
4. Dispatch a narrowly scoped repair agent.
5. Dispatch a fresh independent verifier.
6. Do not proceed until the gate passes.

### Progress reporting

After each stage, report briefly:

- stage completed;
- files changed;
- tests actually run and results;
- verifier verdict;
- remaining blocker, if any.

Do not claim completion based on plans, agent self-reports, compilation alone, or partial test suites.
