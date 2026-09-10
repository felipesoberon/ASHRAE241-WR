# C++2 — Aggregate-Demand-Conditioned Marginal 1D-LHS ECAi Model

C++2 is a C++17 command-line implementation of the Jones et al. (2025) ASHRAE 241 ECAi calculation. It implements a proposed alternative sampling design: **aggregate-demand-conditioned marginal one-dimensional LHS**.

C++2 is not demonstrated to reproduce the publication sampling algorithm. Each complete named parameter bank is marginally stratified. Simulation tuples, within-infector parameter combinations, susceptible–infector pairings, multi-infector sums, and event-selected prefixes are not joint-stratified. The Jones publication describes an LHS matrix but does not specify its dimensions, seed, column mapping, dependence structure, or exact software behavior. C++2 also gives the susceptible and each infector separate person-specific PBR draws, whereas Jones states that occupants share an identical PBR. Results must not be described as joint-LHS, grouped-LHS, or publication-equivalent.

The implementation is isolated from the sibling `C++` tree: it does not modify or link against it, and output paths equal to or inside that tree are rejected. The recorded isolation baseline is [`docs/ISOLATION_BASELINE.sha256`](docs/ISOLATION_BASELINE.sha256).

## Requirements

- WSL with a C++17 compiler
- CMake 3.16 or newer
- C++ standard library only

The build enables `-Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off`.

## Build and test from Windows

Set `CPP2` to the local C++2 directory, then run these commands in WSL:

```bash
wsl.exe -e bash -lc 'cmake -S "$CPP2" -B "$CPP2/build"'
wsl.exe -e bash -lc 'cmake --build "$CPP2/build" -j2'
wsl.exe -e bash -lc 'ctest --test-dir "$CPP2/build" --output-on-failure'
```

The authoritative primary source file is `2025 - Jones et al. - Risk modeling for ASHRAE Standard 241-2023 - control of infectious aerosols.pdf`.

## Run

A category is required for every run. The sole exception is `--help`, which must be the only argument and performs no lookup or simulation.

```bash
"$CPP2/build/ecai_cpp2" --help
"$CPP2/build/ecai_cpp2" --category "Exam room" --simulations 17 --replicates 3
```

The second command is a bounded example, not the default study workload.

Defaults are `N=10000` simulations per replicate, `R=500` replicates, percentile `0.96`, confidence `0.95`, target probability `0.001`, duration `1` hour, and base seed `485`. Community infection rate, occupancy, volume, and mask efficiency come from the selected category. General categories use community rate `0.01` and mask efficiency `0`; the five Health Care categories use `0.03` and `0.30`.

The executable accepts exactly 17 flags. Their domains, provenance rules, category aliases, and supported combinations are defined in [`docs/CLI_REFERENCE.md`](docs/CLI_REFERENCE.md).

## Results and output

A successful run always prints one console-summary line. CSV output is opt-in. `--output-directory DIR` requests `DIR/summary.csv`; `--write-replicates` additionally requests `DIR/replicates.csv`. Literal paths can be selected with `--summary-file` and `--replicates-file`. Existing destinations are never overwritten.

The summary CSV has exactly 89 columns and one data row. The replicate CSV has exactly 19 columns and `R` data rows; it contains replicate percentile estimates and accounting metadata, not all `N` raw ECAi values.

Within each replicate the reported result is the selected percentile of the `N` ECAi values. Across replicates, the raw mean, sample SD, SE, minimum, maximum, and two-sided Student-t CI are calculated from the unrounded replicate percentile estimates. The CI is an approximate randomized-LHS Monte Carlo interval for the mean of that estimator, not a CI for the population percentile, individual ECAi, or total model uncertainty. `5 * ceil(raw_mean / 5)` is reported separately as the publication comparator and never enters the raw statistics.

## Validation status

Deterministic, invariant, regression, writer, and bounded end-to-end checks pass. The bounded Patient room comparator is `FAIL_GT_ONE_INCREMENT`; Museum is `NEAR_WITHIN_ONE_INCREMENT`. Overall publication reproduction is **NOT PASS**. The default `10000 × 500` workload has not been run and remains pending explicit authorization. See [`docs/VALIDATION.md`](docs/VALIDATION.md).

## Documentation and layout

- [`docs/MODEL_SPECIFICATION.md`](docs/MODEL_SPECIFICATION.md) — publication extraction, equations, tables, discrepancies, and implemented conventions
- [`docs/SAMPLING_DESIGN.md`](docs/SAMPLING_DESIGN.md) — RNG, LHS banks, seed inputs, assignment order, and resource checks
- [`docs/CLI_REFERENCE.md`](docs/CLI_REFERENCE.md) — complete current command-line and CSV contract
- [`docs/SEED_POLICY.md`](docs/SEED_POLICY.md) — canonical configuration fingerprint and stream-seed fold
- [`docs/VALIDATION.md`](docs/VALIDATION.md) — bounded validation evidence and limitations
- [`analysis/README.md`](analysis/README.md) — retained-analysis policy and focused validation command
- [`data/reference/publication_cases.csv`](data/reference/publication_cases.csv) — approved small bounded regression reference
- `include/`, `src/` — implementation
- `tests/` — unit, invariant, CLI, output, and end-to-end tests

`build/`, other `build*` directories, `Testing/`, and transient `results/` contents are generated artifacts and are not retained.
