# C++ — Aggregate-Demand-Conditioned Marginal 1D-LHS ECAi Model

This is the current C++ implementation (developed under the working name `C++2`, and still referred to by that name in some historical design documents and in the `ecai_cpp2` executable/target names). It is a C++17 command-line implementation of the Jones et al. (2025) ASHRAE 241 ECAi calculation, implementing a proposed alternative sampling design: **aggregate-demand-conditioned marginal one-dimensional LHS**.

This implementation uses the publication's model equations and parameters with a different sampling design. Each complete named parameter bank is marginally stratified. Simulation tuples, within-infector parameter combinations, susceptible–infector pairings, multi-infector sums, and event-selected prefixes are not joint-stratified. The Jones publication describes an LHS matrix but does not specify its dimensions, seed, column mapping, dependence structure, or exact software behavior, so this implementation's sampling design cannot be compared directly to the publication's. This implementation also gives the susceptible and each infector separate person-specific PBR draws, whereas Jones states that occupants share an identical PBR. Numerical results may differ from the publication depending on these sampling and input choices.

The implementation is isolated from the sibling `C++_legacy` tree (the original C++ port, formerly the sibling `C++` tree before the rename): it does not modify or link against it, and output paths equal to or inside that tree are rejected. The recorded isolation baseline is [`docs/ISOLATION_BASELINE.sha256`](docs/ISOLATION_BASELINE.sha256).

## Requirements

- WSL with a C++17 compiler
- CMake 3.16 or newer
- C++ standard library only

The build enables `-Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off`.

## Build and test from Windows

Set `CPP` to the local C++ directory, then run these commands in WSL:

```bash
wsl.exe -e bash -lc 'cmake -S "$CPP" -B "$CPP/build"'
wsl.exe -e bash -lc 'cmake --build "$CPP/build" -j2'
wsl.exe -e bash -lc 'ctest --test-dir "$CPP/build" --output-on-failure'
```

The authoritative primary source file is `2025 - Jones et al. - Risk modeling for ASHRAE Standard 241-2023 - control of infectious aerosols.pdf`.

## Run

A category is required for every run. The sole exception is `--help`, which must be the only argument and performs no lookup or simulation.

```bash
"$CPP/build/ecai_cpp2" --help
"$CPP/build/ecai_cpp2" --category "Exam room" --simulations 17 --replicates 3
```

The second command is a bounded example, not the default study workload.

Defaults are `N=10000` simulations per replicate, `R=500` replicates, percentile `0.96`, confidence `0.95`, target probability `0.001`, duration `1` hour, and base seed `485`. Community infection rate, occupancy, volume, and mask efficiency come from the selected category. General categories use community rate `0.01` and mask efficiency `0`; the five Health Care categories use `0.03` and `0.30`.

The executable's flags, domains, provenance rules, category aliases, and supported combinations are defined in [`docs/CLI_REFERENCE.md`](docs/CLI_REFERENCE.md).

## Results and output

A successful run always prints one console-summary line. CSV output is opt-in. `--output-directory DIR` requests `DIR/summary.csv`; `--write-replicates` additionally requests `DIR/replicates.csv`. Literal paths can be selected with `--summary-file` and `--replicates-file`. Existing destinations are never overwritten.

The summary CSV has exactly 89 columns and one data row. The replicate CSV has exactly 19 columns and `R` data rows; it contains replicate percentile estimates and accounting metadata, not all `N` raw ECAi values.

Within each replicate the reported result is the selected percentile of the `N` ECAi values. Across replicates, the raw mean, sample SD, SE, minimum, maximum, and two-sided Student-t CI are calculated from the unrounded replicate percentile estimates. The CI is an approximate randomized-LHS Monte Carlo interval for the mean of that estimator, not a CI for the population percentile, individual ECAi, or total model uncertainty. `5 * ceil(raw_mean / 5)` is reported separately as the publication comparator and never enters the raw statistics.

For an opt-in paired run, pass an ordered comma-separated list to
`--target-probabilities` and a new directory to `--raw-output-directory`.
The same sampled scenarios are used for every target; the first target supplies
the seed fingerprint. Each replicate's individual ECAi values are streamed to
`ecai.bin`, with a completed `manifest.txt`. For arbitrary per-replicate
percentiles, use `python3 analysis/read_paired.py DIRECTORY --percentiles
0.75,0.9,0.96`. See the CLI reference for the format and safety details.

## Validation status

Deterministic, invariant, regression, writer, and bounded end-to-end checks pass. For the bounded Patient room case, the computed ECAi differs from the publication's adjusted Table 4 value by more than one 5 L/s/person increment (status `FAIL_GT_ONE_INCREMENT`); for Museum, the difference is within one increment (status `NEAR_WITHIN_ONE_INCREMENT`). These statuses describe numerical distance from the published Table 4 value under this implementation's sampling design and inputs. The default `10000 × 500` workload has not been run and remains pending explicit authorization. See [`docs/VALIDATION.md`](docs/VALIDATION.md).

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
