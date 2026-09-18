# Bounded full-Jones vs. Jones-fitted comparison procedure

## Purpose and scope

This document defines a bounded procedure for comparing `--qer-distribution jones` (full
mechanistic Jones et al.) against `--qer-distribution jones-fitted` (the direct log10-normal fit
to Jones output) under matched controls. It does not cover Mikszewski profiles, does not run the
full default workload, and does not itself contain results: it is a predeclared procedure to be
executed and reported separately, per `docs/VALIDATION.md`'s authorization requirement for any
run at or near the full 10,000 × 500 default workload.

`jones-fitted` is an empirical lognormal approximation to full Jones, not an alternate encoding of
the same calculation. Any comparison between the two families is therefore a **statistical**
comparison of two different emission realizations under matched non-QER controls, never a
bit-identical replay. Full agreement, partial agreement, and disagreement are all possible
outcomes and none is assumed in advance.

## Why the comparison is statistical, not bit-identical

Per `docs/SEED_POLICY.md`, QER family is not a field of `EffectiveConfig` and is not folded into
the seed fingerprint. Both families derive an identical `infector_count` stream and identical
susceptible banks (`susceptible_pbr`, `gamma`, `lambda`) at matched controls, and the shared
`direct_qer` stream (`StreamId` 12) is itself seeded identically regardless of family. But full
Jones computes QER mechanistically from eight sampled infector-side inputs, while `jones-fitted`
transforms a single direct-QER LHS draw through a fitted `(mu_log10, sigma_log10)` pair. These are
two different sampling and transform mechanisms feeding the same downstream ECAi/probability
equations, so per-row and per-replicate outputs are not expected to match exactly even at a shared
base seed. Comparison must therefore rely on across-replicate summary statistics and their
confidence intervals, not on row-level or replicate-level exact equality.

## Matched controls

Every comparison run pair must hold the following fixed and identical between the `jones` and
`jones-fitted` invocation:

- occupancy category (`--category`);
- community infection rate (default: category value; do not override unless explicitly varying
  prevalence as a separate axis);
- mask efficiency (default: category value; do not override unless explicitly varying mask as a
  separate axis);
- percentile (`--percentile`, default `0.96`);
- simulations per replicate `N` (`--simulations`);
- replicate count `R` (`--replicates`);
- base seed (`--base-seed`);
- duration (`--duration-hours`), target probability (`--target-probability`), and confidence level
  (`--confidence`) — all three are fingerprinted (`docs/SEED_POLICY.md`) and must be identical
  across the pair or the two runs sample different streams entirely.

Only `--qer-distribution` differs between the paired invocations. Any deviation from this list
must be documented as an explicit, deliberate second axis of variation (e.g. a mask-on/mask-off
sub-comparison), never as an incidental difference.

### Seed and replicate controls

- Use one fixed base seed per comparison pair; do not vary the base seed within a pair.
- To assess sensitivity to the base seed itself, repeat the full matched pair at two or more
  additional base seeds and report each pair's statistics separately — do not pool across base
  seeds into a single interval.
- Replicate count `R` must be identical within a pair and should be large enough to produce a
  usable across-replicate confidence interval (`dispersion_status` is `UNAVAILABLE_SINGLE_REPLICATE`
  for `R=1`).

## Output modes

Run both output modes where applicable, as separate matched pairs:

- **ECAi mode** (`--mode ecai`, default): compares the raw across-replicate mean of the
  replicate-percentile ECAi estimator between families.
- **Probability mode** (`--mode probability`): compares the raw across-replicate mean of the
  replicate-percentile forward infection probability between families, at a fixed `--input-ecai`
  held identical across the pair (default: the category's published Table 4 ECAi, unless the
  comparison explicitly targets a different fixed input as its own documented variant).

Each mode's pair is independent; do not mix ECAi-mode and probability-mode outputs into one
statistic.

## Command placeholders

Substitute the category, N, R, base seed, and (for probability mode) input ECAi actually used.
`<temporary>` denotes a scratch output directory removed after data collection, per existing
project convention.

```text
ecai_cpp2 --category "<CATEGORY>" --qer-distribution jones \
  --simulations <N> --replicates <R> --base-seed <SEED> --output-directory <temporary>

ecai_cpp2 --category "<CATEGORY>" --qer-distribution jones-fitted \
  --simulations <N> --replicates <R> --base-seed <SEED> --output-directory <temporary>
```

Probability-mode variant of each command above:

```text
ecai_cpp2 --category "<CATEGORY>" --qer-distribution jones --mode probability \
  --input-ecai <ECAI> --simulations <N> --replicates <R> --base-seed <SEED> \
  --output-directory <temporary>

ecai_cpp2 --category "<CATEGORY>" --qer-distribution jones-fitted --mode probability \
  --input-ecai <ECAI> --simulations <N> --replicates <R> --base-seed <SEED> \
  --output-directory <temporary>
```

An explicit `--qer-distribution`, including an explicit `jones`, marks each run
`custom_scenario`; the existing full-Jones publication comparator is not applicable to either run
in a `jones-fitted` pair and must not be reported as if it were.

## Outputs and measurements to collect

For each of the four (or more, if additional base seeds or a mask/prevalence axis are added)
matched runs, collect from `summary.csv`:

- `raw_mean_l_s_person` (ECAi mode) or `raw_mean_probability` (probability mode);
- `raw_sd_l_s_person`/`raw_se_l_s_person` or their probability-mode equivalents;
- `raw_ci_lower_*`/`raw_ci_upper_*` at the run's `--confidence` level;
- `qer_family` (and Mikszewski-profile field, if ever extended to that family) provenance columns;
- `all_replicates_zero_infector_count`, `all_replicates_negative_clamp_count` (ECAi mode only —
  absent in probability-mode schema), and `all_replicates_output_zero_count`, to characterize
  zero-mass structure the same way `docs/VALIDATION.md` does for the existing bounded runs;
- `comparator_status` (expected `NOT_APPLICABLE_CUSTOM` for both runs in the pair, since
  `--qer-distribution` marks the run custom).

Retain each pair's `summary.csv` (and `replicates.csv` if `--write-replicates` is used) as the
frozen record of the comparison, analogous to `data/reference/publication_cases.csv` for the
existing bounded Jones runs.

## Predeclared interpretation

State the interpretation rule before running anything, and do not revise it after seeing results:

- The two families' across-replicate means are compared using their respective confidence
  intervals. Overlapping intervals are evidence consistent with statistical agreement at the
  chosen confidence level; non-overlapping intervals are evidence of a statistically detectable
  difference. Neither outcome is a claim of scientific correctness for either family — both
  remain implementation baselines, and `jones-fitted` remains an approximation to `jones`, not an
  independent validation of it.
- Do not claim agreement, disagreement, or model validity before the measurements above are
  actually collected and recorded. A predicted or assumed outcome must never be written into this
  document or `docs/VALIDATION.md` in place of a measured one.
- Any comparison run at or approaching the full 10,000 × 500 default workload requires the same
  explicit user authorization already required for the existing bounded Jones runs, and must be
  reported separately without replacing existing bounded reference rows.
- Results, once measured, belong in a dedicated results section or a follow-up document — not
  retrofitted into this procedure document, which defines the method only.
