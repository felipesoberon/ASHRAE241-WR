# Paired multi-target ECAi: implementation and verification plan

## Objective

Calculate ECAi for several target infection probabilities using the same sampled
scenario, while retaining the existing single-target behavior. Preserve each
replicate and simulation result so percentiles can be calculated later.

## Design

- Add an explicit opt-in list of unique, finite target probabilities in (0, 1).
  Reject ambiguous combinations with the existing single-target option and
  probability-output mode. Preserve the existing command and output contracts.
- For each sampled row, calculate QER and other scenario inputs once, then
  evaluate the existing ECAi equation for every target, including the existing
  zero-infector and negative-clamping rules. Support both infector policies and
  all supported QER families without changing their existing paths.
- Use the first listed target in the existing configuration fingerprint for
  paired sampling. List order therefore changes the sampled scenarios when it
  changes the first target. Record the original list order and retain all other
  existing fingerprint inputs. A one-target list uses the same seed as a
  corresponding single-target invocation.
- Stream each replicate's indexed results to a versioned, documented binary
  dataset with target order, units, category, policy, model version, seed,
  configuration and replicate metadata. Validate dimensions and resource use;
  avoid overwriting existing results and prevent incomplete files from being
  mistaken for completed datasets.
- Provide a reader/analysis path for arbitrary percentiles per target and
  replicate, using the existing interpolation convention. Keep replicate-level
  results distinct so Monte Carlo uncertainty is not confused with input variation.

## Verification

- Preserve existing single-target regression tests and output behavior.
- Check invalid, duplicate and reordered target lists; verify the sampling is
  unchanged when the first target is unchanged, and that changing the first
  target follows the existing fingerprint policy.
- Compare each paired target to the existing single-target equation under
  identical sampled inputs, including zero-infectors, clamping, both infector
  policies and QER paths. Single-target CLI invocations with different target
  probabilities are not row-wise comparators because the legacy fingerprint
  includes the target.
- Round-trip saved values through the reader; independently calculate several
  percentiles, including extremes and repeated/zero values. Exercise truncated
  files and output failure behavior.
- For a few representative occupancy categories and both infector policies,
  run 50 replicates of 10,000 scenarios, compare against the prior implementation,
  quantify any expected sampling differences, and retain verification evidence
  under the separate `ASHRAE241-WR-Data` directory. Do not run all categories.
- Benchmark time and storage on the intended machine before scaling up.
