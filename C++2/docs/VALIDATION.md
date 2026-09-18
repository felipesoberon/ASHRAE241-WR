# Bounded validation and publication comparison

## Status and scope

The results have separate statuses:

- implementation/invariant/regression gates: **PASS**
- Patient publication comparator: `FAIL_GT_ONE_INCREMENT`
- Museum publication comparator: `NEAR_WITHIN_ONE_INCREMENT`
- overall publication reproduction: **NOT PASS**

The passing gates cover deterministic algebra regression, executable wiring, metadata,
accounting, finite output, and same-seed reproducibility. They do not establish independent
scientific correctness, publication equivalence, calibration, or overall model validity. In
particular, Patient is not a pass and is not relabelled as one.

The required investigations are complete. They covered the equation, inputs, masks, percentile
convention, stochastic transforms, and sampling-design effects. They found no production bug and
did not justify changing the frozen alternative design. The Patient comparator remains explicitly failed,
and overall publication reproduction remains not passed. The implementation may advance without
changing either status.

The full default workload of 10,000 simulations × 500 replicates was **not run** and remains
pending explicit user authorization.

## Authoritative source and audit result

Primary source file:

`2025 - Jones et al. - Risk modeling for ASHRAE Standard 241-2023 - control of infectious aerosols.pdf`

B. Jones, C. Iddon, M. Zaatari, P. Wargocki, and R. Bruns, “Risk Modeling for
ASHRAE Standard 241-2023 – Control of Infectious Aerosols,” *Building and Environment* (2025),
doi:10.1016/j.buildenv.2025.113318.

The independent source audit covered equations on printed pp. 9–12, general inputs and Table 1
on printed pp. 14–19, category Tables 2–3 on printed pp. 20–21, and adjusted Table 4 comparators
on p. 22. Table 4 reports Patient room = 35 and Museum = 30 L/s/person. Only the Table 4 values
and source parameters are independent publication references. The raw C++2 means, confidence
intervals, aggregate counters, and other generated values are deterministic C++2 implementation
baselines for bounded regression; they are not independent scientific-correctness evidence.

Section 4 on printed p. 23 says category cumulative distribution functions are available in the
Supplementary Materials but were not used by the Project Committee. That supplementary file was
unavailable and was not inspected. Its absence limits exact publication-algorithm reconstruction;
it does not alter the printed Table 4 values.

## Architecture interpretation

C++2 uses **aggregate-demand-conditioned marginal one-dimensional LHS** banks. It samples infector counts,
then creates exact-size parameter banks for active susceptible rows and all realized infectors.
Every active susceptible and every individual infector receives person-specific separate PBR and
other sampled values. Zero-infector rows consume no other parameter values. This independently
specified architecture provides deterministic cursor and accounting invariants.

Jones et al. explicitly states on printed pp. 10 and 19 that susceptible and infected occupants
use explicit identical PBR. That person-specific-separate versus identical-PBR rule is the
definitive documented difference. On printed p. 12, the article describes a matrix of uniformly distributed
random numbers produced by Latin Hypercube sampling in MATLAB/Excel, but does not publish its
dimensions, seed, column mapping, dependence structure, or exact software behavior. Because the
matrix design is underspecified, broader LHS equivalence or non-equivalence cannot be established.
C++2 is independently specified and not demonstrated equivalent to the publication algorithm;
it is not claimed to be definitively broadly different from unknown publication details.

Equations and published marginal/source parameters can therefore be audited, but C++2 generated
results must not be called publication-equivalent. No tuning, conditioning, gating, or `I0 - 1`
change was made.

## Deterministic regression fixture and wiring checks

`tests/test_end_to_end.cpp` calls production pure model functions and checks hard-coded
double-precision regression fixtures using absolute `1e-12` or relative `1e-10`:

| Quantity | Frozen regression value |
| --- | ---: |
| V\*drop | 1.1309733552923255e-11 |
| GVL (copies/m³) | 1000000 |
| VER (virions/h) | 2.2619467105846509e-7 |
| QER (quanta/h) | 1.1309733552923255e-8 |
| ECAi fixture (L/s/person) | 1.01560284592605 |

These are deterministic algebra regression checks, not an independent scientific validation.
The end-to-end test also runs a small Exam room case through the executable and CSV writer,
checks the fixed 89-column summary and 19-column replicate schemas, exact row counts and
accounting flags, verifies custom-comparator suppression, and requires byte-identical output for
a same-seed repeat. Its output directories are temporary and removed by the test.

## Bounded publication-condition regression runs

Both bounded runs retained the Jones category source parameters and C++2 science defaults:
96th-percentile replicate estimator, 95% confidence level, target probability 0.001 over 1 hour,
category community rate, and category mask. They used base seed 485. Only replicate count was
bounded below the unauthorized default workload.

```text
ecai_cpp2 --category "Patient room" --simulations 10000 --replicates 10 --output-directory <temporary>
ecai_cpp2 --category Museum --simulations 10000 --replicates 5 --output-directory <temporary>
```

The CI in this section and in `publication_cases.csv` is an approximate randomized-LHS Monte
Carlo Student-t CI for the mean replicate-percentile estimator. It is not a CI for a population
percentile, individual variation, or total model uncertainty, and the publication supplies no
corresponding CI.

| Category | N | R | Mean replicate-p96 estimator | Approximate randomized-LHS Monte Carlo Student-t CI for the mean | Ceiling5 | Jones Table 4 adjusted | Signed / absolute difference | Comparator status |
| --- | ---: | ---: | ---: | --- | ---: | ---: | --- | --- |
| Patient room | 10000 | 10 | 0.66694130367595861 | [0.38942117462747522, 0.94446143272444205] | 5 | 35 | -30 / 30 | `FAIL_GT_ONE_INCREMENT` |
| Museum | 10000 | 5 | 22.437573428908024 | [21.653508104290935, 23.221638753525113] | 25 | 30 | -5 / 5 | `NEAR_WITHIN_ONE_INCREMENT` |

These generated values are frozen at maximum round-trip precision in
`data/reference/publication_cases.csv`. The Patient result remains a failed comparator. Museum is
near within exactly one reporting increment, not an exact reproduction.

## Reproducible aggregate accounting and zero mass

The following values are aggregate metadata from the same bounded runs shown above and are
retained in `data/reference/publication_cases.csv`. Their source columns are
`all_replicates_n`, `all_replicates_zero_infector_count`,
`all_replicates_negative_clamp_count`, and `all_replicates_output_zero_count` in the retained
89-column generated summaries. The focused test reruns Patient at N=10,000, R=10, seed 485 and
Museum at N=10,000, R=5, seed 485, then requires exact equality between every retained counter
and the newly generated summary.

| Category | All rows | Zero-infector count | Negative-clamp count | Zero-output count | Positive-output count |
| --- | ---: | ---: | ---: | ---: | ---: |
| Patient room | 100000 | 91269 (91.269%) | 4653 (4.653%) | 95922 (95.922%) | 4078 (4.078%) |
| Museum | 50000 | 896 (1.792%) | 43261 (86.522%) | 44157 (88.314%) | 5843 (11.686%) |

Every percentage is `100 × count / all_replicates_n`; positive-output count is
`all_replicates_n - all_replicates_output_zero_count`. Under the recorded short-circuit and clamp
policies, the test also requires
`all_replicates_output_zero_count = all_replicates_zero_infector_count +
all_replicates_negative_clamp_count` exactly.

For Patient, the retained zero-output point mass is 95.922%, just 0.078 percentage points below
the requested 96th percentile. Thus the p96 lies immediately above the observed zero point mass,
which makes that estimator sensitive to the sparse positive tail. This count-based observation
does not reconstruct the publication's unpublished sampling matrix or prove which unpublished
algorithm produced Table 4. Museum does not have the same near-96% zero-mass boundary.

The zero-infector counts are consistent with the run's binomial infector-count mechanism and the
retained source facts (`I0` and community rate). The aggregate counters support accounting and
quantile-boundary interpretation only; they do not turn the comparator into a scientific pass.

## Acceptance boundary

The focused test requires finite values, reconciled counters and banks, exhausted cursors,
canonical metadata, raw-statistics/ceiling formulas, the comparator policy including the known
Patient failure, same-seed reproducibility, exact bounded aggregate counters, and removal of
generated test output. These are deterministic regression, wiring, and accounting checks.
Passing them is compatible with a failed publication comparator and does not establish scientific
validity. Any future 10,000 × 500 run requires explicit authorization and must be reported
separately without replacing these bounded reference rows.

## Probability-mode verification

Probability outputs are bounded to [0,1]. A value of 1.0 is permitted numerical saturation, not
an error. Zero-infector outputs normalize to +0.0. Nonfinite intermediates are rejected, and
rejected invocations produce no output. Verification currently passes 12/12 CTest plus 133 output
checks.

## Risk-averse verification

Direct smoke runs verified both risk-averse ECAi and probability output paths: each exited 0, produced its versioned schema (`stage13-ra-v1` or `stage13-prob-ra-v1`), reported `NOT_APPLICABLE_RISK_AVERSE`, and produced no zero-infector rows. An invalid risk-averse `p=0` invocation exited nonzero and produced no files. The simulation-focused executable passed 292 checks and the full CTest suite passed 12/12. These checks verify implementation wiring and invariants, not independent scientific validation.

## QER distribution families

`--qer-distribution` (`docs/CLI_REFERENCE.md`) selects one of three QER calculation families,
resolved to an immutable descriptor once before sampling:

- `jones` (default) — the existing full mechanistic Jones et al. calculation, unchanged: three
  susceptible banks plus the eight infector-side QER banks, the current bank order, seeds,
  cursors, and equations.
- `jones-fitted` — a direct log10-normal fit to full-Jones output (`QER_i = 10 **
  (mu_log10 + sigma_log10 * normal_ppf(u))`), one draw per realized infector from a single
  direct-QER bank, using the 25-record parameter table carried over from the legacy `C++`
  implementation.
- a full canonical Mikszewski profile name (e.g. `mikszewski-sars-cov-2-standing-speaking`) —
  the same direct log10-normal mechanism using literal median-quanta/h and organism-specific
  `sigma_log10` values from the Mikszewski parameter table. Bare `mikszewski` is rejected.

Full Jones and both direct families apply the same susceptible-side inputs and the same ECAi/
probability equations; only the infector-side emission realization differs. Family selection is
orthogonal to output mode: both ECAi and probability modes dispatch through the same
family-resolved QER path (`docs/QER_DISTRIBUTION_IMPLEMENTATION_PLAN.md` Phase 5).

### Matched controls and stream sharing

QER family and Mikszewski profile are not fields of `EffectiveConfig` and are not folded into the
seed fingerprint (`docs/SEED_POLICY.md`). For an identical `(category, overrides, N, R,
percentile, confidence, target, duration, base seed)`, all three families derive an identical
`infector_count` stream and identical susceptible banks, and the shared `direct_qer` stream
(`StreamId` 12) is seeded identically across `jones-fitted` and every Mikszewski profile. This
means a bounded full-Jones-versus-Jones-fitted run at matched category, prevalence, mask,
percentile, N, and base seed compares two different emission realizations under otherwise-matched
sampling, not a bit-identical replay. See
`QER_FULL_JONES_VS_JONES_FITTED_VALIDATION.md` for the bounded comparison procedure. This is a
statistical comparison, not a claim of numerical agreement, and no such agreement is asserted
here.

An explicit `--qer-distribution`, including an explicit `jones`, marks the run `custom_scenario`.

### Schemas, provenance, and comparator suppression

Direct-family runs do not introduce a new summary/replicate schema version; they populate the
existing versioned schema's QER-family provenance columns (family, canonical Mikszewski profile
where applicable, parameter-table/version identifier, transform and units, direct-bank layout,
and direct values generated/consumed/evaluated) alongside the fields already documented above.
The existing full-Jones publication comparator (Table 4 adjusted ECAi) is not valid for
`jones-fitted` or Mikszewski output: alternative-family runs report the comparator as not
applicable, the same suppression already used for custom scenarios and risk-averse runs, rather
than being scored against a Jones-derived reference value.

### Interaction with risk-averse infector policy

`--infector-policy risk-averse` and `--qer-distribution` are independent selectors that compose:
risk-averse conditioning on infector count (`N_active == N`, no zero-infector rows) is unchanged
by QER family, and each realized infector still consumes exactly one direct-QER draw (direct
families) or the eight full-Jones draws (`jones`) regardless of infector policy. Risk-averse
direct-family runs report `NOT_APPLICABLE_RISK_AVERSE` for the population comparator, consistent
with community-policy direct-family runs reporting the QER-family comparator as not applicable.
