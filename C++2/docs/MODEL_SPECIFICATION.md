# C++2 Model Specification

Status: implemented. The current C++2 source is the implementation authority; Jones et al. (2025)
remains the scientific authority. Current sampling mechanics are specified in
[`SAMPLING_DESIGN.md`](SAMPLING_DESIGN.md), the interface in
[`CLI_REFERENCE.md`](CLI_REFERENCE.md), seed derivation in [`SEED_POLICY.md`](SEED_POLICY.md),
and bounded evidence in [`VALIDATION.md`](VALIDATION.md).

This document:

1. Extracts, with exact page/table/equation citations, everything Jones et al. (2025) states
   about the model, its inputs, and its outputs.
2. Records the engineering decisions the publication does not make (RNG, seeds, tolerances,
   rounding direction, resource limits, percentile convention, person-specific sampling). Every such decision is explicitly labelled
   **C++2-specific** and is not attributed to Jones et al.

Discrepancies between sources are recorded in Section 9. D2 (Dwelling unit `PBR_GM`) and D7
(category naming convention) are resolved engineering decisions. D1 and D5 are terminology-only
publication inconsistencies with no numeric impact.

---

## 1. Sources checked

### 1.1 Primary source (authoritative for science)

- **File:** `2025 - Jones et al. - Risk modeling for ASHRAE Standard 241-2023 - control of infectious aerosols.pdf`
- **Citation:** B. Jones, C. Iddon, M. Zaatari, P. Wargocki, R. Bruns, "Risk Modeling for
  ASHRAE Standard 241-2023 – Control of Infectious Aerosols," *Building and Environment*
  (2025), doi:10.1016/j.buildenv.2025.113318. Journal Pre-proof, revised manuscript, 40 PDF
  pages (unnumbered cover + printed pages 2–37).
- **Method of extraction:** every page was read from the PDF (pages 1–20, then 21–40 in full,
  plus pages 10, 11, 20, 21, 22, 23 re-rendered at 4× zoom and, where the source table is
  typeset in landscape/rotated orientation on a portrait page (Table 2), cropped and rotated
  programmatically with PyMuPDF + Pillow for legibility). All numeric table values below were
  read directly off these high-resolution renders, not inferred from the reference code.
- **Pages actually inspected (printed page numbers as they appear in the PDF footer):**
  cover, 2 (title/abstract), 3 (highlights), 3–4 (Nomenclature), 5–8 (§1 Introduction),
  6–8 (§2 Design parameters, items 1–10), 9 (§3 Modeling, §3.1, Eq. 1, Eq. 2),
  10 (Eq. 3, Eq. 4, discussion of `R`), 11 (Eq. 5, Eq. 6), 12 (Eq. 7, MATLAB/Excel note),
  13 (Figure 1, §3.2 start), 14 (Eq. 8, §3.3 General inputs), 15 (Table 1),
  16 (E, RTD discussion), 17 (GVL, VF, DK discussion), 18 (λ, γ discussion, §3.4 start),
  19 (§3.4 continued, PBR discussion), 20 (Table 2), 21 (§3.5 Setting ECAi, Table 3),
  22 (Table 4), 23 (§4 Equivalent clean airflow rates, §5 Discussion start),
  24 (§5.2 The model), 25 (§5.3 Inputs and outputs), 26 (§5.4 Limitations summary),
  27 (§6 Summary and conclusions, Acknowledgements, CRediT), 28–37 (References [1]–[60]),
  plus the trailing highlights/declaration-of-interest pages.
- **Citation page convention:** every page citation in this document is
  the **printed footer page number** (the number typeset at the bottom of the page), never the
  PDF viewer/reader page index. The two numbering systems are offset by exactly one throughout
  the body of the article: PDF viewer page `N` = printed footer page `N − 1`, valid for viewer
  pages 3–39 (footer pages 2–37 inclusive). Viewer pages 1–2 (the Elsevier/production cover
  sheet and the title/abstract/highlights page) carry no footer page number, and the trailing
  highlights/declaration-of-interest pages likewise carry none. This offset was confirmed
  directly against the PDF for this task (footer pages 2, 8–13, 15–22 each spot-checked against
  their corresponding viewer index) and is recorded here once so no individual citation below
  needs to restate it.
- **Supplementary material:** §4 (p. 23) states "Cumulative distribution functions for each
  occupant category are available in the Supplementary Materials, although they were not used
  by the Project Committee." This supplementary file was not present alongside the PDF in the
  working directories available to this task and was not inspected. Per the publication's own
  text it was not used to derive Table 4, so its absence does not affect any value extracted
  here; noted for completeness because published supplementary material is next in the recorded
  source hierarchy when available.

### 1.2 Secondary implementation references (read-only, comparison only, do not override §1.1)

- `C++/src/model.h`, `C++/src/model.cpp` — occupancy parameter table, `sample_parameters`,
  `QER`, `infection_probability`, `compute_ECAi`.
- `C++/src/random_manager.h`, `C++/src/random_manager.cpp` — inverse-CDF implementations and
  the existing per-distribution-family LHS buffer scheme.
- `C++/src/ecai.cpp` — production CLI driver: percentile interpolation, ceiling-to-5 rounding,
  CSV output.
- `C++/analysis/ecai_uncertainty_study.cpp` — replicate/bootstrap/order-statistic/reproduction
  statistical-uncertainty study; source of the existing Student-t CI and percentile-bootstrap
  conventions, and of the "reproduction mode" correspondence settings (which deliberately
  deviate from the Jones defaults — see §9.9).
- `Python/model.py`, `Python/randomManager.py` — original reference implementation; numerically
  identical to the C++ port for every value checked.

### 1.3 Experimental design note (architecture only, not a scientific source)

- `conditional-lhs-by-infector-count.md`
  — proposes pre-sampling `ninf` via a 1D LHS + binomial inverse CDF, then generating
  conditional parameter sets. Used only to inform the sampling-architecture discussion in
  §10; contributes no physical constants.

### 1.4 Governing plan

- `C++2/IMPLEMENTATION_PLAN.md` — provenance for several C++2-specific decisions
  items (ceiling-to-5 formula, SplitMix64 mixing, `R`-vs-`RTD` distinction, canonical
  configuration fingerprint instead of `std::hash`). Where this specification freezes such an
  item, the plan section is cited alongside the (non-)citation from Jones et al.

---

## 2. Nomenclature (verbatim from PDF printed pp. 3–4; PDF viewer page index is one greater — see the citation page convention note in §1.1)

| Symbol | Meaning | Unit |
|---|---|---|
| γ | Surface deposition rate | h⁻¹ |
| λ | Biological decay | h⁻¹ |
| μ | Face mask efficiency | — |
| d̄ | Mean diameter of dehydrated respiratory particles | m |
| φ | Total removal rate | h⁻¹ |
| C_drop | Respiratory particles per unit volume of exhaled air | # m⁻³ |
| D | Duration of exposure | h |
| DK | Dose constant | — |
| E | Evaporation coefficient | — |
| GVL_i | Genomic viral load | RNA copies m⁻³ (as used in Eq. 4; see the `GVL_i` units note below) |
| i | Individual person in a space | — |
| I0 | Number of people in a space | — |
| P(I) | Probability of a susceptible person becoming infected | — |
| PBR | Pulmonary breathing rate | m³ h⁻¹ |
| Q | Quanta of pathogen deposited in the respiratory tract of an exposed person | — (quanta are a defined dose-count unit, not independently dimensioned; see §10) |
| QER_i | Quanta emission rate of a person | quanta h⁻¹ per person |
| R | Fractional reduction of viability after emission | — |
| RTD | Respiratory tract deposition fraction | — |
| T_ECAi | Total equivalent clean airflow rate | m³ h⁻¹ |
| V*_drop | Expelled volume ratio of hydrated respiratory particles to air | — |
| VER_i | Virion emission rate | virions h⁻¹ per person |
| VF | Viable fraction | — |
| VOL | Space volume | m³ |
| ECAi | Equivalent clean airflow rate per person | L s⁻¹ per person |
| IRMM | Infection Risk Management Mode | — |
| RR | Relative risk | — |
| WR | Wells-Riley | — |

**Note on `RTD` naming (PDF p. 4 vs. p. 15):** the Nomenclature (p. 4) calls `RTD` the
"Respiratory tract *deposition* fraction," but Table 1 (p. 15) labels the same symbol/row the
"Respiratory tract *absorption* fraction." Same symbol, same [0.43, 0.65] uniform distribution,
same citation ([30], Darquenne). This is an internal terminology inconsistency in the
publication (deposition vs. absorption), not a numeric discrepancy; recorded as D5 (§9).

**Note on `GVL_i` units (Nomenclature p. 4 vs. Table 1 p. 15):** the Nomenclature gives the
unit `GVL_i` has *when used in Eq. 4* (RNA copies m⁻³ — the post-conversion quantity). Table 1
gives the unit it is *sampled in* (log₁₀ RNA copies mL⁻¹, Normal(μ=7.0, σ=1.4)) together with a
footnote converting it to RNA copies m⁻³ by multiplying by 10⁶. These are consistent, not
contradictory, once sequenced correctly (sample in log₁₀ copies/mL, exponentiate, then ×1e6);
see §10 for the full algebraic check. Recorded here only because the two mentions could
otherwise be mistaken for a unit conflict.

---

## 3. Equations (verbatim, with page citations)

All equation numbers match the publication exactly.

**Eq. 1 (p. 9) — Wells-Riley absolute risk:**
`P(I) = 1 − exp(−Q)`

**Eq. 2 (p. 9) — Quanta dose:**
`Q = [PBR · D · (1 − μ)² / (φ · VOL)] · Σ_{i=1}^{I0} QER_i`
(`QER_i` is presumed zero for uninfected occupants; the sum is effectively over the realized
infectors only.)

**Eq. 3 (p. 10) — Quanta emission rate:**
`QER_i ≡ RTD · VER_i / DK`

**Eq. 4 (p. 10) — Virion emission rate:**
`VER_i = PBR · V*_drop · GVL_i · VF · RD`

The publication's typeset equation literally prints the trailing factor as **`RD`**, but the
prose immediately preceding and following Eq. 4 discusses only the symbol **`R`** ("`R` is the
fraction of viability lost due to rapid decay soon after the initial exhalation... we use a
value of `R = 1` hereon"), and the Nomenclature (p. 4) and Table 1 (p. 15, "Viability loss
fraction after exhalation, `R`" = 1.0, ref. [22]) define only `R`, never `RD`. This is recorded
as Discrepancy D1 (§9.1); C++2 uses the symbol `R` (constant, 1.0) in Eq. 4, and that constant
must never be confused with `RTD` (Eq. 3's respiratory tract deposition/absorption fraction),
per the recorded `R`-versus-`RTD` engineering decision.

**Eq. 5 (p. 11) — Expelled volume ratio:**
`V*_drop = (π / 6) · (d̄ · E)³ · C_drop`

**Eq. 6 (p. 11) — Total removal rate:**
`φ = γ + λ + T_ECAi / VOL`

**Eq. 7 (p. 12) — Equivalent clean airflow rate per person:**
`ECAi = T_ECAi / (I0 × 3.6)`

**Frozen ECAi inversion used by C++2 (algebraic rearrangement of Eqs. 1, 2, 6, and 7; not a separately numbered Jones equation):**

```text
Q_target       = -log(1 - target_probability)
mask_factor    = (1 - μ)^2
QER_sum        = sum(QER_i) over the realized infectors
term1_m3_h     = PBR_s * D * mask_factor * QER_sum / Q_target
term2_m3_h     = VOL * (γ + λ)
T_ECAi_raw     = term1_m3_h - term2_m3_h
ECAi_raw       = T_ECAi_raw / (3.6 * I0)
ECAi_clamped   = max(ECAi_raw, 0)
```

`Q_target` is dimensionless; `mask_factor` is dimensionless; both `term1_m3_h` and
`term2_m3_h` are m³/h; `ECAi_raw` and `ECAi_clamped` are L/s/person. The zero-infector
orchestration shortcut returns `ECAi_clamped = 0` before requesting `PBR_s`, `γ`, `λ`, or any
infector parameter. This shortcut is a C++2 sampling decision; it is algebraically consistent
with `QER_sum = 0` followed by clamping.

**Eq. 8 (p. 14) — Relative risk (out of scope for C++2 v1 output; documented for provenance
only — Section 1 of `IMPLEMENTATION_PLAN.md` does not list an RR output):**
`RR = (φ1 · I0,2) / (φ2 · I0,1)`
(subscript 1 = baseline/starting conditions, subscript 2 = intervention.)

### 3.1 Forward probability evaluation used by probability mode (C++2-specific, algebraic mirror of §3 above)

`ecai_cpp2` supports a second output mode, selected with `--mode probability` (default remains
`--mode ecai`, byte-identical to the mode's absence). Instead of inverting Eqs. 1, 2, 6 and 7 for
the ECAi that reproduces a target probability, probability mode takes a fixed **input** ECAi
(the selected category's published Table 4 value, §8 below, unless `--input-ecai` supplies an
explicit one) and evaluates Eq. 1 **forward** for every sampled simulation row. This is a C++2
engineering addition, not a separate Jones et al. equation; it is the algebraic mirror of the
frozen inversion already documented above, using the same Eqs. 1, 2, 6 and 7.

```text
ECA_total  = 3.6 * I0 * ECAi_in                            [m^3/h]   (Eq. 7, forward)
Lambda     = ECA_total + VOL * (gamma + lambda)            [m^3/h]   (Eq. 6, total removal)
dose_num   = PBR_s * D * (1 - mu)^2 * QER_sum               [m^3 quanta / h]
Q          = dose_num / Lambda                              [quanta, dimensionless]
P_raw      = -expm1(-Q)                                     (Eq. 1, forward)
P          = (P_raw > 0) ? P_raw : +0.0                     (sign normalization)
```

`-expm1(-Q)` is used instead of `1 - exp(-Q)` to retain relative precision for the small `Q`
(~1e-3) the default target implies. The forward evaluation is deliberately decoupled from
`target_probability`: it takes no target, validates none, and reads none. The target is used only
as a comparator at the output stage (below), exactly as Table 4's published ECAi is a comparator
for the ECAi-mode raw mean.

**Row rules, implemented exactly as in ECAi mode's zero-infector shortcut:**

- A zero-infector row is set to `P = +0.0` exactly before any bank value is requested; no
  equation runs. This mirrors the ECAi-mode shortcut (§11.1) and consumes no other sampled value.
- `Lambda > 0` always, given a validated `ECAi_in >= 0`, `VOL > 0` and `gamma + lambda > 0`, so
  there is no division by zero and no clamp is applied; `negative_clamp_count` is structurally
  zero in this mode and is not reported.
- **Range:** `Q >= 0` and finite implies `P` lies in the **closed** interval `[0, 1]`. `P == 1.0`
  exactly is reachable: IEEE-754 double `-expm1(-Q)` rounds to `1.0` once `exp(-Q) < 2^-53`
  (`Q` beyond roughly `36.7`). This is **numerical saturation** of the forward formula in double
  precision, not a model claim that infection is certain, and not an error. It is recorded per
  row as `saturated = (P == 1.0)` and counted per replicate (`saturated_count`).
- Every intermediate (`ECA_total`, `Lambda`, `dose_num`, `Q`, `P_raw`) is checked finite before
  use; a non-finite intermediate (for example an overflowing `ECA_total` from a pathological
  `--input-ecai`) is rejected with a diagnostic naming the intermediate rather than silently
  propagated. `Q >= 0` is additionally asserted; a negative `Q` is impossible after input
  validation and would indicate a defect, so it is reported rather than clamped.
- `P == 0` for an active row occurs only when `dose_num == 0` (mask efficiency `mu == 1`, or an
  underflowed `QER_sum`); `compute_qer` guarantees `QER_i > 0` for in-domain inputs, so this is
  the only active-row zero source. The reported `output_zero_count` therefore equals
  `zero_infector_count` plus these rows.

**Target comparator (probability mode):** the quantity compared against `target_probability` is
`raw_mean_probability`, the across-replicate mean of each replicate's within-replicate percentile
of per-simulation probabilities (default percentile `0.96`, same Type 7 linear convention as
ECAi mode, §12). There is no ceiling, no rounding, and no 5 L/s/person grid in this comparison.
The comparator is eligible only for an unmodified published-category run (`comparator_eligible =
!custom_scenario`); any of `--target-probability`, `--duration-hours`, `--percentile`,
`--community-rate`, `--occupants`, `--volume-m3`, `--mask-efficiency`, or `--input-ecai` makes the
run custom and suppresses it. `comparator_ratio = mean / target` (target is always in `(0, 1)`,
so never a division by zero); status is `AT_TARGET_WITHIN_10PCT` for `0.9 <= ratio <= 1.1`
(closed, evaluated with a few ULPs of rounding tolerance at the boundary), `BELOW_TARGET` for a
smaller ratio, `ABOVE_TARGET` for a larger one. This 10% band is a self-consistency tolerance
against the C++2 sample, not a standard's acceptance criterion, and the comparator does not
validate either implementation (§10 of `VALIDATION.md` and the caveats below).

**Zero-infector fraction and the percentile boundary:** zero-infector rows enter each replicate's
percentile as `P = +0.0`, exactly as ECAi mode's zero-infector rows enter it as `0`. When a
replicate's zero-infector fraction is at or above the requested percentile, that replicate's
percentile is exactly `0`; a mean of `0` across replicates yields `BELOW_TARGET`, which is the
correct reading of the sample, not a defect. Categories with a large zero-infector mass (Patient
room, per the bounded ECAi-mode run in `VALIDATION.md`) are the most sensitive to this boundary.
Probability mode has no negative-clamping step, so its zero-output fraction for a given sample is
not directly comparable to ECAi mode's zero-output fraction for the same category.

**Seed and sampling preservation:** mode and the resolved input ECAi are deliberately **not**
fingerprint fields (see `SEED_POLICY.md`). For identical `(category, overrides, N, R, percentile,
confidence, target, duration, base seed)`, ECAi mode and probability mode derive identical stream
seeds, identical infector counts, identical parameter banks, and identical consumption order —
only the deterministic per-row equation differs. This is intentional: it is what lets a
probability-mode run and an ECAi-mode run at the same settings be compared row-for-row.

**Scientific caveats specific to probability mode** (in addition to every caveat this document
already records for ECAi mode, §9–§11, all of which probability mode inherits unchanged, since it
reuses the same sampling architecture, distributions, and person-specific-PBR deviation D4):

- The reported quantity is the `p`-th percentile (0.96 by default) of per-simulation infection
  probability at a fixed input ECAi, evaluated over the effective duration `D`. It is not `1 - p`
  of anything, and it is not the probability that ASHRAE 241's target is met in general.
- Comparing the mean replicate percentile with `target_probability` is a self-consistency
  diagnostic of the C++2 sample at the published ECAi, not independent validation: Table 4 values
  are themselves ceilinged to multiples of 5 and were produced by an unpublished sampling matrix
  (§8 below), so agreement or disagreement with the frozen default target `0.001` does not
  validate or invalidate either implementation.
- `P == 1.0` is numerical saturation of `-expm1(-Q)` in double precision, not a statement that
  infection is certain; it is reachable only at very small input ECAi or extreme sampled
  emission.
- Duration `D` enters `Q` linearly; the frozen default target is labelled hourly only when the
  effective duration is exactly 1 h, exactly as in ECAi mode.
- `--occupants` overrides rescale total ECA (the per-person input ECAi times the effective `I0`)
  and, like every other scientific override, suppress the comparator.

### 3.2 Direct-QER families (C++2-specific extension, not a Jones et al. equation)

`ecai_cpp2` supports three QER (quanta emission rate) calculation families, selected with
`--qer-distribution` (`docs/CLI_REFERENCE.md`). This is a C++2 engineering extension layered on
top of the Jones et al. framework in §3 above; Jones et al. defines only the mechanistic
calculation (family 1 below) and does not define, endorse, or discuss families 2 or 3.

1. **Full Jones (default, family `jones`).** QER_i is computed mechanistically via Eq. 3
   (`QER_i = RTD · VER_i / DK`) from eight independently sampled infector-side inputs (`PBR_i`,
   `C_drop`, `d̄`, `E`, `GVL_i`, `VF`, `RTD`, `DK`; §5–§7, §11.1). This is the calculation
   documented everywhere else in this specification and is unchanged by the existence of the
   other two families.
2. **Jones-fitted (family `jones-fitted`).** One direct log10-normal distribution per canonical
   occupancy category, fitted to full-Jones QER output, replaces the eight-input mechanistic
   calculation for the infector side only:

   ```text
   QER_i = 10 ** (mu_log10 + sigma_log10 * normal_ppf(u))
   ```

   `QER_i` is in quanta/h per infector, one draw per realized infector. The table has exactly 25
   records, one per canonical publication category (§6 D7 naming), with `(mu_log10, sigma_log10)`
   copied exactly from the legacy `C++/src/qer_fitted.cpp` `qer_fitted_params` map
   (table version `qer_direct::kJonesFittedTableVersion`). Thirteen of the 25 legacy short names
   differ from their canonical category name (`Food`, `Transportation`, `Lecture`, `Sorting`,
   `Exam`, `Group`, `Patient`, `Resident`, `Waiting`, `Place`, `Spectator`, `Common`, `Dwelling` →
   e.g. `Food` → `Food and beverage facilities`, `Transportation` → `Transportation waiting`,
   `Sorting` → `Sorting, packing, light assembly`); the remaining 12 legacy short names are
   identical to their canonical category name. The complete, explicit mapping is
   `qer_direct::jones_fitted_legacy_aliases()`
   (`include/qer_direct.h`), independent of the CLI category-alias table in §6 D7.
3. **Mikszewski (family `mikszewski-<organism>-<activity>`).** 30 published profiles — 10
   organisms (SARS-CoV-1, MERS, TB on treatment, influenza, coxsackievirus, rhinovirus,
   SARS-CoV-2, TB untreated, adenovirus, measles) crossed with 3 activities (resting,
   standing-speaking, light-speaking-loudly) — each store a literal published median quanta/h and
   an organism-specific `sigma_log10` (Mikszewski et al. 2022, Tables 1–2, via the legacy
   `C++/src/qer_fitted.cpp` `mikszewski_params` map; table version
   `qer_direct::kMikszewskiTableVersion`). `mu_log10 = log10(median)` is computed at
   initialization without rounding. Bare `mikszewski` is rejected; a full canonical profile name is
   required (`docs/CLI_REFERENCE.md`).

Both direct families replace only the infector-side QER draw. Every susceptible-side input
(`PBR_s`, `γ`, `λ`), Eq. 1/2/6/7 (ECAi mode, §3 above), and the forward probability evaluation
(§3.1 above) are unchanged and apply identically to all three families and to both output modes.
A direct-family QER value is consumed exactly once per realized infector, in place of the
eight-input Eq. 3/Eq. 4 mechanism, not in addition to it (`docs/QER_DISTRIBUTION_IMPLEMENTATION_
PLAN.md` Phase 2–3). Bank layout and resource accounting for the three families are specified in
[`SAMPLING_DESIGN.md`](SAMPLING_DESIGN.md); output schema and comparator handling are in
[`CLI_REFERENCE.md`](CLI_REFERENCE.md) and `VALIDATION.md`. Selecting a direct family is not a
claim that its output is scientifically equivalent to full Jones; see
`QER_FULL_JONES_VS_JONES_FITTED_VALIDATION.md` for the bounded statistical comparison procedure.

---

## 4. Design parameters (verbatim summary, §2, pp. 6–8)

The Project Committee framed the model with ten numbered design parameters:

1. Occupancy categories already meet baseline IAQ under ASHRAE 62.1/62.2/170; Standard 241
   requirements are supplementary.
2. The standard applies only during Infection Risk Management Mode (IRMM), triggered by public
   health authorities; the standard itself does not define the trigger.
3. Only indoor **long-range** transmission (>1 m from an infector) is modeled; short-range,
   fomite, and close/direct contact are explicitly out of scope.
4. SARS-CoV-2 clinical data (not outbreak data) are used because it is the best-available,
   highest-quality dataset, and clinical data avoid the high-emitter bias of outbreak data.
5. An adapted Wells-Riley model (Wells; Riley et al.) is used.
6. Population-level risk accounts for the probability an infector is present via community
   incidence and occupant count; this is not a conditional/individual risk model and is
   explicitly stated to over-estimate risk in most spaces (§5.2).
7. Risk is evaluated **per hour** and is the same across all occupancy categories; the
   hourly transmission rate is assumed low enough that the linear form of Wells-Riley applies
   (§3.2 discusses the nonlinear/relative-risk simplification separately).
8. All uninfected people are assumed fully susceptible (no immunity) — pandemic-start
   assumption.
9. All respiratory particles are assumed instantaneously and uniformly dispersed (well-mixed
   space; no air-distribution-pattern modeling).
10. Uncertainties are estimated probabilistically (LHS-based Monte Carlo) "to make the results
    more robust."

---

## 5. Table 1 — Inputs common to all occupancy categories (p. 15, verbatim)

### Constants

| Parameter | Value | Source |
|---|---|---|
| Viability loss fraction after exhalation, `R` | 1.0 | [22] |
| Exposure duration, `D` (h) | 1.0 | (uncited — design choice, §2 item 7) |
| Community infection rate (%) — Healthcare | 3.0 | [28] |
| Community infection rate (%) — All other | 1.0 | [28] |
| Probability of infection, `P(I)` (%) | 0.10 | (uncited — design choice) |

### Variables (all probabilistic; "min/max" are hard bounds of the distribution's support)

| Parameter | Distribution | Parameters | Source |
|---|---|---|---|
| Respiratory particle evaporation factor, `E` | Beta | α=5.0, β=2.0, min=2.0, max=5.0 | [15] |
| Genomic viral load, `GVL_i` (log₁₀ RNA copies mL⁻¹) | Normal | μ=7.0, σ=1.4 | [29] |
| Viable fraction, `VF` | Beta | α=5.0, β=2.0, min=1e−4, max=1e−2 | [19] |
| Respiratory tract absorption/deposition fraction, `RTD` | Uniform | min=0.43, max=0.65 | [30] |
| Dose constant, `DK` | Uniform | min=5.0, max=15 | [19] |
| Biological decay rate, `λ` (h⁻¹) | Log-normal | GM=0.52, GSD=1.9 | [31] |
| Surface deposition rate, `γ` (h⁻¹) | Uniform | min=0.42, max=0.61 | [32] |

Footnote (verbatim): "`GVL_i` is converted into RNA copies per m³ by multiplying by 10⁶." GM =
geometric mean; GSD = geometric standard deviation; α, β = Beta shape parameters.

**Cross-check against `C++/src/random_manager.cpp` + `C++/src/model.cpp`:** every distribution
family, shape parameter, and bound above matches the existing implementation exactly
(`beta(5,2,loc=2,scale=3)` for `E`; `beta(5,2,loc=1e-4,scale=1e-2-1e-4)` for `VF`;
`uniform(0.43,0.65)` for `RTD`; `uniform(5,15)` for `DK`; `lognormal(ln 0.52, ln 1.9)` for `λ`;
`uniform(0.42,0.61)` for `γ`; `log10normal(7.0,1.4)` for `GVL`, then `×1e6` to convert
mL⁻¹→m⁻³). No discrepancy found for Table 1.

---

## 6. Table 2 — Occupancy category geometry and baseline airflow (p. 20, verbatim)

All 25 categories, in publication order, exactly as printed (landscape table, extracted via
rotated high-resolution render):

| Occupancy Category | Floor area (m²) | Space volume, VOL (m³) | Default density (people/100 m²) | Max occupants, I0 | P(susceptible present) | People outdoor airflow (L s⁻¹/person) | Area outdoor airflow (L s⁻¹ m⁻²) |
|---|---:|---:|---:|---:|---:|---:|---:|
| **Correctional Facilities** | | | | | | | |
| Cell | 120 | 320 | 25 | 30 | 0.25 | 2.5 | 0.6 |
| Dayroom | 100 | 600 | 30 | 30 | 0.25 | 2.5 | 0.3 |
| **Commercial/Retail** | | | | | | | |
| Food and beverage facilities | 100 | 360 | 50 | 50 | 0.38 | 3.8 | 0.9 |
| Gym | 500 | 1900 | 36 | 180 | 0.83 | 10 | 0.3 |
| Office | 1000 | 2700 | 5.0 | 50 | 0.38 | 2.5 | 0.3 |
| Retail | 1000 | 6000 | 15 | 150 | 0.77 | 3.8 | 0.6 |
| Transportation waiting | 1000 | 10000 | 10 | 100 | 0.62 | 3.8 | 0.3 |
| **Educational Facilities** | | | | | | | |
| Classroom | 120 | 320 | 25 | 30 | 0.25 | 5.0 | 0.6 |
| Lecture hall | 200 | 1800 | 75 | 150 | 0.77 | 2.5 | 0.3 |
| **Industrial** | | | | | | | |
| Manufacturing | 1000 | 12000 | 7.0 | 70 | 0.50 | 5.0 | 0.3 |
| Sorting, packing, light assembly | 400 | 4800 | 5.0 | 20 | 0.17 | 5.0 | 0.3 |
| Warehouse | 100 | 1200 | 20 | 20 | 0.17 | 5.0 | 0.3 |
| **Health Care** | | | | | | | |
| Exam room | 15 | 41 | 20 | 3 | 0.057 | 7.5 | 0.0 |
| Group treatment area | 100 | 270 | 20 | 20 | 0.43 | 7.5 | 0.0 |
| Patient room | 30 | 81 | 10 | 3 | 0.057 | 15 | 0.0 |
| Resident room | 30 | 81 | 10 | 3 | 0.057 | 15 | 0.0 |
| Waiting room | 100 | 270 | 30 | 30 | 0.57 | 5.0 | 0.0 |
| **Public Assembly/Sports/Entertainment** | | | | | | | |
| Auditorium | 500 | 4600 | 30 | 150 | 0.77 | 2.5 | 0.3 |
| Place of religious worship | 200 | 2400 | 90 | 180 | 0.83 | 2.5 | 0.3 |
| Museum | 1000 | 10000 | 40 | 400 | 0.97 | 3.8 | 0.3 |
| Convention | 1000 | 10000 | 40 | 400 | 0.97 | 3.8 | 0.3 |
| Spectator area | 300 | 6300 | 33 | 100 | 0.62 | 3.8 | 0.3 |
| Lobbies | 500 | 4600 | 30 | 150 | 0.77 | 2.5 | 0.3 |
| **Residential** | | | | | | | |
| Common space | 500 | 4600 | 30 | 150 | 0.77 | 2.5 | 0.3 |
| Dwelling unit | 200 | 540 | 3.0 | 6 | 0.049 | 6.8 | 0.0 |

**Cross-check:** `VOL` and `I0` for all 25 rows match `occupancy_params` in
`C++/src/model.cpp` / `Python/model.py` exactly, including the Museum/Convention maximum
occupancy of 400. No discrepancy found for Table 2's `VOL`/`I0` columns. The
"probability of susceptible people present," "people outdoor airflow rate," and "area outdoor
airflow rate" columns are not currently consumed by the reference C++/Python model code
(they feed Table 4's baseline-airflow column and the standard's own comparison, not the ECAi
Monte Carlo itself) but are recorded here for completeness and provenance.

---

## 7. Table 3 — Occupant-specific distributions by category (p. 21, verbatim)

All distributions are log-normal (breathing rate, `C_drop`, `d̄`); "Respiratory activity" is a
fixed breathing:talking:vocalizing time-weighting (%), not itself sampled — it was already
folded into the GM/GSD values below by the publication's authors (per §3.4, p. 19).

| Occupancy Category | PBR GM (m³ h⁻¹) | PBR GSD | Activity breathing:talking:vocalizing (%) | C_drop GM (×10⁴ m⁻³) | C_drop GSD | d̄ GM (×10⁻⁶ m) | d̄ GSD |
|---|---:|---:|---|---:|---:|---:|---:|
| **Correctional Facilities** | | | | | | | |
| Cell | 0.55 | 1.2 | 80:20:0 | 14 | 1.1 | 1.9 | 1.1 |
| Dayroom | 0.55 | 1.2 | 75:25:0 | 19 | 1.1 | 1.9 | 1.1 |
| **Commercial/Retail** | | | | | | | |
| Food and beverage facilities | 0.55 | 1.2 | 75:25:0 | 15 | 1.1 | 1.9 | 1.1 |
| Gym | 0.62 | 1.3 | 100:0:0 | 9.8 | 1.1 | 1.8 | 1.1 |
| Office | 0.55 | 1.2 | 80:20:0 | 14 | 1.1 | 1.9 | 1.1 |
| Retail | 0.55 | 1.2 | 95:5:0 | 11 | 1.1 | 1.8 | 1.1 |
| Transportation waiting | 0.64 | 1.2 | 75:25:0 | 15 | 1.1 | 1.9 | 1.1 |
| **Educational Facilities** | | | | | | | |
| Classroom | 0.55 | 1.2 | 70:30:0 | 17 | 1.1 | 1.9 | 1.1 |
| Lecture hall | 0.55 | 1.2 | 100:0:0 | 9.8 | 1.1 | 1.8 | 1.1 |
| **Industrial** | | | | | | | |
| Manufacturing | 0.73 | 1.3 | 75:25:0 | 15 | 1.1 | 1.9 | 1.1 |
| Sorting, packing, light assembly | 0.77 | 1.3 | 75:25:0 | 19 | 1.1 | 1.9 | 1.1 |
| Warehouse | 0.77 | 1.3 | 75:25:0 | 15 | 1.1 | 1.9 | 1.1 |
| **Health Care** | | | | | | | |
| Exam room | 0.62 | 1.3 | 100:0:0 | 32 | 1.1 | 2.0 | 1.1 |
| Group treatment area | 0.55 | 1.2 | 90:10:0 | 12 | 1.1 | 1.8 | 1.1 |
| Patient room | 0.62 | 1.3 | 100:0:0 | 30 | 1.1 | 2.0 | 1.1 |
| Resident room | 0.62 | 1.3 | 100:0:0 | 21 | 1.1 | 2.0 | 1.1 |
| Waiting room | 0.55 | 1.2 | 90:10:0 | 11 | 1.1 | 1.8 | 1.1 |
| **Public Assembly/Sports/Entertainment** | | | | | | | |
| Auditorium | 0.55 | 1.2 | 90:10:0 | 9.8 | 1.1 | 1.8 | 1.1 |
| Place of religious worship | 0.50 | 1.2 | 100:0:0 | 9.8 | 1.1 | 1.8 | 1.1 |
| Museum | 0.49 | 1.2 | 100:0:0 | 9.8 | 1.1 | 1.8 | 1.1 |
| Convention | 0.49 | 1.2 | 100:0:0 | 9.8 | 1.1 | 1.8 | 1.1 |
| Spectator area | 0.55 | 1.2 | 70:30:0 | 15 | 1.1 | 1.9 | 1.1 |
| Lobbies | 0.55 | 1.2 | 90:10:0 | 9.8 | 1.1 | 1.8 | 1.1 |
| **Residential** | | | | | | | |
| Common space | 0.55 | 1.2 | 90:10:0 | 9.8 | 1.1 | 1.8 | 1.1 |
| Dwelling unit | **0.51** | 4.3 | 15:70:15 | 40 | 1.1 | 2.2 | 1.1 |

Body text, not a footnote (p. 20, continuing §3.4 immediately below Table 2): "standard deviations are set at 10%
of the mean for both `d̄` and `C_drop`" — consistent with every GSD of 1.1 shown (a log-normal
GSD of 1.1 corresponds to
≈10% relative spread), except Dwelling unit's `d̄`/`C_drop` GSD which is also 1.1 (only its
`PBR` GSD of 4.3 is the outlier, reflecting the much wider breathing-rate range across
sleeping/sedentary/active occupants in a home per §3.4, p. 19).

**Cross-check against `occupancy_params` (`C++/src/model.cpp` / `Python/model.py`):** every
`Cdrop_GM`, `Cdrop_GSD`, `d_GM`, `d_GSD`, and `PBR_GSD` value matches for all 25 categories.
**Dwelling unit `PBR_GM`:** Jones et al. Table 3 lists **0.51 m³ h⁻¹**. The Python and C++2
implementations now use this value; this corrects the previous factor-of-10 value of 0.051 m³ h⁻¹.

---

## 8. Table 4 — Proposed ECAi and relative risk (p. 22, verbatim)

Caption (verbatim): "Proposed airflow rates of equivalent clean air per person (ECAi) and
estimates of relative risk (RR), calculated using Equation 8. The values of ECAi appear in
Table 5-1 of ASHRAE Standard 241-2023."

| Occupancy Category | Baseline outdoor airflow rate (L s⁻¹/person) | ECAi (L s⁻¹/person) | RR |
|---|---:|---:|---:|
| **Correctional Facilities** | | | |
| Cell | 4.9 | 15 | 0.44 |
| Dayroom | 3.5 | 20 | 0.35 |
| **Commercial/Retail** | | | |
| Food and beverage facilities | 5.6 | 30 | 0.24 |
| Gym | 11 | 40 | 0.32 |
| Office | 8.5 | 15 | 0.78 |
| Retail | 7.8 | 20 | 0.61 |
| Transportation waiting | 6.8 | 30 | 0.60 |
| **Educational Facilities** | | | |
| Classroom | 7.4 | 20 | 0.45 |
| Lecture hall | 2.9 | 25 | 0.22 |
| **Industrial** | | | |
| Manufacturing | 9.3 | 25 | 0.78 |
| Sorting, packing, light assembly | 11 | 10 | 1.0 |
| Warehouse | 6.5 | 10 | 0.87 |
| **Health Care** | | | |
| Exam room | 7.5 | 20 | 0.47 |
| Group treatment area | 7.5 | 35 | 0.29 |
| Patient room | 15 | 35 | 0.53 |
| Resident room | 15 | 25 | 0.69 |
| Waiting room | 5.0 | 45 | 0.16 |
| **Public Assembly/Sports/Entertainment** | | | |
| Auditorium | 3.5 | 25 | 0.36 |
| Place of religious worship | 2.8 | 25 | 0.23 |
| Museum | 4.6 | 30 | 0.31 |
| Convention | 4.6 | 30 | 0.31 |
| Spectator area | 4.7 | 25 | 0.52 |
| Lobbies | 3.5 | 25 | 0.34 |
| **Residential** | | | |
| Common space | 1.0 | 25 | 0.26 |
| Dwelling unit | 6.8 | 15 | 0.80 |

**Cross-check:** all 25 `ECAi` values match `occupancy_params[...].ECAi` in
`C++/src/model.cpp` / `Python/model.py` exactly. No discrepancy found for Table 4.

**Classification:** Table 4's `ECAi` column
is an **adjusted publication comparator**, not a raw/exact Monte Carlo fixture. §3.5 (p. 21,
verbatim, quoted in full with no elision): "an ECAi was selected so that P(I) = 0.1% is achieved 96%
of the time and was adjusted to be a multiple of 5 L s⁻¹ per person." The adjustment rule is
therefore: run the probabilistic model, take the 96th-percentile estimate of the required ECAi,
then round to a multiple of 5 L s⁻¹ per person. **The publication's own text does not state the
rounding direction** (nearest vs. ceiling) — see Discrepancy/decision D6 (§9.6). The comparison
policy requires exact equality between C++2's ceiling-to-5 value and this table for an unmodified Jones category; a one-increment
(5 L/s/p) difference is `NEAR_WITHIN_ONE_INCREMENT` (documented, not silently accepted); a
larger difference is a failed validation gate.

**Corroborating text, ≤-target wording.** Immediately preceding Table 4,
p. 22 states (verbatim): "The model ensures that the probability of infection for a one-hour
exposure is ≤ P(I) = 0.1% in 96% of spaces." This corroborates the 96th-percentile target
described in §3.5 (p. 21) and additionally makes explicit that the target is an **at-most
(≤)** probability, not an exact-equality target — i.e., the modeled infection probability for
96% of (simulated) spaces is *no greater than* 0.1%, not *equal to* 0.1%. Note the wording
differs slightly between the two pages: p. 21 says the target "is achieved 96% of the time,"
while p. 22 says the exposure probability "is ≤ P(I) = 0.1% in 96% of spaces" — both describe
the same 96th-percentile design target, but "of the time" (p. 21) and "of spaces" (p. 22) are
not verbatim identical phrasings and are not merged into a single paraphrase here.

---

## 9. Discrepancy log

Each entry records the description, sources compared, and status. D2 and D7 are resolved.
D1 and D5 are terminology-only publication inconsistencies. D3, D4, D6, D8, and D9 are standing
C++2-specific decisions, not publication resolutions.

### D1 — `R` vs. `RD` in Equation 4 (UNRESOLVED, low scientific impact — terminology only)

Eq. 4 (p. 10) is typeset as `VER_i = PBR · V*_drop · GVL_i · VF · RD`. The surrounding prose
(p. 10, verbatim, quoted here in full with no elision): "`R` is the fraction of
viability lost due to rapid decay soon after the initial exhalation of respiratory particles
caused by chemical and physical changes to them during evaporation [21, 22, 23, 24]. They then
revert to a more stable biological half-life. There is not enough knowledge to consider this
phenomenon and so we use a value of `R` = 1 hereon." The Nomenclature (p. 4) and Table 1 (p. 15)
define and use only `R` ("Fractional reduction of viability after emission," constant = 1.0,
ref. [22]); `RD` is never separately defined anywhere in the publication. This is an internal
publication inconsistency (almost certainly a typesetting artifact, since `D` is already used
elsewhere for exposure duration in hours and its reappearance here would be dimensionally and
contextually wrong).
`C++/src/model.cpp` already encodes this constant under the local variable name `RD` (`double
RD = 1;`), which is a naming choice borrowed from the equation's typography, not a second
quantity — functionally it is `R`. **C++2 decision (C++2-specific):** use the symbol/name `R`
(not `RD`) for this constant in code and documentation, value 1.0. This does not change any
numeric output.

### D2 — Dwelling unit `PBR_GM`: 0.51 (publication) — RESOLVED

Table 3 (p. 21) lists `0.51` m³ h⁻¹ for the Dwelling unit breathing-rate geometric mean. The Python
and C++2 implementations now use this value. The previous `0.051` m³ h⁻¹ value was a factor-of-10
implementation error and is no longer used.

**Resolution:** Python and C++2 use the publication value `PBR_GM = 0.51` m³ h⁻¹ for the Dwelling
unit category. Only the Dwelling unit category (I0 = 6) was affected.

### D3 — Zero-infector short-circuit changes total draw count relative to the existing implementation (C++2-specific design decision, not a publication discrepancy)

The existing C++/Python implementation calls `sample_parameters()` (which draws `λ`, `γ`, and
susceptible `PBR`) **unconditionally**, before the binomial infector count is known; only the
infector-side `QER` draws are skipped when `n_infected == 0`. `IMPLEMENTATION_PLAN.md` §4.3
requires that C++2, for a zero-infector row, "consume no other sampled values" — i.e., no
susceptible-side draw either. Both approaches yield the same final `ECAi = 0` for a
zero-infector row (the algebra is unreachable either way), but they consume a different number
of LHS values overall, which matters for C++2's exact-size aggregate-demand bank sizing
(`N_active` vs. `N`). **This is a C++2-specific architectural decision, not attributable to
Jones et al.**, who do not discuss LHS consumption bookkeeping at all.

### D4 — Independent per-person `PBR` draws: a deliberate C++2 deviation from the publication's stated assumption

**The publication's explicit statement (p. 10, verbatim, quoted in full — no elision).**
Immediately before Eq. 4, the publication states: "We assume that the breathing rate of each
infected person is the same as that for all susceptible people when `VER_i` is determined to be
[20]." The bracketed `[20]` is the publication's own numbered-reference citation marker (to
reference [20] in its bibliography), **not** a reference to Equation 4; an earlier draft of this
document mis-transcribed that bracket as `[Eq. 4]`, which this revision corrects.

**A second, more explicit statement (p. 19, verbatim, §3.4 Occupancy category inputs).** "The
breathing rate, `PBR`, appears twice in the Wells-Riley equation, first in Equation 2
representing breathing of susceptible people, and again in Equation 4 for infected people.
Herein they are assumed to be identical for all occupants." This is the clearer, unambiguous
form of the same assumption; both statements say the same thing (a single shared
`PBR` value for every occupant in a space, susceptible and infected alike), and neither is
ambiguous on its own once read in full.

**Figure 1 evidence pointing the other way (p. 13, verbatim caption: "A model of the
relationships between infection risk probability and pathogen removal rate. Blue: deterministic
input; Red: biological; Green: unique inputs for each infected person.").** Figure 1's system
diagram encloses "Infected person breathing rate" inside the green dashed box explicitly
labelled "Unique for each infected person," alongside "Ratio of volumes of respiratory fluid to
exhaled air," "Viral load of genomic material," "Viable fraction," "Mean droplet diameter,"
"Evaporation coefficient," and "Respiratory particles per unit volume of exhaled air." The
figure therefore visually depicts the infector's `PBR` as a per-infector unique input, which is
in tension with the identical-for-all-occupants prose quoted above.

**C++2 decision — a labelled deviation, not a claim of publication equivalence.** The existing
C++/Python reference code already draws a fresh `PBR` for the susceptible (in
`sample_parameters`) and a fresh, independent `PBR` for every single infector (inside the
`QER()` call made once per infector in the summation loop), contradicting the literal p.10/p.19
"identical for all occupants" text. Independent per-person `PBR` sampling is **mandatory in
C++2 by explicit user design** (per `IMPLEMENTATION_PLAN.md` §2 "Explicit design decision on
person-specific values"): every infector receives its own complete parameter set including its
own `PBR_i`, and the susceptible receives an independent `PBR_s`, with no realized value ever
shared between roles. This document does **not** present that choice as a faithful reproduction
or reasonable reading of the publication's stated assumption — it is recorded here as a
**deliberate, explicit C++2 deviation** from the publication's own words (p. 10 and p. 19 both),
motivated by (a) the existing code's prior behavior and (b) Figure 1's contrary depiction, but
making no claim that C++2's output is therefore equivalent to what Jones et al. describe in
prose. The *exact-size, one-LHS-list-per-named-parameter* architecture that enforces
independence structurally (§11.1 below) is C++2-specific engineering, layered on top of this
scientific deviation.

### D5 — RTD "deposition" (Nomenclature, p. 4) vs. "absorption" (Table 1, p. 15) (UNRESOLVED, terminology only, no numeric impact)

Same symbol `RTD`, same Uniform(0.43, 0.65) distribution, same citation [30]; the publication
uses "deposition" in one place and "absorption" in another. No action required beyond noting
it; C++2 documentation uses "respiratory tract deposition fraction" (the Nomenclature's term
and the term used throughout `IMPLEMENTATION_PLAN.md`).

### D6 — Rounding direction for the multiple-of-5 adjustment is not stated (C++2-specific freeze)

§3.5 (p. 21) states the ECAi "was adjusted to be a multiple of 5 L s⁻¹ per person" without
specifying ceiling, floor, or nearest. The governing plan's rounding decision uses
`ceiling_5 = 5 * ceil(raw / 5)`, and the existing
`C++/src/ecai.cpp` already implements ceiling (`std::ceil(percentile_96 / 5.0) * 5`). **C++2
uses ceiling-only rounding.** This is a C++2-adopted convention, not a literal statement from
Jones et al.

### D7 — Category naming: publication full names vs. existing code's abbreviations — RESOLVED

Table 2/3/4 use full descriptive names ("Food and beverage facilities," "Transportation
waiting," "Group treatment area," "Patient room," "Resident room," "Place of religious
worship," "Spectator area," "Common space," "Dwelling unit," "Sorting, packing, light
assembly"). The existing C++/Python `occupancy_params` map uses shortened keys ("Food",
"Transportation", "Group", "Patient", "Resident", "Place", "Spectator", "Common", "Dwelling",
"Sorting"). All underlying numeric values are identical under either naming; this is a naming
convention only.

**Resolution:**

- **Canonical identifiers** are the publication's full names exactly as printed in Tables 2/3/4
  (e.g., `"Food and beverage facilities"`, `"Group treatment area"`, `"Place of religious
  worship"`, `"Dwelling unit"`), never the existing code's shortened keys.
- **Normalization for input matching** is case-insensitive, and treats spaces, commas, and
  hyphens as equivalent, interchangeable delimiters (so `"food-and-beverage-facilities"`,
  `"Food, and beverage facilities"`, and `"FOOD AND BEVERAGE FACILITIES"` all normalize to the
  same canonical key). No other punctuation is normalized.
- **Legacy short aliases:** C++2 v1 unconditionally accepts exactly these ten additional input
  aliases at the CLI/API boundary: `"Food"`, `"Transportation"`, `"Group"`, `"Patient"`,
  `"Resident"`, `"Place"`, `"Spectator"`, `"Common"`, `"Dwelling"`, and `"Sorting"`. Each alias
  is an exact, explicitly enumerated string under the same case/delimiter normalization; no alias
  is inferred.
- **Canonical metadata always emits full publication names.** Regardless of which input form
  (full name or legacy alias) was used to select a category, any output that reports the
  category identifier — CSV headers, JSON metadata, log messages, category listings — always
  emits the full publication name, never a legacy alias.
- **No fuzzy fallback.** There is no edit-distance, prefix, or typo-tolerant matching of any
  kind. An input string that does not exactly match a canonical full name (under the
  case/delimiter normalization above) or an explicitly enumerated legacy alias is rejected with
  a diagnostic, never guessed.

`category_data.h` and `category_data.cpp` implement this rule. It is a C++2-specific engineering
decision; Jones et al. does not define a CLI vocabulary.

### D8 — Percentile interpolation convention is not stated by Jones et al. (C++2-specific freeze, per verified current implementation)

Nowhere does the publication specify a finite-sample percentile interpolation formula for the
"96% of the time" statement in §3.5. The existing
`C++/src/ecai.cpp` and `C++/analysis/ecai_uncertainty_study.cpp` both use linear interpolation
matching `numpy.percentile`'s default ("linear") method: `rank = p·(n−1)`, interpolate between
`floor(rank)` and `ceil(rank)`. **C++2 freezes this same linear-interpolation convention**,
inherited from the verified current implementation, not from Jones et al.

### D9 — Reproduction-mode community rate (0.03125) is an out-of-scope correspondence value, not a Jones default

`C++/analysis/ecai_uncertainty_study.cpp`'s `reproduction` mode uses a default healthcare
community rate of `0.03125`, explicitly documented in-file as a deliberate deviation used only
for matching an external correspondence/notebook, distinct from the Jones-published 3.0%
(Table 1, p. 15). **C++2 v1 uses the Jones-published defaults (1% general / 3% healthcare,
Table 1)**, not this analysis-tool value.

---

## 10. Dimensional analysis

Quanta (`Q`, `QER`) are treated as a dimensionless dose-count unit throughout — the publication
states (p. 9) "Quanta cannot be independently measured and [Eq. 1] is the defining
relationship," i.e., quanta have no independent physical dimension; they cancel algebraically
like a pure count. With that convention:

| Quantity | Expression | Unit build-up | Result |
|---|---|---|---|
| `V*_drop` (Eq. 5) | `(π/6)·(d̄·E)³·C_drop` | `m³ (particle volume) × particles·m⁻³ (air)` | dimensionless (m³ fluid per m³ air) |
| `GVL_i` as used in Eq. 4 | `10^Normal(7.0,1.4)` [log₁₀ RNA copies mL⁻¹] `× 1e6` | `RNA copies·mL⁻¹ × mL·m⁻³` (1 m³ = 1e6 mL) | RNA copies·m⁻³ |
| `VER_i` (Eq. 4) | `PBR · V*_drop · GVL_i · VF · R` | `(m³·h⁻¹)·(dimensionless)·(RNA copies·m⁻³)·(virions per RNA copy)·(dimensionless)` | virions·h⁻¹ per person |
| `QER_i` (Eq. 3) | `RTD · VER_i / DK` | `(dimensionless)·(virions·h⁻¹) / (virions per quantum, implicit in `DK`'s definition as "reciprocal of the probability a single virion initiates an infection")` | quanta·h⁻¹ per person |
| `Q` (Eq. 2) | `[PBR·D·(1−μ)²/(φ·VOL)]·ΣQER_i` | `(m³·h⁻¹·h·dimensionless)/(h⁻¹·m³) · (quanta·h⁻¹)` = `(m³/(h⁻¹·m³)) · (quanta·h⁻¹)` = `h · quanta·h⁻¹` | dimensionless (quanta, treated as count) |
| `φ·VOL` (Eq. 6, term2 in inversion) | `(γ+λ)·VOL + T_ECAi` | `h⁻¹·m³ + m³·h⁻¹` | m³·h⁻¹ |
| `T_ECAi` (inversion of Eq. 2/6, i.e. `term1 − term2`) | `PBR·D·mask_factor·QER_sum/Q_target − VOL·(γ+λ)` | `(m³·h⁻¹·h)·(quanta·h⁻¹)/(quanta, treated as count) − m³·h⁻¹·h⁻¹·... ` → both terms reduce to m³·h⁻¹ | m³·h⁻¹ |
| `ECAi` (Eq. 7) | `T_ECAi / (I0 × 3.6)` | `(m³·h⁻¹) / (dimensionless × 3.6)`; `1 m³·h⁻¹ = 1000 L / 3600 s = (1/3.6) L·s⁻¹` | L·s⁻¹ per person |

`DK` is defined (p. 10) as "the reciprocal of the probability that a single virion initiates an
infection." Under a single-hit/Poisson dose-response model, `quanta_dose = virions / DK`, so
`DK` carries implicit units of virions-per-quantum; this is how Eq. 3 converts a virion rate
into a quanta rate. This mechanism is stated in the publication's prose (p. 10) but the
unit-cancellation is not spelled out algebraically there — the derivation above is the C++2
documentation of it, consistent with, not contradicting, the cited text.

The constant `3.6` in Eq. 7 converts m³ h⁻¹ to L s⁻¹ (`1 m³/h = 1/3.6 L/s`), confirmed
algebraically above. The constant `1e6` (Table 1 footnote) converts GVL from RNA copies mL⁻¹ to
RNA copies m⁻³ (`1 m³ = 1e6 mL`), confirmed algebraically above.

---

## 11. C++2-specific decisions frozen by this document

Everything in this section is an engineering decision required by `IMPLEMENTATION_PLAN.md` but
**not stated by Jones et al.** Each is explicitly labelled C++2-specific.

### 11.1 Person-specific sampling (deliberate C++2 deviation from the shared-PBR prose; see D4)

C++2 implements **aggregate-demand-conditioned marginal one-dimensional LHS**. Each complete
named bank is marginally stratified. Simulation tuples, within-infector parameter vectors,
susceptible–infector pairings, multi-infector sums, and event-selected prefixes are not
joint-stratified. The Jones LHS matrix is underspecified, so this design is not demonstrated
publication-equivalent. See [`SAMPLING_DESIGN.md`](SAMPLING_DESIGN.md) for the implemented bank,
cursor, assignment, and resource contracts.

- Susceptible parameters (`PBR_s`, `γ`, `λ`) are drawn from one exact-size LHS list each, sized
  `N_active`.
- Every infector-specific parameter (`PBR_i`, `C_drop`, `d`, `E`, `GVL`, `VF`, `RTD`, `DK`) is
  drawn from its own exact-size LHS list, sized `K_total` (aggregate demand across all
  infectors in the replicate).
- No realized value is ever shared between the susceptible and any infector, or between two
  infectors, even though they may share a distribution family and parameters.
- A zero-infector simulation consumes zero values from every list (susceptible and infector
  alike) — see D3.

### 11.2 RNG, open-interval conversion, unbiased bounded-integer generation, Fisher-Yates shuffle (C++2-specific, reproducibility scope qualified below — option (a) of `IMPLEMENTATION_PLAN.md` §6)

Jones et al. says nothing about RNG mechanics (§3.1, p. 12 mentions only that "a matrix of
uniformly distributed random numbers was created using Latin Hypercube sampling" in MATLAB/
Excel — no algorithm, no seed, no reproducibility claim). C++2 freezes:

- **Base engine:** `std::mt19937_64` (64-bit Mersenne Twister). Its generation algorithm is
  specified by the C++ standard (given a
  seed, every conforming implementation produces the identical 64-bit output sequence). This
  is why `std::uniform_real_distribution` and `std::shuffle` are explicitly avoided for
  anything requiring cross-toolchain bit-for-bit reproducibility: the C++ standard does not mandate their exact algorithms, only their
  statistical properties, so two conforming standard libraries may legally disagree given an
  identical engine state.
- **Open-interval integer→double conversion:** draw one raw `uint64_t` `x` from the engine and
  take the high 53 bits, `m = x >> 11` (`0 ≤ m < 2^53`). Compute the binary64 candidate
  `u = (double(m) + 0.5) / 2^53`. For the maximum `m`, binary64 round-to-nearest can round the
  numerator to `2^53` and therefore produce `u == 1.0`; when the candidate is not less than
  `1.0`, C++2 replaces it with `nextafter(1.0, 0.0)`. The source also guards the lower endpoint:
  when the candidate is not greater than `0.0`, it replaces it with `nextafter(0.0, 1.0)`.
  Stratum construction then computes `lower + u * (upper - lower)`; if binary64 arithmetic rounds
  that point onto or beyond either stratum boundary, C++2 nudges it inside with
  `nextafter(lower, upper)` or `nextafter(upper, lower)`, respectively. The returned unit value and
  every stratum point are strictly open after these corrections.
- **Unbiased bounded-integer generation** (needed for Fisher-Yates), for a nonzero `uint64_t`
  bound `n` representing `[0,n)`: calculate `threshold = uint64_t(0 - n) % n` using defined
  unsigned wraparound; repeatedly draw raw `uint64_t x` values until `x >= threshold`; return
  `x % n`. This avoids trying to represent `2^64` in a `uint64_t`, is fully executable using
  standard unsigned arithmetic, and has no modulo bias.
- **Fisher-Yates (Durstenfeld) shuffle:** for `i` from `n−1` down to `1`, draw
  `j = unbiased_bounded_integer(0, i+1)` (i.e., uniformly in `[0, i]`), swap elements `i` and
  `j`. Fully specified given the bounded-integer generator above.
- **Reproducibility scope:** the integer seed fold, raw `mt19937_64` engine outputs,
  threshold-rejection bounded integers, and Fisher-Yates index permutation use standard-defined
  integer behavior and explicit algorithms. Bit identity for the open-double conversion,
  `nextafter` endpoint/boundary corrections, binary64 stratum arithmetic, or shuffled double
  values additionally requires supported IEEE-754 binary64 semantics and the same relevant
  compiler, standard/mathematical library, target architecture, and floating-point flags. It is
  not claimed for every conforming C++17 implementation. The claim does **not** automatically
  extend to inverse-CDF transforms applied to a uniform draw. Normal, log-normal, log10-normal,
  beta, binomial-PPF, and downstream arithmetic feeding into `QER`, `Q`, and `ECAi` call `erf`/`erfinv`/`log`/`exp`/`pow`/incomplete-beta routines from the
  platform's libm, whose exact bit-level rounding is not standardized across compilers, libm
  implementations, or CPU architectures. Consequently: **end-to-end `ECAi` output is
  bit-for-bit reproducible only on a frozen, explicitly named toolchain** — a specific compiler
  (name and version), a specific `libstdc++`/libm build, a specific target architecture, and a
  specific floating-point flag set — held fixed and documented for a given set of results.
  Across a different supported toolchain (different compiler, different libm, different architecture,
  or different flags), transformed and end-to-end C++2 values are assessed for **numerical
  agreement within the frozen tolerances of
  §11.5**, not bit-identical output. **Fast-math (`-ffast-math`, `/fp:fast`, or equivalent) is
  prohibited anywhere reproducibility is required, and floating-point contraction (FMA fusing of
  multiply-add) must be explicitly disabled** (e.g. `-ffp-contract=off`) rather than left at a
  compiler's default, because both change rounding/associativity in ways that break even the
  frozen-toolchain bit-reproducibility guarantee above. No `std::uniform_real_distribution`,
  `std::discrete_distribution`, or `std::shuffle` may be used anywhere a reproducible sequence
  is required.

### 11.3 Seed architecture: SplitMix64 mixing + canonical configuration fingerprint (C++2-specific, per `IMPLEMENTATION_PLAN.md` §6)

- **Mixing function:** the public-domain SplitMix64 algorithm (Vigna), applied to a running
  64-bit state:
  ```
  state += 0x9E3779B97F4A7C15
  z = state
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
  z = (z ^ (z >> 27)) * 0x94D049BB133111EB
  z = z ^ (z >> 31)
  return z
  ```
- **Stream seed derivation:** define `splitmix64(x)` as the exact value transformation above,
  including the initial addition. Initialize `state = splitmix64(base_seed ^
  0x4350503253454544)` (ASCII domain tag `CPP2SEED`). For each `uint64_t element` of the ordered
  tuple `(design_version, effective_configuration_fingerprint, replicate_index, stream_id,
  list_length)`, set `state = splitmix64(state ^ splitmix64(element))`. The final `state` is the
  seed used to construct that stream's `mt19937_64`. This exact fold is normative; no alternate
  “XOR or additive” choice is permitted. Changing a tuple component is expected to avalanche
  the seed, while uniqueness remains an exhaustively tested supported-domain property rather
  than a mathematical collision-impossibility claim.
- **Canonical `effective_configuration_fingerprint` encoding:** per `IMPLEMENTATION_PLAN.md`
  §6, `std::hash<std::string>` is explicitly prohibited (its output is not guaranteed stable
  across standard library versions/builds). C++2 instead: (1) serializes the effective
  configuration as an ordered, fixed-field ASCII string — category name, `I0`, `VOL`, mask
  efficiency, community rate, target probability, duration, percentile, confidence level,
  simulations `N`, replicates `R` — each numeric field formatted with a fixed, round-trip-safe
  format (`%.17g` for doubles, decimal for integers), fields joined with a fixed delimiter
  (`|`); then (2) hashes that byte string with the 64-bit FNV-1a algorithm (offset basis
  `0xcbf29ce484222325`, prime `0x100000001b3`), a simple, fully specified, portable,
  non-cryptographic hash with a stable, standard-independent definition. This fingerprint is
  itself one input to the SplitMix64-based stream-seed derivation above, not a substitute for
  it.
- **Stream IDs:** a fixed enum, one value per named stochastic parameter (binomial infector
  count, susceptible `PBR`, `γ`, `λ`, infector `PBR`, `C_drop`, `d`, `E`, `GVL`, `VF`, `RTD`,
  `DK`), assigned once and never renumbered/reused for a different parameter.

### 11.4 Resource / memory budget policy (C++2-specific)

Before allocating any parameter bank, C++2 must compute the total number of `double` elements
required (`N` for the binomial list, `N_active` for each of the 3 susceptible-list families,
`K_total` for each of the 8 infector-list families: total = `N + 3·N_active + 8·K_total`) using
checked/overflow-safe arithmetic (reject if any intermediate sum or product would overflow
`size_t`, or if `N·I0` — the maximum possible `K_total` — would overflow first). C++2 freezes a
**default soft memory budget of 4 GiB** for the combined size of all parameter banks in a
single replicate (`total_elements × 8 bytes ≤ 4 GiB`, i.e., `total_elements ≤ 2^29`); a run
that would exceed this budget must be rejected with a diagnostic **before** any allocation
occurs, leaving no partial output. This budget is an engineering default, not derived from
Jones et al. The CLI exposes the `N` ceiling but does not expose a resource-budget override.

### 11.5 Validation acceptance table

All tolerances below are current C++2-specific acceptance criteria. They may be changed only with
documented numerical evidence and independent approval.

| Check | Tolerance | Basis |
|---|---|---|
| Deterministic algebra (QER/ECAi hand-calculated fixtures) | absolute 1e-12 **or** relative 1e-10 | Current bounded regression criterion |
| Normal / lognormal / log10-normal inverse CDF vs. SciPy | absolute 1e-9 **or** relative 1e-9 | Existing Acklam+Halley implementation's advertised ~1.15e-9 max relative error before refinement; refinement tightens further, so this is an achievable, conservative bound |
| Uniform inverse CDF vs. hand-calculated | absolute 1e-14 | Pure affine formula; double-precision rounding floor only |
| Beta inverse CDF (bisection on continued-fraction CDF) vs. SciPy | absolute 1e-9 **or** relative 1e-8 | 60-bisection-iteration interval width is far tighter (~2⁻⁶⁰); bound set by `betacf`'s own `EPS=3e-16` convergence and double precision, not by iteration count |
| Binomial PPF (smallest `k` with `CDF(k) ≥ u`) vs. SciPy | exact integer match required; a ±1 discrepancy is tolerated **only** when `u` is within relative 1e-12 of the exact CDF value at the boundary `k`, and must be documented, not silently accepted, even then | Discrete step function; ties are a floating-point boundary artifact, not a genuine algorithm disagreement |
| Percentile interpolation / mean / sample SD / SE known-answer vectors | absolute 1e-9 | Pure arithmetic, no iterative approximation |
| Student-t 0.975 critical value (incomplete-beta + bisection) vs. SciPy | absolute 1e-6 | Existing bisection-based implementation style; looser than the beta-CDF bound to allow for the additional bisection-on-CDF composition, tightened later only with documented evidence |
| Publication ECAi comparison (ceiling-to-5 vs. Table 4, §8 above) | exact match = pass; exactly one 5 L/s/person increment = `NEAR_WITHIN_ONE_INCREMENT` (documented); more than one increment = failed gate | Current comparison policy |

---

## 12. Implemented conventions summary

- **Defaults (Table 1 + §3.5, all verified above):** category community rate 1% general / 3%
  healthcare; target hourly infection probability 0.1% (0.001); exposure duration 1.0 h;
  percentile 96%; mask efficiency 0.0 general / 0.30 for the five categories under the Table 2
  **Health Care** group (Exam room, Group treatment area, Patient room, Resident room, Waiting
  room). The value `μ = 0.30` is stated in §3.3, p. 16; applying it to those five categories is a
  C++2 mapping inferred from the Table 2 group heading, not a separately enumerated list in the
  prose. `R` (viability loss fraction) is fixed at 1.0.
  25 categories total, in the publication order given in Table 2 (§6 above).
- **Percentile convention:** linear interpolation (`numpy.percentile` "linear" default),
  frozen per D8 from the verified current implementation, not stated by Jones et al.
- **Rounding:** report both the raw (unrounded) 96th-percentile mean, SD, SE, and CI in
  L s⁻¹ per person, and separately `5 × ceil(raw_mean / 5)` as the ceiling-to-5 comparator;
  never feed the ceiling-adjusted value into SD/SE/CI. Ceiling direction frozen per D6.
  Confidence intervals characterize the mean of the replicate percentile estimator only (an
  approximate randomized-LHS Monte Carlo interval), never the population percentile or
  individual ECAi variation.
- **RNG/seeds/resources/tolerances:** frozen per §11.2–11.5 above; all explicitly C++2-specific
  engineering decisions, not attributable to Jones et al.
- **Probability mode (`--mode probability`, §3.1):** a C++2-specific forward evaluation of Eq. 1
  at a fixed input ECAi (published Table 4 value by default, or an explicit `--input-ecai`),
  reporting the configurable within-replicate percentile (default `0.96`) of per-simulation
  infection probability. Reuses the same seeds, banks and person-specific sampling as ECAi mode;
  see `SEED_POLICY.md` for the fingerprint implications and `CLI_REFERENCE.md` for the flags and
  output schema.
- **Adjudicated items:** D2 uses the publication value `PBR_GM = 0.51 m³/h` for Dwelling unit;
  D7 uses full publication category names as canonical identifiers with strict normalized input
  matching and unconditionally accepts exactly the thirteen explicit legacy aliases listed in D7.
  D1 and D5 remain terminology-only publication inconsistencies with no numeric effect.

---

## 13. Validation boundary

The current verification status, bounded reference cases, exact counters, and non-equivalence
limitations are recorded in [`VALIDATION.md`](VALIDATION.md). Passing implementation tests does
not establish independent scientific validity or publication reproduction.

## Risk-averse policy

Risk-averse mode conditions the initial binomial infector count on at least one infector for both ECAi and forward infection probability. Its conditional PPF is inverted in CDF space, with explicit endpoint and limiting branches; it is not an affine remapping of a uniform. `I0==0` and `p<=0` are rejected. For `p==1`, the count is deterministic and `zero_mass=0`; any `zero_mass==1` limiting branch is recorded in provenance.

Risk-averse results use distinct schemas (`stage13-ra-v1` and `stage13-prob-ra-v1`), require zero-infector count 0, active count N, and total infectors at least N. Population comparators are suppressed as `NOT_APPLICABLE_RISK_AVERSE`; this estimand is conditional on an infector and is not comparable with the published community comparator.
