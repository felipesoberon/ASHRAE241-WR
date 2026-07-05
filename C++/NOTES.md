C++ Port -- Design Notes and Validation Record
==============================================

Condensed record of the design decisions and validation behind the C++ port of
the ASHRAE 241 risk model. The original task-by-task implementation plans
(plan.md, plan_gui.md) were removed once the work landed; they remain in git
history if the full detail is ever needed.

Equivalence to the Python model
-------------------------------
The C++ port is statistically equivalent to the frozen Python reference, not
bit-identical: it uses its own RNG to fill the LHS buffers, so individual draws
differ, but the distributions, LHS stratification, formulas, and all 25 x 11
occupancy parameters match exactly.

Validation at N=1,000,000 vs the stored Python 1M reference:
  - Infection probability: 25/25 categories within 1.4% relative (mean ~0.6%).
  - Required ECAi: every meaningful value (> 5 L/s/p) within 3.6%. The only
    larger relative gaps (Patient, Warehouse) are near-zero categories where the
    absolute difference is < 0.25 L/s/p -- Monte Carlo noise on a ~0 quantity,
    not a model discrepancy.
Inverse CDFs are unit-tested against scipy values to 1e-8; zero-infector
fractions match the theoretical (1 - rate)^I0.

Key design decisions
--------------------
- Category order: iterated via the category_order vector (Python/ASHRAE
  insertion order), not the alphabetical std::map, so output rows match Python.
- Percentiles: linear interpolation, rank = q * (N - 1), matching numpy default.
- Raw data: custom .bin format (not .npz), read by percentiles / boxplot.
- Community infection rate (general vs healthcare):
    * --community-rate-general (default 1%) and --community-rate-healthcare
      (default 3%); healthcare = Exam, Group, Patient, Resident, Waiting.
    * is_healthcare_category() classifies by the elevated 3% default
      (threshold 0.02).
    * Each flag defaults to the sentinel -1 (unset), so the category keeps its
      own default -- runs without flags are byte-for-byte the original model.
- GUI (mainGUI_C++.py): Windows Python/Tkinter driving the C++ ELF engine via
  WSL. The WSL distro is auto-detected and paths are converted with wslpath, so
  it is not tied to a specific distro name or drive letter.
- Driver analysis (--save-inputs):
    * probability_ecai --save-inputs writes per-simulation input parameters
      (21 fields: PBR, lambda_bio, gamma, n_infected, phi, QER_sum,
      mask_factor, Q, P, and 12 QER component values) to a binary file.
    * Implemented via QER_with_inputs() and infection_probability_with_inputs()
      -- new functions that mirror the originals exactly (same RNG call
      sequence, same formulas) but also return intermediate values.
    * The original QER() and infection_probability() are untouched; when
      --save-inputs is absent, the original code path runs unchanged.
    * Verified: QER_sum == qer_QER_val for n_infected=1 (diff = 0),
      P == 1 - exp(-Q) (diff = 1.11e-16), P96 from inputs file matches
      terminal output.

Python analysis scripts (analysis/ folder)
------------------------------------------

The following scripts were added to support the presentation analysis plan.
They read the C++ .bin format via bin_reader.py:

  bin_reader.py        Reads --save-all .bin format (numpy arrays)
  exceedance.py        Exceedance fraction table (mean, median, P96, % sims
                       > 0.1/0.5/1.0%)
  ccdf_plot.py         Complementary CDF (exceedance curve) plot
  percentile_heatmap.py 25x11 percentile heatmap with 0.1% contour
  boxplot_bin.py       Matplotlib box-and-whisker on log scale, reads .bin
  summary_table.py     Summary table with compliance grouping (A/B/C)
  inputs_reader.py     Reads --save-inputs .bin format (count x 21 arrays)
