Analysis Scripts
================

Post-processing tools for the raw simulation data produced by the
ASHRAE 241 risk model.

There are two data formats:

- **Python `.npz`** -- written by `probabilityECAi.py --save-all`.
  Contains one `float32` array of raw infection probabilities (0-1)
  per occupancy category. Read by `percentiles.py` and `boxplot.py`.

- **C++ `.bin`** -- written by `C++/build/probability_ecai --save-all`.
  Custom binary format with the same content. Read by `bin_reader.py`
  and all C++-format scripts below.

- **C++ inputs `.bin`** -- written by
  `C++/build/probability_ecai --save-inputs`. Contains 21 fields per
  simulation (sampled parameters + intermediate values for driver
  analysis). Read by `inputs_reader.py`.

Python .npz scripts
-------------------

percentiles.py
--------------

Computes percentile statistics from a raw `.npz` file. Produces two
reports:

- **table** -- per-category percentile table (default: 50th, 75th,
  96th).
- **threshold** -- finds the highest percentile at which *all*
  categories stay below a target probability (default 0.1%) using a
  coarse-to-fine search: scans 25 to 95 in steps of 20, then
  repeatedly refines within the bracketing interval at steps of 10,
  5, and 1. Each scan flags the worst-case category and the maximum
  probability across categories.

### Usage

    # Both reports on the default file
    python analysis/percentiles.py

    # Point at a specific raw file
    python analysis/percentiles.py run10k.npz

    # Only the 0.1% threshold sweep
    python analysis/percentiles.py --report threshold

    # Custom percentiles in the table
    python analysis/percentiles.py --report table --percentiles 50 90 95 99

    # Custom target for the threshold sweep (in %)
    python analysis/percentiles.py --target 0.05

boxplot.py
----------

Draws a horizontal box-and-whisker plot of the raw probabilities,
one box per occupancy category (log x-axis). Boxes show quartiles
with a black median line and 1.5x-IQR whiskers. Each category's 96th
percentile is marked with a red tick, along with a dashed reference
line at the target probability.

### Usage

    python analysis/boxplot.py run10k.npz --save boxplot.png

C++ .bin scripts
----------------

bin_reader.py
-------------

Utility that reads the C++ binary raw-data format. Returns a dict
of {category: numpy float32 array}. Used by all C++-format scripts.

    python analysis/bin_reader.py probECAi_1M_cpp.bin

boxplot_bin.py
--------------

Same as boxplot.py but reads the C++ `.bin` format directly via
bin_reader.py. Horizontal box-and-whisker on log scale, 96th
percentile marked, 0.1% target line.

### Usage

    python analysis/boxplot_bin.py probECAi_1M_cpp.bin --save boxplot.png

### Options

| Option      | Default                    | Description                                       |
|-------------|----------------------------|---------------------------------------------------|
| `infile`    | `probabilityECAi_raw.bin`  | Raw `.bin` file to plot (positional).             |
| `--target`  | `0.1`                      | Target probability (%) drawn as a reference line. |
| `--xmin`    | `0.001`                    | Left edge of the log x-axis in %.                 |
| `--xmax`    | `100`                      | Right edge of the log x-axis in %.                |
| `--save`    | *(show interactively)*     | Save the figure to this path instead of showing.  |

exceedance.py
-------------

Computes exceedance fractions from raw C++ `.bin` data. For each
category, reports mean, median, 96th percentile, and the fraction
of simulations exceeding 0.1%, 0.5%, and 1.0% targets. Optionally
writes results to CSV.

### Usage

    python analysis/exceedance.py probECAi_1M_cpp.bin --csv results.csv

### Options

| Option      | Default                    | Description                                       |
|-------------|----------------------------|---------------------------------------------------|
| `infile`    | `probabilityECAi_raw.bin`  | Raw `.bin` file to analyze (positional).           |
| `--targets` | `0.1 0.5 1.0`              | Target probabilities in % for exceedance.          |
| `--csv`     | *(none)*                   | Write results to this CSV file.                   |

ccdf_plot.py
------------

Plots complementary CDF (exceedance curves) from raw C++ `.bin`
data. One curve per category on a log x-axis, with the 0.1% target
line. Categories are sorted by 96th percentile (worst at top of
legend).

### Usage

    python analysis/ccdf_plot.py probECAi_1M_cpp.bin --save ccdf.png

### Options

| Option      | Default                    | Description                                       |
|-------------|----------------------------|---------------------------------------------------|
| `infile`    | `probabilityECAi_raw.bin`  | Raw `.bin` file to plot (positional).             |
| `--target`  | `0.1`                      | Target probability (%) drawn as a vertical line. |
| `--xmin`    | `0.001`                    | Left edge of the log x-axis in %.                 |
| `--xmax`    | `100`                      | Right edge of the log x-axis in %.                |
| `--save`    | *(show interactively)*     | Save the figure to this path instead of showing.  |

percentile_heatmap.py
---------------------

Draws a 25x11 percentile heatmap (categories x percentiles) from
raw C++ `.bin` data, color-coded by infection probability (log
scale). Categories sorted by 96th percentile. A contour line marks
the 0.1% target. Cell values are annotated.

### Usage

    python analysis/percentile_heatmap.py probECAi_1M_cpp.bin --save heatmap.png

### Options

| Option          | Default                     | Description                                    |
|-----------------|-----------------------------|------------------------------------------------|
| `infile`        | `probabilityECAi_raw.bin`    | Raw `.bin` file to analyze (positional).       |
| `--target`      | `0.1`                       | Target probability (%) drawn as a contour.     |
| `--percentiles` | `50 60 65 70 75 80 85 ...`  | Percentiles for the columns.                   |
| `--save`        | *(show interactively)*      | Save the figure to this path instead of showing.|

summary_table.py
----------------

Produces a comprehensive summary table from raw C++ `.bin` data.
One row per category with: ECAi, median, P75, P90, P96, P99,
exceedance fraction >0.1%, and compliance group.

Compliance groups:
- **A** (compliant): P96 < 0.1%
- **B** (tail risk): median < 0.1% but P96 >= 0.1%
- **C** (systematic): median >= 0.1%

### Usage

    python analysis/summary_table.py probECAi_1M_cpp.bin --csv summary.csv

### Options

| Option     | Default                    | Description                                    |
|------------|----------------------------|------------------------------------------------|
| `infile`   | `probabilityECAi_raw.bin`  | Raw `.bin` file to analyze (positional).       |
| `--target` | `0.1`                      | Target probability (%) for exceedance/grouping.|
| `--csv`    | *(none)*                   | Write results to this CSV file.                |

C++ inputs script
------------------

inputs_reader.py
-----------------

Reads the C++ `--save-inputs` binary format (per-simulation input
parameters for driver/sensitivity analysis). Returns a dict of
{category: 2D numpy array of shape (count, 21)}.

Field order (21 doubles per simulation):

  0: PBR_sample    1: lambda_bio    2: gamma
  3: n_infected    4: phi           5: QER_sum
  6: mask_factor   7: Q             8: P
  9: qer.PBR_qer  10: qer.C_drop   11: qer.d
 12: qer.E        13: qer.Vdrop    14: qer.GVL_ml
 15: qer.GVL_m3   16: qer.VF       17: qer.RTD
 18: qer.DK       19: qer.VER     20: qer.QER_val

### Usage

    python analysis/inputs_reader.py probECAi_inputs.bin

Driver analysis scripts
-----------------------

These scripts read the --save-inputs format to identify which input
parameters drive the high-probability tail.

exceed_vs_nonexceed.py
-----------------------

Splits simulations into exceed (P > 0.1%) and non-exceed (P <= 0.1%)
groups per category. For each input parameter, computes the mean in
each group and the ratio (exceed/non-exceed). Produces a tornado plot
and CSV. Parameters with the largest ratio are the tail drivers.

### Usage

    python analysis/exceed_vs_nonexceed.py probECAi_inputs.bin --csv results.csv --save tornado.png

    # Single category:
    python analysis/exceed_vs_nonexceed.py probECAi_inputs.bin --category Classroom

### Options

| Option       | Default                    | Description                                    |
|--------------|----------------------------|------------------------------------------------|
| `infile`     | `probabilityECAi_inputs.bin` | Inputs .bin file (positional).              |
| `--target`   | `0.1`                      | Target probability (%) for the split.          |
| `--csv`      | *(none)*                   | Write results to this CSV file.               |
| `--save`     | *(none)*                   | Save tornado plot to this file.               |
| `--category` | *(all)*                    | Only analyze this category.                   |

correlation_analysis.py
-----------------------

Computes Spearman rank correlation between each input parameter and
infection probability P. Produces a correlation heatmap (categories
x parameters) and CSV.

### Usage

    python analysis/correlation_analysis.py probECAi_inputs.bin --csv corr.csv --save heatmap.png

### Options

| Option   | Default                    | Description                                    |
|----------|----------------------------|------------------------------------------------|
| `infile` | `probabilityECAi_inputs.bin` | Inputs .bin file (positional).              |
| `--csv`  | *(none)*                   | Write correlation table to CSV.               |
| `--save` | *(none)*                   | Save heatmap to this file.                    |

parameter_scan.py
-----------------

Univariate sensitivity scan: bins P by each parameter value and
computes the median P per bin, showing the marginal effect of each
parameter on infection probability. Produces a multi-panel figure.

### Usage

    python analysis/parameter_scan.py probECAi_inputs.bin --category Classroom --save scan.png

### Options

| Option       | Default                    | Description                                    |
|--------------|----------------------------|------------------------------------------------|
| `infile`     | `probabilityECAi_inputs.bin` | Inputs .bin file (positional).              |
| `--category` | `Classroom`                | Category to analyze.                          |
| `--save`     | *(none)*                   | Save figure to this file.                    |

Generating the raw data
-----------------------

From the `Python/` directory:

    python probabilityECAi.py 1000000 --save-all --outfile run1M.npz

From the repository root via WSL (C++):

    ./C++/build/probability_ecai 1000000 --save-all \
        --outfile probECAi_1M.bin

    # With input parameters for driver analysis:
    ./C++/build/probability_ecai 1000000 --save-all --save-inputs \
        --outfile probECAi_1M.bin \
        --inputs-file probECAi_1M_inputs.bin

All values are stored as raw probability (0-1); the analysis scripts
convert to percent for display.