Analysis Scripts
================

Post-processing tools for the raw simulation data produced by the ASHRAE 241
risk model. They operate on the compressed `.npz` files written by
`probabilityECAi.py --save-all`, which contain one `float32` array of raw
infection probabilities (values in the range 0–1) per occupancy category.

percentiles.py
--------------

Computes percentile statistics from a raw `.npz` file. Produces two reports:

- **table** – per-category percentile table (default percentiles: 50th, 75th, 96th).
- **threshold** – finds the highest percentile at which *all* categories stay
  below a target probability (default 0.1%) using a coarse-to-fine search:
  it scans 25 → 95 in steps of 20, then repeatedly refines within the
  bracketing interval at steps of 10, 5, and 1. Each scan flags the worst-case
  category and the maximum probability across categories.

### Usage

Run from the project root:

    # Both reports on the default file (probabilityECAi_raw.npz)
    python analysis/percentiles.py

    # Point at a specific raw file
    python analysis/percentiles.py run10k.npz

    # Only the 0.1% threshold sweep
    python analysis/percentiles.py --report threshold

    # Custom percentiles in the table
    python analysis/percentiles.py --report table --percentiles 50 90 95 99

    # Custom target for the threshold sweep (in %)
    python analysis/percentiles.py --target 0.05

### Options

| Option           | Default                    | Description                                        |
|------------------|----------------------------|----------------------------------------------------|
| `infile`         | `probabilityECAi_raw.npz`  | Raw `.npz` file to analyze (positional).           |
| `--report`       | `both`                     | `table`, `threshold`, or `both`.                   |
| `--percentiles`  | `50 75 96`                 | Percentiles for the table report.                  |
| `--target`       | `0.1`                      | Target probability (%) for the threshold report.   |

boxplot.py
----------

Draws a horizontal box-and-whisker plot of the raw probabilities, one box per
occupancy category (categories on the y-axis, infection probability on a
**log** x-axis). Boxes show the quartiles with a black median line and
1.5x-IQR whiskers (outlier fliers suppressed for readability). Each category's
96th percentile is marked with a red tick, along with a dashed reference line
at the target probability.

### Usage

Run from the project root:

    # Plot the default file (probabilityECAi_raw.npz)
    python analysis/boxplot.py

    # Plot a specific raw file
    python analysis/boxplot.py run10k.npz

    # Custom target line and save to a file instead of showing
    python analysis/boxplot.py --target 0.05 --save boxplot.png

### Options

| Option      | Default                    | Description                                          |
|-------------|----------------------------|------------------------------------------------------|
| `infile`    | `probabilityECAi_raw.npz`  | Raw `.npz` file to plot (positional).                |
| `--target`  | `0.1`                      | Target probability (%) drawn as a reference line.    |
| `--xmin`    | `0.001`                    | Left edge of the log x-axis in % (clips lower tail). |
| `--xmax`    | `100`                      | Right edge of the log x-axis in % (clips upper tail).|
| `--save`    | *(show interactively)*     | Save the figure to this path instead of displaying.  |

Generating the raw data
-----------------------

From the project root:

    python probabilityECAi.py 10000 --save-all --outfile run10k.npz

All values are stored as raw probability (0–1); the analysis scripts convert to
percent for display.
