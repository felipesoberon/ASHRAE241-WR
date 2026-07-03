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

Generating the raw data
-----------------------

From the project root:

    python probabilityECAi.py 10000 --save-all --outfile run10k.npz

All values are stored as raw probability (0–1); the analysis scripts convert to
percent for display.
