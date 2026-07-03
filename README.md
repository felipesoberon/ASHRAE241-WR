ASHRAE 241 Risk Model Simulation
================================

This repository contains Python tools for running Monte Carlo simulations to estimate airborne infection risk in indoor environments using the ASHRAE Standard 241 framework. The model implemented here is based on the methodology reported by Jones et al. (2025) in *Building and Environment*, which estimates infection risk as a function of equivalent clean airflow rates (ECAi), viral emission characteristics, occupant behavior, and space properties. These tools reproduce and extend the simulations described in that publication.


Features
--------

- Simulates infection probability using a mechanistic Wells-Riley-based approach.
- Supports a wide range of occupancy categories (Classroom, Gym, Office, Warehouse, etc.).
- Community prevalence (infection rate) is incorporated into all simulations.
- Identifies the minimum ECAi where the 96th percentile of infection risk drops below 0.1%.
- Implements **Latin Hypercube Sampling (LHS) per distribution** for improved random coverage and lower variance.
- Modular design for easy validation and extension.
- Generates histograms for input parameter distributions for sanity checking.

Requirements
------------

- Python 3.x
- numpy
- matplotlib
- scipy

To install dependencies:

    pip install numpy matplotlib scipy

Program Overview
----------------

**model.py**  
- Contains the **core functions** and definitions for the simulation model, including occupancy parameters, random sampling of variables, and calculations of quanta emission rate (QER), infection probability, and equivalent clean airflow (ECAi).
- All simulation scripts rely on this module to implement the methodology from the Jones et al. (2025) publication.


**qer.py**  
- Samples and plots the distributions of all parameters used in the Quanta Emission Rate (QER) calculation for a selected occupancy category.
- Default: `Classroom`, or specify a category at the command line.
- *Example:*
    ```
    python qer.py Office
    ```

**distributions.py**  
- Generates and plots random samples for all distribution types (Normal, Beta, Lognormal, Uniform, etc.) for a selected occupancy.
- Useful for sanity-checking the range and distribution of model inputs.
- Default: `Classroom`, or specify a category.

**randomManager.py**  
- Provides the LHS random number engine with independent buffers for each distribution type.
- All other scripts use this module for low-bias, stratified sampling of random variables.

**ecai.py**  
- Main simulation script: For each occupancy, runs 10,000 (default) simulations at a fixed infection probability target (0.1%).
- Determines the required ECAi, reports the 96th percentile, rounded values, and the percentage of simulations with zero infected.
- Results are printed as a compact table and saved to CSV.
- *Example:*
    ```
    python ecai.py
    ```
    (Optionally specify number of runs, e.g. `python ecai.py 5000`)

**probabilityScan.py**  
- Replaces the old `probability.py`.
- Calculates the infection probability over 10,000 runs for a range of ECAi values (5, 10, ... 50 L/s/person) for each occupancy.
- Reports the minimum ECAi at which the 96th percentile of infection probability is below 0.1%.
- Results are saved to a CSV file.

**probabilityECAi.py**  
- Uses the ECAi values specified in the ASHRAE Standard 241 for each occupancy category.
- Determines the 96th percentile infection probability for each occupancy, including **only simulation runs where there is at least one infector present** (i.e., if the random draw results in zero infectors, it redraws until at least one is present).
- This method answers:
  _“If I spend an hour in this location and there is at least one infected person present, what is the probability of catching the infection given the location is compliant with ASHRAE 241?”_
- Results are printed in a table and saved to CSV.
- Optionally saves **every** raw simulation probability (not just the 96th percentile) to a compressed `.npz` file for later analysis (e.g. computing other percentiles or the mean without re-running the model).
- *Options:*
  - `N` (positional) — number of simulations per category (default 10000).
  - `--save-all` — save all raw probabilities to a `.npz` file (one `float32` array per category, values 0–1).
  - `--outfile` — name of the raw-data file (default `probabilityECAi_raw.npz`).
- *Example:*
    ```
    python probabilityECAi.py 1000000 --save-all --outfile run1M.npz
    ```

**singleProbability.py**  
- Runs infection probability simulations for a **single occupancy category**.
- Reports the 96th percentile infection probability and the percentage of simulations with zero infectors.
- Optionally plots linear and log-scale histograms of the infection probability distribution.
- Includes command-line options to specify occupancy category, number of simulations, community infection rate, ECAi override, and plot display.
- For help on available options, run:
    ```
    python singleProbability.py --help
    ```
- *Example:*
    ```
    python singleProbability.py --category Classroom --N 10000 --show_plots
    ```


**mainGUI.py**  
- Provides a **Tkinter-based desktop GUI** for running both ECAi and infection probability simulations.
- Allows selection of calculation type, number of simulations, and whether to include zero-infector scenarios.
- Displays results as horizontal bar plots comparing simulated values to ASHRAE 241 values.

**streamlitGUI.py**  
- Provides a **Streamlit web-based GUI** for interactive simulations.
- Allows running ECAi or infection probability simulations across all occupancy categories.
- Displays results with horizontal bar plots and log-scale options for infection probability.
- To run:
    ```
    streamlit run streamlitGUI.py
    ```

Analysis
--------

The `analysis/` folder contains post-processing tools that operate on the raw
`.npz` files saved by `probabilityECAi.py --save-all`.

**analysis/percentiles.py**
- Computes percentile statistics from a raw `.npz` file.
- Produces a per-category percentile table and a coarse-to-fine threshold search
  that finds the highest percentile at which all categories stay below a target
  probability (default 0.1%).
- *Example:*
    ```
    python analysis/percentiles.py run1M.npz
    ```
- See `analysis/README.md` for full options.

Usage
-----

To simulate for a specific occupancy (e.g., Classroom), run:

    python ecai.py Classroom

To generate and check parameter distributions for a different category:

    python qer.py Office

For other script usage, see the in-line help or comments in each file.

Output
------

- Terminal output summarizing 96th percentile, rounded ECAi, CFM equivalents, and zero-infected percentages.
- Plots (when relevant) for each parameter or outcome.
- Summary tables saved as CSV files.

Simulation Details
------------------

- **Random sampling:** All random draws use Latin Hypercube Sampling with a separate buffer for each distribution type, ensuring independent stratification and even coverage for each random variable.
- **Performance:** Batch size for LHS buffers is configurable (default: 1,000,000).
- **Sanity-checking:** Parameter distributions can be visualized with `qer.py` or `distributions.py`.

Reference
---------

This simulation is based on the methodology in:

Benjamin Jones, Christopher Iddon, Marwa Zaatari, Pawel Wargocki, Richard Bruns.  
**Risk modeling for ASHRAE Standard 241-2023 – Control of infectious aerosols**.  
*Building and Environment*, Volume 248, 2025, 113318.  
https://doi.org/10.1016/j.buildenv.2025.113318
