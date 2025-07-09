ASHRAE 241 Risk Model Simulation
================================

This repository contains Python tools for running Monte Carlo simulations to estimate airborne infection risk in indoor environments using the ASHRAE Standard 241 framework. The models estimate infection risk as a function of equivalent clean airflow rates (ECAi), viral emission characteristics, occupant behavior, and space properties.

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
