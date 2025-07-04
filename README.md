ASHRAE 241 Risk Model Simulation
================================

This repository contains a Monte Carlo simulation to estimate airborne infection risk in indoor environments using the framework from ASHRAE Standard 241. The model evaluates the infection risk as a function of equivalent clean airflow rates (ECAi), viral emission characteristics, occupant behavior, and space properties.

Features
--------

- Simulates infection probability using a mechanistic Wells-Riley-based approach.
- Supports various occupancy categories: Classroom, Gym, Office, Warehouse, etc.
- Community prevalence (infection rate) incorporated into simulations.
- Runs until the 96th percentile of infection risk drops below 0.1%.
- Generates a histogram (log scale) of infection probability from the final iteration.

Requirements
------------

- Python 3.x
- numpy
- matplotlib

To install dependencies:

    pip install numpy matplotlib

Usage
-----

Run the simulation from the command line, providing the occupancy category as an argument:

    python main.py Classroom

Replace `Classroom` with other valid categories such as `Office`, `Gym`, `Warehouse`, etc.

Output
------

- Terminal output showing the mean, median, min, max, and 96th percentile infection probability for each ECAi value.
- Plot showing the histogram of infection probabilities (logarithmic y-axis).
- Simulation automatically stops once the 96th percentile drops below 0.1%.

Reference
---------

The simulation is based on the methodology described in the following publication:

Benjamin Jones, Christopher Iddon, Marwa Zaatari, Pawel Wargocki, Richard Bruns.  
**Risk modeling for ASHRAE Standard 241-2023 – Control of infectious aerosols**.  
*Building and Environment*, Volume 248, 2025, 113318.  
https://doi.org/10.1016/j.buildenv.2025.113318
