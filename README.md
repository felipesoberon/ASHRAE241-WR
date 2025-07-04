ASHRAE 241 Risk Model Simulation
================================

This repository contains a Monte Carlo simulation to estimate airborne infection risk in indoor environments using the ASHRAE Standard 241 framework. It models quanta generation, airflow, deposition, and removal rates to compute infection probabilities across various occupancy categories.

Features
--------

- Simulates infection probability using a mechanistic model
- Supports a wide range of occupancy types (e.g., Classroom, Gym, Office)
- Includes variable parameters: breathing rate, particle sizes, viral load, ventilation
- Automatically stops when the 96th percentile of infection probability drops below 0.1%
- Generates a log-scale histogram of infection probability from the final run

Requirements
------------

- Python 3.x
- numpy
- matplotlib

Install dependencies with:

    pip install numpy matplotlib

Usage
-----

Run the simulation with the desired occupancy category:

    python main.py Classroom

Replace `Classroom` with any valid category such as `Gym`, `Office`, `Warehouse`, etc.

Output
------

- Console output showing infection probability statistics per ECAi (L/s/person)
- Stops iterating once the 96th percentile is < 0.1%
- A plot showing the histogram of infection probabilities (log scale)

Key Variables
-------------

- ECAi: Equivalent clean air per person (L/s/person)
- I0: Number of people in the space
- community_rate: Community infection prevalence (0–1)
- QER: Quanta emission rate (quanta/hour)
- P: Infection probability (dimensionless)

Simulation Logic
----------------

- Each simulation draws random samples from biological and physical distributions
- Number of infectious individuals is simulated using community_rate as a threshold
- Model includes deposition, biological decay, ventilation, and mask effectiveness

Reference
---------

This simulation approach is based on the methodology published in:

Peng, Z. et al., 2025. A mechanistic framework for evaluating air-cleaning technologies under ASHRAE Standard 241. *Building and Environment*, 248, 113318. https://doi.org/10.1016/j.buildenv.2025.113318
