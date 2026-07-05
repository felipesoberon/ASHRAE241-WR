ASHRAE 241 Risk Model Simulation
================================

Monte Carlo simulation of airborne infection risk in indoor
environments using the ASHRAE Standard 241 Wells-Riley framework.
Based on the methodology in Jones et al. (2025).

Implementations
----------------

**Python/** -- Reference implementation (numpy, scipy, matplotlib).
See `Python/README.md` for full documentation.

**C++/** -- High-performance C++17 port (~20x faster).
See `C++/README.md` for build instructions, benchmarks, and
validation results.


Reference
---------

Benjamin Jones, Christopher Iddon, Marwa Zaatari, Pawel Wargocki,
Richard Bruns.
**Risk modeling for ASHRAE Standard 241-2023 -- Control of
infectious aerosols**.
*Building and Environment*, Volume 248, 2025, 113318.
https://doi.org/10.1016/j.buildenv.2025.113318