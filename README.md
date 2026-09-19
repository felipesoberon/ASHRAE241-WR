ASHRAE 241 Risk Model Simulation
================================

Monte Carlo simulation of airborne infection risk in indoor
environments using the ASHRAE Standard 241 Wells-Riley framework.
Based on the methodology in Jones et al. (2025).

Implementations
----------------

**C++/** -- Current C++17 implementation: an aggregate-demand-conditioned
marginal 1D-LHS ECAi model (executable `ecai_cpp2`). See `C++/README.md`
for scope, isolation, and validation status.

**C++_legacy/** -- Original C++17 port of the Python model (~20x faster).
Superseded by `C++/` as the actively maintained implementation but retained
for reference and comparison. See `C++_legacy/README.md` for build
instructions, benchmarks, and validation results.

**Python/** -- Reference implementation (numpy, scipy, matplotlib).
See `Python/README.md` for full documentation.


Reference
---------

Benjamin Jones, Christopher Iddon, Marwa Zaatari, Pawel Wargocki,
Richard Bruns.
**Risk modeling for ASHRAE Standard 241-2023 -- Control of
infectious aerosols**.
*Building and Environment*, Volume 248, 2025, 113318.
https://doi.org/10.1016/j.buildenv.2025.113318