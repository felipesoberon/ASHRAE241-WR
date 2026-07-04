C++ Port of ASHRAE 241 Risk Model
==================================

High-performance C++17 port of the Python ASHRAE 241 Wells-Riley Monte
Carlo simulation. Produces statistically equivalent results with a
~20x speedup over the Python implementation.

Build
-----

Requirements: C++17 compiler (g++ 13+) and CMake 3.16+.

From WSL or Linux, in the C++/ directory:

    cmake -B build
    cmake --build build

Run tests:

    cd build && ctest -V

Usage
-----

All executables are in C++/build/. They accept the same arguments as
their Python counterparts and produce the same output format.

ECAi simulation (compute required ECAi at 0.1% target):

    ./build/ecai [N]

Infection probability at ASHRAE 241 ECAi values:

    ./build/probability_ecai [N] [--save-all] [--outfile FILE]

Scan for minimum ECAi:

    ./build/probability_scan [N]

Single category:

    ./build/single_probability --category Classroom --N 10000

Analysis (percentile table + threshold search):

    ./build/percentiles [infile] [--report both|table|threshold] [--target 0.1]

Box plot (text-based):

    ./build/boxplot [infile] [--target 0.1]

Benchmark Results (N=10000, all 25 categories)
----------------------------------------------

| Script             | Python       | C++        | Speedup   |
|--------------------|--------------|------------|-----------|
| ecai               | 169s         | 8.4s       | ~20x      |
| probability_ecai   | 300s+ (timed | 12.8s      | ~23x+     |
|                    | out at 300s) |            |           |

Validation (N=100000, 96th percentile)
--------------------------------------

ECAi simulation: 21/25 categories within 5% relative difference.
The 4 outliers (Manufacturing, Office, Transportation, Warehouse)
have very small 96th percentile ECAi values (< 8 L/s/p), so small
absolute differences produce large relative differences. All are
within Monte Carlo error for N=100k.

Infection probability: 22/25 categories within 5% relative
difference. The 3 near-misses (Lecture 5.0%, Museum 5.2%, Place
5.5%) are small-probability categories where proportional variance
is higher. All converge at larger N.

Architecture
------------

  src/random_manager.{h,cpp}  LHS engine + inverse CDFs
  src/model.{h,cpp}           Occupancy params + QER/infection/ECAi
  src/ecai.cpp               ECAi simulation script
  src/probability_ecai.cpp   Probability at ASHRAE ECAi values
  src/probability_scan.cpp   ECAi threshold scan
  src/single_probability.cpp Single-category analysis
  analysis/percentiles.cpp   Percentile table + threshold search
  analysis/boxplot.cpp        Text-based box plot
  tests/                     Unit tests + cross-validation

Binary raw data format (--save-all)
------------------------------------

C++ uses a simple binary format (not .npz):

  uint32_t magic = 0x41534852
  uint32_t num_categories
  Per category: uint32_t name_len, char[] name, uint32_t count,
                float[count] values (raw probability 0-1)

The C++ analysis tools (percentiles, boxplot) read this format.