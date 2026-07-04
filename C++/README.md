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

All executables are in C++/build/. They mirror their Python counterparts'
arguments and output format, with a few added flags (below) that the Python
scripts do not expose.

ECAi simulation (compute required ECAi at 0.1% target):

    ./build/ecai [N] [--require-infectors]
                     [--community-rate-general RATE]
                     [--community-rate-healthcare RATE]

Infection probability at ASHRAE 241 ECAi values:

    ./build/probability_ecai [N] [--save-all] [--outfile FILE]
                             [--no-require-infectors]
                             [--community-rate-general RATE]
                             [--community-rate-healthcare RATE]

Scan for minimum ECAi:

    ./build/probability_scan [N] [--community-rate-general RATE]
                                 [--community-rate-healthcare RATE]

Community infection rate options: `--community-rate-general` sets the rate for
general spaces (default 1%) and `--community-rate-healthcare` for healthcare
spaces (Exam, Group, Patient, Resident, Waiting; default 3%). When a flag is
omitted, that group keeps its per-category ASHRAE 241 default, so results are
unchanged from the original model.

Single category:

    ./build/single_probability --category Classroom --N 10000

Analysis (percentile table + threshold search):

    ./build/percentiles [infile] [--report both|table|threshold] [--target 0.1]

Box plot (text-based):

    ./build/boxplot [infile] [--target 0.1]

Graphical interface (C++-backed)
--------------------------------

mainGUI_C++.py is a Python/Tkinter GUI that mirrors the original mainGUI.py but
runs the C++ engine instead of the Python model. It is a Windows Python app that
invokes the ecai / probability_ecai executables through WSL.

    python mainGUI_C++.py      (run from the C++/ directory)

Requirements:
  - WSL with a Linux distro (auto-detected).
  - The C++ engine already built in C++/build/ (see Build above). The GUI shows
    a build hint if the executables are missing.
  - Windows Python with tkinter, matplotlib, and numpy.

Controls mirror mainGUI.py, plus separate community-rate entries for general
(default 1%) and healthcare (default 3%) spaces, which map to the
--community-rate-general / --community-rate-healthcare engine flags.

Benchmark Results (N=10000, all 25 categories)
----------------------------------------------

| Script             | Python       | C++        | Speedup   |
|--------------------|--------------|------------|-----------|
| ecai               | 169s         | 8.4s       | ~20x      |
| probability_ecai   | 300s+ (timed | 12.8s      | ~23x+     |
|                    | out at 300s) |            |           |

Validation (N=1,000,000, 96th percentile, vs stored Python reference)
---------------------------------------------------------------------

Infection probability: 25/25 categories within 1.4% relative
difference (mean ~0.6%). This is the robust comparison and is
essentially exact.

ECAi simulation: every category with a meaningful value (> 5 L/s/p)
agrees within 3.6%. The only larger relative gaps are Patient (33.8%)
and Warehouse (9.5%) -- near-zero categories (82-91% zero-infector
sims) where the 96th percentile is ~0-1 L/s/p, so the absolute
difference is < 0.25 L/s/p and rounds to the same value. This is
Monte Carlo noise on a near-zero quantity, not a model discrepancy.

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
  mainGUI_C++.py             Tkinter GUI (C++ engine via WSL)
  tests/                     Unit tests + cross-validation

Binary raw data format (--save-all)
------------------------------------

C++ uses a simple binary format (not .npz):

  uint32_t magic = 0x41534852
  uint32_t num_categories
  Per category: uint32_t name_len, char[] name, uint32_t count,
                float[count] values (raw probability 0-1)

The C++ analysis tools (percentiles, boxplot) read this format.