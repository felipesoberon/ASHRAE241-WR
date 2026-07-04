# C++ Port of ASHRAE 241 Risk Model — Implementation Plan

> **For Hermes:** Use subagent-driven-development skill to implement this plan task-by-task.

**Goal:** Port the Python ASHRAE 241 Wells-Riley Monte Carlo simulation to C++ for significant speed improvement, producing statistically identical results.

**Architecture:** A standalone C++ project inside `C++/` that mirrors the Python module structure: a LHS random manager, a model module with occupancy parameters and the four core functions (sample_parameters, QER, infection_probability, compute_ECAi), and command-line simulation scripts that produce the same terminal output and CSV files as their Python counterparts. No external dependencies — all inverse CDFs implemented from scratch. Build system: CMake. The Python code is never modified; it serves as the reference implementation for validation.

**Tech Stack:** C++17, CMake, standard library only (no Boost, no scipy). Optional: miniz (single-header ZIP) for .npz output compatibility.

---

## Constraints

1. **Python model is frozen.** No edits to any .py file. It is the reference.
2. **Statistical equivalence, not bit-identical.** C++ uses its own RNG for LHS buffers, so individual draws differ, but distributions, stratification logic, formulas, and parameter values must match exactly. At the same N, 96th percentile values should agree within Monte Carlo error (~1-3% relative at N=10000).
3. **Output compatibility.** C++ scripts write the same CSV column names and terminal table format as Python. For `--save-all`, C++ writes a simple binary format (.bin) that C++ analysis tools read. .npz compatibility is a later optional task.
4. **Platform:** Windows 10 with MSYS/Git Bash. Must build with g++ (MSYS) or MSVC. CMake handles both.

---

## Project Structure

```
C++/
├── CMakeLists.txt
├── README.md
├── src/
│   ├── random_manager.h          # LHS engine class + distribution samplers
│   ├── random_manager.cpp
│   ├── model.h                    # occupancy_params, sample_parameters, QER,
│   │                              #   infection_probability, compute_ECAi
│   ├── model.cpp
│   ├── ecai.cpp                   # port of ecai.py
│   ├── probability_ecai.cpp       # port of probabilityECAi.py
│   ├── probability_scan.cpp       # port of probabilityScan.py
│   └── single_probability.cpp    # port of singleProbability.py
├── analysis/
│   ├── percentiles.cpp            # port of analysis/percentiles.py
│   └── boxplot.cpp                # port of analysis/boxplot.py (text summary;
│                                  #   PNG output optional via gnuplot)
├── tests/
│   ├── test_random_manager.cpp    # unit tests for inverse CDFs + LHS
│   ├── test_model.cpp             # unit tests for QER, infection_probability
│   └── test_validation.cpp        # integration: compare C++ vs Python CSV output
└── third_party/
    └── (miniz.h — optional, for .npz output, added later)
```

---

## Phase 1: Setup

### Task 1: Create directory structure and CMakeLists.txt

**Objective:** Scaffold the C++ project with CMake build system.

**Files:**
- Create: `C++/CMakeLists.txt`
- Create: `C++/README.md` (placeholder, fill at end)

**Step 1: Write CMakeLists.txt**

```cmake
cmake_minimum_required(VERSION 3.16)
project(ASHRAE241_WR LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Warnings as errors for safety
if(MSVC)
    add_compile_options(/W4 /O2)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -O2)
endif()

# Core library (random_manager + model)
add_library(ashrae241 STATIC
    src/random_manager.cpp
    src/model.cpp
)
target_include_directories(ashrae241 PUBLIC src)

# Simulation executables
add_executable(ecai            src/ecai.cpp)
add_executable(probability_ecai  src/probability_ecai.cpp)
add_executable(probability_scan src/probability_scan.cpp)
add_executable(single_probability src/single_probability.cpp)

foreach(tgt ecai probability_ecai probability_scan single_probability)
    target_link_libraries(${tgt} PRIVATE ashrae241)
endforeach()

# Analysis executables
add_executable(percentiles  analysis/percentiles.cpp)
add_executable(boxplot      analysis/boxplot.cpp)

# Tests
enable_testing()
add_executable(test_random_manager tests/test_random_manager.cpp)
target_link_libraries(test_random_manager PRIVATE ashrae241)
add_test(NAME random_manager COMMAND test_random_manager)

add_executable(test_model tests/test_model.cpp)
target_link_libraries(test_model PRIVATE ashrae241)
add_test(NAME model COMMAND test_model)

add_executable(test_validation tests/test_validation.cpp)
target_link_libraries(test_validation PRIVATE ashrae241)
add_test(NAME validation COMMAND test_validation)
```

**Step 2: Create empty README.md placeholder**

**Step 3: Verify CMake configures**

Run: `cd C++ && cmake -B build && cmake --build build`
Expected: Build fails (no source files yet) — that's fine, CMake configures successfully.

**Step 4: Commit**

```bash
git add C++/
git commit -m "feat: scaffold C++ project with CMake"
```

---

### Task 2: Verify compiler availability

**Objective:** Confirm a C++17 compiler and CMake are available on this Windows machine.

**Step 1: Check for compilers**

Run: `g++ --version 2>/dev/null; cmake --version 2>/dev/null`
If g++ is missing, check MSVC: `where cl.exe 2>/dev/null`

**Step 2: If no compiler found, install one**

Options (in order of preference for MSYS/Git Bash environment):
- MSYS2 pacman: `pacman -S mingw-w64-x86_64-gcc cmake`
- Or install Visual Studio Build Tools with C++ workload

**Step 3: Verify a trivial build**

```bash
echo '#include <cstdio>
int main() { printf("ok\n"); return 0; }' > /tmp/t.cpp
g++ -std=c++17 /tmp/t.cpp -o /tmp/t && /tmp/t
```

Expected: prints `ok`

---

## Phase 2: LHS Random Manager

This is the most critical piece — it must match the Python LHS stratification exactly.

### Task 3: Implement inverse normal CDF

**Objective:** Provide `normal_ppf(u)` returning the same values as `scipy.stats.norm.ppf(u)` for u in (0,1).

**Files:**
- Create: `C++/src/random_manager.h`
- Create: `C++/src/random_manager.cpp`
- Create: `C++/tests/test_random_manager.cpp`

**Step 1: Write the header skeleton**

`C++/src/random_manager.h`:
```cpp
#pragma once
#include <vector>
#include <cstdint>

// Inverse CDFs (percent point functions) — match scipy.stats.*.ppf
double normal_ppf(double u);       // scipy.stats.norm.ppf(u)
double lognormal_ppf(double u, double mean, double sigma);  // lognorm.ppf(u, s=sigma, scale=exp(mean))
double uniform_ppf(double u, double low, double high);
double beta_ppf(double u, double a, double b);  // scipy.stats.beta.ppf(u, a, b) on [0,1]
double log10normal_ppf(double u, double mu, double sigma);  // 10^normal_ppf(u, mu, sigma)
int    binomial_ppf(double u, int n, double p);  // inverse CDF of Binomial(n, p)
```

**Step 2: Implement normal_ppf using Acklam's algorithm**

```cpp
// Acklam's inverse normal CDF approximation (max rel error ~1.15e-9)
double normal_ppf(double u) {
    // Constants
    static const double a[] = {
        -3.969683028665376e+01, 2.209460984245205e+02,
        -2.759285188461852e+02, 1.383577518672690e+02,
        -3.066479806614716e+01, 2.506628277453239e+00
    };
    static const double b[] = {
        -5.447609879822406e+01, 1.615858368580409e+02,
        -1.556989798598866e+02, 6.680131188771972e+01,
        -1.328068155288572e+01
    };
    static const double c[] = {
        -7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00,
        4.374664141118049e+00, 2.938163982698783e+00
    };
    static const double d[] = {
        7.784695709041462e-03, 3.224671290700398e-01,
        2.445134137142996e+00, 3.754408661907416e+00
    };
    double plow = 0.02425, phigh = 1 - plow;
    double q, r, x;
    if (u < plow) {
        q = sqrt(-2 * log(u));
        x = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
            ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1);
    } else if (u <= phigh) {
        q = u - 0.5;
        r = q * q;
        x = (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q /
            (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1);
    } else {
        q = sqrt(-2 * log(1 - u));
        x = -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
            ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1);
    }
    // One Halley refinement step using the error function
    double e = 0.5 * erfc(-x / sqrt(2.0)) - u;
    double u1 = e * sqrt(2.0 * M_PI) * exp(x*x / 2.0);
    x = x - u1 / (1 + x * u1 / 2);
    return x;
}
```

**Step 3: Implement the trivial wrappers**

```cpp
double lognormal_ppf(double u, double mean, double sigma) {
    return exp(normal_ppf(u) * sigma + mean);  // lognorm with s=sigma, scale=exp(mean)
}

double uniform_ppf(double u, double low, double high) {
    return low + (high - low) * u;
}

double log10normal_ppf(double u, double mu, double sigma) {
    return pow(10.0, normal_ppf(u) * sigma + mu);
}
```

**Step 4: Write test — verify against known scipy values**

```cpp
// tests/test_random_manager.cpp (initial)
#include "random_manager.h"
#include <cstdio>
#include <cmath>

int main() {
    // scipy.stats.norm.ppf(0.5) == 0.0
    if (fabs(normal_ppf(0.5)) > 1e-9) { printf("FAIL: normal_ppf(0.5)\n"); return 1; }
    // scipy.stats.norm.ppf(0.025) ~= -1.959963985
    if (fabs(normal_ppf(0.025) - (-1.959963985)) > 1e-6) { printf("FAIL: normal_ppf(0.025) = %.10f\n", normal_ppf(0.025)); return 1; }
    // scipy.stats.norm.ppf(0.975) ~= 1.959963985
    if (fabs(normal_ppf(0.975) - 1.959963985) > 1e-6) { printf("FAIL: normal_ppf(0.975) = %.10f\n", normal_ppf(0.975)); return 1; }
    // scipy.stats.norm.ppf(0.001) ~= -3.090232306
    if (fabs(normal_ppf(0.001) - (-3.090232306)) > 1e-6) { printf("FAIL: normal_ppf(0.001)\n"); return 1; }

    // lognormal: lognorm.ppf(0.5, s=0.1, scale=exp(0.3)) = exp(0.3) = 1.3498...
    double ln = lognormal_ppf(0.5, 0.3, 0.1);
    if (fabs(ln - exp(0.3)) > 1e-6) { printf("FAIL: lognormal_ppf\n"); return 1; }

    // uniform: uniform_ppf(0.3, 2, 5) = 2 + 0.3*3 = 2.9
    if (fabs(uniform_ppf(0.3, 2.0, 5.0) - 2.9) > 1e-12) { printf("FAIL: uniform_ppf\n"); return 1; }

    // log10normal: 10^normal_ppf(0.5, 7.0, 1.4) = 10^7.0 = 1e7
    double l1n = log10normal_ppf(0.5, 7.0, 1.4);
    if (fabs(l1n - 1e7) / 1e7 > 1e-9) { printf("FAIL: log10normal_ppf = %f\n", l1n); return 1; }

    printf("PASS: inverse CDFs\n");
    return 0;
}
```

**Step 5: Run test to verify pass**

Run: `cd C++ && cmake -B build && cmake --build build && cd build && ctest -V -R random_manager`
Expected: PASS

**Step 6: Commit**

```bash
git add C++/
git commit -m "feat: implement inverse normal CDF and trivial wrappers"
```

---

### Task 4: Implement inverse beta CDF

**Objective:** Provide `beta_ppf(u, a, b)` matching `scipy.stats.beta.ppf(u, a, b)`.

This is the hardest inverse CDF. We need the incomplete beta function and its inverse.

**Files:**
- Modify: `C++/src/random_manager.cpp`

**Step 1: Implement the incomplete beta function (continued fraction)**

```cpp
// Incomplete beta function I_x(a,b) via continued fraction (Lentz's algorithm)
// Based on Numerical Recipes betacf
static double betacf(double a, double b, double x) {
    const int MAXIT = 200;
    const double EPS = 3e-16;
    const double FPMIN = 1e-300;
    double qab = a + b, qap = a + 1, qam = a - 1;
    double c = 1, d = 1 - qab * x / qap;
    if (fabs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    double h = d;
    for (int m = 1, m2 = 2; m <= MAXIT; m++, m2 += 2) {
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1 + aa * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1 + aa / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1 + aa * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1 + aa / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (fabs(del - 1) < EPS) break;
    }
    return h;
}

// Regularized incomplete beta function I_x(a,b) = P(X <= x) for Beta(a,b)
static double beta_cdf(double x, double a, double b) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    double bt = exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * log(x) + b * log(1 - x));
    if (x < (a + 1) / (a + b + 2))
        return bt * betacf(a, b, x) / a;
    else
        return 1 - bt * betacf(b, a, 1 - x) / b;
}
```

**Step 2: Implement beta_ppf via bisection**

```cpp
double beta_ppf(double u, double a, double b) {
    if (u <= 0) return 0;
    if (u >= 1) return 1;
    // Bisection on beta_cdf — robust, ~50 iterations for double precision
    double lo = 0, hi = 1, mid;
    for (int i = 0; i < 60; i++) {
        mid = 0.5 * (lo + hi);
        double cdf = beta_cdf(mid, a, b);
        if (cdf < u) lo = mid;
        else hi = mid;
    }
    return 0.5 * (lo + hi);
}
```

**Step 3: Add tests for beta_ppf**

Append to `tests/test_random_manager.cpp`:
```cpp
// scipy.stats.beta.ppf(0.5, 5, 2) ~= 0.73557
double bp = beta_ppf(0.5, 5.0, 2.0);
if (fabs(bp - 0.73557) > 1e-5) { printf("FAIL: beta_ppf(0.5,5,2) = %.8f\n", bp); return 1; }
// scipy.stats.beta.ppf(0.1, 5, 2) ~= 0.49015
bp = beta_ppf(0.1, 5.0, 2.0);
if (fabs(bp - 0.49015) > 1e-5) { printf("FAIL: beta_ppf(0.1,5,2) = %.8f\n", bp); return 1; }
// beta_ppf with loc/scale: beta(a,b,loc=L,scale=S) = beta_ppf(u,a,b)*S + L
// Python: random_beta_lhs(rng, a=5.0, b=2.0, loc=2.0, scale=3.0) => range [2, 5]
// median = beta_ppf(0.5,5,2)*3 + 2 = 0.73557*3 + 2 = 4.2067
double be = beta_ppf(0.5, 5.0, 2.0) * 3.0 + 2.0;
if (fabs(be - 4.20671) > 1e-4) { printf("FAIL: beta with loc/scale = %.8f\n", be); return 1; }
```

**Step 4: Run test**

Run: `cd C++/build && cmake --build . && ctest -V -R random_manager`
Expected: PASS

**Step 5: Commit**

```bash
git add C++/
git commit -m "feat: implement inverse beta CDF via incomplete beta function"
```

---

### Task 5: Implement binomial inverse CDF

**Objective:** Provide `binomial_ppf(u, n, p)` matching the Python `random_binomial_lhs` behavior (sequential CDF search).

**Files:**
- Modify: `C++/src/random_manager.cpp`

**Step 1: Implement binomial_ppf**

The Python code sums PMF values sequentially. We do the same but in C++ (much faster). We use `log(lgamma)` for the PMF to avoid overflow with large n.

```cpp
int binomial_ppf(double u, int n, double p) {
    if (u <= 0) return 0;
    if (u >= 1) return n;
    // Match Python: sum PMF from k=0 to n, return first k where cumulative >= u
    double cdf = 0.0;
    double log_p = log(p);
    double log_1mp = log(1.0 - p);
    for (int k = 0; k <= n; k++) {
        double log_pmf = lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
                         + k * log_p + (n - k) * log_1mp;
        cdf += exp(log_pmf);
        if (u <= cdf)
            return k;
    }
    return n;
}
```

**Step 2: Add tests**

```cpp
// Binomial(10, 0.3): P(X<=2) = 0.38278, P(X<=3) = 0.64961
// So binomial_ppf(0.4, 10, 0.3) = 3 (first k where cdf >= 0.4)
int bk = binomial_ppf(0.4, 10, 0.3);
if (bk != 3) { printf("FAIL: binomial_ppf(0.4,10,0.3) = %d\n", bk); return 1; }
// Edge: u=0 => 0
if (binomial_ppf(0.0, 10, 0.3) != 0) { printf("FAIL: binomial_ppf(0,...)\n"); return 1; }
// Edge: u=1 => n
if (binomial_ppf(1.0, 10, 0.3) != 10) { printf("FAIL: binomial_ppf(1,...)\n"); return 1; }
```

**Step 3: Run test**

Expected: PASS

**Step 4: Commit**

```bash
git add C++/
git commit -m "feat: implement binomial inverse CDF"
```

---

### Task 6: Implement LHS buffer manager

**Objective:** Port `RandomNumberManager` class with per-distribution stratified buffers, matching the Python LHS logic exactly.

**Files:**
- Modify: `C++/src/random_manager.h`
- Modify: `C++/src/random_manager.cpp`

**Step 1: Add the class to the header**

```cpp
class RandomNumberManager {
public:
    explicit RandomNumberManager(size_t batch_size = 1000000);

    // Returns one stratified uniform from the buffer for `dist`
    double get(const std::string& dist);

private:
    size_t batch_size_;
    std::map<std::string, std::vector<double>> buffers_;
    std::map<std::string, size_t> indices_;
    std::mt19937_64 rng_;  // Mersenne Twister for buffer generation

    void refill(const std::string& dist);
};

// Distribution sampling functions (match Python API)
double random_beta_lhs(RandomNumberManager& rng, double a, double b,
                       double loc = 0, double scale = 1);
double random_normal_lhs(RandomNumberManager& rng, double mu, double sigma);
double random_lognormal_lhs(RandomNumberManager& rng, double mean, double sigma);
double random_uniform_lhs(RandomNumberManager& rng, double low, double high);
int    random_binomial_lhs(RandomNumberManager& rng, int n, double p);
double random_log10normal_lhs(RandomNumberManager& rng, double mu, double sigma);
```

**Step 2: Implement the class**

The key: match Python's `_refill` exactly:
1. Create `batch_size + 1` evenly spaced edges on [0,1] (`np.linspace(0, 1, batch_size + 1)`)
2. Draw one uniform from each segment `[edge[i], edge[i+1])`
3. Shuffle the buffer
4. Serve values sequentially; refill when exhausted

```cpp
#include <random>
#include <algorithm>
#include <map>
#include <string>

RandomNumberManager::RandomNumberManager(size_t batch_size)
    : batch_size_(batch_size), rng_(std::random_device{}())
{
    for (const auto& dist : {"beta", "normal", "lognormal", "uniform", "binomial", "log10normal"}) {
        refill(dist);
    }
}

void RandomNumberManager::refill(const std::string& dist) {
    std::vector<double> buffer(batch_size_);
    // Segment edges: np.linspace(0, 1, batch_size+1)
    // Draw one uniform per segment, then shuffle — matches Python exactly
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    for (size_t i = 0; i < batch_size_; i++) {
        double lo = static_cast<double>(i) / batch_size_;
        double hi = static_cast<double>(i + 1) / batch_size_;
        buffer[i] = lo + uni(rng_) * (hi - lo);
    }
    std::shuffle(buffer.begin(), buffer.end(), rng_);
    buffers_[dist] = std::move(buffer);
    indices_[dist] = 0;
}

double RandomNumberManager::get(const std::string& dist) {
    if (indices_[dist] >= batch_size_)
        refill(dist);
    return buffers_[dist][indices_[dist]++];
}
```

**Step 3: Implement sampling functions**

```cpp
double random_beta_lhs(RandomNumberManager& rng, double a, double b,
                       double loc, double scale) {
    double u = rng.get("beta");
    return beta_ppf(u, a, b) * scale + loc;
}

double random_normal_lhs(RandomNumberManager& rng, double mu, double sigma) {
    double u = rng.get("normal");
    return normal_ppf(u) * sigma + mu;
}

double random_lognormal_lhs(RandomNumberManager& rng, double mean, double sigma) {
    double u = rng.get("lognormal");
    // scipy lognorm.ppf(u, s=sigma, scale=exp(mean)) = exp(normal_ppf(u)*sigma + mean)
    return exp(normal_ppf(u) * sigma + mean);
}

double random_uniform_lhs(RandomNumberManager& rng, double low, double high) {
    double u = rng.get("uniform");
    return low + (high - low) * u;
}

int random_binomial_lhs(RandomNumberManager& rng, int n, double p) {
    double u = rng.get("binomial");
    return binomial_ppf(u, n, p);
}

double random_log10normal_lhs(RandomNumberManager& rng, double mu, double sigma) {
    double u = rng.get("log10normal");
    return pow(10.0, normal_ppf(u) * sigma + mu);
}
```

**Step 4: Add test — verify LHS stratification**

```cpp
// Verify LHS: draw 10000 from 'uniform' buffer, check stratification
// Each of 10000 segments should have exactly 1 value
{
    RandomNumberManager mgr(10000);
    std::vector<int> bins(10, 0);
    for (int i = 0; i < 10000; i++) {
        double v = mgr.get("uniform");
        int b = std::min(int(v * 10), 9);
        bins[b]++;
    }
    // With LHS, each decile should have ~1000 values (±100 for tolerance)
    for (int b = 0; b < 10; b++) {
        if (bins[b] < 800 || bins[b] > 1200) {
            printf("FAIL: LHS bin %d has %d (expected ~1000)\n", b, bins[b]);
            return 1;
        }
    }
}
printf("PASS: LHS stratification\n");
```

**Step 5: Run test**

Expected: PASS

**Step 6: Commit**

```bash
git add C++/
git commit -m "feat: implement LHS random manager with stratified buffers"
```

---

## Phase 3: Model

### Task 7: Port occupancy parameters

**Objective:** Port all 25 occupancy categories with exact parameter values from `model.py`.

**Files:**
- Create: `C++/src/model.h`
- Create: `C++/src/model.cpp`

**Step 1: Define the OccupancyParams struct and the map**

```cpp
// model.h
#pragma once
#include "random_manager.h"
#include <string>
#include <map>

struct OccupancyParams {
    double PBR_GM;
    double PBR_GSD;
    double Cdrop_GM;
    double Cdrop_GSD;
    double d_GM;
    double d_GSD;
    double mask_eff;
    int    I0;
    double volume_m3;
    double community_rate;
    double ECAi;
};

extern const std::map<std::string, OccupancyParams> occupancy_params;

const OccupancyParams& get_occupancy_parameters(const std::string& category);

struct SimParameters {
    double D;
    double VOL;
    int    I0;
    double lambda_bio;
    double gamma;
    double mask_efficiency;
    double community_rate;
    double PBR;
};

SimParameters sample_parameters(RandomNumberManager& rng, const std::string& category);
double QER(RandomNumberManager& rng, const std::string& category);
std::pair<double, int> infection_probability(
    double ECAi, const SimParameters& par, RandomNumberManager& rng,
    const std::string& category, bool require_infectors = false,
    double override_community_rate = 0);
std::pair<double, int> compute_ECAi(
    const SimParameters& par, double target_P, RandomNumberManager& rng,
    const std::string& category, bool require_infectors = false,
    double override_community_rate = 0);
```

**Step 2: Implement occupancy_params with all 25 categories**

Copy every value from `model.py` lines 17-41 exactly. This is mechanical but critical — any typo changes results.

**Step 3: Write a test that verifies the count and a few values**

```cpp
// tests/test_model.cpp (initial)
#include "model.h"
#include <cstdio>
#include <cmath>

int main() {
    if (occupancy_params.size() != 25) { printf("FAIL: expected 25 categories, got %zu\n", occupancy_params.size()); return 1; }
    const auto& c = occupancy_params.at("Classroom");
    if (c.PBR_GM != 0.55 || c.I0 != 30 || c.volume_m3 != 320 || c.ECAi != 20) {
        printf("FAIL: Classroom params mismatch\n"); return 1;
    }
    const auto& g = occupancy_params.at("Gym");
    if (g.PBR_GM != 0.62 || g.I0 != 180 || g.volume_m3 != 1900 || g.ECAi != 40) {
        printf("FAIL: Gym params mismatch\n"); return 1;
    }
    const auto& d = occupancy_params.at("Dwelling");
    if (d.PBR_GSD != 4.3 || d.Cdrop_GM != 40e4 || d.I0 != 6 || d.ECAi != 15) {
        printf("FAIL: Dwelling params mismatch\n"); return 1;
    }
    printf("PASS: occupancy_params\n");
    return 0;
}
```

**Step 4: Run test**

Expected: PASS

**Step 5: Commit**

```bash
git add C++/
git commit -m "feat: port occupancy parameters to C++"
```

---

### Task 8: Implement sample_parameters, QER, infection_probability, compute_ECAi

**Objective:** Port the four core functions, matching `model.py` exactly.

**Files:**
- Modify: `C++/src/model.cpp`

**Step 1: Implement sample_parameters**

```cpp
SimParameters sample_parameters(RandomNumberManager& rng, const std::string& category) {
    const auto& occ = get_occupancy_parameters(category);
    SimParameters par;
    par.D = 1.0;
    par.VOL = occ.volume_m3;
    par.I0 = occ.I0;
    par.mask_efficiency = occ.mask_eff;
    par.community_rate = occ.community_rate;
    par.lambda_bio = random_lognormal_lhs(rng, log(0.52), log(1.9));
    par.gamma = random_uniform_lhs(rng, 0.42, 0.61);
    par.PBR = random_lognormal_lhs(rng, log(occ.PBR_GM), log(occ.PBR_GSD));
    return par;
}
```

**Step 2: Implement QER**

Port from `model.py:84-100`:
```cpp
double QER(RandomNumberManager& rng, const std::string& category) {
    const auto& occ = get_occupancy_parameters(category);
    double PBR = random_lognormal_lhs(rng, log(occ.PBR_GM), log(occ.PBR_GSD));
    double C_drop = random_lognormal_lhs(rng, log(occ.Cdrop_GM), log(occ.Cdrop_GSD));
    double d = random_lognormal_lhs(rng, log(occ.d_GM), log(occ.d_GSD));
    double E = random_beta_lhs(rng, 5.0, 2.0, 2.0, 3.0);
    double Vdrop = (M_PI / 6.0) * pow(d * E, 3) * C_drop;
    double GVL_ml = random_log10normal_lhs(rng, 7.0, 1.4);
    double GVL_m3 = GVL_ml * 1e6;
    double VF = random_beta_lhs(rng, 5.0, 2.0, 1e-4, (1e-2 - 1e-4));
    double RD = 1;
    double RTD = random_uniform_lhs(rng, 0.43, 0.65);
    double DK = random_uniform_lhs(rng, 5, 15);
    double VER = PBR * Vdrop * GVL_m3 * VF * RD;
    return RTD * VER / DK;
}
```

**Step 3: Implement infection_probability**

Port from `model.py:103-126`:
```cpp
std::pair<double, int> infection_probability(
    double ECAi, const SimParameters& par, RandomNumberManager& rng,
    const std::string& category, bool require_infectors,
    double override_community_rate)
{
    double TECAi = ECAi * par.I0 * 3.6;
    double phi = par.gamma + par.lambda_bio + TECAi / par.VOL;
    double comm_rate = (override_community_rate > 0) ? override_community_rate : par.community_rate;
    int n_infected = random_binomial_lhs(rng, par.I0, comm_rate);
    if (require_infectors) {
        while (n_infected == 0)
            n_infected = random_binomial_lhs(rng, par.I0, comm_rate);
    }
    int infected_flag = (n_infected > 0) ? 1 : 0;
    double P = 0;
    if (n_infected > 0) {
        double QER_sum = 0;
        for (int i = 0; i < n_infected; i++)
            QER_sum += QER(rng, category);
        double mask_factor = pow(1 - par.mask_efficiency, 2);
        double Q = (par.PBR * par.D * mask_factor / (phi * par.VOL)) * QER_sum;
        P = 1 - exp(-Q);
    }
    return {P, infected_flag};
}
```

**Step 4: Implement compute_ECAi**

Port from `model.py:129-161`:
```cpp
std::pair<double, int> compute_ECAi(
    const SimParameters& par, double target_P, RandomNumberManager& rng,
    const std::string& category, bool require_infectors,
    double override_community_rate)
{
    double comm_rate = (override_community_rate > 0) ? override_community_rate : par.community_rate;
    int n_infected = random_binomial_lhs(rng, par.I0, comm_rate);
    if (require_infectors) {
        while (n_infected == 0)
            n_infected = random_binomial_lhs(rng, par.I0, comm_rate);
    }
    int infected_flag = (n_infected > 0) ? 1 : 0;
    double ECAi_val = 0;
    if (n_infected > 0) {
        double Q = -log(1 - target_P);
        double QER_sum = 0;
        for (int i = 0; i < n_infected; i++)
            QER_sum += QER(rng, category);
        double mask_factor_sq = pow(1 - par.mask_efficiency, 2);
        double term1 = (par.PBR * par.D * mask_factor_sq / Q) * QER_sum;
        double term2 = par.VOL * (par.gamma + par.lambda_bio);
        ECAi_val = (1.0 / (3.6 * par.I0)) * (term1 - term2);
        if (ECAi_val < 0) ECAi_val = 0;
    }
    return {ECAi_val, infected_flag};
}
```

**Step 5: Add tests — smoke test for each function**

```cpp
// tests/test_model.cpp (append)
#include "random_manager.h"

int main() {
    // ... existing occupancy_params tests ...

    RandomNumberManager rng;

    // sample_parameters: check ranges
    auto par = sample_parameters(rng, "Classroom");
    if (par.D != 1.0 || par.VOL != 320 || par.I0 != 30) { printf("FAIL: sample_parameters\n"); return 1; }
    if (par.lambda_bio <= 0) { printf("FAIL: lambda_bio > 0\n"); return 1; }
    if (par.gamma < 0.42 || par.gamma > 0.61) { printf("FAIL: gamma range\n"); return 1; }

    // QER: should be positive
    double q = QER(rng, "Classroom");
    if (q <= 0 || !std::isfinite(q)) { printf("FAIL: QER = %f\n", q); return 1; }

    // infection_probability: should be in [0, 1)
    auto [prob, flag] = infection_probability(20, par, rng, "Classroom", false, 0);
    if (prob < 0 || prob >= 1 || !std::isfinite(prob)) { printf("FAIL: P = %f\n", prob); return 1; }

    // compute_ECAi: should be non-negative
    auto [ecai, flag2] = compute_ECAi(par, 0.001, rng, "Classroom", false, 0);
    if (ecai < 0 || !std::isfinite(ecai)) { printf("FAIL: ECAi = %f\n", ecai); return 1; }

    printf("PASS: model functions\n");
    return 0;
}
```

**Step 6: Run test**

Expected: PASS

**Step 7: Commit**

```bash
git add C++/
git commit -m "feat: implement sample_parameters, QER, infection_probability, compute_ECAi"
```

---

## Phase 4: Simulation Scripts

### Task 9: Port ecai.cpp

**Objective:** Port `ecai.py` — the main ECAi simulation script.

**Files:**
- Create: `C++/src/ecai.cpp`

**Step 1: Implement ecai.cpp**

Port from `ecai.py`:
- Read N from command line (default 10000)
- Loop over all 25 categories, N simulations each
- Compute 96th percentile ECAi
- Round to nearest 5 L/s/p and 10 CFM (LPS_TO_CFM = 2.11888)
- Track zero-infector count
- Print the same table format
- Write CSV with same column names: Category, Percentile_96_Lps, Rounded_Lps, Rounded_CFM, ZeroInfectedPercent

**Step 2: Verify output format matches Python**

Run Python: `python ecai.py 1000 > /tmp/py_ecai.txt`
Run C++: `./build/ecai 1000 > /tmp/cpp_ecai.txt`

Compare table structure (headers, column widths). Values won't match exactly (different RNG) but format must be identical.

**Step 3: Commit**

---

### Task 10: Port probability_ecai.cpp

**Objective:** Port `probabilityECAi.py` — probability at ASHRAE ECAi values with `--save-all` option.

**Files:**
- Create: `C++/src/probability_ecai.cpp`

**Step 1: Implement**

Port from `probabilityECAi.py`:
- Read N (positional, default 10000), `--save-all`, `--outfile`
- Loop categories, use `require_infectors=True`
- Compute 96th percentile probability
- Print same table format: Category, ECAi (L/s/p), 96th per P (%)
- Write CSV: Category, ECAi_Lps, P_96th_percentile
- For `--save-all`: write binary file (one double array per category, with a simple header listing category names and sizes). C++ analysis tools will read this. .npz compatibility is deferred.

**Step 2: Verify output format**

Run Python: `python probabilityECAi.py 1000`
Run C++: `./build/probability_ecai 1000`

Compare table headers and CSV columns.

**Step 3: Commit**

---

### Task 11: Port probability_scan.cpp

**Objective:** Port `probabilityScan.py`.

**Files:**
- Create: `C++/src/probability_scan.cpp`

**Step 1: Implement**

Port from `probabilityScan.py`:
- Scan ECAi from 5 to 100 L/s/p in steps of 5
- Find minimum ECAi where 96th percentile < 0.1%
- Print same table, write CSV: Category, ECAi_Lps_for_P_lt_0.1pct

**Step 2: Verify and commit**

---

### Task 12: Port single_probability.cpp

**Objective:** Port `singleProbability.py`.

**Files:**
- Create: `C++/src/single_probability.cpp`

**Step 1: Implement**

Port from `singleProbability.py`:
- Args: --N, --category, --community_rate, --ecai, --no_zero_infectors
- Run N simulations for one category
- Report 96th percentile and zero-infector percentage
- Skip plotting (no matplotlib in C++) — print a text histogram summary instead (optional: could pipe to gnuplot)

**Step 2: Verify and commit**

---

## Phase 5: Analysis Tools

### Task 13: Port percentiles.cpp

**Objective:** Port `analysis/percentiles.py` to read binary raw data files.

**Files:**
- Create: `C++/analysis/percentiles.cpp`

**Step 1: Implement**

Port the two reports (table + threshold) from `analysis/percentiles.py`. Read the binary format written by `probability_ecai --save-all`. Same coarse-to-fine search logic.

**Step 2: Verify and commit**

---

### Task 14: Port boxplot.cpp

**Objective:** Port `analysis/boxplot.py`. Since C++ has no matplotlib, produce either a text-based box plot or save data for external plotting.

**Files:**
- Create: `C++/analysis/boxplot.cpp`

**Step 1: Implement**

Options (pick one):
- (a) Text-based horizontal box plot in terminal (quartiles, whiskers, 96th percentile marker)
- (b) Write a CSV/TSV that can be piped to gnuplot
- (c) Use a header-only plotting library (plotutils, etc.)

Recommend (a) for simplicity — a text box plot that shows the median, quartiles, whiskers, and 96th percentile per category.

**Step 2: Verify and commit**

---

## Phase 6: Validation

### Task 15: Cross-validate C++ vs Python

**Objective:** Prove C++ produces statistically equivalent results to Python.

**Files:**
- Create: `C++/tests/test_validation.cpp`

**Step 1: Write validation test**

Strategy:
1. Run Python `ecai.py 100000` → save output CSV
2. Run C++ `ecai 100000` → save output CSV
3. Compare 96th percentile ECAi per category: relative difference should be < 5% for most categories (Monte Carlo error at N=100k)
4. Some high-variance categories (Dwelling, Exam) may have larger error — allow up to 10%

The test reads both CSVs and compares. It does NOT need to run the Python — it expects the Python CSV to already exist (or the test runner generates it as a prerequisite).

**Step 2: Run validation**

```bash
# Generate reference data
python ecai.py 100000
mv ecai_results.csv ecai_results_py.csv
./build/ecai 100000
# Run comparison test
./build/test_validation
```

Expected: PASS (all categories within tolerance)

**Step 3: Repeat for probabilityECAi**

```bash
python probabilityECAi.py 100000
mv ecai_ashrae241_96th_percentile.csv prob_py.csv
./build/probability_ecai 100000
./build/test_validation --mode prob
```

**Step 4: Commit**

```bash
git add C++/
git commit -m "test: cross-validate C++ vs Python at N=100000"
```

---

## Phase 7: Benchmarking

### Task 16: Benchmark C++ vs Python

**Objective:** Measure the speedup.

**Step 1: Time Python**

```bash
time python ecai.py 100000
```

**Step 2: Time C++**

```bash
time ./build/ecai 100000
```

**Step 3: Record results in README.md**

Update `C++/README.md` with:
- Build instructions
- Usage (matching Python CLI)
- Validation results (C++ vs Python 96th percentile comparison table)
- Benchmark results (speedup factor)

**Step 4: Commit**

```bash
git add C++/
git commit -m "docs: add C++ README with benchmarks and validation results"
```

---

## Risks and Tradeoffs

1. **Inverse beta CDF accuracy.** The bisection approach is robust but slower than Newton's method. If it proves too slow (unlikely — it's called once per QER, not per simulation iteration), switch to Newton's method with the beta PDF as derivative.

2. **Binomial with large n.** For Museum/Convention (I0=400) with community_rate=0.01, the binomial CDF loop runs up to ~401 iterations per draw. In C++ this is fast (~microseconds), but if it becomes a bottleneck, a normal approximation could be used for large n (but this would change results slightly — keep the exact version unless profiling shows it's the bottleneck).

3. **LHS buffer memory.** 6 distributions × 1M doubles × 8 bytes = 48 MB. Same as Python. Can reduce batch_size if needed.

4. **No plotting.** C++ scripts that produced histograms in Python (qer.py, distributions.py, singleProbability.py --show_plots) will not produce plots. The computational results are identical; visualization is deferred to Python or external tools.

5. **.npz compatibility.** Deferred. If needed later, add miniz (single-header ZIP library) to write .npz files readable by Python analysis tools. For now, C++ uses its own binary format.

6. **Platform variance.** `M_PI` may not be defined in MSVC — use `constexpr double PI = 3.14159265358979323846;` if needed. `std::lgamma` requires C++11 (available). `erfc` is in `<cmath>` (C++11).

---

## Open Questions

1. Should we also port the GUIs (mainGUI.py Tkinter, streamlitGUI.py Streamlit)? These would require a GUI framework in C++ (Qt, ImGui, web server). Recommend: skip for now, C++ is for batch computation speed.

2. Should C++ support the same `--save-all` .npz format for interchange with Python analysis tools? Recommend: add later if needed via miniz.

3. Build with MSVC or MinGW g++? Recommend: support both via CMake; test whichever is available on Felipe's machine first.