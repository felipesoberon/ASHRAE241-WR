# C++-Backed GUI Implementation Plan

> **For Hermes:** Use subagent-driven-development skill to implement this plan task-by-task.

**Goal:** Create `C++/mainGUI_C++.py` — a Python Tkinter GUI that is visually and functionally identical to `mainGUI.py`, but delegates all Monte Carlo computation to the C++ engine instead of the Python model.

**Architecture:** The GUI is a thin Python/Tkinter layer that constructs and runs C++ simulation executables as subprocesses, parses their CSV output, and displays the results using matplotlib — exactly as the original GUI does. The C++ engine (`ecai` or `probability_ecai`) runs as a background process; the GUI reads the resulting CSV and plots it. No C++ code changes are needed. The Python model is not imported or called.

**Tech Stack:** Python 3, Tkinter, matplotlib, subprocess, CSV parsing. C++ executables from `C++/build/`.

---

## Constraints

1. **Do not modify `mainGUI.py`.** It stays as the reference.
2. **Do not modify any C++ source files.** The C++ engine is complete and validated.
3. **The GUI must look and behave the same as `mainGUI.py`:**
   - Same window title, size, layout (left control panel, right chart).
   - Same controls: Calculation Type (ECAi / Infection Probability), Allow Zero Infector checkbox, Use ASHRAE 241 Community Infection Rate checkbox, Community Infection Rate entry, Number of Simulations entry, Run button, progress bar, status label.
   - Same chart: horizontal bar plot, simulated vs ASHRAE values, log scale for probability mode, same colors, labels, legend, grid.
4. **The C++ engine does NOT support all the Python GUI's runtime options directly.** Specifically:
   - The C++ `ecai` executable does not accept `--require_infectors` or `--override_community_rate` flags — it always uses `require_infectors=False` and the category default community rate.
   - The C++ `probability_ecai` executable always uses `require_infectors=True` and the category default community rate.
   - The "Allow Zero Infector" and "Community Infection Rate" GUI controls affect the Python model at runtime, but the C++ executables don't expose these as CLI flags.
   - **Resolution:** We add minimal CLI flags to the C++ executables (`--require-infectors`, `--no-require-infectors`, `--community-rate`) so the GUI can pass them through. This is a small, backward-compatible addition to the C++ code that does not change default behavior.
5. **File location:** `C++/mainGUI_C++.py` — inside the C++ folder, as requested.
6. **Build artifacts:** The GUI assumes the C++ executables are already built (in `C++/build/`). It does not build them.

---

## Project Structure

```
C++/
├── mainGUI_C++.py          # NEW: Python Tkinter GUI using C++ engine
├── src/
│   ├── ecai.cpp            # MODIFY: add --require-infectors, --community-rate flags
│   ├── probability_ecai.cpp # MODIFY: add --no-require-infectors, --community-rate flags
│   └── ... (all other files unchanged)
├── build/
│   ├── ecai                # rebuilt after adding flags
│   └── probability_ecai   # rebuilt after adding flags
└── ...
```

---

## Phase 1: Add CLI Flags to C++ Executables

### Task 1: Add runtime flags to `ecai.cpp`

**Objective:** Add `--require-infectors` and `--community-rate RATE` flags to the C++ `ecai` executable, matching the Python GUI's runtime options. Default behavior unchanged.

**Files:**
- Modify: `C++/src/ecai.cpp`

**Step 1: Add argument parsing**

After the existing N parsing, add:

```cpp
bool require_infectors = false;  // default: match current behavior
double community_rate = -1;       // -1 = use category default

for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];
    if (arg == "--require-infectors") {
        require_infectors = true;
    } else if (arg == "--community-rate" && i + 1 < argc) {
        community_rate = std::atof(argv[++i]);
    } else if (arg[0] != '-' && i == 1) {
        // First positional = N (existing behavior)
        N = std::atoi(arg.c_str());
    }
}
```

**Step 2: Pass flags to compute_ECAi**

Change the compute_ECAi call:
```cpp
auto [ecai_val, infected_flag] = compute_ECAi(par, target_P, rng, category,
    require_infectors, community_rate);
```

**Step 3: Rebuild and verify default behavior unchanged**

Run: `./build/ecai 1000` (no flags) — output must match previous behavior.
Run: `./build/ecai 1000 --require-infectors --community-rate 0.05` — must use those options.

**Step 4: Commit**

---

### Task 2: Add runtime flags to `probability_ecai.cpp`

**Objective:** Add `--no-require-infectors` and `--community-rate RATE` flags. Default behavior unchanged (require_infectors=True).

**Files:**
- Modify: `C++/src/probability_ecai.cpp`

**Step 1: Add argument parsing**

```cpp
bool require_infectors = true;   // default: match current behavior
double community_rate = -1;      // -1 = use category default

// In the existing arg loop, add:
if (arg == "--no-require-infectors") {
    require_infectors = false;
} else if (arg == "--community-rate" && i + 1 < argc) {
    community_rate = std::atof(argv[++i]);
}
```

**Step 2: Pass flags to infection_probability**

```cpp
auto [prob, _] = infection_probability(ECAi, par, rng, category,
    require_infectors, community_rate);
```

**Step 3: Rebuild and verify**

**Step 4: Commit**

---

## Phase 2: GUI Implementation

### Task 3: Write `mainGUI_C++.py`

**Objective:** Create the GUI that mirrors `mainGUI.py` exactly, but calls C++ executables instead of the Python model.

**Files:**
- Create: `C++/mainGUI_C++.py`

**Step 1: Copy the structure from `mainGUI.py`**

The file will have the same layout:
- Same imports (tkinter, matplotlib, numpy)
- Remove `from model import ...` and `from randomManager import ...`
- Add `import subprocess, csv, os, tempfile`
- Same window setup, same controls, same layout
- Replace `generate_simulation_data()` with a function that:
  1. Builds the C++ command line with the appropriate flags
  2. Runs the C++ executable as a subprocess
  3. Reads the output CSV
  4. Returns the results in the same format as the original function

**Step 2: Implement the C++ subprocess runner**

```python
def run_cpp_simulation(calculation_type, N, allow_zero_infector, use_ashrae_cir, community_rate):
    """Run the C++ simulation and return (results, ashrae_results) lists."""
    # Determine executable and flags
    if calculation_type == "ECAi":
        exe = os.path.join(BUILD_DIR, "ecai")
        args = [str(N)]
        if not allow_zero_infector:
            args.append("--require-infectors")
        if not use_ashrae_cir:
            args.append("--community-rate")
            args.append(str(community_rate))
        csv_file = "ecai_results.csv"
    else:
        exe = os.path.join(BUILD_DIR, "probability_ecai")
        args = [str(N)]
        if allow_zero_infector:
            args.append("--no-require-infectors")
        if not use_ashrae_cir:
            args.append("--community-rate")
            args.append(str(community_rate))
        csv_file = "ecai_ashrae241_96th_percentile.csv"

    # Run the C++ executable
    subprocess.run([exe] + args, capture_output=True, text=True, check=True)

    # Read the CSV output
    results = []
    ashrae_results = []
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            ...
    return results, ashrae_results
```

**Step 3: Map the results to the GUI's display format**

For ECAi mode: read `Percentile_96_Lps`, round to nearest 5, read ASHRAE ECAi from the params (hardcoded or read from a small config).
For Probability mode: read `P_96th_percentile`.

**Step 4: Keep the chart code identical to mainGUI.py**

Same matplotlib FigureCanvasTkAgg, same barh plot, same colors, same log scale for probability, same legend, same axis labels.

**Step 5: Handle the progress bar and status label**

The C++ subprocess runs synchronously (blocks until done). During the run:
- Set progress to 0, status to "Running simulation..."
- Run the subprocess
- Set progress to 100, status to "Calculation complete."
- Update the chart

Note: The original GUI updates progress per-category during simulation. The C++ version runs all categories in one subprocess call, so we can only show 0% and 100%. This is an acceptable difference — the subprocess is fast (seconds, not minutes).

**Step 6: Test the GUI**

Run: `cd C++ && python mainGUI_C++.py`
Verify: window appears, controls work, Run button produces a chart.

**Step 7: Commit**

---

### Task 4: Rebuild C++ and end-to-end test

**Objective:** Rebuild the C++ executables with the new flags and test the GUI end-to-end.

**Step 1: Rebuild C++**

Run: `cd C++ && cmake --build build`

**Step 2: Run the GUI**

Run: `cd C++ && python mainGUI_C++.py`

**Step 3: Test each mode**

- ECAi mode, default settings, N=1000 → chart appears with simulated vs ASHRAE bars
- Infection Probability mode, default settings, N=1000 → log-scale chart
- Toggle "Allow Zero Infector" off → results change
- Toggle "Use ASHRAE 241 Community Infection Rate" off, enter custom rate → results change
- Increase N to 10000 → chart updates (slower but completes)

**Step 4: Commit**

---

## Risks and Tradeoffs

1. **Progress bar granularity.** The original GUI updates progress per-category. The C++ version runs all categories in one subprocess, so progress goes 0% → 100%. This is acceptable because the C++ run is ~20x faster.

2. **Subprocess blocking.** The C++ subprocess blocks the GUI thread during computation. For N=10000 this takes ~8 seconds (vs ~169 seconds in Python). If this is a problem, the subprocess could be run in a background thread with periodic polling, but this adds complexity. Start with the simple synchronous approach.

3. **CSV file location.** The C++ executables write CSVs to the current working directory. The GUI must `os.chdir()` to a known location or use absolute paths. Simplest: run the subprocess with `cwd` set to the project root.

4. **ASHRAE ECAi values for the chart.** The Python GUI reads these from `occupancy_params`. The C++ GUI can either: (a) hardcode them, (b) read them from the C++ CSV output (which includes them), or (c) import them from the Python model (but we want to avoid importing the Python model). Best: read from the C++ CSV output, which includes `ECAi_Lps` for probability mode, or hardcode a small dict for ECAi mode.

5. **Category order.** The C++ executables now output in Python insertion order (after the cosmetic fixes), so the chart will display categories in the same order as the Python GUI.

---

## Open Questions

1. Should the GUI offer to build the C++ executables if they're missing? Recommend: no — show an error message instead, pointing to the build instructions.

2. Should the subprocess run in a background thread to avoid blocking? Recommend: start synchronous, add threading only if the user finds the blocking noticeable.