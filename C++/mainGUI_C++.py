import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import subprocess
import csv
import os
import threading
import shlex
import re
import numpy as np

# ---- C++ engine configuration ----

# Directory containing the C++ executables (build output).
# The GUI runs from the C++ folder, so build/ is a sibling.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BUILD_DIR = os.path.join(SCRIPT_DIR, "build")

# Project root (parent of C++/) — C++ executables write CSVs to their cwd.
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

# ASHRAE 241 ECAi values per category — needed for the ECAi-mode comparison bars.
# Same values as model.py occupancy_params[...]["ECAi"].
ASHRAE_ECAI = {
    "Cell": 15, "Dayroom": 20, "Food": 30, "Gym": 40, "Office": 15,
    "Retail": 20, "Transportation": 30, "Classroom": 20, "Lecture": 25,
    "Manufacturing": 25, "Sorting": 10, "Warehouse": 10, "Exam": 20,
    "Group": 35, "Patient": 35, "Resident": 25, "Waiting": 45,
    "Auditorium": 25, "Place": 25, "Museum": 30, "Convention": 30,
    "Spectator": 25, "Lobbies": 25, "Common": 25, "Dwelling": 15,
}

# Category display order (matches Python model.py insertion order).
categories = [
    "Cell", "Dayroom", "Food", "Gym", "Office", "Retail", "Transportation",
    "Classroom", "Lecture", "Manufacturing", "Sorting", "Warehouse",
    "Exam", "Group", "Patient", "Resident", "Waiting",
    "Auditorium", "Place", "Museum", "Convention", "Spectator",
    "Lobbies", "Common", "Dwelling",
]


# ---- WSL plumbing (auto-detected, not hardcoded) ----
# These helpers only locate the distro and translate paths so the C++ engine
# can be launched. They do not touch simulation arguments, flags, or results.

_WSL_DISTRO = None
_WSL_DISTRO_RESOLVED = False


def default_wsl_distro():
    """Return the default WSL distro (the one marked '*' by `wsl -l -v`),
    skipping utility distros like docker-desktop. Returns None if detection
    fails, in which case callers omit -d and let WSL use its own default.
    Cached after the first lookup. wsl.exe emits UTF-16 output."""
    global _WSL_DISTRO, _WSL_DISTRO_RESOLVED
    if _WSL_DISTRO_RESOLVED:
        return _WSL_DISTRO
    _WSL_DISTRO_RESOLVED = True
    distro = None
    try:
        out = subprocess.run(["wsl", "-l", "-v"], capture_output=True).stdout
        text = out.decode("utf-16-le", "ignore")
        for raw in text.splitlines():
            line = raw.strip()
            if line.startswith("*"):                     # default distro row
                m = re.search(r"[A-Za-z0-9._-]+", line[1:])
                if m and "docker" not in m.group(0).lower():
                    distro = m.group(0)
                break
    except Exception:
        distro = None
    _WSL_DISTRO = distro
    return _WSL_DISTRO


def _wsl_base():
    """`wsl` command prefix, pinned to the detected distro when available."""
    distro = default_wsl_distro()
    return ["wsl", "-d", distro] if distro else ["wsl"]


def to_wsl_path(win_path):
    """Convert a Windows path to its WSL mount path via `wslpath` (handles any
    drive letter and spaces); fall back to a C: string replace if wslpath fails."""
    try:
        out = subprocess.run(_wsl_base() + ["--", "wslpath", "-a", win_path],
                             capture_output=True).stdout
        p = out.decode("utf-8", "ignore").strip()
        if p:
            return p
    except Exception:
        pass
    return win_path.replace("C:\\", "/mnt/c/").replace("\\", "/")


def run_cpp_simulation(calculation_type, N, allow_zero_infector,
                       use_ashrae_cir, general_rate, healthcare_rate,
                       progress_callback=None, use_fitted_qer=False):
    """Run the C++ simulation executable and return (results, ashrae_results).

    The C++ executables are Linux ELF binaries built in WSL, so they must be
    invoked through WSL. The CSV output is written to the project root.

    If progress_callback is provided, it is called as:
        progress_callback(category_name, completed_count, total_count)

    Returns:
        results: list of float — 96th percentile values per category
        ashrae_results: list of float — ASHRAE ECAi values per category
                        (empty for Infection Probability mode)
    """
    if calculation_type == "ECAi":
        exe_rel = "C++/build/ecai"
        args = [str(N)]
        if not allow_zero_infector:
            args.append("--require-infectors")
        if not use_ashrae_cir:
            args += ["--community-rate-general", str(general_rate),
                     "--community-rate-healthcare", str(healthcare_rate)]
        if use_fitted_qer:
            args.append("--use-fitted-qer")
        csv_file = os.path.join(PROJECT_ROOT, "ecai_results.csv")
    else:
        exe_rel = "C++/build/probability_ecai"
        args = [str(N)]
        if allow_zero_infector:
            args.append("--no-require-infectors")
        if not use_ashrae_cir:
            args += ["--community-rate-general", str(general_rate),
                     "--community-rate-healthcare", str(healthcare_rate)]
        if use_fitted_qer:
            args.append("--use-fitted-qer")
        csv_file = os.path.join(PROJECT_ROOT, "ecai_ashrae241_96th_percentile.csv")

    # Convert the Windows project path to its WSL path (handles any drive
    # letter and spaces via wslpath) and locate the executable inside it.
    wsl_root = to_wsl_path(PROJECT_ROOT)
    wsl_exe = wsl_root + "/" + exe_rel

    # Fail early with a helpful message if the C++ engine has not been built.
    check = subprocess.run(_wsl_base() + ["--", "test", "-x", wsl_exe])
    if check.returncode != 0:
        raise RuntimeError(
            "C++ executable not found or not built: " + exe_rel + "\n"
            "Build it in WSL first:\n"
            "    cd " + wsl_root + "/C++ && cmake -B build && cmake --build build")

    # Run the C++ executable through WSL with line-buffered output
    # so we can read each category's result line as it's produced.
    # Quote the paths so directories containing spaces still work; the
    # simulation arguments (args) are passed through unchanged.
    shell_cmd = "stdbuf -oL " + shlex.quote(wsl_exe) + " " + " ".join(args)
    cmd = _wsl_base() + ["--", "bash", "-c",
                         "cd " + shlex.quote(wsl_root) + " && " + shell_cmd]

    total = len(categories)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)
    completed = 0
    for line in proc.stdout:
        line = line.strip()
        # Data rows look like: "| Cell            |         11.60 | ..."
        if line.startswith("|") and not line.startswith("|-") and not line.startswith("| Category"):
            completed += 1
            # Extract category name (first field after "|")
            parts = line.split("|")
            if len(parts) >= 3:
                cat_name = parts[1].strip()
                if progress_callback:
                    progress_callback(cat_name, completed, total)
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"C++ executable failed with code {proc.returncode}")

    # Read the CSV output
    results = []
    ashrae_results = []
    with open(csv_file, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            cat = row["Category"]
            if calculation_type == "ECAi":
                p96 = float(row["Percentile_96_Lps"])
                rounded = int(np.ceil(p96 / 5.0)) * 5
                results.append(rounded)
                ashrae_results.append(ASHRAE_ECAI.get(cat, 0))
            else:
                p96 = float(row["P_96th_percentile"])
                results.append(p96)
                ashrae_results.append(0)

    return results, ashrae_results


def generate_simulation_data(N=10000, mode="ECAi", progress_callback=None):
    # Build flags from GUI state. When not using the ASHRAE 241 defaults,
    # the general (default 1%) and healthcare (default 3%) community infection
    # rates are read separately from their entry fields.
    allow_zero = allow_zero_infector.get()
    use_ashrae = use_ashrae_cir.get()
    if use_ashrae:
        general_rate, healthcare_rate = 0.01, 0.03
    else:
        general_rate = float(community_rate_general_var.get())
        healthcare_rate = float(community_rate_healthcare_var.get())

    # Run C++ simulation
    values, ashrae_vals = run_cpp_simulation(
        mode, N, allow_zero, use_ashrae, general_rate, healthcare_rate,
        progress_callback, use_fitted_qer.get())

    if mode == "ECAi":
        xmax = max(max(values), max(ashrae_vals), 50)
        xlabel = "L/s/person"
    else:
        xlabel = "Infection Probability (%)"

    return values, ashrae_vals, xlabel, xmax if mode == "ECAi" else 100


def update_option_states(*args):
    use_ashrae_ecai.set(True)
    ashrae_ecai_check.state(["disabled", "selected"])
    entry_state = "disabled" if use_ashrae_cir.get() else "normal"
    community_rate_general_entry.config(state=entry_state)
    community_rate_healthcare_entry.config(state=entry_state)


def run_simulation():
    progress.config(mode="determinate")
    progress["value"] = 0
    status_label.config(text="Running simulation...")
    run_button.config(state="disabled")

    try:
        N = int(simulation_count_var.get())
        if N < 1:
            raise ValueError
    except ValueError:
        N = 10000
        simulation_count_var.set("10000")
        status_label.config(text="Invalid N, using default = 10000")

    mode = calculation_type.get()

    def on_progress(cat_name, completed, total):
        # Called from the worker thread — marshal to GUI thread
        root.after(0, lambda: _update_progress(cat_name, completed, total))

    def worker():
        try:
            values, ashrae_vals, xlabel, xmax = generate_simulation_data(
                N=N, mode=mode, progress_callback=on_progress)
            root.after(0, lambda: finish_simulation(values, ashrae_vals, xlabel, xmax, mode))
        except Exception as e:
            root.after(0, lambda: simulation_error(str(e)))

    thread = threading.Thread(target=worker, daemon=True)
    thread.start()


def _update_progress(cat_name, completed, total):
    pct = completed / total * 100
    progress["value"] = pct
    status_label.config(text=f"Simulating: {cat_name} ({completed}/{total})")


def finish_simulation(values, ashrae_vals, xlabel, xmax, mode):
    progress.stop()
    progress.config(mode="determinate")
    progress["value"] = 100

    ax.clear()
    y_pos = np.arange(len(categories))
    bar_width = 0.4

    # Avoid zero in log scale: set minimum displayable value
    if mode == "Infection Probability":
        values = [max(v, 0.001) for v in values]  # Avoid log(0)
        ax.set_xscale('log')
        ax.set_xlim(0.01, 100)
        ax.set_xticks([0.01, 0.1, 1, 10, 100])
        ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())

    ax.barh(y_pos - bar_width / 2, values, height=bar_width, label="Simulated", color="steelblue")

    if mode == "ECAi":
        ax.barh(y_pos + bar_width / 2, ashrae_vals, height=bar_width, label="ASHRAE 241", color="gray")
        for i, (v1, v2) in enumerate(zip(values, ashrae_vals)):
            ax.text(v1 + 1, i - bar_width / 2, f"{int(v1)}", va='center', fontsize=8)
            ax.text(v2 + 1, i + bar_width / 2, f"{int(v2)}", va='center', fontsize=8)
        ax.set_xlim(0, xmax)
    else:
        for i, v in enumerate(values):
            ax.text(v + 0.01, i - bar_width / 2, f"{v:.2f}", va='center', fontsize=8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(categories)
    ax.set_xlabel(xlabel)
    ax.invert_yaxis()
    ax.grid(axis='x', linestyle='--', color='gray', linewidth=0.7)
    ax.legend(loc='lower right')
    fig.tight_layout()
    canvas.draw()
    status_label.config(text="Calculation complete.")
    run_button.config(state="normal")


def simulation_error(msg):
    progress.stop()
    progress.config(mode="determinate")
    progress["value"] = 0
    status_label.config(text=f"Error: {msg}")
    run_button.config(state="normal")


# ---- GUI setup (identical layout to mainGUI.py) ----

root = tk.Tk()
root.title("Simulation Selector")
root.geometry("1000x740")

left_frame = ttk.Frame(root)
left_frame.pack(side=tk.LEFT, padx=20, pady=20, anchor="n", fill="y")

ttk.Label(left_frame, text="Calculation Type:").pack(anchor="w", pady=(0, 5))
calculation_type = ttk.Combobox(left_frame, values=["ECAi", "Infection Probability"], state="readonly")
calculation_type.set("ECAi")
calculation_type.pack(anchor="w", fill="x", pady=(0, 10))
calculation_type.bind("<<ComboboxSelected>>", update_option_states)

use_ashrae_ecai = tk.BooleanVar(value=True)
ashrae_ecai_check = ttk.Checkbutton(left_frame, text="Use ASHRAE 241 ECAi values", variable=use_ashrae_ecai)
ashrae_ecai_check.state(["disabled", "selected"])
ashrae_ecai_check.pack(anchor="w", pady=(0, 10))

allow_zero_infector = tk.BooleanVar(value=True)
zero_infector_check = ttk.Checkbutton(left_frame, text="Allow Zero Infector Simulations", variable=allow_zero_infector)
zero_infector_check.pack(anchor="w", pady=(0, 10))

use_ashrae_cir = tk.BooleanVar(value=True)
ashrae_cir_check = ttk.Checkbutton(left_frame, text="Use ASHRAE 241 Community Infection Rate", variable=use_ashrae_cir, command=update_option_states)
ashrae_cir_check.pack(anchor="w", pady=(0, 5))

ttk.Label(left_frame, text="General spaces rate (0-1):").pack(anchor="w")
community_rate_general_var = tk.StringVar(value="0.01")
community_rate_general_entry = ttk.Entry(left_frame, textvariable=community_rate_general_var, width=10)
community_rate_general_entry.pack(anchor="w", pady=(0, 5))

ttk.Label(left_frame, text="Healthcare spaces rate (0-1):").pack(anchor="w")
community_rate_healthcare_var = tk.StringVar(value="0.03")
community_rate_healthcare_entry = ttk.Entry(left_frame, textvariable=community_rate_healthcare_var, width=10)
community_rate_healthcare_entry.pack(anchor="w", pady=(0, 10))

use_fitted_qer = tk.BooleanVar(value=False)
fitted_qer_check = ttk.Checkbutton(left_frame, text="Use fitted log-normal QER (faster)", variable=use_fitted_qer)
fitted_qer_check.pack(anchor="w", pady=(0, 10))

ttk.Label(left_frame, text="Number of Simulations (N):").pack(anchor="w")
simulation_count_var = tk.StringVar(value="10000")
simulation_count_entry = ttk.Entry(left_frame, textvariable=simulation_count_var, width=10)
simulation_count_entry.pack(anchor="w", pady=(0, 10))

run_button = ttk.Button(left_frame, text="Run", command=run_simulation)
run_button.pack(anchor="w", pady=(0, 10))

progress = ttk.Progressbar(left_frame, orient="horizontal", length=200, mode="determinate")
progress.pack(anchor="w", pady=(0, 10))

status_label = ttk.Label(left_frame, text="")
status_label.pack(anchor="w")

right_frame = ttk.Frame(root)
right_frame.pack(side=tk.RIGHT, padx=10, pady=10, fill=tk.BOTH, expand=True)

fig, ax = plt.subplots(figsize=(6, 10))
initial_values = [0 for _ in categories]
y_pos = np.arange(len(categories))
bar_width = 0.4
ax.barh(y_pos - bar_width / 2, initial_values, height=bar_width, label="Simulated", color="steelblue")
ax.barh(y_pos + bar_width / 2, initial_values, height=bar_width, label="ASHRAE 241", color="gray")
ax.set_yticks(y_pos)
ax.set_yticklabels(categories)
ax.set_xlabel('L/s/person')
ax.set_xlim(0, 50)
ax.invert_yaxis()
ax.grid(axis='x', linestyle='--', color='gray', linewidth=0.7)
ax.legend(loc='lower right')
fig.tight_layout()

canvas = FigureCanvasTkAgg(fig, master=right_frame)
canvas.draw()
canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

update_option_states()
root.mainloop()