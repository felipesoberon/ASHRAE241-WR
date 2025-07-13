import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import random
import numpy as np

from model import occupancy_params, sample_parameters, compute_ECAi, infection_probability
from randomManager import RandomNumberManager

categories = list(occupancy_params.keys())

def generate_simulation_data(N=10000, mode="ECAi"):
    rng = RandomNumberManager()
    results = []

    if mode == "ECAi":
        ashrae_results = []
        for category in categories:
            ECAi_list = []
            for _ in range(N):
                par = sample_parameters(rng, category=category)
                comm_rate = 0 if use_ashrae_cir.get() else float(community_rate_var.get())
                ECAi_val, _ = compute_ECAi(par, target_P=0.001, rng=rng, category=category,
                                           require_infectors=not allow_zero_infector.get(),
                                           override_community_rate=comm_rate)
                ECAi_list.append(ECAi_val)
            percentile_96 = np.percentile(ECAi_list, 96)
            rounded_lps = int(np.ceil(percentile_96 / 5.0)) * 5
            results.append(rounded_lps)
            ashrae_results.append(occupancy_params[category].get("ECAi", 0))
            progress["value"] = (len(results)) / len(categories) * 100
            status_label.config(text=f"Simulating: {category} ({len(results)}/{len(categories)})")
            root.update_idletasks()
        xmax = max(max(results), max(ashrae_results), 50)
        return results, ashrae_results, "L/s/person", xmax

    else:  # Infection Probability
        for category in categories:
            probs = []
            for _ in range(N):
                par = sample_parameters(rng, category=category)
                comm_rate = 0 if use_ashrae_cir.get() else float(community_rate_var.get())
                ECAi_val = occupancy_params[category].get("ECAi", 0)
                prob, _ = infection_probability(ECAi_val, par, rng, category=category,
                                                require_infectors=not allow_zero_infector.get(),
                                                override_community_rate=comm_rate)
                probs.append(prob * 100)
            percentile_96 = np.percentile(probs, 96)
            results.append(percentile_96)
            progress["value"] = (len(results)) / len(categories) * 100
            status_label.config(text=f"Simulating: {category} ({len(results)}/{len(categories)})")
            root.update_idletasks()
        xmax = max(max(results), 100)
        return results, [], "Infection Probability (%)", xmax


def update_option_states(*args):
    use_ashrae_ecai.set(True)
    ashrae_ecai_check.state(["disabled", "selected"])
    community_rate_entry.config(state="disabled" if use_ashrae_cir.get() else "normal")


def run_simulation():
    progress["value"] = 0
    status_label.config(text="Running simulation...")

    try:
        N = int(simulation_count_var.get())
        if N < 1:
            raise ValueError
    except ValueError:
        N = 10000
        simulation_count_var.set("10000")
        status_label.config(text="Invalid N, using default = 10000")

    mode = calculation_type.get()
    values, ashrae_vals, xlabel, xmax = generate_simulation_data(N=N, mode=mode)

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
            ax.text(v1 + 1, i - bar_width / 2, f"{v1:.1f}", va='center', fontsize=8)
            ax.text(v2 + 1, i + bar_width / 2, f"{v2:.1f}", va='center', fontsize=8)
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
    progress["value"] = 100
    status_label.config(text="Calculation complete.")



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

ttk.Label(left_frame, text="Community Infection Rate (0–1):").pack(anchor="w")
community_rate_var = tk.StringVar(value="0.01")
community_rate_entry = ttk.Entry(left_frame, textvariable=community_rate_var, width=10)
community_rate_entry.pack(anchor="w", pady=(0, 10))

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
