import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from model import occupancy_params, sample_parameters, compute_ECAi, infection_probability
from randomManager import RandomNumberManager

st.set_page_config(layout="wide")
st.title("Simulation Selector")

categories = list(occupancy_params.keys())

# Sidebar for user inputs
st.sidebar.header("Simulation Settings")
calculation_type = st.sidebar.selectbox("Calculation Type:", ["ECAi", "Infection Probability"])

#use_ashrae_ecai = st.sidebar.checkbox("Use ASHRAE 241 ECAi values", value=True)
st.sidebar.markdown("#### Model uses ASHRAE 241 ECAi values (Probability Calculation)")

use_ashrae_ecai = True

allow_zero_infector = st.sidebar.checkbox("Allow Zero Infector Simulations", value=True)
use_ashrae_cir = st.sidebar.checkbox("Use ASHRAE 241 Community Infection Rate", value=True)

community_rate = st.sidebar.text_input("Community Infection Rate (0–1):", value="0.01")
simulation_count = st.sidebar.number_input("Number of Simulations (N):", min_value=1, value=10000, step=1000)

run = st.sidebar.button("Run Simulation")

if run:
    st.write(f"Running {simulation_count} simulations for each occupancy category...")
    progress_bar = st.progress(0)
    status_text = st.empty()

    rng = RandomNumberManager()
    N = simulation_count
    results = []
    ashrae_results = []

    for idx, category in enumerate(categories):
        data = []
        for _ in range(N):
            par = sample_parameters(rng, category=category)
            comm_rate = 0 if use_ashrae_cir else float(community_rate)
            if calculation_type == "ECAi":
                val, _ = compute_ECAi(par, 0.001, rng, category, not allow_zero_infector, comm_rate)
            else:
                val, _ = infection_probability(
                    occupancy_params[category].get("ECAi", 0),
                    par,
                    rng,
                    category,
                    not allow_zero_infector,
                    comm_rate
                )
                val *= 100
            data.append(val)
        percentile_96 = np.percentile(data, 96)
        rounded_val = int(np.ceil(percentile_96 / 5.0)) * 5 if calculation_type == "ECAi" else percentile_96
        results.append(rounded_val)
        ashrae_results.append(occupancy_params[category].get("ECAi", 0))
        progress_bar.progress((idx + 1) / len(categories))
        status_text.write(f"Simulating {category} ({idx + 1}/{len(categories)})")

    fig, ax = plt.subplots(figsize=(8, 10))
    y_pos = np.arange(len(categories))
    bar_width = 0.4

    if calculation_type == "Infection Probability":
        results = [max(v, 0.001) for v in results]
        ax.set_xscale('log')
        ax.set_xlim(0.01, 100)
        ax.set_xticks([0.01, 0.1, 1, 10, 100])
        ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        xlabel = "Infection Probability (%)"
    else:
        xlabel = "L/s/person"

    ax.barh(y_pos - bar_width / 2, results, height=bar_width, label="Simulated", color="steelblue")
    if calculation_type == "ECAi":
        ax.barh(y_pos + bar_width / 2, ashrae_results, height=bar_width, label="ASHRAE 241", color="gray")

    for i, v in enumerate(results):
        offset = 1 if calculation_type == "ECAi" else 0.01
        ax.text(v + offset, i - bar_width / 2, f"{v:.1f}", va='center', fontsize=8)
    if calculation_type == "ECAi":
        for i, v in enumerate(ashrae_results):
            ax.text(v + 1, i + bar_width / 2, f"{v:.1f}", va='center', fontsize=8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(categories)
    ax.set_xlabel(xlabel)
    ax.invert_yaxis()
    ax.grid(axis='x', linestyle='--', color='gray')
    ax.legend()

    st.pyplot(fig)
    st.success("Calculation complete.")
