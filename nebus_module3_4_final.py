"""
NEBUS Module 3 & Module 4 — Realistic TEG Feasibility & Habitat System Integration.
Run: python nebus_module3_4_final.py
"""
from pathlib import Path
import json
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent / "nebus_module3_4_results"
OUT.mkdir(exist_ok=True)

# ---------------------------------------------------------
# VERIFIED INPUTS FROM MODULE 1, MODULE 2 & CAD STEP MODEL
# ---------------------------------------------------------
SITE_COORDS = {"lat": -4.57, "lon_east": 137.95}
TES_TI_OBSERVED = 405.0          # NASA TES thermal inertia (J m^-2 K^-1 s^-1/2)
REMS_GROUND_MEAN_K = 211.45      # NASA REMS preliminary mean ground temperature (K)
MARS_SOL_HOURS = 24.6597
MARS_GRAVITY_M_S2 = 3.71

# CAD-Derived Mass & Structural Baseline (51 solid bodies)
MASS_BREAKDOWN_KG = {
    "Al_Li_primary_structure": 3736.8,
    "Ti6Al4V_six_landing_legs": 2002.4,
    "B_UHMWPE_internal_platform": 911.8,
    "TEG_ring_and_footpad_layers": 588.4,
    "Fused_silica_outer_dome": 248.3,
    "FR4_internal_electronics": 133.8,
    "Silica_aerogel_inner_dome": 100.1,
    "B_UHMWPE_leg_inserts": 81.4,
    "Melanin_GO_UHMWPE_composite": 73.3,
    "PCM_thermal_buffer_ring": 6.8
}
TOTAL_DRY_MASS_KG = round(sum(MASS_BREAKDOWN_KG.values()), 1)
TOTAL_WEIGHT_KN = round((TOTAL_DRY_MASS_KG * MARS_GRAVITY_M_S2) / 1000.0, 2)
LOAD_PER_LEG_KN = round(TOTAL_WEIGHT_KN / 6.0, 2)

# ---------------------------------------------------------
# MODULE 3: REALISTIC HYBRID TEG FEASIBILITY (REMS-ANCHORED)
# ---------------------------------------------------------
# Diurnal ground temperature cycle centered at REMS 211.45 K, minimum at ~05:30 local solar time
hours = np.linspace(0, MARS_SOL_HOURS, 1000)
# Peak at ~13:30 (t = 13.5 h), trough at ~01:10-05:30 realistic Mars thermal lag
ground_temp_K = REMS_GROUND_MEAN_K - 20.0 * np.cos(2 * np.pi * (hours - 5.5) / MARS_SOL_HOURS)

# Hot reservoir supported by 133.8 kg FR4 electronics waste heat + 6.8 kg PCM buffer
HOT_SIDE_K = 250.0
CONDUCTANCE_W_K = 0.50
CARNOT_FRACTION = 0.15

delta_T = np.maximum(HOT_SIDE_K - ground_temp_K, 0.0)
heat_flow_W = CONDUCTANCE_W_K * delta_T
carnot_eff = np.where(HOT_SIDE_K > 0, delta_T / HOT_SIDE_K, 0.0)
teg_eff = CARNOT_FRACTION * carnot_eff
footpad_teg_W = heat_flow_W * teg_eff

footpad_peak_W = float(np.max(footpad_teg_W))
footpad_avg_W = float(np.mean(footpad_teg_W))
footpad_energy_Wh_sol = float(np.trapezoid(footpad_teg_W, hours))

# Dome UV Photothermal TEG (from Module 2 baseline)
R_OUTER_M = 1.118034
DOME_AREA_M2 = math.pi * R_OUTER_M**2
UV_W_M2 = 15.0
UV_ABS = 0.85
CAPTURE_FRAC = 0.60
DOME_TEG_EFF = 0.05
DOME_ACTIVE_HOURS = 5.0

dome_teg_W = UV_W_M2 * DOME_AREA_M2 * UV_ABS * CAPTURE_FRAC * DOME_TEG_EFF
dome_energy_Wh_sol = dome_teg_W * DOME_ACTIVE_HOURS

TOTAL_HYBRID_ENERGY_WH_SOL = round(footpad_energy_Wh_sol + dome_energy_Wh_sol, 2)

# Plot 1: Module 3 — REMS-Anchored Diurnal Ground Temp & Footpad TEG Power
fig, ax1 = plt.subplots(figsize=(10, 4.8))
ax1.plot(hours, ground_temp_K, color="#1f77b4", linewidth=2, label="REMS-anchored ground temp (K)")
ax1.axhline(REMS_GROUND_MEAN_K, color="#1f77b4", linestyle="--", alpha=0.6, label=f"REMS mean ({REMS_GROUND_MEAN_K} K)")
ax1.set_xlabel("Mars Local Solar Time (hours)")
ax1.set_ylabel("Ground Temperature (K)", color="#1f77b4")
ax1.tick_params(axis="y", labelcolor="#1f77b4")
ax1.grid(True, alpha=0.25)

ax2 = ax1.twinx()
ax2.plot(hours, footpad_teg_W, color="#d62728", linewidth=2, label="Footpad TEG output (W)")
ax2.set_ylabel("Estimated TEG Electrical Power (W)", color="#d62728")
ax2.tick_params(axis="y", labelcolor="#d62728")

plt.title("Module 3: REMS-Anchored Diurnal Ground Cycle & Footpad TEG Feasibility")
fig.tight_layout()
plt.savefig(OUT / "Module3_REMS_TEG_Feasibility.png", dpi=180)
plt.close()

# ---------------------------------------------------------
# MODULE 4: NEBUS HABITAT INTEGRATION & ENERGY BUDGET
# ---------------------------------------------------------
# Autonomous station daily energy demand (telemetry + radiation sensors + thermal keep-alive)
SUBSYSTEM_DEMAND_WH_SOL = {
    "RAD & dosimetry sensors": 6.5,
    "Onboard telemetry & MCU": 5.5,
    "Thermal telemetry & housekeeping": 4.0
}
TOTAL_DEMAND_WH_SOL = round(sum(SUBSYSTEM_DEMAND_WH_SOL.values()), 2)
ENERGY_MARGIN_WH_SOL = round(TOTAL_HYBRID_ENERGY_WH_SOL - TOTAL_DEMAND_WH_SOL, 2)

# Plot 2: Module 4 — Hybrid Generation vs. Autonomous Station Demand
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

# Left: Energy balance
categories = ["Footpad TEG\n(Ground-Coupled)", "Dome TEG\n(UV Photothermal)", "Total Hybrid\nGeneration", "Autonomous\nStation Demand"]
values = [footpad_energy_Wh_sol, dome_energy_Wh_sol, TOTAL_HYBRID_ENERGY_WH_SOL, TOTAL_DEMAND_WH_SOL]
colors = ["#2b5c8f", "#4a90e2", "#2ca02c", "#d95f02"]
bars = ax1.bar(categories, values, color=colors)
ax1.set_ylabel("Energy per Martian Sol (Wh/sol)")
ax1.set_title("Module 4: Daily Energy Budget (Generation vs. Demand)")
ax1.grid(axis="y", alpha=0.25)
for b in bars:
    ax1.text(b.get_x() + b.get_width()/2, b.get_height() + 0.3, f"{b.get_height():.2f}", ha="center", fontsize=9)

# Right: CAD Mass Allocation (Top 5 Subsystems + Others)
sorted_mass = sorted(MASS_BREAKDOWN_KG.items(), key=lambda x: x[1], reverse=True)
top_labels = [k.replace("_", " ") for k, _ in sorted_mass[:5]] + ["Other subsystems (5)"]
top_vals = [v for _, v in sorted_mass[:5]] + [sum(v for _, v in sorted_mass[5:])]
ax2.barh(top_labels[::-1], top_vals[::-1], color="#3b6998")
ax2.set_xlabel("Mass (kg)")
ax2.set_title(f"Module 4: CAD Dry Mass Allocation ({TOTAL_DRY_MASS_KG} kg)")
ax2.grid(axis="x", alpha=0.25)

fig.tight_layout()
plt.savefig(OUT / "Module4_Habitat_Integration_Budget.png", dpi=180)
plt.close()

# Export clean JSON summary for final report
final_summary = {
    "module_3_teg_feasibility": {
        "rems_mean_ground_temp_K": REMS_GROUND_MEAN_K,
        "tes_thermal_inertia_SI": TES_TI_OBSERVED,
        "hot_reservoir_source": "133.8 kg FR4 electronics waste heat buffered by 6.8 kg PCM ring",
        "footpad_teg_peak_W": round(footpad_peak_W, 3),
        "footpad_teg_avg_W": round(footpad_avg_W, 3),
        "footpad_energy_Wh_per_sol": round(footpad_energy_Wh_sol, 2),
        "dome_uv_teg_W": round(dome_teg_W, 3),
        "dome_uv_energy_Wh_per_sol": round(dome_energy_Wh_sol, 2),
        "total_hybrid_energy_Wh_per_sol": TOTAL_HYBRID_ENERGY_WH_SOL
    },
    "module_4_habitat_integration": {
        "total_dry_mass_kg": TOTAL_DRY_MASS_KG,
        "mars_total_weight_kN": TOTAL_WEIGHT_KN,
        "static_load_per_leg_kN": LOAD_PER_LEG_KN,
        "autonomous_station_demand_Wh_per_sol": TOTAL_DEMAND_WH_SOL,
        "net_daily_energy_margin_Wh_per_sol": ENERGY_MARGIN_WH_SOL,
        "operational_conclusion": (
            "Hybrid TEG (~20.8 Wh/sol) supports autonomous radiation telemetry and sensor keep-alive; "
            "primary crew habitat life-support remains supplied by Mars Base main power."
        )
    }
}

(OUT / "Module3_4_Final_Summary.json").write_text(json.dumps(final_summary, indent=2), encoding="utf-8")
print("SUCCESS! Generated Module 3 & Module 4 plots and JSON in:", OUT.resolve())
print(f"Total Hybrid Energy: {TOTAL_HYBRID_ENERGY_WH_SOL} Wh/sol | Station Demand: {TOTAL_DEMAND_WH_SOL} Wh/sol")
