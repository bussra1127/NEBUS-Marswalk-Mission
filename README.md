# NEBUS — Integrated Martian Terrain Assessment and Autonomous Surface Station

**NEBUS** is an integrated planetary engineering and preliminary design (PDR) concept developed to support future human surface operations in **Gale Crater, Mars (`4.57° S, 137.95° E`)**.

## Project Overview

* **Module 1 — Martian Terrain & EVA Route Assessment:** Evaluates Gale Crater topography, surface slope (`<15°` safe corridor), surface roughness, and MGS TES thermal inertia (`405 SI`, neighborhood median `373 SI`) alongside NASA Curiosity RAD radiation observations (`601 records`) and REMS thermal baselines (`211.45 K`).
* **Module 2 — Autonomous Surface Station & Hybrid Energy Feasibility:** Presents a 51-body, six-legged autonomous station (`7,883.1 kg` dry mass; `~4.88 kN/leg` static load under `3.71 m/s²` Martian gravity) featuring a multilayer protective dome (`Fused Silica` / `12 wt% Melanin + 8 wt% GO + 80 wt% H-UHMWPE` / `Silica Aerogel`) and Ti-6Al-4V/B-UHMWPE/aerogel landing supports.
* **Structural & Energy Verification:** SimScale dynamic FEA confirms low stress levels (`0.6324 MPa` max Von Mises stress, `0.27 mm` max displacement vs. `~880 MPa` Ti-6Al-4V yield strength). Parametric hybrid TEG modeling yields `19.98 Wh/sol` (`12.47 Wh/sol` footpad ground-coupled + `7.51 Wh/sol` dome UV-photothermal) against a `16.00 Wh/sol` autonomous telemetry demand (`+3.98 Wh/sol` margin).

## Repository Contents

* `nebus_module2_nasa_data.py` — NASA RAD & TES data ingestion, composite hydrogen comparison, and parametric sensitivity analysis.
* `nebus_tes_to_teg.py` — MGS TES thermal inertia (`405 SI`) extraction and diurnal ground-coupled TEG modeling.
* `nebus_module3_4_final.py` — REMS-anchored (`211.45 K`) thermal cycle, hybrid TEG feasibility, and CAD mass/energy budget integration.
* `NEBUS_module2_results.json` & `Module3_4_Final_Summary.json` — Structured numerical outputs separating empirical NASA data from parametric engineering models.
* `README_LIMITATIONS.txt` — Explicit verification boundaries and model assumptions.
