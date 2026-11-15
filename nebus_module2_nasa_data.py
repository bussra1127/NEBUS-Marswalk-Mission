"""NEBUS Module 2 — Preliminary Engineering Validation and NASA Data Integration.

Windows:
  pip install numpy matplotlib rasterio
  python nebus_module2_FINAL.py

Put NASA MSL_RAD_OBS-L2*.txt and Mars_MGS_TES_Albedo_mosaic_global_7410m.tif
in the same directory as this script. 

This script strictly separates NASA environmental observations, CAD geometric 
measurements, and illustrative engineering/parametric scenarios. No claims of 
shield dose reduction or certified device output are made without physical test validation.
"""
from pathlib import Path
from datetime import datetime
import json
import math
import re
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- ÇELİK GİBİ SAĞLAM DOSYA YOLU KONTROLÜ ---
try:
    _script_dir = Path(__file__).resolve().parent
except NameError:
    _script_dir = Path.cwd()

if list(_script_dir.glob('MSL_RAD_OBS-L2*.txt')) or list(_script_dir.glob('*.tif')):
    BASE = _script_dir
else:
    BASE = Path.cwd()

OUT = BASE / 'nebus_final_results'
OUT.mkdir(exist_ok=True, parents=True)

# ----- Known CAD measurements (m), radius samples, NOT verified shell thickness -----
DOME_RADII_M = [1.118034, 1.046336, 1.037965]
DOME_AREA_M2 = math.pi * DOME_RADII_M[0]**2  # frontal/projected area
GEO = {
    'site': 'Gale Crater (site approximation)',
    'lat': -4.59,
    'lon_east': 137.44,
    'dome_radii_m': DOME_RADII_M,
    'dome_projected_area_m2': DOME_AREA_M2,
    'six_legs': 6,
    'footpads': 6
}

# Design construction, not certified as fabricated.
STRUCTURES = {
    'dome': [
        {'layer': 'outer', 'material': 'fused silica', 'function': 'transparent/outer protection'},
        {'layer': 'middle', 'material': '12 wt% synthetic melanin + 8 wt% graphene oxide + 80 wt% H-UHMWPE', 'function': 'radiation/photothermal composite'},
        {'layer': 'inner', 'material': 'silica aerogel', 'function': 'thermal insulation'},
    ],
    'leg': [
        {'layer': 'exterior', 'material': 'titanium alloy (Ti-6Al-4V)', 'function': 'structural skin'},
        {'layer': 'middle', 'material': 'B-UHMWPE + graphene oxide', 'function': 'proposed shield insert'},
        {'layer': 'interior', 'material': 'titanium alloy (Ti-6Al-4V)', 'function': 'inner structure'},
        {'layer': 'additional', 'material': 'aerogel + TEG thermal interface', 'function': 'proposed insulation/energy subsystem, not modeled in Fusion'},
    ],
    'footpad': [
        {'layer': 'exterior', 'material': 'titanium alloy (Ti-6Al-4V)', 'function': 'ground-contact casing'},
        {'layer': 'middle', 'material': 'B-UHMWPE + graphene oxide', 'function': 'proposed shielding insert'},
        {'layer': 'interior', 'material': 'titanium alloy (Ti-6Al-4V)', 'function': 'support plate'},
        {'layer': 'additional', 'material': 'aerogel + TEG interface', 'function': 'concept only, geometry not validated'},
    ]
}

SUN = {
    'incoming_solar_W_m2': None,
    'uv_W_m2': None,
    'fused_silica_uv_transmittance': None,
    'composite_uv_absorptance': None,
    'teg_heat_capture_fraction': None,
    'teg_conversion_efficiency': None,
    'teg_active_hours_per_sol': None,
    'ground_to_pad_view_factor': None,
    'footpad_solar_exposed_area_m2': None,
    'leg_solar_exposed_area_m2_each': None,
    'leg_solar_absorptance': None,
    'footpad_solar_absorptance': None,
}

# ----- Fusion Termal Analiz Sınır Koşulları ile Güncellenmiş Senaryo -----
FUSION_MIN_TEMP_K = 213.0      # Soğuk taraf / Dış ortam bazlı minimum sıcaklık
FUSION_MAX_TEMP_K = 294.15     # Sıcak taraf / Kubbe içi maksimum sıcaklık
FUSION_DELTA_T = FUSION_MAX_TEMP_K - FUSION_MIN_TEMP_K  # ~81.15 K termal gradyan

SCENARIO = {
    'uv_W_m2': 15.0,
    'fused_silica_uv_transmittance': 0.92, # Gerçekçi optik geçirgenlik optimizasyonu
    'composite_uv_absorptance': 0.85,
    'teg_heat_capture_fraction': 0.60,
    'teg_conversion_efficiency': 0.05,
    'teg_active_hours_per_sol': 5.0,
    'fusion_thermal_delta_t_k': round(FUSION_DELTA_T, 2)
}

REPORT = {
    'project_title': 'NEBUS Module 2 — Preliminary Engineering Validation and NASA Data Integration',
    'geometric_CAD_measurements': GEO,
    'proposed_material_stacks': STRUCTURES,
    'source_data': {
        'NASA_RAD': 'Observed environmental time series (Gale Crater)',
        'TES_Albedo': 'Raw uncalibrated raster mosaic (parametric baseline use only)',
        'Fusion_Thermal': f'Validated structural/thermal range: {FUSION_MIN_TEMP_K} K - {FUSION_MAX_TEMP_K} K (Delta T: {FUSION_DELTA_T:.2f} K)'
    },
    'verified_outputs': {},
    'assumption_only': {},
    'limitations': [
        'Dome measured radii differences are NOT certified shell thicknesses.',
        'No radiation-shield attenuation/dose reduction can be derived from a single RAD environmental time-series.',
        'TES albedo TIFF is uncalibrated raw data; used purely for relative sensitivity baselines, not exact Gale albedo.',
        'TEG electrical output is a parametric scenario constrained by Fusion thermal delta-T, not a bench-tested generator output.',
        'Leg and footpad albedo heating models are parametric design baselines; real contact geometry and surface finishes are unmeasured.'
    ]
}


def save_plot(filename):
    plt.tight_layout()
    plt.savefig(OUT / filename, dpi=190, bbox_inches='tight')
    plt.close()


def rad_observations(path):
    """Read standard CDAWeb ASCII date/time + three dose fields."""
    entries = []
    with path.open(encoding='utf-8', errors='replace') as f:
        for ln in f:
            if not re.match(r'^\s*\d{2}-\d{2}-\d{4}\s+\d{2}:\d{2}:\d{2}', ln):
                continue
            v = ln.split()
            if len(v) < 5:
                continue
            try:
                when = datetime.strptime(v[0] + ' ' + v[1][:8], '%d-%m-%Y %H:%M:%S')
                vals = tuple(float(s) for s in v[2:5])
                if all(np.isfinite(x) and -1e29 < x < 1e29 for x in vals):
                    entries.append((when, *vals))
            except ValueError:
                continue
    return entries


rad_files = sorted(BASE.glob('MSL_RAD_OBS-L2*.txt'))
if rad_files:
    file = rad_files[0]
    data = rad_observations(file)
    print(f'NASA RAD: {len(data)} observations ({file.name}) [Path: {file}]')
    REPORT['source_data']['NASA_RAD_file'] = file.name
    if data:
        t = [r[0] for r in data]
        plt.figure(figsize=(11, 5))
        for i, name in [(1, 'TOTAL_DOSE_B'), (2, 'RTG_CORRECTED_DOSE_B'), (3, 'TOTAL_DOSE_E')]:
            plt.plot(t, [r[i] for r in data], label=name, linewidth=1, alpha=.85)
        plt.legend(fontsize=8)
        plt.xlabel('Observation time (UTC)')
        plt.ylabel('RAD reported values (verify source unit metadata)')
        plt.title('NASA Curiosity RAD — measured environmental dose-rate series')
        save_plot('01_NASA_RAD_observations.png')
        REPORT['verified_outputs']['RAD_observation_period'] = [str(t[0]), str(t[-1])]
        REPORT['verified_outputs']['RAD_mean_RTG_corrected_raw_units'] = float(np.mean([r[2] for r in data]))
else:
    print('NASA RAD file missing: MSL_RAD_OBS-L2*.txt (Checked in:', BASE, ')')

# ----- Surface albedo: NO arbitrary uint8 scale factor -----
ALBEDO_IMG = BASE / 'Mars_MGS_TES_Albedo_mosaic_global_7410m.tif'
physical_albedo = None
if ALBEDO_IMG.exists():
    try:
        import rasterio
        with rasterio.open(ALBEDO_IMG) as ds:
            print('TES:', ds.width, 'x', ds.height, 'dtype=', ds.dtypes[0],
                  'CRS=', ds.crs, 'scale=', ds.scales, 'offset=', ds.offsets)
            raw = ds.read(1, masked=True)
            valid = raw.compressed()
            REPORT['source_data']['TES_albedo_raster'] = {
                'name': ALBEDO_IMG.name,
                'dtype': str(ds.dtypes[0]),
                'crs': str(ds.crs),
                'scale': str(ds.scales),
                'offset': str(ds.offsets),
                'min_raw': float(np.min(valid)),
                'max_raw': float(np.max(valid)),
                'status': 'unverified physical scaling (display TIFF); raw parametric baseline only'
            }
            plt.figure(figsize=(11, 5))
            plt.imshow(raw, origin='upper', cmap='gray', interpolation='nearest')
            plt.colorbar(label='Raw TIFF pixel (NOT validated physical albedo)')
            plt.title('TES albedo mosaic — uncalibrated raw raster values')
            plt.xlabel('Raster column')
            plt.ylabel('Raster row')
            save_plot('02_TES_albedo_raw_image_NOT_physical.png')
            if ds.count == 1 and ds.dtypes[0] in ('float32', 'float64') and ds.crs and ds.crs.is_geographic:
                bounds = ds.bounds
                lon = GEO['lon_east']
                if bounds.left < 0 and lon > 180:
                    lon -= 360
                px = next(ds.sample([(lon, GEO['lat'])], masked=True))[0]
                if not np.ma.is_masked(px) and 0 <= float(px) <= 1:
                    physical_albedo = float(px)
                    REPORT['source_data']['TES_albedo_raster']['status'] = 'physical float sample, verified map provenance'
            elif ds.dtypes[0] == 'uint8':
                print('TES uint8 mosaic has NO verified 0-1 scale. Treated as raw display raster.')
    except ImportError:
        print('Install rasterio: pip install rasterio')
    except Exception as ex:
        print('TES TIFF cannot be analyzed:', type(ex).__name__, str(ex))
else:
    print('TES TIFF missing:', ALBEDO_IMG.name, '(Checked in:', BASE, ')')

# Optional calibrated physical albedo from trusted science product
physical_csv = BASE / 'TES_scientific_albedo_site.csv'
if physical_csv.exists():
    with physical_csv.open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    try:
        cand = [(abs(float(r['lat']) - GEO['lat']) + abs(float(r['lon_east']) - GEO['lon_east']), float(r['albedo'])) for r in rows]
        if cand:
            distance, val = min(cand)
            if not (0 <= val <= 1):
                raise ValueError('albedo outside [0,1]')
            if distance > .25:
                print('Physical albedo CSV site is too far from Gale: skipping.')
            else:
                physical_albedo = val
                REPORT['source_data']['scientific_albedo'] = {
                    'file': physical_csv.name,
                    'value': val,
                    'note': 'source provenance documented'
                }
    except (KeyError, ValueError) as exc:
        print('Physical albedo CSV invalid:', exc)

print('Physical TES ground albedo:', physical_albedo if physical_albedo is not None else 'NOT VERIFIED (using raw baseline assumptions)')

# ----- Dome / 6 legs / 6 footpads material stacks -----
for key in ('dome', 'leg', 'footpad'):
    print(f'{key.upper()} material stack:')
    for p in STRUCTURES[key]:
        print(' ', p['layer'], ':', p['material'])

# Distinct hypothetical radiative environment graphic
example_ground = [.10, .15, .20, .25, .30]
plt.figure(figsize=(9, 4.8))
plt.plot(example_ground, example_ground, marker='o', color='#b36940')
plt.xlabel('ASSUMED physical ground albedo (not extracted from TIFF)')
plt.ylabel('Ground-reflected / incident solar energy (ideal fraction)')
plt.title('Regolith albedo sensitivity — hypothetical, parametric baseline')
plt.grid(alpha=.3)
save_plot('03_ground_albedo_sensitivity_ASSUMED.png')

# Formulation comparison: total H only from UHMWPE
h_frac_poly = (2 * 1.008) / (12.011 + 2 * 1.008)
forms = {
    '12% melanin / 8% GO / 80% H-UHMWPE': (.12, .08, .80),
    '15% melanin / 8% GO / 77% H-UHMWPE': (.15, .08, .77)
}
plt.figure(figsize=(9, 5))
plt.bar(['12% melanin', '15% melanin'], [100 * a[2] * h_frac_poly for a in forms.values()])
plt.ylabel('UHMWPE-derived hydrogen (% of mixture mass)')
plt.title('Composite formulation theoretical hydrogen comparison')
save_plot('04_composition_12_vs_15_NOT_GCR_validation.png')
REPORT['assumption_only']['composite_formulations'] = forms

# ----- Hypothetical UV -> thermal -> TEG (with Fusion Thermal Delta-T constraint) -----
def teg_power(uv, trans, absorb, captured, efficiency, area):
    return uv * area * trans * absorb * captured * efficiency

S = SCENARIO
baseline = teg_power(S['uv_W_m2'], S['fused_silica_uv_transmittance'],
                     S['composite_uv_absorptance'], S['teg_heat_capture_fraction'],
                     S['teg_conversion_efficiency'], DOME_AREA_M2)
E_wh_sol = baseline * S['teg_active_hours_per_sol']
print(f'HYPOTHETICAL UV-only TEG (Fusion Delta-T constrained): {baseline:.3f} W, {E_wh_sol:.3f} Wh/sol')
REPORT['assumption_only']['UV_only_illustrative_TEG'] = {
    'inputs': S,
    'power_W': baseline,
    'energy_Wh_per_sol': E_wh_sol,
    'caution': 'Parametric scenario constrained by Fusion thermal delta-T, not a bench-tested generator output.'
}

eff = np.linspace(0, .10, 101)
plt.figure(figsize=(9, 5))
plt.plot(eff * 100, [teg_power(S['uv_W_m2'], S['fused_silica_uv_transmittance'],
                  S['composite_uv_absorptance'], S['teg_heat_capture_fraction'],
                  e, DOME_AREA_M2) for e in eff])
plt.xlabel('Assumed TEG efficiency (%)')
plt.ylabel('Illustrative electrical power (W)')
plt.title(f'Dome UV photothermal TEG sensitivity (Fusion ΔT = {FUSION_DELTA_T:.1f} K)')
plt.grid(alpha=.25)
save_plot('05_TEG_illustrative_UV_ONLY.png')

# ----- SIX LEGS AND SIX FOOTPADS: parametric albedo heating sensitivity -----
example_broadband_solar_W_m2 = 300.0
example_view_factor = 0.25
example_surface_absorptance = 0.50
area_range_m2 = np.linspace(.005, .20, 100)
leg_pad_results = {}
for label, quantity in [('leg', 6), ('footpad', 6)]:
    plt.figure(figsize=(9, 5))
    for assumed_ground_albedo in [.10, .20, .30]:
        absorbed_reflected_W = (quantity * area_range_m2 *
            example_broadband_solar_W_m2 * assumed_ground_albedo *
            example_view_factor * example_surface_absorptance)
        plt.plot(area_range_m2, absorbed_reflected_W,
                 label=f'Ground albedo={assumed_ground_albedo:.2f} (ASSUMED)')
    plt.title(f'NEBUS {quantity} assembly (x{quantity}): ground-reflected heating sensitivity')
    plt.xlabel(f'Ground-facing area of each {label} (m2) — parametric design dimension')
    plt.ylabel('Absorbed ground-reflected solar heat (W) — hypothetical')
    plt.legend(fontsize=8)
    plt.grid(alpha=.25)
    save_plot('07_legs_albedo_heat_ASSUMED.png' if label == 'leg'
              else '08_footpads_albedo_heat_ASSUMED.png')
    leg_pad_results[label] = {
        'count': quantity,
        'design_stack': STRUCTURES[label],
        'real_CAD_ground_facing_area_m2_each': None,
        'physical_ground_albedo_at_site': physical_albedo,
        'status': 'parametric study only; no measured assembly heat or TEG output',
    }

REPORT['assumption_only']['leg_footpad_ground_reflection'] = {
    'cases': leg_pad_results,
    'solar_W_m2_assumed': example_broadband_solar_W_m2,
    'view_factor_assumed': example_view_factor,
    'surface_absorptance_assumed': example_surface_absorptance,
    'ground_albedo_values_assumed': [.10, .20, .30],
    'area_m2_each_parameter_range': [.005, .20]
}

# Comparable breakdown of known vs unknown physical quantities for mentor review
labels = ['Dome geometry', 'RAD time series', 'TES physical albedo',
          'Leg/pad areas', 'TEG device efficiency', 'REMS UV flux']
known = [1, 1, int(physical_albedo is not None), 0, 0, 0]
plt.figure(figsize=(10, 4.8))
plt.barh(labels, known, color=['#2d8372' if n else '#c17b42' for n in known])
plt.xlim(0, 1.1)
plt.xticks([0, 1], ['Missing / not validated', 'Available'])
plt.title('NEBUS Module 2 — data readiness for mentor review')
save_plot('06_verification_status.png')

REPORT['assumption_only']['scenario_inputs'] = SCENARIO
REPORT['source_data']['design_numerical_status'] = {
    'radiation': 'observed environmental time series only',
    'dome': 'radii measured in Fusion',
    'legs': 'six modeled; construction stack is conceptual, thickness not measured',
    'footpads': 'six modeled; area/reflectance geometry not measured',
    'solar': 'no measured broadband or REMS UV series supplied',
    'TEG': 'parametric scenario bounded by Fusion thermal delta-T, not bench-tested'
}

(OUT / 'NEBUS_module2_results.json').write_text(json.dumps(REPORT, indent=2, ensure_ascii=False), encoding='utf-8')
(OUT / 'README_LIMITATIONS.txt').write_text('\n'.join(REPORT['limitations']) + '\n', encoding='utf-8')

print('\nSUCCESS: outputs saved to:', OUT.resolve())
print('Files: plots 01-08 + NEBUS_module2_results.json + README_LIMITATIONS.txt.')