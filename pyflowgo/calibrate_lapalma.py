"""
Re-calibration with literature-corrected parameters for La Palma 1971.
Fixed: vesicle_fraction=0.15, einstein_coef=2.5, max_packing=0.62
Varies: effusion_rate, width, depth, eruption_temp, core_base_distance, crust_cover_fraction
"""
import json, subprocess, re, os, sys

TEMPLATE = "resource/LaPalma1971.json"
PYTHON = sys.executable
TARGET = 2000

def run_flowgo(cfg):
    r = subprocess.run([PYTHON, "main_flowgo.py", cfg], capture_output=True, text=True,
                       cwd=os.path.dirname(os.path.abspath(__file__)))
    dists = re.findall(r"distance from vent \(m\) = ([\d.]+)", r.stdout + r.stderr)
    return float(dists[-1]) if dists else 0.0

def make_config(p, out):
    with open(TEMPLATE, 'r') as f: c = json.load(f)
    c["effusion_rate_init"] = p["er"]
    c["terrain_conditions"]["width"] = p["w"]
    c["terrain_conditions"]["depth"] = p["d"]
    c["eruption_condition"]["eruption_temperature"] = p["t"]
    c["eruption_condition"]["viscosity_eruption"] = p["v"]
    c["lava_state"]["density_dre"] = p["dre"]
    c["lava_state"]["vesicle_fraction"] = 0.15  # FIXED per literature
    c["lava_state"]["crystal_fraction"] = p["xf"]
    c["conduction_parameters"]["core_base_distance"] = p["cbd"]
    c["thermal_parameters"]["crust_cover_fraction"] = p["ccf"]
    c["relative_viscosity_parameters"]["max_packing"] = 0.62   # FIXED
    c["relative_viscosity_parameters"]["einstein_coef"] = 2.5   # FIXED
    with open(out, 'w') as f: json.dump(c, f, indent=2)

tmp = "resource/_cal_tmp.json"
best_diff = float('inf')
best_dist = 0
best_p = None

print("=" * 70)
print("CALIBRATION — Literature-corrected parameters")
print(f"Fixed: vesicle_fraction=0.15, einstein_coef=2.5, max_packing=0.62")
print(f"Target: {TARGET}m")
print("=" * 70)

# Coarse grid
for er in [5, 10, 15, 19.3, 25]:
    for w in [5, 10, 15, 20, 30, 41]:
        for t in [1373.15, 1413.15, 1448.15]:
            for d in [1.5, 2.5]:
                p = {"er":er, "w":w, "d":d, "t":t, "v":100, "dre":2600, "xf":0.10, "cbd":10.0, "ccf":0.8}
                make_config(p, tmp)
                dist = run_flowgo(tmp)
                diff = abs(dist - TARGET)
                if diff < best_diff:
                    best_diff = diff; best_dist = dist; best_p = p.copy()
                if dist > 1500 or diff < 500:
                    tag = " ***BEST***" if diff == best_diff else ""
                    print(f"  ER={er:5.1f} W={w:4.0f} D={d:.1f} T={t:.0f} -> {dist:6.0f}m{tag}")

print(f"\nCoarse best: {best_dist:.0f}m (diff {best_diff:.0f}m)")
print(f"  Params: {best_p}")

# Fine-tune around best
print("\nFine-tuning...")
eb = best_p["er"]; wb = best_p["w"]; db = best_p["d"]; tb = best_p["t"]
for er in [eb*0.7, eb*0.8, eb*0.9, eb, eb*1.1, eb*1.2, eb*1.3]:
    for w in [wb*0.7, wb*0.8, wb*0.9, wb, wb*1.1, wb*1.2, wb*1.3]:
        for cbd in [5, 10, 15, 19]:
            for ccf in [0.5, 0.7, 0.8, 0.9, 1.0]:
                for v in [100, 300]:
                    p = {**best_p, "er":round(er,1), "w":round(w,1), "cbd":cbd, "ccf":ccf, "v":v}
                    make_config(p, tmp)
                    dist = run_flowgo(tmp)
                    diff = abs(dist - TARGET)
                    if diff < best_diff:
                        best_diff = diff; best_dist = dist; best_p = p.copy()
                        print(f"  NEW BEST: ER={er:.1f} W={w:.1f} CBD={cbd} CCF={ccf} V={v} -> {dist:.0f}m (diff={diff:.0f}m)")

print(f"\n{'='*70}")
print(f"FINAL: {best_dist:.0f}m (target {TARGET}m)")
for k,v in sorted(best_p.items()): print(f"  {k}: {v}")
print(f"{'='*70}")

# Save
make_config(best_p, "resource/LaPalma1971_calibrated.json")
print(f"\nSaved to resource/LaPalma1971_calibrated.json")
if os.path.exists(tmp): os.remove(tmp)
