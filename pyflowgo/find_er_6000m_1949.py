"""
Find the effusion rate for the La Palma 1949 flow to reach ~6000m.
Uses the calibrated config (LaPalma1949_calibrated.json) and sweeps ER.
"""
import json, subprocess, re, os
import sys

TEMPLATE = "resource/LaPalma1949_calibrated.json"
PYTHON = sys.executable
TMP = "resource/_find_er_6000_tmp.json"
TARGET = 6000
TOLERANCE = 200

def run_flowgo(cfg):
    r = subprocess.run([PYTHON, "main_flowgo.py", cfg], capture_output=True, text=True,
                       cwd=os.path.dirname(os.path.abspath(__file__)))
    dists = re.findall(r"distance from vent \(m\) = ([\d.]+)", r.stdout + r.stderr)
    return float(dists[-1]) if dists else 0.0

def make_config(er, out=TMP):
    with open(TEMPLATE, 'r') as f:
        c = json.load(f)
    c["effusion_rate_init"] = er
    c["terrain_conditions"]["max_channel_length"] = 8000.0
    with open(out, 'w') as f:
        json.dump(c, f, indent=2)
    return out

# Scan over a range of ER values
er_values = [22, 28, 35, 50, 75, 100, 150, 200, 300, 400, 500, 600, 700, 750, 800, 850, 875, 900]

print(f"{'ER':>8} {'Dist(m)':>10} {'Note'}")
print("-" * 35)

best_er = None
best_dist = None
best_diff = float('inf')

for er in er_values:
    cfg = make_config(er)
    dist = run_flowgo(cfg)
    diff = abs(dist - TARGET)
    note = ""
    if dist >= TARGET - TOLERANCE and dist <= TARGET + TOLERANCE:
        note = " <-- MATCH"
    if diff < best_diff:
        best_diff = diff
        best_er = er
        best_dist = dist
    print(f"{er:>8.1f} {dist:>10.0f}{note}")

print(f"\nBest ER = {best_er} m³/s -> {best_dist:.0f} m (target {TARGET}m, diff {best_diff:.0f}m)")

if os.path.exists(TMP):
    os.remove(TMP)
