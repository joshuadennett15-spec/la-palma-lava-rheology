"""
Fine-grained ER sweep for LaPalma1949 calibrated config at 12deg slope.
Targets ~6000m.
"""
import json, subprocess, re, os
import sys

TEMPLATE = "resource/LaPalma1949_calibrated.json"
PYTHON = sys.executable
TMP = "resource/_fine_er_tmp.json"
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

er_values = [36, 38, 40, 42, 44, 46, 48]

print(f"{'ER':>8} {'Dist(m)':>10} {'Note'}")
print("-" * 35)

best_er, best_dist, best_diff = None, None, float('inf')
for er in er_values:
    cfg = make_config(er)
    dist = run_flowgo(cfg)
    diff = abs(dist - TARGET)
    note = " <-- MATCH" if diff <= TOLERANCE else ""
    if diff < best_diff:
        best_diff, best_er, best_dist = diff, er, dist
    print(f"{er:>8.0f} {dist:>10.0f}{note}")

print(f"\nBest: ER={best_er} -> {best_dist:.0f}m (diff {best_diff:.0f}m)")
if os.path.exists(TMP):
    os.remove(TMP)
