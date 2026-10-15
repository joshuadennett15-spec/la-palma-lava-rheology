"""
Binary search for the effusion rate that produces a ~980m flow length
for the La Palma 2021 Dec 13 terminal phase configuration.
"""
import json, subprocess, re, os
import sys

TEMPLATE = "resource/LaPalma2021_Dec13.json"
PYTHON = sys.executable
TARGET = 980
TOLERANCE = 10  # accept within +/- 10m

def run_flowgo(cfg):
    r = subprocess.run([PYTHON, "main_flowgo.py", cfg], capture_output=True, text=True,
                       cwd=os.path.dirname(os.path.abspath(__file__)))
    dists = re.findall(r"distance from vent \(m\) = ([\d.]+)", r.stdout + r.stderr)
    return float(dists[-1]) if dists else 0.0

def set_effusion_rate(er, out="resource/_cal_2021_tmp.json"):
    with open(TEMPLATE, 'r') as f:
        c = json.load(f)
    c["effusion_rate_init"] = round(er, 2)
    with open(out, 'w') as f:
        json.dump(c, f, indent=2)
    return out

print("=" * 60)
print(f"TARGET: {TARGET}m flow length")
print("=" * 60)

# Binary search between known bounds: 5 m3/s -> 166m, 28 m3/s -> 1487m
lo, hi = 5.0, 28.0
best_er, best_dist, best_diff = None, 0, float('inf')

for iteration in range(30):
    mid = (lo + hi) / 2.0
    cfg = set_effusion_rate(mid)
    dist = run_flowgo(cfg)
    diff = abs(dist - TARGET)
    
    print(f"  [{iteration+1:2d}] ER={mid:6.2f} m3/s -> {dist:6.0f}m (diff={diff:.0f}m)")
    
    if diff < best_diff:
        best_diff = diff
        best_dist = dist
        best_er = mid
    
    if diff <= TOLERANCE:
        print(f"\n  *** CONVERGED ***")
        break
    
    if dist < TARGET:
        lo = mid
    else:
        hi = mid

print(f"\n{'=' * 60}")
print(f"BEST: effusion_rate = {best_er:.2f} m3/s -> {best_dist:.0f}m (target {TARGET}m)")
print(f"{'=' * 60}")

# Save the best config
set_effusion_rate(best_er, TEMPLATE)
print(f"\nUpdated {TEMPLATE} with effusion_rate_init = {best_er:.2f}")

# Cleanup
tmp = "resource/_cal_2021_tmp.json"
if os.path.exists(tmp):
    os.remove(tmp)
