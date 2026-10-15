"""
Multi-parameter calibration for La Palma 2021 Sept 29 — target ~5880m (coast).
Varies effusion rate, width, crust cover, and core-base distance.
"""
import json, subprocess, re, os, itertools
import sys

TEMPLATE = "resource/LaPalma2021_Sept29.json"
PYTHON = sys.executable
TMP = "resource/_cal_sept29_tmp.json"
TARGET = 5880
TOLERANCE = 100

def run_flowgo(cfg):
    r = subprocess.run([PYTHON, "main_flowgo.py", cfg], capture_output=True, text=True,
                       cwd=os.path.dirname(os.path.abspath(__file__)))
    dists = re.findall(r"distance from vent \(m\) = ([\d.]+)", r.stdout + r.stderr)
    return float(dists[-1]) if dists else 0.0

def make_config(p, out=TMP):
    with open(TEMPLATE, 'r') as f:
        c = json.load(f)
    c["effusion_rate_init"] = p["er"]
    c["terrain_conditions"]["width"] = p["w"]
    c["thermal_parameters"]["crust_cover_fraction"] = p["ccf"]
    c["conduction_parameters"]["core_base_distance"] = p["cbd"]
    with open(out, 'w') as f:
        json.dump(c, f, indent=2)
    return out

# Realistic-ish parameter grid for peak 2021 eruption
# Sept 29 was Day 10, very high effusion, tube-fed in places
param_grid = {
    "er":  [30, 50, 75, 100],
    "w":   [5, 8, 12],
    "ccf": [0.9, 0.95, 1.0],
    "cbd": [8, 12, 15, 20],
}

keys = list(param_grid.keys())
combos = list(itertools.product(*[param_grid[k] for k in keys]))
total = len(combos)

print(f"TARGET: {TARGET}m | Combos: {total}")
print("=" * 70)

results = []
for i, vals in enumerate(combos):
    p = dict(zip(keys, vals))
    cfg = make_config(p)
    dist = run_flowgo(cfg)
    diff = abs(dist - TARGET)
    results.append((dist, diff, p.copy()))
    
    if diff <= TOLERANCE:
        print(f"  [{i+1:3d}/{total}] ER={p['er']:4.0f} W={p['w']:2.0f} CCF={p['ccf']:.2f} CBD={p['cbd']:2.0f} -> {dist:6.0f}m  *** MATCH ***")
    elif dist > 2000:
        print(f"  [{i+1:3d}/{total}] ER={p['er']:4.0f} W={p['w']:2.0f} CCF={p['ccf']:.2f} CBD={p['cbd']:2.0f} -> {dist:6.0f}m")

results.sort(key=lambda x: x[1])

print(f"\n{'=' * 70}")
print(f"TOP 15 CLOSEST TO {TARGET}m")
print(f"{'=' * 70}")
print(f"{'Rank':>4} {'Dist':>7} {'Diff':>6} {'ER':>5} {'W':>4} {'CCF':>5} {'CBD':>4}")
print("-" * 45)
for rank, (dist, diff, p) in enumerate(results[:15], 1):
    marker = " <-- BEST" if rank == 1 else ""
    print(f"{rank:4d} {dist:7.0f} {diff:6.0f} {p['er']:5.0f} {p['w']:4.0f} {p['ccf']:5.2f} {p['cbd']:4.0f}{marker}")

best_dist, best_diff, best_p = results[0]
make_config(best_p, TEMPLATE)
print(f"\nUpdated {TEMPLATE}")
print(f"Best: ER={best_p['er']}, W={best_p['w']}, CCF={best_p['ccf']}, CBD={best_p['cbd']}")
print(f"Flow: {best_dist:.0f}m (target {TARGET}m)")

if os.path.exists(TMP):
    os.remove(TMP)
