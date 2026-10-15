"""
Multi-parameter calibration for La Palma 1949 — target 1200m.
Varies effusion rate, width, crust cover, and core-base distance
within realistic bounds. Width fixed at 12m per user.
"""
import json, subprocess, re, os, itertools
import sys

TEMPLATE = "resource/LaPalma1949.json"
PYTHON = sys.executable
TMP = "resource/_cal_1949_tmp.json"
TARGET = 1200
TOLERANCE = 30

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

# Realistic parameter grid (width fixed at 12)
param_grid = {
    "er":  [10, 14, 18, 22, 28],
    "w":   [12],
    "ccf": [0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    "cbd": [2, 5, 8, 10, 15],
}

keys = list(param_grid.keys())
combos = list(itertools.product(*[param_grid[k] for k in keys]))
total = len(combos)

print("=" * 70)
print(f"TARGET: {TARGET}m | Width: 12m | Combos: {total}")
print("=" * 70)

results = []
for i, vals in enumerate(combos):
    p = dict(zip(keys, vals))
    cfg = make_config(p)
    dist = run_flowgo(cfg)
    diff = abs(dist - TARGET)
    results.append((dist, diff, p.copy()))
    
    if diff <= TOLERANCE:
        print(f"  [{i+1:3d}/{total}] ER={p['er']:4.0f} CCF={p['ccf']:.1f} CBD={p['cbd']:2.0f} -> {dist:6.0f}m  *** MATCH ***")
    elif dist > 500:
        print(f"  [{i+1:3d}/{total}] ER={p['er']:4.0f} CCF={p['ccf']:.1f} CBD={p['cbd']:2.0f} -> {dist:6.0f}m")

# Sort by closeness to target
results.sort(key=lambda x: x[1])

print(f"\n{'=' * 70}")
print(f"TOP 15 CLOSEST TO {TARGET}m")
print(f"{'=' * 70}")
print(f"{'Rank':>4} {'Dist':>7} {'Diff':>6} {'ER':>5} {'CCF':>5} {'CBD':>4}")
print("-" * 40)
for rank, (dist, diff, p) in enumerate(results[:15], 1):
    marker = " <-- BEST" if rank == 1 else ""
    print(f"{rank:4d} {dist:7.0f} {diff:6.0f} {p['er']:5.0f} {p['ccf']:5.1f} {p['cbd']:4.0f}{marker}")

# Save best
best_dist, best_diff, best_p = results[0]
make_config(best_p, TEMPLATE)
print(f"\nUpdated {TEMPLATE}: ER={best_p['er']}, CCF={best_p['ccf']}, CBD={best_p['cbd']}")
print(f"Flow length: {best_dist:.0f}m (target {TARGET}m, diff {best_diff:.0f}m)")

if os.path.exists(TMP):
    os.remove(TMP)
