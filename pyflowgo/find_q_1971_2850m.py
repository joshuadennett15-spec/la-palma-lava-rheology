"""
Binary search for the effusion rate Q that makes the 1971 Teneguia flow reach ~2850m,
using 2m step size.
"""
import json
import sys
import subprocess
import os
import pandas as pd

TEMPLATE = 'resource/LaPalma1971_calibrated.json'
TMP = 'resource/tmp_1971_2850_search.json'
PYTHON = sys.executable
RESULTS_DIR = 'results_flowgo'
LAVA_NAME = 'La Palma 1971 Teneguia'

TARGET = 2850
TOLERANCE = 20  # accept within ±20m


def run_simulation(q):
    with open(TEMPLATE, 'r') as f:
        data = json.load(f)

    data['effusion_rate_init'] = q
    data['step_size'] = 2.0
    data['terrain_conditions']['max_channel_length'] = 3000.0

    with open(TMP, 'w') as f:
        json.dump(data, f, indent=4)

    subprocess.run(
        [PYTHON, 'main_flowgo.py', TMP],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        cwd=os.path.dirname(os.path.abspath(__file__))
    )

    output_file = os.path.join(RESULTS_DIR, f'results_main_flowgo_{LAVA_NAME}_{q}m3s.csv')
    if os.path.exists(output_file):
        df = pd.read_csv(output_file)
        return df['position'].max()
    return 0


low = 13.8
high = 80.0

best_q = low
best_dist = 0

print(f"Binary search for Q that gives ~{TARGET}m runout at 2m step size")
print(f"Search range: Q=[{low}, {high}] m³/s\n")

for i in range(20):
    mid = (low + high) / 2
    dist = run_simulation(mid)
    print(f"  Iter {i+1:2d}: Q={mid:.4f} m³/s  ->  runout={dist:.0f} m")

    if abs(dist - TARGET) < TOLERANCE:
        best_q = mid
        best_dist = dist
        print(f"\nConverged! Q={best_q:.4f} m³/s gives {best_dist:.0f}m (target {TARGET}m)")
        break

    if dist < TARGET:
        low = mid
    else:
        high = mid
else:
    # Take the closest result
    best_q = (low + high) / 2
    best_dist = run_simulation(best_q)
    print(f"\nMax iterations reached. Best Q={best_q:.4f} m³/s -> {best_dist:.0f}m")

print(f"\nFinal: Q = {best_q:.4f} m³/s  (runout = {best_dist:.0f} m, target = {TARGET} m)")

if os.path.exists(TMP):
    os.remove(TMP)
