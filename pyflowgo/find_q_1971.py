import json
import subprocess
import os
import pandas as pd

json_file = 'resource/LaPalma1971_calibrated.json'

def run_simulation(q):
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    data['effusion_rate_init'] = q
    
    with open('resource/tmp_1971_search.json', 'w') as f:
        json.dump(data, f, indent=4)
        
    # Run simulation
    subprocess.run(['python', 'main_flowgo.py', 'resource/tmp_1971_search.json'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # Read output
    output_file = os.path.join('results_flowgo', f"results_main_flowgo_La Palma 1971 Teneguia_{q}m3s.csv")
    if os.path.exists(output_file):
        df = pd.read_csv(output_file)
        return df['position'].max()
    return 0

# Binary search for Q that reaches ~2500m (between 2450 and 2550)
target = 2500
low = 12.0
high = 50.0

best_q = low
best_dist = 0

print("Searching for Q to reach ~2500m")
for _ in range(15):
    mid = (low + high) / 2
    dist = run_simulation(mid)
    print(f"Q: {mid:.2f} m3/s -> Distance: {dist} m")
    
    if abs(dist - target) < 20: # Within 20m of 2500 is good enough
        best_q = mid
        best_dist = dist
        break
        
    if dist < target:
        low = mid
    else:
        high = mid

print(f"\\nFound optimal Q: {best_q:.2f} m3/s (Distance: {best_dist} m)")
