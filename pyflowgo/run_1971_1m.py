import json, subprocess, os, pandas as pd
import sys

TEMPLATE = 'resource/LaPalma1971_calibrated.json'
TMP = 'resource/tmp_1971_1m.json'
PYTHON = sys.executable

with open(TEMPLATE) as f:
    data = json.load(f)

data['effusion_rate_init'] = 15.15
data['step_size'] = 1.0
data['terrain_conditions']['max_channel_length'] = 3000.0

with open(TMP, 'w') as f:
    json.dump(data, f, indent=4)

subprocess.run([PYTHON, 'main_flowgo.py', TMP],
               cwd=os.path.dirname(os.path.abspath(__file__)))

out = 'results_flowgo/results_main_flowgo_La Palma 1971 Teneguia_15.15m3s.csv'
if os.path.exists(out):
    df = pd.read_csv(out)
    last = df.tail(3)[['position', 'viscosity', 'tho_0', 'core_temperature']]
    print("\nLast 3 rows:")
    print(last.to_string(index=False))
    print(f"\nFinal viscosity: {df['viscosity'].iloc[-1]:.3e} Pa·s  @ {df['position'].iloc[-1]:.0f} m")

if os.path.exists(TMP):
    os.remove(TMP)
