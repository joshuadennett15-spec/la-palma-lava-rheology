"""
Run La Palma 1971 Teneguia PyFLOWGO with 2m step size to capture viscosity blow-up at ~2490m.
"""
import json, subprocess, os
import sys

TEMPLATE = "resource/LaPalma1971_calibrated.json"
PYTHON = sys.executable
TMP = "resource/_1971_2m_run.json"

with open(TEMPLATE, 'r') as f:
    c = json.load(f)

c["step_size"] = 2.0
# Keep ER=13.8 (as calibrated) and max_channel_length=3000

with open(TMP, 'w') as f:
    json.dump(c, f, indent=2)

print(f"Running PyFLOWGO: ER={c['effusion_rate_init']} m³/s, step=2m, 1971 Teneguia config...")
r = subprocess.run([PYTHON, "main_flowgo.py", TMP],
                   capture_output=False, text=True,
                   cwd=os.path.dirname(os.path.abspath(__file__)))
print(f"\nDone. Exit code: {r.returncode}")

if os.path.exists(TMP):
    os.remove(TMP)
