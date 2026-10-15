"""
Run La Palma 1949 PyFLOWGO with ER=46 m³/s (calibrated config, 12° slope extension)
"""
import json, subprocess, os
import sys

TEMPLATE = "resource/LaPalma1949_calibrated.json"
PYTHON = sys.executable
TMP = "resource/_1949_6000m_run.json"

with open(TEMPLATE, 'r') as f:
    c = json.load(f)

c["effusion_rate_init"] = 46
c["step_size"] = 2.0
c["terrain_conditions"]["max_channel_length"] = 8000.0

with open(TMP, 'w') as f:
    json.dump(c, f, indent=2)

print("Running PyFLOWGO: ER=46 m³/s, step=3m, 12° slope extension, calibrated 1949 config...")
r = subprocess.run([PYTHON, "main_flowgo.py", TMP],
                   capture_output=False, text=True,
                   cwd=os.path.dirname(os.path.abspath(__file__)))
print(f"\nDone. Exit code: {r.returncode}")

if os.path.exists(TMP):
    os.remove(TMP)

