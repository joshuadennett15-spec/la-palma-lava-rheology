import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# 1. Calculate distances for SRA data
df = pd.read_csv('LaPalmaData_MC_Results.csv')
sra_dists = {}
sra_visc_med = {}
sra_visc_low = {}
sra_visc_high = {}
sra_ys_med = {}
sra_ys_low = {}
sra_ys_high = {}

for eruption in ['1949_San_Juan', '1971_Teneguia']:
    sub = df[df['Eruption'] == eruption].sort_values('Site')
    x = sub['Easting (X)'].values
    y = sub['Northing (Y)'].values
    dists = [0]
    for i in range(1, len(x)):
        dists.append(dists[-1] + np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2))
    
    sra_dists[eruption] = np.array(dists)
    sra_visc_med[eruption] = sub['MC_Eta_Median'].values
    sra_visc_low[eruption] = sub['MC_Eta_2.5%'].values
    sra_visc_high[eruption] = sub['MC_Eta_97.5%'].values
    sra_ys_med[eruption] = sub['MC_Ty_Median'].values
    sra_ys_low[eruption] = sub['MC_Ty_2.5%'].values
    sra_ys_high[eruption] = sub['MC_Ty_97.5%'].values

# 2. Load PyFLOWGO data
def load_pf_data(csv_file, max_visc, max_ys):
    pos_list, visc_list, ys_list = [], [], []
    import csv
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pos = float(row['position'])
            visc = float(row['viscosity'])
            ys = float(row['tho_0'])
            if visc < max_visc and ys < max_ys:
                pos_list.append(pos)
                visc_list.append(visc)
                ys_list.append(ys)
    return np.array(pos_list), np.array(visc_list), np.array(ys_list)

pf_1949_csv = '../pyflowgo/results_flowgo/results_main_flowgo_La Palma 1949 (San Juan - Llano del Banco)_22m3s.csv'
pf_1971_csv = '../pyflowgo/results_flowgo/results_main_flowgo_La Palma 1971 Teneguia_12.0m3s.csv'

pos_49, visc_49, ys_49 = load_pf_data(pf_1949_csv, 1e8, 1e6)
pos_71, visc_71, ys_71 = load_pf_data(pf_1971_csv, 1e8, 1e6)

# 3. Plot Combined Viscosity
fig, ax = plt.subplots(figsize=(12, 7))

# 1949
color_49 = '#D32F2F' # Red
ax.errorbar(sra_dists['1949_San_Juan'], sra_visc_med['1949_San_Juan'], 
            yerr=[sra_visc_med['1949_San_Juan'] - sra_visc_low['1949_San_Juan'], 
                  sra_visc_high['1949_San_Juan'] - sra_visc_med['1949_San_Juan']], 
            fmt='o', color=color_49, label='1949 San Juan (Manual SRA)', capsize=4, markersize=7)
ax.plot(pos_49, visc_49, '-', color=color_49, linewidth=2.5, label='1949 San Juan (PyFLOWGO)')

# 1971
color_71 = '#1565C0' # Blue
ax.errorbar(sra_dists['1971_Teneguia'], sra_visc_med['1971_Teneguia'], 
            yerr=[sra_visc_med['1971_Teneguia'] - sra_visc_low['1971_Teneguia'], 
                  sra_visc_high['1971_Teneguia'] - sra_visc_med['1971_Teneguia']], 
            fmt='s', color=color_71, label='1971 Teneguía (Manual SRA)', capsize=4, markersize=7)
ax.plot(pos_71, visc_71, '-', color=color_71, linewidth=2.5, label='1971 Teneguía (PyFLOWGO)')

ax.set_yscale('log')
ax.set_xlabel('Distance from Vent (m)', fontsize=14, fontweight='bold')
ax.set_ylabel(r'Apparent Viscosity $\eta$ (Pa·s)', fontsize=14, fontweight='bold')
ax.set_title('Apparent Viscosity: 1949 San Juan vs 1971 Teneguía', fontsize=16, pad=15, fontweight='bold')
ax.grid(True, which='both', linestyle='--', alpha=0.5)
ax.set_xlim(left=0)
ax.legend(fontsize=12, loc='upper left', framealpha=0.9, fancybox=True)

plt.tight_layout()
plt.savefig('combined_viscosity_1949_1971.png', dpi=300)
plt.close()

# 4. Plot Combined Yield Strength
fig, ax = plt.subplots(figsize=(12, 7))

# 1949
ax.errorbar(sra_dists['1949_San_Juan'], sra_ys_med['1949_San_Juan'], 
            yerr=[sra_ys_med['1949_San_Juan'] - sra_ys_low['1949_San_Juan'], 
                  sra_ys_high['1949_San_Juan'] - sra_ys_med['1949_San_Juan']], 
            fmt='o', color=color_49, label='1949 San Juan (Manual SRA)', capsize=4, markersize=7)
ax.plot(pos_49, ys_49, '-', color=color_49, linewidth=2.5, label='1949 San Juan (PyFLOWGO)')

# 1971
ax.errorbar(sra_dists['1971_Teneguia'], sra_ys_med['1971_Teneguia'], 
            yerr=[sra_ys_med['1971_Teneguia'] - sra_ys_low['1971_Teneguia'], 
                  sra_ys_high['1971_Teneguia'] - sra_ys_med['1971_Teneguia']], 
            fmt='s', color=color_71, label='1971 Teneguía (Manual SRA)', capsize=4, markersize=7)
ax.plot(pos_71, ys_71, '-', color=color_71, linewidth=2.5, label='1971 Teneguía (PyFLOWGO)')

ax.set_yscale('log')
ax.set_xlabel('Distance from Vent (m)', fontsize=14, fontweight='bold')
ax.set_ylabel(r'Yield Strength $\tau_y$ (Pa)', fontsize=14, fontweight='bold')
ax.set_title('Yield Strength: 1949 San Juan vs 1971 Teneguía', fontsize=16, pad=15, fontweight='bold')
ax.grid(True, which='both', linestyle='--', alpha=0.5)
ax.set_xlim(left=0)
ax.legend(fontsize=12, loc='upper left', framealpha=0.9, fancybox=True)

plt.tight_layout()
plt.savefig('combined_yield_strength_1949_1971.png', dpi=300)
plt.close()

print('Saved combined_viscosity_1949_1971.png and combined_yield_strength_1949_1971.png')
