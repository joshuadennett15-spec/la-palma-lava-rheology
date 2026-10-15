import os
from pathlib import Path
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyflowgo.run_flowgo as run_flowgo
import glob

os.chdir(Path(__file__).resolve().parent)
path_to_folder = os.path.join(os.path.abspath(''), "results_flowgo")

resolutions = [2.0, 10.0, 30.0]
base_jsons = {
    '1949_San_Juan': 'resource/LaPalma1949_calibrated.json',
    '1971_Teneguia': 'resource/LaPalma1971_calibrated.json'
}

csv_outputs = {'1949_San_Juan': {}, '1971_Teneguia': {}}

# --- 1. RUN SIMULATIONS ---
for eruption, json_path in base_jsons.items():
    with open(json_path, "r") as f:
        base_data = json.load(f)
    
    orig_name = base_data['lava_name']
    
    for res in resolutions:
        # Create new config
        config = base_data.copy()
        config['step_size'] = res
        new_name = f"{orig_name} (res {int(res)}m)"
        config['lava_name'] = new_name
        
        tmp_json = f'resource/tmp_{eruption}_{int(res)}m.json'
        with open(tmp_json, 'w') as f:
            json.dump(config, f, indent=4)
        
        # Run PyFlowGo
        print(f"Running {eruption} at {res}m resolution...")
        flowgo = run_flowgo.RunFlowgo()
        flowgo.run(tmp_json, path_to_folder)
        
        # Get exact filename
        csv_file = flowgo.get_file_name_results(path_to_folder, tmp_json)
        csv_outputs[eruption][res] = csv_file
        
        # Cleanup temp JSON
        if os.path.exists(tmp_json):
            os.remove(tmp_json)

# --- 2. PREPARE PLOTTING ---
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 11,
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'figure.dpi': 150,
})

def calculate_sra_distances(df, eruption):
    sub = df[df['Eruption'] == eruption].sort_values('Site')
    x = sub['Easting (X)'].values
    y = sub['Northing (Y)'].values
    dists = [0]
    for i in range(1, len(x)):
        dists.append(dists[-1] + np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2))
    return np.array(dists), sub

def load_pf_data(csv_file, max_visc=1e8, max_ys=1e6):
    import csv
    pos_list, visc_list, ys_list = [], [], []
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

sra_df = pd.read_csv('../sra/LaPalmaData_MC_Results.csv')

def plot_resolution_comparison(eruption_name, title, color_base, out_file):
    sra_dists, sub = calculate_sra_distances(sra_df, eruption_name)
    
    sra_visc_med = sub['MC_Eta_Median'].values
    sra_visc_low = sub['MC_Eta_2.5%'].values
    sra_visc_high = sub['MC_Eta_97.5%'].values
    sra_ys_med = sub['MC_Ty_Median'].values
    sra_ys_low = sub['MC_Ty_2.5%'].values
    sra_ys_high = sub['MC_Ty_97.5%'].values

    # Determine scale factor to align SRA final point with 10m resolution PyFLOWGO max distance
    pos_pf_ref, _, _ = load_pf_data(csv_outputs[eruption_name][10.0])
    max_pf_dist = pos_pf_ref[-1]
    scale_factor = max_pf_dist / sra_dists[-1]
    sra_dists_scaled = sra_dists * scale_factor

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Colors for resolutions: 2m (dark), 10m (medium), 30m (light/dashed)
    res_styles = {
        2.0:  {'color': '#1a1a1a', 'ls': '-', 'lw': 2.5, 'label': '2m Resolution Mdl'},
        10.0: {'color': '#666666', 'ls': '--', 'lw': 2.5, 'label': '10m Resolution Mdl'},
        30.0: {'color': '#b3b3b3', 'ls': ':', 'lw': 2.5, 'label': '30m Resolution Mdl'}
    }

    # --- Top Panel: Viscosity ---
    for res in resolutions:
        csv_file = csv_outputs[eruption_name][res]
        pos_pf, visc_pf, ys_pf = load_pf_data(csv_file)
        ax1.plot(pos_pf, visc_pf, color=res_styles[res]['color'], 
                 linestyle=res_styles[res]['ls'], linewidth=res_styles[res]['lw'], 
                 label=res_styles[res]['label'], zorder=2)
                 
    yerr_visc = [sra_visc_med - sra_visc_low, sra_visc_high - sra_visc_med]
    ax1.errorbar(sra_dists_scaled, sra_visc_med, yerr=yerr_visc, fmt='o', color=color_base, 
                 label='Manual SRA Data', capsize=4, markersize=7, zorder=3,
                 markerfacecolor='white', markeredgewidth=1.5)

    ax1.set_yscale('log')
    ax1.set_ylabel(r'Apparent Viscosity ($\eta$, Pa·s)', fontsize=13, fontweight='bold')
    ax1.grid(True, which='both', linestyle='--', alpha=0.4)
    ax1.legend(loc='lower right', framealpha=0.9, fancybox=True)
    ax1.set_title(f'{title} - Resolution Uncertainty Comparison', fontsize=15, fontweight='bold', pad=15)
    ax1.set_xlim(left=0, right=max(sra_dists_scaled[-1], max_pf_dist) * 1.05)
    
    # Annotate sites (Viscosity)
    for i, txt in enumerate(sub['Site']):
        ax1.annotate('S'+str(txt), (sra_dists_scaled[i], sra_visc_med[i]), 
                     xytext=(0, 10), textcoords='offset points', ha='center',
                     fontsize=9, color=color_base, fontweight='bold')

    # --- Bottom Panel: Yield Strength ---
    for res in resolutions:
        csv_file = csv_outputs[eruption_name][res]
        pos_pf, visc_pf, ys_pf = load_pf_data(csv_file)
        ax2.plot(pos_pf, ys_pf, color=res_styles[res]['color'], 
                 linestyle=res_styles[res]['ls'], linewidth=res_styles[res]['lw'], 
                 label=res_styles[res]['label'], zorder=2)
                 
    yerr_ys = [sra_ys_med - sra_ys_low, sra_ys_high - sra_ys_med]
    ax2.errorbar(sra_dists_scaled, sra_ys_med, yerr=yerr_ys, fmt='s', color=color_base, 
                 label='Manual SRA Data', capsize=4, markersize=7, zorder=3,
                 markerfacecolor='white', markeredgewidth=1.5)

    ax2.set_yscale('log')
    ax2.set_xlabel('Distance from Vent (m)', fontsize=13, fontweight='bold')
    ax2.set_ylabel(r'Yield Strength ($\tau_y$, Pa)', fontsize=13, fontweight='bold')
    ax2.grid(True, which='both', linestyle='--', alpha=0.4)
    ax2.legend(loc='lower right', framealpha=0.9, fancybox=True)

    # Annotate sites (Yield Strength)
    for i, txt in enumerate(sub['Site']):
        ax2.annotate('S'+str(txt), (sra_dists_scaled[i], sra_ys_med[i]), 
                     xytext=(0, 10), textcoords='offset points', ha='center',
                     fontsize=9, color=color_base, fontweight='bold')

    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()

plot_resolution_comparison('1949_San_Juan', '1949 San Juan (Llano del Banco)', '#D32F2F', '../sra/resolution_comparison_1949.png')
plot_resolution_comparison('1971_Teneguia', '1971 Teneguía', '#1565C0', '../sra/resolution_comparison_1971.png')

print("Process completed successfully.")
