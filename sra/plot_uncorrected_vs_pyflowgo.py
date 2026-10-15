import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# 1. Load MC Results
mc_df = pd.read_csv('LaPalmaData_MC_Results.csv')
sanjuan_mc = mc_df[mc_df['Eruption'] == '1949_San_Juan'].sort_values(by='Site', ascending=True)
teneguia_mc = mc_df[mc_df['Eruption'] == '1971_Teneguia'].sort_values(by='Site', ascending=True)

# 2. Load PyFlowGo Results
sanjuan_pfg = pd.read_csv(r'../pyflowgo/results_flowgo/results_flowgo_La Palma 1949 (San Juan - Llano del Banco)_22m3s.csv')
teneguia_pfg = pd.read_csv(r'../pyflowgo/results_flowgo/results_flowgo_La Palma 1971 Teneguia (res 10m)_19.3m3s.csv')

def plot_mc_vs_pfg(mc_data, pfg_data, eruption_name, color_mc, color_pfg, 
                   y_med_mc, y_col_pfg, ylabel, title, filename):
    
    fig, ax1 = plt.subplots(figsize=(10, 6))

    if mc_data.empty or pfg_data.empty:
        return
        
    # X-Axis 1: Elevation (for MC data)
    ax1.plot(mc_data['Elevation (m)'], mc_data[y_med_mc], 's--', label=f'{eruption_name} (SRA Uncorrected)', color=color_mc,
             markersize=6, markerfacecolor='white', markeredgewidth=1.5)
    
    ax1.set_xlabel('Elevation (m) [For SRA Data]', fontsize=12, color=color_mc)
    ax1.set_ylabel(ylabel, fontsize=12)
    ax1.tick_params(axis='x', labelcolor=color_mc)
    
    # Invert MC elevation axis so vent (high elevation) is on the left
    ax1.invert_xaxis()

    # X-Axis 2: Distance from Vent (for PyFlowGo data)
    ax2 = ax1.twiny()
    
    # Filter out wild PyFlowGo values at the very tip where it stops
    pfg_valid = pfg_data[(pfg_data['viscosity'] < 1e8) & (pfg_data['tho_0'] < 1e6)]
    
    ax2.plot(pfg_valid['position'], pfg_valid[y_col_pfg], '-', label=f'{eruption_name} (PyFlowGo)', color=color_pfg, linewidth=2.5)
    
    ax2.set_xlabel('Distance from Vent (m) [For PyFlowGo Data]', fontsize=12, color=color_pfg)
    ax2.tick_params(axis='x', labelcolor=color_pfg)

    # Set log scale
    ax1.set_yscale('log')
    ax1.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))
    ax1.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*0.1, numticks=15))
    ax1.yaxis.set_minor_formatter(ticker.NullFormatter())
    
    ax1.set_title(title, fontsize=14, pad=25)
    
    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='best', fontsize=11)
    
    ax1.grid(True, which="major", ls="-", alpha=0.5)
    ax1.grid(True, which="minor", ls=":", alpha=0.4)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

# Viscosity San Juan
plot_mc_vs_pfg(sanjuan_mc, sanjuan_pfg, '1949 San Juan', 
               color_mc='#e67e22', color_pfg='#2c3e50',
               y_med_mc='MC_Eta_Median_Unc', y_col_pfg='viscosity',
               ylabel='Apparent Viscosity $\eta$ (Pa·s)', 
               title='San Juan: Uncorrected SRA vs PyFlowGo Viscosity', 
               filename='viscosity_pfg_vs_uncorrected_sanjuan.png')

# Viscosity Teneguia
plot_mc_vs_pfg(teneguia_mc, teneguia_pfg, '1971 Teneguía', 
               color_mc='#e67e22', color_pfg='#2c3e50',
               y_med_mc='MC_Eta_Median_Unc', y_col_pfg='viscosity',
               ylabel='Apparent Viscosity $\eta$ (Pa·s)', 
               title='Teneguía: Uncorrected SRA vs PyFlowGo Viscosity', 
               filename='viscosity_pfg_vs_uncorrected_teneguia.png')

# Yield Strength San Juan
plot_mc_vs_pfg(sanjuan_mc, sanjuan_pfg, '1949 San Juan', 
               color_mc='#8e44ad', color_pfg='#16a085',
               y_med_mc='MC_Ty_Median_Unc', y_col_pfg='tho_0',
               ylabel='Yield Strength $\\tau_y$ (Pa)', 
               title='San Juan: Uncorrected SRA vs PyFlowGo Yield Strength', 
               filename='yield_strength_pfg_vs_uncorrected_sanjuan.png')

# Yield Strength Teneguia
plot_mc_vs_pfg(teneguia_mc, teneguia_pfg, '1971 Teneguía', 
               color_mc='#8e44ad', color_pfg='#16a085',
               y_med_mc='MC_Ty_Median_Unc', y_col_pfg='tho_0',
               ylabel='Yield Strength $\\tau_y$ (Pa)', 
               title='Teneguía: Uncorrected SRA vs PyFlowGo Yield Strength', 
               filename='yield_strength_pfg_vs_uncorrected_teneguia.png')

print("PyFlowGo vs Uncorrected SRA plots generated successfully.")
