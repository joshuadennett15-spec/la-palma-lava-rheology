import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# Load the MC results
df = pd.read_csv('LaPalmaData_MC_Results.csv')

# Separate the dataset by eruption
sanjuan = df[df['Eruption'] == '1949_San_Juan'].copy()
teneguia = df[df['Eruption'] == '1971_Teneguia'].copy()

sanjuan = sanjuan.sort_values(by='Site', ascending=True)
teneguia = teneguia.sort_values(by='Site', ascending=True)

def create_compare_plot(data, eruption_name, color_corrected, color_uncorrected, 
                        y_col, 
                        y_med_corr, y_lower_corr, y_upper_corr, 
                        y_med_unc, y_lower_unc, y_upper_unc, 
                        ylabel, title, filename):
    fig, ax = plt.subplots(figsize=(10, 6))

    if data.empty:
        return
        
    # --- Uncorrected Data (Constant Q, No Stiffening) ---
    ax.fill_between(data['Site'], data[y_lower_unc], data[y_upper_unc], color=color_uncorrected, alpha=0.15)
    ax.plot(data['Site'], data[y_med_unc], 's--', label=f'{eruption_name} (Uncorrected)', color=color_uncorrected,
            markersize=5, markerfacecolor='white', markeredgewidth=1.2)

    # --- Corrected Data (Decay + Stiffening) ---
    ax.fill_between(data['Site'], data[y_lower_corr], data[y_upper_corr], color=color_corrected, alpha=0.2)
    ax.plot(data['Site'], data[y_med_corr], 'o-', label=f'{eruption_name} (Corrected)', color=color_corrected,
            markersize=6, markerfacecolor='white', markeredgewidth=1.5)

    ax.set_yscale('log')
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*0.1, numticks=15))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    
    ax.set_xlabel('Station Number (Site)', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, pad=15)
    
    ax.legend(fontsize=11)
    ax.grid(True, which="major", ls="-", alpha=0.5)
    ax.grid(True, which="minor", ls=":", alpha=0.4)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

# 1. Viscosity Comparison - San Juan
create_compare_plot(data=sanjuan, eruption_name='1949 San Juan', 
                    color_corrected='#E04C4C', color_uncorrected='#888888',
                    y_col='Viscosity', 
                    y_med_corr='MC_Eta_Median', y_lower_corr='MC_Eta_2.5%', y_upper_corr='MC_Eta_97.5%', 
                    y_med_unc='MC_Eta_Median_Unc', y_lower_unc='MC_Eta_2.5%_Unc', y_upper_unc='MC_Eta_97.5%_Unc',
                    ylabel='Apparent Viscosity $\eta$ (Pa·s)', 
                    title='1949 San Juan: Viscosity Comparison (Corrected vs Constant Q)', 
                    filename='viscosity_compare_sanjuan.png')

# 2. Viscosity Comparison - Teneguia
create_compare_plot(data=teneguia, eruption_name='1971 Teneguía', 
                    color_corrected='#4A80D6', color_uncorrected='#888888',
                    y_col='Viscosity', 
                    y_med_corr='MC_Eta_Median', y_lower_corr='MC_Eta_2.5%', y_upper_corr='MC_Eta_97.5%', 
                    y_med_unc='MC_Eta_Median_Unc', y_lower_unc='MC_Eta_2.5%_Unc', y_upper_unc='MC_Eta_97.5%_Unc',
                    ylabel='Apparent Viscosity $\eta$ (Pa·s)', 
                    title='1971 Teneguía: Viscosity Comparison (Corrected vs Constant Q)', 
                    filename='viscosity_compare_teneguia.png')

# 3. Yield Strength Comparison - San Juan
create_compare_plot(data=sanjuan, eruption_name='1949 San Juan', 
                    color_corrected='#E04C4C', color_uncorrected='#888888',
                    y_col='Yield Strength', 
                    y_med_corr='MC_Ty_Median', y_lower_corr='MC_Ty_2.5%', y_upper_corr='MC_Ty_97.5%', 
                    y_med_unc='MC_Ty_Median_Unc', y_lower_unc='MC_Ty_2.5%_Unc', y_upper_unc='MC_Ty_97.5%_Unc',
                    ylabel='Yield Strength $\\tau_y$ (Pa)', 
                    title='1949 San Juan: Yield Strength Comparison (Corrected vs Constant)', 
                    filename='yield_strength_compare_sanjuan.png')

# 4. Yield Strength Comparison - Teneguia
create_compare_plot(data=teneguia, eruption_name='1971 Teneguía', 
                    color_corrected='#4A80D6', color_uncorrected='#888888',
                    y_col='Yield Strength', 
                    y_med_corr='MC_Ty_Median', y_lower_corr='MC_Ty_2.5%', y_upper_corr='MC_Ty_97.5%', 
                    y_med_unc='MC_Ty_Median_Unc', y_lower_unc='MC_Ty_2.5%_Unc', y_upper_unc='MC_Ty_97.5%_Unc',
                    ylabel='Yield Strength $\\tau_y$ (Pa)', 
                    title='1971 Teneguía: Yield Strength Comparison (Corrected vs Constant)', 
                    filename='yield_strength_compare_teneguia.png')

print("Comparison plots generated.")
