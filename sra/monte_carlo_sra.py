import numpy as np
import pandas as pd
import sys

# Parameters
N = 100000 # Number of MC trials
g = 9.81
rho_mean = 2407.5 # Midpoint of 2403 to 2412
rho_std = 175.5   # Midpoint of 170 to 181

# Downstream Loss Factor (exponential decay constant)
# k = 0.025 overrides the H^3 thinning effect downstream.
k_decay = 0.025


# File paths
input_file = 'sheet2_out.csv'
output_file = 'LaPalmaData_MC_Results.csv'

try:
    df = pd.read_csv(input_file)
except FileNotFoundError:
    print(f"Error: {input_file} not found.")
    sys.exit(1)

# lists to store results
results = []

print(f"Running Monte Carlo simulation with N={N} trials per site...")
print(f"Applying Elevation-Dependent Effusion Rate decay (k={k_decay})...")
# Group by eruption to find the max elevation (Vent elevation)
vent_elevations = df.groupby('Eruption')['Elevation (m)'].max().to_dict()

for idx, row in df.iterrows():
    H_mean, H_sem = row['Height_Mean'], row['Height_SEM']
    W_mean, W_sem = row['Width_Mean'], row['Width_SEM']
    S_mean, S_sem = row['Slope_Mean'], row['Slope_SEM']
    
    E_site = row['Elevation (m)']
    E_vent = vent_elevations[row['Eruption']]
    
    # Calculate the local Effusion Rate mean using exponential decay from the Vent
    Q_vent_mean = row['Q (m³/s)'] # Assuming the CSV column holds the vent Q
    Q_local_mean = Q_vent_mean * np.exp(-k_decay * (E_vent - E_site))
    
    
    # Draw random samples
    H_samp = np.clip(np.random.normal(H_mean, H_sem, size=N), a_min=0.01, a_max=None)
    W_samp = np.clip(np.random.normal(W_mean, W_sem, size=N), a_min=0.01, a_max=None)
    S_samp_deg = np.random.normal(S_mean, S_sem, size=N)
    S_samp_rad = np.radians(S_samp_deg)
    
    rho_samp = np.random.normal(rho_mean, rho_std, size=N)
    
    # Sample Q around the decayed local mean (keeping 10% relative uncertainty)
    Q_samp = np.clip(np.random.normal(Q_local_mean, 0.10 * Q_local_mean, size=N), a_min=0.01, a_max=None)
    
    # Calculate Yield Strength
    ty_samp = rho_samp * g * H_samp * np.sin(S_samp_rad)
    
    # Calculate Shape Factor
    n_samp = 3.0 * (1.0 + H_samp / W_samp)**2
    
    # Calculate Apparent Viscosity
    eta_samp = (rho_samp * g * H_samp**2 * W_samp * np.sin(S_samp_rad)) / (n_samp * Q_samp)
    
    # --- Uncorrected Calculations ---
    Q_local_mean_uncorrected = Q_vent_mean
    
    Q_samp_uncorrected = np.clip(np.random.normal(Q_local_mean_uncorrected, 0.10 * Q_local_mean_uncorrected, size=N), a_min=0.01, a_max=None)
    ty_samp_uncorrected = rho_samp * g * H_samp * np.sin(S_samp_rad)
    eta_samp_uncorrected = (rho_samp * g * H_samp**2 * W_samp * np.sin(S_samp_rad)) / (n_samp * Q_samp_uncorrected)
    
    # Extract statistics
    ty_med = np.median(ty_samp)
    ty_std = np.std(ty_samp)
    ty_2_5 = np.percentile(ty_samp, 2.5)
    ty_97_5 = np.percentile(ty_samp, 97.5)
    
    eta_med = np.median(eta_samp)
    eta_std = np.std(eta_samp)
    eta_2_5 = np.percentile(eta_samp, 2.5)
    eta_97_5 = np.percentile(eta_samp, 97.5)
    
    # Uncorrected statistics
    ty_med_unc = np.median(ty_samp_uncorrected)
    ty_std_unc = np.std(ty_samp_uncorrected)
    ty_2_5_unc = np.percentile(ty_samp_uncorrected, 2.5)
    ty_97_5_unc = np.percentile(ty_samp_uncorrected, 97.5)
    
    eta_med_unc = np.median(eta_samp_uncorrected)
    eta_std_unc = np.std(eta_samp_uncorrected)
    eta_2_5_unc = np.percentile(eta_samp_uncorrected, 2.5)
    eta_97_5_unc = np.percentile(eta_samp_uncorrected, 97.5)
    
    # Store
    results.append({
        'Site': row['Site'],
        'Eruption': row['Eruption'],
        'Q_Local_Mean': Q_local_mean, # Store this to see the decay
        'MC_Ty_Median': ty_med,
        'MC_Ty_Std': ty_std,
        'MC_Ty_2.5%': ty_2_5,
        'MC_Ty_97.5%': ty_97_5,
        'MC_Eta_Median': eta_med,
        'MC_Eta_Std': eta_std,
        'MC_Eta_2.5%': eta_2_5,
        'MC_Eta_97.5%': eta_97_5,
        'MC_Ty_Median_Unc': ty_med_unc,
        'MC_Ty_Std_Unc': ty_std_unc,
        'MC_Ty_2.5%_Unc': ty_2_5_unc,
        'MC_Ty_97.5%_Unc': ty_97_5_unc,
        'MC_Eta_Median_Unc': eta_med_unc,
        'MC_Eta_Std_Unc': eta_std_unc,
        'MC_Eta_2.5%_Unc': eta_2_5_unc,
        'MC_Eta_97.5%_Unc': eta_97_5_unc
    })

results_df = pd.DataFrame(results)
final_df = pd.merge(df, results_df, on=['Eruption', 'Site'])

final_df.to_csv(output_file, index=False)
print(f"Simulation complete. Results saved to {output_file}")
