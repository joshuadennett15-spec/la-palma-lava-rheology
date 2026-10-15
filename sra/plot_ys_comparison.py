"""
plot_ys_comparison.py
Compares PyFLOWGO yield strength (tho_0) output against SRA-derived
yield strength values for the 1949 San Juan and 1971 Teneguia eruptions.
Includes a log-linear LOBF through the SRA data, extended to the full
x-range so any intersection with the PyFLOWGO model is visible.
"""

import os
from pathlib import Path
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats

# ------------------------------------------------------------------ #
# Paths                                                               #
# ------------------------------------------------------------------ #
BASE   = Path(__file__).resolve().parent.parent
PF_DIR = os.path.join(BASE, 'pyflowgo', 'results_flowgo')
SRA_CSV = os.path.join(BASE, 'sra', 'LaPalmaData_MC_Results.csv')
OUT_DIR = os.path.join(BASE, 'sra')

# Best-calibrated 10m-resolution PyFLOWGO CSVs that reach the correct flow length
PF_FILES = {
    '1949_San_Juan': os.path.join(
        PF_DIR,
        'results_flowgo_La Palma 1949 (San Juan - Llano del Banco) (res 10m)_22.0m3s.csv'
    ),
    '1971_Teneguia': os.path.join(
        PF_DIR,
        'results_flowgo_La Palma 1971 Teneguia (res 10m)_12.0m3s.csv'
    ),
}

ERUPTION_META = {
    '1949_San_Juan': {'label': '1949 San Juan (Llano del Banco)', 'color': '#D32F2F'},
    '1971_Teneguia': {'label': '1971 Teneguia',                   'color': '#1565C0'},
}

# ------------------------------------------------------------------ #
# Helpers                                                             #
# ------------------------------------------------------------------ #
def load_pyflowgo(csv_path, max_ys=1e6):
    """Load position and tho_0 (yield strength) from a PyFLOWGO CSV."""
    pos, ys = [], []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            ys_val = float(row['tho_0'])
            if ys_val < max_ys:
                pos.append(float(row['position']))
                ys.append(ys_val)
    return np.array(pos), np.array(ys)


def calculate_sra_distances(sub):
    """Cumulative distance (m) between SRA stations using Easting/Northing."""
    x = sub['Easting (X)'].values
    y = sub['Northing (Y)'].values
    dists = [0.0]
    for i in range(1, len(x)):
        dists.append(dists[-1] + np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2))
    return np.array(dists)


# ------------------------------------------------------------------ #
# Global plot style                                                   #
# ------------------------------------------------------------------ #
plt.rcParams.update({
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':          11,
    'axes.linewidth':     1.2,
    'xtick.major.width':  1.0,
    'ytick.major.width':  1.0,
    'figure.dpi':         150,
})

# ------------------------------------------------------------------ #
# SRA data                                                            #
# ------------------------------------------------------------------ #
sra_df = pd.read_csv(SRA_CSV)

# ------------------------------------------------------------------ #
# One figure per eruption                                             #
# ------------------------------------------------------------------ #
for eruption, meta in ERUPTION_META.items():
    color  = meta['color']
    label  = meta['label']

    # --- SRA data ---
    sub = sra_df[sra_df['Eruption'] == eruption].sort_values('Site')
    sra_dists = calculate_sra_distances(sub)
    ys_med  = sub['MC_Ty_Median'].values
    ys_low  = sub['MC_Ty_2.5%'].values
    ys_high = sub['MC_Ty_97.5%'].values

    # --- PyFLOWGO data ---
    pf_path = PF_FILES[eruption]
    pos_pf, ys_pf = load_pyflowgo(pf_path)

    # Scale SRA distances so the last station aligns with PyFLOWGO's max distance
    scale = pos_pf[-1] / sra_dists[-1] if sra_dists[-1] > 0 else 1.0
    sra_dists_sc = sra_dists * scale

    # ----------------------------------------------------------------
    # Log-linear LOBF through SRA medians
    # Fit: log10(tau_y) = slope * x + intercept  (in log space)
    # Extended across the full plot x-range so intersection is visible
    # ----------------------------------------------------------------
    valid = ys_med > 0
    x_fit = sra_dists_sc[valid]
    y_log = np.log10(ys_med[valid])

    slope, intercept, r_value, p_value, std_err = stats.linregress(x_fit, y_log)
    r_sq = r_value ** 2

    x_max_plot = max(sra_dists_sc[-1], pos_pf[-1]) * 1.05
    x_lobf = np.linspace(0, x_max_plot, 500)
    y_lobf = 10 ** (intercept + slope * x_lobf)

    p_str = "p < 0.001" if p_value < 0.001 else f"p = {p_value:.3f}"
    sign  = '+' if intercept >= 0 else '-'
    eq_str = f"log10(ty) = {slope:.4f}x {sign} {abs(intercept):.2f}"

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(10, 6))

    # PyFLOWGO model line
    ax.plot(pos_pf, ys_pf, color='#333333', linewidth=2.5, linestyle='-',
            label='PyFLOWGO Model (10m res)', zorder=2)

    # SRA error bars
    yerr = [ys_med - ys_low, ys_high - ys_med]
    ax.errorbar(sra_dists_sc, ys_med, yerr=yerr,
                fmt='s', color=color, capsize=4, markersize=7, zorder=3,
                markerfacecolor='white', markeredgewidth=1.8,
                label='SRA Estimate (95% CI)')

    # LOBF — extended across full plot width
    ax.plot(x_lobf, y_lobf, linestyle='--', linewidth=1.8,
            color=color, alpha=0.75,
            label=f'SRA LOBF: {eq_str}', zorder=2)

    # Stats annotation
    stats_text = (f"R2 = {r_sq:.4f}\n"
                  f"{p_str}\n"
                  f"SE_slope = {std_err:.5f}")
    ax.text(0.02, 0.97, stats_text,
            transform=ax.transAxes, fontsize=9.5,
            verticalalignment='top', horizontalalignment='left',
            family='monospace',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor=color, linewidth=1.1, alpha=0.9))

    # Site labels
    for i, site in enumerate(sub['Site']):
        ax.annotate(f'S{int(site)}',
                    xy=(sra_dists_sc[i], ys_med[i]),
                    xytext=(0, 10), textcoords='offset points',
                    ha='center', fontsize=9, color=color, fontweight='bold')

    # Axes
    ax.set_yscale('log')
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10)*0.1, numticks=15))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())

    ax.set_xlim(left=0, right=x_max_plot)
    ax.set_xlabel('Distance from Vent (m)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Yield Strength (tau_y, Pa)', fontsize=13, fontweight='bold')
    ax.set_title(f'{label} - PyFLOWGO vs SRA Yield Strength (with LOBF)',
                 fontsize=14, fontweight='bold', pad=15)

    ax.grid(True, which='both', linestyle='--', alpha=0.4)
    ax.legend(loc='lower right', framealpha=0.9, fancybox=True, fontsize=10)

    plt.tight_layout()
    out_name = f'ys_comparison_{eruption}.png'
    out_path = os.path.join(OUT_DIR, out_name)
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved: {out_path}")

print("Done.")
