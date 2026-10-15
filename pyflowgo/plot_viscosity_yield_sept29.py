"""
Viscosity and Yield Strength plot for La Palma 2021 Sept 29 flow.
"""
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

csv_file = r"results_flowgo/results_main_flowgo_La Palma 2021 (Sept 29 Peak Phase)_100m3s.csv"

position, viscosity, yield_strength = [], [], []
with open(csv_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        pos = float(row['position'])
        visc = float(row['viscosity'])
        ys = float(row['tho_0'])
        if visc < 1e8 and ys < 1e6:
            position.append(pos)
            viscosity.append(visc)
            yield_strength.append(ys)

position = np.array(position)
viscosity = np.array(viscosity)
yield_strength = np.array(yield_strength)

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 11, 'axes.linewidth': 1.2, 'figure.dpi': 150,
})

fig, ax1 = plt.subplots(figsize=(14, 6))

color_visc = '#D32F2F'
ax1.semilogy(position, viscosity, color=color_visc, linewidth=2.2, label='Apparent Viscosity', zorder=3)
ax1.fill_between(position, viscosity, alpha=0.08, color=color_visc, zorder=1)
ax1.set_xlabel('Distance from Vent (m)', fontsize=13, fontweight='bold', labelpad=10)
ax1.set_ylabel('Apparent Viscosity (Pa·s)', fontsize=13, fontweight='bold', color=color_visc, labelpad=10)
ax1.tick_params(axis='y', labelcolor=color_visc, colors=color_visc, labelsize=11)

color_ys = '#1565C0'
ax2 = ax1.twinx()
ax2.plot(position, yield_strength, color=color_ys, linewidth=2.2,
         linestyle='--', dashes=(6, 3), label='Yield Strength', zorder=3)
ax2.fill_between(position, yield_strength, alpha=0.06, color=color_ys, zorder=1)
ax2.set_ylabel('Yield Strength (Pa)', fontsize=13, fontweight='bold', color=color_ys, labelpad=10)
ax2.tick_params(axis='y', labelcolor=color_ys, colors=color_ys, labelsize=11)

ax1.grid(True, alpha=0.25, linestyle='-', linewidth=0.5)
ax1.set_xlim(0, max(position) + 50)

ax1.axvline(x=max(position), color='#424242', linestyle=':', linewidth=1.2, alpha=0.6)
ax1.annotate('Flow front\n(near coast)', xy=(max(position), ax1.get_ylim()[1]),
             xytext=(-70, -20), textcoords='offset points',
             fontsize=10, fontstyle='italic', color='#424242', ha='right',
             arrowprops=dict(arrowstyle='->', color='#424242', lw=1.0))

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=11,
           framealpha=0.9, edgecolor='#bdbdbd', fancybox=True)

fig.suptitle('La Palma 2021 — Sept 29 Peak Phase (Vent to Coast)', fontsize=16, fontweight='bold', y=0.98)
ax1.set_title('Viscosity and Yield Strength Evolution with Distance', fontsize=12, color='#555555', pad=15)

info_text = (f'Effusion rate: 100 m³/s\n'
             f'Channel width: 12 m\n'
             f'T$_{{erupt}}$: 1413 K (1140°C)\n'
             f'Crystal frac.: 0.15\n'
             f'Crust cover: 0.90\n'
             f'Flow length: {max(position):.0f} m')
props = dict(boxstyle='round,pad=0.5', facecolor='#f5f5f5', edgecolor='#bdbdbd', alpha=0.9)
ax1.text(0.55, 0.35, info_text, transform=ax1.transAxes, fontsize=9.5,
         verticalalignment='top', horizontalalignment='right', bbox=props, family='monospace')

plt.tight_layout()
plt.savefig('results_flowgo/LaPalma2021_Sept29_viscosity_yield.png', dpi=200, bbox_inches='tight', facecolor='white')
print(f"Plot saved. Distance: {position[0]:.0f}–{position[-1]:.0f}m")
print(f"Viscosity: {viscosity.min():.0f}–{viscosity.max():.0f} Pa·s")
print(f"Yield strength: {yield_strength.min():.0f}–{yield_strength.max():.0f} Pa")
