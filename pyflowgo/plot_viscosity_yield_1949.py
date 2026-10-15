"""
Publication-quality plot: Viscosity and Yield Strength vs Distance from Vent
La Palma 1949 San Juan (Llano del Banco) — PyFlowGo simulation results
"""
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ── Load data ──────────────────────────────────────────────────────
csv_file = r"results_flowgo/results_main_flowgo_La Palma 1949 (San Juan - Llano del Banco)_46m3s.csv"

position = []
viscosity = []
yield_strength = []

with open(csv_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        pos = float(row['position'])
        visc = float(row['viscosity'])
        ys = float(row['tho_0'])
        if ys < 1e6:  # keep all viscosity values to show rheological cessation
            position.append(pos)
            viscosity.append(visc)
            yield_strength.append(ys)

position = np.array(position)
viscosity = np.array(viscosity)
yield_strength = np.array(yield_strength)

# ── Figure setup ───────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 11,
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'xtick.minor.width': 0.6,
    'ytick.minor.width': 0.6,
    'figure.dpi': 150,
})

fig, ax1 = plt.subplots(figsize=(12, 6))

# ── Viscosity (left axis, log scale) ──────────────────────────────
color_visc = '#D32F2F'
ax1.semilogy(position, viscosity, color=color_visc, linewidth=2.2,
             label='Apparent Viscosity', zorder=3)
ax1.fill_between(position, viscosity, alpha=0.08, color=color_visc, zorder=1)

ax1.set_xlabel('Distance from Vent (m)', fontsize=13, fontweight='bold', labelpad=10)
ax1.set_ylabel('Apparent Viscosity (Pa·s)', fontsize=13, fontweight='bold',
               color=color_visc, labelpad=10)
ax1.tick_params(axis='y', labelcolor=color_visc, colors=color_visc, labelsize=11)
ax1.tick_params(axis='x', labelsize=11)
ax1.yaxis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(2, 10)*0.1, numticks=20))
ax1.yaxis.set_minor_formatter(ticker.NullFormatter())

# ── Yield Strength (right axis, linear scale) ─────────────────────
color_ys = '#1565C0'
ax2 = ax1.twinx()
ax2.plot(position, yield_strength, color=color_ys, linewidth=2.2,
         linestyle='--', dashes=(6, 3), label='Yield Strength', zorder=3)
ax2.fill_between(position, yield_strength, alpha=0.06, color=color_ys, zorder=1)

ax2.set_ylabel('Yield Strength (Pa)', fontsize=13, fontweight='bold',
               color=color_ys, labelpad=10)
ax2.tick_params(axis='y', labelcolor=color_ys, colors=color_ys, labelsize=11)

# ── Grid and styling ──────────────────────────────────────────────
ax1.grid(True, alpha=0.25, linestyle='-', linewidth=0.5)
ax1.set_xlim(0, max(position) + 20)

# ── Flow front annotation ────────────────────────────────────────
ax1.axvline(x=max(position), color='#424242', linestyle=':', linewidth=1.2, alpha=0.6)
ax1.annotate('Flow front', xy=(max(position), ax1.get_ylim()[1]),
             xytext=(-50, -20), textcoords='offset points',
             fontsize=10, fontstyle='italic', color='#424242',
             ha='right',
             arrowprops=dict(arrowstyle='->', color='#424242', lw=1.0))

# ── Combined legend ───────────────────────────────────────────────
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
legend = ax1.legend(lines1 + lines2, labels1 + labels2,
                    loc='upper left', fontsize=11, framealpha=0.9,
                    edgecolor='#bdbdbd', fancybox=True)

# ── Title ─────────────────────────────────────────────────────────
fig.suptitle('La Palma 1949 San Juan (Llano del Banco)', fontsize=16, fontweight='bold', y=0.98)
ax1.set_title('Viscosity and Yield Strength Evolution with Distance',
              fontsize=12, color='#555555', pad=15)

# ── Info box ──────────────────────────────────────────────────────
info_text = (f'Effusion rate: 46 m³/s\n'
             f'Channel width: 5 m\n'
             f'T$_{{erupt}}$: 1413 K (1140°C)\n'
             f'Crystal frac.: 0.12\n'
             f'Crust cover: 0.9\n'
             f'Slope ext.: 12°\n'
             f'Flow length: {max(position):.0f} m')
props = dict(boxstyle='round,pad=0.5', facecolor='#f5f5f5',
             edgecolor='#bdbdbd', alpha=0.9)
ax1.text(0.98, 0.35, info_text, transform=ax1.transAxes,
         fontsize=9.5, verticalalignment='top', horizontalalignment='right',
         bbox=props, family='monospace')

plt.tight_layout()
plt.savefig('results_flowgo/LaPalma1949_viscosity_yield_strength.png',
            dpi=200, bbox_inches='tight', facecolor='white')
print(f"Plot saved to results_flowgo/LaPalma1949_viscosity_yield_strength.png")
print(f"Distance range: {position[0]:.0f} – {position[-1]:.0f} m")
print(f"Viscosity range: {viscosity.min():.0f} – {viscosity.max():.0f} Pa·s")
print(f"Yield strength range: {yield_strength.min():.0f} – {yield_strength.max():.0f} Pa")
