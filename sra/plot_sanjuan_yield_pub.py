import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy import stats

# ------------------------------------------------------------------ #
# Data Loading
# ------------------------------------------------------------------ #
data_list = []
with open('LaPalmaData_MC_Results.csv', 'r') as f:
    for row in csv.DictReader(f):
        if row['Eruption'] == '1949_San_Juan':
            data_list.append(row)

data_list.sort(key=lambda x: int(x['Site']))

# Compute cumulative distance from vent using Easting/Northing
eastings = np.array([float(r['Easting (X)']) for r in data_list])
northings = np.array([float(r['Northing (Y)']) for r in data_list])
cum_dists = [0.0]
for i in range(1, len(eastings)):
    cum_dists.append(cum_dists[-1] + np.sqrt(
        (eastings[i] - eastings[i-1])**2 + (northings[i] - northings[i-1])**2))
cum_dists = np.array(cum_dists)

# Scale to vent range (1000–2500 m) as in resolution comparison
vent_range = (1000, 2500)
data_distances = vent_range[0] + (cum_dists / cum_dists[-1]) * (vent_range[1] - vent_range[0])
data_y_med = np.array([float(r['MC_Ty_Median']) for r in data_list])
data_y_lower = np.array([float(r['MC_Ty_2.5%']) for r in data_list])
data_y_upper = np.array([float(r['MC_Ty_97.5%']) for r in data_list])

# ------------------------------------------------------------------ #
# Global Plot Style
# ------------------------------------------------------------------ #
plt.rcParams.update({
    'font.family': 'serif',
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'figure.dpi': 300
})

# Create figure with single-column journal width
fig, ax = plt.subplots(figsize=(3.35, 2.8))

# ------------------------------------------------------------------ #
# Data & Line of Best Fit
# ------------------------------------------------------------------ #
x = data_distances
y = data_y_med

valid = y > 0
x_fit = x[valid]
y_log = np.log10(y[valid])

if len(x_fit) >= 3:
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_fit, y_log)
    r_sq = r_value ** 2

    x_smooth = np.linspace(x_fit.min(), x_fit.max(), 300)
    y_smooth = 10 ** (intercept + slope * x_smooth)

    # Calculate 95% Confidence Interval for the fit
    n = len(x_fit)
    y_pred = intercept + slope * x_fit
    s_err_fit = np.sqrt(np.sum((y_log - y_pred)**2) / (n - 2))
    t_val = stats.t.ppf(0.975, n - 2)
    mean_x = np.mean(x_fit)
    sum_sq_diff_x = np.sum((x_fit - mean_x)**2)
    
    y_line_log = intercept + slope * x_smooth
    confs = t_val * s_err_fit * np.sqrt(1/n + (x_smooth - mean_x)**2 / sum_sq_diff_x)
    
    lower_bound = 10 ** (y_line_log - confs)
    upper_bound = 10 ** (y_line_log + confs)

    # Plot 1: 95% CI Band
    ax.fill_between(x_smooth, lower_bound, upper_bound, 
                    color='#a63d40', alpha=0.12, label='95% CI', zorder=1)

    # Plot 2: Best-fit line
    ax.plot(x_smooth, y_smooth, '--', color='#a63d40', linewidth=1.5, 
            label='Log-linear fit', zorder=2)

    # Plot 3: Scatter points & error bars
    yerr = [data_y_med - data_y_lower, data_y_upper - data_y_med]
    ax.errorbar(data_distances, data_y_med, yerr=yerr, fmt='none', 
                ecolor='#2c3e6b', elinewidth=0.8, capsize=2, capthick=0.8, zorder=2.5)
                
    ax.scatter(data_distances, data_y_med, s=40, c='#2c3e6b', 
               edgecolors='white', linewidths=0.5, label='Data', zorder=3)

    # ------------------------------------------------------------------ #
    # Annotations & Formatting
    # ------------------------------------------------------------------ #
    # Axes limits (optional padding if needed)
    ax.set_ylim(bottom=data_y_lower.min() * 0.8)

    # Y-axis log formatting
    ax.set_yscale('log')
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=15))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())

    # Labels
    ax.set_xlabel('Distance from Vent (m)')
    ax.set_ylabel('Yield Strength, $\\tau_y$ (Pa)')

    # Legend
    ax.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.8, 
              edgecolor='lightgrey', fontsize=8)

    # Statistics text
    if p_value < 0.001:
        p_str = "$p$ < 0.001"
    else:
        p_str = f"$p$ = {p_value:.3f}"

    stats_text = (
        f"$R^2$ = {r_sq:.4f}\n"
        f"{p_str}\n"
        f"SE$_{{\\mathrm{{slope}}}}$ = {std_err:.4f}"
    )
    
    # Semi-transparent white bbox with no visible border
    bbox_props = dict(boxstyle='square,pad=0.3', facecolor='white', alpha=0.75, edgecolor='none')
    ax.text(0.95, 0.05, stats_text, transform=ax.transAxes, fontsize=7.5,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=bbox_props, zorder=4)

# Gridlines (Subtle horizontal only)
ax.grid(False, axis='x')
ax.grid(True, axis='y', which='major', color='#dddddd', linewidth=0.5, zorder=0)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Tight layout
plt.tight_layout(pad=0.5)
plt.savefig('yield_strength_plot_sanjuan_pub.png', dpi=300)
plt.savefig('yield_strength_plot_sanjuan_pub.pdf')
plt.close()

print("Publication quality plot generated: yield_strength_plot_sanjuan_pub.png")
