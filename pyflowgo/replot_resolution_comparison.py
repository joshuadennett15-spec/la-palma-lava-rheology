"""
replot_resolution_comparison.py
Rebuilds the resolution comparison plots using existing PyFLOWGO CSVs
and adds log-linear LOBFs through the SRA data in both panels.
Does NOT re-run PyFLOWGO simulations.
"""
import os
from pathlib import Path
import csv as csv_mod
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats

os.chdir(Path(__file__).resolve().parent)

# ------------------------------------------------------------------ #
# Config                                                              #
# ------------------------------------------------------------------ #
resolutions = ['new_13.8', 'new_46']

csv_outputs = {
    '1949_San_Juan': {
        'new_46': 'results_flowgo/results_main_flowgo_La Palma 1949 (San Juan - Llano del Banco)_46m3s.csv',
    },
    '1971_Teneguia': {
        'new_13.8': 'results_flowgo/results_main_flowgo_La Palma 1971 Teneguia_13.8m3s.csv',
    },
}

lobf_config = {
    '1949_San_Juan': {'vent_range': (1000, 2500)},
    '1971_Teneguia': {'vent_offset': 50.0},
}

# ------------------------------------------------------------------ #
# Plot style                                                          #
# ------------------------------------------------------------------ #
plt.rcParams.update({
    'font.family': 'serif',
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'figure.dpi': 300,
})

# ------------------------------------------------------------------ #
# Helpers                                                             #
# ------------------------------------------------------------------ #
def calculate_sra_distances(data_list, eruption):
    sub = [r for r in data_list if r['Eruption'] == eruption]
    sub.sort(key=lambda r: int(r['Site']))
    x = np.array([float(r['Easting (X)']) for r in sub])
    y = np.array([float(r['Northing (Y)']) for r in sub])
    dists = [0]
    for i in range(1, len(x)):
        dists.append(dists[-1] + np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2))
    return np.array(dists), sub


def load_pf_data(csv_file, max_ys=1e6):
    pos_list, visc_list, ys_list = [], [], []
    with open(csv_file, 'r') as f:
        reader = csv_mod.DictReader(f)
        for row in reader:
            pos  = float(row['position'])
            visc = float(row['viscosity'])
            ys   = float(row['tho_0'])
            if ys < max_ys: # remove visc ceiling to see blow-up
                pos_list.append(pos)
                visc_list.append(visc)
                ys_list.append(ys)
    return np.array(pos_list), np.array(visc_list), np.array(ys_list)


def add_lobf(ax, x_data, y_data, color, x_max, measured_range=None):
    valid = y_data > 0
    if valid.sum() < 3:
        return None
    xv = x_data[valid]
    yv = np.log10(y_data[valid])
    slope, intercept, r_value, p_value, std_err = stats.linregress(xv, yv)
    r_sq = r_value ** 2

    n = len(xv)
    y_pred = intercept + slope * xv
    s_err = np.sqrt(np.sum((yv - y_pred)**2) / (n - 2))
    t_val = stats.t.ppf(0.975, n - 2)
    mean_x = np.mean(xv)
    sum_sq_diff_x = np.sum((xv - mean_x)**2)

    x_line = np.linspace(0, x_max, 1000)
    y_line_log = intercept + slope * x_line
    confs = t_val * s_err * np.sqrt(1/n + (x_line - mean_x)**2 / sum_sq_diff_x)

    y_line = 10 ** y_line_log
    lower_bound = 10 ** (y_line_log - confs)
    upper_bound = 10 ** (y_line_log + confs)

    if measured_range is None:
        ax.plot(x_line, y_line, linestyle='--', linewidth=1.8, color=color, alpha=0.75,
                label='SRA trend', zorder=3)
        ax.fill_between(x_line, lower_bound, upper_bound, color=color, alpha=0.15,
                        label='95% CI trend', zorder=1)
    else:
        x0, x1 = measured_range
        mask_fit = (x_line >= x0) & (x_line <= x1)
        mask_ext = (x_line < x0) | (x_line > x1)
        
        ax.plot(x_line[mask_fit], y_line[mask_fit], linestyle='-', linewidth=2.0, color=color, alpha=0.9,
                label='SRA trend', zorder=3)
        ax.plot(x_line[mask_ext], y_line[mask_ext], linestyle=':', linewidth=1.6, color=color, alpha=0.5,
                label='trend extrapolation', zorder=2)
        ax.fill_between(x_line[mask_fit], lower_bound[mask_fit], upper_bound[mask_fit], color=color, alpha=0.15,
                        label='95% CI trend', zorder=1)

    p_str = "$p < 0.001$" if p_value < 0.001 else f"$p$ = {p_value:.3f}"
    return f"$R^2$ = {r_sq:.4f}\n{p_str}\nSE$_{{\\mathrm{{slope}}}}$ = {std_err:.4f}"


sra_data = []
with open('../sra/LaPalmaData_MC_Results.csv', 'r') as f:
    for row in csv_mod.DictReader(f):
        sra_data.append(row)

# ------------------------------------------------------------------ #
# Plotting                                                            #
# ------------------------------------------------------------------ #
def plot_resolution_comparison(eruption_name, title, color_base, out_file, lobf_cfg=None):
    if lobf_cfg is None:
        lobf_cfg = {}
    sra_dists, sub = calculate_sra_distances(sra_data, eruption_name)

    sra_visc_med  = np.array([float(r['MC_Eta_Median']) for r in sub])
    sra_visc_low  = np.array([float(r['MC_Eta_2.5%']) for r in sub])
    sra_visc_high = np.array([float(r['MC_Eta_97.5%']) for r in sub])
    sra_ys_med    = np.array([float(r['MC_Ty_Median']) for r in sub])
    sra_ys_low    = np.array([float(r['MC_Ty_2.5%']) for r in sub])
    sra_ys_high   = np.array([float(r['MC_Ty_97.5%']) for r in sub])

    max_pf_dist = 0
    for res in csv_outputs[eruption_name]:
        pos_pf_ref, _, _ = load_pf_data(csv_outputs[eruption_name][res])
        if len(pos_pf_ref) > 0:
            max_pf_dist = max(max_pf_dist, pos_pf_ref[-1])

    vent_range  = lobf_cfg.get('vent_range')
    vent_offset = lobf_cfg.get('vent_offset')

    if vent_range is not None:
        r0, r1 = vent_range
        sra_dists_scaled = r0 + (sra_dists / sra_dists[-1]) * (r1 - r0)
        measured_range   = (float(r0), float(r1))
    elif vent_offset is not None:
        sra_dists_scaled = sra_dists + vent_offset
        measured_range   = (float(sra_dists_scaled[0]), float(sra_dists_scaled[-1]))
    else:
        scale_factor     = max_pf_dist / sra_dists[-1]
        sra_dists_scaled = sra_dists * scale_factor
        measured_range   = (float(sra_dists_scaled[0]), float(sra_dists_scaled[-1]))

    if eruption_name == '1971_Teneguia':
        x_max_plot = 2600
    else:
        x_max_plot = max(sra_dists_scaled[-1], max_pf_dist) * 1.05

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9), sharex=True)

    res_styles = {
        'new_13.8': {'color': '#555555', 'ls': '-', 'lw': 1.2, 'label': 'PyFLOWGO'},
        'new_46': {'color': '#555555', 'ls': '-', 'lw': 1.2, 'label': 'PyFLOWGO'},
    }

    # Add panel labels
    ax1.text(-0.08, 1.02, '(a)', transform=ax1.transAxes, size=14, weight='bold')
    ax2.text(-0.08, 1.02, '(b)', transform=ax2.transAxes, size=14, weight='bold')

    # ---- Top panel: Viscosity ----
    for res in resolutions:
        if res not in csv_outputs[eruption_name]: continue
        pos_pf, visc_pf, _ = load_pf_data(csv_outputs[eruption_name][res])
        ax1.plot(pos_pf, visc_pf,
                 color=res_styles[res]['color'],
                 linestyle=res_styles[res]['ls'],
                 linewidth=res_styles[res]['lw'],
                 label=res_styles[res]['label'], zorder=2)

    yerr_visc = [sra_visc_med - sra_visc_low, sra_visc_high - sra_visc_med]
    ax1.errorbar(sra_dists_scaled, sra_visc_med, yerr=yerr_visc,
                 fmt='none', ecolor=color_base, elinewidth=1.0, capsize=3, capthick=1.0, zorder=2.5)
    
    ax1.scatter(sra_dists_scaled, sra_visc_med, marker='s', s=45, color='none',
                edgecolors=color_base, linewidths=1.8, label='SRA field data', zorder=3)

    visc_stats = add_lobf(ax1, sra_dists_scaled, sra_visc_med, color_base,
                          x_max_plot, measured_range=measured_range)

    ax1.set_yscale('log')
    ax1.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))
    ax1.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10)*0.1, numticks=15))
    ax1.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax1.set_ylabel(r'Apparent Viscosity ($\eta_A$, Pa·s)', fontsize=13)
    ax1.grid(True, which='major', color='lightgrey', alpha=0.3, linestyle='-')
    # Removed title for publication styling
    # ax1.set_title(f'{title}{title_suffix}', fontsize=15, pad=15)
    
    # Clip the y-axis to prevent the PyFLOWGO terminal blow-up from crushing the field data
    # We set the top to ~5 times the highest field data median to give it breathing room 
    # without showing the full 10^7 spike that wastes 50% of the plot area
    ax1.set_ylim(top=np.max(sra_visc_med) * 5)
    ax1.set_xlim(left=0, right=x_max_plot)
    
    if visc_stats:
        bbox_props = dict(boxstyle='square,pad=0.3', facecolor='white', alpha=0.75, edgecolor='none')
        ax1.text(0.95, 0.05, visc_stats, transform=ax1.transAxes, fontsize=8,
                 verticalalignment='bottom', horizontalalignment='right', 
                 bbox=bbox_props, zorder=4)

    # ---- Bottom panel: Yield Strength ----
    for res in resolutions:
        if res not in csv_outputs[eruption_name]: continue
        pos_pf, _, ys_pf = load_pf_data(csv_outputs[eruption_name][res])
        ax2.plot(pos_pf, ys_pf,
                 color=res_styles[res]['color'],
                 linestyle=res_styles[res]['ls'],
                 linewidth=res_styles[res]['lw'],
                 label=res_styles[res]['label'], zorder=2)

    yerr_ys = [sra_ys_med - sra_ys_low, sra_ys_high - sra_ys_med]
    ax2.errorbar(sra_dists_scaled, sra_ys_med, yerr=yerr_ys,
                 fmt='none', ecolor=color_base, elinewidth=1.0, capsize=3, capthick=1.0, zorder=2.5)
    
    ax2.scatter(sra_dists_scaled, sra_ys_med, marker='s', s=45, color='none',
                edgecolors=color_base, linewidths=1.8, label='SRA field data', zorder=3)

    ys_stats = add_lobf(ax2, sra_dists_scaled, sra_ys_med, color_base,
                        x_max_plot, measured_range=measured_range)

    ax2.set_yscale('log')
    ax2.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))
    ax2.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10)*0.1, numticks=15))
    ax2.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax2.set_xlabel('Distance from Vent (m)', fontsize=13)
    ax2.set_ylabel(r'Yield Strength ($\tau_y$, Pa)', fontsize=13)
    # ax2.grid(True, which='major', color='lightgrey', alpha=0.3, linestyle='-')
    
    # Clip the y-axis for the same reason
    ax2.set_ylim(top=np.max(sra_ys_med) * 3)
    
    if ys_stats:
        bbox_props = dict(boxstyle='square,pad=0.3', facecolor='white', alpha=0.75, edgecolor='none')
        ax2.text(0.95, 0.05, ys_stats, transform=ax2.transAxes, fontsize=8,
                 verticalalignment='bottom', horizontalalignment='right', 
                 bbox=bbox_props, zorder=4)

    hands, labs = [], []
    for ax in [ax1, ax2]:
        h, l = ax.get_legend_handles_labels()
        hands.extend(h)
        labs.extend(l)
    by_label = dict(zip(labs, hands))
    
    order_labels = ['2 m step (13.8 m³/s)', '2 m step (46.0 m³/s)', 'SRA field data', 'SRA trend', '95% CI trend', 'trend extrapolation']
    plot_handles = [by_label[lbl] for lbl in order_labels if lbl in by_label]
    plot_labels = [lbl for lbl in order_labels if lbl in by_label]
    
    # --- Spines, grids, and tight layout ---
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(False, axis='x')
        # We replace the y grid with just horizontal lines behind
        ax.grid(True, axis='y', which='major', color='#dddddd', linewidth=0.5, zorder=0)

    fig.subplots_adjust(top=0.92, bottom=0.15)
    # Legend style updates
    fig.legend(plot_handles, plot_labels, loc='upper center', ncol=3, 
               frameon=True, facecolor='white', framealpha=0.8, edgecolor='lightgrey', 
               fontsize=8, bbox_to_anchor=(0.5, 0.08))

    plt.tight_layout(pad=0.5)
    plt.savefig(out_file, dpi=300)
    plt.savefig(out_file.replace('.png', '.pdf'))
    plt.close()
    print(f"Saved: {out_file} (and .pdf)")


plot_resolution_comparison('1949_San_Juan', '1949 San Juan (Llano del Banco)', '#D32F2F',
                           '../sra/resolution_comparison_1949.png',
                           lobf_cfg=lobf_config['1949_San_Juan'])
plot_resolution_comparison('1971_Teneguia', '1971 Teneguia', '#1565C0',
                           '../sra/resolution_comparison_1971.png',
                           lobf_cfg=lobf_config['1971_Teneguia'])

print("Done.")
