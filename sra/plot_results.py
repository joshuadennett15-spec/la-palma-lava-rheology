import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy import stats
import seaborn as sns

# Set publication quality styling
sns.set_theme(style="ticks", context="paper")
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "axes.linewidth": 1.0,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "axes.edgecolor": "black"
})

# Load the MC results
df = pd.read_csv('LaPalmaData_MC_Results.csv')

# Separate the dataset by eruption
sanjuan = df[df['Eruption'] == '1949_San_Juan'].copy()
teneguia = df[df['Eruption'] == '1971_Teneguia'].copy()

# Sort by Site ascending, so the plot flows sequentially from Vent (Site 1) to Front
sanjuan = sanjuan.sort_values(by='Site', ascending=True)
teneguia = teneguia.sort_values(by='Site', ascending=True)


def create_plot(data, eruption_name, color, y_col, y_med, y_lower, y_upper, ylabel, title, filename):
    fig, ax = plt.subplots(figsize=(10, 6))

    if data.empty:
        return

    # --- 95% uncertainty bars and Median data points ---
    y_err_lower = data[y_med] - data[y_lower]
    y_err_upper = data[y_upper] - data[y_med]
    
    ax.errorbar(data['Site'], data[y_med], yerr=[y_err_lower, y_err_upper],
                fmt='o', label=eruption_name, color=color, ecolor=color,
                elinewidth=1.5, capsize=4, markersize=8, markerfacecolor=color,
                markeredgecolor='black', markeredgewidth=1.0)

    # ---------------------------------------------------------------
    # Line of Best Fit — Exponential regression (log-linear)
    #
    # Because the y-axis is log-scale, we fit:
    #   log10(y) = intercept + slope * Site
    # which is equivalent to: y = 10^intercept * 10^(slope * Site)
    #
    # This gives a straight line on the semi-log plot and is the most
    # statistically meaningful fit for monotonically growing rheology.
    # Stats reported: R, R², p-value (two-tailed), SE of slope.
    # ---------------------------------------------------------------
    x = data['Site'].values
    y = data[y_med].values

    valid = y > 0
    x_fit = x[valid]
    y_log = np.log10(y[valid])

    if len(x_fit) >= 3:
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_fit, y_log)
        r_sq = r_value ** 2

        # Smooth best-fit curve
        x_smooth = np.linspace(x_fit.min(), x_fit.max(), 300)
        y_smooth = 10 ** (intercept + slope * x_smooth)

        # Equation label: y = 10^a * 10^(bx)
        eq_label = (
            f"Best fit: $\\log_{{10}}(y) = {slope:.3f}x {'+' if intercept >= 0 else '-'} {abs(intercept):.2f}$"
        )
        ax.plot(x_smooth, y_smooth, '--', color=color, linewidth=1.8, alpha=0.85, label=eq_label)

        # Format p-value
        if p_value < 0.001:
            p_str = r"$p < 0.001$"
        elif p_value < 0.01:
            p_str = f"$p = {p_value:.4f}$"
        else:
            p_str = f"$p = {p_value:.3f}$"

        # Statistics annotation box (bottom-right corner)
        stats_text = (
            f"$R = {r_value:+.4f}$\n"
            f"$R^2 = {r_sq:.4f}$\n"
            f"{p_str}\n"
            f"$\\mathrm{{SE}}_{{\\mathrm{{slope}}}} = {std_err:.4f}$"
        )
        ax.text(0.97, 0.05, stats_text,
                transform=ax.transAxes,
                fontsize=11,
                verticalalignment='bottom',
                horizontalalignment='right',
                family='sans-serif',
                bbox=dict(boxstyle='round,pad=0.45', facecolor='#FFFFE0', # pale yellow to match reference
                          edgecolor='black', linewidth=1.0, alpha=0.9))

    # --- Axes formatting ---
    ax.set_yscale('log')
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=15))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())

    ax.set_xlabel('Station Number (Site)', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, pad=15, fontweight='bold')

    ax.legend(fontsize=11, frameon=True, edgecolor='black', facecolor='white')
    
    # Light gray grid lines
    ax.grid(True, which="major", ls="-", alpha=0.4, color="#E0E0E0")
    ax.grid(True, which="minor", ls=":", alpha=0.2, color="#E0E0E0")
    
    # Ensure all box spines are black and visible
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1.0)
        
    ax.set_axisbelow(True) # Put grid lines behind plot elements

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()


# 1. Plot Apparent Viscosity for San Juan
create_plot(data=sanjuan, eruption_name='1949 San Juan', color='#E04C4C',
            y_col='Viscosity', y_med='MC_Eta_Median', y_lower='MC_Eta_2.5%', y_upper='MC_Eta_97.5%',
            ylabel='Apparent Viscosity $\\eta$ (Pa·s)',
            title='1949 San Juan: Apparent Viscosity vs Station Number',
            filename='viscosity_plot_sanjuan.png')

# 2. Plot Apparent Viscosity for Teneguia
create_plot(data=teneguia, eruption_name='1971 Teneguía', color='#4A80D6',
            y_col='Viscosity', y_med='MC_Eta_Median', y_lower='MC_Eta_2.5%', y_upper='MC_Eta_97.5%',
            ylabel='Apparent Viscosity $\\eta$ (Pa·s)',
            title='1971 Teneguía: Apparent Viscosity vs Station Number',
            filename='viscosity_plot_teneguia.png')

# 3. Plot Yield Strength for San Juan
create_plot(data=sanjuan, eruption_name='1949 San Juan', color='#E04C4C',
            y_col='Yield Strength', y_med='MC_Ty_Median', y_lower='MC_Ty_2.5%', y_upper='MC_Ty_97.5%',
            ylabel='Yield Strength $\\tau_y$ (Pa)',
            title='1949 San Juan: Yield Strength vs Station Number',
            filename='yield_strength_plot_sanjuan.png')

# 4. Plot Yield Strength for Teneguia
create_plot(data=teneguia, eruption_name='1971 Teneguía', color='#4A80D6',
            y_col='Yield Strength', y_med='MC_Ty_Median', y_lower='MC_Ty_2.5%', y_upper='MC_Ty_97.5%',
            ylabel='Yield Strength $\\tau_y$ (Pa)',
            title='1971 Teneguía: Yield Strength vs Station Number',
            filename='yield_strength_plot_teneguia.png')

print("Plots generated: viscosity_plot_sanjuan.png, viscosity_plot_teneguia.png, "
      "yield_strength_plot_sanjuan.png, yield_strength_plot_teneguia.png")
