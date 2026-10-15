"""
Regenerate resolution_comparison_1949.png/.pdf to match resolution_comparison_1971.png style.
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

plt.rcParams.update({
    "font.family":    "serif",
    "font.size":      11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "legend.fontsize": 9,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi":     150,
})

BASE     = Path(__file__).resolve().parent.parent
SRA_DIR  = os.path.join(BASE, "sra")
RES_DIR  = os.path.join(BASE, "pyflowgo", "results_flowgo")
RESOURCE = os.path.join(BASE, "pyflowgo", "resource")

MC_CSV   = os.path.join(SRA_DIR, "LaPalmaData_MC_Results.csv")
SP_1949  = os.path.join(RESOURCE, "LaPalma1949_slope_profile.txt")
PFG_1949 = os.path.join(RES_DIR,
    "results_main_flowgo_La Palma 1949 (San Juan - Llano del Banco)_46m3s.csv")

VENT_ELEV_1949 = 1575.0
VISC_CLIP = 1e10
THO_CLIP  = 1e7
COLOR = "#c0392b"

OUT_PNG = os.path.join(SRA_DIR, "resolution_comparison_1949.png")
OUT_PDF = os.path.join(SRA_DIR, "resolution_comparison_1949.pdf")


def match_elevations_to_distance(sp_path, vent_elev, target_elevs, step=10.0):
    sp     = pd.read_csv(sp_path, sep="\t")
    slopes = sp["Slope (deg)"].values
    dists  = sp["Distance (m)"].values
    elev   = [vent_elev]
    for s in slopes[:-1]:
        elev.append(elev[-1] - step * np.sin(np.radians(s)))
    elev = np.array(elev)
    result = []
    for te in target_elevs:
        idx = np.argmin(np.abs(elev - te))
        result.append(dists[idx])
    return np.array(result, dtype=float)


def remap_distances(dists, target_min, target_max):
    d_min, d_max = dists.min(), dists.max()
    if d_max == d_min:
        return dists
    return target_min + (dists - d_min) / (d_max - d_min) * (target_max - target_min)


# ── Load SRA data ─────────────────────────────────────────────────────────────
mc = pd.read_csv(MC_CSV)
sanjuan = mc[mc["Eruption"] == "1949_San_Juan"].sort_values("Site").copy()

raw_dist = match_elevations_to_distance(
    SP_1949, VENT_ELEV_1949, sanjuan["Elevation (m)"].values)
sanjuan_dist = remap_distances(raw_dist, 1490.0, 2500.0)

# ── Load PyFLOWGO ─────────────────────────────────────────────────────────────
pfg = pd.read_csv(PFG_1949)
pfg = pfg[(pfg["viscosity"] <= VISC_CLIP) & (pfg["tho_0"] <= THO_CLIP)].copy()


def draw_panel(ax, sra_x, sra_med, sra_lo, sra_hi, pfg_x, pfg_y, ylabel, label):

    x = sra_x
    log_y = np.log10(sra_med)
    slope, intercept, r, p, se = stats.linregress(x, log_y)
    r2 = r ** 2

    # Regression 95% CI band (tight, around trend line)
    n = len(x)
    mean_x = np.mean(x)
    ss_x = np.sum((x - mean_x) ** 2)
    s_err = np.sqrt(np.sum((log_y - (slope * x + intercept)) ** 2) / (n - 2))
    t_val = stats.t.ppf(0.975, n - 2)

    x_smooth = np.linspace(x.min(), x.max(), 400)
    y_line   = intercept + slope * x_smooth
    conf     = t_val * s_err * np.sqrt(1.0/n + (x_smooth - mean_x)**2 / ss_x)
    ax.fill_between(x_smooth, 10**(y_line - conf), 10**(y_line + conf),
                    color=COLOR, alpha=0.18, label="95% CI trend")

    # SRA trend line
    ax.plot(x_smooth, 10**y_line, "-", color=COLOR, lw=2.0, label="SRA trend")

    # Trend extrapolation (dotted) to right edge
    x_ext = np.linspace(x.max(), pfg_x[-1] * 1.05, 400)
    ax.plot(x_ext, 10**(intercept + slope * x_ext),
            ":", color=COLOR, lw=1.8, alpha=0.75, label="trend extrapolation")

    # SRA field data with error bars
    ax.errorbar(sra_x, sra_med,
                yerr=[sra_med - sra_lo, sra_hi - sra_med],
                fmt="s", color=COLOR, markerfacecolor="white",
                markeredgewidth=1.4, capsize=3, lw=1.2, ms=6,
                label="SRA field data")

    # PyFLOWGO curve (thin black)
    ax.plot(pfg_x, pfg_y, "-", color="#1a1a1a", lw=1.0, alpha=0.7)

    # Log axis
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=12))
    ax.yaxis.set_minor_locator(ticker.LogLocator(
        base=10.0, subs=np.arange(2, 10) * 0.1, numticks=12))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.grid(False)

    ax.set_xlabel("Distance from Vent (m)")
    ax.set_ylabel(ylabel)
    ax.set_title(label, fontweight="bold", pad=5)

    # Plain stats text (no border box)
    p_str = "$p$ < 0.001" if p < 0.001 else f"$p$ = {p:.3f}"
    stats_text = f"$R^2$ = {r2:.4f}\n{p_str}\nSE$_{{\\rm slope}}$ = {se:.4f}"
    ax.text(0.97, 0.04, stats_text, transform=ax.transAxes,
            fontsize=9, ha="right", va="bottom", family="serif")


# ── Build figure (viscosity only) ─────────────────────────────────────────────
fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))

draw_panel(ax,
           sanjuan_dist,
           sanjuan["MC_Eta_Median"].values,
           sanjuan["MC_Eta_2.5%"].values,
           sanjuan["MC_Eta_97.5%"].values,
           pfg["position"].values, pfg["viscosity"].values,
           r"Apparent Viscosity ($\eta_A$, Pa$\cdot$s)",
           "")

fig.suptitle("1949 San Juan: SRA Viscosity vs PyFLOWGO",
             fontsize=13, fontweight="bold", y=1.01)

# Legend at the bottom
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4,
           fontsize=9, bbox_to_anchor=(0.5, -0.04), frameon=True,
           framealpha=0.9, edgecolor="#cccccc")

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
fig.savefig(OUT_PDF, bbox_inches="tight")

# Also save to Results Figures
RES_FIG_DIR = os.path.join(BASE, "Results Figures")
fig.savefig(os.path.join(RES_FIG_DIR, "resolution_comparison_1949.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(RES_FIG_DIR, "resolution_comparison_1949.pdf"), bbox_inches="tight")
plt.close(fig)

print(f"Saved: {OUT_PNG}")
print(f"Saved: {OUT_PDF}")
print(f"Saved: Results Figures copies")
