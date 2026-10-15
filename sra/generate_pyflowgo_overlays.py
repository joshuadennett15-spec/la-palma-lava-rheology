"""
PyFlowGo Overlay Plot Generator — Resolution Comparison Update
John - PyFlowGo Plotter Agent
Updated: 2026-03-04

Loads PyFlowGo simulation CSVs and overlays viscosity
and yield-strength profiles for all three DEM resolutions (2 m, 10 m, 30 m)
plus the NEW calibrated best-fit run (Q=15.15 m³/s, 2m step, reaches 2852 m)
onto the reference resolution comparison images for:
  - 1971 Teneguia (Fuencaliente, La Palma)
  - 1949 San Juan / Llano del Banco (La Palma)

Resolution comparison runs (ER = 12.0 m³/s, three DEMs):
  1971: results_flowgo_La Palma 1971 Teneguia (res 2m)_12.0m3s.csv
        results_flowgo_La Palma 1971 Teneguia (res 10m)_12.0m3s.csv
        results_flowgo_La Palma 1971 Teneguia (res 30m)_12.0m3s.csv

NEW calibrated best-fit run (2026-03-04, Q=15.15 m³/s, 2m step, 2852 m):
  1971: results_main_flowgo_La Palma 1971 Teneguia_15.15m3s.csv
        pos 0-2852 m, true flow length 2850 m

1949 San Juan — resolution comparison runs (ER = 22.0 m³/s):
  results_flowgo_La Palma 1949 (San Juan - Llano del Banco) (res 2m)_22.0m3s.csv
  results_flowgo_La Palma 1949 (San Juan - Llano del Banco) (res 10m)_22.0m3s.csv
  results_flowgo_La Palma 1949 (San Juan - Llano del Banco) (res 30m)_22.0m3s.csv

DEM resolution line styles (matching existing reference image convention):
  2 m DSM  — solid black   (#1a1a1a)
  10 m DEM — dashed gray   (#555555)
  30 m DEM — dotted pink   (#e07b7b)
  NEW best-fit (15.15 m³/s) — solid orange (#e67e22), lw=2.4
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
from matplotlib.image import imread
from scipy import stats
from datetime import datetime

# ── Matplotlib style settings ────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.size":        12,
    "axes.labelsize":   12,
    "axes.titlesize":   13,
    "legend.fontsize":  9.5,
    "xtick.labelsize":  11,
    "ytick.labelsize":  11,
    "figure.dpi":       150,
})

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE        = Path(__file__).resolve().parent.parent
SRA_DIR     = os.path.join(BASE, "sra")
RES_DIR     = os.path.join(BASE, "pyflowgo", "results_flowgo")
RESOURCE    = os.path.join(BASE, "pyflowgo", "resource")

REF_1971    = os.path.join(SRA_DIR, "resolution_comparison_1971.png")
REF_1949    = os.path.join(SRA_DIR, "resolution_comparison_1949.png")
MC_CSV      = os.path.join(SRA_DIR, "LaPalmaData_MC_Results.csv")
SP_1971     = os.path.join(RESOURCE, "LaPalma1971_slope_profile.txt")
SP_1949     = os.path.join(RESOURCE, "LaPalma1949_slope_profile.txt")

# Calibrated effusion rates used for resolution comparison runs
ER_1971     = 12.0   # m3/s
ER_1949     = 22.0   # m3/s

# NEW calibrated best-fit run parameters
ER_1971_NEW = 15.15   # m3/s  — new run, Q=15.15 m³/s, 2m step, reaches 2852 m
TRUE_LEN_1971 = 2850  # m     — true observed flow length for 1971 Teneguia

# PyFlowGo CSV files — three resolutions for each event
def pfg_path(event_label, res_label, er):
    """Construct result CSV path from event/resolution/ER labels."""
    er_str = f"{er:.1f}" if er != int(er) else str(er)
    return os.path.join(RES_DIR, f"results_flowgo_{event_label} (res {res_label})_{er_str}m3s.csv")

LABEL_1971 = "La Palma 1971 Teneguia"
LABEL_1949 = "La Palma 1949 (San Juan - Llano del Banco)"

CSV_1971 = {
    "2 m DSM":  pfg_path(LABEL_1971, "2m",  ER_1971),
    "10 m DEM": pfg_path(LABEL_1971, "10m", ER_1971),
    "30 m DEM": pfg_path(LABEL_1971, "30m", ER_1971),
}
CSV_1949 = {
    "2 m DSM":  pfg_path(LABEL_1949, "2m",  ER_1949),
    "10 m DEM": pfg_path(LABEL_1949, "10m", ER_1949),
    "30 m DEM": pfg_path(LABEL_1949, "30m", ER_1949),
}

# New best-fit run CSV path (no resolution label in filename)
CSV_1971_NEW = os.path.join(RES_DIR,
    "results_main_flowgo_La Palma 1971 Teneguia_15.15m3s.csv")

# Line style palette — matches existing reference image convention
RES_STYLES = {
    "2 m DSM":  dict(color="#1a1a1a", lw=1.8, ls="-",   zorder=5),
    "10 m DEM": dict(color="#555555", lw=1.8, ls="--",  zorder=4),
    "30 m DEM": dict(color="#e07b7b", lw=1.8, ls=":",   zorder=3),
}

# New best-fit run style — prominent orange solid line
STYLE_NEW = dict(color="#e67e22", lw=2.4, ls="-", zorder=6)

# SRA field data colours
COLOR_SRA_1971 = "#2471a3"   # blue
COLOR_SRA_1949 = "#c0392b"   # red

TIMESTAMP   = datetime.now().strftime("%Y%m%d_%H%M%S")
OUT_1971    = os.path.join(SRA_DIR, f"pyflowgo_overlay_1971_{TIMESTAMP}.png")
OUT_1949    = os.path.join(SRA_DIR, f"pyflowgo_overlay_1949_{TIMESTAMP}.png")

# Vent elevations (m asl)
VENT_ELEV_1971 = 430.0
VENT_ELEV_1949 = 1575.0

# ── Helper: load and lightly clean a PyFlowGo CSV ────────────────────────────
VISC_CLIP = 1e10   # Pa·s — cap numerical blow-up artefacts above this threshold
THO_CLIP  = 1e7    # Pa

def load_pfg(path, label):
    """
    Load PyFlowGo CSV; clip extreme numerical blow-up values (> VISC_CLIP / THO_CLIP)
    that arise from a single runaway step at the flow terminus.  These artefacts
    compress the y-axis and hide the rest of the curve.  Report the raw max
    before and after clipping.
    """
    df = pd.read_csv(path)
    visc_raw_max = df["viscosity"].max()
    tho_raw_max  = df["tho_0"].max()
    # Clip only the runaway terminus rows — keep all other data intact
    df = df[(df["viscosity"] <= VISC_CLIP) & (df["tho_0"] <= THO_CLIP)].copy()
    clipped_visc = visc_raw_max > VISC_CLIP
    clipped_tho  = tho_raw_max  > THO_CLIP
    clip_note = ""
    if clipped_visc or clipped_tho:
        clip_note = (f"  [clipped visc>{VISC_CLIP:.0e}: {visc_raw_max:.2e}; "
                     f"tho0>{THO_CLIP:.0e}: {tho_raw_max:.2e}]")
    print(f"  {label}: {len(df)} rows | "
          f"pos {df['position'].min():.0f}-{df['position'].max():.0f} m | "
          f"visc {df['viscosity'].min():.2e}-{df['viscosity'].max():.2e} Pa.s | "
          f"tho0 {df['tho_0'].min():.2e}-{df['tho_0'].max():.2e} Pa"
          + clip_note)
    return df

# ── Helper: elevation-based distance-from-vent mapping ───────────────────────
def match_elevations_to_distance(sp_path, vent_elev, target_elevs, step=2.0):
    """
    Reconstruct elevation profile by integrating slope from vent.
    For each target elevation, return the nearest matching distance (m).
    """
    sp     = pd.read_csv(sp_path, sep="\t")
    slopes = sp["Slope (deg)"].values
    dists  = sp["Distance (m)"].values
    elev   = [vent_elev]
    for s in slopes[:-1]:
        dz = -step * np.sin(np.radians(s))
        elev.append(elev[-1] + dz)
    elev = np.array(elev)
    result = []
    for te in target_elevs:
        idx = np.argmin(np.abs(elev - te))
        result.append(dists[idx])
    return np.array(result, dtype=float)

# ── Helper: log-linear trend fit ─────────────────────────────────────────────
def fit_loglinear(x, y):
    """Fit log10(y) ~ a*x + b. Returns x_smooth, y_smooth, stats_string."""
    valid = y > 0
    if valid.sum() < 3:
        return None, None, None, None, ""
    xv, yv          = x[valid], np.log10(y[valid])
    slope, ic, r, p, se = stats.linregress(xv, yv)
    x_s  = np.linspace(xv.min(), xv.max(), 400)
    y_s  = 10 ** (ic + slope * x_s)
    p_str = "p<0.001" if p < 0.001 else f"p={p:.3f}"
    return x_s, y_s, slope, ic, f"R\u00b2={r**2:.4f}, {p_str}, SE={se:.5f}"

# ── Helper: add trend + extrapolation to an axis ──────────────────────────────
def add_trend(ax, x, y, color, x_max_data):
    """Fit + plot log-linear trend on ax; return stats_string (or '')."""
    result = fit_loglinear(x, y)
    if result[0] is None:
        return ""
    x_s, y_s, slope, ic, stats_str = result
    x_ext = np.linspace(x_s.min(), x_max_data * 1.35, 400)
    y_ext = 10 ** (ic + slope * x_ext)
    mask  = x_ext <= x_max_data
    ax.plot(x_ext[mask],  y_ext[mask],  "-",  color=color, lw=2.0, label="SRA trend")
    ax.plot(x_ext[~mask], y_ext[~mask], ":",  color=color, lw=1.8, alpha=0.7,
            label="trend extrapolation")
    return stats_str

# ── Helper: configure log-scale axis ─────────────────────────────────────────
def configure_log_axis(ax):
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=12))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0,
        subs=np.arange(2, 10) * 0.1, numticks=12))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.grid(True, which="major", ls="-",  alpha=0.40)
    ax.grid(True, which="minor", ls=":",  alpha=0.25)

# ── Helper: add regression stats text box ────────────────────────────────────
def add_stats_box(ax, stats_str, color):
    if stats_str:
        ax.text(0.97, 0.04, stats_str,
                transform=ax.transAxes, fontsize=9, ha="right", va="bottom",
                family="monospace",
                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=color, alpha=0.88))

# ── Helper: add vertical true flow length marker ──────────────────────────────
def add_true_length_marker(ax, true_len, label="True terminus\n(2850 m)"):
    ax.axvline(true_len, color="#e67e22", lw=1.4, ls="--", alpha=0.75,
               label=label, zorder=7)

# ══════════════════════════════════════════════════════════════════════════════
# Load SRA field data
# ══════════════════════════════════════════════════════════════════════════════
mc       = pd.read_csv(MC_CSV)
sanjuan  = mc[mc["Eruption"] == "1949_San_Juan"].sort_values("Site").copy()
teneguia = mc[mc["Eruption"] == "1971_Teneguia"].sort_values("Site").copy()

print("Assigning distance-from-vent via elevation matching...")
teneguia_dist = match_elevations_to_distance(
    SP_1971, VENT_ELEV_1971, teneguia["Elevation (m)"].values, step=10.0)
sanjuan_dist  = match_elevations_to_distance(
    SP_1949, VENT_ELEV_1949, sanjuan["Elevation (m)"].values,  step=10.0)

print(f"  1971 Teneguia site distances: {teneguia_dist.min():.0f}-{teneguia_dist.max():.0f} m")
print(f"  1949 San Juan  site distances: {sanjuan_dist.min():.0f}-{sanjuan_dist.max():.0f} m")

# ══════════════════════════════════════════════════════════════════════════════
# Load PyFlowGo CSVs — resolution comparison runs (all three resolutions)
# ══════════════════════════════════════════════════════════════════════════════
print("\nLoading 1971 Teneguia PyFlowGo resolution-comparison runs...")
pfg_1971 = {}
for res_label, csv_path in CSV_1971.items():
    pfg_1971[res_label] = load_pfg(csv_path, res_label)

print("\nLoading NEW 1971 calibrated best-fit run (Q=15.15 m³/s, 2852 m)...")
pfg_1971_new = load_pfg(CSV_1971_NEW, f"NEW best-fit {ER_1971_NEW} m3/s")

print("\nLoading 1949 San Juan PyFlowGo runs...")
pfg_1949 = {}
for res_label, csv_path in CSV_1949.items():
    pfg_1949[res_label] = load_pfg(csv_path, res_label)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 1 — 1971 Teneguia: all-resolution overlay + new best-fit run
# ══════════════════════════════════════════════════════════════════════════════
print(f"\nGenerating 1971 overlay figure...")

ref_img_1971 = imread(REF_1971)
fig1 = plt.figure(figsize=(22, 11))
gs1  = GridSpec(2, 2, figure=fig1, width_ratios=[1.15, 1],
                wspace=0.06, hspace=0.38)
ax_ref  = fig1.add_subplot(gs1[:, 0])
ax_vis  = fig1.add_subplot(gs1[0, 1])
ax_ys   = fig1.add_subplot(gs1[1, 1])

# -- Reference image --
ax_ref.imshow(ref_img_1971)
ax_ref.set_title(
    "Reference: 1971 Teneguia — Resolution Uncertainty Comparison\n"
    "(Updated: New best-fit run Q=15.15 m\u00b3/s, 2 m step, terminus 2852 m)",
    fontsize=11, fontweight="bold", pad=8)
ax_ref.axis("off")

# -- Viscosity panel --
# SRA field data
ax_vis.fill_between(teneguia_dist,
    teneguia["MC_Eta_2.5%"], teneguia["MC_Eta_97.5%"],
    color=COLOR_SRA_1971, alpha=0.12, label="_nolegend_")
ax_vis.errorbar(teneguia_dist, teneguia["MC_Eta_Median"],
    yerr=[teneguia["MC_Eta_Median"] - teneguia["MC_Eta_2.5%"],
          teneguia["MC_Eta_97.5%"]  - teneguia["MC_Eta_Median"]],
    fmt="o", color=COLOR_SRA_1971, markerfacecolor="white", markeredgewidth=1.5,
    capsize=3, lw=1.2, ms=6, label="SRA field data (95% CI)")

ss_vis = add_trend(ax_vis, teneguia_dist,
                   teneguia["MC_Eta_Median"].values,
                   COLOR_SRA_1971, teneguia_dist.max())

# PyFlowGo resolution-comparison curves
for res_label, df in pfg_1971.items():
    style = RES_STYLES[res_label]
    ax_vis.plot(df["position"], df["viscosity"],
                label=f"{res_label} ({ER_1971} m\u00b3/s)",
                **style)

# NEW best-fit run — prominent orange line
ax_vis.plot(pfg_1971_new["position"], pfg_1971_new["viscosity"],
            label=f"NEW best-fit ({ER_1971_NEW} m\u00b3/s, 2 m step)",
            **STYLE_NEW)

# True flow length dashed marker
add_true_length_marker(ax_vis, TRUE_LEN_1971, label=f"True terminus ({TRUE_LEN_1971} m)")

configure_log_axis(ax_vis)
ax_vis.set_ylabel(r"Apparent Viscosity $\eta_A$ (Pa$\cdot$s)", fontsize=12)
ax_vis.set_xlabel("Distance from Vent (m)", fontsize=12)
ax_vis.set_title("(a)  Viscosity — PyFlowGo (all resolutions + new run) vs SRA", fontsize=12, pad=5)
ax_vis.legend(fontsize=9, loc="upper left", framealpha=0.92)
add_stats_box(ax_vis, ss_vis, COLOR_SRA_1971)

# -- Yield Strength panel --
ax_ys.fill_between(teneguia_dist,
    teneguia["MC_Ty_2.5%"], teneguia["MC_Ty_97.5%"],
    color=COLOR_SRA_1971, alpha=0.12, label="_nolegend_")
ax_ys.errorbar(teneguia_dist, teneguia["MC_Ty_Median"],
    yerr=[teneguia["MC_Ty_Median"] - teneguia["MC_Ty_2.5%"],
          teneguia["MC_Ty_97.5%"]  - teneguia["MC_Ty_Median"]],
    fmt="s", color=COLOR_SRA_1971, markerfacecolor="white", markeredgewidth=1.5,
    capsize=3, lw=1.2, ms=6, label="SRA field data (95% CI)")

ss_ys = add_trend(ax_ys, teneguia_dist,
                  teneguia["MC_Ty_Median"].values,
                  COLOR_SRA_1971, teneguia_dist.max())

# Resolution-comparison curves
for res_label, df in pfg_1971.items():
    style = RES_STYLES[res_label]
    ax_ys.plot(df["position"], df["tho_0"],
               label=f"{res_label} ({ER_1971} m\u00b3/s)",
               **style)

# NEW best-fit run
ax_ys.plot(pfg_1971_new["position"], pfg_1971_new["tho_0"],
           label=f"NEW best-fit ({ER_1971_NEW} m\u00b3/s, 2 m step)",
           **STYLE_NEW)

# True flow length marker
add_true_length_marker(ax_ys, TRUE_LEN_1971, label=f"True terminus ({TRUE_LEN_1971} m)")

configure_log_axis(ax_ys)
ax_ys.set_ylabel(r"Yield Strength $\tau_y$ (Pa)", fontsize=12)
ax_ys.set_xlabel("Distance from Vent (m)", fontsize=12)
ax_ys.set_title("(b)  Yield Strength — PyFlowGo (all resolutions + new run) vs SRA", fontsize=12, pad=5)
ax_ys.legend(fontsize=9, loc="upper left", framealpha=0.92)
add_stats_box(ax_ys, ss_ys, COLOR_SRA_1971)

fig1.suptitle(
    f"1971 Teneguia — PyFlowGo Resolution Comparison + NEW Calibrated Run  |  "
    f"NEW: Q = {ER_1971_NEW} m\u00b3/s, 2 m step, terminus 2852 m (true: {TRUE_LEN_1971} m)\n"
    f"Resolution runs: Q = {ER_1971} m\u00b3/s  |  "
    f"Source: results_main_flowgo_La Palma 1971 Teneguia_15.15m3s.csv  |  "
    f"Generated: {TIMESTAMP}",
    fontsize=12, fontweight="bold", y=1.01)

fig1.savefig(OUT_1971, dpi=300, bbox_inches="tight")
plt.close(fig1)
print(f"Saved: {OUT_1971}")


# ══════════════════════════════════════════════════════════════════════════════
# PLOT 2 — 1949 San Juan: all-resolution overlay (unchanged)
# ══════════════════════════════════════════════════════════════════════════════
print(f"\nGenerating 1949 overlay figure...")

ref_img_1949 = imread(REF_1949)
fig2 = plt.figure(figsize=(22, 11))
gs2  = GridSpec(2, 2, figure=fig2, width_ratios=[1.15, 1],
                wspace=0.06, hspace=0.38)
ax_ref2 = fig2.add_subplot(gs2[:, 0])
ax_vis2 = fig2.add_subplot(gs2[0, 1])
ax_ys2  = fig2.add_subplot(gs2[1, 1])

# -- Reference image --
ax_ref2.imshow(ref_img_1949)
ax_ref2.set_title(
    "Reference: 1949 San Juan (Llano del Banco)\n"
    "Resolution Uncertainty Comparison (2 m / 10 m / 30 m DEM)",
    fontsize=11, fontweight="bold", pad=8)
ax_ref2.axis("off")

# -- Viscosity panel --
ax_vis2.fill_between(sanjuan_dist,
    sanjuan["MC_Eta_2.5%"], sanjuan["MC_Eta_97.5%"],
    color=COLOR_SRA_1949, alpha=0.12, label="_nolegend_")
ax_vis2.errorbar(sanjuan_dist, sanjuan["MC_Eta_Median"],
    yerr=[sanjuan["MC_Eta_Median"] - sanjuan["MC_Eta_2.5%"],
          sanjuan["MC_Eta_97.5%"]  - sanjuan["MC_Eta_Median"]],
    fmt="o", color=COLOR_SRA_1949, markerfacecolor="white", markeredgewidth=1.5,
    capsize=3, lw=1.2, ms=6, label="SRA field data (95% CI)")

ss_vis2 = add_trend(ax_vis2, sanjuan_dist,
                    sanjuan["MC_Eta_Median"].values,
                    COLOR_SRA_1949, sanjuan_dist.max())

for res_label, df in pfg_1949.items():
    style = RES_STYLES[res_label]
    ax_vis2.plot(df["position"], df["viscosity"],
                 label=f"{res_label} ({ER_1949} m\u00b3/s)",
                 **style)

configure_log_axis(ax_vis2)
ax_vis2.set_ylabel(r"Apparent Viscosity $\eta_A$ (Pa$\cdot$s)", fontsize=12)
ax_vis2.set_xlabel("Distance from Vent (m)", fontsize=12)
ax_vis2.set_title("(a)  Viscosity — PyFlowGo (all resolutions) vs SRA", fontsize=12, pad=5)
ax_vis2.legend(fontsize=9, loc="upper left", framealpha=0.92)
add_stats_box(ax_vis2, ss_vis2, COLOR_SRA_1949)

# -- Yield Strength panel --
ax_ys2.fill_between(sanjuan_dist,
    sanjuan["MC_Ty_2.5%"], sanjuan["MC_Ty_97.5%"],
    color=COLOR_SRA_1949, alpha=0.12, label="_nolegend_")
ax_ys2.errorbar(sanjuan_dist, sanjuan["MC_Ty_Median"],
    yerr=[sanjuan["MC_Ty_Median"] - sanjuan["MC_Ty_2.5%"],
          sanjuan["MC_Ty_97.5%"]  - sanjuan["MC_Ty_Median"]],
    fmt="s", color=COLOR_SRA_1949, markerfacecolor="white", markeredgewidth=1.5,
    capsize=3, lw=1.2, ms=6, label="SRA field data (95% CI)")

ss_ys2 = add_trend(ax_ys2, sanjuan_dist,
                   sanjuan["MC_Ty_Median"].values,
                   COLOR_SRA_1949, sanjuan_dist.max())

for res_label, df in pfg_1949.items():
    style = RES_STYLES[res_label]
    ax_ys2.plot(df["position"], df["tho_0"],
                label=f"{res_label} ({ER_1949} m\u00b3/s)",
                **style)

configure_log_axis(ax_ys2)
ax_ys2.set_ylabel(r"Yield Strength $\tau_y$ (Pa)", fontsize=12)
ax_ys2.set_xlabel("Distance from Vent (m)", fontsize=12)
ax_ys2.set_title("(b)  Yield Strength — PyFlowGo (all resolutions) vs SRA", fontsize=12, pad=5)
ax_ys2.legend(fontsize=9, loc="upper right", framealpha=0.92)
add_stats_box(ax_ys2, ss_ys2, COLOR_SRA_1949)

fig2.suptitle(
    f"1949 San Juan (Llano del Banco) — PyFlowGo Resolution Comparison Overlay  |  ER = {ER_1949} m\u00b3/s\n"
    f"Most Recent 2 m Run: results_flowgo_La Palma 1949 (San Juan - Llano del Banco) (res 2m)_22.0m3s.csv  |  Generated: {TIMESTAMP}",
    fontsize=13, fontweight="bold", y=1.01)

fig2.savefig(OUT_1949, dpi=300, bbox_inches="tight")
plt.close(fig2)
print(f"Saved: {OUT_1949}")

print("\nAll overlay plots generated successfully.")
print(f"  1971 overlay : {OUT_1971}")
print(f"  1949 overlay : {OUT_1949}")
