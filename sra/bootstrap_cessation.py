"""
Bootstrap Cessation Prediction Analysis
========================================
Tests whether SRA field-data log-linear trends can predict the rheological
cessation conditions found by PyFLOWGO.

Four independent analyses:
  1. 1971 Teneguía — Apparent Viscosity
  2. 1971 Teneguía — Yield Strength
  3. 1949 San Juan  — Apparent Viscosity
  4. 1949 San Juan  — Yield Strength

Method (per analysis):
  1. Fit OLS: log10(y) = a·x + b  to SRA field stations.
  2. Compute residuals r_i.
  3. For B=10 000 bootstrap iterations:
       - resample residuals with replacement → r*
       - construct pseudo-obs:  log10(y*_i) = (a·x_i + b) + r*_i
       - re-fit OLS → a*, b*
       - extrapolate:  d_cessation* = (log10(cessation_value) - b*) / a*
  4. Report median, 95 % CI, and whether the true PyFLOWGO cessation
     distance falls inside that interval.

Author: Antigravity Agent  |  2026-03-04
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from datetime import datetime

# ── Matplotlib style ─────────────────────────────────────────────────────────
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

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE     = Path(__file__).resolve().parent.parent
SRA_DIR  = os.path.join(BASE, "sra")
RES_DIR  = os.path.join(BASE, "pyflowgo", "results_flowgo")
RESOURCE = os.path.join(BASE, "pyflowgo", "resource")

MC_CSV   = os.path.join(SRA_DIR, "LaPalmaData_MC_Results.csv")
SP_1971  = os.path.join(RESOURCE, "LaPalma1971_slope_profile.txt")
SP_1949  = os.path.join(RESOURCE, "LaPalma1949_slope_profile.txt")

PFG_1971 = os.path.join(RES_DIR,
    "results_main_flowgo_La Palma 1971 Teneguia_15.15m3s.csv")
PFG_1949 = os.path.join(RES_DIR,
    "results_main_flowgo_La Palma 1949 (San Juan - Llano del Banco)_46m3s.csv")

VENT_ELEV_1971 = 430.0   # m asl
VENT_ELEV_1949 = 1575.0  # m asl

# ── Clipping thresholds (same as generate_pyflowgo_overlays.py) ──────────────
VISC_CLIP = 1e10   # Pa·s
THO_CLIP  = 1e7    # Pa

# ── Bootstrap parameters ─────────────────────────────────────────────────────
B   = 10_000       # number of bootstrap iterations
RNG = np.random.default_rng(42)

# ── Output timestamp ─────────────────────────────────────────────────────────
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")


# ══════════════════════════════════════════════════════════════════════════════
# Helper functions
# ══════════════════════════════════════════════════════════════════════════════

def match_elevations_to_distance(sp_path, vent_elev, target_elevs, step=10.0):
    """
    Reconstruct elevation profile by integrating slope from vent.
    For each target elevation, return the nearest matching distance (m).
    (Reused from generate_pyflowgo_overlays.py)
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

def remap_distances(dists, target_min, target_max):
    """Linearly stretch an array of distances to fit exactly within [target_min, target_max]."""
    d_min = dists.min()
    d_max = dists.max()
    if d_max == d_min: return dists
    return target_min + (dists - d_min) / (d_max - d_min) * (target_max - target_min)

def load_pyflowgo(path, visc_clip=VISC_CLIP, tho_clip=THO_CLIP):
    """Load and lightly clip a PyFLOWGO CSV; return DataFrame."""
    df = pd.read_csv(path)
    df = df[(df["viscosity"] <= visc_clip) & (df["tho_0"] <= tho_clip)].copy()
    return df


def interpolate_pfg(df, col, step=1.0):
    """Interpolate a PyFLOWGO column to 1-m resolution; return (pos, vals)."""
    pos = df["position"].values
    vals = df[col].values
    f = interp1d(pos, vals, kind="linear", fill_value="extrapolate")
    pos_fine = np.arange(pos.min(), pos.max() + step, step)
    return pos_fine, f(pos_fine)


def cessation_point(df, col):
    """
    Return (cessation_distance, cessation_value) from clipped PyFLOWGO.
    Cessation = last position in the simulation (flow has stopped).
    """
    pos_fine, vals_fine = interpolate_pfg(df, col)
    return pos_fine[-1], vals_fine[-1]


def bootstrap_cessation(distance, log_y, cessation_value, B=B):
    """
    Paired-residual bootstrap.

    Parameters
    ----------
    distance : array-like   — SRA station distances (m)
    log_y    : array-like   — log10 of observed rheological variable
    cessation_value : float — PyFLOWGO value at cessation (natural units)
    B : int                 — number of bootstrap iterations

    Returns
    -------
    dict with keys:
        slope, intercept, r2, p, se          — observed OLS fit
        boot_distances   — array of predicted cessation distances (may have NaN)
        n_valid          — count of valid (positive-slope) iterations
    """
    x = np.asarray(distance, dtype=float)
    y = np.asarray(log_y, dtype=float)

    # Observed OLS fit
    slope, intercept, r, p, se = stats.linregress(x, y)
    residuals = y - (slope * x + intercept)

    log_cess = np.log10(cessation_value)

    boot_d = np.full(B, np.nan)
    for i in range(B):
        # Resample residuals with replacement
        r_star = RNG.choice(residuals, size=len(residuals), replace=True)
        y_star = (slope * x + intercept) + r_star
        # Re-fit
        a_star, b_star, *_ = stats.linregress(x, y_star)
        if a_star > 0:
            boot_d[i] = (log_cess - b_star) / a_star

    return {
        "slope": slope,
        "intercept": intercept,
        "r2": r ** 2,
        "p": p,
        "se": se,
        "boot_distances": boot_d,
        "n_valid": int(np.isfinite(boot_d).sum()),
    }


def print_results(label, result, true_cessation_dist):
    """Pretty-print one analysis block."""
    bd = result["boot_distances"]
    valid = bd[np.isfinite(bd)]

    print(f"\n{'=' * 72}")
    print(f"  {label}")
    print(f"{'=' * 72}")
    print(f"  Observed OLS fit:")
    print(f"    slope     = {result['slope']:.8f}")
    print(f"    intercept = {result['intercept']:.4f}")
    print(f"    R²        = {result['r2']:.6f}")
    p = result["p"]
    print(f"    p-value   = {'<0.001' if p < 0.001 else f'{p:.4f}'}")
    print(f"    SE(slope) = {result['se']:.8f}")
    print()
    print(f"  PyFLOWGO cessation distance : {true_cessation_dist:.0f} m")
    print()

    if len(valid) == 0:
        print("  Bootstrap: NO valid iterations (all slopes <= 0)")
        return

    med   = np.median(valid)
    lo    = np.percentile(valid, 2.5)
    hi    = np.percentile(valid, 97.5)
    frac  = result["n_valid"] / len(bd) * 100

    inside = lo <= true_cessation_dist <= hi

    print(f"  Bootstrap results (B = {len(bd):,}):")
    print(f"    Valid iterations  : {result['n_valid']:,} / {len(bd):,}  ({frac:.1f} %)")
    print(f"    Median predicted  : {med:.0f} m")
    print(f"    95 % CI           : [{lo:.0f}, {hi:.0f}] m")
    print(f"    True cessation    : {true_cessation_dist:.0f} m")
    print(f"    Inside 95 % CI?   : {'YES' if inside else 'NO'}")
    print()


# ══════════════════════════════════════════════════════════════════════════════
# Load SRA field data
# ══════════════════════════════════════════════════════════════════════════════
print("Loading SRA field data...")
mc = pd.read_csv(MC_CSV)

sanjuan  = mc[mc["Eruption"] == "1949_San_Juan"].sort_values("Site").copy()
teneguia = mc[mc["Eruption"] == "1971_Teneguia"].sort_values("Site").copy()

# Map elevation → distance from vent
print("Computing station distances via elevation matching...")
raw_teneguia_dist = match_elevations_to_distance(
    SP_1971, VENT_ELEV_1971, teneguia["Elevation (m)"].values, step=10.0)
raw_sanjuan_dist  = match_elevations_to_distance(
    SP_1949, VENT_ELEV_1949, sanjuan["Elevation (m)"].values,  step=10.0)

# Remap to user-specified bounds
print("Remapping distances to user-specified bounds...")
teneguia_dist = remap_distances(raw_teneguia_dist, 250.0, 2000.0)
sanjuan_dist  = remap_distances(raw_sanjuan_dist,  1490.0, 2500.0)

print(f"  1971 Teneguia: {len(teneguia)} stations, "
      f"remapped dist {teneguia_dist.min():.0f}–{teneguia_dist.max():.0f} m")
print(f"  1949 San Juan: {len(sanjuan)} stations, "
      f"remapped dist {sanjuan_dist.min():.0f}–{sanjuan_dist.max():.0f} m")

# Log10 of corrected MC medians
ten_log_visc = np.log10(teneguia["MC_Eta_Median"].values)
ten_log_ys   = np.log10(teneguia["MC_Ty_Median"].values)
sj_log_visc  = np.log10(sanjuan["MC_Eta_Median"].values)
sj_log_ys    = np.log10(sanjuan["MC_Ty_Median"].values)


# ══════════════════════════════════════════════════════════════════════════════
# Load PyFLOWGO results and determine cessation points
# ══════════════════════════════════════════════════════════════════════════════
print("\nLoading PyFLOWGO 2 m results...")
pfg_1971 = load_pyflowgo(PFG_1971)
pfg_1949 = load_pyflowgo(PFG_1949)

print(f"  1971 Teneguia: {len(pfg_1971)} rows, "
      f"pos {pfg_1971['position'].min():.0f}–{pfg_1971['position'].max():.0f} m")
print(f"  1949 San Juan: {len(pfg_1949)} rows, "
      f"pos {pfg_1949['position'].min():.0f}–{pfg_1949['position'].max():.0f} m")

# Cessation points
cess_1971_visc_dist, cess_1971_visc_val = cessation_point(pfg_1971, "viscosity")
cess_1971_ys_dist,   cess_1971_ys_val   = cessation_point(pfg_1971, "tho_0")
cess_1949_visc_dist, cess_1949_visc_val = cessation_point(pfg_1949, "viscosity")
cess_1949_ys_dist,   cess_1949_ys_val   = cessation_point(pfg_1949, "tho_0")

print(f"\n  Cessation values (from clipped 2 m runs):")
print(f"    1971 viscosity : {cess_1971_visc_val:.2e} Pa·s  @ {cess_1971_visc_dist:.0f} m")
print(f"    1971 yield str : {cess_1971_ys_val:.2e} Pa     @ {cess_1971_ys_dist:.0f} m")
print(f"    1949 viscosity : {cess_1949_visc_val:.2e} Pa·s  @ {cess_1949_visc_dist:.0f} m")
print(f"    1949 yield str : {cess_1949_ys_val:.2e} Pa     @ {cess_1949_ys_dist:.0f} m")


# ══════════════════════════════════════════════════════════════════════════════
# Run bootstrap analyses
# ══════════════════════════════════════════════════════════════════════════════
print(f"\nRunning bootstrap analyses (B = {B:,}) ...")

# 1. 1971 Teneguía — Viscosity
res_1971_visc = bootstrap_cessation(
    teneguia_dist, ten_log_visc, cess_1971_visc_val)
print_results("1971 Teneguía — Apparent Viscosity",
              res_1971_visc, cess_1971_visc_dist)

# 2. 1971 Teneguía — Yield Strength
res_1971_ys = bootstrap_cessation(
    teneguia_dist, ten_log_ys, cess_1971_ys_val)
print_results("1971 Teneguía — Yield Strength",
              res_1971_ys, cess_1971_ys_dist)

# 3. 1949 San Juan — Viscosity
res_1949_visc = bootstrap_cessation(
    sanjuan_dist, sj_log_visc, cess_1949_visc_val)
print_results("1949 San Juan — Apparent Viscosity",
              res_1949_visc, cess_1949_visc_dist)

# 4. 1949 San Juan — Yield Strength
res_1949_ys = bootstrap_cessation(
    sanjuan_dist, sj_log_ys, cess_1949_ys_val)
print_results("1949 San Juan — Yield Strength",
              res_1949_ys, cess_1949_ys_dist)

print("\nAll four bootstrap analyses complete.")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Trend + Extrapolation panels (2 × 2)
# ══════════════════════════════════════════════════════════════════════════════
print("\nGenerating figures...")

COLOR_1971 = "#2471a3"
COLOR_1949 = "#c0392b"
COLOR_PFG  = "#1a1a1a"

analyses = [
    # (title, sra_dist, sra_log_y, pfg_df, pfg_col, cess_dist, cess_val, result, color, ylabel)
    ("(a) 1971 Teneguía — Viscosity",
     teneguia_dist, ten_log_visc, pfg_1971, "viscosity",
     cess_1971_visc_dist, cess_1971_visc_val, res_1971_visc, COLOR_1971,
     r"Apparent Viscosity $\eta_A$ (Pa$\cdot$s)"),
    ("(b) 1971 Teneguía — Yield Strength",
     teneguia_dist, ten_log_ys, pfg_1971, "tho_0",
     cess_1971_ys_dist, cess_1971_ys_val, res_1971_ys, COLOR_1971,
     r"Yield Strength $\tau_y$ (Pa)"),
    ("(c) 1949 San Juan — Viscosity",
     sanjuan_dist, sj_log_visc, pfg_1949, "viscosity",
     cess_1949_visc_dist, cess_1949_visc_val, res_1949_visc, COLOR_1949,
     r"Apparent Viscosity $\eta_A$ (Pa$\cdot$s)"),
    ("(d) 1949 San Juan — Yield Strength",
     sanjuan_dist, sj_log_ys, pfg_1949, "tho_0",
     cess_1949_ys_dist, cess_1949_ys_val, res_1949_ys, COLOR_1949,
     r"Yield Strength $\tau_y$ (Pa)"),
]

fig1, axes1 = plt.subplots(2, 2, figsize=(14, 10))
axes1 = axes1.ravel()

for idx, (title, sra_d, sra_ly, pfg_df, pfg_col,
          cd, cv, res, color, ylabel) in enumerate(analyses):
    ax = axes1[idx]

    # PyFLOWGO curve (background)
    ax.plot(pfg_df["position"], pfg_df[pfg_col],
            color=COLOR_PFG, lw=1.5, ls="-", alpha=0.5, label="PyFLOWGO (2 m)", zorder=2)

    # SRA data points
    ax.scatter(sra_d, 10**sra_ly, c=color, s=40, zorder=5,
               edgecolors="white", linewidths=0.6, label="SRA field data")

    # OLS fit line through data range
    slope, intercept = res["slope"], res["intercept"]
    x_fit = np.linspace(sra_d.min(), sra_d.max(), 300)
    ax.plot(x_fit, 10**(slope * x_fit + intercept),
            color=color, lw=2.0, ls="-", zorder=4, label="OLS trend")

    # Extrapolation to cessation (if positive slope)
    if slope > 0:
        x_ext = np.linspace(sra_d.max(), cd * 1.15, 300)
        ax.plot(x_ext, 10**(slope * x_ext + intercept),
                color=color, lw=1.8, ls=":", alpha=0.7, zorder=4,
                label="Trend extrapolation")

    # Bootstrap 95% CI band on cessation distance (vertical shaded region)
    bd = res["boot_distances"]
    valid = bd[np.isfinite(bd)]
    if len(valid) > 10:
        lo = np.percentile(valid, 2.5)
        hi = np.percentile(valid, 97.5)
        # Only draw if the range is physically meaningful
        if lo > 0 and hi < cd * 3:
            ax.axvspan(lo, hi, alpha=0.12, color=color, zorder=1,
                       label=f"Bootstrap 95% CI")

    # True PyFLOWGO cessation line
    ax.axvline(cd, color=COLOR_PFG, ls="--", lw=1.2, alpha=0.8, zorder=3,
               label=f"PyFLOWGO cessation ({cd:.0f} m)")

    # Cessation value horizontal line
    ax.axhline(cv, color="#888888", ls=":", lw=0.8, alpha=0.5, zorder=1)

    # Formatting
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=12))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0,
        subs=np.arange(2, 10) * 0.1, numticks=12))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.grid(True, which="major", ls="-", alpha=0.3)
    ax.grid(True, which="minor", ls=":", alpha=0.15)

    ax.set_xlabel("Distance from Vent (m)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11, fontweight="bold", pad=6)
    ax.legend(fontsize=7.5, loc="upper left", framealpha=0.9)

    # Stats annotation box
    p_str = "$p$ < 0.001" if res["p"] < 0.001 else f"$p$ = {res['p']:.3f}"
    stats_text = f"$R^2$ = {res['r2']:.3f}\n{p_str}\nSE = {res['se']:.2e}"
    ax.text(0.97, 0.05, stats_text, transform=ax.transAxes,
            fontsize=8, ha="right", va="bottom", family="serif",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=color, alpha=0.88))

fig1.suptitle("Bootstrap Cessation Prediction: SRA Trend Extrapolation vs PyFLOWGO",
              fontsize=13, fontweight="bold", y=0.98)
fig1.tight_layout(rect=[0, 0, 1, 0.96])

out_trend = os.path.join(SRA_DIR, f"bootstrap_trend_{TIMESTAMP}.png")
fig1.savefig(out_trend, dpi=300, bbox_inches="tight")
plt.close(fig1)
print(f"  Saved: {out_trend}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Bootstrap Histograms (2 × 2)
# ══════════════════════════════════════════════════════════════════════════════
print("\nGenerating bootstrap histogram figure...")

from scipy.stats import gaussian_kde

visc_analyses = [analyses[0], analyses[2]]  # viscosity panels only

fig2, axes2 = plt.subplots(1, 2, figsize=(14, 5))

for idx, (title, sra_d, sra_ly, pfg_df, pfg_col,
          cd, cv, res, color, ylabel) in enumerate(visc_analyses):
    ax = axes2[idx]
    bd  = res["boot_distances"]
    valid = bd[np.isfinite(bd)]

    if len(valid) < 10:
        ax.text(0.5, 0.5,
                f"Insufficient valid iterations\n({len(valid)} / {len(bd)})\n"
                "OLS slope is negative:\ntrend cannot reach cessation value.",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=10, style="italic",
                bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="#cccccc"))
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlabel("Predicted Cessation Distance (m)")
        ax.set_ylabel("Frequency")
        continue

    lo  = np.percentile(valid, 2.5)
    hi  = np.percentile(valid, 97.5)
    med = np.median(valid)
    inside = lo <= cd <= hi

    # Histogram
    ax.hist(valid, bins=60, density=True, color=color, alpha=0.45)

    # KDE
    kde  = gaussian_kde(valid)
    x_kd = np.linspace(valid.min(), valid.max(), 500)
    ax.plot(x_kd, kde(x_kd), color=color, lw=2.0)

    # 95% CI shaded under KDE
    mask = (x_kd >= lo) & (x_kd <= hi)
    ax.fill_between(x_kd, kde(x_kd), where=mask, alpha=0.25, color=color,
                    label=f"95% CI [{lo:.0f}, {hi:.0f}]")

    # PyFLOWGO cessation vertical line
    ax.axvline(cd, color=COLOR_PFG, ls="--", lw=1.8, zorder=5,
               label=f"PyFLOWGO cessation ({cd:.0f} m)")

    # Bootstrap median vertical line
    ax.axvline(med, color=color, ls="-", lw=1.5, alpha=0.8,
               label=f"Bootstrap median ({med:.0f} m)")

    # Within / Outside annotation
    ci_label = "Within 95% CI" if inside else "Outside 95% CI"
    ci_color = "#1a7a3c"      if inside else "#b03a2e"
    ax.text(0.97, 0.97, ci_label, transform=ax.transAxes,
            ha="right", va="top", fontsize=10, fontweight="bold",
            color=ci_color)

    ax.set_xlabel("Predicted Cessation Distance (m)")
    ax.set_ylabel("Density")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.legend(fontsize=7.5, loc="upper left", framealpha=0.88)

    # Cap x-axis for 1971 panel
    if idx == 0:
        ax.set_xlim(right=6000)

fig2.suptitle(f"Bootstrap Distribution of Predicted Cessation Distances (B = {B:,})",
              fontsize=13, fontweight="bold", y=0.98)
fig2.tight_layout(rect=[0, 0, 1, 0.96])

out_hist = os.path.join(SRA_DIR, f"bootstrap_histogram_{TIMESTAMP}.png")
fig2.savefig(out_hist, dpi=300, bbox_inches="tight")
plt.close(fig2)
print(f"  Saved: {out_hist}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Publication Figure (2 × 1 Viscosity focus)
# ══════════════════════════════════════════════════════════════════════════════
print("\nGenerating final publication focus figure...")

# Filter to just viscosity
pub_analyses = [analyses[0], analyses[2]] 

fig3, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharey=False)
pub_axes = [ax1, ax2]

for idx, (title, sra_d, sra_ly, pfg_df, pfg_col,
          cd, cv, res, color, ylabel) in enumerate(pub_analyses):
    ax = pub_axes[idx]

    # PyFLOWGO curve (background)
    ax.plot(pfg_df["position"], pfg_df[pfg_col],
            color=COLOR_PFG, lw=1.5, ls="-", alpha=0.4, label="PyFLOWGO (2 m)", zorder=2)

    # SRA data points
    ax.scatter(sra_d, 10**sra_ly, c=color, s=45, zorder=5,
               edgecolors="white", linewidths=0.7, label="SRA field data")

    # OLS fit line through data range
    slope, intercept = res["slope"], res["intercept"]
    x_fit = np.linspace(sra_d.min(), sra_d.max(), 300)
    ax.plot(x_fit, 10**(slope * x_fit + intercept),
            color=color, lw=2.2, ls="-", zorder=4, label="OLS trend")

    # Extrapolation to cessation
    x_ext = np.linspace(sra_d.max(), max(cd, sra_d.max()) * 1.1, 300)
    ax.plot(x_ext, 10**(slope * x_ext + intercept),
            color=color, lw=1.8, ls=":", alpha=0.7, zorder=4,
            label="Trend extrapolation")

    # Bootstrap 95% CI band on cessation distance (vertical shaded region)
    bd = res["boot_distances"]
    valid = bd[np.isfinite(bd)]
    if len(valid) > 10:
        lo = np.percentile(valid, 2.5)
        hi = np.percentile(valid, 97.5)
        if 0 < lo < cd * 5:
            ax.axvspan(lo, hi, alpha=0.12, color=color, zorder=1, label=f"Bootstrap 95% CI")

    # True PyFLOWGO cessation line (removed)

    # Formatting
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=12))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=12))
    ax.grid(True, which="major", ls="-", alpha=0.25)
    ax.grid(True, which="minor", ls=":", alpha=0.1)

    ax.set_xlabel("Distance from Vent (m)", fontweight="bold")
    ax.set_ylabel(ylabel, fontweight="bold")
    ax.set_title(title.replace("(a) ", "").replace("(c) ", ""), fontsize=12, fontweight="bold", pad=10)
    
    # Text box for R2
    p_str = "p < 0.001" if res["p"] < 0.001 else f"p = {res['p']:.3f}"
    stats_text = f"$R^2$ = {res['r2']:.3f}\n{p_str}"
    ax.text(0.96, 0.06, stats_text, transform=ax.transAxes,
            fontsize=10, ha="right", va="bottom", family="serif", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.5", fc="white", ec=color, alpha=0.9))

# Shared Legend
handles, labels = ax1.get_legend_handles_labels()
fig3.legend(handles, labels, loc='lower center', ncol=6, fontsize=9, 
           bbox_to_anchor=(0.5, -0.05), frameon=True)

fig3.suptitle("SRA Viscosity Trend Extrapolation vs PyFLOWGO Cessation",
              fontsize=14, fontweight="bold", y=1.02)
fig3.tight_layout()

out_pub = os.path.join(SRA_DIR, f"bootstrap_viscosity_publication_{TIMESTAMP}.png")
fig3.savefig(out_pub, dpi=300, bbox_inches="tight")
plt.close(fig3)
print(f"  Saved: {out_pub}")

print("\nFinal plots saved.")
