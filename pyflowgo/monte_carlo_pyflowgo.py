"""
monte_carlo_pyflowgo.py
PyFLOWGO Monte Carlo uncertainty analysis for La Palma 1949 and 1971.

Perturbs 7 input parameters using physically motivated distributions:
  - Depth H:           Normal(nominal, σ=0.75 m)
  - Width w:           Normal(nominal, σ=4.0 m)
  - Slope profile:     Normal offset(0, σ=1.25°) added to every row
  - Effusion rate Q:   Normal(nominal, σ=15% of nominal)
  - Bulk density ρ:    Normal(nominal, σ=100 kg/m³)
  - Eruptive temp T:   Triangular(mode=nominal, ±10 K)
  - Init. crystallinity Cᵢ: Uniform[10%, 20%]

Outputs 95% CI envelopes on viscosity and yield strength.
"""

import os
from pathlib import Path
import json
import csv as csv_mod
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats

os.chdir(Path(__file__).resolve().parent)
import pyflowgo.run_flowgo as run_flowgo

# ------------------------------------------------------------------ #
# Config                                                              #
# ------------------------------------------------------------------ #
N_TRIALS    = 500
RESULTS_DIR = os.path.join(os.path.abspath(''), 'results_flowgo')
SRA_OUT     = '../sra'

ERUPTIONS = {
    '1949_San_Juan': {
        'json':  'resource/LaPalma1949_calibrated.json',
        'color': '#D32F2F',
        'label': '1949 San Juan (Llano del Banco)',
        'out':   f'{SRA_OUT}/pyflowgo_mc_1949_San_Juan.png',
    },
    '1971_Teneguia': {
        'json':  'resource/LaPalma1971_calibrated.json',
        'color': '#1565C0',
        'label': '1971 Teneguía',
        'out':   f'{SRA_OUT}/pyflowgo_mc_1971_Teneguia.png',
    },
}

# Parameter uncertainties
SIGMA_DEPTH  = 0.75      # m       Normal
SIGMA_WIDTH  = 4.0       # m       Normal
SIGMA_SLOPE  = 1.25      # °       Normal offset applied uniformly to profile
SIGMA_Q_REL  = 0.15      # frac    Normal (15 % of nominal)
SIGMA_RHO    = 100.0     # kg/m³   Normal
T_HALF_RANGE = 10.0      # K       Triangular ± this value around nominal
CI_MIN, CI_MAX = 0.10, 0.20    # Uniform

GRID_POINTS  = 300       # interpolation resolution

# ------------------------------------------------------------------ #
# Helpers                                                             #
# ------------------------------------------------------------------ #
def read_slope_file(path):
    """Return (header_line, list_of_row_lists) from a slope profile txt."""
    header = None
    rows = []
    with open(path, 'r') as f:
        for line in f:
            if header is None:
                header = line
                continue
            parts = line.rstrip('\n\r').split('\t')
            if len(parts) >= 5:
                rows.append(parts)
    return header, rows


def write_perturbed_slope(src_path, dst_path, slope_offset_deg):
    """Write a slope file with every slope value shifted by slope_offset_deg."""
    header, rows = read_slope_file(src_path)
    with open(dst_path, 'w') as f:
        f.write(header)
        for parts in rows:
            s = float(parts[4]) + slope_offset_deg
            s = max(0.1, min(60.0, s))
            parts[4] = f'{s:.7f}'
            f.write('\t'.join(parts) + '\n')


def load_pf_csv(csv_path, max_visc=1e9, max_ys=1e7):
    """Read (position, viscosity, tho_0) from a PyFLOWGO results CSV."""
    pos, visc, ys = [], [], []
    try:
        with open(csv_path, 'r') as f:
            reader = csv_mod.DictReader(f)
            for row in reader:
                v = float(row['viscosity'])
                y = float(row['tho_0'])
                if v < max_visc and y < max_ys:
                    pos.append(float(row['position']))
                    visc.append(v)
                    ys.append(y)
    except Exception:
        pass
    return np.array(pos), np.array(visc), np.array(ys)


def sample_and_run(base_cfg, slope_src, trial_idx, nominal_length):
    """
    Draw one set of perturbed parameters, run PyFLOWGO, return
    (pos, visc, ys) or None if the run failed / stopped too early.
    """
    cfg = json.loads(json.dumps(base_cfg))

    # Read nominals
    d_nom = float(cfg['terrain_conditions']['depth'])
    w_nom = float(cfg['terrain_conditions']['width'])
    Q_nom = float(cfg['effusion_rate_init'])
    T_nom = float(cfg['eruption_condition']['eruption_temperature'])
    r_nom = float(cfg['lava_state']['density_dre'])

    # Sample
    d_s  = float(np.clip(np.random.normal(d_nom, SIGMA_DEPTH),  0.1,   10.0))
    w_s  = float(np.clip(np.random.normal(w_nom, SIGMA_WIDTH),  1.0,  100.0))
    Q_s  = float(np.clip(np.random.normal(Q_nom, SIGMA_Q_REL * Q_nom), 1.0, 200.0))
    r_s  = float(np.clip(np.random.normal(r_nom, SIGMA_RHO),  2100.0, 2800.0))
    T_s  = float(stats.triang.rvs(c=0.5,
                                   loc=T_nom - T_HALF_RANGE,
                                   scale=2 * T_HALF_RANGE))
    ci_s = float(np.random.uniform(CI_MIN, CI_MAX))
    sl_s = float(np.random.normal(0.0, SIGMA_SLOPE))

    # Apply to config
    cfg['terrain_conditions']['depth']              = round(d_s,  3)
    cfg['terrain_conditions']['width']              = round(w_s,  3)
    cfg['effusion_rate_init']                       = round(Q_s,  3)
    cfg['eruption_condition']['eruption_temperature'] = round(T_s, 3)
    cfg['lava_state']['density_dre']                = round(r_s,  1)
    cfg['lava_state']['crystal_fraction']           = round(ci_s, 4)

    # Temp slope file
    tmp_slope = f'resource/_mc_slope_{trial_idx}.txt'
    write_perturbed_slope(slope_src, tmp_slope, sl_s)
    cfg['slope_file'] = tmp_slope
    cfg['lava_name']  = f"MC_trial_{trial_idx}"

    # Temp JSON
    tmp_json = f'resource/_mc_cfg_{trial_idx}.json'
    with open(tmp_json, 'w') as f:
        json.dump(cfg, f)

    pos, visc, ys = np.array([]), np.array([]), np.array([])
    try:
        fg = run_flowgo.RunFlowgo()
        fg.run(tmp_json, RESULTS_DIR)
        csv_path = fg.get_file_name_results(RESULTS_DIR, tmp_json)
        pos, visc, ys = load_pf_csv(csv_path)
        if os.path.exists(csv_path):
            os.remove(csv_path)
    except Exception as e:
        pass
    finally:
        for fpath in [tmp_json, tmp_slope]:
            if os.path.exists(fpath):
                os.remove(fpath)

    params = {
        'Depth H': d_s,
        'Width w': w_s,
        'Effusion Q': Q_s,
        'Density rho': r_s,
        'Erupt Temp T': T_s,
        'Cryst C_i': ci_s,
        'Slope offset': sl_s
    }
    # Discard very short runs
    if len(pos) == 0 or pos[-1] < 0.05 * nominal_length:
        return None
    return pos, visc, ys, params


# ------------------------------------------------------------------ #
# Main loop                                                           #
# ------------------------------------------------------------------ #
plt.rcParams.update({
    'font.family':     'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':       11,
    'axes.linewidth':  1.2,
    'figure.dpi':      150,
})

for eruption, cfg_meta in ERUPTIONS.items():
    print(f"\n{'='*62}")
    print(f"  MC — {cfg_meta['label']}  (N={N_TRIALS})")
    print(f"{'='*62}")

    with open(cfg_meta['json'], 'r') as f:
        base_cfg = json.load(f)

    slope_src = base_cfg['slope_file']
    color     = cfg_meta['color']

    # --- Nominal run ---
    fg_nom = run_flowgo.RunFlowgo()
    fg_nom.run(cfg_meta['json'], RESULTS_DIR)
    nom_csv  = fg_nom.get_file_name_results(RESULTS_DIR, cfg_meta['json'])
    nom_pos, nom_visc, nom_ys = load_pf_csv(nom_csv)
    nominal_length = float(nom_pos[-1]) if len(nom_pos) > 0 else 2000.0
    print(f"  Nominal run length: {nominal_length:.0f} m")

    grid = np.linspace(0, nominal_length, GRID_POINTS)

    all_visc, all_ys = [], []
    all_params = []
    n_failed = 0

    for i in range(N_TRIALS):
        result = sample_and_run(base_cfg, slope_src, i, nominal_length)
        if result is None:
            n_failed += 1
            continue
        p, v, y, params = result
        all_visc.append(np.interp(grid, p, v, right=float(v[-1])))
        all_ys.append(  np.interp(grid, p, y, right=float(y[-1])))
        all_params.append(params)
        if (i + 1) % 10 == 0:
            print(f"  [{i+1:3d}/{N_TRIALS}]  Success: {i+1-n_failed}  Failed: {n_failed}")

    n_ok = N_TRIALS - n_failed
    print(f"\n  Completed: {n_ok}/{N_TRIALS}  ({n_failed} failed)")

    if n_ok < 5:
        print("  Too few successful runs — skipping plot.")
        continue

    AV  = np.array(all_visc)
    AY  = np.array(all_ys)
    vp2, vp50, vp97 = (np.percentile(AV, p, axis=0) for p in [12.5, 50, 87.5])
    yp2, yp50, yp97 = (np.percentile(AY, p, axis=0) for p in [12.5, 50, 87.5])

    # ---- Plot ----
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 9), sharex=True)

    # Viscosity panel
    ax1.fill_between(grid, vp2, vp97, color=color, alpha=0.20, label='75% CI')
    ax1.plot(grid, vp50, '--', color=color,   lw=1.8, alpha=0.85, label='MC Median')
    ax1.plot(nom_pos, nom_visc, '-', color='#111111', lw=2.2,     label='Nominal run')
    ax1.set_yscale('log')
    ax1.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
    ax1.yaxis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(2,10)*0.1, numticks=15))
    ax1.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax1.set_ylabel(r'Apparent Viscosity ($\eta$, Pa·s)', fontsize=13, fontweight='bold')
    ax1.set_title(f"{cfg_meta['label']} — PyFLOWGO Monte Carlo (N={n_ok})",
                  fontsize=14, fontweight='bold', pad=15)
    ax1.grid(True, which='both', ls='--', alpha=0.4)
    ax1.legend(fontsize=10, loc='lower right', framealpha=0.9)
    ax1.set_xlim(0, nominal_length * 1.02)

    ann = (f"σ(H)=±0.75 m   σ(w)=±4.0 m\n"
           f"σ(Q)=±15%    σ(ρ)=±100 kg/m³\n"
           f"T_erupt: Triangular ±10 K\n"
           f"C_i: Uniform [10–20%]\n"
           f"Slope: Normal offset ±1.25°\n"
           f"{n_ok} successful / {N_TRIALS} total")
    ax1.text(0.02, 0.98, ann, transform=ax1.transAxes, fontsize=8.5,
             va='top', ha='left', family='monospace',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                       edgecolor=color, lw=1.1, alpha=0.93))

    # Yield strength panel
    ax2.fill_between(grid, yp2, yp97, color=color, alpha=0.20, label='75% CI')
    ax2.plot(grid, yp50, '--', color=color,   lw=1.8, alpha=0.85, label='MC Median')
    ax2.plot(nom_pos, nom_ys, '-', color='#111111', lw=2.2,       label='Nominal run')
    ax2.set_yscale('log')
    ax2.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
    ax2.yaxis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(2,10)*0.1, numticks=15))
    ax2.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax2.set_xlabel('Distance from Vent (m)', fontsize=13, fontweight='bold')
    ax2.set_ylabel(r'Yield Strength ($\tau_y$, Pa)', fontsize=13, fontweight='bold')
    ax2.grid(True, which='both', ls='--', alpha=0.4)
    ax2.legend(fontsize=10, loc='lower right', framealpha=0.9)

    plt.tight_layout()
    plt.savefig(cfg_meta['out'], dpi=300)
    plt.close()
    print(f"  Saved: {cfg_meta['out']}")

    # ---- Spearman Rank Correlation ----
    final_visc = AV[:, -1]
    final_ys = AY[:, -1]
    
    param_names = list(all_params[0].keys())
    spearman_visc = []
    spearman_ys = []
    
    for p_name in param_names:
        p_vals = [prm[p_name] for prm in all_params]
        r_v, _ = stats.spearmanr(p_vals, final_visc)
        r_y, _ = stats.spearmanr(p_vals, final_ys)
        spearman_visc.append(r_v)
        spearman_ys.append(r_y)
        
    fig_s, (ax_sv, ax_sy) = plt.subplots(1, 2, figsize=(12, 5))
    
    y_pos = np.arange(len(param_names))
    
    # Sort for Viscosity
    sorted_idx_v = np.argsort(spearman_visc)
    sorted_params_v = [param_names[i] for i in sorted_idx_v]
    sorted_corr_v = [spearman_visc[i] for i in sorted_idx_v]
    
    ax_sv.barh(y_pos, sorted_corr_v, color=color, alpha=0.7)
    ax_sv.set_yticks(y_pos)
    ax_sv.set_yticklabels(sorted_params_v)
    ax_sv.set_xlim(-1, 1)
    ax_sv.grid(True, ls='--', alpha=0.4)
    ax_sv.set_title(f"Spearman's r with Final Viscosity")
    ax_sv.set_xlabel("Spearman's rank correlation coefficient")
    
    # Sort for Yield Strength
    sorted_idx_y = np.argsort(spearman_ys)
    sorted_params_y = [param_names[i] for i in sorted_idx_y]
    sorted_corr_y = [spearman_ys[i] for i in sorted_idx_y]
    
    ax_sy.barh(y_pos, sorted_corr_y, color=color, alpha=0.7)
    ax_sy.set_yticks(y_pos)
    ax_sy.set_yticklabels(sorted_params_y)
    ax_sy.set_xlim(-1, 1)
    ax_sy.grid(True, ls='--', alpha=0.4)
    ax_sy.set_title(f"Spearman's r with Final Yield Strength")
    ax_sy.set_xlabel("Spearman's rank correlation coefficient")
    
    plt.suptitle(f"{cfg_meta['label']} — Sensitivity Analysis (N={n_ok})", fontsize=14, fontweight='bold')
    plt.tight_layout()
    spearman_out = cfg_meta['out'].replace('.png', '_spearman.png')
    plt.savefig(spearman_out, dpi=300)
    plt.close()
    print(f"  Saved Spearman plot: {spearman_out}")

print("\nAll done.")
