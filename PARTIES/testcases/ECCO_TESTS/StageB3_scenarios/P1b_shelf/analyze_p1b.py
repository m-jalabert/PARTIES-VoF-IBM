#!/usr/bin/env python3
"""Post-run analysis of B.3 P1b (job 20977443) against its matched control.

Inputs (scratch): P1b default-operator run, the matched default-operator
control and the Yang-operator control.  Uses ecco_profiles.csv (every 10 steps),
mobile.dat, budgets.json and mid-plane slices of a few Data frames.  Never
submits jobs.  Writes figures/CSVs/summary.json to --out.

Conventions: z = y_int - y is the depth below the initial ice interface (d);
the tracer integral is a meltwater-equivalent volume (d^3); the grain volume
is pi/6 d^3.  Run with the Anvil anaconda 2025.06 python.
"""
import argparse
import csv
import glob
import json
import re
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path('/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf')
RUNS = {'P1b': ROOT / 'run_prod_20977443', 'ctrl': ROOT / 'run_ctrl_20977444',
        'ctrl_yang': ROOT / 'run_ctrl_20973341'}
Y_INT, VG, WS_UREF = 20.0, np.pi / 6, None

# dataviz reference palette (light): categorical slots 1-3, ordinal blue ramp,
# blue<->red diverging with a gray midpoint, text/surface tokens.
SURF, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
CAT = {'P1b': '#2a78d6', 'ctrl': '#eb6834', 'ctrl_yang': '#1baf7a'}
LABEL = {'P1b': 'P1b (grain, default operator)', 'ctrl': 'control (no grain, default operator)',
         'ctrl_yang': 'control (no grain, Yang operator)'}
ORD = ['#86b6ef', '#3987e5', '#256abf', '#184f95', '#0d366b']
SEQ = LinearSegmentedColormap.from_list('seq', ['#fcfcfb', '#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b'])
DIV = LinearSegmentedColormap.from_list('div', ['#104281', '#3987e5', '#b7d3f6', '#f0efec',
                                                '#f4b8b5', '#e34948', '#8f1f1f'])


def style():
    plt.rcParams.update({'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF,
                         'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                         'ytick.color': INK2, 'text.color': INK, 'axes.grid': True, 'grid.color': GRID,
                         'grid.linewidth': 0.6, 'axes.spines.top': False, 'axes.spines.right': False,
                         'font.size': 10, 'lines.linewidth': 2, 'legend.frameon': False})


def profiles(run):
    """Return times, y and per-(time,row) arrays from ecco_profiles.csv (complete blocks only)."""
    a = np.genfromtxt(run / 'ecco_profiles.csv', delimiter=',', names=True)
    y = np.unique(a['y'])
    ny = len(y)
    times, blocks = [], []
    for t in np.unique(a['time']):
        b = a[a['time'] == t][-ny:]
        if len(b) == ny:
            times.append(t)
            blocks.append(b[np.argsort(b['y'])])
    stack = {k: np.array([b[k] for b in blocks]) for k in a.dtype.names}
    return np.array(times), y, stack


def at_times(t_src, field, t_new):
    """Linear interpolation in time of a (time, row) array."""
    out = np.empty((len(t_new), field.shape[1]))
    for j in range(field.shape[1]):
        out[:, j] = np.interp(t_new, t_src, field[:, j])
    return out


def frame(run, t):
    files = sorted(glob.glob(str(run / 'Data_*.h5')), key=lambda p: int(re.findall(r'_(\d+)\.h5', p)[0]))
    return min(files, key=lambda p: abs(float(h5py.File(p)['time'][0]) - t))


def midplane(path, name, k):
    with h5py.File(path) as f:
        return float(f['time'][0]), f[name][k, :-1, :-1], f['grid/xc'][...].ravel()[:-1], f['grid/yc'][...].ravel()[:-1]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    out = a.out
    out.mkdir(parents=True, exist_ok=True)
    style()
    case = json.loads((RUNS['P1b'] / 'case.json').read_text())
    ws = case['scaling']['w_s_over_Uref']
    st = case['scaling']['St']
    tref = case['scaling']['t_ref_s']
    summary = dict(runs={k: str(v) for k, v in RUNS.items()}, w_s_SN=ws, t_ref_s=tref)

    # ------------------------------------------------------------- grain
    m = np.loadtxt(RUNS['P1b'] / 'mobile.dat', delimiter=',')
    t_g, y_g, v_g = m[:, 0], m[:, 3], m[:, 6]
    z_g = Y_INT - y_g
    bottom_gap = y_g - 0.5
    plateau = (t_g > 4) & (bottom_gap > 3)
    v_term = float(np.mean(v_g[plateau]))
    land = int(np.argmax(bottom_gap < 0.05)) if np.any(bottom_gap < 0.05) else None
    summary['grain'] = dict(
        terminal_velocity=v_term, terminal_over_ws_SN=-v_term / ws,
        terminal_m_s=-v_term * case['scaling']['U_ref_m_s'],
        lateral_drift_max=float(np.max(np.hypot(m[:, 2] - 4, m[:, 4] - 4))),
        time_to_0p9_terminal=float(t_g[np.argmax(-v_g >= 0.9 * -v_term)]),
        depth_at_0p9_terminal=float(z_g[np.argmax(-v_g >= 0.9 * -v_term)]),
        contact_time=float(t_g[land]) if land is not None else None,
        final=dict(t=float(t_g[-1]), y=float(y_g[-1]), v=float(v_g[-1])),
        min_bottom_gap=float(bottom_gap.min()),
        speed_at_depth={f'{zz}d': float(-np.interp(zz, z_g, v_g) / -v_term) for zz in (1, 2, 3, 4, 6, 8)},
        start_depth=float(z_g[0]),
        decel_onset_bottom_gap=float(bottom_gap[np.argmax((t_g > 5) & (-v_g < 0.95 * -v_term))]),
        gap_at_max_deceleration=float(bottom_gap[1:][np.argmin(np.where(t_g[1:] > 9, np.diff(-v_g) / np.diff(t_g), np.inf))]),
        note='Approach to terminal mixes ceiling (no-slip ice) hindrance and stratification; not separable without a no-ice control.')
    np.savetxt(out / 'grain_trajectory.csv', np.c_[t_g, y_g, z_g, v_g, -v_g / -v_term],
               delimiter=',', header='t,y,z_below_ice,v,v_over_terminal', comments='')

    # ------------------------------------------------------------- profiles
    P = {k: profiles(v) for k, v in RUNS.items()}
    tP, y, sP = P['P1b']
    tC, _, sC = P['ctrl']
    z = Y_INT - y
    tcom = tP[(tP >= tC[0]) & (tP <= tC[-1])]
    trP = at_times(tP, sP['tracer_integral'], tcom)
    trC = at_times(tC, sC['tracer_integral'], tcom)
    dtr = trP - trC                                  # meltwater-equivalent excess per row
    depths = [2, 4, 6, 8, 12]
    E = {zz: dtr[:, z > zz].sum(axis=1) for zz in depths}
    np.savetxt(out / 'meltwater_excess_vs_time.csv',
               np.c_[tcom, *[E[zz] for zz in depths]], delimiter=',',
               header='t,' + ','.join(f'excess_below_{zz}d' for zz in depths), comments='')
    i_end = len(tcom) - 1
    summary['meltwater'] = dict(
        common_final_time=float(tcom[-1]),
        excess_below_depth_d3={f'{zz}d': float(E[zz][-1]) for zz in depths},
        excess_below_depth_grain_volumes={f'{zz}d': float(E[zz][-1] / VG) for zz in depths},
        peak_excess_below_depth_grain_volumes={f'{zz}d': float(E[zz].max() / VG) for zz in depths},
        deficit_top_2d_d3=float(dtr[i_end, z <= 2].sum()),
        total_change_d3=float(dtr[i_end].sum()),
        centroid_depth_of_excess_d=float(np.sum(np.clip(dtr[i_end], 0, None) * z) / np.sum(np.clip(dtr[i_end], 0, None))))
    np.savetxt(out / 'meltwater_excess_profile_final.csv', np.c_[z, dtr[i_end]], delimiter=',',
               header='z_below_ice,excess_tracer_per_row_d3', comments='')

    # Liquid-weighted mean S and theta, P1b - control, at the final common time
    def meanfield(s, k, t):
        return np.array([np.interp(t, s[0], s[2][k][:, j]) for j in range(len(y))])
    dS = meanfield(P['P1b'], 'mean_S', tcom[-1]) - meanfield(P['ctrl'], 'mean_S', tcom[-1])
    dT = meanfield(P['P1b'], 'mean_T', tcom[-1]) - meanfield(P['ctrl'], 'mean_T', tcom[-1])
    lq = z > 0.3                                     # liquid, below the diffuse band
    summary['stratification_change_final'] = dict(
        max_abs_dS_mean=float(np.nanmax(np.abs(dS[lq]))), depth_of_max_dS=float(z[lq][np.nanargmax(np.abs(dS[lq]))]),
        max_abs_dT_mean=float(np.nanmax(np.abs(dT[lq]))), depth_of_max_dT=float(z[lq][np.nanargmax(np.abs(dT[lq]))]))
    np.savetxt(out / 'mean_profile_change_final.csv', np.c_[z, dS, dT], delimiter=',',
               header='z_below_ice,dS_mean_P1b_minus_ctrl,dtheta_mean_P1b_minus_ctrl', comments='')

    # Dissipation (integrated over the domain) and its time integral
    epsP = at_times(tP, sP['epsilon_integral'], tcom).sum(axis=1)
    epsC = at_times(tC, sC['epsilon_integral'], tcom).sum(axis=1)
    summary['dissipation'] = dict(time_integral_P1b=float(np.trapezoid(epsP, tcom)),
                                  time_integral_ctrl=float(np.trapezoid(epsC, tcom)),
                                  peak_rate_P1b=float(epsP.max()), t_peak=float(tcom[np.argmax(epsP)]))

    # ------------------------------------------------------------- melt / ice
    melt = {}
    for k, (tt, yy, ss) in P.items():
        ice = ss['ice_volume']
        melt[k] = (tt, ice[0].sum() - ice.sum(axis=1), ice[:, (Y_INT - yy) > 0.5].sum(axis=1))
    summary['melt'] = {k: dict(final_t=float(v[0][-1]), net_melt_d3=float(v[1][-1]),
                               ice_more_than_0p5d_below_interface_max_d3=float(v[2].max()))
                       for k, v in melt.items()}
    np.savetxt(out / 'net_melt_series_P1b.csv', np.c_[melt['P1b'][0], melt['P1b'][1], melt['P1b'][2]],
               delimiter=',', header='t,net_melt_d3,ice_below_0p5d_d3', comments='')

    # ------------------------------------------------------------- operator comparison (no grain)
    tY, _, sY = P['ctrl_yang']
    op = {}
    for tt in (1, 3, 7, 12, 15):
        if tt > min(tC[-1], tY[-1]):
            continue
        dth = meanfield(P['ctrl'], 'mean_T', tt) - meanfield(P['ctrl_yang'], 'mean_T', tt)
        dsl = meanfield(P['ctrl'], 'mean_S', tt) - meanfield(P['ctrl_yang'], 'mean_S', tt)
        liq = z > 0.3
        reach = float(z[liq][np.abs(dth[liq]) > 0.01].max()) if np.any(np.abs(dth[liq]) > 0.01) else 0.0
        reach_s = float(z[liq][np.abs(dsl[liq]) > 0.01].max()) if np.any(np.abs(dsl[liq]) > 0.01) else 0.0
        op[tt] = (dth, dsl)
        summary.setdefault('operator_comparison', {})[f't{tt}'] = dict(
            deepest_depth_dtheta_gt_0p01_d=reach, deepest_depth_dS_gt_0p01_d=reach_s, dtheta_at_0p5d=float(np.interp(0.5, z[::-1], dth[::-1])),
            dS_at_0p5d=float(np.interp(0.5, z[::-1], dsl[::-1])))

    # ------------------------------------------------------------- budgets
    summary['budgets'] = {}
    for k, run in RUNS.items():
        b = json.loads((run / 'budgets.json').read_text()) if (run / 'budgets.json').exists() else {}
        tt, yy, ss = P[k]
        Sint = ss['S_integral'].sum(axis=1)
        SF = np.nansum(ss['mean_S'] * ss['liquid_volume'], axis=1)
        H = ss['T_integral'].sum(axis=1) - ss['ice_volume'].sum(axis=1) / st
        mt = ss['ice_volume'][0].sum() - ss['ice_volume'].sum(axis=1)
        trc = ss['tracer_integral'].sum(axis=1) - ss['tracer_integral'][0].sum()
        summary['budgets'][k] = dict(
            analyze_profiles=b, salt_int_s_rel_drift_max=float(np.max(np.abs(Sint / Sint[0] - 1))),
            salt_int_FS_rel_drift_final=float(SF[-1] / SF[0] - 1),
            conserved_salt=('int (F+delta) S (Yang)' if k == 'ctrl_yang' else 'int s (default operator)'),
            heat_drift_max_abs=float(np.max(np.abs(H - H[0]))),
            tracer_minus_melt_max_abs=float(np.max(np.abs(trc - mt))))

    # ------------------------------------------------------------- figures
    # F1 grain
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(t_g, z_g, color=CAT['P1b'])
    ax[0].invert_yaxis(); ax[0].set_xlabel('t (code units)'); ax[0].set_ylabel('grain centre depth below ice (d)')
    ax[0].axhspan(0, 2 * case['three_equation']['delta_S_d'], color='#e9e8e3', zorder=0)
    ax[0].text(t_g[-1] * 0.98, 2.2, 'salt sublayer (2δ_S)', ha='right', color=INK2, fontsize=9)
    ax[0].set_title('a  Grain depth', loc='left')
    sel = bottom_gap > -1
    ax[1].plot(z_g[sel], -v_g[sel] / -v_term, color=CAT['P1b'])
    ax[1].axhline(1, color=INK2, lw=1, ls='--'); ax[1].text(0.3, 1.02, f'terminal = {-v_term:.3f} U_ref = {-v_term/ws:.3f} w_s(SN)', color=INK2, fontsize=9)
    ax[1].axvspan(0, 2 * case['three_equation']['delta_S_d'], color='#e9e8e3', zorder=0)
    ax[1].set_xlabel('grain centre depth below ice (d)'); ax[1].set_ylabel('settling speed / terminal')
    ax[1].set_title('b  Settling speed through the meltwater layer', loc='left')
    fig.tight_layout(); fig.savefig(out / 'F1_grain.png', dpi=160); plt.close(fig)

    # F2 meltwater excess
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for c, zz in zip(ORD, depths):
        ax[0].plot(tcom, E[zz] / VG, color=c, label=f'deeper than {zz} d')
        ax[0].text(tcom[-1], E[zz][-1] / VG, f' {zz} d', color=INK2, va='center', fontsize=9)
    ax[0].set_xlabel('t (code units)'); ax[0].set_ylabel('meltwater carried below depth (grain volumes)')
    ax[0].set_title('a  Meltwater excess, P1b − control', loc='left'); ax[0].legend(loc='upper left', fontsize=8)
    ax[1].axvline(0, color=INK2, lw=1)
    zc = summary.get('operator_comparison', {}).get('t15', {}).get('deepest_depth_dtheta_gt_0p01_d')
    if zc:
        ax[1].axhspan(-0.3, zc, color='#e9e8e3', zorder=0)
        ax[1].text(0.98, zc - 0.2, 'band-kernel affected (t ≈ 15)', transform=ax[1].get_yaxis_transform(),
                   ha='right', va='bottom', color=INK2, fontsize=8)
    ax[1].plot(dtr[i_end] / VG * 24, z, color=CAT['P1b'])
    ax[1].invert_yaxis(); ax[1].set_xlabel('excess per unit depth (grain volumes per d)')
    ax[1].set_ylabel('depth below ice (d)'); ax[1].set_title(f'b  Excess profile at t = {tcom[-1]:.1f}', loc='left')
    fig.tight_layout(); fig.savefig(out / 'F2_meltwater_excess.png', dpi=160); plt.close(fig)

    # F3 mid-plane slices (tracer and P1b - control), equal aspect so the grain stays round
    times = [t for t in (1.0, 2.5, 5.0, 10.0, 18.0) if t <= tcom[-1] + 1e-6]
    fig, axes = plt.subplots(2, len(times), figsize=(1.9 * len(times) + 1.6, 10.2), sharex=True, sharey=True,
                             constrained_layout=True)
    k_mid = 96
    vmax_d = 0.0
    slices = []
    for t in times:
        fp = frame(RUNS['P1b'], t)
        tp, cP, x, yy = midplane(fp, 'Conc/2', k_mid)
        _, csP, _, _ = midplane(fp, 'VOF/C_S', k_mid)
        tc, cC, _, _ = midplane(frame(RUNS['ctrl'], t), 'Conc/2', k_mid)
        solid = csP >= 0.05                          # sediment support: extension values, not fluid
        cP, dC = np.where(solid, np.nan, cP), np.where(solid, np.nan, cP - cC)
        slices.append((tp, tc, cP, dC, x, yy))
        vmax_d = max(vmax_d, np.nanmax(np.abs(dC)))
    ring = np.linspace(0, 2 * np.pi, 60)
    for i, (tp, tc, cP, dC, x, yy) in enumerate(slices):
        im0 = axes[0, i].pcolormesh(x, Y_INT - yy, cP, cmap=SEQ, vmin=0, vmax=0.4, shading='auto')
        im1 = axes[1, i].pcolormesh(x, Y_INT - yy, dC, cmap=DIV, vmin=-vmax_d, vmax=vmax_d, shading='auto')
        tg = np.interp(tp, t_g, z_g)
        for r in (0, 1):
            axes[r, i].set_aspect('equal')
            axes[r, i].grid(False)
            axes[r, i].plot(4 + 0.5 * np.cos(ring), tg + 0.5 * np.sin(ring), color=INK, lw=1)
        axes[0, i].set_title(f't = {tp:.1f}', loc='left')
        axes[1, i].set_xlabel('x (d)')
    axes[0, 0].set_ylim(20.2, -0.6)
    axes[0, 0].set_ylabel('meltwater tracer, P1b\ndepth below ice (d)')
    axes[1, 0].set_ylabel('P1b minus control\ndepth below ice (d)')
    fig.colorbar(im0, ax=axes[0, :].tolist(), shrink=0.7, label='tracer (mid-plane, z = 4 d)')
    fig.colorbar(im1, ax=axes[1, :].tolist(), shrink=0.7, label='tracer difference')
    fig.savefig(out / 'F3_tracer_slices.png', dpi=150); plt.close(fig)

    # F4 operator comparison
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for c, (tt, (dth, dsl)) in zip(ORD, op.items()):
        ax[0].plot(dth, z, color=c, label=f't = {tt}')
        ax[1].plot(dsl, z, color=c, label=f't = {tt}')
    for a_ in ax:
        a_.axvline(0, color=INK2, lw=1); a_.set_ylim(4, -0.3)
        a_.axhspan(-0.3, 0.3, color='#e9e8e3', zorder=0)
    ax[1].text(0.02, 0.32, 'diffuse band: salt definitions differ (volume-averaged vs liquid)', transform=ax[1].get_yaxis_transform(),
               va='top', color=INK2, fontsize=8)
    ax[0].set_xlabel('Δθ (default − Yang), liquid-weighted mean'); ax[0].set_ylabel('depth below ice (d)')
    ax[1].set_xlabel('Δ mean S (default − Yang)')
    ax[0].set_title('a  Temperature', loc='left'); ax[1].set_title('b  Salinity', loc='left')
    ax[0].legend(fontsize=8, loc='lower right')
    fig.tight_layout(); fig.savefig(out / 'F4_operator_comparison.png', dpi=160); plt.close(fig)

    # F5 melt
    fig, ax = plt.subplots(figsize=(6, 4))
    for k, (tt, mlt, _) in melt.items():
        ax.plot(tt, mlt, color=CAT[k], label=LABEL[k])
        ax.text(tt[-1], mlt[-1], f' {k}', color=INK2, va='center', fontsize=9)
    ax.axhline(0, color=INK2, lw=1)
    ax.set_xlabel('t (code units)'); ax.set_ylabel('net melt (d³), negative = net freezing')
    ax.set_title('Net melt: band-kernel sign artifact', loc='left'); ax.legend(fontsize=8, loc='lower left')
    fig.tight_layout(); fig.savefig(out / 'F5_net_melt.png', dpi=160); plt.close(fig)

    (out / 'summary.json').write_text(json.dumps(summary, indent=2, default=float) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'budgets'}, indent=1, default=float))
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != 'analyze_profiles'} for k, v in summary['budgets'].items()}, indent=1))


if __name__ == '__main__':
    main()
