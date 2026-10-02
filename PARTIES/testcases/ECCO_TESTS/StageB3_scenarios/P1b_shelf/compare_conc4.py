#!/usr/bin/env python3
"""Concentration variant: P1b (one grain per 8 x 8 d) vs conc4 (one per 4 x 4 d).

Both use the default operator and build_v2.  The no-grain control is horizontally
uniform, so its tracer is compared per unit area.  Results are per grain:
meltwater-equivalent excess below depth z, in grain volumes.  Writes
conc4_summary.json and F6_concentration.png to --out.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import analyze_p1b as A

ROOT = A.ROOT
RUNS = {'P1b': (ROOT / 'run_prod_20977443', 64.0), 'conc4': (ROOT / 'run_conc4_21006061', 16.0)}
CTRL = (ROOT / 'run_ctrl_20977444', 64.0)
COL = {'P1b': A.CAT['P1b'], 'conc4': '#eda100'}       # categorical slots 1 and 4
LAB = {'P1b': 'one grain per 8×8 d (P1b)', 'conc4': 'one grain per 4×4 d (conc4)'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    A.style()
    tC, y, sC = A.profiles(CTRL[0])
    z = A.Y_INT - y
    res, series = {}, {}
    for name, (run, area) in RUNS.items():
        t, _, s = A.profiles(run)
        tc = t[(t >= tC[0]) & (t <= tC[-1])]
        tr = A.at_times(t, s['tracer_integral'], tc)
        trc = A.at_times(tC, sC['tracer_integral'], tc) * area / CTRL[1]
        d = tr - trc
        m = np.loadtxt(run / 'mobile.dat', delimiter=',')
        tg, yg, vg = m[:, 0], m[:, 3], m[:, 6]
        gap = yg - 0.5
        vt = float(np.mean(vg[(tg > 4) & (gap > 3)]))
        E = {zz: d[:, z > zz].sum(axis=1) / A.VG for zz in (4, 6, 8, 12)}
        series[name] = (tc, E, tg, A.Y_INT - yg, vg / vt)
        tcmp = min(tc[-1], 17.5)
        i = np.argmin(abs(tc - tcmp))
        res[name] = dict(area_per_grain_d2=area, debris_volume_fraction_one_grain_sheet=float(np.pi / 6 / area),
                         terminal_velocity=vt, terminal_over_P1b=None,
                         contact_time=float(tg[np.argmax(gap < 0.05)]) if np.any(gap < 0.05) else None,
                         speed_at_depth={f'{zz}d': float(-np.interp(zz, A.Y_INT - yg, vg) / -vt) for zz in (1, 2, 3, 4)},
                         compare_time=float(tc[i]),
                         excess_grain_volumes_at_compare={f'{zz}d': float(E[zz][i]) for zz in E},
                         peak_excess_grain_volumes={f'{zz}d': float(E[zz].max()) for zz in E},
                         final_time=float(tg[-1]))
    for name in res:
        res[name]['terminal_over_P1b'] = res[name]['terminal_velocity'] / res['P1b']['terminal_velocity']

    from matplotlib.lines import Line2D
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for name, (tc, E, tg, zg, vr) in series.items():
        for zz, ls in zip((4, 8, 12), ('-', '--', ':')):
            ax[0].plot(tc, E[zz], color=COL[name], ls=ls)
        ax[0].text(tc[-1], E[4][-1], f' {name}', color=A.INK2, va='center', fontsize=9)
        sel = zg < 6
        ax[1].plot(zg[sel], vr[sel], color=COL[name], label=LAB[name])
    ax[0].set_xlabel('t (code units)'); ax[0].set_ylabel('meltwater below depth, per grain\n(grain volumes)')
    ax[0].set_title('a  Per-grain meltwater excess', loc='left')
    runs = [Line2D([], [], color=COL[n], label=LAB[n]) for n in series]
    deps = [Line2D([], [], color=A.INK2, ls=ls, label=f'below {zz} d') for zz, ls in zip((4, 8, 12), ('-', '--', ':'))]
    ax[0].set_xlim(0, tc[-1] * 1.62)                  # empty band right of the curve ends holds the legends
    leg = ax[0].legend(handles=runs, fontsize=8, loc='upper right')
    ax[0].add_artist(leg)
    ax[0].legend(handles=deps, fontsize=8, loc='lower right')
    ax[1].axhline(1, color=A.INK2, lw=1, ls='--')
    ax[1].set_xlabel('grain centre depth below ice (d)'); ax[1].set_ylabel('speed / own terminal speed')
    ax[1].set_title('b  Approach to terminal speed', loc='left'); ax[1].legend(fontsize=8, loc='lower right')
    fig.tight_layout(); fig.savefig(a.out / 'F6_concentration.png', dpi=160); plt.close(fig)
    (a.out / 'conc4_summary.json').write_text(json.dumps(res, indent=2) + '\n')
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
