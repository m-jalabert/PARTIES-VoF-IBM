#!/usr/bin/env python3
"""Figure 11 — where the meltwater goes, and one thing to check before B.3.

Panels (a)-(c) come from ``ecco_profiles.csv``, the in-solver ECCO_PROFILES
output (default-off; written every ten steps and at every field output).  Its
horizontal means and covariance fluxes are liquid-weighted; the conserved
scalar integrals it also carries are over the whole physical mesh.

Panel (d) is an observation made while preparing these figures, not a recorded
gate: with the hard ice-penalization threshold active, downward ice tongues
grow after release and melt away again; the threshold-disabled twin, run with
the same deck to the same time, never grows any.  It is measured with the
project's own ``real_ice`` definition (ice masked to C_S < 0.05) so it is not
the IBM support being miscounted.  The scan is cached in data/.

    python fig11_mixing.py --rescan    # redo the field scan for panel (d)
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stageb_common as C

plt = C.use_style()
CM = C.cmaps()

ap = argparse.ArgumentParser()
ap.add_argument('--rescan', action='store_true')
args = ap.parse_args()

RUN = C.run('coupled16')
rel = C.release_time(RUN)
tm, Xm, _ = C.mobile(RUN)
wall_hit = float(tm[np.argmax(Xm[:, 1] < 0.51)])
times, yy, F = C.profile_blocks(RUN)

CACHE = C.DATADIR / 'low_ice_scan.json'


def scan_low_ice(key, y_cut=3.0):
    """Volume of ACTUAL ice (C_S < 0.05) below y_cut, per saved output."""
    out = []
    for _, path in C.snapshots(C.run(key)):
        f = C.load_fields(path, names=('VOF/C_L', 'VOF/C_S'))
        y = f['y']
        dv = (f['x'][1] - f['x'][0]) * (y[1] - y[0]) * (f['z'][1] - f['z'][0])
        real = np.where(f['VOF/C_S'] < 0.05, f['ice'], 0.0)[:, y < y_cut, :]
        m = real > 0.9
        out.append([f['t'], float(real[m].sum() * dv) if m.any() else 0.0,
                    int(m.sum())])
    return out


if args.rescan or not CACHE.exists():
    CACHE.write_text(json.dumps({k: scan_low_ice(k) for k in
                                 ('coupled16', 'coupled16_smooth')}, indent=1) + '\n')
low = {k: np.array(v) for k, v in json.loads(CACHE.read_text()).items()}

fig = plt.figure(figsize=(12.6, 7.4))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.24,
                      left=0.068, right=0.975, top=0.815, bottom=0.095)


def overlay(ax):
    ax.plot(tm, Xm[:, 1], color=C.INK, lw=1.6, label='grain centre')
    ax.axvline(rel, color=C.STATUS['critical'], lw=1.2, ls='--')
    ax.text(rel - 0.8, 5.8, 'release', color=C.STATUS['critical'], fontsize=8,
            ha='right', va='top')
    ax.set_xlim(0, times.max())
    ax.set_ylim(0, 6)
    ax.set_xlabel('time  $t$')
    ax.set_ylabel('height  $y$')
    ax.grid(False)


# --- (a) horizontally averaged meltwater --------------------------------
ax = C.finish(fig.add_subplot(gs[0, 0]))
mw = F['mean_tracer']
pm = ax.pcolormesh(times, yy, mw.T, cmap=CM['blue'], shading='auto',
                   vmin=0, vmax=np.nanmax(mw))
cb = fig.colorbar(pm, ax=ax, pad=0.02)
cb.set_label(r'$\langle C_{mw}\rangle$', fontsize=9, color=C.INK2)
cb.ax.tick_params(labelsize=8, colors=C.MUTED)
cb.outline.set_edgecolor(C.AXIS)
overlay(ax)
ax.set_title('a   Meltwater fills the box from the ice down')
ax.legend(loc='lower left')

# --- (b) vertical meltwater flux ----------------------------------------
ax = C.finish(fig.add_subplot(gs[0, 1]))
fl = F['flux_vtracer']
lim = float(np.nanpercentile(np.abs(fl), 99.5))
pm = ax.pcolormesh(times, yy, fl.T, cmap=CM['diverging'], shading='auto',
                   vmin=-lim, vmax=lim)
cb = fig.colorbar(pm, ax=ax, pad=0.02, extend='both')
cb.set_label(r"liquid-weighted  $\langle v' C_{mw}'\rangle$", fontsize=9, color=C.INK2)
cb.ax.tick_params(labelsize=8, colors=C.MUTED)
cb.outline.set_edgecolor(C.AXIS)
overlay(ax)
ax.set_title('b   The grain drags meltwater downward as it falls')
ax.legend(loc='lower left')

# --- (c) integrated diagnostics -----------------------------------------
# two stacked panels rather than one twin axis: the two quantities have
# different units and a shared y-scale would be meaningless.
sub = gs[1, 0].subgridspec(2, 1, hspace=0.12, height_ratios=[1, 1])
axA = C.finish(fig.add_subplot(sub[0]))
axB = C.finish(fig.add_subplot(sub[1]))
ice = np.nansum(F['ice_volume'], axis=1)
eps = np.nansum(F['epsilon_integral'], axis=1)
axA.plot(times, ice, color=C.BLUE, lw=2.0)
axA.set_ylabel('ice volume', fontsize=9)
axA.set_xticklabels([])
axA.set_title('c   Melt rate and dissipation, from the in-solver profiles')
axB.plot(times, eps, color=C.ORANGE, lw=2.0)
axB.set_ylabel(r'dissipation  $\varepsilon$', fontsize=9)
axB.set_xlabel('time  $t$')
for a in (axA, axB):
    a.set_xlim(0, times.max())
    a.axvline(rel, color=C.STATUS['critical'], lw=1.2, ls='--')
axA.text(rel - 0.8, ice.max(), 'release', color=C.STATUS['critical'],
         fontsize=8, ha='right', va='top')
axA.text(1, ice.min() + 0.2, 'the profile reconstruction of all 32 coupled snapshots\n'
         'agrees with the saved fields to 7.0e-15 (scaled)',
         fontsize=8.2, color=C.MUTED, linespacing=1.5, va='bottom')

# --- (d) the threshold observation --------------------------------------
ax = C.finish(fig.add_subplot(gs[1, 1]))
a = low['coupled16']
b = low['coupled16_smooth']
ax.fill_between(a[:, 0], 0, a[:, 1], color='#fdf1e9', lw=0)
ax.plot(a[:, 0], a[:, 1], 'o-', ms=4, color=C.STATUS['warning'],
        label='hard ice-penalization threshold ACTIVE  (20492853)')
ax.plot(b[:, 0], b[:, 1], 's-', ms=4, color=C.STATUS['good'],
        label='threshold DISABLED  (20500720)')
ax.axvline(rel, color=C.STATUS['critical'], lw=1.2, ls='--')
ax.text(rel + 0.8, 0.005, 'release', color=C.STATUS['critical'], fontsize=8,
        ha='left', va='bottom')
k = int(np.argmax(a[:, 1]))
ax.annotate(f'peak {a[k,1]:.3f}\n'
            f'({100*a[k,1]/F["ice_volume"][0].sum():.1f} % of the initial ice),\n'
            'then it melts away again',
            xy=(a[k, 0] - 3.2, 0.245), xytext=(a[k, 0] - 16, 0.21),
            fontsize=8.5, color=C.STATUS['warning'], ha='center',
            arrowprops=dict(arrowstyle='-', color=C.STATUS['warning'], lw=1.0))
ax.set_xlim(0, times.max())
ax.set_ylim(-0.005, 0.38)
ax.set_xlabel('time  $t$')
ax.set_ylabel('volume of actual ice below  $y$ = 3')
ax.set_title('d   Flag for B.3 — an observation, not a recorded gate')
ax.legend(loc='upper left', fontsize=8)

fig.suptitle('Stage B · where the meltwater goes, and one thing to check before B.3',
             x=0.068, ha='left', fontsize=15, fontweight='semibold', color=C.INK, y=0.965)
fig.text(0.068, 0.930,
         'Panels a–c: the default-off ECCO_PROFILES output of the accepted coupled gate.\n'
         'No mixing efficiency is inferred — an open melting box needs its source and boundary '
         'BPE accounting, and the actual EOS, before η means anything.',
         ha='left', fontsize=9.5, color=C.INK2, va='top', linespacing=1.5)

C.save(fig, 'fig11_mixing',
       'Sources: ecco_profiles.csv of B2closure_coupled16_rampfix_20492853 (a–c); '
       'field scan of that run and of B2closure_coupled16_smooth_20500720, cached in '
       'data/low_ice_scan.json (d). Panel d was measured for these figures and has not '
       'been through the campaign audit chain.')
