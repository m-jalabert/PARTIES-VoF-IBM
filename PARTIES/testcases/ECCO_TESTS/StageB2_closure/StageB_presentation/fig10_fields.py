#!/usr/bin/env python3
"""Figure 10 — what the simulation actually looks like: mid-plane montage.

Four saved snapshots of the accepted B.2 coupled gate, on the vertical plane
through the grain centre (z = Lz/2).

  top row     temperature theta — what drives the melt and the convection
  bottom row  meltwater anomaly C_mw - <C_mw>_liquid(t) — where the melt goes

The box is closed, so the mean meltwater concentration rises monotonically;
subtracting the instantaneous liquid mean is what makes the plume structure
visible at every time on one fixed colour scale.  Tracer values inside the IBM's
excluded support are extension values, not physical concentrations, and are
masked out rather than plotted.
"""
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stageb_common as C

plt = C.use_style()
CM = C.cmaps()

RUN = C.run('coupled16')
rel = C.release_time(RUN)
snaps = dict(C.snapshots(RUN))

TIMES = [10.0, 34.0, 40.0, 64.0]
CAPTIONS = ['ice melting from below,\ngrain locked inside it',
            f'release, t = {rel:.2f}',
            'falling through\nits own meltwater',
            'at rest on the wall,\nplume still rising']


def nearest(t):
    best = None
    for idx, path in C.snapshots(RUN):
        with h5py.File(path.with_name(path.name.replace('Data_', 'Particle_')), 'r') as p:
            tt = float(p['time'][0])
        if best is None or abs(tt - t) < abs(best[1] - t):
            best = (idx, tt)
    return best


fig = plt.figure(figsize=(13.4, 9.0))
gs = fig.add_gridspec(2, 4, hspace=0.12, wspace=0.10,
                      left=0.048, right=0.895, top=0.845, bottom=0.055)

TMAX = 0.08   # meltwater anomaly range, fixed across all panels
im_T = im_M = None

for k, tt in enumerate(TIMES):
    idx, t_actual = nearest(tt)
    f = C.load_fields(snaps[idx])
    X, U, R, _ = C.particle_state(str(snaps[idx]).replace('Data_', 'Particle_'))
    x, y, z = f['x'], f['y'], f['z']
    jz = int(np.argmin(np.abs(z - X[2])))
    ext = [x[0], x[-1], y[0], y[-1]]

    # liquid-volume mean of the tracer at this time, over the physical mesh
    liquid = (f['VOF/C_S'] < 0.05) & (f['ice'] < 0.5)
    mw_mean = float(f['Conc/2'][liquid].mean())

    # field arrays are [z, y, x]; a z-slice is [y, x], which is what imshow wants
    cs, ice = f['VOF/C_S'][jz], f['ice'][jz]
    solid = (cs >= 0.05) | (ice > 0.5)
    panels = [
        (np.ma.masked_where(solid, f['Conc/0'][jz]), CM['orange'], 0.0, 1.0),
        (np.ma.masked_where(solid, f['Conc/2'][jz] - mw_mean),
         CM['diverging'], -TMAX, TMAX),
    ]

    for row, (field, cmap, lo, hi) in enumerate(panels):
        ax = fig.add_subplot(gs[row, k])
        ax.set_facecolor('#eceae5')
        ax.grid(False)
        h = ax.imshow(field, origin='lower', extent=ext, aspect='equal',
                      cmap=cmap, vmin=lo, vmax=hi, interpolation='bilinear')
        if row == 0:
            im_T = h
        else:
            im_M = h
        # ice block and its melt front
        ax.imshow(np.ma.masked_where(ice <= 0.5, np.ones_like(ice)),
                  origin='lower', extent=ext, aspect='equal', zorder=3,
                  cmap=plt.matplotlib.colors.ListedColormap(['#e3edf7']),
                  vmin=0, vmax=1, interpolation='nearest')
        ax.contour(x, y, ice, levels=[0.5], colors=[C.AQUA], linewidths=1.2, zorder=4)
        # grain, and the IBM support that is excluded from the scalar fields
        ax.contourf(x, y, cs, levels=[0.5, 1.01], colors=['#4f4840'], zorder=5)
        ax.contour(x, y, cs, levels=[0.05], colors=[C.ORANGE], linewidths=0.9,
                   linestyles=[(0, (3, 2))], zorder=6)
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(y[0], y[-1])
        ax.tick_params(labelsize=8)
        ax.set_xticks([0, 1, 2, 3, 4])
        ax.set_yticks([0, 2, 4, 6])
        for s in ax.spines.values():
            s.set_color(C.AXIS)
        if k:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel('$y$')
        if row == 0:
            ax.set_xticklabels([])
            ax.set_title(f'$t$ = {t_actual:.1f}', fontsize=11.5)
            ax.text(0.5, -0.045, CAPTIONS[k], transform=ax.transAxes, fontsize=8.5,
                    color=C.INK2, ha='center', va='top')
        else:
            ax.set_xlabel('$x$')

fig.text(0.012, 0.655, 'temperature  θ', fontsize=10.5, color=C.INK,
         fontweight='semibold', rotation=90, va='center', ha='left')
fig.text(0.012, 0.235, 'meltwater anomaly', fontsize=10.5, color=C.INK,
         fontweight='semibold', rotation=90, va='center', ha='left')

cax1 = fig.add_axes([0.908, 0.505, 0.013, 0.30])
cb1 = fig.colorbar(im_T, cax=cax1)
cb1.set_label('θ   (0 = ice, 1 = warm wall)', color=C.INK2, fontsize=9)
cax2 = fig.add_axes([0.908, 0.085, 0.013, 0.30])
cb2 = fig.colorbar(im_M, cax=cax2, extend='both')
cb2.set_label(r'$C_{mw} - \langle C_{mw}\rangle_{liquid}(t)$', color=C.INK2, fontsize=9)
for cb in (cb1, cb2):
    cb.outline.set_edgecolor(C.AXIS)
    cb.ax.tick_params(labelsize=8, colors=C.MUTED)

fig.text(0.908, 0.44, 'melt front', color=C.AQUA, fontsize=9)
fig.text(0.908, 0.420, 'grain', color='#4f4840', fontsize=9)
fig.text(0.908, 0.400, 'IBM support', color=C.ORANGE, fontsize=9)

fig.suptitle('Stage B · the melt, the release and the plume, on the plane through the grain',
             x=0.048, ha='left', fontsize=15, fontweight='semibold', color=C.INK, y=0.968)
fig.text(0.048, 0.905,
         'Accepted B.2 coupled gate, job 20492853 — vertical slice at z = 2.0 of a '
         '4 × 6 × 4 box; grain diameter d = 1 at d/Δx = 16; cold top wall, warm bottom wall',
         ha='left', fontsize=9.5, color=C.INK2)

C.save(fig, 'fig10_fields',
       'Source: B2closure_coupled16_rampfix_20492853 Data_*.h5 / Particle_*.h5. '
       'Synthetic debugging parameters (Pe = {10, 70, 10}, St = 0.25, passive salt), '
       'not a physical ECCO scenario.')
