#!/usr/bin/env python3
"""Diagnostic — the four "corner" ice tongues in video 1 are ONE structure.

x and z are periodic, so the box drawn from 0 to 4 is an arbitrary window onto a
tiled plane.  The descending ice column sits at (x, z) = (0, 0), which is a
CORNER of that window, so it is sliced into four quarter-columns, one at each
corner.  Rolling the fields by half a period puts it in the middle of the window
and it becomes one object.

Left   the window as video 1 draws it
Right  the same instant, rolled by (Lx/2, Lz/2)
Below  what actually happened, from the corner column itself

    python diag_periodic_shift.py            # t = 44
    python diag_periodic_shift.py --time 40
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stageb_common as C

plt = C.use_style()

ap = argparse.ArgumentParser()
ap.add_argument('--time', type=float, default=44.0)
args = ap.parse_args()

from matplotlib.animation import FFMpegWriter  # noqa: F401,E402  (style parity)
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402
from skimage import measure  # noqa: E402

RUN = C.run('coupled16')
D = C.COUPLED16
ICE_RGB = np.array([0.42, 0.62, 0.84])
GRAIN = '#4f4840'

snap = min(C.snapshots(RUN),
           key=lambda p: abs(C.load_fields(p[1], names=('VOF/C_L', 'VOF/C_S'))['t'] - args.time))
f = C.load_fields(snap[1], names=('VOF/C_L', 'VOF/C_S', 'v'))
X, U, R, _ = C.particle_state(str(snap[1]).replace('Data_', 'Particle_'))
x, y, z = f['x'], f['y'], f['z']
real = np.where(f['VOF/C_S'] < 0.05, f['ice'], 0.0)          # [z, y, x]
nz, ny, nx = real.shape


def draw(ax, field, grain_xz, title):
    ax.set_facecolor(C.PAGE)
    ax.set_axis_off()
    pts = np.array([[0, 0, 0], [D['Lx'], 0, 0], [D['Lx'], D['Lz'], 0], [0, D['Lz'], 0],
                    [0, 0, D['Ly']], [D['Lx'], 0, D['Ly']],
                    [D['Lx'], D['Lz'], D['Ly']], [0, D['Lz'], D['Ly']]], float)
    for a, b in [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
                 (0, 4), (1, 5), (2, 6), (3, 7)]:
        ax.plot(*zip(pts[a], pts[b]), color=C.AXIS, lw=0.8, zorder=0)
    dx, dy, dz = x[1] - x[0], y[1] - y[0], z[1] - z[0]
    v, faces, nrm, _ = measure.marching_cubes(field, 0.5, spacing=(dz, dy, dx),
                                              step_size=1)
    tri = np.stack([v[:, 2] + x[0], v[:, 0] + z[0], v[:, 1] + y[0]], 1)[faces]
    nn = np.stack([nrm[:, 2], nrm[:, 0], nrm[:, 1]], 1)[faces].mean(1)
    nn /= np.linalg.norm(nn, axis=1)[:, None] + 1e-12
    ld = np.array([0.35, -0.55, 0.76])
    ld /= np.linalg.norm(ld)
    diff = np.abs(nn @ ld)
    cols = np.clip(ICE_RGB[None, :] * (0.45 + 0.55 * diff)[:, None]
                   + 0.30 * (diff ** 14)[:, None], 0, 1)
    fc = np.c_[cols, np.full(len(cols), 1.0)]
    ax.add_collection3d(Poly3DCollection(tri, facecolors=fc, edgecolors=fc,
                                         linewidths=0.0, antialiased=False, zorder=4))
    u_ = np.linspace(0, 2 * np.pi, 40)
    v_ = np.linspace(0, np.pi, 20)
    ax.plot_surface(grain_xz[0] + R * np.outer(np.cos(u_), np.sin(v_)),
                    grain_xz[1] + R * np.outer(np.sin(u_), np.sin(v_)),
                    X[1] + R * np.outer(np.ones_like(u_), np.cos(v_)),
                    color=GRAIN, linewidth=0, shade=True, zorder=8)
    ax.set_xlim(0, D['Lx'])
    ax.set_ylim(0, D['Lz'])
    ax.set_zlim(0, D['Ly'])
    ax.set_box_aspect((D['Lx'], D['Lz'], D['Ly']), zoom=1.28)
    ax.view_init(elev=13, azim=-62)
    ax.set_title(title, fontsize=10.5, color=C.INK, y=0.99)


fig = plt.figure(figsize=(13.0, 8.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.62], hspace=0.22, wspace=0.05,
                      left=0.04, right=0.975, top=0.845, bottom=0.085)

ax = fig.add_subplot(gs[0, 0], projection='3d', computed_zorder=False)
draw(ax, real, (X[0], X[2]), 'as video 1 draws it — window [0, 4] × [0, 4]')

# roll by half a period; the grain moves with the field
sx, sz = nx // 2, nz // 2
rolled = np.roll(np.roll(real, sz, axis=0), sx, axis=2)
gx = (X[0] + D['Lx'] / 2) % D['Lx']
gz = (X[2] + D['Lz'] / 2) % D['Lz']
ax = fig.add_subplot(gs[0, 1], projection='3d', computed_zorder=False)
draw(ax, rolled, (gx, gz), 'the same instant, window shifted by half a period')

# --- what happened in that column, from the fields ----------------------
ax = C.finish(fig.add_subplot(gs[1, :]))
times, corner_ice, corner_v, corner_front = [], [], [], []
for _, path in C.snapshots(RUN):
    g = C.load_fields(path, names=('VOF/C_L', 'VOF/C_S', 'v'))
    yy = g['y']
    ri = np.where(g['VOF/C_S'] < 0.05, g['ice'], 0.0)
    band = (yy > 1.0) & (yy < 3.4)
    times.append(g['t'])
    corner_ice.append(float(ri[:2, band, :2].max()))
    corner_v.append(float(g['v'][:2, band, :2].mean()))
    solid = ri + g['VOF/C_S']
    rev = solid[:, ::-1, :] < 0.5
    ii = np.argmax(rev, axis=1)
    ii[~rev.any(axis=1)] = len(yy) - 1
    corner_front.append(float(yy[len(yy) - 1 - ii][0, 0]))
times = np.array(times)

ax.plot(times, corner_ice, 'o-', color=C.BLUE, ms=4,
        label='ice fraction in the descending column  (max over 1 < y < 3.4)')
ax.plot(times, -np.array(corner_v), 's-', color=C.ORANGE, ms=4,
        label='downward speed there  $-\\langle v\\rangle$')
ax.axhline(0.9, color=C.STATUS['critical'], lw=1.4, ls='--')
ax.text(1, 0.925, 'VOF_DIFFUSE_ICE_PENAL_THRESHOLD fires at ice fraction 0.9 — '
        'the Darcy mask flips 0 → 1 and the velocity in the cell is clamped to zero',
        color=C.STATUS['critical'], fontsize=8.6)
k = int(np.argmax(np.array(corner_ice) > 0.9))
ax.axvline(times[k], color=C.STATUS['critical'], lw=1.2, ls=':')
ax.annotate('the flow in the column stops here,\n'
            'and it can no longer be flushed by warm water',
            xy=(times[k], 0.5), xytext=(times[k] - 13, 0.62), fontsize=8.6,
            color=C.STATUS['critical'],
            arrowprops=dict(arrowstyle='-|>', color=C.STATUS['critical'], lw=1.1))
ax.set_xlim(0, times.max())
ax.set_ylim(-0.02, 1.05)
ax.set_xlabel('time  $t$')
ax.set_ylabel('ice fraction  /  downward speed')
ax.set_title('a runaway, not a boundary effect: the plume seeds ice, '
             'the threshold freezes the flow, the column locks solid')
ax.legend(loc='upper left')

fig.suptitle('Diagnostic · the four “corner” ice tongues are one column at the '
             'periodic antipode of the grain',
             x=0.04, ha='left', fontsize=15, fontweight='semibold', color=C.INK, y=0.965)
fig.text(0.04, 0.895,
         f"Accepted B.2 coupled gate, job 20492853, at t = {f['t']:.1f}. x and z are "
         'periodic, so nothing physical happens at x = 0 or z = 0 — the drawing '
         'window is arbitrary.',
         ha='left', fontsize=9.5, color=C.INK2)

C.save(fig, 'diag_periodic_shift',
       'Source: B2closure_coupled16_rampfix_20492853 Data_*.h5. Ice is the project’s '
       'real_ice: max(0, 1 − C_L − C_S) masked to C_S < 0.05. This diagnostic was made '
       'to answer a question about video 1 and has not been through the campaign audit '
       'chain — see README, figure 11d.')
