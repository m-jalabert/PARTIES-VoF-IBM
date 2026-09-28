#!/usr/bin/env python3
"""Video 2 — four fields on the vertical plane through the grain.

temperature θ · meltwater anomaly · vertical velocity v · speed |u|

Same run and the same time convention as video 1: Eulerian fields are linearly
blended between the 33 saved outputs (Δt = 2), while the grain is drawn at its
exact position from mobile.dat.  The ice is overlaid as a pale block with its
melt front; the grain's IBM support is masked, because scalar values inside it
are extension values rather than physical concentrations.

    python video2_slices.py            # full movie
    python video2_slices.py --sub 1    # one frame per saved output (fast)
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stageb_common as C

plt = C.use_style()
C.use_ffmpeg()
CM = C.cmaps()

ap = argparse.ArgumentParser()
ap.add_argument('--sub', type=int, default=4)
ap.add_argument('--fps', type=int, default=12)
ap.add_argument('--dpi', type=int, default=130)
ap.add_argument('--limit', type=int, default=0)
ap.add_argument('--out', default=str(C.VIDDIR / 'video2_slices.mp4'))
args = ap.parse_args()

RUN = C.run('coupled16')
D = C.COUPLED16
rel = C.release_time(RUN)
tm, Xm, Um = C.mobile(RUN)
snaps = C.snapshots(RUN)
wall_hit = float(tm[np.argmax(Xm[:, 1] < 0.51)])
GRAIN = '#4f4840'

NAMES = ('VOF/C_L', 'VOF/C_S', 'Conc/0', 'Conc/2', 'u', 'v', 'w')


def load(path):
    f = C.load_fields(path, names=NAMES)
    f['real_ice'] = np.where(f['VOF/C_S'] < 0.05, f['ice'], 0.0)
    f['speed'] = np.sqrt(f['u'] ** 2 + f['v'] ** 2 + f['w'] ** 2)
    return f


from matplotlib.animation import FFMpegWriter    # noqa: E402

fig = plt.figure(figsize=(13.0, 6.3))
fig.patch.set_facecolor(C.PAGE)
# equal-aspect panels shrink inside their gridspec cell, so place the colour
# bars from the height the 4:6 domain actually takes rather than guessing.
LEFT, RIGHT, TOP, WSPACE = 0.035, 0.985, 0.815, 0.30
FW, FH = fig.get_size_inches()
PW = (RIGHT - LEFT) / (4 + 3 * WSPACE)
PH = (D['Ly'] / D['Lx']) * PW * FW / FH
gs = fig.add_gridspec(1, 4, left=LEFT, right=RIGHT, top=TOP, bottom=TOP - PH,
                      wspace=WSPACE)
axes = [fig.add_subplot(gs[0, i]) for i in range(4)]
for axx in axes:
    axx.set_anchor('N')
CB_Y = TOP - PH - 0.125
caxes = [fig.add_axes([LEFT + i * (1 + WSPACE) * PW, CB_Y, PW, 0.016])
         for i in range(4)]

TMAX, VMAX, SMAX = 0.08, 0.9, 1.1
SPECS = [
    ('temperature  θ', CM['orange'], 0.0, 1.0, 'Conc/0', None),
    (r'meltwater anomaly  $C_{mw}-\langle C_{mw}\rangle$', CM['diverging'],
     -TMAX, TMAX, 'Conc/2', 'anomaly'),
    ('vertical velocity  $v$', CM['diverging'], -VMAX, VMAX, 'v', None),
    ('speed  $|u|$', CM['blue'], 0.0, SMAX, 'speed', None),
]

fig.suptitle('Stage B · four fields on the plane through the grain',
             x=0.035, ha='left', fontsize=16, fontweight='semibold',
             color=C.INK, y=0.982)
head = fig.text(0.035, 0.934,
                'Accepted B.2 coupled gate, job 20492853 — 64 × 96 × 64, d/Δx = 16, '
                'cold top wall, warm bottom wall. Synthetic debugging parameters, '
                'not a physical ECCO scenario.',
                ha='left', fontsize=9, color=C.INK2)
state = fig.text(0.035, 0.876, '', ha='left', fontsize=12.5, color=C.BLUE,
                 fontweight='semibold')
fig.text(0.035, 0.018,
         'Eulerian fields blended between the 33 saved outputs (Δt = 2); grain position '
         'from mobile.dat at full step resolution. Values inside the IBM support are '
         'extension values and are masked.',
         ha='left', fontsize=7.5, color=C.MUTED)

writer = FFMpegWriter(fps=args.fps, codec='h264', bitrate=7000,
                      metadata=dict(title='PARTIES Stage B — mid-plane fields'))

t0, nframe = time.time(), 0
Path(args.out).parent.mkdir(exist_ok=True)
with writer.saving(fig, args.out, args.dpi):
    prev = None
    for k in range(len(snaps) - 1):
        fa = prev if prev is not None else load(snaps[k][1])
        fb = load(snaps[k + 1][1])
        prev = fb
        for sidx in range(args.sub):
            w = sidx / args.sub
            t = (1 - w) * fa['t'] + w * fb['t']
            blend = {key: (1 - w) * fa[key] + w * fb[key]
                     for key in ('real_ice', 'VOF/C_S', 'Conc/0', 'Conc/2', 'v', 'speed')}
            x, y, z = fa['x'], fa['y'], fa['z']

            i = min(max(int(np.searchsorted(tm, t)), 1), len(tm) - 1)
            f2 = (t - tm[i - 1]) / max(tm[i] - tm[i - 1], 1e-12)
            Xg = (1 - f2) * Xm[i - 1] + f2 * Xm[i]
            Ug = (1 - f2) * Um[i - 1] + f2 * Um[i]

            jz = int(np.argmin(np.abs(z - Xg[2])))
            ice2 = blend['real_ice'][jz]
            cs2 = blend['VOF/C_S'][jz]
            solid2 = (cs2 >= 0.05) | (ice2 > 0.5)
            liquid = (blend['VOF/C_S'] < 0.05) & (blend['real_ice'] < 0.5)
            ext = [x[0], x[-1], y[0], y[-1]]

            for axx, cax, (title, cmap, lo, hi, key, mode) in zip(axes, caxes, SPECS):
                axx.clear()
                axx.grid(False)
                axx.set_facecolor('#eceae5')
                field = blend[key]
                if mode == 'anomaly':
                    field = field - (float(field[liquid].mean()) if liquid.any() else 0.0)
                h = axx.imshow(np.ma.masked_where(solid2, field[jz]), origin='lower',
                               extent=ext, aspect='equal', cmap=cmap,
                               vmin=lo, vmax=hi, interpolation='bilinear')
                axx.imshow(np.ma.masked_where(ice2 <= 0.5, np.ones_like(ice2)),
                           origin='lower', extent=ext, aspect='equal', zorder=3,
                           cmap=plt.matplotlib.colors.ListedColormap(['#e3edf7']),
                           vmin=0, vmax=1, interpolation='nearest')
                axx.contour(x, y, ice2, levels=[0.5], colors=[C.AQUA],
                            linewidths=1.1, zorder=4)
                axx.contourf(x, y, cs2, levels=[0.5, 1.01], colors=[GRAIN], zorder=5)
                axx.set_xlim(x[0], x[-1])
                axx.set_ylim(y[0], y[-1])
                axx.set_xticks([0, 2, 4])
                axx.set_yticks([0, 2, 4, 6])
                axx.tick_params(labelsize=8)
                axx.set_xlabel('$x$')
                if axx is axes[0]:
                    axx.set_ylabel('$y$')
                else:
                    axx.set_yticklabels([])
                for sp in axx.spines.values():
                    sp.set_color(C.AXIS)
                axx.set_title(title, fontsize=10.5, color=C.INK)
                cax.clear()
                cb = fig.colorbar(h, cax=cax, orientation='horizontal')
                cb.outline.set_edgecolor(C.AXIS)
                cb.ax.tick_params(labelsize=7.5, colors=C.MUTED)

            phase = ('locked in the ice' if t < rel else
                     'falling' if t < wall_hit else 'at rest on the wall')
            col = (C.BLUE if t < rel else C.ORANGE if t < wall_hit else C.AQUA)
            state.set_text(f't = {t:6.2f}   ·   {phase}   ·   '
                           f'y = {Xg[1]:.3f}   ·   $U_y$ = {Ug[1]:+.3f}')
            state.set_color(col)

            writer.grab_frame(facecolor=C.PAGE)
            nframe += 1
            if nframe % 20 == 0:
                print(f'  {nframe} frames, t = {t:.1f}, '
                      f'{(time.time()-t0)/nframe:.2f} s/frame', flush=True)
            if args.limit and nframe >= args.limit:
                break
        if args.limit and nframe >= args.limit:
            break

print(f'wrote {args.out}  ({nframe} frames, {time.time()-t0:.0f} s)')
