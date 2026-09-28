#!/usr/bin/env python3
"""Video 1 — the 3-D lifecycle of the accepted B.2 coupled gate.

Left   a 3-D view of the melt front (marching-cubes isosurface of the ACTUAL ice
       field, i.e. 1 - C_L - C_S masked to C_S < 0.05 so the grain's own IBM
       support is never drawn as ice) and the sediment grain.
Right  the vertical plane through the grain: meltwater anomaly, ice, grain.
Below  a timeline of the grain's centre height, with the phase it is in.

Time resolution.  PARTIES wrote 33 field outputs, one every 2.0 time units.
Eulerian fields are linearly blended between consecutive outputs to give a
watchable frame rate; the grain's position and velocity come from mobile.dat at
full step resolution, so the grain moves exactly as it did in the run.  The
overlay states this on every frame.

    python video1_3d_lifecycle.py             # full movie
    python video1_3d_lifecycle.py --sub 1     # one frame per saved output (fast)
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stageb_common as C

plt = C.use_style()
FFMPEG = C.use_ffmpeg()
CM = C.cmaps()

ap = argparse.ArgumentParser()
ap.add_argument('--sub', type=int, default=4,
                help='interpolated frames per saved-output interval')
ap.add_argument('--fps', type=int, default=12)
ap.add_argument('--dpi', type=int, default=120)
ap.add_argument('--limit', type=int, default=0, help='stop after N frames (testing)')
ap.add_argument('--out', default=str(C.VIDDIR / 'video1_3d_lifecycle.mp4'))
args = ap.parse_args()

RUN = C.run('coupled16')
D = C.COUPLED16
rel = C.release_time(RUN)
tm, Xm, Um = C.mobile(RUN)
snaps = C.snapshots(RUN)
wall_hit = float(tm[np.argmax(Xm[:, 1] < 0.51)])

GRAIN = '#4f4840'
ICE_RGB = np.array([0.42, 0.62, 0.84])


def load(path):
    f = C.load_fields(path, names=('VOF/C_L', 'VOF/C_S', 'Conc/2'))
    f['real_ice'] = np.where(f['VOF/C_S'] < 0.05, f['ice'], 0.0)
    return f


def sphere(Xc, R, n=40):
    u = np.linspace(0, 2 * np.pi, n)
    v = np.linspace(0, np.pi, n // 2)
    sx = Xc[0] + R * np.outer(np.cos(u), np.sin(v))
    sz = Xc[2] + R * np.outer(np.sin(u), np.sin(v))
    sy = Xc[1] + R * np.outer(np.ones_like(u), np.cos(v))
    return sx, sz, sy


def box_edges(ax, Lx, Ly, Lz):
    pts = np.array([[0, 0, 0], [Lx, 0, 0], [Lx, Lz, 0], [0, Lz, 0],
                    [0, 0, Ly], [Lx, 0, Ly], [Lx, Lz, Ly], [0, Lz, Ly]], float)
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7)]
    for a, b in edges:
        ax.plot(*zip(pts[a], pts[b]), color=C.AXIS, lw=0.8, alpha=0.8, zorder=0)


from matplotlib.animation import FFMpegWriter          # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402
from skimage import measure                            # noqa: E402

fig = plt.figure(figsize=(12.8, 7.2))
fig.patch.set_facecolor(C.PAGE)
gs = fig.add_gridspec(2, 2, width_ratios=[1.18, 1.0], height_ratios=[1.0, 0.20],
                      left=0.042, right=0.955, top=0.855, bottom=0.085,
                      hspace=0.20, wspace=0.02)
ax3 = fig.add_subplot(gs[0, 0], projection='3d', computed_zorder=False)
ax2 = fig.add_subplot(gs[0, 1])
axt = fig.add_subplot(gs[1, :])

# fixed timeline
axt.plot(tm, Xm[:, 1], color=C.BLUE, lw=1.8)
axt.axvspan(0, rel, color='#eef3fa', lw=0)
axt.axvspan(rel, wall_hit, color='#fdf1e9', lw=0)
axt.axvspan(wall_hit, tm[-1], color='#eef7f3', lw=0)
axt.set_xlim(0, tm[-1])
axt.set_ylim(0, 6)
axt.set_yticks([0, 3, 6])
axt.set_xlabel('time  $t$', labelpad=1)
axt.set_ylabel('$y$', labelpad=1)
C.finish(axt)
cursor = axt.axvline(0, color=C.STATUS['critical'], lw=1.6)
dot, = axt.plot([0], [D['y_grain0']], 'o', ms=6, color=C.STATUS['critical'],
                mec=C.SURFACE, mew=1.2)

fig.suptitle('Stage B · a sediment grain melts out of sea ice, settles, and lands',
             x=0.02, ha='left', fontsize=16, fontweight='semibold', color=C.INK, y=0.975)
sub = fig.text(0.02, 0.918,
               'Accepted B.2 coupled gate, job 20492853 — 64 × 96 × 64, d/Δx = 16, '
               'max_dt = 0.01, cold top wall, warm bottom wall. Synthetic debugging '
               'parameters, not a physical ECCO scenario.',
               ha='left', fontsize=9, color=C.INK2)
fig.text(0.245, 0.845, 'melt front (actual ice) and the grain, in 3-D',
         ha='center', fontsize=10.5, color=C.INK, fontweight='semibold')
fig.text(0.012, 0.010,
         'Eulerian fields blended between the 33 saved outputs (Δt = 2); '
         'grain position from mobile.dat at full step resolution.',
         ha='left', fontsize=7.5, color=C.MUTED)

writer = FFMpegWriter(fps=args.fps, codec='h264', bitrate=6000,
                      metadata=dict(title='PARTIES Stage B — sediment from ice',
                                    artist='ECCO Stage B campaign'))

TMAX = 0.08
t0 = time.time()
nframe = 0
Path(args.out).parent.mkdir(exist_ok=True)
with writer.saving(fig, args.out, args.dpi):
    prev = None
    for k in range(len(snaps) - 1):
        fa = prev if prev is not None else load(snaps[k][1])
        fb = load(snaps[k + 1][1])
        prev = fb
        for s in range(args.sub):
            w = s / args.sub
            t = (1 - w) * fa['t'] + w * fb['t']
            ice = (1 - w) * fa['real_ice'] + w * fb['real_ice']
            mw = (1 - w) * fa['Conc/2'] + w * fb['Conc/2']
            cs = (1 - w) * fa['VOF/C_S'] + w * fb['VOF/C_S']
            x, y, z = fa['x'], fa['y'], fa['z']
            dx, dy, dz = x[1] - x[0], y[1] - y[0], z[1] - z[0]

            i = int(np.searchsorted(tm, t))
            i = min(max(i, 1), len(tm) - 1)
            f2 = (t - tm[i - 1]) / max(tm[i] - tm[i - 1], 1e-12)
            Xg = (1 - f2) * Xm[i - 1] + f2 * Xm[i]
            Ug = (1 - f2) * Um[i - 1] + f2 * Um[i]

            # ---- 3-D panel -------------------------------------------------
            ax3.clear()
            ax3.set_facecolor(C.PAGE)
            ax3.set_axis_off()
            box_edges(ax3, D['Lx'], D['Ly'], D['Lz'])
            if ice.max() > 0.5:
                v, faces, nrm, _ = measure.marching_cubes(
                    ice, 0.5, spacing=(dz, dy, dx), step_size=1)
                verts = np.stack([v[:, 2] + x[0], v[:, 0] + z[0], v[:, 1] + y[0]], 1)
                tri = verts[faces]
                nn = np.stack([nrm[:, 2], nrm[:, 0], nrm[:, 1]], 1)[faces].mean(1)
                nn /= np.linalg.norm(nn, axis=1)[:, None] + 1e-12
                # headlight-ish: above and in front of the camera.  |n.l| because
                # the isosurface normals are not consistently outward-oriented.
                lightdir = np.array([0.35, -0.55, 0.76])
                lightdir /= np.linalg.norm(lightdir)
                diff = np.abs(nn @ lightdir)
                cols = np.clip(ICE_RGB[None, :] * (0.45 + 0.55 * diff)[:, None]
                               + 0.30 * (diff ** 14)[:, None], 0, 1)
                fc = np.c_[cols, np.full(len(cols), 1.0)]
                pc = Poly3DCollection(tri, facecolors=fc, edgecolors=fc,
                                      linewidths=0.0, antialiased=False, zorder=4)
                ax3.add_collection3d(pc)
            sx, sz, sy = sphere(Xg, D['R'])
            ax3.plot_surface(sx, sz, sy, color=GRAIN, linewidth=0, shade=True,
                             antialiased=True, zorder=8)
            ax3.set_xlim(0, D['Lx'])
            ax3.set_ylim(0, D['Lz'])
            ax3.set_zlim(0, D['Ly'])
            ax3.set_box_aspect((D['Lx'], D['Lz'], D['Ly']), zoom=1.30)
            ax3.view_init(elev=13, azim=-62 + 10 * np.sin(2 * np.pi * nframe / 260))


            # ---- slice panel -----------------------------------------------
            ax2.clear()
            ax2.grid(False)
            ax2.set_facecolor('#eceae5')
            jz = int(np.argmin(np.abs(z - Xg[2])))
            ice2, cs2 = ice[jz], cs[jz]
            solid = (cs2 >= 0.05) | (ice2 > 0.5)
            liquid3 = (cs < 0.05) & (ice < 0.5)
            mw_mean = float(mw[liquid3].mean()) if liquid3.any() else 0.0
            ext = [x[0], x[-1], y[0], y[-1]]
            ax2.imshow(np.ma.masked_where(solid, mw[jz] - mw_mean), origin='lower',
                       extent=ext, aspect='equal', cmap=CM['diverging'],
                       vmin=-TMAX, vmax=TMAX, interpolation='bilinear')
            ax2.imshow(np.ma.masked_where(ice2 <= 0.5, np.ones_like(ice2)),
                       origin='lower', extent=ext, aspect='equal', zorder=3,
                       cmap=plt.matplotlib.colors.ListedColormap(['#e3edf7']),
                       vmin=0, vmax=1, interpolation='nearest')
            ax2.contour(x, y, ice2, levels=[0.5], colors=[C.AQUA], linewidths=1.2,
                        zorder=4)
            ax2.contourf(x, y, cs2, levels=[0.5, 1.01], colors=[GRAIN], zorder=5)
            ax2.set_xlim(x[0], x[-1])
            ax2.set_ylim(y[0], y[-1])
            ax2.set_xticks([0, 2, 4])
            ax2.set_yticks([0, 2, 4, 6])
            ax2.tick_params(labelsize=8)
            ax2.set_xlabel('$x$')
            ax2.set_ylabel('$y$')
            for sp in ax2.spines.values():
                sp.set_color(C.AXIS)
            ax2.set_title(r'meltwater anomaly on $z = z_{grain}$',
                          fontsize=10.5, color=C.INK)

            phase = ('locked in the ice' if t < rel else
                     'falling' if t < wall_hit else 'at rest on the wall')
            col = (C.BLUE if t < rel else C.ORANGE if t < wall_hit else C.AQUA)
            ax2.text(0.03, 0.975, f't = {t:6.2f}\n{phase}\n$U_y$ = {Ug[1]:+.3f}',
                     transform=ax2.transAxes, fontsize=10, color=col, va='top',
                     zorder=9, linespacing=1.5,
                     bbox=dict(boxstyle='round,pad=0.35', fc='#ffffffd9', ec='none'))

            cursor.set_xdata([t, t])
            dot.set_data([t], [Xg[1]])

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
