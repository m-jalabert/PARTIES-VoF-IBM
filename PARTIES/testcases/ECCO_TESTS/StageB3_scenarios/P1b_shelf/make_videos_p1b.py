#!/usr/bin/env python3
"""Videos of P1b (job 20977443): a fine-sand grain released beneath an Antarctic
ice-shelf base, falling through the meltwater sublayer.

  python make_videos_p1b.py 2d      --out DIR   full-depth mid-plane, 3 fields + grain speed
  python make_videos_p1b.py closeup --out DIR   window following the grain (wake, streamlines)
  python make_videos_p1b.py 3d      --out DIR   rotating 3-D view: ice, meltwater layer, dragged meltwater, grain
  options: --frames N (render only the first N, for a test)  --procs P  --fps F

Meltwater "excess" = P1b tracer minus the horizontally uniform no-grain control
(20977444, interpolated in time), so it isolates what the grain transports.
Videos stop at t = 18.3, the end of the control.  Run with the Anvil anaconda python.
"""
import argparse
import glob
import re
import subprocess
from multiprocessing import Pool
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

ROOT = Path('/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf')
RUN, CTRL = ROOT / 'run_prod_20977443', ROOT / 'run_ctrl_20977444'
TREF, UREF, D_MM = 0.027382, 7.304, 0.2        # s, mm/s, mm
DX, Y_INT, K_MID, T_MAX = 1 / 24, 20.0, 96, 18.3

SURF, INK, INK2 = '#fcfcfb', '#0b0b0b', '#52514e'
SEQ = LinearSegmentedColormap.from_list('seq', ['#fcfcfb', '#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b'])
DIV = LinearSegmentedColormap.from_list('div', ['#104281', '#3987e5', '#b7d3f6', '#f0efec',
                                                '#f4b8b5', '#e34948', '#8f1f1f'])
ICE_RGBA = (0.80, 0.89, 0.98, 0.85)
GRAIN = '#2f2f2d'


def frames():
    fs = sorted(glob.glob(str(RUN / 'Data_*.h5')), key=lambda p: int(re.findall(r'_(\d+)\.h5', p)[0]))
    out = []
    for p in fs:
        with h5py.File(p) as f:
            t = float(f['time'][0])
        if t <= T_MAX + 1e-6:
            out.append((p, t))
    return out


def control_profile():
    a = np.genfromtxt(CTRL / 'ecco_profiles.csv', delimiter=',', names=True)
    y = np.unique(a['y'])
    T, P = [], []
    for t in np.unique(a['time']):
        b = a[a['time'] == t][-len(y):]
        if len(b) == len(y):
            T.append(t)
            P.append(b[np.argsort(b['y'])]['tracer_integral'] / (64.0 * DX))
    return np.array(T), np.array(P)


CT, CP = control_profile()
M = np.loadtxt(RUN / 'mobile.dat', delimiter=',')


def cmean(t):
    return np.array([np.interp(t, CT, CP[:, j]) for j in range(CP.shape[1])])


def grain(t):
    return np.interp(t, M[:, 0], M[:, 3]), np.interp(t, M[:, 0], M[:, 6])


def style():
    plt.rcParams.update({'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF,
                         'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                         'ytick.color': INK2, 'text.color': INK, 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})


def read_slice(p, k=K_MID):
    with h5py.File(p) as f:
        g = lambda name: f[name][k, :-1, :-1]
        v = f['v'][k, :, :-1]
        u = f['u'][k, :-1, :]
        return dict(t=float(f['time'][0]), tr=g('Conc/2'), cl=g('VOF/C_L'), cs=g('VOF/C_S'),
                    v=0.5 * (v[:-1] + v[1:]), u=0.5 * (u[:, :-1] + u[:, 1:]),
                    x=f['grid/xc'][...].ravel()[:-1], y=f['grid/yc'][...].ravel()[:-1])


def grain_disc(ax, zg):
    from matplotlib.patches import Circle
    ax.add_patch(Circle((4.0, zg), 0.53, facecolor=GRAIN, edgecolor=GRAIN, lw=0.8, zorder=5))  # covers the cell-mask staircase


def overlay_solids(ax, s, X, Z):
    ice = np.where((s['cl'] < 0.5) & (s['cs'] < 0.05), 1.0, np.nan)
    ax.pcolormesh(X, Z, ice, cmap=ListedColormap([ICE_RGBA]), shading='auto')
    sed = np.where(s['cs'] >= 0.05, 1.0, np.nan)
    ax.pcolormesh(X, Z, sed, cmap=ListedColormap([GRAIN]), shading='auto')


# ------------------------------------------------------------------ 2D full depth
def render_2d(args):
    i, p, t, out = args
    s = read_slice(p)
    X, Z = s['x'], Y_INT - s['y']
    ex = s['tr'] - cmean(s['t'])[:, None]
    vmm = s['v'] * UREF
    yg, vg = grain(s['t'])
    style()
    fig = plt.figure(figsize=(16, 9), dpi=120)
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.25], wspace=0.35, left=0.04, right=0.97, top=0.86, bottom=0.08)
    panels = [(s['tr'], SEQ, 0, 0.35, 'Meltwater fraction'),
              (ex, DIV, -0.25, 0.25, 'Meltwater carried by the grain\n(minus no-grain control)'),
              (vmm, DIV, -15, 15, 'Vertical velocity (mm/s)')]
    for c, (field, cmap, lo, hi, title) in enumerate(panels):
        ax = fig.add_subplot(gs[0, c])
        im = ax.pcolormesh(X, Z, np.where(s['cs'] >= 0.05, np.nan, field), cmap=cmap, vmin=lo, vmax=hi, shading='auto')
        overlay_solids(ax, s, X, Z)
        grain_disc(ax, Y_INT - yg)
        ax.set_aspect('equal'); ax.set_ylim(20.0, -2.0); ax.set_xlim(0, 8)
        ax.set_title(title, loc='left', fontsize=11)
        ax.set_xlabel('x (grain diameters)')
        if c == 0:
            ax.set_ylabel('depth below the ice base (grain diameters)')
        cb = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.035, pad=0.07)
        cb.ax.tick_params(labelsize=9)
        if c == 2:
            cb.set_label('blue = downward', fontsize=9, color=INK2)
    ax = fig.add_subplot(gs[0, 3])
    sel = M[:, 0] <= s['t']
    ax.plot(M[:, 0] * TREF * 1e3, -M[:, 6] * UREF, color='#e9e8e3', lw=2)
    ax.plot(M[sel, 0] * TREF * 1e3, -M[sel, 6] * UREF, color='#2a78d6', lw=2.2)
    ax.plot([s['t'] * TREF * 1e3], [-vg * UREF], 'o', color='#2a78d6', ms=8, mec=SURF, mew=2)
    ax.set_xlabel('time (ms)'); ax.set_ylabel('grain settling speed (mm/s)')
    ax.set_title('Grain speed', loc='left', fontsize=11)
    ax.grid(True, color='#e4e3df', lw=0.6)
    ax.set_xlim(0, T_MAX * TREF * 1e3); ax.set_ylim(0, 16)
    fig.suptitle('Fine-sand grain (0.2 mm) released beneath an Antarctic ice-shelf base', x=0.04, ha='left',
                 fontsize=16, y=0.975)
    fig.text(0.04, 0.915, f't = {s["t"] * TREF * 1e3:6.1f} ms   (code time {s["t"]:5.2f})    '
             f'grain depth {Y_INT - yg:5.2f} d    speed {-vg * UREF:5.2f} mm/s    '
             'mid-plane slice through the grain; ice in pale blue, grain in charcoal',
             color=INK2, fontsize=11, ha='left')
    fig.savefig(out / f'frame_{i:04d}.png'); plt.close(fig)


# ------------------------------------------------------------------ close-up following the grain
def render_closeup(args):
    i, p, t, out = args
    s = read_slice(p)
    yg, vg = grain(s['t'])
    zg = Y_INT - yg
    X, Z = s['x'], Y_INT - s['y']
    ex = s['tr'] - cmean(s['t'])[:, None]
    z0, z1 = zg - 4.4, zg + 1.6                       # wake trails above the falling grain
    if z0 < -0.6:
        z0, z1 = -0.6, 5.4
    if z1 > 20.0:
        z0, z1 = 14.0, 20.0
    xs = (X > 1.9) & (X < 6.1)
    zs = (Z > z0 - 0.2) & (Z < z1 + 0.2)
    style()
    fig, ax = plt.subplots(figsize=(8, 10), dpi=120)
    fig.subplots_adjust(left=0.13, right=0.97, top=0.88, bottom=0.14)
    im = ax.pcolormesh(X, Z, np.where(s['cs'] >= 0.05, np.nan, ex), cmap=DIV, vmin=-0.25, vmax=0.25, shading='auto')
    overlay_solids(ax, s, X, Z)
    grain_disc(ax, zg)
    # streamlines of the in-plane velocity relative to the grain (depth axis points down: dz/dt = -v)
    uu = s['u'][np.ix_(zs, xs)]
    ww = -(s['v'][np.ix_(zs, xs)] - vg)
    solid = s['cs'][np.ix_(zs, xs)] >= 0.05
    uu, ww = np.where(solid, np.nan, uu), np.where(solid, np.nan, ww)
    Xs, Zs = X[xs], Z[zs]
    order = np.argsort(Zs)
    ax.streamplot(Xs, Zs[order], uu[order], ww[order], color=INK2, density=1.3, linewidth=0.6, arrowsize=0.7)
    ax.set_aspect('equal'); ax.set_xlim(2.0, 6.0); ax.set_ylim(z1, z0)
    ax.set_xlabel('x (grain diameters)'); ax.set_ylabel('depth below the ice base (grain diameters)')
    cb = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.04, pad=0.08)
    cb.set_label('meltwater carried by the grain (minus no-grain control)')
    fig.suptitle('Following the grain: its wake and the meltwater it drags down', x=0.13, ha='left', fontsize=14, y=0.975)
    fig.text(0.13, 0.925, f't = {s["t"] * TREF * 1e3:6.1f} ms    depth {zg:5.2f} d    speed {-vg * UREF:5.2f} mm/s    '
             'streamlines: flow seen from the grain', color=INK2, fontsize=11, ha='left')
    fig.savefig(out / f'frame_{i:04d}.png'); plt.close(fig)


# ------------------------------------------------------------------ 3D
def shade(verts, faces, normals, rgb, light=(0.35, -0.55, 0.75)):
    L = np.array(light) / np.linalg.norm(light)
    n = normals[faces].mean(axis=1)
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    inten = 0.30 + 0.70 * np.abs(n @ L)
    return np.clip(np.outer(inten, rgb), 0, 1)


def iso(field, level, step):
    from skimage.measure import marching_cubes
    if not (field.min() < level < field.max()):
        return None
    v, f_, n, _ = marching_cubes(field, level=level, spacing=(step * DX,) * 3)
    v += 0.5 * DX                                          # cell centres
    return v, f_, n


def render_3d(args):
    i, p, t, out = args
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    step = 2
    with h5py.File(p) as f:
        tt = float(f['time'][0])
        tr = f['Conc/2'][...][:-1, :-1, :-1]
    exc = tr - cmean(tt)[None, :, None]
    yg, vg = grain(tt)
    style()
    fig = plt.figure(figsize=(10.8, 10.8), dpi=100)
    ax = fig.add_axes([0.0, 0.07, 1.0, 0.82], projection='3d')
    ax.set_box_aspect((8, 8, 22), zoom=1.0)
    surfaces = [(exc[::step, ::step, ::step], 0.04, np.array([0.89, 0.29, 0.28]), 0.92, 'excess'),
                (tr[::step, ::step, ::step], 0.15, np.array([0.22, 0.53, 0.90]), 0.16, 'layer')]
    for fld, lev, rgb, alpha, _ in surfaces:
        r = iso(fld, lev, step)
        if r is None:
            continue
        v, f_, n = r
        # array axes are (z, y, x) -> plot (x, z, y) so height y is vertical
        P = v[:, [2, 0, 1]]
        Nn = n[:, [2, 0, 1]]
        col = shade(P, f_, Nn, rgb)
        pc = Poly3DCollection(P[f_], facecolors=np.c_[col, np.full(len(col), alpha)], edgecolors='none')
        ax.add_collection3d(pc)
    # ice slab (static) and grain
    xx, zz = np.meshgrid([0, 8], [0, 8])
    for yy in (Y_INT, 22.0):
        ax.plot_surface(xx, zz, np.full_like(xx, yy, dtype=float), color=ICE_RGBA[:3], alpha=0.35, shade=False)
    for a, b in [((0, 0), (8, 0)), ((0, 8), (8, 8)), ((0, 0), (0, 8)), ((8, 0), (8, 8))]:
        X_ = np.array([[a[0], b[0]], [a[0], b[0]]]); Z_ = np.array([[a[1], b[1]], [a[1], b[1]]])
        Y_ = np.array([[Y_INT, Y_INT], [22, 22]])
        ax.plot_surface(X_, Z_, Y_, color=ICE_RGBA[:3], alpha=0.25, shade=False)
    uu, vv = np.mgrid[0:2 * np.pi:28j, 0:np.pi:16j]
    ax.plot_surface(4 + 0.5 * np.cos(uu) * np.sin(vv), 4 + 0.5 * np.sin(uu) * np.sin(vv), yg + 0.5 * np.cos(vv),
                    color=GRAIN, shade=True, linewidth=0)
    for xe in (0, 8):
        for ze in (0, 8):
            ax.plot([xe, xe], [ze, ze], [0, 22], color='#c3c2b7', lw=0.6)
    for ye in (0, 22):
        ax.plot([0, 8, 8, 0, 0], [0, 0, 8, 8, 0], [ye] * 5, color='#c3c2b7', lw=0.6)
    ax.set_xlim(0, 8); ax.set_ylim(0, 8); ax.set_zlim(0, 22)
    ax.view_init(elev=18, azim=-55 + 1.0 * i)
    ax.set_axis_off()
    fig.text(0.04, 0.965, 'Fine-sand grain beneath an Antarctic ice-shelf base (3-D)', fontsize=17, ha='left', va='top')
    fig.text(0.04, 0.925, f't = {tt * TREF * 1e3:6.1f} ms    grain depth {Y_INT - yg:5.2f} d    speed {-vg * UREF:5.2f} mm/s',
             fontsize=12, color=INK2, ha='left', va='top')
    fig.text(0.04, 0.03, 'pale blue: ice    translucent blue: lower edge of the meltwater layer (meltwater fraction 0.15)\n'
             'red: meltwater dragged down by the grain (excess over the no-grain control > 0.04)    '
             'charcoal: grain (0.2 mm)    box 1.6 × 4.4 × 1.6 mm',
             fontsize=10.5, color=INK2, ha='left')
    fig.savefig(out / f'frame_{i:04d}.png'); plt.close(fig)


def encode(frames_dir, mp4, fps):
    import imageio_ffmpeg
    exe = imageio_ffmpeg.get_ffmpeg_exe()
    subprocess.run([exe, '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', str(frames_dir / 'frame_%04d.png'),
                    '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2', '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '20',
                    '-movflags', '+faststart', str(mp4)], check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['2d', 'closeup', '3d'])
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--frames', type=int, default=0)
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--fps', type=int, default=15)
    a = ap.parse_args()
    fr = frames()
    if a.mode == 'closeup':
        fr = [(p, t) for p, t in fr if t <= 11.0]
    if a.frames:
        fr = fr[:a.frames]
    work = a.out / f'frames_{a.mode}'
    work.mkdir(parents=True, exist_ok=True)
    fn = {'2d': render_2d, 'closeup': render_closeup, '3d': render_3d}[a.mode]
    jobs = [(i, p, t, work) for i, (p, t) in enumerate(fr)]
    with Pool(a.procs) as pool:
        pool.map(fn, jobs, chunksize=1)
    if not a.frames or a.frames > 5:
        name = {'2d': 'P1b_2D_midplane.mp4', 'closeup': 'P1b_2D_grain_closeup.mp4', '3d': 'P1b_3D_meltwater.mp4'}[a.mode]
        encode(work, a.out / name, a.fps)
        print('wrote', a.out / name)


if __name__ == '__main__':
    main()
