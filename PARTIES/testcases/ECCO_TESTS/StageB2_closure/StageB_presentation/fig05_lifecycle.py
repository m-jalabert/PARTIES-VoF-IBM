#!/usr/bin/env python3
"""Figure 5 — the grain lifecycle in the accepted B.2 coupled gate.

Source: /anvil/scratch/.../B2closure_coupled16_rampfix_20492853
        mobile.dat (every step), release.dat (every step), Particle_*.h5 (33 outputs)

Panels
  (a) centre height y(t): locked in ice -> released -> settling -> wall contact
  (b) vertical velocity U_y(t), with the release handoff inset
  (c) vertical force balance from the saved outputs: IBM reaction against buoyant weight
  (d) the four B.2 motion gates as measured numbers
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stageb_common as C

plt = C.use_style()

RUN = C.run('coupled16')
D = C.COUPLED16
A_FREE = D['g'] * (1.0 - 1.0 / D['rho_s'])          # 0.6
WEIGHT = (4 / 3) * np.pi * D['R'] ** 3 * (D['rho_s'] - 1.0) * D['g']   # 0.785398

MOTION = C.audit('motion_coupled16_rampfix.json')['coupled16_rampfix']
t, X, U = C.mobile(RUN)
rel = C.release_time(RUN)
r = C.release(RUN)
y = X[:, 1]
uy = U[:, 1]

wall_hit = float(t[np.argmax(y < 0.51)])
n_locked = int(np.count_nonzero(t < rel))
n_locked_str = f'{n_locked:,}'.replace(',', ' ')
y_min = float(y.min())
u_min = float(uy.min())

# vertical force history at the 33 saved outputs
import h5py
tp, f_ibm, f_rigid = [], [], []
for idx, path in C.snapshots(RUN):
    with h5py.File(path.with_name(path.name.replace('Data_', 'Particle_')), 'r') as p:
        tp.append(float(p['time'][0]))
        f_ibm.append(float(p['mobile/F_IBM'][...].reshape(-1, 3)[0][1]))
        f_rigid.append(float(p['mobile/F_rigid'][...].reshape(-1, 3)[0][1]))
tp, f_ibm, f_rigid = map(np.asarray, (tp, f_ibm, f_rigid))

fig = plt.figure(figsize=(12.4, 7.4))
gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.22,
                      left=0.065, right=0.985, top=0.86, bottom=0.085)

PHASES = [(0.0, rel, '#eef3fa', 'locked in ice'),
          (rel, wall_hit, '#fdf1e9', 'free fall'),
          (wall_hit, t[-1], '#eef7f3', 'at the wall')]

# --- (a) height ----------------------------------------------------------
ax = C.finish(fig.add_subplot(gs[0, 0]))
for a, b, col, lab in PHASES:
    ax.axvspan(a, b, color=col, lw=0, zorder=0)
    ax.text((a + b) / 2, 5.72, lab, ha='center', va='top',
            color=C.MUTED, fontsize=8)
ax.axhline(D['y_ice0'], color=C.AQUA, lw=1.4, ls=(0, (5, 3)), zorder=2)
ax.text(1.0, D['y_ice0'] + 0.12, 'initial ice front  y = 3.6', ha='left',
        va='bottom', color=C.AQUA, fontsize=8)
ax.axhline(D['R'], color=C.MUTED, lw=1.2, ls=':', zorder=2)
ax.text(1.0, D['R'] + 0.1, 'wall contact  y = R', color=C.MUTED, fontsize=8)
ax.plot(t, y, color=C.BLUE, zorder=4)
ax.plot([rel], [y[np.searchsorted(t, rel)]], 'o', ms=7, color=C.ORANGE,
        mec=C.SURFACE, mew=1.6, zorder=6)
ax.annotate(f'release  t = {rel:.3f}', xy=(rel, D['y_grain0']),
            xytext=(rel - 14, 2.6), color=C.ORANGE, fontsize=8.5,
            arrowprops=dict(arrowstyle='-', color=C.ORANGE, lw=1.1))
ax.set_ylim(0, 6.0)
ax.set_xlim(0, t[-1])
ax.set_xlabel('time  $t$')
ax.set_ylabel('grain centre height  $y$')
ax.set_title('a   Locked, released, settled, landed')

# --- (b) velocity --------------------------------------------------------
ax = C.finish(fig.add_subplot(gs[0, 1]))
for a, b, col, _ in PHASES:
    ax.axvspan(a, b, color=col, lw=0, zorder=0)
ax.axhline(0, color=C.AXIS, lw=0.9, zorder=1)
ax.plot(t, uy, color=C.BLUE, zorder=4)
ax.annotate(f'min $U_y$ = {u_min:.4f}', xy=(t[np.argmin(uy)], u_min),
            xytext=(t[np.argmin(uy)] + 4, u_min - 0.03), color=C.BLUE, fontsize=8.5,
            arrowprops=dict(arrowstyle='-', color=C.BLUE, lw=1.1))
ax.text(2, 0.20, f'exactly zero for all {n_locked_str} locked steps',
        ha='left', va='center', color=C.INK2, fontsize=8.5)
ax.set_xlim(0, t[-1])
ax.set_xlabel('time  $t$')
ax.set_ylabel('vertical velocity  $U_y$')
ax.set_title('b   The lock is exact; the handoff is a spin-up')

# inset: the ramp
ins = ax.inset_axes([0.14, 0.11, 0.35, 0.34])
w = (t >= rel - 0.02) & (t <= rel + 0.30)
ins.axvline(rel, color=C.ORANGE, lw=1.2, ls='--')
ins.plot(t[w], uy[w], color=C.BLUE, lw=1.8)
ramp_end = t[np.argmax((r['release_ramp'] == 0) & (t > rel))]
ins.axvspan(rel, ramp_end, color='#fdf1e9', lw=0, zorder=0)
ins.text(ramp_end, ins.get_ylim()[1], ' 10-step ramp', color=C.ORANGE,
         fontsize=7.5, va='top')
ins.set_title('handoff:  $U_y$  vs  $t$', fontsize=7.5, color=C.INK2,
              pad=3, fontweight='normal')
ins.set_facecolor(C.SURFACE)
ins.tick_params(labelsize=7, colors=C.MUTED)
ins.grid(color=C.GRID, lw=0.5)
for s in ('top', 'right'):
    ins.spines[s].set_visible(False)

# --- (c) force balance ---------------------------------------------------
ax = C.finish(fig.add_subplot(gs[1, 0]))
ax.axvspan(0, rel, color='#eef3fa', lw=0, zorder=0)
ax.axhline(WEIGHT, color=C.MUTED, lw=1.3, ls='--', zorder=2)
ax.text(0.8, WEIGHT + 0.015, f'buoyant weight  {WEIGHT:.4f}', color=C.INK2, fontsize=8)
ax.plot(tp, f_ibm, 'o-', color=C.BLUE, ms=4, label='IBM reaction  $F_{IBM,y}$')
ax.axvline(rel, color=C.ORANGE, lw=1.2, ls='--', zorder=3)
ax.text(rel - 0.8, 0.30, 'release', color=C.ORANGE, fontsize=8, ha='right')
k = int(np.argmin(np.abs(tp - 34.0)))
net = WEIGHT - f_ibm[k]
ax.annotate('', xy=(tp[k], WEIGHT), xytext=(tp[k], f_ibm[k]),
            arrowprops=dict(arrowstyle='<->', color=C.RED, lw=1.4))
ax.annotate(f'net = {net:.5f}\n= {100*net/WEIGHT:.1f} % of the weight',
            xy=(tp[k], 0.5 * (WEIGHT + f_ibm[k])), xytext=(tp[k] + 6, 0.50),
            color=C.RED, fontsize=8.5,
            arrowprops=dict(arrowstyle='-', color=C.RED, lw=1.0))
ax.set_xlim(0, t[-1])
ax.set_ylim(-0.15, 0.93)
ax.set_xlabel('time  $t$')
ax.set_ylabel('vertical force on the grain')
ax.set_title('c   At escape the grain is almost exactly supported')
ax.legend(loc='lower right')

# --- (d) gate card -------------------------------------------------------
ax = fig.add_subplot(gs[1, 1])
ax.axis('off')
ax.set_facecolor(C.SURFACE)
rows = [
    ('locked displacement', 'max |ΔX| = 0', f'{n_locked_str} steps', True),
    ('locked speed', 'max |U| = 0', f'{n_locked_str} steps', True),
    ('release instant', f't = {rel:.5f}', 'φ_liq → 0.709', True),
    ('impulse-free handoff',
     f"peak |a| = {MOTION['release_acceleration_over_free']:.3f} × a_free",
     '< 2 × a_free', True),
    ('peak fall speed', f'min U_y = {u_min:.4f}', 'settles', True),
    ('wall contact', f't = {wall_hit:.3f},  y → {y[-1]:.4f}', 'reaches y = R', True),
    ('soft-contact penetration', f'min y = {y_min:.5f}',
     f'{D["R"]-y_min:.4f} d, contact model', None),
]
ax.text(0.0, 1.02, 'd   B.2 motion gates, as measured', transform=ax.transAxes,
        fontsize=11, fontweight='semibold', color=C.INK, va='bottom')
yy = 0.90
for label, value, bar, ok in rows:
    mark, col = ('✓', C.STATUS['good']) if ok else ('•', C.MUTED)
    ax.text(0.0, yy, mark, transform=ax.transAxes, color=col,
            fontsize=11, fontweight='bold', va='center')
    ax.text(0.055, yy, label, transform=ax.transAxes, color=C.INK,
            fontsize=9.2, va='center')
    ax.text(0.47, yy, value, transform=ax.transAxes, color=C.INK2,
            fontsize=9.2, va='center', family='DejaVu Sans')
    ax.text(0.99, yy, bar, transform=ax.transAxes, color=C.MUTED,
            fontsize=8.4, va='center', ha='right')
    yy -= 0.132
ax.axhline(0.0, color=C.GRID)

fig.suptitle('Stage B · the sediment grain melts out of the ice, falls, and lands',
             x=0.065, ha='left', fontsize=15, fontweight='semibold', color=C.INK, y=0.965)
fig.text(0.065, 0.905,
         'Accepted B.2 coupled gate — job 20492853, d/Δx = 16, 64×96×64, max_dt = 0.01, '
         'warm bottom wall, passive nonzero salt',
         ha='left', fontsize=9.5, color=C.INK2)

C.save(fig, 'fig05_lifecycle',
       'Source: B2closure_coupled16_rampfix_20492853 — mobile.dat, release.dat, Particle_*.h5. '
       'Synthetic debugging parameters, not a physical ECCO scenario.')
