#!/usr/bin/env python3
"""Figure 8 — the one thing B.2 could NOT converge, and why.

The matched-ramp timestep triple (max_dt = 0.01 / 0.005 / 0.0025, ramp duration
held fixed) fails its debugging tolerances.  Investigating it found the cause:
while the grain creeps out of the ice it is almost exactly supported, so the net
driving force is a near-cancellation and its timing drifts instead of converging.

Panels
  (a) release-relative fall: the trajectories do not collapse as dt halves
  (b) escape arrival time at y = 3.30 against dt — a drift, not a limit
  (c) the cancellation itself: IBM reaction against the buoyant weight
  (d) free-fall convergence, isolated three ways: the settling dynamics ARE
      converged; the hard ice-penalization threshold is what costs the accuracy

Everything is recomputed with B2_closure/refine_timesteps.py on mobile.dat.
"""
import importlib.util
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stageb_common as C

plt = C.use_style()
N = C.ROADMAP_NUMBERS

_spec = importlib.util.spec_from_file_location('rt', C.B2 / 'refine_timesteps.py')
rt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rt)

WEIGHT = N['escape']['buoyant_weight']
DTS = [0.01, 0.005, 0.0025]
TRIPLE = [C.run(k) for k in ('coupled16', 'coupled16_halfdt', 'coupled16_quarterdt')]
data = [rt.read(r) for r in TRIPLE]
pairs = [rt.pair(data[i], data[i + 1]) for i in range(2)]
tau = [float(np.interp(rt.ESCAPE_Y, d['y'], d['tau'])) for d in data]
p_obs = rt.order(pairs[0]['escape_lag'], pairs[1]['escape_lag'])


def f_ibm_at(rundir, t_target=34.0):
    best = None
    for _, path in C.snapshots(rundir):
        ppath = path.with_name(path.name.replace('Data_', 'Particle_'))
        with h5py.File(ppath, 'r') as p:
            t = float(p['time'][0])
            if best is None or abs(t - t_target) < abs(best[0] - t_target):
                best = (t, float(p['mobile/F_IBM'][...].reshape(-1, 3)[0][1]))
    return best


fibm = [f_ibm_at(r)[1] for r in TRIPLE]
net = [WEIGHT - f for f in fibm]


def uy_of_height(rundir):
    a = np.loadtxt(Path(rundir) / 'mobile.dat', delimiter=',')
    return a[:, 3][::-1], a[:, 6][::-1]


def peak_relative(run_a, run_b, y0=None, y1=None):
    """max |ΔU_y(y)| / peak |U_y| over the common free-fall height window."""
    ya, va = uy_of_height(run_a)
    yb, vb = uy_of_height(run_b)
    lo = max(ya.min(), yb.min()) + 1e-6 if y0 is None else y0
    hi = min(ya.max(), yb.max()) - 1e-6 if y1 is None else y1
    y = np.linspace(lo, hi, 3001)
    fa, fb = np.interp(y, ya, va), np.interp(y, yb, vb)
    return float(np.max(np.abs(fa - fb)) / max(abs(fa).max(), abs(fb).max()))


CONFIGS = [
    ('ice-free settling\nthreshold never fires',
     peak_relative(C.run('settle24_dt10'), C.run('settle24_dt05')),
     C.STATUS['good'], 'jobs 20500848 / 20500849'),
    ('coupled\nthreshold DISABLED',
     rt.pair(rt.read(C.run('coupled16_smooth')),
             rt.read(C.run('coupled16_smooth_halfdt')))['open_water_Uy_max_relative'],
     C.BLUE, 'jobs 20500720 / 20500721'),
    ('coupled\nhard threshold ACTIVE',
     pairs[0]['open_water_Uy_max_relative'],
     C.STATUS['warning'], 'jobs 20492853 / 20492854'),
]

fig = plt.figure(figsize=(12.6, 7.6))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.30,
                      left=0.068, right=0.985, top=0.845, bottom=0.10)

# --- (a) release-relative trajectories -----------------------------------
ax = C.finish(fig.add_subplot(gs[0, 0]))
for d, dt, col in zip(data, DTS, [C.BLUE, C.ORANGE, C.VIOLET]):
    ax.plot(d['tau'][::-1], d['y'][::-1], color=col, lw=2.0,
            label=f'max_dt = {dt}')
ax.axhline(rt.ESCAPE_Y, color=C.MUTED, lw=1.1, ls=':')
ax.text(11.6, rt.ESCAPE_Y + 0.1, 'reference height  y = 3.30', color=C.MUTED,
        fontsize=8, ha='right')
for d, tt, col in zip(data, tau, [C.BLUE, C.ORANGE, C.VIOLET]):
    ax.plot([tt], [rt.ESCAPE_Y], 'o', ms=6, color=col, mec=C.SURFACE, mew=1.4, zorder=5)
ax.annotate('', xy=(tau[0], 2.55), xytext=(tau[2], 2.55),
            arrowprops=dict(arrowstyle='<->', color=C.STATUS['critical'], lw=1.3))
ax.text(np.mean(tau), 2.40, f'{tau[2]-tau[0]:.2f} time units apart —\nand still moving',
        color=C.STATUS['critical'], fontsize=8.5, ha='center', va='top')
ax.set_xlim(0, 11.6)
ax.set_ylim(0.4, 4.3)
ax.set_xlabel(r'time since release  $\tau = t - t_{release}$')
ax.set_ylabel('grain centre height  $y$')
ax.set_title('a   Halving the timestep does not collapse the escape')
ax.legend(loc='lower left')

# --- (b) escape arrival vs dt --------------------------------------------
ax = C.finish(fig.add_subplot(gs[0, 1]))
ax.set_xscale('log')
ax.plot(DTS, tau, 'o-', color=C.STATUS['critical'], ms=8, mec=C.SURFACE, mew=1.5,
        label='measured')
for dt, tt in zip(DTS, tau):
    ax.annotate(f'{tt:.2f}', (dt, tt), textcoords='offset points', xytext=(0, 11),
                ha='center', fontsize=8.5, color=C.INK2)
# what a first-order scheme would look like, anchored on the finest point
first = tau[2] + (np.array(DTS) - DTS[2]) * (tau[1] - tau[2]) / (DTS[1] - DTS[2])
ax.plot(DTS, first, ls=(0, (4, 3)), color=C.MUTED, lw=1.6,
        label='what first order would look like')
ax.set_xticks(DTS)
ax.set_xticklabels([str(d) for d in DTS])
ax.minorticks_off()
ax.invert_xaxis()
ax.set_ylim(5.3, 7.95)
ax.set_xlabel('max_dt')
ax.set_ylabel(r'time from release to $y = 3.30$')
ax.set_title('b   The escape time drifts; it does not converge')
ax.text(0.0099, 7.88,
        f'successive increments  {pairs[0]["escape_lag"]:.3f} → {pairs[1]["escape_lag"]:.3f}\n'
        f'ratio {pairs[1]["escape_lag"]/pairs[0]["escape_lag"]:.2f};  '
        f'observed order {p_obs:.2f}  (healthy = 1)',
        color=C.STATUS['critical'], fontsize=8.5, ha='left', va='top')
ax.legend(loc='lower right')

# --- (c) the near-cancellation -------------------------------------------
ax = C.finish(fig.add_subplot(gs[1, 0]))
xs = np.arange(3)
ax.bar(xs, fibm, width=0.46, color=C.BLUE, zorder=3)
ax.bar(xs, net, width=0.46, bottom=fibm, color=C.STATUS['critical'], zorder=3)
ax.axhline(WEIGHT, color=C.MUTED, lw=1.3, ls='--', zorder=4)
ax.text(0.5, 0.762, f'buoyant weight {WEIGHT:.4f}', color=C.INK2,
        fontsize=8, ha='center', va='top')
for x, f, n in zip(xs, fibm, net):
    ax.text(x, WEIGHT + 0.055, f'net {n:.5f}\n{100*n/WEIGHT:.1f} % of the weight',
            ha='center', fontsize=8.3, color=C.STATUS['critical'])
    ax.text(x, f / 2, f'IBM reaction\n{f:.5f}', ha='center', va='center',
            color='white', fontsize=8.8, linespacing=1.5)
ax.set_xticks(xs)
ax.set_xticklabels([f'max_dt = {d}' for d in DTS])
ax.set_ylim(0, 1.02)
ax.set_ylabel('vertical force at  $t = 34$')
ax.set_title('c   A 2 % change in the reaction is a 20× change in the net')
ax.grid(axis='x', visible=False)

# --- (d) free-fall convergence, isolated ---------------------------------
ax = C.finish(fig.add_subplot(gs[1, 1]))
vals = [100 * v for _, v, _, _ in CONFIGS]
cols = [c for _, _, c, _ in CONFIGS]
labs = [l for l, _, _, _ in CONFIGS]
bars = ax.barh(np.arange(3)[::-1], vals, height=0.5, color=cols, zorder=3)
for i, (v, (lab, _, col, src)) in enumerate(zip(vals, CONFIGS)):
    ax.text(v + 0.09, 2 - i, f'{v:.2f} %', va='center', fontsize=10,
            color=col, fontweight='semibold')
    ax.text(v + 0.09, 2 - i - 0.28, src, va='center', fontsize=7.5, color=C.MUTED)
ax.set_yticks(np.arange(3)[::-1])
ax.set_yticklabels(labs, fontsize=9)
ax.set_xlim(0, 4.4)
ax.set_xlabel(r'max $|\Delta U_y|$ / peak $|U_y|$  between max_dt 0.01 and 0.005  (%)')
ax.set_title('d   The settling itself is converged — the threshold is the cost')
ax.grid(axis='y', visible=False)
ax.annotate('', xy=(vals[2], 0.38), xytext=(vals[1], 0.38),
            arrowprops=dict(arrowstyle='<->', color=C.INK2, lw=1.2))
ax.text(0.5 * (vals[1] + vals[2]), 0.50, '≈ 8×', ha='center', color=C.INK2, fontsize=9.5)

fig.suptitle('Stage B · the one quantity that does not converge, and the measurements that explain it',
             x=0.068, ha='left', fontsize=15, fontweight='semibold', color=C.INK, y=0.965)
fig.text(0.068, 0.895,
         'Matched-ramp timestep triple at d/Δx = 16. Chaos, Darcy leakage and a dt-dependent '
         'melt rate were each tested and refuted before conditioning was accepted as the cause.',
         ha='left', fontsize=9.5, color=C.INK2)

C.save(fig, 'fig08_conditioning',
       'Recomputed with B2_closure/refine_timesteps.py from mobile.dat / release.dat; forces from '
       'Particle_*.h5 at t ≈ 34. Panel d uses one metric for all three arms (peak-normalised '
       'U_y(y)); the closure report quotes 0.28 % for the ice-free arm under a pointwise-relative '
       'metric from t ≥ 0.5 — the ordering and the ≈8× threshold penalty are the same either way.')
