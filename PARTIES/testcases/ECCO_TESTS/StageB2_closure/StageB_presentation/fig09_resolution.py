#!/usr/bin/env python3
"""Figure 9 — settling: grid refinement, and isolating the sediment/ice coupling.

Panels
  (a) ice-free settling transient at d/Δx = 16, 24, 32
  (b) the difference between successive resolutions — it shrinks, which is what
      supports the d/Δx >= 24 choice for the production grid
  (c) the Darcy-penalization control: a stiff tau and an effectively absent one
      give bit-identical particle and fluid fields, and a conventional IBM run
      agrees to roundoff — so the ice penalization is not driving the settling
  (d) what is and is not claimed about terminal velocity

Sources: mobile.dat of settle16/24/32 and the recorded settle16_controls.json.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stageb_common as C

plt = C.use_style()

ARMS = [('16', 'settle16', C.BLUE, '20460632'),
        ('24', 'settle24', C.ORANGE, '20477406'),
        ('32', 'settle32', C.AQUA, '20477407')]
traces = {}
for lab, key, col, job in ARMS:
    t, X, U = C.mobile(C.run(key))
    traces[lab] = (t, X[:, 1], U[:, 1], col, job)

ref = C.audit('settling_refinement.json')
ctrl = C.audit('settle16_controls.json')

fig = plt.figure(figsize=(12.4, 7.2))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.24,
                      left=0.068, right=0.985, top=0.845, bottom=0.135)

TEND = 6.0

# --- (a) settling transient ---------------------------------------------
ax = C.finish(fig.add_subplot(gs[0, 0]))
for lab, _, _, job in ARMS:
    t, y, uy, col, _ = traces[lab]
    m = t <= TEND
    ax.plot(t[m], -uy[m], color=col, lw=2.0, label=f'd/Δx = {lab}')
ax.plot([6], [-ref['16']['t6_Uy']], 'o', ms=5, color=C.BLUE, mec=C.SURFACE, mew=1.2)
ax.plot([6], [-ref['24']['t6_Uy']], 'o', ms=5, color=C.ORANGE, mec=C.SURFACE, mew=1.2)
ax.plot([6], [-ref['32']['t6_Uy']], 'o', ms=5, color=C.AQUA, mec=C.SURFACE, mew=1.2)
ax.set_xlim(0, TEND)
ax.set_xlabel('time  $t$')
ax.set_ylabel('downward speed  $-U_y$')
ax.set_title('a   Ice-free settling, resolved-IBM sphere')
ax.legend(loc='lower right')
ax.text(0.15, 1.18,
        'the short column reaches the bottom before a\n'
        'clean terminal plateau — these are transients',
        color=C.MUTED, fontsize=8, va='top')

# --- (b) refinement differences -----------------------------------------
ax = C.finish(fig.add_subplot(gs[0, 1]))
tt = np.linspace(1.0, TEND, 1200)
curves = {lab: np.interp(tt, traces[lab][0], -traces[lab][2]) for lab, _, _, _ in ARMS}
for ((a, b), col, dy) in [(('16', '24'), C.VIOLET, -16), (('24', '32'), C.STATUS['good'], 9)]:
    d = 100 * np.abs(curves[a] - curves[b]) / np.abs(curves[b])
    ax.plot(tt, d, color=col, lw=2.0, label=f'{a} → {b}')
    ax.annotate(f'{d[-1]:.2f} % at t = 6', (tt[-1], d[-1]),
                textcoords='offset points', xytext=(-8, dy), ha='right',
                fontsize=8.5, color=col)
ax.set_xlim(1.0, TEND)
ax.set_xlabel('time  $t$')
ax.set_ylabel('relative difference in $U_y$  (%)')
ax.set_title('b   Refinement shrinks the difference — the basis for d/Δx ≥ 24')
ax.legend(loc='upper right', title='grid pair', title_fontsize=8.5)

# --- (c) isolation controls ---------------------------------------------
ax = C.finish(fig.add_subplot(gs[1, 0]))
ax.set_xscale('log')
keys = ['X', 'U', 'u', 'v', 'w', 'p']
labels = ['particle position', 'particle velocity', 'fluid $u$', 'fluid $v$',
          'fluid $w$', 'pressure']
FLOOR = 1e-16
refr = [max(ctrl['reference']['max_abs_diff'][k], FLOOR) for k in keys]
ypos = np.arange(len(keys))[::-1]
ax.barh(ypos, refr, height=0.5, color=C.BLUE, zorder=3)
for yv, k in zip(ypos, keys):
    ax.text(max(ctrl['reference']['max_abs_diff'][k], FLOOR) * 1.7, yv,
            f'{ctrl["reference"]["max_abs_diff"][k]:.1e}',
            va='center', fontsize=8.5, color=C.BLUE)
ax.set_yticks(ypos)
ax.set_yticklabels(labels)
ax.set_xlim(FLOOR, 3e-11)
ax.set_ylim(-2.4, len(keys) - 0.4)
ax.set_xlabel('max |difference|  vs  a conventional-IBM run, at the common saved times')
ax.set_title('c   The ice penalization is not what is driving the settling')
ax.grid(axis='y', visible=False)
ax.text(1.4e-16, -2.30,
        'And the Darcy stiffness itself does nothing here: the $\\tau = 10^{-3}$ arm and an\n'
        'arm with an effectively absent Darcy force are BIT-IDENTICAL in all six\n'
        'quantities — difference exactly 0  (jobs 20460632 / 20460634).',
        fontsize=8.4, color=C.STATUS['good'], va='bottom', linespacing=1.6)

# --- (d) what is and is not claimed -------------------------------------
ax = fig.add_subplot(gs[1, 1])
ax.axis('off')
ax.text(0.0, 1.02, 'd   What this establishes, and what it does not',
        transform=ax.transAxes, fontsize=11, fontweight='semibold',
        color=C.INK, va='bottom')
claims = [
    (True, 'Transient settling converges under refinement',
     f"{100*ref['24_to_32']['relative_Uy_t6']:.3f} % at t = 6; "
     f"≤ {100*ref['24_to_32']['max_relative_Uy_t1_to_t6']:.3f} % over t = 1–6"),
    (True, 'The ice-penalization path does not perturb it',
     'bit-identical to the negligible-Darcy arm; roundoff vs conventional IBM'),
    (True, 'Time integration of the settling is converged',
     '0.08 % in $U_y$ between max_dt 0.01 and 0.005 (figure 8)'),
    (False, 'NOT a terminal-velocity measurement',
     'the column reaches the bottom before a clean plateau'),
    (False, 'NOT validated against Schiller–Naumann',
     'the unbounded estimate is 1.478; evaluating Cd at the\nmeasured confined speed is not an independent reference'),
    (False, 'NOT a confinement correction',
     'a Stokes tube-wall factor cannot validate this finite-Re periodic box'),
]
yy = 0.90
for ok, head, body in claims:
    mark, col = ('✓', C.STATUS['good']) if ok else ('✕', C.STATUS['critical'])
    ax.text(0.0, yy, mark, transform=ax.transAxes, color=col, fontsize=11,
            fontweight='bold', va='center')
    ax.text(0.05, yy, head, transform=ax.transAxes, color=C.INK, fontsize=9.4,
            va='center')
    ax.text(0.05, yy - 0.058, body, transform=ax.transAxes, color=C.INK2,
            fontsize=8.5, va='top', linespacing=1.4)
    yy -= 0.165 + 0.075 * body.count('\n')

fig.suptitle('Stage B · how fast does a freed grain actually fall, and how sure are we',
             x=0.068, ha='left', fontsize=15, fontweight='semibold', color=C.INK, y=0.965)
fig.text(0.068, 0.895,
         'Genuine ice-free settling decks (ice layer moved outside the box, both EOS '
         'coefficients zero) — the legacy “ice-free” deck was not ice-free and is not used here.',
         ha='left', fontsize=9.5, color=C.INK2)

C.save(fig, 'fig09_resolution',
       'Sources: mobile.dat of B2closure_settle16_20460632 / settle24_20477406 / '
       'settle32_20477407; B2_closure/settling_refinement.json and settle16_controls.json '
       '(jobs 20460632 / 20460634 / 20460708, all recorded TIMEOUT after these observations).')
