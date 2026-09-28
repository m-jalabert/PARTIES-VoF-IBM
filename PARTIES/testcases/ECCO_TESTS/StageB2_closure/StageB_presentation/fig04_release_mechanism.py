#!/usr/bin/env python3
"""Figure 4 — the hold-in-ice / interface-triggered release mechanism (B.1.3, B.1.4).

Panels
  (a) the mechanism, drawn: shell-averaged liquid fraction phi_liq over the grain
  (b) phi_liq(t) measured in the accepted 3-D gate, with the F_release threshold,
      the saturation ceiling, and the two traps the calibration found
  (c) F_release calibration: phi_liq against how far the melt front has passed
  (d) the meltwater-tracer budget before and after the sediment-clip fix

Sources: coupled16_rampfix release.dat (a, b);  roadmap B.1.3-fn calibration
table and B.1.4 discriminator runs (c, d) via stageb_common.ROADMAP_NUMBERS.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stageb_common as C

plt = C.use_style()
N = C.ROADMAP_NUMBERS

RUN = C.run('coupled16')
r = C.release(RUN)
rel = C.release_time(RUN)
F_REL = C.deck(RUN)['F_release'][0]
CEIL = N['phi_liq_ceiling']

fig = plt.figure(figsize=(12.4, 7.6))
gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.24,
                      left=0.065, right=0.985, top=0.855, bottom=0.085)

# --- (a) the mechanism, drawn -------------------------------------------
ax = fig.add_subplot(gs[0, 0])
ax.set_aspect('equal')
ax.axis('off')
ax.set_xlim(-0.3, 10.2)
ax.set_ylim(-0.2, 3.6)
ax.set_title('a   How the grain knows it is free', y=1.04)

stages = [(1.3, 0.12, 'front below the grain\nφ_liq small — locked'),
          (4.9, 0.62, 'front at the centre\nφ_liq ≈ 0.57 — locked'),
          (8.5, 1.35, 'front just cleared it\nφ_liq = 0.709 — release')]
for cx, front, label in stages:
    # ice block above the front, water below
    ax.add_patch(plt.Rectangle((cx - 1.45, front), 2.9, 2.45 - front,
                               facecolor='#dcebf8', edgecolor=C.AQUA, lw=1.2))
    ax.add_patch(plt.Rectangle((cx - 1.45, 0.0), 2.9, front,
                               facecolor='#fbf4ec', edgecolor=C.AXIS, lw=0.8))
    # release shell (2 cells thick) then the grain
    ax.add_patch(plt.Circle((cx, 1.25), 0.82, facecolor='none',
                            edgecolor=C.ORANGE, lw=1.3, ls=(0, (3, 2))))
    ax.add_patch(plt.Circle((cx, 1.25), 0.6, facecolor='#6b6259',
                            edgecolor='#443e38', lw=1.0))
    ax.text(cx, -0.12, label, ha='center', va='top', fontsize=8, color=C.INK2)
ax.text(1.3, 2.26, 'ice', color=C.AQUA, fontsize=9, ha='center')
ax.text(8.5, 0.33, 'warm water', color='#a07c55', fontsize=8.5, ha='center', va='center')
ax.annotate('shell = grain support + 2 cells;\nφ_liq is the liquid fraction averaged over it',
            xy=(8.5 + 0.82, 1.9), xytext=(5.4, 3.25), fontsize=8, color=C.ORANGE,
            ha='center', arrowprops=dict(arrowstyle='-', color=C.ORANGE, lw=1.0))

# --- (b) phi_liq(t) ------------------------------------------------------
ax = C.finish(fig.add_subplot(gs[0, 1]))
ax.axhspan(CEIL, 1.0, color='#fbe7e7', lw=0, zorder=0)
ax.text(1, 0.995, f'F_release ≥ 0.98 never fires — φ_liq saturates at ≈ {CEIL}',
        color=C.STATUS['critical'], fontsize=8, va='top')
ax.axhline(F_REL, color=C.ORANGE, lw=1.4, ls='--', zorder=2)
ax.text(1, F_REL - 0.03, f'F_release = {F_REL}  (calibrated)', color=C.ORANGE, fontsize=8.5)
ax.axhline(0.9, color=C.MUTED, lw=1.1, ls=':', zorder=2)
ax.text(64, 0.885, 'the original 0.9 fires 0.67 d late', color=C.MUTED,
        fontsize=8, ha='right', va='top')
ax.plot(r['time'], r['phi_liq'], color=C.BLUE, zorder=4)
ax.axvline(rel, color=C.ORANGE, lw=1.1, ls='--', zorder=3)
ax.plot([rel], [F_REL], 'o', ms=7, color=C.ORANGE, mec=C.SURFACE, mew=1.6, zorder=6)
ax.annotate(f'release  t = {rel:.2f}', xy=(rel, F_REL), xytext=(rel + 4, 0.40),
            color=C.ORANGE, fontsize=8.5,
            arrowprops=dict(arrowstyle='-', color=C.ORANGE, lw=1.1))
ax.set_ylim(0, 1.03)
ax.set_xlim(0, r['time'].max())
ax.set_xlabel('time  $t$')
ax.set_ylabel(r'shell liquid fraction  $\varphi_{liq}$')
ax.set_title('b   The trigger, measured in 3-D')

# --- (c) F_release calibration ------------------------------------------
ax = C.finish(fig.add_subplot(gs[1, 0]))
cal = N['F_release_cal']
xs = np.array(sorted(cal))
ys = np.array([cal[k] for k in xs])
ax.axhspan(CEIL, 1.02, color='#fbe7e7', lw=0, zorder=0)
ax.plot(xs, ys, 'o-', color=C.BLUE, ms=7, mec=C.SURFACE, mew=1.5, zorder=4)
for x, y in zip(xs, ys):
    ax.annotate(f'{y:.3f}', (x, y), textcoords='offset points', xytext=(0, 9),
                ha='center', fontsize=8, color=C.INK2)
ax.axvline(1.0, color=C.ORANGE, lw=1.3, ls='--', zorder=2)
ax.axhline(F_REL, color=C.ORANGE, lw=1.3, ls='--', zorder=2)
ax.text(1.04, 0.545, 'mechanically free:\nthe front has just\ncleared the grain',
        color=C.ORANGE, fontsize=8.5, va='center')
ax.text(1.80, 1.045, f'ceiling {CEIL}: the shell always holds part of the\n'
        'diffuse band and the non-melting sediment rim',
        color=C.STATUS['critical'], fontsize=7.8, ha='right', va='top')
ax.set_xlabel('melt-front exposure  (grain diameters past the grain underside)')
ax.set_ylabel(r'$\varphi_{liq}$ at that exposure')
ax.set_ylim(0.45, 1.09)
ax.set_xlim(0.35, 1.82)
ax.set_title('c   Calibrating F_release, so the grain leaves when it is free')

# --- (d) tracer budget, before and after the clip fix --------------------
ax = C.finish(fig.add_subplot(gs[1, 1]))
ax.set_yscale('log')
leak = N['tracer_leak']
noP = N['tracer_noparticle']
tt = np.array(sorted(leak))
ax.plot(tt, [abs(leak[k]) for k in tt], 'o-', color=C.STATUS['critical'], ms=6,
        label='with a resolved grain — before the fix')
tn = np.array(sorted(noP))
ax.plot(tn, [abs(noP[k]) for k in tn], 's-', color=C.BLUE, ms=6,
        label='no grain (the transport itself)')
ax.axhline(1e-3, color=C.MUTED, lw=1.2, ls=':')
ax.text(30, 1.25e-3, 'budget gate 1e-3', color=C.MUTED, fontsize=8, ha='right')
final = C.audit('coupled16_rampfix_profiles.json')['tracer_max_relative']
ax.axhline(final, color=C.STATUS['good'], lw=2.0)
ax.plot([30], [final], '*', ms=13, color=C.STATUS['good'], mec=C.SURFACE, mew=0.8)
ax.text(2, final * 1.4, f'after the fix, 3-D coupled gate:  {final:.2e}',
        color=C.STATUS['good'], fontsize=8.5)
ax.set_xlabel('time  $t$')
ax.set_ylabel(r'$|\int C_{mw} - \int dF| \;/\; \int dF$')
ax.set_ylim(1e-8, 3e-1)
ax.set_title('d   The 9.3 % tracer error was a liquid-mass leak in the sediment cells')
ax.legend(loc='center right')

fig.suptitle('Stage B · release mechanism and the coupling defect it exposed',
             x=0.065, ha='left', fontsize=15, fontweight='semibold', color=C.INK, y=0.965)
fig.text(0.065, 0.90,
         'B.1.3 hold-in-ice → interface-triggered release;  B.1.4 the uncompensated '
         'C_L clip inside the IBM support, and its fix (CH no-flux sediment mask)',
         ha='left', fontsize=9.5, color=C.INK2)

C.save(fig, 'fig04_release_mechanism',
       'Sources: (a,b) B2closure_coupled16_rampfix_20492853/release.dat; '
       '(c) roadmap B.1.3-fn calibration table; (d) roadmap B.1.4 discriminator runs '
       '19761565 / StageB_release2d, and coupled16_rampfix_profiles.json for the final value.')
