#!/usr/bin/env python3
"""Figure 1 — Stage B at a glance: what was built, what passed, what B.3 inherits.

Top      the four Stage-B blocks and their status
Bottom   left, the six B.2 acceptance gates with their measured values;
         right, the constraint the closure audit hands to B.3 — which
         quantities are trustworthy at max_dt = 0.01 and which are not.

Every value shown here is quoted from the closure report / roadmap through
``stageb_common.ROADMAP_NUMBERS``; the figures that follow re-derive them from
the run directories.
"""
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stageb_common as C

plt = C.use_style()
N = C.ROADMAP_NUMBERS

# gate 1's step count is measured here rather than quoted, so it cannot drift
_t, _X, _U = C.mobile(C.run('coupled16'))
_locked = _t < C.release_time(C.run('coupled16'))
GATES = list(N['gates'])
_nsteps = f'{int(_locked.sum()):,}'.replace(',', ' ')
_su = next(j['cpu_hours'] for j in C.campaign_jobs()['jobs']
           if j['JobID'] == '20492853')
GATES[0] = (GATES[0][0],
            f'max |ΔX| = 0, max |U| = 0 over all {_nsteps} locked steps',
            GATES[0][2], GATES[0][3])

fig = plt.figure(figsize=(13.0, 7.6))
gs = fig.add_gridspec(2, 2, height_ratios=[0.78, 1.0], hspace=0.30, wspace=0.16,
                      left=0.035, right=0.975, top=0.845, bottom=0.075)

# --- top: the four blocks ------------------------------------------------
ax = fig.add_subplot(gs[0, :])
ax.axis('off')
ax.set_xlim(0, 100)
ax.set_ylim(0, 26)

BLOCKS = [
    ('B.0', 'Non-invasiveness gate',
     'Does turning VOF_IBM and LAG_PARTICLE_RESOLVED on disturb the '
     'validated Stage-A physics?',
     'PASSED — worst residual 9.2e-9 in ice volume at t = 10, six orders '
     'inside the already-accepted binary-lineage bar',
     C.STATUS['good'], 'done'),
    ('B.1', 'Code and inputs',
     'Polar EOS re-fit · resolved IBM beside VOF · hold-in-ice release · '
     'meltwater tracer · ECCO inits · mixing diagnostics',
     'DONE — the sediment-clip leak took four rounds; resolved by a no-flux '
     'CH mask on the grain. F_release calibrated to 0.709',
     C.STATUS['good'], 'done'),
    ('B.2', '3-D shake-out',
     'Exercise the whole coupled path cheaply: melt, lock, release, settle, '
     'collide, conserve — on a synthetic 4 × 6 × 4 box',
     'CLOSED 2026-09-08 — all six gates pass on job 20492853; the escape '
     'TIMING is documented as not converging',
     C.STATUS['good'], 'done'),
    ('B.3', 'ECCO production',
     'Choose the physical scenario, fit its EOS, set the resolution, cost it, '
     'and run the two precursors',
     'NEXT — held for the user’s physical choices and approval. No real ECCO '
     'scenario has been submitted or fitted',
     C.ORANGE, 'next'),
]

W, GAP = 22.5, 3.3
WRAP = 38
for i, (tag, title, what, result, col, state) in enumerate(BLOCKS):
    x = i * (W + GAP)
    fc = '#f2f7f2' if state == 'done' else '#fdf3e9'
    ax.add_patch(plt.Rectangle((x, 1.0), W, 23.5, facecolor=fc, edgecolor=col,
                               lw=1.6, joinstyle='round'))
    ax.add_patch(plt.Rectangle((x, 21.0), W, 3.5, facecolor=col, edgecolor=col, lw=0))
    ax.text(x + 1.0, 22.7, f'{tag}   {title}', color='white', fontsize=10.5,
            fontweight='semibold', va='center')
    ax.text(x + 1.0, 19.6, textwrap.fill(what, WRAP), color=C.INK2, fontsize=8.0,
            va='top', linespacing=1.55)
    ax.text(x + 1.0, 11.4, textwrap.fill(result, WRAP),
            color=col if state == 'next' else C.INK,
            fontsize=8.1, va='top', linespacing=1.55, fontweight='semibold')
    if i < len(BLOCKS) - 1:
        ax.annotate('', xy=(x + W + GAP - 0.5, 13), xytext=(x + W + 0.5, 13),
                    arrowprops=dict(arrowstyle='-|>', color=C.MUTED, lw=1.6))

# --- bottom left: the six gates -----------------------------------------
ax = fig.add_subplot(gs[1, 0])
ax.axis('off')
ax.text(0.0, 1.0, 'The six B.2 acceptance gates, as measured',
        transform=ax.transAxes, fontsize=11.5, fontweight='semibold',
        color=C.INK, va='bottom')
ax.text(0.0, 0.955, f'job 20492853 · d/Δx = 16 · max_dt = 0.01 · 64 ranks · {_su:.0f} SU',
        transform=ax.transAxes, fontsize=8.5, color=C.MUTED, va='bottom')
yy = 0.86
for name, measured, bar, _ in GATES:
    ax.text(0.0, yy, '✓', transform=ax.transAxes, color=C.STATUS['good'],
            fontsize=12, fontweight='bold', va='center')
    ax.text(0.048, yy, name, transform=ax.transAxes, color=C.INK, fontsize=9.6,
            va='center')
    ax.text(0.99, yy, bar, transform=ax.transAxes, color=C.MUTED, fontsize=8.3,
            va='center', ha='right')
    ax.text(0.048, yy - 0.062, measured, transform=ax.transAxes, color=C.INK2,
            fontsize=8.8, va='center')
    yy -= 0.155

# --- bottom right: what B.3 inherits ------------------------------------
ax = fig.add_subplot(gs[1, 1])
ax.axis('off')
ax.text(0.0, 1.0, 'What B.3 inherits from the closure audit',
        transform=ax.transAxes, fontsize=11.5, fontweight='semibold',
        color=C.INK, va='bottom')
ax.text(0.0, 0.955, 'the constraint is on WHICH quantities are converged, not on whether '
        'the solver is',
        transform=ax.transAxes, fontsize=8.5, color=C.MUTED, va='bottom')

GOOD = [('conservation budgets', 'machine precision'),
        ('pre-release fields', 'order ≈ 0.9 in u, θ and $C_L$'),
        ('MPI reproducibility', 'bit-identical, 64 vs 32 ranks'),
        ('ice-free settling velocity', '0.08–0.28 % at max_dt = 0.01'),
        ('stability at max_dt = 0.01', 'established by completion')]
BAD = [('absolute release time', 'drifts 5.75 → 6.63 → 7.42 as dt halves'),
       ('escape duration', 'observed order 0.15, where healthy is 1'),
       ('release-to-contact timing', 'inherits the same one-signed bias')]

ax.add_patch(plt.Rectangle((0.0, 0.50), 1.0, 0.40, transform=ax.transAxes,
                           facecolor='#f2f7f2', edgecolor=C.STATUS['good'], lw=1.2))
ax.text(0.018, 0.865, 'TRUSTWORTHY at max_dt = 0.01', transform=ax.transAxes,
        color=C.STATUS['good'], fontsize=9, fontweight='semibold', va='center')
yy = 0.805
for k, v in GOOD:
    ax.text(0.03, yy, '✓', transform=ax.transAxes, color=C.STATUS['good'],
            fontsize=10, va='center')
    ax.text(0.068, yy, k, transform=ax.transAxes, color=C.INK, fontsize=9, va='center')
    ax.text(0.98, yy, v, transform=ax.transAxes, color=C.INK2, fontsize=8.4,
            va='center', ha='right')
    yy -= 0.064

ax.add_patch(plt.Rectangle((0.0, 0.13), 1.0, 0.31, transform=ax.transAxes,
                           facecolor='#fbe7e7', edgecolor=C.STATUS['critical'], lw=1.2))
ax.text(0.018, 0.405, 'NOT TRUSTWORTHY at any affordable timestep',
        transform=ax.transAxes, color=C.STATUS['critical'], fontsize=9,
        fontweight='semibold', va='center')
yy = 0.345
for k, v in BAD:
    ax.text(0.03, yy, '✕', transform=ax.transAxes, color=C.STATUS['critical'],
            fontsize=10, va='center')
    ax.text(0.068, yy, k, transform=ax.transAxes, color=C.INK, fontsize=9, va='center')
    ax.text(0.98, yy, v, transform=ax.transAxes, color=C.INK2, fontsize=8.4,
            va='center', ha='right')
    yy -= 0.064

ax.text(0.0, 0.085,
        'So B.3 must not report or depend on release timing, and should prefer grains '
        'initialized already free.\n'
        'It still needs the user’s physical choices, then its own EOS fit, '
        'dimensionless parameters, resolution and cost.',
        transform=ax.transAxes, color=C.INK2, fontsize=8.6, va='top', linespacing=1.6)

fig.suptitle('STAGE B — ECCO sediment-from-ice study',
             x=0.035, ha='left', fontsize=18, fontweight='semibold',
             color=C.INK, y=0.972)
fig.text(0.035, 0.895,
         'Can PARTIES follow a sediment grain from frozen inside sea ice, through the melt '
         'that frees it, to the sea floor — conservatively, reproducibly, and without '
         'disturbing the Stage-A physics?',
         ha='left', fontsize=10, color=C.INK2)

C.save(fig, 'fig01_stageb_map',
       'Values quoted from YANG_ECCO_B2_CLOSURE_REPORT.md and the Stage-B sections of '
       'YANG_ECCO_IMPLEMENTATION_ROADMAP.md; the following figures re-derive them from '
       'the run directories on /anvil/scratch/x-mjalabert/ECCO_StageB.')
