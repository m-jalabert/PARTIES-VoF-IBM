#!/usr/bin/env python3
"""Figure 2 — the configuration that every Stage-B number in this deck comes from.

(a) the box, drawn to scale, with the boundary conditions
(b) the initial profiles, read back from the t = 0 snapshot
(c) the deck, read back from parties.inp
(d) which of those numbers are synthetic debugging choices and which are not

The point of (d) is that nothing in B.2 is a physical ECCO scenario.  The
parameters were chosen so the melt front clears the grain inside an affordable
run, and the salt field is deliberately passive (eos_betaS = 0) so it exercises
conservation without pretending to be the polar equation of state.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stageb_common as C

plt = C.use_style()

RUN = C.run('coupled16')
D = C.COUPLED16
deck = C.deck(RUN)
f0 = C.load_fields(C.snapshots(RUN)[0][1])
X0, _, R, _ = C.particle_state(str(C.snapshots(RUN)[0][1]).replace('Data_', 'Particle_'))

fig = plt.figure(figsize=(12.8, 7.4))
gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 0.85, 1.25],
                      hspace=0.42, wspace=0.30,
                      left=0.045, right=0.985, top=0.845, bottom=0.085)

# --- (a) the box ---------------------------------------------------------
ax = fig.add_subplot(gs[:, 0])
ax.set_aspect('equal')
ax.axis('off')
ax.set_xlim(-1.5, 5.3)
ax.set_ylim(-0.9, 7.0)
ax.set_title('a   The box, to scale', y=0.985)
ax.add_patch(plt.Rectangle((0, 0), D['Lx'], D['y_ice0'], facecolor='#fbf4ec',
                           edgecolor='none'))
ax.add_patch(plt.Rectangle((0, D['y_ice0']), D['Lx'], D['Ly'] - D['y_ice0'],
                           facecolor='#dcebf8', edgecolor='none'))
ax.add_patch(plt.Rectangle((0, 0), D['Lx'], D['Ly'], facecolor='none',
                           edgecolor=C.AXIS, lw=1.2))
ax.plot([0, D['Lx']], [D['y_ice0']] * 2, color=C.AQUA, lw=1.8)
ax.add_patch(plt.Circle((X0[0], X0[1]), R, facecolor='#4f4840', edgecolor='#2f2b26'))
ax.annotate('sediment grain\n$d$ = 1,  $\\rho_s/\\rho_f$ = 2.5\n'
            'frozen 0.5 $d$ above the front',
            xy=(X0[0], X0[1] - R * 1.05), xytext=(0.62, 2.85),
            fontsize=8.2, ha='left', va='top', color=C.INK2, linespacing=1.5,
            arrowprops=dict(arrowstyle='-', color=C.MUTED, lw=0.9))
ax.text(0.15, 5.4, 'ice', color=C.AQUA, fontsize=11, fontweight='semibold')
ax.text(3.85, 1.15, 'water', color='#a07c55', fontsize=11,
        fontweight='semibold', ha='right')
ax.annotate('', xy=(-0.35, 0), xytext=(-0.35, D['Ly']),
            arrowprops=dict(arrowstyle='<->', color=C.MUTED, lw=1.0))
ax.text(-0.5, D['Ly'] / 2, f"$L_y$ = {D['Ly']:.0f}", rotation=90, ha='right',
        va='center', color=C.MUTED, fontsize=8.5)
ax.annotate('', xy=(0, -0.35), xytext=(D['Lx'], -0.35),
            arrowprops=dict(arrowstyle='<->', color=C.MUTED, lw=1.0))
ax.text(D['Lx'] / 2, -0.5, f"$L_x = L_z$ = {D['Lx']:.0f}", ha='center', va='top',
        color=C.MUTED, fontsize=8.5)
ax.text(D['Lx'] / 2, D['Ly'] + 0.22, 'cold no-slip wall,  θ = 0',
        ha='center', color=C.BLUE, fontsize=9)
ax.text(D['Lx'] / 2, -0.78, 'warm no-slip wall,  θ = 1', ha='center', va='top',
        color=C.ORANGE, fontsize=9)
for xx in (0, D['Lx']):
    ax.text(xx, D['Ly'] * 0.52, 'periodic', rotation=90, ha='center', va='center',
            color=C.MUTED, fontsize=8, zorder=0)
ax.text(D['Lx'] / 2, 0.45, 'x and z periodic,  y walls',
        ha='center', color=C.MUTED, fontsize=8.2)
ax.annotate('', xy=(4.65, 0.5), xytext=(4.65, 4.1),
            arrowprops=dict(arrowstyle='-|>', color=C.STATUS['critical'], lw=1.6))
ax.text(4.8, 2.3, 'the grain falls', rotation=90, va='center',
        color=C.STATUS['critical'], fontsize=8.5)

# --- (b) initial profiles ------------------------------------------------
ax = C.finish(fig.add_subplot(gs[0, 1]))
y = f0['y']
prof = [
    ('liquid fraction $C_L$', f0['VOF/C_L'].mean(axis=(0, 2)), C.BLUE, '-', 3.2),
    ('temperature θ', f0['Conc/0'].mean(axis=(0, 2)), C.ORANGE, (0, (4, 3)), 2.0),
    ('salinity  $s$ ×50', 50 * f0['Conc/1'].mean(axis=(0, 2)), C.AQUA, '-', 2.0),
]
for label, v, col, ls, lw in prof:
    ax.plot(v, y, color=col, lw=lw, ls=ls, label=label)
ax.axhline(D['y_ice0'], color=C.AQUA, lw=1.1, ls='--')
ax.text(1.02, D['y_ice0'] + 0.08, 'ice front', color=C.AQUA, fontsize=8, ha='right')
ax.set_ylim(0, D['Ly'])
ax.set_xlim(-0.05, 1.05)
ax.set_ylabel('$y$')
ax.set_xlabel('value at  $t$ = 0')
ax.set_title('b   Initial state')
ax.legend(loc='lower left', fontsize=8)

# --- (c) what is synthetic ----------------------------------------------
ax = fig.add_subplot(gs[1, 1])
ax.axis('off')
ax.text(0.0, 1.06, 'c   What is synthetic here', transform=ax.transAxes,
        fontsize=11, fontweight='semibold', color=C.INK, va='bottom')
notes = [
    ('Pe and St', 'chosen so the front clears the grain inside\nan affordable run'),
    ('salinity', 'passive: eos_betaS = 0, so it tests conservation,\nnot the polar equation of state'),
    ('d/Δx = 16', 'debug resolution; production wants ≥ 24'),
    ('Cn, τ, max_dt', 'diagnostic values, not dimensional\nproduction defaults'),
]
yy = 0.90
for head, body in notes:
    ax.text(0.0, yy, '•', transform=ax.transAxes, color=C.STATUS['warning'],
            fontsize=13, va='center')
    ax.text(0.05, yy, head, transform=ax.transAxes, color=C.INK, fontsize=9.2,
            va='center', fontweight='semibold')
    ax.text(0.05, yy - 0.085, body, transform=ax.transAxes, color=C.INK2,
            fontsize=8.4, va='top', linespacing=1.5)
    yy -= 0.175 + 0.075 * body.count('\n')

# --- (d) the deck --------------------------------------------------------
ax = fig.add_subplot(gs[:, 2])
ax.axis('off')
ax.text(0.0, 1.0, 'd   The deck, read back from parties.inp',
        transform=ax.transAxes, fontsize=11, fontweight='semibold',
        color=C.INK, va='bottom')


def g(key, i=0, fmt='{:g}'):
    return fmt.format(deck[key][i])


table = [
    ('geometry and grid', None),
    ('domain', f"{g('xmax')} × {g('ymax')} × {g('zmax')}   (x, z periodic)"),
    ('grid', f"{g('NXM')} × {g('NYM')} × {g('NZM')}   →  d/Δx = 16"),
    ('flow', None),
    ('Re', g('Re')),
    ('gravity', 'ĝ = −ŷ,  |g| = 1'),
    ('phase change and interface', None),
    ('Stefan number St', g('stefan')),
    ('Cahn number Cn', g('Cn')),
    ('CH Péclet', g('Pe_CH')),
    ('melt band ε', g('melt_band_eps')),
    ('liquidus slope', g('liquidus_slope')),
    ('Darcy τ (ice rigidity)', g('darcy_tau', fmt='{:g}')),
    ('scalars', None),
    ('Pe  (θ, s, meltwater)',
     f"{g('Pe', 0)}, {g('Pe', 1)}, {g('Pe', 2)}   →  Sc = 7"),
    ('inits', 'conc_init_type = 34 (θ), 35 (s), 0 (tracer)'),
    ('meltwater tracer', 'on  (NConc = 3)'),
    ('equation of state', None),
    ('eos_betaT / eos_betaS', f"{g('eos_betaT')} / {g('eos_betaS')}   — salt is passive"),
    ('eos_q, Tmd0, slope',
     f"{g('eos_q')}, {g('eos_Tmd0')}, {g('eos_Tmd_slope')}"),
    ('particle', None),
    ('density ratio', g('rho_s')),
    ('F_release', f"{g('F_release')}   (calibrated, §B.1.3-fn)"),
    ('release shell / ramp',
     f"{g('release_shell_cells')} cells / {g('release_ramp_steps')} steps"),
    ('time', None),
    ('max_dt', g('max_dt')),
    ('CFL / time_max', f"{g('cfl')} / {g('time_max')}"),
]
yy = 0.985
for key, val in table:
    if val is None:
        yy -= 0.010
        ax.text(0.0, yy, key.upper(), transform=ax.transAxes, color=C.ORANGE,
                fontsize=8.2, fontweight='semibold', va='top')
        yy -= 0.036
        continue
    ax.text(0.02, yy, key, transform=ax.transAxes, color=C.INK2, fontsize=8.8,
            va='top')
    ax.text(0.40, yy, val, transform=ax.transAxes, color=C.INK, fontsize=8.8,
            va='top')
    yy -= 0.0345

fig.suptitle('Stage B · the configuration behind every number in this deck',
             x=0.045, ha='left', fontsize=15, fontweight='semibold', color=C.INK, y=0.965)
fig.text(0.045, 0.895,
         'Accepted B.2 coupled gate — job 20492853, binary built from the promoted '
         'shared source with both new switches off by default.',
         ha='left', fontsize=9.5, color=C.INK2)

C.save(fig, 'fig02_configuration',
       'Sources: parties.inp and Data_0.h5 / Particle_0.h5 of '
       'B2closure_coupled16_rampfix_20492853. Reproducible inputs are archived in '
       'ECCO_TESTS/StageB2_closure/B2_closure/coupled16_rampfix/.')
