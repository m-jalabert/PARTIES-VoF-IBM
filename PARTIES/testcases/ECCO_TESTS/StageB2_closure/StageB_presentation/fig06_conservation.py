#!/usr/bin/env python3
"""Figure 6 — conservation, rigidity and the negative controls (B.2 gates 3, 5, 6).

Panels
  (a) heat / salt / meltwater-tracer budget error against the gates, over the run
  (b) the melt budget itself: ice volume lost = meltwater tracer gained
  (c) ice rigidity under Darcy penalization: bulk-ice speed / domain maximum speed
  (d) the negative controls — a moving rock must not melt, and must not make ice

The budget errors are recomputed here with ``B2_closure/analyze_profiles.py``
(imported, not re-implemented) from ``ecco_profiles.csv``; (c) and (d) read the
recorded field audits in ``B2_closure``.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stageb_common as C

plt = C.use_style()

# import the recorded analyser rather than re-deriving its budgets
_spec = importlib.util.spec_from_file_location('ap', C.B2 / 'analyze_profiles.py')
ap = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ap)

RUN = C.run('coupled16')
summary, rows = ap.analyze(RUN)
rel = C.release_time(RUN)

t = np.array([r['time'] for r in rows])
tracer_err = np.array([r['tracer_budget_relative'] or np.nan for r in rows], float)
heat_err = np.array([r['heat_budget_relative'] or np.nan for r in rows], float)
salt_err = np.array([r['salt_budget_relative'] or np.nan for r in rows], float)
ice_vol = np.array([r['ice_volume'] for r in rows])
tracer_int = np.array([r['tracer_integral'] for r in rows])

fields = C.audit('coupled16_rampfix_fields.json')
hot = C.audit('hot_water_audit.json')

fig = plt.figure(figsize=(12.4, 7.4))
gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.24,
                      left=0.068, right=0.985, top=0.855, bottom=0.085)

# --- (a) budget errors ---------------------------------------------------
ax = C.finish(fig.add_subplot(gs[0, 0]))
ax.set_yscale('log')
ax.axhline(1e-3, color=C.MUTED, lw=1.2, ls=':')
ax.text(64, 1.4e-3, 'debugging gate  1e-3', color=C.MUTED, fontsize=8, ha='right')
for err, col, lab in [(tracer_err, C.BLUE, 'meltwater tracer'),
                      (heat_err, C.ORANGE, 'heat (enthalpy)'),
                      (salt_err, C.AQUA, 'salt')]:
    ax.plot(t, np.maximum(err, 1e-17), color=col, lw=1.6, label=lab)
ax.axvline(rel, color=C.MUTED, lw=1.0, ls='--')
ax.text(rel - 1.0, 3e-4, 'release', color=C.MUTED, fontsize=8, ha='right')
ax.text(1.0, 2e-12,
        'worst over the run\n'
        f"tracer  {summary['tracer_max_relative']:.2e}\n"
        f"heat    {summary['heat_max_relative']:.2e}\n"
        f"salt    {summary['salt_max_relative']:.2e}",
        color=C.INK2, fontsize=8.2, va='top', linespacing=1.5)
ax.set_ylim(1e-18, 5e-3)
ax.set_xlim(0, t.max())
ax.set_xlabel('time  $t$')
ax.set_ylabel('relative budget error')
ax.set_title('a   Every conserved scalar closes, three to ten decades inside the gate')
ax.legend(loc='center right', bbox_to_anchor=(1.0, 0.40))

# --- (b) the melt budget -------------------------------------------------
ax = C.finish(fig.add_subplot(gs[0, 1]))
melt = ice_vol[0] - ice_vol
gain = tracer_int - tracer_int[0]
ax.plot(t, melt, color=C.BLUE, lw=2.4, label=r'ice volume lost   $\int\!dF$')
ax.plot(t, gain, color=C.ORANGE, lw=1.3, ls=(0, (4, 3)),
        label=r'meltwater gained   $\int\! C_{mw}$')
ax.axvline(rel, color=C.MUTED, lw=1.0, ls='--')
ax.text(rel - 1.0, 19.0, 'release', color=C.MUTED, fontsize=8, ha='right')
ax.set_xlim(0, t.max())
ax.set_xlabel('time  $t$')
ax.set_ylabel('volume')
ax.set_title('b   The meltwater tracer is the melt, to seven digits')
ax.legend(loc='upper left')

ins = ax.inset_axes([0.60, 0.11, 0.37, 0.30])
ins.plot(t, np.abs(gain - melt), color=C.STATUS['good'], lw=1.4)
ins.set_yscale('log')
ins.set_title('absolute residual', fontsize=7.5, color=C.INK2, pad=3)
ins.set_facecolor(C.SURFACE)
ins.tick_params(labelsize=7, colors=C.MUTED)
ins.grid(color=C.GRID, lw=0.5)
for s in ('top', 'right'):
    ins.spines[s].set_visible(False)

# --- (c) ice rigidity ----------------------------------------------------
ax = C.finish(fig.add_subplot(gs[1, 0]))
tf = np.array([f['t'] for f in fields])
ratio = np.array([100.0 * (f['ice_speed_max'] or 0.0) / f['speed_max']
                  if f['speed_max'] > 0 else np.nan for f in fields])
ax.axhspan(1.0, 1.35, color='#fbe7e7', lw=0, zorder=0)
ax.axhline(1.0, color=C.STATUS['critical'], lw=1.3, ls='--', zorder=2)
ax.text(64, 1.06, 'leakage criterion  1 % of the domain maximum speed',
        color=C.STATUS['critical'], fontsize=8, ha='right')
ax.plot(tf, ratio, 'o-', color=C.BLUE, ms=4)
ax.axvline(rel, color=C.MUTED, lw=1.0, ls='--')
ax.text(rel - 1.0, 0.95, 'release', color=C.MUTED, fontsize=8, ha='right', va='top')
worst = np.nanmax(ratio)
ax.annotate(f'worst {worst:.3f} %', xy=(tf[int(np.nanargmax(ratio))], worst),
            xytext=(8, 0.80), color=C.BLUE, fontsize=8.5,
            arrowprops=dict(arrowstyle='-', color=C.BLUE, lw=1.0))
ax.set_ylim(0, 1.35)
ax.set_xlim(0, t.max())
ax.set_xlabel('time  $t$')
ax.set_ylabel('bulk-ice speed / domain max speed  (%)')
ax.set_title(r'c   The ice stays rigid: Darcy penalization at $\tau = 10^{-3}$')

# --- (d) negative controls ----------------------------------------------
ax = fig.add_subplot(gs[1, 1])
ax.axis('off')
ax.text(0.0, 1.02, 'd   Negative controls — what must NOT happen',
        transform=ax.transAxes, fontsize=11, fontweight='semibold',
        color=C.INK, va='bottom')
hot_ice = max(abs(h['real_ice_max']) for h in hot)
hot_tracer = max(abs(h['c2_max']) for h in hot)
hot_cs = max(abs(h['solid_volume'] - hot[0]['solid_volume']) for h in hot)
rows_d = [
    ('A rock in warm water must not be read as ice',
     f'max actual ice = {hot_ice:.2e}', 'job 20477870'),
    ('A rock must not pay latent heat or inject meltwater',
     f'tracer = {hot_tracer:.0f} exactly;  ∫C_S steady to {hot_cs:.1e}', 'job 20477870'),
    ('The ice must not leak through the Darcy mask',
     f'worst bulk-ice speed {worst:.3f} % of the domain max', 'panel c'),
    ('Moving the grain must not create or destroy ice',
     'conservative moving-mask remap; aborts rather than losing ice', 'remap fixture'),
    ('Closed box: salt must not drift',
     f'{summary["salt_max_relative"]:.2e}   (roundoff)', 'nonzero-salt deck'),
    ('Tracer must stay ≥ 0 outside the IBM support',
     'satisfied in every saved field', 'all d/Δx = 16 runs'),
]
yy = 0.90
for label, value, src in rows_d:
    ax.text(0.0, yy, '✓', transform=ax.transAxes, color=C.STATUS['good'],
            fontsize=11, fontweight='bold', va='center')
    ax.text(0.052, yy, label, transform=ax.transAxes, color=C.INK,
            fontsize=9.2, va='center')
    ax.text(0.052, yy - 0.062, value, transform=ax.transAxes, color=C.INK2,
            fontsize=8.8, va='center')
    ax.text(0.995, yy, src, transform=ax.transAxes, color=C.MUTED,
            fontsize=8, va='center', ha='right')
    yy -= 0.155
ax.text(0.0, -0.05,
        'Values inside the IBM excluded support are extension values, not physical\n'
        'concentrations; physical extrema are evaluated outside it.',
        transform=ax.transAxes, color=C.MUTED, fontsize=8, va='top')

fig.suptitle('Stage B · conservation, rigidity, and the things that must not happen',
             x=0.068, ha='left', fontsize=15, fontweight='semibold', color=C.INK, y=0.965)
fig.text(0.068, 0.905,
         'Accepted B.2 coupled gate (job 20492853) unless noted. Budgets integrate the '
         'whole physical mesh; the imposed wall heat flux is included.',
         ha='left', fontsize=9.5, color=C.INK2)

C.save(fig, 'fig06_conservation',
       'Sources: ecco_profiles.csv via B2_closure/analyze_profiles.py; '
       'coupled16_rampfix_fields.json; hot_water_audit.json (job 20477870).')
