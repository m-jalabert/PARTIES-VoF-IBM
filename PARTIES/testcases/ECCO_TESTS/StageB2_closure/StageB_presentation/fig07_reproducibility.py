#!/usr/bin/env python3
"""Figure 7 — reproducibility, and the non-invasiveness contract.

Panels
  (a) every reproducibility check in the campaign on one difference ladder
  (b) B.0 part 2: enabling VOF_IBM + LAG_PARTICLE_RESOLVED with no particle
      present perturbs the Stage-A physics only at roundoff, and that roundoff
      grows at a Lyapunov rate rather than as a bias
  (c) rank-count agreement of the coupled d/Δx = 8 pair, before and after the
      release-ramp fix — recomputed here from mobile.dat
  (d) the non-invasiveness contract: how each ECCO change is kept unreachable
      from a pre-ECCO run

The point of the panel-(a) ladder is that "bit-identical" and "agrees to
roundoff" are different claims, and the campaign makes the stronger one only
where it is true.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stageb_common as C

plt = C.use_style()
N = C.ROADMAP_NUMBERS

FLOOR = 1e-17          # where "exactly 0" is drawn on the log ladder

fig = plt.figure(figsize=(12.8, 7.6))
gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.26,
                      left=0.115, right=0.982, top=0.845, bottom=0.105)

# --- (a) the difference ladder ------------------------------------------
ax = C.finish(fig.add_subplot(gs[0, 0]))
ax.set_xscale('log')
rows = N['reproducibility']
ypos = np.arange(len(rows))[::-1]
for yv, (label, val, note) in zip(ypos, rows):
    exact = (val == 0.0)
    x = FLOOR if exact else val
    col = C.STATUS['good'] if exact else C.BLUE
    ax.plot([FLOOR, x], [yv, yv], color=C.GRID, lw=1.2, zorder=1)
    ax.plot([x], [yv], 'o', ms=8, color=col, mec=C.SURFACE, mew=1.2, zorder=4)
    ax.text(x * 2.2, yv, ('bit-identical' if exact else f'{val:.2e}') + f'   ·  {note}',
            va='center', fontsize=8.2, color=C.INK2)
ax.axvspan(FLOOR / 2, FLOOR * 1.7, color='#e8f5e8', lw=0, zorder=0)
ax.set_yticks(ypos)
ax.set_yticklabels([r[0] for r in rows], fontsize=8.6)
ax.set_xlim(FLOOR / 2, 1e-1)
ax.set_ylim(ypos[-1] - 0.8, ypos[0] + 0.8)
ax.set_xticks([1e-17, 1e-14, 1e-11, 1e-8, 1e-5, 1e-2])
ax.set_xticklabels(['exactly 0', '$10^{-14}$', '$10^{-11}$', '$10^{-8}$',
                    '$10^{-5}$', '$10^{-2}$'])
ax.set_xlabel('largest difference against the reference')
ax.set_title('a   Where the campaign claims bit-identity, and where it does not')
ax.grid(axis='y', visible=False)

# --- (b) B.0 roundoff growth --------------------------------------------
ax = C.finish(fig.add_subplot(gs[0, 1]))
ax.set_yscale('log')
ni, lin = N['noninvasive'], N['lineage_bar']
t_ni = np.array(sorted(ni))
ax.axhspan(N['discriminating_size'], 1, color='#fbe7e7', lw=0, zorder=0)
ax.axhline(N['discriminating_size'], color=C.STATUS['critical'], lw=1.3, ls='--')
ax.text(10, N['discriminating_size'] * 1.5,
        f"discriminating size  {N['discriminating_size']:.1e}", ha='right',
        color=C.STATUS['critical'], fontsize=8.5)
ax.plot(t_ni, [ni[k] for k in t_ni], 'o-', color=C.BLUE, ms=7, mec=C.SURFACE,
        mew=1.2, label='Stage-B build vs Stage-A reference')
t_l = np.array(sorted(lin))
ax.plot(t_l, [lin[k] for k in t_l], 's--', color=C.MUTED, ms=6,
        label='already-accepted Stage-A binary lineage test')
for k in t_ni:
    ax.annotate(f'{ni[k]:.1e}', (k, ni[k]), textcoords='offset points',
                xytext=(0, 10), ha='center', fontsize=8.2, color=C.BLUE)
ax.annotate('six orders tighter than the\naccepted bar at the same time',
            xy=(6, lin[6]), xytext=(4.0, 1.4e-7), fontsize=8.5, color=C.INK2,
            ha='center', linespacing=1.5,
            arrowprops=dict(arrowstyle='-', color=C.MUTED, lw=1.0))
ax.set_xlim(1, 11)
ax.set_ylim(1e-13, 3e-2)
ax.set_xlabel('time  $t$  of the freshwater 288² comparison')
ax.set_ylabel('relative difference in ice volume')
ax.set_title('b   Turning the ECCO flags on perturbs Stage A only at roundoff')
ax.legend(loc='lower right')

# --- (c) rank-count agreement -------------------------------------------
ax = C.finish(fig.add_subplot(gs[1, 0]))
ax.set_yscale('log')


def rank_pair(a, b):
    ta, Xa, _ = C.mobile(C.run(a))
    tb, Xb, _ = C.mobile(C.run(b))
    t = np.linspace(0, min(ta[-1], tb[-1]), 4000)
    return t, np.abs(np.interp(t, ta, Xa[:, 1]) - np.interp(t, tb, Xb[:, 1]))


t_new, d_new = rank_pair('v9_rampfix', 'v9_rampfix_mpi8')
old = C.audit('coupled_comparisons.json')['coupled_v9_mpi8']
rel8 = C.release_time(C.run('v9_rampfix'))
CONTACT = 44.0
ax.plot(t_new, np.maximum(d_new, 1e-18), color=C.STATUS['good'], lw=1.8,
        label='after the release-ramp fix  (20492915 / 20492916)')
ko = np.array(sorted(float(k) for k in old))
ax.plot(ko, [old[str(int(k))]['max_y_diff'] for k in ko], 's--', color=C.MUTED,
        ms=6, label='before it  (20477868 / 20478207), recorded')
ax.axvline(rel8, color=C.STATUS['critical'], lw=1.2, ls='--')
ax.text(rel8 - 0.7, 2e-9, 'release', color=C.STATUS['critical'], fontsize=8,
        ha='right', va='bottom')
ax.axvline(CONTACT, color=C.MUTED, lw=1.1, ls=':')
ax.text(CONTACT + 0.8, 1.5e-9, 'wall contact', color=C.MUTED, fontsize=8,
        rotation=90, va='bottom')
m = t_new <= CONTACT
ax.annotate(f'{d_new[m].max():.2e} before contact,\n'
            f"{old['44']['max_y_diff']/d_new[m].max():.0f}× better than before the fix",
            xy=(CONTACT, d_new[m].max()), xytext=(19, 1.0e-7),
            fontsize=8.5, color=C.STATUS['good'],
            arrowprops=dict(arrowstyle='-', color=C.STATUS['good'], lw=1.0))
ax.set_xlim(0, 62)
ax.set_ylim(1e-9, 3e-2)
ax.set_xlabel('time  $t$')
ax.set_ylabel(r'$\max\,|\Delta y|$  between 8 and 16 ranks')
ax.set_title('c   Same physics on a different rank count')
ax.legend(loc='upper left')

# --- (d) the contract ----------------------------------------------------
ax = fig.add_subplot(gs[1, 1])
ax.axis('off')
ax.text(0.0, 1.05, 'd   Non-invasiveness: every ECCO change is flag-gated',
        transform=ax.transAxes, fontsize=11, fontweight='semibold',
        color=C.INK, va='bottom')
contract = [
    ('release machinery, motion lock, ramp, release.dat',
     'compiled out', 'LAG_PARTICLE_RESOLVED && VOF_IBM'),
    ('sediment ice transport, moving-mask remap',
     'compiled out', 'VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT'),
    ('in-solver profile output',
     'compiled out', 'ECCO_PROFILES'),
    ('meltwater tracer source',
     'runtime-gated off', 'meltwater_tracer = 0 and NConc ≥ 3'),
    ('ECCO scalar inits 34 / 35',
     'never selected', 'a pre-ECCO deck does not name them'),
    ('release threshold F_release',
     'off by default', 'default −1.0; the calibrated 0.709 lives in the decks'),
]
yy = 0.86
for change, how, gate in contract:
    ax.text(0.0, yy, '✓', transform=ax.transAxes, color=C.STATUS['good'],
            fontsize=11, fontweight='bold', va='center')
    ax.text(0.05, yy, change, transform=ax.transAxes, color=C.INK,
            fontsize=9.2, va='center')
    ax.text(0.99, yy, how, transform=ax.transAxes, color=C.STATUS['good'],
            fontsize=8.8, va='center', ha='right')
    ax.text(0.05, yy - 0.055, gate, transform=ax.transAxes, color=C.MUTED,
            fontsize=8.2, va='top')
    yy -= 0.145
ax.text(0.0, -0.02,
        'Stated honestly: the builds are not byte-reproducible, so the contract is\n'
        'enforced functionally. The Stage-A gate has been re-run on every round and\n'
        'reproduces λ = 0.2184495545645199 and enthalpy drift 4.975e-10 exactly.',
        transform=ax.transAxes, color=C.INK2, fontsize=8.3, va='top', linespacing=1.5)

fig.suptitle('Stage B · does it reproduce, and does it disturb anything that came before',
             x=0.02, ha='left', fontsize=15, fontweight='semibold', color=C.INK, y=0.965)
fig.text(0.02, 0.895,
         'Reproducibility is claimed per quantity, not globally: phase fields and '
         'pre-release trajectories are bit-identical; the post-contact convective flow is not, '
         'and is not claimed to be.',
         ha='left', fontsize=9.5, color=C.INK2)

C.save(fig, 'fig07_reproducibility',
       'Sources: recorded values from the closure report and roadmap B.0 / B.1.0 '
       '(a, b, d); mobile.dat of B2closure_coupled_v9_rampfix_20492915 / _mpi8_20492916 '
       'and B2_closure/coupled_comparisons.json (c).')
