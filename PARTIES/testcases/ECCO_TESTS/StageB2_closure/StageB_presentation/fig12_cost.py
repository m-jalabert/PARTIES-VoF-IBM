#!/usr/bin/env python3
"""Figure 12 — what Stage B cost, and what B.3 needs.

(a) allocation phy250167 as recorded in the roadmap
(b) where the B.2 closure campaign's SU went, from Slurm accounting
    (B2_closure/campaign_jobs.json, refreshed by audit_campaign.py)
(c) the forward requirement, against what is left and the two request levels
(d) the cost basis, and the single highest-leverage thing to settle first

SU is charged per allocated core-hour; campaign_jobs.json records cpu_hours per
job, including the runs that timed out, failed or were cancelled.
"""
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stageb_common as C

plt = C.use_style()
N = C.ROADMAP_NUMBERS
ALLOC = N['allocation']

jobs = C.campaign_jobs()['jobs']

BUCKETS = OrderedDict([
    ('ice-free settling controls\nand resolution ladder', ('B2settle',)),
    ('coupled gate + timestep triple\n(with and without the ramp fix)',
     ('B2coupled16', 'B2coupled_v9')),
    ('MPI, restart, remap, budgets,\nStage-A regression',
     ('B2remap', 'B2restart', 'B2stagea', 'B2seed', 'B2profile', 'B2hot_water',
      'B2melt_closed')),
    ('rejected hypotheses\n(refill, v6, superseded coupled)',
     ('B2tiny', 'B2coupled_v6', 'B2coupled_2')),
])


def bucket_of(name):
    for label, prefixes in BUCKETS.items():
        if any(name.startswith(p) for p in prefixes):
            return label
    return 'rejected hypotheses\n(refill, v6, superseded coupled)'


su = OrderedDict((k, 0.0) for k in BUCKETS)
timeout = OrderedDict((k, 0.0) for k in BUCKETS)
for j in jobs:
    b = bucket_of(j['JobName'])
    su[b] += j['cpu_hours']
    if j['State'].split()[0] != 'COMPLETED':
        timeout[b] += j['cpu_hours']
total = sum(su.values())

fig = plt.figure(figsize=(12.6, 7.2))
gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.26,
                      left=0.068, right=0.975, top=0.845, bottom=0.10)

# --- (a) allocation ------------------------------------------------------
# The roadmap's 309,843-SU snapshot is 2026-08-31 and already contains the
# 3,425 SU Stage B had spent by then; the closure campaign ran afterwards, so
# it comes out of what the snapshot called "remaining".
ax = C.finish(fig.add_subplot(gs[0, 0]))
stage_b_before = 3425
before_stage_b = ALLOC['used'] - stage_b_before
remaining_now = ALLOC['remaining'] - total
stage_b_total = stage_b_before + total
spent_now = ALLOC['used'] + total
segments = [('Stage A / Yang, and everything before', before_stage_b, C.AXIS),
            ('Stage B to 2026-08-31', stage_b_before, C.BLUE),
            ('B.2 closure campaign, September', total, C.ORANGE),
            ('remaining', remaining_now, C.STATUS['good'])]
left = 0.0
handles = []
for label, val, col in segments:
    h = ax.barh([0], [val], left=left, height=0.55, color=col, zorder=3,
                label=f'{label}   {val:,.0f} SU')
    handles.append(h)
    left += val
ax.set_xlim(0, ALLOC['limit'])
ax.set_ylim(-2.6, 0.6)
ax.set_yticks([])
ax.set_xlabel('SU of the phy250167 allocation')
ax.set_title(f'a   Stage B is {100*stage_b_total/spent_now:.1f} % of what the '
             'project has spent')
ax.grid(axis='y', visible=False)
ax.legend(loc='lower left', ncol=1, handlelength=1.1, handleheight=0.9,
          borderaxespad=0.2)

# --- (b) where the closure SU went --------------------------------------
ax = C.finish(fig.add_subplot(gs[0, 1]))
labels = list(su)
vals = np.array([su[k] for k in labels])
tos = np.array([timeout[k] for k in labels])
ypos = np.arange(len(labels))[::-1]
ax.barh(ypos, vals, height=0.5, color=C.BLUE, zorder=3, label='completed')
ax.barh(ypos, tos, height=0.5, color='#9ec5f4', zorder=4,
        label='timed out, failed or cancelled')
for yv, v, t in zip(ypos, vals, tos):
    ax.text(v + 9, yv, f'{v:,.0f} SU', va='center', fontsize=9, color=C.INK2)
ax.set_yticks(ypos)
ax.set_yticklabels(labels, fontsize=8.6)
ax.set_xlim(0, max(vals) * 1.32)
ax.set_xlabel('SU (core-hours), from Slurm accounting')
ax.set_title(f'b   The B.2 closure campaign: {total:,.0f} SU in {len(jobs)} jobs')
ax.legend(loc='lower right')
ax.grid(axis='y', visible=False)

# --- (c) forward requirement --------------------------------------------
ax = C.finish(fig.add_subplot(gs[1, 0]))
rows = N['forward']
ypos = np.arange(len(rows))[::-1]
for yv, (label, lo, hi) in zip(ypos, rows):
    ax.plot([lo, hi], [yv, yv], color=C.AXIS, lw=3, solid_capstyle='round', zorder=2)
    ax.plot([lo], [yv], 'o', ms=8, color=C.STATUS['good'], mec=C.SURFACE, mew=1.2,
            zorder=4, label='tuned throughput' if yv == ypos[0] else None)
    ax.plot([hi], [yv], 'o', ms=8, color=C.STATUS['warning'], mec=C.SURFACE, mew=1.2,
            zorder=4, label='measured throughput' if yv == ypos[0] else None)
ax.axvline(ALLOC['remaining'], color=C.STATUS['critical'], lw=1.6)
ax.text(ALLOC['remaining'] * 0.90, ypos[0] + 0.80,
        f"remaining  {ALLOC['remaining'] - total:,.0f} SU",
        color=C.STATUS['critical'], fontsize=8.2, rotation=90, va='top', ha='right')
for key, col, ls in [('minimum', C.BLUE, '--'), ('recommended', C.VIOLET, ':')]:
    ax.axvline(N['request'][key], color=col, lw=1.4, ls=ls)
    ax.text(N['request'][key] * 1.06, ypos[0] + 0.80,
            f"{key} request  {N['request'][key]:,} SU", color=col, fontsize=8.2,
            rotation=90, va='top', ha='left')
ax.set_xscale('log')
ax.set_yticks(ypos)
ax.set_yticklabels([r[0] for r in rows], fontsize=8.8)
ax.set_xlim(5e3, 1.2e6)
ax.set_ylim(ypos[-1] - 0.9, ypos[0] + 0.9)
ax.set_xlabel('SU required')
ax.set_title('c   What B.3 costs, against what is left')
ax.legend(loc='lower left')
ax.grid(axis='y', visible=False)

# --- (d) the cost basis --------------------------------------------------
ax = fig.add_subplot(gs[1, 1])
ax.axis('off')
ax.text(0.0, 1.06, 'd   The cost basis, and what to settle first',
        transform=ax.transAxes, fontsize=11, fontweight='semibold',
        color=C.INK, va='bottom')
items = [
    ('Measured throughput', '2 763 cell-steps per core-second',
     'job 20218407: 655 k cells × 5 600 steps on 64 ranks in 5.76 h'),
    ('It is a pessimistic floor', '≈ 3× recoverable',
     'the shakeout runs a 22³ per-rank subdomain and is communication-bound'),
    ('The timestep is the leverage', 'a factor of 6',
     'CFL 0.3 at Δx = d/24 allows dt = 1.25e-2; the deck asks for 2e-3, and the\n'
     '3-D shakeout ran stably at 1e-2. Confirm the true limit before requesting time.'),
    ('The least certain line item', 'the Sc = 70 arms',
     'the salinity Batchelor scale is 8.4× thinner on the same grid — the cell count\n'
     'does not change, but the scalar-solver iteration count does'),
    ('Not affordable as specified', 'the 36-run full factorial',
     'replace it in the write-up by the one-at-a-time design costed in panel c'),
]
yy = 0.90
for head, value, body in items:
    ax.text(0.0, yy, '›', transform=ax.transAxes, color=C.ORANGE, fontsize=12,
            fontweight='bold', va='center')
    ax.text(0.035, yy, head, transform=ax.transAxes, color=C.INK, fontsize=9.4,
            va='center')
    ax.text(0.99, yy, value, transform=ax.transAxes, color=C.ORANGE, fontsize=9.4,
            va='center', ha='right', fontweight='semibold')
    ax.text(0.035, yy - 0.055, body, transform=ax.transAxes, color=C.INK2,
            fontsize=8.3, va='top', linespacing=1.4)
    yy -= 0.175 + 0.05 * body.count('\n')

fig.suptitle('Stage B · what it cost, and what B.3 will need',
             x=0.068, ha='left', fontsize=15, fontweight='semibold', color=C.INK, y=0.965)
fig.text(0.068, 0.895,
         'Allocation figures are the roadmap’s 2026-08-31 snapshot; the campaign '
         'breakdown is live Slurm accounting and includes every rejected trial.',
         ha='left', fontsize=9.5, color=C.INK2)

C.save(fig, 'fig12_cost',
       'Sources: B2_closure/campaign_jobs.json (b, refreshed by audit_campaign.py); '
       'roadmap "RESOURCE STATUS AND FORWARD REQUIREMENT" for the allocation, the '
       'forward table and the request levels. Production lines multiply by 6 if '
       'max_dt = 2e-3 turns out to be required.')
