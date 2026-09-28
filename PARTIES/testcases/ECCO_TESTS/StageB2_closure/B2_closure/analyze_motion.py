#!/usr/bin/env python3
"""Particle-motion gate quantities for the B.2 coupled release/settle arms.

Reconstructs the B.2 gate table (checks 1, 2 and 4) from `mobile.dat` and
`release.dat`.  The original summary was generated inline and not preserved, so
this script re-derives it; `--check` validates it against a recorded
motion_summary.json entry before it is trusted on new runs.

a_free = |g| (1 - rho_f/rho_s) is the buoyancy-corrected free acceleration; the
gate compares the largest acceleration in the time unit after release with it.
"""
import argparse
import json
from pathlib import Path
import numpy as np


def deck_value(run, key, default=None):
    for line in (Path(run)/'parties.inp').read_text().splitlines():
        line = line.split('#')[0].strip()
        if line.startswith(key) and '=' in line:
            return line.split('=', 1)[1].strip()
    return default


def summarize(run, wall_y=0.51, window=1.0):
    run = Path(run)
    a = np.loadtxt(run/'mobile.dat', delimiter=',')
    t, X, U = a[:, 0], a[:, 2:5], a[:, 5:8]

    rho_s = float(deck_value(run, 'rho_s', '2.5'))
    grav = deck_value(run, 'grav', '{0.0, -1.0, 0.0}').strip('{}')
    g = abs(float(grav.split(',')[1]))
    a_free = g*(1.0 - 1.0/rho_s)

    r = np.genfromtxt(run/'release.dat', delimiter=',', names=True)
    got = r[r['t_released'] >= 0]
    out = {'a_free': a_free}
    if not len(got):
        return dict(out, release_time=None, released=False)
    rel = float(got['t_released'][0])

    locked = t < rel
    out['release_time'] = rel
    out['released'] = True
    out['locked_displacement_max'] = float(np.max(np.linalg.norm(X[locked]-X[0], axis=1))) if locked.any() else 0.
    out['locked_speed_max'] = float(np.max(np.linalg.norm(U[locked], axis=1))) if locked.any() else 0.

    acc = np.diff(U, axis=0)/np.diff(t)[:, None]
    idx = np.where((t >= rel) & (t <= rel+window))[0]
    out['release_acceleration_over_free'] = float(
        np.max(np.linalg.norm(acc[idx[0]:idx[-1]], axis=1))/a_free)

    fall = t >= rel
    out['min_y'] = float(X[fall, 1].min())
    out['final_y'] = float(X[-1, 1])
    out['final_U'] = [float(v) for v in U[-1]]
    out['min_Uy'] = float(U[fall, 1].min())
    hit = np.where(X[:, 1] < wall_y)[0]
    out['first_wall_time'] = float(t[hit[0]]) if len(hit) else None
    return out


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('runs', type=Path, nargs='+')
    p.add_argument('--name', nargs='*', help='label per run (default: directory name)')
    p.add_argument('--check', type=Path, help='recorded motion_summary.json to validate against')
    p.add_argument('--output', type=Path)
    args = p.parse_args()

    names = args.name if args.name else [r.name for r in args.runs]
    result = {n: summarize(r) for n, r in zip(names, args.runs)}

    if args.check:
        recorded = json.loads(args.check.read_text())
        worst = 0.0
        for n, got in result.items():
            if n not in recorded:
                continue
            for k, v in recorded[n].items():
                if k not in got or v is None or got[k] is None:
                    continue
                d = np.max(np.abs(np.asarray(v, float) - np.asarray(got[k], float)))
                worst = max(worst, float(d))
        result['_validation_max_difference_vs_recorded'] = worst
        result['_validation_reproduces_recorded'] = bool(worst < 1e-12)

    text = json.dumps(result, indent=2)+'\n'
    print(text)
    if args.output:
        args.output.write_text(text)
