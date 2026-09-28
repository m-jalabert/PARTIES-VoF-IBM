#!/usr/bin/env python3
"""Timestep refinement for the coupled release/settle diagnostic.

`compare_timesteps.py` scores one pair against fixed debugging tolerances on
release-relative trajectories.  That single score mixes two errors with very
different character, so this script separates them and, given three timesteps,
measures the observed order of each.

  escape_lag        The grain is released while still Darcy-penalized inside the
                    ice.  Break-out is a threshold event, so its timing carries
                    an O(dt) error.  Everything after it inherits a pure time
                    offset, measured here as the difference in release-relative
                    arrival time at a reference height BELOW the ice.

  open_water_Uy     After break-out the trajectories are parallel.  Comparing
                    Uy as a function of HEIGHT removes the inherited offset and
                    leaves the temporal error of the settling dynamics itself.

Neither is a production accuracy specification; both are convergence diagnostics.
"""
import argparse
import json
from pathlib import Path
import numpy as np

# Reference heights, in grain diameters (d = 1).  ESCAPE_Y is below the melting
# ice and below the height at which the Uy(y) curves rejoin; WATER_Y0/1 bound the
# open-water fall, staying clear of the wall-collision layer.
ESCAPE_Y = 3.30
WATER_Y0, WATER_Y1 = 0.80, 3.30


def read(run):
    run = Path(run)
    a = np.loadtxt(run/'mobile.dat', delimiter=',')
    r = np.genfromtxt(run/'release.dat', delimiter=',', names=True)
    released = r[r['t_released'] >= 0]
    if not len(released):
        raise ValueError(f'{run}: no release recorded')
    release = float(released['t_released'][0])
    fall = a[a[:, 0] >= release]
    if fall[:, 3].min() > WATER_Y0:
        raise ValueError(f'{run}: fall does not reach y={WATER_Y0}')
    # After release the centre height decreases monotonically; reverse it so it
    # is an increasing coordinate for np.interp.
    y = fall[:, 3][::-1]
    return dict(release=release, y=y, Uy=fall[:, 6][::-1],
                tau=(fall[:, 0]-release)[::-1],
                contact=float(fall[fall[:, 3] < .51, 0][0]))


def pair(a, b):
    tau_a = float(np.interp(ESCAPE_Y, a['y'], a['tau']))
    tau_b = float(np.interp(ESCAPE_Y, b['y'], b['tau']))
    y = np.linspace(WATER_Y0, WATER_Y1, 2001)
    va = np.interp(y, a['y'], a['Uy'])
    vb = np.interp(y, b['y'], b['Uy'])
    peak = max(abs(va).max(), abs(vb).max())
    return dict(
        release_difference=abs(a['release']-b['release']),
        escape_lag=abs(tau_a-tau_b),
        escape_arrival_a=tau_a, escape_arrival_b=tau_b,
        open_water_Uy_max_relative=float(np.max(abs(va-vb))/peak),
        open_water_Uy_rms_relative=float(np.sqrt(((va-vb)**2).mean())/peak),
        contact_difference=abs((a['contact']-a['release'])-(b['contact']-b['release'])))


def order(coarse_fine, fine_finest, ratio=2.0):
    """Observed order p from two successive differences at refinement `ratio`."""
    if fine_finest <= 0 or coarse_fine <= 0:
        return None
    return float(np.log(coarse_fine/fine_finest)/np.log(ratio))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('runs', type=Path, nargs='+',
                   help='run directories, coarsest timestep first')
    p.add_argument('--dt', type=float, nargs='+', required=True,
                   help='max_dt of each run, same order')
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    if len(args.runs) != len(args.dt):
        p.error('one --dt per run')

    data = [read(r) for r in args.runs]
    result = dict(runs=[str(r) for r in args.runs], max_dt=list(args.dt), pairs=[])
    for i in range(len(data)-1):
        d = pair(data[i], data[i+1])
        d['max_dt'] = [args.dt[i], args.dt[i+1]]
        result['pairs'].append(d)

    if len(result['pairs']) >= 2:
        c, f = result['pairs'][-2], result['pairs'][-1]
        ratio = args.dt[-3]/args.dt[-2]
        result['observed_order'] = {
            k: order(c[k], f[k], ratio)
            for k in ('escape_lag', 'open_water_Uy_max_relative',
                      'open_water_Uy_rms_relative', 'contact_difference')}
        # Richardson estimate of the finest run's own error, using the observed
        # order where it is meaningful (p > 0.2) and first order otherwise.
        result['richardson_error_of_finest'] = {}
        for k in result['observed_order']:
            pk = result['observed_order'][k]
            pk = pk if (pk is not None and pk > 0.2) else 1.0
            result['richardson_error_of_finest'][k] = float(f[k]/(ratio**pk - 1.0))

    text = json.dumps(result, indent=2)+'\n'
    print(text)
    if args.output:
        args.output.write_text(text)
