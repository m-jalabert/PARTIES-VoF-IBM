#!/usr/bin/env python3
"""Compare release-relative falling trajectories before wall impact.

Acceptance for this debugging comparison: <=2% of fall distance in y,
<=5% of peak falling speed in Uy, <=2% in release-to-contact duration,
and <=0.5% in release time. These are debugging tolerances, not a DNS
production timestep-accuracy specification. Contact uses centre y < 0.51 d.
"""
import argparse
import json
from pathlib import Path
import numpy as np


def read(run):
    a = np.loadtxt(run/'mobile.dat', delimiter=',')
    r = np.genfromtxt(run/'release.dat', delimiter=',', names=True)
    released = r[r['t_released'] >= 0]
    if not len(released) or not np.any(a[:, 3] < .51):
        raise ValueError('Both release and wall contact must be observed')
    release = float(released['t_released'][0])
    contact = float(a[a[:, 3] < .51, 0][0])
    return a, release, contact


def compare(a, b):
    ma, ra, ca = read(a)
    mb, rb, cb = read(b)
    end = min(ca-ra, cb-rb)-.5
    assert end > .1
    t = np.linspace(.1, end, 3001)
    ya = np.interp(t, ma[:, 0]-ra, ma[:, 3])
    yb = np.interp(t, mb[:, 0]-rb, mb[:, 3])
    va = np.interp(t, ma[:, 0]-ra, ma[:, 6])
    vb = np.interp(t, mb[:, 0]-rb, mb[:, 6])
    fall = max(abs(ma[0, 3]-.5), abs(mb[0, 3]-.5))
    result = dict(release_a=ra, release_b=rb, contact_a=ca, contact_b=cb,
        release_relative_difference=abs(ra-rb)/max(abs(ra), abs(rb)),
        fall_duration_relative_difference=abs((ca-ra)-(cb-rb))/max(ca-ra, cb-rb),
        max_y_difference_over_fall=float(np.max(abs(ya-yb))/fall),
        max_Uy_difference_over_peak=float(np.max(abs(va-vb))/max(abs(va).max(), abs(vb).max())),
        comparison_ends_before_contact_by=.5)
    result['debug_timestep_pass'] = bool(
        result['release_relative_difference'] <= .005 and
        result['fall_duration_relative_difference'] <= .02 and
        result['max_y_difference_over_fall'] <= .02 and
        result['max_Uy_difference_over_peak'] <= .05)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('a', type=Path)
    p.add_argument('b', type=Path)
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    text = json.dumps(compare(args.a, args.b), indent=2)+'\n'
    print(text)
    if args.output:
        args.output.write_text(text)
