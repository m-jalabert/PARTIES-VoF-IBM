#!/usr/bin/env python3
"""Confirm the release-ramp fix changes nothing before the grain is released.

The fix touches only the `release_ramp > 0` branch, which cannot execute while
t_released < 0, so every field up to the release time should be unchanged.  That
is an argument about the source, not about the binary: the two executables were
compiled from different text, so this checks the claim on the saved fields.

Reports the maximum absolute field difference per snapshot.  Snapshots strictly
before the recorded release time are expected to be BIT-IDENTICAL; a nonzero but
roundoff-scale difference there would mean the compiler reordered arithmetic
outside the changed branch, which is still acceptable but must not be called
bit-identical.
"""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np

FIELDS = ['u', 'v', 'w', 'p', 'VOF/C_L', 'VOF/C_S', 'VOF/F',
          'Conc/0', 'Conc/1', 'Conc/2']


def release_time(run):
    r = np.genfromtxt(Path(run)/'release.dat', delimiter=',', names=True)
    got = r[r['t_released'] >= 0]
    return float(got['t_released'][0]) if len(got) else float('inf')


def snapshots(run):
    return sorted(Path(run).glob('Data_*.h5'),
                  key=lambda p: int(p.stem.split('_')[1]))


def compare(a, b):
    rel = min(release_time(a), release_time(b))
    out = {'release_time': rel, 'snapshots': []}
    for fa, fb in zip(snapshots(a), snapshots(b)):
        with h5py.File(fa, 'r') as ha, h5py.File(fb, 'r') as hb:
            ta, tb = float(ha['time'][0]), float(hb['time'][0])
            # Exclude the duplicated high boundary planes, as the other audits do.
            s = (slice(0, 64), slice(0, 96), slice(0, 64))
            worst, where = 0.0, None
            for name in FIELDS:
                if name not in ha or name not in hb:
                    continue
                d = float(np.max(np.abs(ha[name][:][s] - hb[name][:][s])))
                if d > worst:
                    worst, where = d, name
            out['snapshots'].append(dict(
                file=fa.name, time_a=ta, time_b=tb,
                before_release=bool(ta < rel),
                max_abs_difference=worst, worst_field=where,
                bit_identical=bool(worst == 0.0)))
    pre = [s for s in out['snapshots'] if s['before_release']]
    out['pre_release_snapshots'] = len(pre)
    out['pre_release_all_bit_identical'] = bool(pre) and all(s['bit_identical'] for s in pre)
    out['pre_release_max_difference'] = max([s['max_abs_difference'] for s in pre], default=None)
    return out


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
