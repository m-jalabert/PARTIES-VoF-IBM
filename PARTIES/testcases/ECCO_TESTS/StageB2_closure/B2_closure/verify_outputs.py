#!/usr/bin/env python3
"""Independent checks of saved physical fields; fail on a violated tolerance."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from analyze_profiles import deck


def compare(a, b, groups):
    differences = {}
    with h5py.File(a) as fa, h5py.File(b) as fb:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset) and any(name == g or name.startswith(g + '/') for g in groups):
                x, y = obj[...], fb[name][...]
                assert x.shape == y.shape, name
                differences[name] = float(np.max(np.abs(x - y))) if x.size else 0.
        fa.visititems(visit)
    assert differences, 'No fields compared'
    return differences


def profiles(run):
    d = deck(run / 'parties.inp')
    rows = np.genfromtxt(run / 'ecco_profiles.csv', delimiter=',', names=True)
    errors = {}
    count = 0
    for path in sorted(run.glob('Data_*.h5')):
        with h5py.File(path) as f:
            time = f['time'][0]
            b = rows[rows['time'] == time]
            ny = int(d['NYM'][0])
            assert len(b) == ny, (path, len(b), ny)
            def field(name):
                return f[name][...][:-1, :-1, :-1]
            x, y, z = [f['grid/' + c][...].ravel()[:-1] for c in ['xc', 'yc', 'zc']]
            dx, dy, dz = x[1]-x[0], y[1]-y[0], z[1]-z[0]
            dv = dx*dy*dz
            cl, cs = field('VOF/C_L'), field('VOF/C_S')
            weight = np.where(cs < .05, cl, 0.) * dv
            den = weight.sum(axis=(0, 2))
            ice = np.where(cs < .05, np.maximum(0., 1.-cs-cl), 0.)
            T, S, tr = [field('Conc/' + str(c)) for c in range(3)]
            u, v, w = [f[c][...] for c in ['u', 'v', 'w']]
            uc = .5*(u[:-1, :-1, :-1]+u[:-1, :-1, 1:])
            vc = .5*(v[:-1, :-1, :-1]+v[:-1, 1:, :-1])
            wc = .5*(w[:-1, :-1, :-1]+w[1:, :-1, :-1])
            expected = dict(liquid_volume=den, ice_volume=ice.sum(axis=(0, 2))*dv)
            def mean(a):
                return np.divide((weight*a).sum(axis=(0, 2)), den,
                                 out=np.full_like(den, np.nan), where=den > 0)
            expected['mean_v'] = mean(vc)
            for name, scalar, integral in [('T', T, 'T'), ('S', S, 'S'), ('tracer', tr, 'tracer')]:
                expected['mean_'+name] = mean(scalar)
                expected['flux_v'+name] = mean(vc*scalar)-mean(vc)*mean(scalar)
                expected[integral+'_integral'] = scalar.sum(axis=(0, 2))*dv
            # Independent centered differences: periodic x,z and interior y.
            ux = (u[:-1, :-1, 1:]-u[:-1, :-1, :-1])/dx
            vy = (v[:-1, 1:, :-1]-v[:-1, :-1, :-1])/dy
            wz = (w[1:, :-1, :-1]-w[:-1, :-1, :-1])/dz
            def dd(a, axis, spacing):
                return (np.roll(a, -1, axis)-np.roll(a, 1, axis))/(2*spacing)
            eps = 2/d['Re'][0]*(ux*ux+vy*vy+wz*wz + .5*(
                (dd(uc, 1, dy)+dd(vc, 2, dx))**2 +
                (dd(uc, 0, dz)+dd(wc, 2, dx))**2 +
                (dd(vc, 0, dz)+dd(wc, 1, dy))**2))
            expected['epsilon_integral'] = (weight*eps).sum(axis=(0, 2))
            for name, a in expected.items():
                observed = b[name]
                if name == 'epsilon_integral':
                    a, observed = a[1:-1], observed[1:-1]
                assert np.array_equal(np.isnan(a), np.isnan(observed)), name
                err = float(np.nanmax(np.abs(a-observed)))
                scale = max(1., float(np.nanmax(np.abs(a))))
                errors[name] = max(errors.get(name, 0.), err/scale)
                assert err/scale < 1.e-11, (path, name, err, scale)
        count += 1
    assert count
    return dict(snapshots=count, max_scaled_errors=errors,
                epsilon_scope='all x,z; interior y rows; wall stencil not checked')


def remap(a, b):
    diff = compare(a/'Data_0.h5', b/'Data_0.h5', ['VOF'])
    assert max(diff.values()) == 0., diff
    d = deck(a/'parties.inp')
    with h5py.File(a/'Data_0.h5') as f, h5py.File(a/'Particle_0.h5') as p:
        X = p['mobile/X'][...].reshape(-1, 3)[0]
        R = p['mobile/R'][...].ravel()[0]
        x, y, z = [f['grid/'+c][...].ravel()[:-1] for c in ['xc', 'yc', 'zc']]
        delta = [x[None, None, :]-X[0], y[None, :, None]-X[1], z[:, None, None]-X[2]]
        for axis, key in [(0, 'x'), (2, 'z')]:
            length = d[key+'max'][0]-d[key+'min'][0]
            delta[axis] -= length*np.rint(delta[axis]/length)
        radius = np.sqrt(sum(c*c for c in delta))
        cn = d['Cn'][0]
        expected = .5-.5*np.tanh((radius-(R-np.sqrt(2)*np.log(19)*cn))/(2*np.sqrt(2)*cn))
        actual = f['VOF/C_S'][...][:-1, :-1, :-1]
        # The solver deliberately truncates negligible far tails.
        err = float(np.max(np.abs(expected-actual)))
        assert err < 1.e-10, err
    mass = []
    import re
    for run in [a, b]:
        m = re.search(r'REMAP_FIXTURE before=(\S+) after=(\S+) delta=(\S+)', (run/'run.log').read_text())
        assert m, run
        before, after, change = map(float, m.groups())
        assert abs(change) < 1.e-12*abs(before), change
        mass.append(dict(before=before, after=after, change=change))
    return dict(field_max_absolute_differences=diff, analytic_periodic_CS_max_error=err, mass=mass)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=['profiles', 'remap', 'restart'])
    p.add_argument('a', type=Path)
    p.add_argument('b', type=Path, nargs='?')
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    if args.mode == 'profiles':
        result = profiles(args.a)
    elif args.mode == 'remap':
        result = remap(args.a, args.b)
    else:
        result = {}
        for number in [3, 4]:
            for prefix, groups in [('Data', ['u', 'v', 'w', 'p', 'VOF', 'Conc']), ('Particle', ['mobile'])]:
                file = f'{prefix}_{number}.h5'
                diff = compare(args.a/file, args.b/file, groups)
                assert max(diff.values()) == 0., (file, diff)
                result[file] = diff
    text = json.dumps(result, indent=2)+'\n'
    print(text)
    if args.output:
        args.output.write_text(text)
