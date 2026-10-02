#!/usr/bin/env python3
"""Score dev3: Yang repair at resolved sediment (build_v3) and v2/v3 non-invasiveness.

Pass criteria (declared before the result):
  R1  Y_fix completes to t=1.2 with liquid S in [s_i - 0.02, 1 + 1e-6] at every saved frame
      (Y_orig is expected to leave that range, reproducing job 20973340)
  R2  G2r (1-D Yang gate, no grain) lambda equals round-2 G2_yang_n24 (0.019413) to 1e-6
  R3  NI_v2 and NI_v3 Data files bit-identical in every field (default operator untouched)
"""
import argparse
import glob
import json
import re
from pathlib import Path

import h5py
import numpy as np


def frames(run):
    return sorted(glob.glob(str(run / 'Data_*.h5')), key=lambda p: int(re.findall(r'_(\d+)\.h5', p)[0]))


def salt_ranges(run):
    out = []
    for p in frames(run):
        with h5py.File(p) as f:
            g = lambda k: f[k][...][:-1, :-1, :-1]
            cs, cl, s = g('VOF/C_S'), g('VOF/C_L'), g('Conc/1')
            liq = (cs < 0.05) & (cl > 0.95)
            fl = cs < 0.05
            out.append(dict(t=float(f['time'][0]), S_liquid=[float(s[liq].min()), float(s[liq].max())],
                            S_fluid=[float(s[fl].min()), float(s[fl].max())]))
    return out


def identical(a, b):
    diffs = {}
    for pa in frames(a):
        pb = b / Path(pa).name
        with h5py.File(pa) as fa, h5py.File(pb) as fb:
            def walk(name, obj):
                if isinstance(obj, h5py.Dataset) and name in fb:
                    x, y = obj[...], fb[name][...]
                    if x.shape != y.shape or not np.array_equal(x, y):
                        diffs[f'{Path(pa).name}:{name}'] = float(np.max(np.abs(x.astype(float) - y.astype(float)))) if x.shape == y.shape else 'shape'
            fa.visititems(walk)
    return diffs


def main():
    p = argparse.ArgumentParser()
    p.add_argument('dev3', type=Path)
    a = p.parse_args()
    d = a.dev3
    case = json.loads((d / 'case.json').read_text())
    s_i = case['initial_state']['s_interface']
    res = dict(status=(d / 'status.txt').read_text().split('\n'))
    for r in ('Y_orig', 'Y_fix'):
        rng = salt_ranges(d / r)
        log = (d / r / 'run.log').read_text(errors='replace')
        last = re.findall(r'Iteration (\d+), dt = ([\deE.+-]+), time = ([\deE.+-]+)', log)
        m = np.atleast_2d(np.loadtxt(d / r / 'mobile.dat', delimiter=','))
        ok = all(x['S_liquid'][0] >= s_i - 0.02 and x['S_liquid'][1] <= 1 + 1e-6 for x in rng)
        res[r] = dict(frames=rng, last_iteration=last[-1] if last else None,
                      poisson_failure='did not converge' in log, grain_final=dict(t=float(m[-1, 0]), y=float(m[-1, 3]), v=float(m[-1, 6])),
                      bounded=ok)
    res['R1'] = res['Y_fix']['bounded'] and float(res['Y_fix']['last_iteration'][2]) >= 1.19
    import importlib.util
    spec = importlib.util.spec_from_file_location('ag', d / 'analyze_gate.py')
    ag = importlib.util.module_from_spec(spec); spec.loader.exec_module(ag)
    ref = json.loads((d / 'G2r' / '..' / 'case.json').read_text())['gate']['reference']
    ref.update(Pe_T=case['scaling']['Pe_T'])
    g = ag.score(d / 'G2r', ref)
    res['G2r'] = {k: g[k] for k in ('lambda_fit', 'lambda_rel_error', 'salt_rel_drift', 'theta_err', 's_err')}
    res['R2'] = abs(g['lambda_fit'] - 0.019412720023085504) < 1e-6
    diffs = identical(d / 'NI_v2', d / 'NI_v3')
    res['NI_differences'] = diffs
    res['R3'] = len(diffs) == 0 and len(frames(d / 'NI_v2')) > 1
    (d / 'dev3_result.json').write_text(json.dumps(res, indent=2) + '\n')
    print(json.dumps({k: v for k, v in res.items() if k not in ('Y_orig', 'Y_fix')}, indent=1))
    for r in ('Y_orig', 'Y_fix'):
        print(r, res[r]['last_iteration'], 'poisson_failure', res[r]['poisson_failure'], 'bounded', res[r]['bounded'])
        for x in res[r]['frames']:
            print('   t=%.2f S_liquid [%.4f, %.4f]  S_fluid [%.3f, %.3f]' % (x['t'], *x['S_liquid'], *x['S_fluid']))


if __name__ == '__main__':
    main()
