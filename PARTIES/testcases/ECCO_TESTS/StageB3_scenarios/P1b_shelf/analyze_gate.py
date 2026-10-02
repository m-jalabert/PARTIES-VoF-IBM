#!/usr/bin/env python3
"""Score the P1b 1-D kernel gate against the exact binary-Stefan solution.

Criteria (declared 2026-09-28, before any gate result existed):
  C1  |lambda_fit/lambda_exact - 1| <= 0.05, lambda from X(t) = a + b sqrt(t),
      t in [20, 80], X = ice lost per unit area (ecco_profiles.csv)
  C2  net melt never decreases (no net freezing at any profile time)
  C3  at the last field frame, at the F=0.5 point: |theta - theta_i| <= 0.05
      and |s/F - s_i| <= 0.05
  C4  enthalpy |dH|/(X/St) <= 1e-3; salt relative drift <= 1e-8;
      tracer-minus-melt <= 1e-3 of the melt
  C5  (consistency) n32 lambda within 3% and half-dt lambda within 1% of n24
Decision: P1b uses liquid_referenced_salinity=1 only if G1_sliq1_n24 passes
C1-C4; if only G0 passes use 0; if neither passes, P1b is NOT submitted.

Round 2 (added 2026-09-29 after round 1 failed, before any round-2 result):
G2_yang_* use yang_salt_transport=1, where the salt field IS the liquid
salinity and the conserved salt is the integral of (F+delta)*S (from the
profiles: sum of mean_S*liquid_volume).  Same C1-C5.  P1b uses the Yang
operator if G2_yang_n24 passes C1-C4 and G2_yang_n32 passes C5.
"""
import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np


def deck(path):
    d = {}
    for line in path.read_text().splitlines():
        line = line.split('#')[0].strip()
        if '=' in line:
            k, v = line.split('=', 1)
            d[k.strip()] = v.strip()
    return d


def profiles(run):
    rows = list(csv.DictReader(open(run / 'ecco_profiles.csv')))
    out = {}
    for r in rows:
        t = float(r['time'])
        acc = out.setdefault(t, dict(ice=0., T=0., S=0., SF=0., C=0., n=0))
        acc['ice'] += float(r['ice_volume'])
        acc['T'] += float(r['T_integral'])
        acc['S'] += float(r['S_integral'])
        lv = float(r['liquid_volume'])
        acc['SF'] += float(r['mean_S']) * lv if lv > 0 else 0.0
        acc['C'] += float(r['tracer_integral'])
        acc['n'] += 1
    ny = max(v['n'] for v in out.values())
    times = sorted(t for t, v in out.items() if v['n'] == ny)
    return times, [out[t] for t in times]


def interface_state(run, yang):
    files = sorted(run.glob('Data_*.h5'), key=lambda p: int(p.stem.split('_')[1]))
    with h5py.File(files[-1]) as f:
        cl = f['VOF/C_L'][...][:-1, :-1, :-1].mean(axis=(0, 2))
        th = f['Conc/0'][...][:-1, :-1, :-1].mean(axis=(0, 2))
        s = f['Conc/1'][...][:-1, :-1, :-1].mean(axis=(0, 2))
        y = f['grid/yc'][...].ravel()[:-1]
        t = float(f['time'][0])
    j = int(np.where((cl[:-1] >= 0.5) & (cl[1:] < 0.5))[0][0])
    w = (cl[j] - 0.5) / (cl[j] - cl[j + 1])
    lerp = lambda a: float(a[j] + w * (a[j + 1] - a[j]))
    sl = s if yang else s / np.maximum(cl, 1e-12)
    return dict(time=t, y=lerp(y), theta=lerp(th), s_liquid=lerp(s) / (1.0 if yang else 0.5),
                s_liquid_ice_side_max=float(np.max(np.where(cl > 1e-3, sl, 0))))


def score(run, ref):
    d = deck(run / 'parties.inp')
    st, pe_t = float(d['stefan']), ref['Pe_T']
    area = float(d['xmax']) * float(d['zmax'])
    times, acc = profiles(run)
    t = np.array(times)
    X = np.array([(acc[0]['ice'] - a['ice']) / area for a in acc])
    sel = (t >= 20) & (t <= 80.0 + 1e-9)
    b, a0 = np.polyfit(np.sqrt(t[sel]), X[sel], 1)
    lam = b * np.sqrt(pe_t) / 2
    H = np.array([a['T'] - a['ice'] / st for a in acc])
    melt = (acc[0]['ice'] - acc[-1]['ice'])
    yang = int(float(d.get('yang_salt_transport', '0'))) == 1
    key = 'SF' if yang else 'S'
    salt = abs(acc[-1][key] - acc[0][key]) / abs(acc[0][key])
    tracer = abs((acc[-1]['C'] - acc[0]['C']) - melt) / max(melt, 1e-30)
    dH = float(np.max(np.abs(H - H[0])))
    iface = interface_state(run, yang)
    res = dict(lambda_fit=float(lam), lambda_exact=ref['lambda_'],
               lambda_rel_error=float(lam / ref['lambda_'] - 1), fit_offset=float(a0),
               X_final=float(X[-1]), X_exact_final=float(2 * ref['lambda_'] * np.sqrt(t[-1] / pe_t)),
               min_dX=float(np.min(np.diff(X))), enthalpy_over_latent=dH / max(X[-1] * area / st, 1e-30),
               salt_rel_drift=float(salt), tracer_minus_melt_rel=float(tracer), interface=iface,
               theta_err=iface['theta'] - ref['theta_interface'],
               s_err=iface['s_liquid'] - ref['s_interface'], final_time=float(t[-1]))
    res['C1'] = abs(res['lambda_rel_error']) <= 0.05
    res['C2'] = res['min_dX'] >= -1e-9
    res['C3'] = abs(res['theta_err']) <= 0.05 and abs(res['s_err']) <= 0.05
    res['C4'] = res['enthalpy_over_latent'] <= 1e-3 and salt <= 1e-8 and tracer <= 1e-3
    res['pass_C1_C4'] = all(res[k] for k in ['C1', 'C2', 'C3', 'C4'])
    res['series'] = dict(t=t[::10].tolist(), X=X[::10].tolist())
    return res


def main():
    p = argparse.ArgumentParser()
    p.add_argument('gate_dir', type=Path)
    a = p.parse_args()
    ref = json.loads((a.gate_dir / 'reference.json').read_text())
    out = dict(reference=ref)
    for run in sorted(a.gate_dir.iterdir()):
        if (run / 'ecco_profiles.csv').exists():
            try:
                out[run.name] = score(run, ref)
            except Exception as exc:              # report, never hide
                out[run.name] = dict(error=repr(exc))
    g1, g0, g2 = (out.get(k, {}) for k in ('G1_sliq1_n24', 'G0_sliq0_n24', 'G2_yang_n24'))
    for base, pairs in [('G1_sliq1_n24', [('G1_sliq1_n32', 0.03), ('G1_sliq1_n24_halfdt', 0.01)]),
                        ('G2_yang_n24', [('G2_yang_n32', 0.03)])]:
        ref24 = out.get(base, {})
        for name, tol in pairs:
            if 'lambda_fit' in ref24 and 'lambda_fit' in out.get(name, {}):
                out[name]['C5_vs_n24'] = abs(out[name]['lambda_fit'] / ref24['lambda_fit'] - 1) <= tol
    yang_ok = g2.get('pass_C1_C4') and out.get('G2_yang_n32', {}).get('C5_vs_n24', False)
    out['decision'] = ('liquid_referenced_salinity=1' if g1.get('pass_C1_C4') else
                       'liquid_referenced_salinity=0' if g0.get('pass_C1_C4') else
                       'yang_salt_transport=1' if yang_ok else
                       'STOP: no band variant passes; do not submit P1b')
    (a.gate_dir / 'gate_result.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({k: ({kk: v[kk] for kk in v if kk not in ('series',)} if isinstance(v, dict) else v)
                      for k, v in out.items()}, indent=1))


if __name__ == '__main__':
    main()
