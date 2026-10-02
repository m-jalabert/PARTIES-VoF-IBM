#!/usr/bin/env python3
"""Check the P1b smoke run: initial profiles against the prescribed formulas,
geometry, first-step melt direction, grain acceleration, and cost per step."""
import argparse
import json
import re
from pathlib import Path

import h5py
import numpy as np
from scipy.special import erf


def main():
    p = argparse.ArgumentParser()
    p.add_argument('run', type=Path)
    p.add_argument('--case', type=Path, required=True)
    a = p.parse_args()
    case = json.loads(a.case.read_text())
    ini = case['initial_state']
    y_int, cn = case['geometry']['y_interface'], case['numerics']['Cn']
    out = {}
    with h5py.File(a.run / 'Data_0.h5') as f:
        y = f['grid/yc'][...].ravel()[:-1]
        cl = f['VOF/C_L'][...][:-1, :-1, :-1]
        cs = f['VOF/C_S'][...][:-1, :-1, :-1]
        k, i = 0, 0                                  # far from the grain (periodic corner column)
        col = lambda name: f[name][...][:-1, :-1, :-1][k, :, i]
        th, s, c, F = col('Conc/0'), col('Conc/1'), col('Conc/2'), cl[k, :, i]
    zeta = np.maximum(y_int - y, 0.0)
    fl = 0.5 * (1 - np.tanh((y - y_int) / (2 * np.sqrt(2) * cn)))
    s_liq = ini['s_interface'] + (1 - ini['s_interface']) * erf(zeta / ini['ell_S_d'])
    yang = int(case['numerics'].get('yang_salt_transport', 0)) == 1
    s_expected = s_liq if yang else fl * s_liq     # Yang: salt field is the liquid salinity
    out['yang_salt_transport'] = yang
    out['ic_max_error'] = dict(theta=float(np.abs(th - erf(zeta / ini['ell_T_d'])).max()),
                               salinity=float(np.abs(s - s_expected).max()),
                               tracer=float(np.abs(c - fl * (1 - s_liq)).max()),
                               liquid_fraction_vs_tanh=float(np.abs(F - fl).max()))
    out['tracer_min'] = float(c.min())
    ice = np.where(cs < 0.05, np.clip(1 - cs - cl, 0, 1), 0)
    out['ice_volume_fraction_top2d'] = float(ice[:, y > y_int + 0.3, :].mean())
    solid = cs >= 0.05
    yy = np.broadcast_to(y[None, :, None], cs.shape)
    out['grain_support_y_range'] = [float(yy[solid].min()), float(yy[solid].max())]
    out['grain_top_liquid_fraction_min'] = float(cl[solid | (cs > 0.001)].min()) if solid.any() else None
    log = (a.run / 'run.log').read_text(errors='replace')
    it = re.findall(r'Iteration (\d+), dt = ([\deE.+-]+), time = ([\deE.+-]+)', log)
    out['iterations'] = len(it)
    out['last'] = it[-1] if it else None
    m = np.atleast_2d(np.loadtxt(a.run / 'mobile.dat', delimiter=','))
    out['grain'] = dict(t=float(m[-1, 0]), y0=float(m[0, 3]), y=float(m[-1, 3]),
                        v=float(m[-1, 6]), v_over_ws=float(-m[-1, 6] / case['scaling']['w_s_over_Uref']))
    prof = np.genfromtxt(a.run / 'ecco_profiles.csv', delimiter=',', names=True)
    t0, t1 = prof['time'].min(), prof['time'].max()
    out['ice_volume_change'] = float(prof['ice_volume'][prof['time'] == t1].sum()
                                     - prof['ice_volume'][prof['time'] == t0].sum())
    out['nan_in_log'] = 'nan' in log.lower()
    (a.run / 'smoke_check.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
