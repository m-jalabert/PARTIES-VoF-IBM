#!/usr/bin/env python3
"""Audit available P1 output, including an unfinished feasibility segment.

Read complete HDF5 frames only. Budget summaries use the archived B2 audit;
partial profile writes should be retried after the next completed block.
Never infer completed precursor acceptance from an unfinished/no-release run.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import re
import subprocess

import h5py
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument('run', type=Path)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    run = a.run
    prof = Path(__file__).parent.parent/'B2_closure/analyze_profiles.py'
    if not prof.exists():
        prof = run/'analyze_profiles.py'
    spec = importlib.util.spec_from_file_location('profiles', prof)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    summary, rows = module.analyze(run)
    out = dict(run=str(run), budgets=summary,
               net_melt_volume=rows[0]['ice_volume']-rows[-1]['ice_volume'],
               initial_ice_volume=rows[0]['ice_volume'],
               final_ice_volume=rows[-1]['ice_volume'])
    r = np.atleast_1d(np.genfromtxt(run/'release.dat', delimiter=',', names=True))
    m = np.atleast_2d(np.loadtxt(run/'mobile.dat', delimiter=','))
    out['particle'] = dict(released=bool(np.any(r['t_released']>=0)),
                           max_displacement=float(np.max(np.linalg.norm(m[:,2:5]-m[0,2:5],axis=1))),
                           max_speed=float(np.max(np.linalg.norm(m[:,5:8],axis=1))),
                           initial_y=float(m[0,3]),final_y=float(m[-1,3]),
                           initial_logged_shell_fraction=float(r['phi_liq'][0]) if 'phi_liq' in r.dtype.names else None,
                           first_computed_shell_fraction=float(r['phi_liq'][r['time']>0][0]) if 'phi_liq' in r.dtype.names and np.any(r['time']>0) else None,
                           shell_fraction_note='The t=0 logged value is an initialization placeholder; first positive-time sample is the first computed shell fraction.',
                           release_log_columns=list(r.dtype.names))
    paths = sorted(run.glob('Data_*.h5'), key=lambda x:int(x.stem.split('_')[-1]))
    frames = []
    for path in dict.fromkeys([paths[0], paths[-1]]):
        try:
            with h5py.File(path) as f:
                cs=f['VOF/C_S'][...][:-1,:-1,:-1]
                cl=f['VOF/C_L'][...][:-1,:-1,:-1]
                ice=np.where(cs<.05,np.clip(1-cs-cl,0,1),0)
                speed2=np.zeros_like(cl)
                for key in ['u','v','w']:
                    v=f[key][...][:-1,:-1,:-1]
                    speed2+=v*v
                bulk=(ice>.95)&(cs<.001)
                physical=(cs<.05)
                liquid=physical&(cl>.95)
                vmax=float(np.sqrt(speed2.max()))
                imax=float(np.sqrt(speed2[bulk].max())) if bulk.any() else None
                frame=dict(file=path.name,time=float(f['time'][0]),speed_max=vmax,
                           bulk_ice_speed_max=imax,bulk_ice_speed_ratio=imax/vmax if vmax>1e-14 and imax is not None else None)
                for i in range(3):
                    v=f[f'Conc/{i}'][...][:-1,:-1,:-1]
                    frame[f'c{i}_physical_range']=[float(v[physical].min()),float(v[physical].max())]
                    frame[f'c{i}_liquid_range']=[float(v[liquid].min()),float(v[liquid].max())]
                y=f['grid/yc'][...].ravel()[:-1]
                below=(y[None,:,None]<7)&(cs<.001)
                frame['ice_max_below_y7']=float(ice[below].max())
                frames.append(frame)
        except (OSError,KeyError) as exc:
            frames.append(dict(file=path.name,incomplete=str(exc)))
    out['frames']=frames
    log=(run/'run.log').read_text(errors='replace')
    iterations=re.findall(r'Iteration (\d+), dt = ([\deE.+-]+), time = ([\deE.+-]+)',log)
    out['last_iteration']=iterations[-1] if iterations else None
    out['source_diagnostics']=[line for line in log.splitlines() if any(x in line.lower() for x in ['unplaced','residual ice','nan','failed','error'])][-30:]
    jobid=run.name.split('_')[-1]
    out['slurm']=subprocess.check_output(['sacct','-j',jobid,'--format=JobID,State,ElapsedRaw,AllocCPUS,ExitCode','-P'],text=True)
    out['interpretation']='Feasibility segment only; P1 post-escape acceptance and P2 comparison remain unevaluated.'
    a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2,allow_nan=False))


if __name__ == '__main__':
    main()
