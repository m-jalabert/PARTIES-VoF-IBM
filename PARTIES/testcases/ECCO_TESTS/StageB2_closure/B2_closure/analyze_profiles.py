#!/usr/bin/env python3
"""Audit physical-cell budgets and integrate measured dissipation/heat flux.
The energy check supports the tested uniform grid with kappa_ice_ratio[0]=1.
No mixing efficiency is inferred from BPE changes in an open melting box.
"""
import argparse,csv,json
from pathlib import Path
import numpy as np

def deck(path):
    d={}
    for line in path.read_text().splitlines():
        line=line.split('#')[0].split('//')[0].strip()
        if '=' not in line:continue
        key,val=line.split('=',1);val=val.strip().strip('{}')
        try:d[key.strip()]=[float(x) for x in val.split(',')]
        except ValueError:pass
    return d

def analyze(run):
    d=deck(run/'parties.inp');a=np.genfromtxt(run/'ecco_profiles.csv',delimiter=',',names=True)
    ny=int(d['NYM'][0]);dy=(d['ymax'][0]-d['ymin'][0])/ny
    area=(d['xmax'][0]-d['xmin'][0])*(d['zmax'][0]-d['zmin'][0]);Pe=d['Pe'][0];St=d['stefan'][0]
    if abs(d['kappa_ice_ratio'][0]-1)>1e-14:raise ValueError('Heat audit requires unit ice/water thermal diffusivity ratio')
    rows=[]
    for time in np.unique(a['time']):
        b=a[a['time']==time][-ny:]
        if len(b)!=ny or len(np.unique(b['y']))!=ny:raise ValueError('Incomplete profile block')
        sums={k:float(b[k].sum()) for k in ['ice_volume','T_integral','S_integral','tracer_integral','epsilon_integral']}
        heat=0.0
        for side,row in [('S',b[0]),('N',b[-1])]:
            A=d['BC_A'+side][0];B=d['BC_B'+side][0];C=d['BC_C'+side][0]
            meanT=row['T_integral']/(area*dy)
            heat+=area/Pe*(C-B*meanT)/(A+B*dy/2)
        rows.append(dict(time=float(time),**sums,heat_in_rate=heat))
    t=np.array([r['time'] for r in rows]);n=len(rows)
    for field,target in [('epsilon_integral','dissipation_integrated'),('heat_in_rate','heat_in_integrated')]:
        v=np.array([r[field] for r in rows]);c=np.r_[0.,np.cumsum(.5*(v[1:]+v[:-1])*np.diff(t))]
        for r,val in zip(rows,c):r[target]=float(val)
    first=rows[0];max_tr=max_h=max_s=0.0
    for r in rows:
        melt=first['ice_volume']-r['ice_volume'];tr=r['tracer_integral']-first['tracer_integral']
        r['tracer_budget_absolute']=tr-melt
        r['tracer_budget_relative']=abs(tr-melt)/abs(melt) if abs(melt)>1e-8 else None
        if r['tracer_budget_relative'] is not None:max_tr=max(max_tr,r['tracer_budget_relative'])
        r['salt_budget_absolute']=r['S_integral']-first['S_integral']
        r['salt_budget_relative']=abs(r['salt_budget_absolute'])/abs(first['S_integral']) if abs(first['S_integral'])>1e-10 else None
        if r['salt_budget_relative'] is not None:max_s=max(max_s,r['salt_budget_relative'])
        if St>0:
            H=r['T_integral']-r['ice_volume']/St;H0=first['T_integral']-first['ice_volume']/St
            r['heat_budget_absolute']=H-H0-r['heat_in_integrated']
            r['heat_budget_relative']=abs(r['heat_budget_absolute'])/max(abs(H0),abs(r['heat_in_integrated']),1e-20)
            max_h=max(max_h,r['heat_budget_relative'])
    result={'start_time':float(t[0]),'end_time':float(t[-1]),'profile_times':n,'max_time_gap':float(np.diff(t).max()) if n>1 else 0.,
      'tracer_max_relative':max_tr,'tracer_max_absolute':max(abs(r['tracer_budget_absolute']) for r in rows),
      'heat_max_relative':max_h if St>0 else None,'salt_nonzero_tested':abs(first['S_integral'])>1e-10,
      'salt_max_relative':max_s if abs(first['S_integral'])>1e-10 else None,
      'salt_max_absolute':max(abs(r['salt_budget_absolute']) for r in rows),
      'dissipation_integrated':rows[-1]['dissipation_integrated'],
      'minimum_dissipation':min(r['epsilon_integral'] for r in rows)}
    return result,rows

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('--output',type=Path);args=p.parse_args()
    result,rows=analyze(args.run)
    print(json.dumps(result,indent=2))
    if args.output:
        args.output.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
        with args.output.with_suffix('.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[-1]));w.writeheader();w.writerows(rows)
