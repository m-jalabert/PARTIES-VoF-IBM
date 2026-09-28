#!/usr/bin/env python3
"""Audit recorded fields, excluding duplicated high boundary planes.
No settling correlation or pass verdict is inferred from repeatability.
"""
import argparse, json
from pathlib import Path
import h5py
import numpy as np

def audit(path):
    out=[]
    for file in sorted(path.glob('Data_*.h5'),key=lambda p:int(p.stem.split('_')[-1])):
        with h5py.File(file,'r') as f:
            t=float(f['time'][0]); cl=f['VOF/C_L'][...][:-1,:-1,:-1];cs=f['VOF/C_S'][...][:-1,:-1,:-1]
            x=f['grid/xc'][...].ravel()[:-1];y=f['grid/yc'][...].ravel()[:-1];z=f['grid/zc'][...].ravel()[:-1]
            dx=x[1]-x[0];dy=y[1]-y[0];dz=z[1]-z[0];dv=float(dx*dy*dz)
            ice=np.clip(1-cl-cs,0,1)
            vel=[f[c][...][:-1,:-1,:-1] for c in ['u','v','w']]
            speed=np.sqrt(sum(a*a for a in vel))
            real_ice=np.where(cs<0.05,ice,0.0)
            row=dict(real_ice_volume=float(real_ice.sum()*dv),real_ice_max=float(real_ice.max()),t=t,liquid_volume=float(cl.sum()*dv),solid_volume=float(cs.sum()*dv),ice_volume=float(ice.sum()*dv),speed_max=float(speed.max()))
            for c in range(3):
                a=f[f'Conc/{c}'][...][:-1,:-1,:-1];row[f'c{c}_volume']=float(a.sum()*dv);row[f'c{c}_min']=float(a.min());row[f'c{c}_max']=float(a.max())
            bulk=(ice>.95)&(cs<.001)
            row['ice_speed_max']=float(speed[bulk].max()) if bulk.any() else None
            pf=file.with_name(file.name.replace('Data_','Particle_'))
            if pf.exists():
                with h5py.File(pf) as p:
                    X=p['mobile/X'][...].reshape(-1,3)[0];U=p['mobile/U'][...].reshape(-1,3)[0];R=float(p['mobile/R'][...].ravel()[0])
                row.update(y=float(X[1]),uy=float(U[1]))
                xx=x[None,None,:]-X[0];yy=y[None,:,None]-X[1];zz=z[:,None,None]-X[2]
                r=np.sqrt(xx**2+yy**2+zz**2)
                water_region=(y[None,:,None]<9.0)&(cs<.001)
                wake=water_region&(yy>0)&(r>R+3*dx)
                below=water_region&(yy<0)&(r>R+3*dx)
                for name,mask in [('wake',wake),('below',below),('water',water_region)]:
                    row[name+'_ice_mean']=float(ice[mask].mean()) if mask.any() else None
                    row[name+'_ice_max']=float(ice[mask].max()) if mask.any() else None
                    row[name+'_ice_gt09']=int(np.count_nonzero((ice>.9)&mask))
            out.append(row)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--output',type=Path);a=p.parse_args()
    rows=audit(a.directory)
    s=json.dumps(rows,indent=2)
    if a.output:a.output.write_text(s+'\n')
    else:print(s)
