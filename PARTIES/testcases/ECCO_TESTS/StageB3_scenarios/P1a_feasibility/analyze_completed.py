#!/usr/bin/env python3
"""Read every P1 field frame and preserve budgets, accounting and figures.

No solver execution or raw-output modification. All arrays exclude duplicated
high boundary planes. The tracer is signed under freezing; it is not assumed
to be a nonnegative water-origin fraction. Source reconstruction is diagnostic
at saved times, not an RK-integrated source audit.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


def write_csv(path, rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def crossing(y, v, level):
    hits=np.where((v[:-1]>=level)&(v[1:]<level))[0]
    if not len(hits): return None
    i=hits[np.argmin(abs(y[hits]-8.))]
    return float(y[i]+(level-v[i])*(y[i+1]-y[i])/(v[i+1]-v[i]))


def plot_final_slice(out, x, y, last, time):
    fig,axes=plt.subplots(1,3,figsize=(12,5),constrained_layout=True)
    for ax,arr,title in zip(axes,[last[2],34*last[3],last[4]],['Temperature θ','Model salinity 34s (g/kg)','Signed phase-change tracer']):
        data=np.ma.masked_where(last[1]>=.05,arr)
        signed=title.startswith('Signed')
        norm=TwoSlopeNorm(vcenter=0.,vmin=min(float(data.min()),-1e-15),vmax=max(float(data.max()),1e-15)) if signed else None
        im=ax.pcolormesh(x,y,data,shading='nearest',cmap='coolwarm' if signed else 'viridis',norm=norm);fig.colorbar(im,ax=ax)
        ax.contour(x,y,last[0],levels=[.5],colors='black',linewidths=.8);ax.set_ylim(6.8,9.2);ax.set_aspect('equal');ax.set_xlabel('x/d');ax.set_title(title)
    axes[0].set_ylabel('y/d');fig.suptitle(f'P1 final saved field, t={time:.6f}; sediment masked')
    fig.savefig(out/'P1_final_slice.png',dpi=180);fig.savefig(out/'P1_final_slice.pdf');plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path);ap.add_argument('output',type=Path)
    args=ap.parse_args();run=args.run.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    case=json.loads((run/'case.json').read_text())
    spec=importlib.util.spec_from_file_location('budgets',run/'analyze_profiles.py')
    bmod=importlib.util.module_from_spec(spec);spec.loader.exec_module(bmod)
    budget,brows=bmod.analyze(run)
    write_csv(out/'budgets.csv',brows)
    (out/'budgets.json').write_text(json.dumps(budget,indent=2)+'\n')
    d=bmod.deck(run/'parties.inp');dx=1/24;dv=dx**3
    files=sorted(run.glob('Data_*.h5'),key=lambda p:int(p.stem.split('_')[-1]))
    assert [int(p.stem.split('_')[-1]) for p in files]==list(range(len(files)))
    frames=[];profiles=[];slices={}
    for n,path in enumerate(files):
        with h5py.File(path) as f:
            t=float(f['time'][0]);x=f['grid/xc'][...].ravel()[:-1];y=f['grid/yc'][...].ravel()[:-1];z=f['grid/zc'][...].ravel()[:-1]
            cs=f['VOF/C_S'][...][:-1,:-1,:-1];cl=f['VOF/C_L'][...][:-1,:-1,:-1]
            physical=cs<.05
            # The C_S tanh is shifted so C_S=.05 marks the grain radius.
            # Its integral is not grain volume. As in ECCO_write_profiles,
            # exclude C_S>=.05 before integrating actual ice.
            ice=np.where(physical,np.clip(1-cs-cl,0,1),0)
            liquid=physical&(cl>.95);bulk=(ice>.95)&(cs<.001)
            speed2=np.zeros_like(cl)
            for key in ['u','v','w']:
                v=f[key][...][:-1,:-1,:-1];assert np.isfinite(v).all(),(path,key);speed2+=v*v
            assert np.isfinite(cs).all() and np.isfinite(cl).all()
            vmax=float(np.sqrt(speed2.max()));imax=float(np.sqrt(speed2[bulk].max()))
            row=dict(frame=n,time=t,ice_volume=float(ice.sum()*dv),solid_indicator_integral=float(cs.sum()*dv),excluded_sediment_volume=float(np.count_nonzero(~physical)*dv),
                     speed_max=vmax,bulk_ice_speed_max=imax,ice_speed_ratio=imax/vmax if vmax>0 else 0.,
                     bulk_ice_max_y=float(y[np.unravel_index(np.argmax(np.where(bulk,speed2,-1)),cl.shape)[1]]))
            scalar=[]
            for i in range(3):
                v=f[f'Conc/{i}'][...][:-1,:-1,:-1];assert np.isfinite(v).all(),(path,i);scalar.append(v)
                row[f'c{i}_physical_min']=float(v[physical].min());row[f'c{i}_physical_max']=float(v[physical].max())
                row[f'c{i}_liquid_min']=float(v[liquid].min());row[f'c{i}_liquid_max']=float(v[liquid].max())
            theta,s,tr=scalar
            superheat=theta-d['T_melt'][0]+d['liquidus_slope'][0]*s
            pore=1-cs;weight=cl*np.maximum(pore-cl,0)/(pore+1e-30)/(np.sqrt(2)*d['Cn'][0])
            m=d['stefan'][0]/(d['Pe'][0]*d['melt_band_eps'][0])*superheat*weight
            m[~physical]=0
            row.update(source_positive=float(np.maximum(m,0).sum()*dv),
                       source_negative=float(np.minimum(m,0).sum()*dv),
                       tracer_negative_integral=float(np.minimum(tr,0).sum()*dv),
                       tracer_positive_integral=float(np.maximum(tr,0).sum()*dv),
                       negative_tracer_liquid_cells=int(np.count_nonzero(liquid&(tr< -1e-10))))
            band=physical&(cl>.05)&(cl<.95)
            row['band_superheat_min']=float(superheat[band].min());row['band_superheat_max']=float(superheat[band].max())
            for bound in [6,7]:
                mask=np.broadcast_to(y[None,:,None]<bound,cl.shape)&(cs<.001)
                row[f'ice_max_below_y{bound}']=float(ice[mask].max())
            # Far from the stationary grain: horizontal distance > 1.25d.
            far=((z[:,None]-2)**2+(x[None,:]-2)**2)>1.25**2
            layer=np.moveaxis(cl,1,0)[:,far].mean(axis=1)
            yi=crossing(y,layer,.5);y95=crossing(y,layer,.95);y05=crossing(y,layer,.05)
            row.update(far_interface_y=yi,far_interface_width_cells=(y05-y95)/dx)
            if n in [0,40,200,400]:
                far_scalars=[np.moveaxis(v,1,0)[:,far].mean(axis=1) for v in scalar]
                for j,yj in enumerate(y):
                    profiles.append(dict(frame=n,time=t,y=float(yj),liquid_fraction=float(layer[j]),
                                         theta=float(far_scalars[0][j]),salinity=float(far_scalars[1][j]),tracer=float(far_scalars[2][j])))
            if n in [0,len(files)-1]:
                iz=np.argmin(abs(z-2));slices[n]=[a[iz].copy() for a in [cl,cs,theta,s,tr,np.sqrt(speed2)]]
                row['tracer_negative_cell_location']=list(map(float,[x[np.unravel_index(np.argmin(np.where(physical,tr,np.inf)),tr.shape)[2]],y[np.unravel_index(np.argmin(np.where(physical,tr,np.inf)),tr.shape)[1]],z[np.unravel_index(np.argmin(np.where(physical,tr,np.inf)),tr.shape)[0]]]))
            else:
                row['tracer_negative_cell_location']=None
            frames.append(row)
        if n%40==0: print(f'Fields audited: {n+1}/{len(files)}',flush=True)
    write_csv(out/'fields.csv',frames);write_csv(out/'far_profiles.csv',profiles)
    r=np.atleast_1d(np.genfromtxt(run/'release.dat',delimiter=',',names=True));mob=np.atleast_2d(np.loadtxt(run/'mobile.dat',delimiter=','))
    log=(run/'run.log').read_text();it=re.findall(r'Iteration (\d+), dt = ([\deE.+-]+), time = ([\deE.+-]+)',log)
    steps=int(it[-1][0]);job=run.name.split('_')[-1]
    accounting=subprocess.check_output(['sacct','-j',job,'--format=JobID,State,Start,End,ElapsedRaw,AllocCPUS,CPUTimeRAW,AllocTRES,ExitCode','-P'],text=True)
    (out/'accounting.txt').write_text(accounting)
    acct=list(csv.DictReader(accounting.splitlines(),delimiter='|'))[0]
    elapsed=int(acct['ElapsedRaw']);cores=int(acct['AllocCPUS']);su=int(acct['CPUTimeRAW'])/3600
    balance=subprocess.check_output(['mybalance'],text=True);(out/'allocation.txt').write_text(balance)
    inventory=[]
    for p in sorted(run.iterdir()):
        if p.is_file(): inventory.append(dict(name=p.name,bytes=p.stat().st_size,symlink=p.is_symlink(),target=str(p.readlink()) if p.is_symlink() else ''))
    write_csv(out/'raw_inventory.csv',inventory)
    sha=[]
    for p in sorted(run.iterdir()):
        if p.is_file() and not p.is_symlink() and (p.stat().st_size<40_000_000 or p.name in ['Data_0.h5','Data_400.h5','Resume_400.h5']):
            h=hashlib.sha256()
            with p.open('rb') as f:
                for chunk in iter(lambda:f.read(8*1024**2),b''):h.update(chunk)
            sha.append(f'{h.hexdigest()}  {p.name}')
    (out/'raw_selected.sha256').write_text('\n'.join(sha)+'\n')
    with h5py.File(run/'Resume.h5') as f:
        checkpoint=dict(target=str((run/'Resume.h5').readlink()),time=float(f['time'][0]),step=int(f['ntime'][0]),transport_version=int(f['sediment_ice_transport_version'][0]))
    nonzero=[f for f in frames if f['time']>0];late=[f for f in frames if f['time']>=1]
    worst=max(nonzero,key=lambda f:f['ice_speed_ratio']);worstlate=max(late,key=lambda f:f['ice_speed_ratio'])
    times=np.array([v['time'] for v in frames]);mp=np.array([v['source_positive'] for v in frames]);mn=np.array([v['source_negative'] for v in frames])
    unplaced=[float(x) for x in re.findall(r'unplaced=([\deE.+-]+)',log)]
    eosT=case['eos']['CT_range_C'];eosS=case['eos']['SA_range_g_kg'];phys=case['physical'];scale=case['scaling']
    ctmin=min(v['c0_physical_min'] for v in frames)*scale['DeltaT_K']+phys['CT_ice_C'];ctmax=max(v['c0_physical_max'] for v in frames)*scale['DeltaT_K']+phys['CT_ice_C']
    samin=min(v['c1_physical_min'] for v in frames)*34;samax=max(v['c1_physical_max'] for v in frames)*34
    profile_times=np.array([b['time'] for b in brows])
    matches=[brows[int(np.argmin(abs(profile_times-v['time'])))] for v in frames]
    time_error=max(abs(v['time']-b['time']) for v,b in zip(frames,matches))
    ice_error=max(abs(v['ice_volume']-b['ice_volume']) for v,b in zip(frames,matches))
    assert time_error<1e-10 and ice_error<1e-10,(time_error,ice_error)
    result=dict(job_id=job,status=acct['State'],exit_code=acct['ExitCode'],elapsed_seconds=elapsed,allocated_cores=cores,
        cost_SU=su,cost_basis='CPUTimeRAW/3600 on shared (billing CPU weight=1); job row only, no double-counting steps.',
        steps=steps,throughput_cell_steps_per_core_second=2211840*steps/(cores*elapsed),SU_per_nondimensional_time=su/budget['end_time'],
        end_time=budget['end_time'],physical_duration_s=budget['end_time']*scale['t_ref_s'],budget=budget,
        net_ice_gain=brows[-1]['ice_volume']-brows[0]['ice_volume'],net_ice_gain_fraction=(brows[-1]['ice_volume']/brows[0]['ice_volume']-1),
        mean_equivalent_ice_growth_d=(brows[-1]['ice_volume']-brows[0]['ice_volume'])/16,
        released=bool(np.any(r['t_released']>=0)),particle_displacement_max=float(np.max(np.linalg.norm(mob[:,2:5]-mob[0,2:5],axis=1))),particle_speed_max=float(np.max(np.linalg.norm(mob[:,5:8],axis=1))),
        first_computed_shell_fraction=float(r['phi_liq'][r['time']>0][0]),first_computed_shell_time=float(r['time'][r['time']>0][0]),last_shell_fraction=float(r['phi_liq'][-1]),release_threshold=.709,
        rigidity_worst=worst,rigidity_worst_t_ge_1=worstlate,rigidity_exceedances=sum(f['ice_speed_ratio']>.01 for f in nonzero),nonzero_frames=len(nonzero),
        maximum_bulk_ice_speed_m_s=max(v['bulk_ice_speed_max'] for v in frames)*scale['U_ref_m_s'],
        maximum_unplaced_residual=max(map(abs,unplaced)),all_field_arrays_finite=True,
        observed_CT_range_C=[ctmin,ctmax],observed_SA_range_g_kg=[samin,samax],observed_states_within_fitted_box=ctmin>=eosT[0]-1e-12 and ctmax<=eosT[1]+1e-12 and samin>=eosS[0]-1e-12 and samax<=eosS[1]+1e-12,
        negative_tracer_liquid_cells_max=max(v['negative_tracer_liquid_cells'] for v in frames),
        tracer_physical_min=min(v['c2_physical_min'] for v in frames),tracer_physical_max=max(v['c2_physical_max'] for v in frames),
        saved_time_source_positive_integral=float(np.trapezoid(mp,times)),saved_time_source_negative_integral=float(np.trapezoid(mn,times)),
        initial_frame=frames[0],final_frame=frames[-1],checkpoint=checkpoint,raw_file_count=len(inventory),raw_bytes_excluding_symlinks=sum(v['bytes'] for v in inventory if not v['symlink']),
        field_frame_count=len(files),resume_frame_count=len(list(run.glob('Resume_*.h5'))),particle_frame_count=len(list(run.glob('Particle_*.h5'))),
        field_profile_max_time_difference=time_error,field_profile_ice_max_absolute_difference=ice_error,
        next_decision='Do not extend blindly: no net melt-out, signed tracer under freezing, rigidity criterion exceeded. Resolve these model/numerical interpretation gates before long P1 or P2.')
    (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    end=brows[-1];first=brows[0]
    interpretation=dict(heat_error_final_absolute=end['heat_budget_absolute'],
        heat_error_final_over_boundary_heat=abs(end['heat_budget_absolute']/end['heat_in_integrated']),
        heat_error_final_over_temperature_change=abs(end['heat_budget_absolute']/(end['T_integral']-first['T_integral'])),
        cumulative_net_boundary_heat=end['heat_in_integrated'],
        net_ice_gain_monotonic_at_profile_times=bool(np.all(np.diff([b['ice_volume'] for b in brows])>=0)),
        P2_same_short_window_SU_estimate=su*2*5242880/2211840,
        one_diameter_conduction_scenario_SU_at_observed_rate=result['SU_per_nondimensional_time']*case['heat_planning']['one_d_retreat_conduction_time'],
        wall_minutes_with_30pct_margin=elapsed/60*1.3)
    (out/'interpretation_metrics.json').write_text(json.dumps(interpretation,indent=2)+'\n')
    # Reproducible, standalone figures with explicit signed/source conventions.
    bt=np.array([v['time'] for v in brows]);melt=np.array([brows[0]['ice_volume']-v['ice_volume'] for v in brows])
    fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    axes[0,0].plot(bt,melt,label='Net ice lost');axes[0,0].plot(bt,[v['tracer_integral'] for v in brows],'--',label='Signed tracer');axes[0,0].set_ylabel('Volume (d³)');axes[0,0].legend();axes[0,0].set_title('Net freezing throughout the segment')
    axes[0,1].plot(r['time'][1:],r['phi_liq'][1:]);axes[0,1].axhline(.709,color='r',ls='--',label='Release trigger');axes[0,1].set_ylabel('Shell liquid fraction');axes[0,1].legend()
    axes[1,0].plot(times[1:],[100*v['ice_speed_ratio'] for v in frames[1:]]);axes[1,0].axhline(1,color='r',ls='--');axes[1,0].set_ylabel('Bulk-ice / domain speed (%)')
    for key,label in [('salt_budget_relative','Salt'),('tracer_budget_relative','Signed tracer'),('heat_budget_relative','Heat')]:
        vals=np.array([v[key] if v[key] is not None else np.nan for v in brows]);axes[1,1].semilogy(bt,np.where(vals>0,vals,np.nan),label=label)
    axes[1,1].legend();axes[1,1].set_ylabel('Relative budget error')
    for ax in axes.flat:ax.set_xlabel('Nondimensional time');ax.grid(alpha=.2)
    fig.suptitle(f'P1 job {job}: {su:.2f} SU, no release')
    fig.savefig(out/'P1_diagnostics.png',dpi=180);fig.savefig(out/'P1_diagnostics.pdf');plt.close(fig)
    plot_final_slice(out,x,y,slices[len(files)-1],frames[-1]['time'])
    for name in ['analyze_profiles.py','case.json','parties.inp','p_mobile.inp','p_fixed.inp','job.sh','manifest.json','run.sha256']:
        shutil.copy2(run/name,out/name)
    if Path(__file__).resolve() != (out/Path(__file__).name).resolve():
        shutil.copy2(Path(__file__),out/Path(__file__).name)
    print(json.dumps({k:result[k] for k in ['cost_SU','end_time','net_ice_gain','released','rigidity_exceedances','tracer_physical_min','raw_bytes_excluding_symlinks']},indent=2))


if __name__=='__main__':main()
