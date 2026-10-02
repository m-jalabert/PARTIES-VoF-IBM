#!/usr/bin/env python3
"""Capture small, read-only extracts from P1b and render supervisor figures.

Default: render the saved extracts, capturing them on the first invocation.
--refresh: replace the extracts with the latest available simulation outputs.
No simulation files are modified. Dependencies: numpy, scipy, matplotlib, h5py.
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch
import numpy as np
from scipy.special import erf

HERE = Path(__file__).resolve().parent
ROOT = Path('/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf')
RUNS = {'grain': 'run_prod_20977443', 'default': 'run_ctrl_20977444',
        'yang': 'run_ctrl_20973341'}
NAVY, TEAL, ORANGE, GREY = '#18334b', '#087f8c', '#c65c2c', '#667888'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
    'axes.titlesize': 12, 'axes.labelsize': 11, 'axes.edgecolor': '#bac6ce',
    'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.color': GREY, 'ytick.color': GREY, 'text.color': NAVY,
    'axes.labelcolor': NAVY, 'figure.facecolor': 'white',
    'savefig.facecolor': 'white', 'pdf.fonttype': 42})


def available(run):
    result = []
    for p in sorted(run.glob('Data_*.h5'), key=lambda p: int(p.stem.split('_')[1])):
        try:
            with h5py.File(p, 'r') as f:
                t = float(f['time'][0])
                if all(key in f for key in ['Conc/0', 'Conc/1', 'Conc/2', 'VOF/C_L', 'VOF/C_S']):
                    result.append((t, p))
        except (OSError, KeyError):
            continue
    return result


def extract(path, centre=True):
    with h5py.File(path, 'r') as f:
        x, y, z = (f['grid/' + a][:-1] for a in ['xc', 'yc', 'zc'])
        k = int(np.argmin(abs(z - 4))) if centre else 0
        out = {'x': x, 'y': y, 'z_slice': z[k], 'time': float(f['time'][0])}
        for name, key in [('theta','Conc/0'), ('salt','Conc/1'), ('tracer','Conc/2'),
                          ('cl','VOF/C_L'), ('cs','VOF/C_S')]:
            # Slice directly from disk; never load a full 19-million-cell field.
            a = f[key][k, :-1, :-1]
            if not np.isfinite(a).all():
                raise ValueError(f'Nonfinite snapshot: {path}: {key}')
            out[name] = a if centre else a[:, 0]
    return out


def capture():
    data = HERE / 'data'
    data.mkdir(exist_ok=True)
    case = json.loads((ROOT / RUNS['grain'] / 'case.json').read_text())
    (data / 'case.json').write_text(json.dumps(case, indent=2) + '\n')
    grain_files = available(ROOT / RUNS['grain'])
    # A running job can be writing the newest file. Fall back if extraction fails.
    for _, p in reversed(grain_files):
        try:
            field = extract(p)
            break
        except (OSError, KeyError, ValueError):
            continue
    else:
        raise RuntimeError('No readable grain snapshot')
    np.savez_compressed(data / 'grain_slice.npz', **field)
    grain_path = p
    defaults, yangs = available(ROOT / RUNS['default']), available(ROOT / RUNS['yang'])
    pairs = [(td, pd, py) for td, pd in defaults for ty, py in yangs
             if td > 0 and abs(td-ty) < 1e-9]
    if not pairs:
        raise RuntimeError('No matched positive-time control snapshots')
    _, pd, py = max(pairs, key=lambda row: row[0])
    np.savez_compressed(data / 'control_default.npz', **extract(pd, False))
    np.savez_compressed(data / 'control_yang.npz', **extract(py, False))
    rows = []
    # Ignore a trailing partial line in an actively appended text file.
    for line in (ROOT / RUNS['grain'] / 'mobile.dat').read_text().splitlines():
        try:
            row = [float(v) for v in line.split(',')]
            if len(row) == 8 and np.isfinite(row).all():
                rows.append(row)
        except ValueError:
            pass
    np.savetxt(data / 'particle.csv', rows, delimiter=',',
               header='time,particle_id,x,y,z,u,v,w', comments='')
    manifest = {'captured_utc': datetime.now(timezone.utc).isoformat(),
        'runs': RUNS, 'grain_field_source': str(grain_path),
        'grain_field_time': float(field['time']),
        'particle_source': str(ROOT / RUNS['grain'] / 'mobile.dat'),
        'particle_last_time': rows[-1][0],
        'default_control_source': str(pd), 'yang_control_source': str(py),
        'control_time': float(extract(pd, False)['time']),
        'note': 'Read-only mid-run snapshot; times are nondimensional. No acceptance verdict.'}
    (data / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


def page(title, subtitle):
    fig = plt.figure(figsize=(13.33, 7.5))
    fig.text(.055, .94, title, fontsize=23, weight='bold')
    fig.text(.055, .895, subtitle, fontsize=11, color=GREY)
    return fig


def foot(fig, text):
    fig.text(.055, .035, text, fontsize=9, color=GREY)


def polish(ax):
    ax.grid(alpha=.15)
    ax.set_axisbelow(True)


def render():
    data = HERE / 'data'
    case = json.loads((data / 'case.json').read_text())
    meta = json.loads((data / 'manifest.json').read_text())
    g = dict(np.load(data / 'grain_slice.npz'))
    controls = {k: dict(np.load(data / ('control_' + k + '.npz'))) for k in ['default', 'yang']}
    particle = np.atleast_2d(np.loadtxt(data / 'particle.csv', delimiter=',', skiprows=1))
    scale, ini, num = case['scaling'], case['initial_state'], case['numerics']
    tref = scale['t_ref_s']; uref = scale['U_ref_m_s']; si = ini['s_interface']
    ellS, ellT = ini['ell_S_d'], ini['ell_T_d']
    figures = []

    fig = page('P1b | A grain settling beneath an ice shelf',
        'Scenario schematic · initial condition · one resolved 0.2 mm quartz sphere · melting active')
    ax = fig.add_axes([.065,.15,.32,.69])
    y = np.linspace(0,22,800); depth = np.maximum(20-y,0)
    c = (1-si)*(1-erf(depth/ellS))
    ax.imshow(np.repeat(c[:,None],120,axis=1), origin='lower',extent=[0,8,0,22],
              cmap='Blues', vmin=0, vmax=.4, aspect='equal')
    ax.add_patch(Rectangle((0,20),8,2,facecolor='#dce8eb',edgecolor=NAVY))
    ax.text(4,21,'Fresh ice',ha='center',va='center',weight='bold')
    ax.axhline(20,color=NAVY,lw=1)
    ax.add_patch(Circle((4,19.15),.5,facecolor=ORANGE,edgecolor='white',lw=1.3))
    ax.annotate('',xy=(4,15.5),xytext=(4,18.45),arrowprops={'arrowstyle':'->','color':ORANGE,'lw':2})
    ax.text(4,9,'Water initially\nat rest',ha='center',color=NAVY)
    ax.set(xlim=(0,8),ylim=(0,22),xlabel=r'$x/d$',ylabel=r'$y/d$  (upward)')
    ax.set_xticks([0,4,8]); ax.set_yticks([0,5,10,15,20,22])
    ax.text(4,-3.0,'Periodic in x and z',ha='center',fontsize=10)
    ax.text(4,1,'No-slip, insulated bottom',ha='center',fontsize=8)
    fig.text(.43,.79,'A local snapshot of the meltwater layer',fontsize=17,weight='bold')
    fig.text(.43,.735,'700 dbar  •  far-field CT = 0 °C  •  SA = 34.85 g/kg',fontsize=12)
    fig.text(.43,.64,'BOX',color=TEAL,fontsize=10,weight='bold')
    fig.text(.43,.594,'8d × 22d × 8d = 1.6 × 4.4 × 1.6 mm',fontsize=15)
    fig.text(.43,.52,'INITIAL STATE',color=TEAL,fontsize=10,weight='bold')
    fig.text(.43,.473,'Ice interface: y/d = 20\nGrain centre: (4, 19.15, 4)d; top clearance: 0.35d\nInterface water: −1.670 °C, 21.30 g/kg',fontsize=12,linespacing=1.65,va='top')
    fig.text(.43,.285,'WHAT WE MEASURE',color=TEAL,fontsize=10,weight='bold')
    fig.text(.43,.24,'Settling speed • meltwater carried downward\nWake evolution • scalar transport • continued melting',fontsize=12,linespacing=1.7,va='top')
    foot(fig,'Schematic is to scale in x–y; shading indicates prescribed meltwater. Release is the initial condition; no background current is imposed.')
    figures.append(('01_domain_and_question',fig))

    fig = page('Initial sublayer | Cold and fresh beneath the ice',
        'Prescribed pure-liquid profiles · depth measured downward from the initial ice interface · no imposed turbulent forcing')
    gs=fig.add_gridspec(1,4,left=.07,right=.97,bottom=.23,top=.80,wspace=.35)
    depth=np.linspace(0,20,1000); th=erf(depth/ellT); sl=si+(1-si)*erf(depth/ellS)
    temp=scale['T_ref_C']+scale['DeltaT_K']*th; salt=34.85*sl; tr=100*(1-sl)
    eos=case['eos']['coefficients']
    b=lambda t,s: -eos['eos_betaT']*(t-eos['eos_Tmd0']-eos['eos_Tmd_slope']*s)**2+s
    drho=(b(th,sl)-b(1.,1.))*case['eos']['dimensional_a_b_Tmd0_TmdSlope_h'][-1]*34.85
    for i,(values,label,color) in enumerate([(temp,'Temperature (°C)',ORANGE),
             (salt,'Liquid salinity (g/kg)',TEAL),(tr,'Meltwater (%)','#416caa'),
             (drho,r'$\rho-\rho_{far}$ (kg/m³)',NAVY)]):
        ax=fig.add_subplot(gs[i]); ax.plot(values,depth,color=color,lw=2.5)
        ax.axhline(case['three_equation']['delta_S_d'],color=TEAL,ls=':',lw=1)
        ax.axhline(case['three_equation']['delta_T_d'],color=ORANGE,ls=':',lw=1)
        ax.set(xlabel=label,ylim=(20,0)); polish(ax)
        if i==0: ax.set_ylabel(r'Depth below ice / $d$')
        else: ax.set_yticklabels([])
    fig.text(.08,.135,'Salt-layer scale: 3.23d',color=TEAL,fontsize=12,weight='bold')
    fig.text(.39,.135,'Thermal-layer scale: 12.73d',color=ORANGE,fontsize=12,weight='bold')
    fig.text(.73,.135,'1d = 0.2 mm',color=GREY,fontsize=12)
    foot(fig,'Analytical initialization, before diffuse-band weighting. Density uses the configured nonlinear EOS; far-field values are references, not wall forcing.')
    figures.append(('02_initial_profiles',fig))

    t=float(g['time']); fig=page('Simulation fields | Latest captured grain snapshot',
        f'Job 20977443 · default salt operator · t = {t:.5f} = {t*tref*1000:.3f} ms · z/d = {float(g["z_slice"]):.4f}')
    gs=fig.add_gridspec(1,3,left=.07,right=.96,bottom=.23,top=.79,wspace=.32)
    excluded=(g['cl']<.999)|(g['cs']>.001)
    arrays=[(scale['T_ref_C']+scale['DeltaT_K']*g['theta'],'Temperature (°C)','magma',-1.67,0),
        (34.85*g['salt'],'Salinity in pure liquid (g/kg)','viridis',21.3,34.85),
        (100*g['tracer'],'Meltwater tracer in pure liquid (%)','Blues',0,40)]
    for i,(a,label,cmap,vmin,vmax) in enumerate(arrays):
        ax=fig.add_subplot(gs[i]); ax.set_facecolor('#e4e8eb')
        pcm=ax.pcolormesh(g['x'],g['y'],np.ma.array(a,mask=excluded),shading='nearest',cmap=cmap,vmin=vmin,vmax=vmax,rasterized=True)
        ax.contour(g['x'],g['y'],g['cs'],levels=[.5],colors=[ORANGE],linewidths=1.5)
        ax.axhline(20,color=NAVY,ls='--',lw=.8)
        ax.set(xlim=(2.5,5.5),ylim=(16,20.5),xlabel=r'$x/d$'); ax.set_aspect('equal')
        if i==0: ax.set_ylabel(r'$y/d$')
        cb=fig.colorbar(pcm,ax=ax,orientation='horizontal',pad=.13,fraction=.045)
        cb.set_label(label,fontsize=9)
    fig.text(.08,.125,('Initial output only: a developed wake is not yet available.' if t==0 else 'An early transient snapshot; this is not a converged transport result.'),fontsize=12,weight='bold')
    foot(fig,r'Actual HDF5 slices. Grey excludes ice, grain and diffuse support ($C_L<0.999$ or $C_S>0.001$); orange outlines $C_S=0.5$.')
    figures.append(('03_grain_field_snapshot',fig))

    fig=page('Grain motion | Initial acceleration recorded so far',
        f'Job 20977443 · actual particle output · captured through t = {particle[-1,0]:.5f} ({particle[-1,0]*tref*1000:.3f} ms)')
    gs=fig.add_gridspec(1,2,left=.085,right=.95,bottom=.28,top=.79,wspace=.30)
    tm=particle[:,0]*tref*1000
    ax=fig.add_subplot(gs[0]);ax.plot(tm,(particle[0,3]-particle[:,3])*.2*1000,color=TEAL,lw=2.4)
    ax.set(xlabel='Elapsed physical time (ms)',ylabel='Downward displacement (µm)');polish(ax)
    ax=fig.add_subplot(gs[1]);ax.plot(tm,-particle[:,6]*uref*1000,color=ORANGE,lw=2.4)
    ax.set(xlabel='Elapsed physical time (ms)',ylabel='Downward grain speed (mm/s)');polish(ax)
    fig.text(.085,.17,'Planned observation window: 821 ms',fontsize=13,weight='bold')
    fig.text(.085,.12,'Only the recorded interval is plotted. These curves do not establish terminal settling or precursor acceptance.',fontsize=11,color=GREY)
    foot(fig,'Source: mobile.dat, columns time / y / vertical velocity. Downward displacement and speed are plotted as positive.')
    figures.append(('04_early_grain_motion',fig))

    fig=page('From precursor to larger-domain production',
        'P1b establishes operability and signal size; P2 screens resolution sensitivity; production geometry remains to be selected')
    ax=fig.add_axes([.04,.13,.92,.69]);ax.set_axis_off()
    boxes=[(.015,'P1b + control','24 cells / d\n192 × 528 × 192\n8d × 22d × 8d','CURRENT'),
           (.355,'P2 refinement','32 cells / d\nSame physical scenario\nFiner time and interface','PLANNED'),
           (.695,'Larger domain','More lateral separation\nLonger useful fall window\nGeometry and cost: pending','CONDITIONAL')]
    for x,title,body,status in boxes:
        ax.add_patch(FancyBboxPatch((x,.47),.285,.44,boxstyle='round,pad=.012,rounding_size=.015',facecolor='#eef5f7',edgecolor='#c1d6de'))
        ax.text(x+.018,.85,status,color=TEAL,fontsize=9,weight='bold')
        ax.text(x+.018,.765,title,fontsize=16,weight='bold')
        ax.text(x+.018,.66,body,fontsize=11,va='top',linespacing=1.8)
    for x in [.307,.647]: ax.annotate('',xy=(x+.04,.68),xytext=(x,.68),arrowprops={'arrowstyle':'->','lw':2,'color':TEAL})
    ax.text(.02,.33,'HEALTH & ATTRIBUTION',color=TEAL,fontsize=10,weight='bold')
    ax.text(.02,.265,'Bounded fields; budgets\nGrain minus matched control\nOperator sensitivity',fontsize=11,va='top',linespacing=1.7)
    ax.text(.36,.33,'PROPOSED SCREENING TARGETS',color=TEAL,fontsize=10,weight='bold')
    ax.text(.36,.265,'≤ 2% change in settling speed\n≤ 5% in selected transport metrics\nCompare at matched physical states',fontsize=11,va='top',linespacing=1.7)
    ax.text(.70,.33,'PRODUCTION DECISION',color=TEAL,fontsize=10,weight='bold')
    ax.text(.70,.265,'Accept resolution and model limits\nReassess walls and periodic images\nConfirm runtime and storage',fontsize=11,va='top',linespacing=1.7)
    foot(fig,'Acceptance remains open. Grid agreement cannot remove interface-model error or the selected Sc = 70 diffusivity approximation.')
    figures.append(('05_precursor_to_production',fig))

    ct=float(controls['default']['time']);fig=page('Ice-only controls | Early salt-operator comparison',
        f'Matched saved time t = {ct:.6f} = {ct*tref*1000:.3f} ms · default job 20977444 vs Yang job 20973341')
    gs=fig.add_gridspec(1,3,left=.075,right=.965,bottom=.28,top=.79,wspace=.33)
    specs=[('theta','Temperature (°C)'),('salt',r'Salt inventory per cell volume ($S_{ref}$ units)'),('tracer','Meltwater tracer (%)')]
    for i,(key,label) in enumerate(specs):
        ax=fig.add_subplot(gs[i])
        for name,color in [('default',ORANGE),('yang',TEAL)]:
            q=controls[name]; dep=20-q['y']; a=q[key].copy()
            if key=='theta': a=scale['T_ref_C']+scale['DeltaT_K']*a
            if key=='salt' and name=='yang': a=q['cl']*a
            if key=='tracer': a=100*a
            selected=(dep>=-.3)&(dep<=1.5)
            ax.plot(a[selected],dep[selected],color=color,lw=2.1,label=name.capitalize())
        ax.axhline(0,color=GREY,ls=':',lw=1)
        ax.set(ylim=(1.5,-.3),xlabel=label);polish(ax)
        if i==0: ax.set_ylabel(r'Depth below initial interface / $d$');ax.legend(frameon=False,loc='lower left')
        else: ax.set_yticklabels([])
    fig.text(.075,.17,'The comparison isolates background-layer sensitivity to the salt operator.',fontsize=12,weight='bold')
    fig.text(.075,.115,'It does not measure sensitivity of the grain–layer interaction. This early snapshot is not an acceptance test.',fontsize=11,color=GREY)
    foot(fig,r'Actual profiles at x,z near the periodic corner. Salt compares default $s$ with Yang $C_L S$; it is not interface liquid salinity.')
    figures.append(('06_matched_control_profiles',fig))

    with PdfPages(HERE/'P1b_supervisor_figures.pdf') as pdf:
        pdf.infodict()['Title']='P1b: setup, captured fields and path to production'
        pdf.infodict()['Subject']='Read-only extracts from active simulations; preliminary results'
        for name,fig in figures:
            fig.savefig(HERE/(name+'.png'),dpi=200)
            fig.savefig(HERE/(name+'.pdf'))
            pdf.savefig(fig)
            plt.close(fig)
    from PIL import Image
    sheet=Image.new('RGB',(1600,1380),'#dde4e8')
    for i,(name,_) in enumerate(figures):
        with Image.open(HERE/(name+'.png')) as original:
            thumb=original.convert('RGB')
        thumb.thumbnail((790,444))
        sheet.paste(thumb,((i%2)*800+5,(i//2)*460+5))
    sheet.save(HERE/'figure_contact_sheet.png')
    print(json.dumps(meta,indent=2))
    print(f'Saved {len(figures)} PNG/PDF figures and combined PDF to {HERE}')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refresh',action='store_true')
    args=parser.parse_args()
    if args.refresh or not (HERE/'data/manifest.json').exists(): capture()
    render()
