#!/usr/bin/env python3
"""Reproduce the September 25 P1 physical/EOS/cost preflight; never submit jobs.

Requires numpy, scipy, gsw. All temperatures are Conservative Temperature.
Use exact Gibbs density after CT->in-situ conversion, including diluted water
outside the oceanographic funnel of the efficient density polynomial.
"""
import json
from pathlib import Path
import re

import gsw
import numpy as np
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
B2 = HERE.parents[1] / 'StageB2_closure' / 'B2_closure'


def rho(s, t):
    return gsw.rho_t_exact(s, gsw.t_from_CT(s, t, 0.), 0.)


def metrics(e):
    return dict(max_absolute=float(np.max(np.abs(e))), rms=float(np.sqrt(np.mean(e*e))))


def main():
    d, si, tb, gamma, rhop = .001, 34., 0., .01, 2650.
    nu, kt, cp, latent, ds_phys, sc = 1.8e-6, 1.4e-7, 3991.86795711963, 333500., 1e-9, 70.
    tr = float(gsw.CT_freezing(si, 0., 0.))
    dt = tb-tr
    # Includes production's deeper bottom, complete dilution, and the slightly
    # positive freshwater CT freezing point. Also covers diffuse-band states.
    smax = si+gamma*22*d
    sa, ct = np.meshgrid(np.linspace(0., smax, 201), np.linspace(tr, .03, 101))
    sa, ct = sa.ravel(), ct.ravel()
    density = rho(sa, ct)
    linear = np.linalg.lstsq(np.array([np.ones(sa.size), ct, sa]).T, density, rcond=None)[0]
    # First fit the density-maximum locus; then thermal response; finally
    # intercept and haline amplitude. Density-only LS suppresses thermal
    # expansion by ~55%, despite small error on the total haline scale.
    smd = np.linspace(0., smax, 101)
    tmd = np.array([minimize_scalar(lambda t: -float(rho(s, t)),
                    bounds=(-10., 6.), method='bounded').x for s in smd])
    k, t0 = np.polyfit(smd, tmd, 1)
    eps = 1e-4
    derivative = (rho(sa, ct+eps)-rho(sa, ct-eps))/(2*eps)
    x = -2*(ct-t0-k*sa)
    b = np.dot(x, derivative)/np.dot(x, x)
    a, h = np.linalg.lstsq(np.array([np.ones(sa.size), sa]).T,
                          density+b*(ct-t0-k*sa)**2, rcond=None)[0]
    assert b > 0 and h > 0
    # Independent, denser grid; report local thermal errors separately.
    sv, tv = np.meshgrid(np.linspace(0., smax, 302), np.linspace(tr, .03, 152))
    exact = rho(sv, tv)
    fitted = a-b*(tv-t0-k*sv)**2+h*sv
    drdt = (rho(sv, tv+eps)-rho(sv, tv-eps))/(2*eps)
    drds = (rho(sv+eps, tv)-rho(np.maximum(sv-eps, 0), tv))/(sv+eps-np.maximum(sv-eps, 0))
    thermal_error = -2*b*(tv-t0-k*sv)-drdt
    haline_error = (2*b*k*(tv-t0-k*sv)+h)-drds
    rho0 = float(rho(si, tb))
    # beta is the fitted explicit haline coefficient / rho0. The nonlinear
    # cross term contributes an additional local haline derivative.
    beta = h/rho0
    uref = float(np.sqrt(9.81*beta*si*d))
    tref, re_, pr = d/uref, uref*d/nu, nu/kt
    pe, st = re_*pr, cp*dt/latent
    gstar = 1/(beta*si)
    # Secant passes through both exact freshwater and ambient freezing points.
    af = float(gsw.CT_freezing(0., 0., 0.))
    mf = (af-tr)/si
    freeze_error = af-mf*sv[0]-gsw.CT_freezing(sv[0], 0., 0.)
    eos = dict(eos_q=2., eos_betaT=float(b*dt**2/(h*si)), eos_betaS=1.,
               eos_Tmd0=float((t0-tr)/dt), eos_Tmd_slope=float(k*si/dt))
    assert np.all(2*b*k*(tv-t0-k*sv)+h > 0)
    assert -2*b*(tb-t0-k*si) < 0  # warm ambient seawater is lighter
    assert -2*b*(tb-t0) > 0       # warm fresh water near freezing is denser
    assert pr < sc
    cap, start, tend, fields = .002, .0002, 10., .025
    source_scale = 4*pe*.03125*.0625
    assert cap < .1*source_scale
    # Heat-only estimate: retreat by one diameter, reservoir used first,
    # mean water path 8.5d. Neglects cold-top extraction, dilution, convection.
    reservoir = 8*st
    conduction_t = (1-reservoir)*pe*8.5/st
    steps = tend/cap
    throughput = [31300., 4630.]
    cost = lambda cells, n: [float(cells*n/(3600*q)) for q in throughput]
    record = dict(
        selection_date='2026-09-25', selection='User accepted proposed idealized case',
        physical=dict(d_m=d, rho_p_kg_m3=rhop, CT_bottom_C=tb, SA_interface_g_kg=si,
                      gamma_S_g_kg_m=gamma, p0_dbar=0., CT_ice_C=tr,
                      rho0_kg_m3=rho0, nu_m2_s=nu, kappa_T_m2_s=kt, cp_J_kg_K=cp,
                      latent_J_kg=latent, salt_diffusivity_reference_m2_s=ds_phys,
                      properties_note='Fixed idealized constants, not a site calibration. cp is TEOS-10 cp0 for CT; constant latent heat and unit ice/water diffusivity are model approximations.'),
        eos=dict(choice='Existing EOS_NONLINEAR, q=2; thermal-response fit with pinned density-maximum locus',
                 gsw_version=gsw.__version__, pressure_dbar=0., density_routine='rho_t_exact(SA,t_from_CT(SA,CT,p),p)',
                 CT_range_C=[tr,.03], SA_range_g_kg=[0.,smax],
                 linear_coefficients_rho_c0_cT_cS=linear.tolist(),
                 linear_rejected='Positive d(rho)/dCT gives wrong sign in ambient seawater; a total-density least-squares objective hides local thermal errors.',
                 dimensional_rho_a_b_Tmd0_TmdSlope_h=[float(v) for v in [a,b,t0,k,h]],
                 density_error_kg_m3=metrics(fitted-exact),
                 thermal_derivative_error_kg_m3_K=metrics(thermal_error),
                 haline_derivative_error_kg_m3_per_g_kg=metrics(haline_error),
                 density_maximum_locus_error_K=metrics(t0+k*smd-tmd),
                 ambient_thermal_relative_error_max=float(np.max(np.abs((thermal_error/np.where(np.abs(drdt)>1e-10,drdt,np.nan))[sv>=28]))),
                 coefficients=eos),
        liquidus=dict(a_f_C=af,m_f_K_per_g_kg=mf, error_K=metrics(freeze_error),
                      T_melt=(af-tr)/dt,liquidus_slope=mf*si/dt,
                      note='Secant exact at SA=0 and 34; includes pressure-dependent freshwater CT freezing. Interior residual is reported, not hidden.'),
        scaling=dict(U_ref_m_s=uref,t_ref_s=tref,DeltaT_K=dt,beta_per_g_kg=beta,
                     Re=re_,Pr=pr,Sc_sim=sc,Le=sc/pr,Pe_T=pe,Pe_S=re_*sc,
                     St=st,G_star=gstar,rho_s=rhop/rho0,
                     Ga=float(re_*np.sqrt((rhop/rho0-1)*gstar)),
                     salt_diffusivity_sim_m2_s=nu/sc,
                     N_squared_s2=float(9.81*(h+2*b*k*(tb-t0-k*si))/rho0*gamma)),
        bounded_P1=dict(time_max=tend,physical_duration_s=tend*tref,max_dt=cap,default_dt=start,
                        intended_ramp_duration=.1,release_ramp_steps=50,field_interval=fields,
                        nominal_steps_at_cap=steps, source_relaxation_scale=source_scale,
                        speed_planning_bound=5., advective_cfl_at_bound=5*cap*24,
                        nominal_SU=cost(2211840,steps),SU_with_30pct_contingency=[1.3*c for c in cost(2211840,steps)],
                        allocated_cores=128,wall_hours_hard_cap=3,SU_hard_cap=384,
                        note='Bounded feasibility segment, not an assertion of completed P1 or a post-escape window; stop/checkpoint at 170 min, scheduler limit 180 min. No automatic continuation.',
                        storage_estimate_GiB=250),
        heat_planning=dict(initial_reservoir_melt_depth_d=reservoir,
                           one_d_retreat_conduction_time=conduction_t,
                           one_d_retreat_physical_s=conduction_t*tref,
                           one_d_retreat_SU=cost(2211840,conduction_t/cap),
                           caveat='Order-of-magnitude heat estimate, not release prediction or rigorous cost bound; cold top reduces heat, convection can increase it, dilution changes liquidus. Actual detachment timing is unvalidated.'),
        P2_planning=dict(cells=5242880,max_dt=cap/2,default_dt=start/2,
                         same_duration_SU=cost(5242880,2*steps),authorized=False),
        references=['https://teos-10.github.io/GSW-Python/density.html',
                    'https://teos-10.github.io/GSW-Python/gsw_flat.html#gsw.CT_freezing'])
    record['eos']['density_error_over_haline_scale'] = metrics((fitted-exact)/(h*si))
    (HERE/'case.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    overrides = dict(ymax=10.,NXM=96,NYM=240,NZM=96,time_max=tend,
                     output_time_interval=fields,output_time_interval_2d=fields,
                     default_dt=start,max_dt=cap,Re=re_,vof_slab_x0=.8,Cn=.03125,
                     Pe_CH=28.8,stefan=st,T_melt=(af-tr)/dt,liquidus_slope=mf*si/dt,
                     melt_band_eps=.0625,Pe='{'+', '.join(map(str,[pe,re_*sc,re_*sc]))+'}',
                     cbd2=1-gamma*d/si*2,cbd5=1+gamma*d/si*8,
                     rho_s=rhop/rho0,rho_prImp=rhop/rho0,grav='{0.0, '+str(-gstar)+', 0.0}',
                     release_ramp_steps=50,**eos)
    text = (B2/'coupled16_rampfix/parties.inp').read_text()
    text = '\n'.join(line for line in text.splitlines() if not line.startswith('#'))+'\n'
    for key,value in overrides.items():
        text,n = re.subn(r'^'+re.escape(key)+r'\s*=.*$',key+' = '+str(value),text,flags=re.M)
        assert n == 1,(key,n)
    text = text.replace('[phase_change]','[phase_change]\nyang_salt_transport = 0\nliquid_referenced_salinity = 0')
    (HERE/'parties.inp').write_text('# P1 bounded feasibility segment. See case.json and README.md.\n'+text)
    (HERE/'p_mobile.inp').write_text('1\n2.0 8.5 2.0 0.5\n')
    (HERE/'p_fixed.inp').write_text('0\n')
    (HERE/'stop.inp').write_text('0\n')
    (HERE/'Boundary.scenario.h').write_bytes((B2/'Boundary.validation.h').read_bytes())
    print(json.dumps({k:record[k] for k in ['scaling','bounded_P1','heat_planning']},indent=2))


if __name__ == '__main__':
    main()
